#include "Modules/Hydrodynamics/HydroFittedMomentum.h"
#include "Modules/Hydrodynamics/Detail/HydroMomentumAssembly.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/Detail/HydroFittedBoundaryMomentum.h"
#include "Modules/Hydrodynamics/Detail/HydroNativeRZRowPreparation.h"
#include "Modules/Hydrodynamics/HydroBoundary.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
using Type = HydroFittedDomain::BoundaryType;
using Face = HydroFittedDomain::NativeFace;
using Triplet = Eigen::Triplet<double>;

void require(bool value, const char* reason)
{
    if (!value) throw std::invalid_argument(std::string("[HydroFittedMomentum] ") + reason);
}

bool pressure(Type type)
{
    return type == Type::PressureReservoir || type == Type::PressureOutlet;
}

class Assembly
{
public:
    Assembly(const QuadtreeMesh& mesh, const HydroFittedDomain& domain, double dt,
        const CellField<double>& rho, const CellField<double>& mu,
        const VertField<double>& vertexRho, const VertField<double>& vertexMu,
        const Eigen::VectorXd& current, const Eigen::VectorXd& previous,
        const CellField<double>& bodyX, const CellField<double>& bodyY,
        const std::vector<double>& tolerances, const HydroFittedNativeInterior* nativeInterior)
        : mesh_(mesh), domain_(domain), dt_(dt), rho_(rho), mu_(mu),
          vertexRho_(vertexRho), vertexMu_(vertexMu), current_(current),
          previous_(previous), bodyX_(bodyX), bodyY_(bodyY), tolerances_(tolerances),
          nr_(mesh.vEdgeIndexNum()), count_(nr_ + mesh.hEdgeIndexNum())
          , nativeInterior_(nativeInterior)
    {
        domain_.validateMeshIdentity(mesh_);
        require(domain_.axisymmetric(), "this native fitted momentum assembly requires RZ geometry");
        require(std::isfinite(dt_) && dt_ > 0.0, "dt must be finite and positive");
        require(static_cast<int>(domain_.faces().size()) == count_
            && static_cast<int>(domain_.velocityRowByNative().size()) == count_,
            "domain does not cover the complete native index space");
        require(current_.size() == count_ && previous_.size() == count_,
            "velocity histories must use complete native indexing");
        require(rho_.size() == mesh_.cellIndexNum() && mu_.size() == mesh_.cellIndexNum()
            && bodyX_.size() == mesh_.cellIndexNum() && bodyY_.size() == mesh_.cellIndexNum()
            && vertexRho_.size() == mesh_.vertIndexNum() && vertexMu_.size() == mesh_.vertIndexNum()
            && rho_.value && mu_.value && bodyX_.value && bodyY_.value
            && vertexRho_.value && vertexMu_.value, "field storage does not match the mesh");
        const int freeRows = static_cast<int>(domain_.nativeByVelocityRow().size());
        require(domain_.uniformLeaves() || nativeInterior_, "mixed-level momentum requires prepared native interior rows");
        if (nativeInterior_)
        {
            require(nativeInterior_->layout != nullptr
                && nativeInterior_->matrix.rows()==count_ && nativeInterior_->matrix.cols()==count_
                && nativeInterior_->rhs.size()==count_ && nativeInterior_->density.size()==count_
                && nativeInterior_->viscosity.size()==count_, "invalid prepared native interior layout");
            nativeInterior_->layout->validate(mesh_, domain_);
        }
        require(tolerances_.empty() || static_cast<int>(tolerances_.size()) == freeRows,
            "outlet tolerances must be empty or indexed by compact free row");
        for (double tolerance : tolerances_)
            require(std::isfinite(tolerance) && tolerance >= 0.0, "invalid outlet backflow tolerance");
        assembly_.rhs = Eigen::VectorXd::Zero(freeRows);
        result_.densityByVelocityRow.resize(freeRows);
        assembly_.entries.reserve(12 * freeRows);
        // Never inspect a solid-only field/history value, even for validation.
        for (int index : domain_.cellIndexByRow())
        {
            require(std::isfinite(rho_.value[index]) && rho_.value[index] > 0.0
                && std::isfinite(mu_.value[index]) && mu_.value[index] >= 0.0,
                "active-fluid cell density/viscosity is invalid");
            require(std::isfinite(bodyX_.value[index]) && std::isfinite(bodyY_.value[index]),
                "active-fluid body force is nonfinite");
        }
        for (const auto& f : domain_.faces())
        {
            if (!active(f)) continue;
            const double u = value(f.nativeRow, false), pre = value(f.nativeRow, true);
            if (f.boundary.type == Type::PressureReservoir)
                require(std::isfinite(f.boundary.incomingDensity) && f.boundary.incomingDensity > 0.0
                    && std::isfinite(f.boundary.incomingViscosity) && f.boundary.incomingViscosity >= 0.0,
                    "reservoir incoming material is invalid");
            if (pressure(f.boundary.type))
                require(std::isfinite(f.boundary.pressure), "pressure traction is nonfinite");
            if (f.boundary.type == Type::PressureOutlet)
            {
                const double tolerance = tolerances_.empty() ? 0.0 : tolerances_[f.velocityRow];
                if (f.outwardSign * u < -tolerance || f.outwardSign * pre < -tolerance)
                {
                    std::ostringstream message;
                    message << std::setprecision(17) << "[HydroFittedMomentum] pressure outlet backflow"
                        << " requires incoming state; nativeRow=" << f.nativeRow
                        << " outward=" << f.outwardSign * u << " previous=" << f.outwardSign * pre
                        << " tolerance=" << tolerance;
                    throw std::runtime_error(message.str());
                }
            }
        }
    }

    HydroFittedMomentumSystem run()
    {
        const auto* boundary = nativeInterior_ ? nativeInterior_->layout->boundary() : nullptr;
        if (boundary) boundary->append(mesh_, domain_, dt_, rho_, mu_, vertexRho_, vertexMu_,
            current_, previous_, bodyX_, bodyY_, result_, assembly_);
        for (int native : domain_.nativeByVelocityRow())
        {
            if (boundary && boundary->contains(native)) continue;
            const auto& f = face(native);
            HydroNativeRZRow::Row row;
            const bool adaptiveNative = nativeInterior_ && nativeInterior_->layout->nativeInteriorRows()[native];
            const bool regularNative = !adaptiveNative && regular(f, row);
            if (adaptiveNative)
            {
                for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator entry(nativeInterior_->matrix,native); entry; ++entry)
                    add(native,static_cast<int>(entry.col()),entry.value());
                assembly_.rhs[f.velocityRow]=nativeInterior_->rhs[native];
                const double density=nativeInterior_->density[native];
                require(std::isfinite(density) && density>0.0, "invalid native interpolated face density");
                result_.densityByVelocityRow[f.velocityRow]=density;
                ++result_.adaptiveNativeRows;
            }
            else if (regularNative)
            {
                const auto emit = [&](int r, int c, double a) { add(r, c, a); };
                const auto now = HydroNativeRZRow::values(row.indices, [&](int k) { return value(k, false); });
                const auto pre = HydroNativeRZRow::values(row.indices, [&](int k) { return value(k, true); });
                if (f.horizontal)
                {
                    HydroNativeRZRow::appendAxialMatrix(row, dt_, emit);
                    assembly_.rhs[f.velocityRow] = HydroNativeRZRow::axialRhs(row, now, pre,
                        bodyY_.value[f.positiveCell], bodyY_.value[f.negativeCell], dt_);
                }
                else
                {
                    HydroNativeRZRow::appendRadialMatrix(row, dt_, emit);
                    assembly_.rhs[f.velocityRow] = HydroNativeRZRow::radialRhs(row, now, pre,
                        bodyX_.value[f.positiveCell], bodyX_.value[f.negativeCell], dt_);
                }
                result_.densityByVelocityRow[f.velocityRow] = row.materials.rho;
                ++result_.regularNativeRows;
            }
            else
            {
                pieces(f);
                ++result_.boundaryPieceRows;
            }
            // Preserve the original full-native volume arithmetic on regular
            // rows; boundary pieces use the view's actual dual bounds.
            if (!f.horizontal) hoop(f, regularNative ? row.metrics.volume : f.dualVolume);
        }
        require(assembly_.rhs.allFinite(), "assembled compact momentum RHS is nonfinite");
        result_.rhs = std::move(assembly_).finish(result_.matrix);
        for (int k = 0; k < result_.matrix.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(result_.matrix, k); it; ++it)
                require(std::isfinite(it.value()), "assembled compact momentum matrix is nonfinite");
        return std::move(result_);
    }

private:
    const QuadtreeMesh& mesh_;
    const HydroFittedDomain& domain_;
    double dt_;
    const CellField<double>& rho_;
    const CellField<double>& mu_;
    const VertField<double>& vertexRho_;
    const VertField<double>& vertexMu_;
    const Eigen::VectorXd& current_;
    const Eigen::VectorXd& previous_;
    const CellField<double>& bodyX_;
    const CellField<double>& bodyY_;
    const std::vector<double>& tolerances_;
    int nr_, count_;
    const HydroFittedNativeInterior* nativeInterior_;
    HydroFittedMomentumSystem result_;
    HydroMomentumAssembly assembly_;

    static bool active(const Face& f) { return f.negativeCell >= 0 || f.positiveCell >= 0; }
    const Face& face(int native) const
    {
        require(native >= 0 && native < count_, "missing native face in a required fluid stencil");
        const auto& f = domain_.faces()[native];
        require(f.nativeRow == native, "domain native-face index is inconsistent");
        return f;
    }
    const Cell& cell(int index) const
    {
        require(index >= 0 && index < mesh_.cellIndexNum()
            && domain_.cellRowByIndex()[index] >= 0, "attempt to access a solid-only cell");
        return mesh_.cells(mesh_.cellIndexToLabel[index]);
    }
    const Edge& edge(const Face& f) const
    {
        return f.horizontal ? mesh_.hEdges(f.edgeLabel) : mesh_.vEdges(f.edgeLabel);
    }
    int native(bool horizontal, int label) const
    {
        require(label >= 0, "required native neighbor is absent");
        return (horizontal ? nr_ : 0)
            + (horizontal ? mesh_.hEdges(label).index() : mesh_.vEdges(label).index());
    }
    int normalFace(const Face& f, int owner, int side) const
    {
        const auto& c = cell(owner);
        return native(f.horizontal, f.horizontal ? c.cellToHEdge[side][0] : c.cellToVEdge[side][0]);
    }
    int transverseFace(const Face& f, int owner, int side) const
    {
        const auto& c = cell(owner);
        return native(!f.horizontal, f.horizontal ? c.cellToVEdge[side][0] : c.cellToHEdge[side][0]);
    }
    int parallelFace(const Face& f, int side) const
    {
        const auto& v = mesh_.verts(edge(f).edgeToVert[side]);
        return native(f.horizontal, f.horizontal ? v.vertToHEdge[side] : v.vertToVEdge[side]);
    }
    double value(int nativeRow, bool previous) const
    {
        const auto& f = face(nativeRow);
        require(active(f), "solid-only native velocity has no physical trace");
        const double u = previous ? previous_[nativeRow] : current_[nativeRow];
        require(std::isfinite(u), "active native velocity/history is nonfinite");
        if (f.velocityRow < 0)
            require((f.boundary.type == Type::NoSlipWall || f.boundary.type == Type::Symmetry) && u == 0.0,
                "wall/symmetry native velocity must be exactly zero");
        return u;
    }
    double vertexValue(const VertField<double>& field, int index, bool positive) const
    {
        const double v = field.value[index];
        require(std::isfinite(v) && (positive ? v > 0.0 : v >= 0.0),
            "active-fluid vertex material is invalid");
        return v;
    }
    void add(int r, int c, double a)
    {
        const auto& target = face(r);
        const auto& source = face(c);
        if (target.velocityRow < 0 || !active(source) || !std::isfinite(a))
        {
            std::ostringstream message;
            message << std::setprecision(17)
                << "[HydroFittedMomentum] invalid momentum coefficient or solid-only column"
                << "; nativeRow=" << r << ", nativeColumn=" << c << ", coefficient=" << a
                << "; targetFace=(" << target.center.x << ',' << target.center.y << ')'
                << ", sourceFace=(" << source.center.x << ',' << source.center.y << ')'
                << "; targetVelocityRow=" << target.velocityRow
                << ", sourceCells=(" << source.negativeCell << ',' << source.positiveCell << ')'
                << ", adaptiveRow=" << bool(nativeInterior_ && nativeInterior_->layout->nativeInteriorRows()[r]);
            throw std::invalid_argument(message.str());
        }
        // Exact homogeneous constraints: no ghost value or solid history read.
        if (source.velocityRow >= 0) assembly_.entries.emplace_back(target.velocityRow, source.velocityRow, a);
        else require(source.boundary.type == Type::NoSlipWall || source.boundary.type == Type::Symmetry,
            "unsupported eliminated native constraint");
    }

    bool regular(const Face& f, HydroNativeRZRow::Row& row) const
    {
        if (f.negativeCell < 0 || f.positiveCell < 0 || f.boundary.type != Type::Interior) return false;
        const auto& e = edge(f);
        const auto& start = mesh_.verts(e.edgeToVert[0]);
        const auto& end = mesh_.verts(e.edgeToVert[1]);
        const int lowLabel = f.horizontal ? start.vertToHEdge[0] : start.vertToVEdge[0];
        const int highLabel = f.horizontal ? end.vertToHEdge[1] : end.vertToVEdge[1];
        if (lowLabel < 0 || highLabel < 0) return false;
        auto& i = row.indices;
        i.center = f.nativeRow;
        i.west = f.horizontal ? native(true, lowLabel) : normalFace(f, f.negativeCell, 0);
        i.east = f.horizontal ? native(true, highLabel) : normalFace(f, f.positiveCell, 1);
        i.south = f.horizontal ? normalFace(f, f.negativeCell, 0) : native(false, lowLabel);
        i.north = f.horizontal ? normalFace(f, f.positiveCell, 1) : native(false, highLabel);
        i.crossSW = transverseFace(f, f.negativeCell, 0);
        i.crossSE = transverseFace(f, f.horizontal ? f.negativeCell : f.positiveCell, f.horizontal ? 1 : 0);
        i.crossNW = transverseFace(f, f.horizontal ? f.positiveCell : f.negativeCell, f.horizontal ? 0 : 1);
        i.crossNE = transverseFace(f, f.positiveCell, 1);
        for (int k : {i.west, i.east, i.south, i.north, i.crossSW, i.crossSE, i.crossNW, i.crossNE})
        {
            const auto& neighbor = face(k);
            if (neighbor.boundary.type != Type::Interior || neighbor.negativeCell < 0
                || neighbor.positiveCell < 0 || neighbor.velocityRow < 0) return false;
        }
        const auto& master = cell(f.positiveCell);
        const auto& slave = cell(f.negativeCell);
        HydroNativeRZRowDetail::prepareRegularData(row, !f.horizontal, e, master, slave, start, end, rho_, mu_,
            f.horizontal ? f.center.x - face(i.west).center.x : f.center.y - face(i.south).center.y,
            f.horizontal ? face(i.east).center.x - f.center.x : face(i.north).center.y - f.center.y,
            [&](int index) { return vertexValue(vertexMu_, index, false); });
        HydroNativeRZRowDetail::prepareTransverseDensity(row.materials, !f.horizontal,
            start.index(), end.index(), [&](int index) { return vertexValue(vertexRho_, index, true); });
        return true;
    }

    double boundaryDensity(const Face& boundary, int owner, double outward) const
    {
        return boundary.boundary.type == Type::PressureReservoir && outward < 0.0
            ? boundary.boundary.incomingDensity : rho_.value[owner];
    }

    double incomingTangent(const Face& f, int owner, int tangentSide, bool previous) const
    {
        const double sign = tangentSide ? 1.0 : -1.0;
        double derivative = 0.0;
        if (f.negativeCell >= 0 && f.positiveCell >= 0)
            derivative = (value(transverseFace(f, f.positiveCell, tangentSide), previous)
                - value(transverseFace(f, f.negativeCell, tangentSide), previous)) / f.normalSpacing;
        else
        {
            // At a physical pressure/pressure corner use TWO actual inward
            // native boundary faces, as in the existing corner zero-shear trace.
            const int inner = normalFace(f, owner, f.outwardSign > 0.0 ? 0 : 1);
            const auto& innerFace = face(inner);
            const int innerOwner = f.outwardSign > 0.0 ? innerFace.negativeCell : innerFace.positiveCell;
            require(innerOwner >= 0, "pressure-corner trace lacks its second inward fluid cell");
            const int here = transverseFace(f, owner, tangentSide);
            const int inside = transverseFace(f, innerOwner, tangentSide);
            const auto& a = face(here); const auto& b = face(inside);
            require(active(b) && (f.horizontal ? a.center.x == b.center.x : a.center.y == b.center.y),
                "pressure-corner trace leaves its physical boundary plane");
            derivative = f.outwardSign * (value(here, previous) - value(inside, previous)) / f.normalSpacing;
        }
        return HydroBoundaryDetail::zeroShearTrace(value(f.nativeRow, previous),
            f.tangentSpacing, sign, derivative);
    }

    bool wholeFluidSide(const Face& f, int tangentSide, double* flux, double* fluxPre)
    {
        // Decide from this actual shared side, NOT from a row-wide near-wall
        // label. A complete free/free dual side must have the same flux on
        // BOTH rows, including when its neighbor uses the far-native kernel.
        if (f.negativeCell < 0 || f.positiveCell < 0) return false;
        const int high = transverseFace(f, f.positiveCell, tangentSide);
        const int low = transverseFace(f, f.negativeCell, tangentSide);
        for (int k : {high, low})
        {
            const auto& adjacent = face(k);
            if (adjacent.boundary.type != Type::Interior || adjacent.negativeCell < 0
                || adjacent.positiveCell < 0 || adjacent.velocityRow < 0) return false;
        }
        const int neighbor = parallelFace(f, tangentSide);
        const auto& adjacent = face(neighbor);
        if (adjacent.boundary.type != Type::Interior || adjacent.negativeCell < 0
            || adjacent.positiveCell < 0 || adjacent.velocityRow < 0) return false;
        const auto& vertex = mesh_.verts(edge(f).edgeToVert[tangentSide]);
        const auto& master = cell(f.positiveCell);
        const double area = f.horizontal ? 2.0 * Pi * vertex.coord().x * master.size().y
            : 2.0 * Pi * f.center.x * master.size().x;
        const double sign = tangentSide ? 1.0 : -1.0;
        const int slot = f.horizontal ? tangentSide : 2 + tangentSide;
        const double density = vertexValue(vertexRho_, vertex.index(), true);
        const double viscosity = vertexValue(vertexMu_, vertex.index(), false);
        const double v = 0.5 * (value(high, false) + value(low, false));
        const double vPre = 0.5 * (value(high, true) + value(low, true));
        const double q = density * v * area;
        const double qPre = density * vPre * area;
        const double u = value(f.nativeRow, false), pre = value(f.nativeRow, true);
        const double un = value(neighbor, false), unPre = value(neighbor, true);
        flux[slot] += q * (sign * q > 0.0 ? u : un);
        fluxPre[slot] += qPre * (sign * qPre > 0.0 ? pre : unPre);
        const double areaDivSpacing = area / f.tangentSpacing;
        const double coefficient = 0.5 * viscosity * areaDivSpacing;
        const double crossCoefficient = 0.5 * viscosity * area / f.normalSpacing;
        add(f.nativeRow, f.nativeRow, coefficient);
        add(f.nativeRow, neighbor, -coefficient);
        add(f.nativeRow, high, -sign * crossCoefficient);
        add(f.nativeRow, low, sign * crossCoefficient);
        const double distance = f.horizontal
            ? sign * (adjacent.center.x - f.center.x) : sign * (adjacent.center.y - f.center.y);
        require(distance > 0.0 && std::isfinite(distance), "invalid shared native shear distance");
        assembly_.rhs[f.velocityRow] += 0.5 * viscosity * area * ((un - u) / distance
            + sign * (value(high, false) - value(low, false)) / f.normalSpacing);
        return true;
    }

    void pieces(const Face& f)
    {
        const int r = f.nativeRow, compact = f.velocityRow;
        const double u = value(r, false), pre = value(r, true);
        const bool full = f.negativeCell >= 0 && f.positiveCell >= 0;
        const int onlyOwner = f.negativeCell >= 0 ? f.negativeCell : f.positiveCell;
        require(full || pressure(f.boundary.type), "free boundary row lacks full pressure traction");
        const double density = full ? 0.5 * (rho_.value[f.positiveCell] + rho_.value[f.negativeCell])
            : rho_.value[onlyOwner];
        result_.densityByVelocityRow[compact] = density;
        const double mass = f.dualVolume * density / dt_;
        add(r, r, mass);
        const auto& body = f.horizontal ? bodyY_ : bodyX_;
        assembly_.rhs[compact] = mass * u + (full
            ? f.dualVolume * 0.5 * (body.value[f.positiveCell] + body.value[f.negativeCell])
            : f.dualVolume * body.value[onlyOwner]);
        // Positive-coordinate flux arrays: W,E,S,N. Aggregate actual side
        // pieces, then retain the E-W+N-S AB2 order (never use a wall ghost donor).
        double flux[4] = {}, fluxPre[4] = {};
        bool wholeSide[2] = {};
        for (int side = 0; side < 2; ++side)
            wholeSide[side] = wholeFluidSide(f, side, flux, fluxPre);
        for (int normalSide = 0; normalSide < 2; ++normalSide)
        {
            const int owner = normalSide ? f.positiveCell : f.negativeCell;
            if (owner < 0) continue;
            const auto& c = cell(owner);
            const double normalSign = normalSide ? 1.0 : -1.0;
            const int inner = normalFace(f, owner, normalSide);
            const double ui = value(inner, false), uiPre = value(inner, true);
            const double innerArea = f.horizontal ? 2.0 * Pi * f.center.x * c.size().x
                : 2.0 * Pi * c.centerCoord().x * c.size().y;
            const double viscous = mu_.value[owner] * innerArea / f.normalSpacing;
            add(r, r, viscous); add(r, inner, -viscous);
            assembly_.rhs[compact] -= viscous * (u - ui);
            const int normalSlot = f.horizontal ? 2 + normalSide : normalSide;
            const double q = rho_.value[owner] * 0.5 * (u + ui) * innerArea;
            const double qPre = rho_.value[owner] * 0.5 * (pre + uiPre) * innerArea;
            flux[normalSlot] += q * (normalSign * q > 0.0 ? u : ui);
            fluxPre[normalSlot] += qPre * (normalSign * qPre > 0.0 ? pre : uiPre);

            for (int tangentSide = 0; tangentSide < 2; ++tangentSide)
            {
                if (wholeSide[tangentSide]) continue;
                const double sign = tangentSide ? 1.0 : -1.0;
                const int slot = f.horizontal ? tangentSide : 2 + tangentSide;
                const auto& vertex = mesh_.verts(edge(f).edgeToVert[tangentSide]);
                // Each owner contributes only its REAL half of this dual side.
                // A radial row's horizontal piece is an exact annular area;
                // an axial row's radial piece has its actual half axial height.
                const double a = std::min(f.center.x, c.centerCoord().x);
                const double b = std::max(f.center.x, c.centerCoord().x);
                const double area = f.horizontal
                    ? 2.0 * Pi * vertex.coord().x * std::abs(f.center.y - c.centerCoord().y)
                    : Pi * (b - a) * (b + a);
                require(std::isfinite(area) && area >= 0.0, "invalid physical dual-side piece area");
                const int other = transverseFace(f, owner, tangentSide);
                const auto& boundary = face(other);
                require(active(boundary), "fluid half-dual meets a solid-only transverse face");
                if (boundary.boundary.type == Type::NoSlipWall)
                {
                    // True no-slip distance h_t/2. The factor 1/2 is CN,
                    // not a geometric doubling/odd ghost used by advection.
                    const double wallCoefficient = 0.5 * mu_.value[owner] * area / (0.5 * f.tangentSpacing);
                    add(r, r, wallCoefficient);
                    assembly_.rhs[compact] -= wallCoefficient * u;
                    continue; // Both history advective wall fluxes are exactly zero.
                }
                if (boundary.boundary.type == Type::Symmetry)
                {
                    // Stationary slip/symmetry: its normal native constraint
                    // is exactly zero and tau_nt=0, at any physical area.
                    // The radial axis is the additional zero-measure case.
                    continue;
                }
                if (pressure(boundary.boundary.type))
                {
                    // Prescribed full traction has zero tangential component.
                    // Its incoming transported tangential velocity still has
                    // the explicit zero-shear trace; it is NOT zero momentum.
                    for (int history = 0; history < 2; ++history)
                    {
                        const bool old = history != 0;
                        const double v = value(other, old);
                        const double outward = sign * v;
                        const double densityAtBoundary = boundaryDensity(boundary, owner, outward);
                        const double donor = outward < 0.0 ? incomingTangent(f, owner, tangentSide, old)
                            : value(r, old);
                        (old ? fluxPre : flux)[slot] += densityAtBoundary * v * area * donor;
                    }
                    continue;
                }
                require(boundary.boundary.type == Type::Interior
                    && boundary.negativeCell >= 0 && boundary.positiveCell >= 0,
                    "unclassified transverse fluid side piece");
                const int neighbor = parallelFace(f, tangentSide);
                const auto& parallel = face(neighbor);
                // The half-dual zero-shear continuation below is paired only
                // with another pressure half dual (or a real constrained wall
                // trace). At a pressure/mask staircase a full free neighbor
                // would retain this shared shear while this half discarded it.
                // That unmatched traction needs a separate junction closure.
                require(full || parallel.velocityRow < 0 || parallel.negativeCell < 0
                    || parallel.positiveCell < 0,
                    "pressure half/full transverse junction has no conservative shear closure");
                const double un = value(neighbor, false), unPre = value(neighbor, true);
                const double vertexDensity = vertexValue(vertexRho_, vertex.index(), true);
                const double v = value(other, false), vPre = value(other, true);
                // Boundary-local closure: each primal face is constant on its
                // own half-piece. Never average a normal flux across a wall.
                const double qSide = vertexDensity * v * area;
                const double qSidePre = vertexDensity * vPre * area;
                flux[slot] += qSide * (sign * qSide > 0.0 ? u : un);
                fluxPre[slot] += qSidePre * (sign * qSidePre > 0.0 ? pre : unPre);
                if (!full)
                {
                    // A pressure-normal row has only an inward half dual.
                    // Linear zero-shear extension of the OTHER component
                    // gives d_n(u_t)=-d_t(u_n), so these two terms in tau_nt
                    // cancel on this FLUID piece. This does not cancel a
                    // no-slip piece, handled above with its true h_t/2.
                    continue;
                }
                const double viscosity = vertexValue(vertexMu_, vertex.index(), false);
                const int high = transverseFace(f, f.positiveCell, tangentSide);
                const int low = transverseFace(f, f.negativeCell, tangentSide);
                const double coefficient = 0.5 * viscosity * area / f.tangentSpacing;
                const double crossCoefficient = 0.5 * sign * viscosity * area / f.normalSpacing;
                add(r, r, coefficient); add(r, neighbor, -coefficient);
                add(r, high, -crossCoefficient); add(r, low, crossCoefficient);
                assembly_.rhs[compact] += 0.5 * viscosity * area * ((un - u) / f.tangentSpacing
                    + sign * (value(high, false) - value(low, false)) / f.normalSpacing);
            }
        }
        if (!full)
        {
            const int inner = normalFace(f, onlyOwner, f.outwardSign > 0.0 ? 0 : 1);
            const int slot = (f.horizontal ? 2 : 0) + (f.outwardSign > 0.0 ? 1 : 0);
            for (int history = 0; history < 2; ++history)
            {
                const bool old = history != 0;
                const double ub = value(r, old), outward = f.outwardSign * ub;
                const double densityAtBoundary = boundaryDensity(f, onlyOwner, outward);
                // Explicit reservoir continuation is a transported-state model
                // only: the solved boundary ub remains the mass-flux velocity.
                const double donor = f.boundary.type == Type::PressureReservoir && outward < 0.0
                    ? value(inner, old) : ub;
                (old ? fluxPre : flux)[slot] += densityAtBoundary * ub * f.area * donor;
            }
            // External total normal stress is prescribed. Incoming viscosity
            // is validated reservoir data, but cannot be added a second time
            // to a boundary whose full traction is already -p_res*n.
            assembly_.rhs[compact] -= f.outwardSign * f.area * f.boundary.pressure;
        }
        const double n = flux[1] - flux[0] + flux[3] - flux[2];
        const double nPre = fluxPre[1] - fluxPre[0] + fluxPre[3] - fluxPre[2];
        assembly_.rhs[compact] += -1.5 * n + 0.5 * nPre;
    }

    void hoop(const Face& f, double fullVolume)
    {
        const int r = f.nativeRow, compact = f.velocityRow;
        if (f.negativeCell >= 0 && f.positiveCell >= 0)
        {
            const double viscosity = nativeInterior_ && nativeInterior_->layout->nativeInteriorRows()[f.nativeRow]
                ? nativeInterior_->viscosity[f.nativeRow]
                : 0.5 * (mu_.value[f.positiveCell] + mu_.value[f.negativeCell]);
            const double coefficient = viscosity * fullVolume / (f.center.x * f.center.x);
            require(f.center.x > 0.0 && std::isfinite(coefficient) && coefficient >= 0.0,
                "invalid full-dual RZ hoop stress");
            add(r, r, coefficient);
            assembly_.rhs[compact] -= coefficient * value(r, false);
            return;
        }
        require(pressure(f.boundary.type), "unsupported free radial half-dual hoop stress");
        const int owner = f.negativeCell >= 0 ? f.negativeCell : f.positiveCell;
        const int inner = normalFace(f, owner, f.outwardSign > 0.0 ? 0 : 1);
        const long double radius = f.center.x;
        const long double delta = -f.outwardSign * 0.5L * f.normalSpacing;
        const long double dr = -f.outwardSign * f.normalSpacing;
        const long double logarithm = std::log1p(delta / radius);
        const long double innerWeight = -f.outwardSign * (delta - radius * logarithm) / dr;
        const long double boundaryWeight = -f.outwardSign * logarithm - innerWeight;
        const long double multiplier = 2.0L * Pi * edge(f).length() * mu_.value[owner];
        const double diagonal = static_cast<double>(multiplier * boundaryWeight);
        const double neighbor = static_cast<double>(multiplier * innerWeight);
        require(radius > 0.0L && radius + delta > 0.0L && std::isfinite(diagonal)
            && std::isfinite(neighbor) && boundaryWeight >= 0.0L && innerWeight >= 0.0L,
            "invalid integrated radial pressure-half-dual hoop stress");
        add(r, r, diagonal); add(r, inner, neighbor);
        assembly_.rhs[compact] -= diagonal * value(r, false) + neighbor * value(inner, false);
    }
};
}


HydroFittedMomentumSystem assembleHydroFittedMomentum(
    const QuadtreeMesh& mesh, const HydroFittedDomain& domain, double dt,
    const CellField<double>& cellRho, const CellField<double>& cellMu,
    const VertField<double>& vertexRho, const VertField<double>& vertexMu,
    const Eigen::VectorXd& currentNative, const Eigen::VectorXd& previousNative,
    const CellField<double>& bodyX, const CellField<double>& bodyY,
    const std::vector<double>& outletBackflowToleranceByVelocityRow,
    const HydroFittedNativeInterior* nativeInterior)
{
    return Assembly(mesh, domain, dt, cellRho, cellMu, vertexRho, vertexMu,
        currentNative, previousNative, bodyX, bodyY, outletBackflowToleranceByVelocityRow, nativeInterior).run();
}
