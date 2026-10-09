#include "Modules/Electricity/DielectricElectrostatics.h"
// S08_FROZEN_BASELINE_BEGIN
// Frozen experimental baseline: sharpHermiteP3, including its fitted entry.
// Retained only for explicit A/B comparisons and legacy case regression.
// Do not grow this P3 correction/stencil into sharpHHO2 or add SI-HHO-Q features.
// Electricity V2 is a separate HHO implementation; no implicit legacy fallback.
// Scope and source identities: docs/task/S08_Frozen_Baselines.json.
// S08_FROZEN_BASELINE_END
#include "Modules/Electricity/DielectricAssemblyDiagnostics.h"
#include "Modules/Electricity/SurfaceChargeSource.h"
#include "Modules/Electricity/FittedSharpDielectricElectrostatics.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Core/Geometry/InterfaceGeometry.h"
#include "Core/Geometry/CellChartPoint.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Solvers/Operators/QuadtreeFaceGradient.h"
#include "Modules/Electricity/Detail/DielectricCubicTrace.h"
#include "Modules/Electricity/Detail/DielectricAxisTransmission.h"
#include <variant>
#include "Modules/Electricity/Detail/DielectricInterfaceFlux.h"
#include "Modules/Electricity/Detail/DielectricTraceSampling.h"

// This TU is compiled with precise floating-point semantics and without LTO.
// No reference shape, test fixture, or Detail type crosses its public boundary.
namespace
{
namespace P3 = Quantum::Electricity::Detail::DielectricCubicTrace;
namespace Flux = Quantum::Electricity::Detail::DielectricInterfaceFlux;
namespace Sampling = Quantum::Electricity::Detail::DielectricTraceSampling;
namespace AxisTransmission = Quantum::Electricity::Detail::DielectricAxisTransmission;
class NativeTrace
{
    std::variant<P3::Reconstruction,AxisTransmission::Reconstruction> value_;
public:
    NativeTrace(P3::Reconstruction value):value_(std::move(value)){}
    NativeTrace(AxisTransmission::Reconstruction value):value_(std::move(value)){}
    NativeTrace(Vector2D p,Vector2D n,double k,double ks,double em,double ep,
        const std::vector<P3::SamplePoint>& samples,double h,bool rz):value_(std::in_place_type<P3::Reconstruction>,p,n,k,ks,em,ep,samples,h,rz){}
    const AxisTransmission::Reconstruction* axis() const{return std::get_if<AxisTransmission::Reconstruction>(&value_);}
    double conditionEstimate() const{return std::visit([](const auto& v){return v.conditionEstimate();},value_);}
    double normalizedReproductionError() const{return std::visit([](const auto& v){return v.normalizedReproductionError();},value_);}
    template<class... A>auto evaluate(A&&... a) const{return std::visit([&](const auto& v){return v.evaluate(std::forward<A>(a)...);},value_);}
    template<class... A>auto particularPotential(A&&... a) const{return std::visit([&](const auto& v){return v.particularPotential(std::forward<A>(a)...);},value_);}
    template<class... A>auto particularDisplacementNormal(A&&... a) const{return std::visit([&](const auto& v){return v.particularDisplacementNormal(std::forward<A>(a)...);},value_);}
    template<class... A>auto displacementNormalWeights(A&&... a) const{return std::visit([&](const auto& v){return v.displacementNormalWeights(std::forward<A>(a)...);},value_);}
    template<class... A>auto correctionWeights(A&&... a) const{return std::visit([&](const auto& v){return v.correctionWeights(std::forward<A>(a)...);},value_);}
    template<class... A>auto affineCorrection(A&&... a) const{return std::visit([&](const auto& v){return v.affineCorrection(std::forward<A>(a)...);},value_);}
};
using Solver = PolarizationPoissonFVSolverCellCentered;
using Correction = Solver::SegmentFluxCorrection;
using FittedBoundary = Solver::FittedBoundary;
using BCType = Solver::FittedBoundaryType;

// Most checks run per cell/face. Materialize their diagnostic only on failure;
// the string overload below retains the existing dynamically composed errors.
void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(std::string("[SharpDielectricElectrostatics] ") + message);
}

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error("[SharpDielectricElectrostatics] " + message);
}

double checked(long double value, const char* message)
{
    const double result = static_cast<double>(value);
    require(std::isfinite(value) && std::isfinite(result), message);
    return result;
}

double retainedCurvatureDerivative(const LevelSetHermiteGeometry::Snapshot& snapshot,
    const LevelSetHermiteGeometry::Frame& frame)
{
    require(frame.hasCellLocalPoint != frame.hasCellChartPoint,
        "retained electric frame requires exactly one local or boundary-chart provenance");
    if (snapshot.usesNativeLeaves())
    {
        require(frame.hasCellChartPoint && !frame.hasCellLocalPoint && frame.cellI == -1 && frame.cellJ == -1,
            "native electric curvature derivative requires actual leaf/chart provenance");
        return snapshot.curvatureTangentialDerivativeAtNativeCellChart(frame.nativeCellIndex, frame.cellChartPoint);
    }
    return frame.hasCellChartPoint
        ? snapshot.curvatureTangentialDerivativeAtCellChart(frame.cellI, frame.cellJ, frame.cellChartPoint)
        : snapshot.curvatureTangentialDerivativeAtCellLocal(frame.cellI, frame.cellJ, frame.cellLocalPoint);
}

bool lowerElectrode(Boundary boundary, bool xAxis)
{
    return boundary == (xAxis ? Boundary::west : Boundary::south);
}

bool upperElectrode(Boundary boundary, bool xAxis)
{
    return boundary == (xAxis ? Boundary::east : Boundary::north);
}

P3::Phase phase(const LevelSetHermiteGeometry::Snapshot& snapshot, Vector2D point)
{
    return snapshot.negative(point) ? P3::Phase::Negative : P3::Phase::Positive;
}

P3::Jump knownJump(const SurfaceChargeSource& source, const LevelSetHermiteGeometry::Frame& frame)
{
    const auto jet = source.evaluate(frame);
    require(std::isfinite(jet.density) && std::isfinite(jet.firstDerivative)
        && std::isfinite(jet.secondDerivative), "prescribed surface-charge jet is not finite");
    return {jet.density, jet.firstDerivative, jet.secondDerivative};
}

std::pair<Vector2D, Vector2D> faceEndpoints(const QuadtreeMesh& mesh, const FaceSegment& face)
{
    const auto& edge = face.isHorizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
    // Use the actual shared vertices, including exact domain endpoints. Do not
    // manufacture/clamp an endpoint by subtracting a rounded half edge length.
    Vector2D a = mesh.verts(edge.edgeToVert[0]).coord();
    Vector2D b = mesh.verts(edge.edgeToVert[1]).coord();
    const Vector2D tangent(-face.n.y, face.n.x);
    if ((b.x - a.x) * tangent.x + (b.y - a.y) * tangent.y < 0.0) std::swap(a, b);
    return {a, b};
}

void requireInClosedCell(const QuadtreeMesh& mesh, const Cell& cell, Vector2D point)
{
    Vector2D minimum(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
    Vector2D maximum(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
    for (int vertex : cell.cellToVert)
    {
        const auto p = mesh.verts(vertex).coord();
        minimum.x = std::min(minimum.x, p.x); minimum.y = std::min(minimum.y, p.y);
        maximum.x = std::max(maximum.x, p.x); maximum.y = std::max(maximum.y, p.y);
    }
    require(std::isfinite(point.x) && std::isfinite(point.y)
        && point.x >= minimum.x && point.x <= maximum.x
        && point.y >= minimum.y && point.y <= maximum.y,
        "boundary face or half-cell link is outside its certified closed owner cell");
}

bool certifyBoundary(const QuadtreeMesh& mesh, const LevelSetHermiteGeometry::Snapshot& snapshot,
    const FaceSegment& face, const std::pair<Vector2D, Vector2D>& endpoints)
{
    const auto& owner = mesh.cells(face.K);
    const auto center = owner.centerCoord();
    // A nonzero rigorous Bernstein closed-cell sign excludes every zero,
    // including tangency, a zero endpoint and an identically-zero wall edge.
    // Zero means UNPROVED, not a point value; reject instead of guessing phase.
    const int sign = snapshot.closedCellSign(center);
    require(sign == -1 || sign == 1,
        "boundary owner cell is not certified strictly single-phase; interface/wall contact or unresolved separation is unsupported");
    for (const auto& point : {center, face.faceCenter, endpoints.first, endpoints.second})
        requireInClosedCell(mesh, owner, point);
    const bool negative = sign < 0;
    const auto facePieces = snapshot.partitionSegment(endpoints.first, endpoints.second);
    const auto halfLinkPieces = snapshot.partitionSegment(center, face.faceCenter);
    require(!facePieces.empty() && !halfLinkPieces.empty(), "empty certified boundary partition");
    for (const auto* pieces : {&facePieces, &halfLinkPieces})
        for (const auto& piece : *pieces)
        {
            requireInClosedCell(mesh, owner, piece.begin);
            requireInClosedCell(mesh, owner, piece.end);
            require(piece.negative == negative, "boundary partition disagrees with the strict closed-cell sign certificate");
        }
    return negative;
}

struct PhaseElectricFlux
{
    bool mixed = false, negative = false, zeroAreaAxis = false, zeroNormalBoundary = false;
    double negativeMeasure = 0.0, positiveMeasure = 0.0; // ds in planar, 2*pi*r*ds in RZ.
    std::vector<int> indices;
    std::array<std::vector<double>, 2> weights;
    std::array<double, 2> constant{{0.0, 0.0}};
};

// Integrate each phase's E.n functional using the SAME global-face P3 fit.
// Its gradient is quadratic, hence two-point Gauss is exact for that local
// polynomial in real arithmetic. No extra QR or independently fitted trace.
PhaseElectricFlux phaseElectricFlux(const NativeTrace& trace,
    const Sampling::Window& window, const std::vector<P3::SegmentPiece>& pieces,
    Vector2D normal, const DielectricConfig& config, const P3::Jump* jump)
{
    PhaseElectricFlux result;
    result.mixed = true;
    result.indices = window.indices;
    const std::size_t count = window.indices.size();
    std::array<std::vector<long double>, 2> weights{{std::vector<long double>(count), std::vector<long double>(count)}};
    std::array<long double, 2> particular{{0.0L, 0.0L}}, lengths{{0.0L, 0.0L}};
    const long double offsetScale = 0.5L / std::sqrt(3.0L);
    for (const auto& piece : pieces)
    {
        const int side = piece.phase == P3::Phase::Negative ? 0 : 1;
        const double epsilon = side == 0 ? config.permittivityNegative : config.permittivityPositive;
        const double wholeLength = std::hypot(piece.end.x - piece.begin.x, piece.end.y - piece.begin.y);
        const double length = piece.retainedRange
            ? checked(wholeLength * piece.fraction,"retained phase electric-flux length overflow") : wholeLength;
        require(std::isfinite(length) && length > 0.0, "invalid phase electric-flux piece length");
        const auto nativePoint = [&](long double fraction) {
            const auto component = [&](double a, double b) {
                const long double delta = static_cast<long double>(b)-a;
                return checked(piece.fromEnd
                    ? static_cast<long double>(b)-delta*piece.offset-delta*piece.fraction*fraction
                    : static_cast<long double>(a)+delta*piece.offset+delta*piece.fraction*fraction,
                    "retained phase Gauss point overflow");
            };
            return Vector2D(component(piece.begin.x,piece.end.x),component(piece.begin.y,piece.end.y));
        };
        lengths[side] += CartesianMeasure::face(length,
            piece.retainedRange ? nativePoint(.5L).x
                : checked(.5L * piece.begin.x + .5L * piece.end.x, "phase face midpoint radius overflow"),
            config.axisymmetric);
        if(const auto* axis=trace.axis()) {
            const auto functional=axis->integratedDisplacementWeights(piece,normal);
            require(functional.size()==count,"harmonic phase electric-flux functional size mismatch");
            for(std::size_t i=1;i<count;++i)weights[side][i]+=functional[i]/epsilon;
            if(jump)particular[side]+=axis->particularIntegratedDisplacement(piece,normal,*jump)/epsilon;
            continue;
        }
        const long double scale = 0.5L * length / epsilon;
        require(std::isfinite(scale) && scale > 0.0L, "phase electric-flux integration scale overflow or underflow");
        for (int sign : {-1, 1})
        {
            const Vector2D point = piece.retainedRange ? nativePoint(.5L+sign*offsetScale) : Vector2D(
                checked(0.5L * piece.begin.x + 0.5L * piece.end.x
                    + sign * offsetScale * (static_cast<long double>(piece.end.x) - piece.begin.x), "phase Gauss x overflow"),
                checked(0.5L * piece.begin.y + 0.5L * piece.end.y
                    + sign * offsetScale * (static_cast<long double>(piece.end.y) - piece.begin.y), "phase Gauss y overflow"));
            const auto functional = trace.displacementNormalWeights(point, piece.phase, normal);
            require(functional.size() == count, "phase electric-flux functional size mismatch");
            const long double physicalScale = config.axisymmetric ? scale * CartesianMeasure::radialWeight(point.x) : scale;
            require(std::isfinite(physicalScale) && physicalScale > 0.0L,
                "phase physical electric-flux integration scale overflow or underflow");
            for (std::size_t i = 1; i < count; ++i) weights[side][i] += physicalScale * functional[i];
            if (jump) particular[side] += physicalScale * trace.particularDisplacementNormal(point, piece.phase, normal, *jump);
        }
    }
    result.negativeMeasure = checked(lengths[0], "negative phase measure overflow");
    result.positiveMeasure = checked(lengths[1], "positive phase measure overflow");
    require(result.negativeMeasure > 0.0 && result.positiveMeasure > 0.0, "mixed electric flux requires both positive phase measures");
    std::vector<double> sampleParticular(count, 0.0);
    if (jump)
        for (std::size_t i = 0; i < count; ++i)
            sampleParticular[i] = trace.particularPotential(window.samples[i].point, window.samples[i].phase, *jump);
    for (int side = 0; side < 2; ++side)
    {
        result.weights[side].assign(count, 0.0);
        long double tail = 0.0L, constant = particular[side];
        for (std::size_t i = 1; i < count; ++i)
        {
            result.weights[side][i] = checked(weights[side][i], "phase electric-flux weight overflow");
            tail += result.weights[side][i];
            constant -= result.weights[side][i]
                * (static_cast<long double>(sampleParticular[i]) - sampleParticular[0]);
        }
        result.weights[side][0] = checked(-tail, "phase electric-flux constant-mode weight overflow");
        result.constant[side] = checked(constant, "phase electric-flux particular constant overflow");
    }
    return result;
}

struct Stencil
{
    std::vector<double> epsilon;
    std::vector<Correction> correction;
    std::vector<std::size_t> physicalOwner;
    std::vector<double> boundaryRhs;
    std::vector<PhaseElectricFlux> electricFlux; // Allocated only for a conduction request; first owners only.
};

struct FittedContext
{
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot;
    const LevelSetInterfacePanels& panels;
    const std::vector<FittedBoundary>& boundaries;
    CartesianGridScope scope;
    bool native = false;
    const std::vector<const LevelSetInterfacePanels*>* surfaces = nullptr;
    bool completeNativeCollection = false;

    LevelSetHermiteGeometry::Frame projectedFrame(Vector2D point) const {
        auto closest=panels.projectedFrame(point);
        if(!surfaces || surfaces->size()==1)return closest.frame;
        double second=std::numeric_limits<double>::infinity();
        for(std::size_t i=1;i<surfaces->size();++i) {
            auto candidate=(*surfaces)[i]->projectedFrame(point);
            if(candidate.distance<closest.distance){second=closest.distance;closest=std::move(candidate);}
            else second=std::min(second,candidate.distance);
        }
        const auto spacing=panels.nativeChain().minimumActiveSpacing();
        const double roundoff=64*std::numeric_limits<double>::epsilon()
            *std::max({spacing.x,spacing.y,closest.distance,second});
        require(second-closest.distance>roundoff,"electric interface projection is tied between components");
        return closest.frame;
    }
};

bool axisNeckWindow(const QuadtreeMesh& mesh,const LevelSetHermiteGeometry::Snapshot& snapshot,
    const std::vector<const LevelSetInterfacePanels*>& surfaces,const LevelSetHermiteGeometry::Frame& frame,
    const Sampling::Window& window,const DielectricConfig& config)
{
    if(!config.nativeAxisNeckHarmonics || !config.axisymmetric || !snapshot.usesNativeLeaves()
        || !snapshot.radialAxisSymmetry() || snapshot.minCoord().x!=0 || !frame.hasCellChartPoint
        || snapshot.nativeCellBounds(frame.nativeCellIndex)[0].x!=0 || frame.normal.x<=0)return false;
    if(std::any_of(window.samples.begin(),window.samples.end(),[](const auto& s){return s.phase==P3::Phase::Negative;}))return false;
    double lower=std::numeric_limits<double>::infinity(),upper=-lower,radial=0;
    for(int cell:window.indices) {
        const auto bounds=snapshot.nativeCellBounds(cell);lower=std::min(lower,bounds[0].y);
        upper=std::max(upper,bounds[1].y);radial=std::max(radial,bounds[1].x);
    }
    const auto source=surfaces.front()->nativeChain().source();
    require(source && source.get()==&snapshot,"axis-neck model lost its retained scalar");
    const LevelSetMappedSegment axis(source,{0,lower},{0,upper});
    if(!axis.roots().empty() || axis.ranges().empty()
        || std::any_of(axis.ranges().begin(),axis.ranges().end(),[](const auto& r){return !r.negative;}))return false;
    bool present=false;
    for(const auto* surface:surfaces) {
        const auto& chain=surface->nativeChain();
        for(std::size_t g=0;g<chain.graphs().size();++g) {
            const auto bounds=snapshot.nativeCellBounds(chain.graphs()[g].nativeCellIndex);
            if(bounds[1].y<=lower || bounds[0].y>=upper || bounds[0].x>=radial)continue;
            present=true;
            if(chain.frameAt(static_cast<int>(g),.5).normal.x<=0
                || !chain.zeroSetPartialResolved(static_cast<int>(g),0,0,1))return false;
        }
    }
    return present;
}

void requireCompleteNativeCollection(const QuadtreeMesh& mesh,const FittedContext& fitted)
{
    std::vector<Vector2D> contacts;
    const auto& symmetry=fitted.panels.nativeChain().symmetryPlanes();
    for(const auto* surface:*fitted.surfaces) {
        require(surface->nativeChain().symmetryPlanes()==symmetry,
            "fitted components disagree on the scalar reflection planes");
        for(const auto point:surface->boundaryContacts())
            if(std::none_of(contacts.begin(),contacts.end(),[&](auto p){return p.x==point.x&&p.y==point.y;}))
                contacts.push_back(point);
    }
    const auto discovered=LevelSetInterfacePanels::discoverNativeComponents(mesh,fitted.snapshot,
        fitted.panels.axisymmetric(),contacts,symmetry);
    require(discovered.curves.size()==fitted.surfaces->size(),"fitted collection omits a scalar interface component");
    std::vector<bool> used(discovered.curves.size(),false);
    const auto sameGraph=[](const auto& a,const auto& b) {
        return a.nativeCellIndex==b.nativeCellIndex && a.coordinateIsX==b.coordinateIsX
            && a.localLower.x==b.localLower.x && a.localLower.y==b.localLower.y
            && a.localUpper.x==b.localUpper.x && a.localUpper.y==b.localUpper.y
            && a.begin.side==b.begin.side && a.begin.rootOrdinal==b.begin.rootOrdinal
            && a.end.side==b.end.side && a.end.rootOrdinal==b.end.rootOrdinal
            && a.controlledIntegration==b.controlledIntegration;
    };
    for(const auto* surface:*fitted.surfaces) {
        const auto& actual=surface->nativeChain();bool found=false;
        for(std::size_t k=0;k<discovered.curves.size();++k) {
            const auto& expected=discovered.curves[k];
            if(used[k] || actual.closed()!=expected.closed()
                || actual.beginEndpointKind()!=expected.beginEndpointKind()
                || actual.endEndpointKind()!=expected.endEndpointKind()
                || actual.graphs().size()!=expected.graphs().size())continue;
            if(std::equal(actual.graphs().begin(),actual.graphs().end(),expected.graphs().begin(),sameGraph)) {
                used[k]=true;found=true;break;
            }
        }
        require(found,"fitted collection duplicates or substitutes an original scalar curve");
    }
}

// Only geometry-dependent local work is retained. Charge jets are deliberately
// absent: every electric solve must evaluate the current trial's callback.
struct ChargedFaceFit
{
    std::size_t segment;
    LevelSetHermiteGeometry::Frame frame;
    Sampling::Window window;
    NativeTrace trace;
    std::vector<P3::SegmentPiece> pieces;
    Vector2D owner, neighbor;
    P3::Phase ownerPhase, neighborPhase;
    bool splitCurrent;
    P3::LinkGeometry linkGeometry = P3::LinkGeometry::Orthogonal;
};

struct NativePartition
{
    double negativeFraction = 0;
    bool mixed = false, negative = false;
    std::vector<P3::SegmentPiece> pieces;
};

NativePartition nativePartition(const FittedContext &fitted, Vector2D begin, Vector2D end)
{
    const LevelSetMappedSegment mapped(fitted.snapshot, begin, end);
    NativePartition result;
    long double negative = 0, positive = 0;
    for (const auto &range : mapped.ranges())
    {
        require(std::isfinite(range.length) && range.length > 0, "invalid retained native phase width");
        (range.negative ? negative : positive) += range.length;
        result.pieces.push_back({begin, end, range.negative ? P3::Phase::Negative : P3::Phase::Positive, true,
                                 range.fromEnd, range.offset, range.length});
    }
    require(std::abs(negative + positive - 1) <= 256 * std::numeric_limits<double>::epsilon(),
            "native phases do not cover the complete query segment");
    result.mixed = negative > 0 && positive > 0;
    result.negative = positive == 0;
    result.negativeFraction = checked(negative / (negative + positive), "native phase fraction overflow");
    return result;
}

struct PreparedFaces
{
    Stencil stencil;
    DielectricReport assemblyReport;
    std::vector<ChargedFaceFit> fits;
};

bool isInterior(const FaceSegment& face, std::size_t segment, const FittedContext* fitted)
{
    return fitted ? fitted->boundaries[segment].type == BCType::Interior : face.boundary == Boundary::non;
}

double electricDistance(const QuadtreeMesh& mesh, const FaceSegment& face, bool interior)
{
    if (interior) return face.df;
    const auto point = mesh.cells(face.K).centerCoord();
    const double distance = (face.faceCenter.x-point.x)*face.n.x + (face.faceCenter.y-point.y)*face.n.y;
    require(std::isfinite(distance) && distance > 0.0, "invalid actual electrode half-cell distance");
    return distance;
}

bool certifyFittedBoundary(const QuadtreeMesh& mesh, const FaceSegment& face,
    const std::pair<Vector2D, Vector2D>& endpoints, const FittedContext& fitted)
{
    const auto& owner = mesh.cells(face.K);
    const auto center = owner.centerCoord();
    int sign = fitted.snapshot->closedCellSign(center);
    if (sign == 0)
    {
        // Only a declared exact contact vertex may be a zero of this otherwise
        // strictly pure cell. Never exempt a near vertex, zero edge or hidden
        // interior root. Panels already certify the declared open-chain roots.
        Vector2D low(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
        Vector2D high(-low.x, -low.y);
        for (int label : owner.cellToVert)
        {
            const auto p = mesh.verts(label).coord();
            low.x=std::min(low.x,p.x); low.y=std::min(low.y,p.y);
            high.x=std::max(high.x,p.x); high.y=std::max(high.y,p.y);
        }
        unsigned corners = 0;
        for(const auto* surface:*fitted.surfaces) for (const auto& pin : surface->boundaryContacts())
        {
            if (pin.x==low.x && pin.y==low.y) corners |= 1u;
            if (pin.x==high.x && pin.y==low.y) corners |= 2u;
            if (pin.x == high.x && pin.y == high.y)
                corners |= 4u;
            if (pin.x == low.x && pin.y == high.y)
                corners |= 8u;
        }
        if (corners)
            sign = fitted.snapshot->closedCellSignExceptVertices(center, corners);
    }
    if (sign == 0 && fitted.native)
    {
        // Native wall cells may contain interface away from the actual
        // electric face/half-link. Certify THESE full closed paths instead:
        // every positive range has the center's phase; every zero must be an
        // exact declared contact at an input endpoint. Coincident, interior,
        // unresolved and undeclared roots still reject. No root is snapped.
        const bool negative = fitted.snapshot->negative(center);
        for (const auto &line : {endpoints, std::make_pair(center, face.faceCenter)})
        {
            const LevelSetMappedSegment query(fitted.snapshot, line.first, line.second);
            require(!query.ranges().empty(), "empty native boundary phase certificate");
            for (const auto &range : query.ranges())
                require(range.negative == negative, "native boundary face/half-link crosses material phases");
            for (std::size_t root = 0; root < query.roots().size(); ++root)
            {
                bool declared = false;
                for (bool atEnd : {false, true})
                    if (query.rootAtInputEndpoint(root, atEnd))
                    {
                        const auto point = atEnd ? line.second : line.first;
                        for(const auto* surface:*fitted.surfaces) for (const auto &pin : surface->boundaryContacts())
                            declared |= point.x == pin.x && point.y == pin.y;
                    }
                require(declared, "native electric boundary contains an interior or undeclared zero");
            }
        }
        sign = negative ? -1 : 1;
    }
    require(sign == -1 || sign == 1,
            "fitted boundary owner is not certified single-phase apart from declared zero contact vertices");
    for (auto point : {center, face.faceCenter, endpoints.first, endpoints.second})
        requireInClosedCell(mesh, owner, point);
    return sign < 0;
}

// Read-only native interior-face preparation on ONE retained scalar. Charge
// callbacks, physical-mate reversal, matrix writes and report accumulation
// stay in assemble(), in their original segment order.
struct NativeInteriorGeometry
{
    NativePartition link, face;
    std::map<int, long double> gradient;
    bool coarseFine = false, supportNeedsTrace = false;
    LevelSetHermiteGeometry::Frame frame;
};

NativeInteriorGeometry prepareNativeInteriorGeometry(const QuadtreeMesh& mesh,
    const FaceSegment& face, const FittedContext& fitted, const DielectricConfig& config,
    bool hasCharge, bool hasConductivity)
{
    const auto& snapshot = *fitted.snapshot;
    const auto endpoints = faceEndpoints(mesh, face);
    const auto owner = mesh.cells(face.K).centerCoord(), neighbor = mesh.cells(face.L).centerCoord();
    NativeInteriorGeometry result;
    result.link = nativePartition(fitted, owner, neighbor);
    result.face = nativePartition(fitted, endpoints.first, endpoints.second);
    result.coarseFine = mesh.cells(face.K).level() != mesh.cells(face.L).level();
    if (result.coarseFine)
    {
        const auto& edge = face.isHorizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
        // Per-face owned support: no mutable cross-face map enters workers.
        std::map<std::pair<int, int>, supCoeff> childCenters;
        const auto centers = [&](const Cell& cell, int child) -> const supCoeff& {
            const auto key = std::make_pair(cell.index(), child);
            auto found = childCenters.find(key);
            if (found == childCenters.end())
                found = childCenters.emplace(key, QuadtreeFaceGradient::childCenter(mesh, cell, child)).first;
            return found->second;
        };
        QuadtreeFaceGradient::visitGradient(mesh, edge, face.isHorizontal, centers,
            [&](int cell, double coefficient) {
                require(cell >= 0 && cell < mesh.cellIndexNum() && snapshot.nativeCellActive(cell),
                        "native electric virtual-child support crosses inactive cells");
                result.gradient[cell] += coefficient;
                // Point-potential samples, not phase volume averages.
                result.supportNeedsTrace |=
                    snapshot.negative(mesh.cells(mesh.cellIndexToLabel[cell]).centerCoord()) != result.face.negative;
            });
    }
    const bool cut = (config.permittivityNegative != config.permittivityPositive || hasCharge) &&
                     (result.link.mixed || result.face.mixed || result.supportNeedsTrace);
    const bool splitCurrent = hasConductivity && result.face.mixed;
    result.frame.normal = Vector2D(1.0, 0.0); // Original inactive mixture-helper argument.
    if (cut || splitCurrent)
        result.frame = fitted.projectedFrame(face.faceCenter);
    return result;
}

struct NativeFacePreparation
{
    std::vector<std::unique_ptr<NativeInteriorGeometry>> faces;
    std::vector<std::exception_ptr> errors;

    const NativeInteriorGeometry& at(std::size_t segment) const
    {
        if (errors.at(segment)) std::rethrow_exception(errors[segment]);
        require(static_cast<bool>(faces.at(segment)), "missing native interior-face preparation");
        return *faces[segment];
    }
};

NativeFacePreparation prepareNativeFaces(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const FittedContext* fitted, const DielectricConfig& config, bool hasCharge, bool hasConductivity)
{
    NativeFacePreparation result;
    if (!fitted || !fitted->native) return result;
    const auto& segments = graph.segments();
    result.faces.resize(segments.size());
    result.errors.resize(segments.size());
    std::vector<std::size_t> work;
    std::map<std::pair<bool, int>, std::size_t> physical;
    for (std::size_t s = 0; s < segments.size(); ++s)
        if (physical.emplace(std::make_pair(segments[s].isHorizontal, segments[s].edgeLabel), s).second &&
            isInterior(segments[s], s, fitted)) work.push_back(s);
    const auto prepare = [&](int index) {
        const auto s = work[index];
        try
        {
            result.faces[s] = std::make_unique<NativeInteriorGeometry>(
                prepareNativeInteriorGeometry(mesh, segments[s], *fitted, config, hasCharge, hasConductivity));
        }
        catch (...) { result.errors[s] = std::current_exception(); }
    };
    // First real query observes the caller's FP environment. Failures are
    // consumed only at their original segment, after all preceding boundary,
    // mate, charge and assembly checks; no exception crosses an OpenMP region.
    if (!work.empty()) prepare(0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) if(work.size() > 32)
#endif
    for (int i = 1; i < static_cast<int>(work.size()); ++i) prepare(i);
    return result;
}

Stencil assemble(const QuadtreeMesh &mesh, const FaceSegmentGraph &graph, const DielectricConfig &config,
                 const LevelSetHermiteGeometry::Snapshot &snapshot, DielectricReport &report,
                 DielectricSharpAssemblyDiagnostics *diagnostics, const SurfaceChargeSource *surfaceCharge,
                 const DielectricConductivity *conductivity, const FittedContext *fitted = nullptr,
                 PreparedFaces *prepared = nullptr)
{
    const auto &segments = graph.segments();
    const bool xAxis = config.electrodeAxis == "x";
    const bool symmetryAxis = config.axisymmetric && mesh.minCoord.x == 0.0 && snapshot.radialAxisSymmetry();
    const double voltage = config.upperPotential - config.lowerPotential;
    const bool native = fitted && fitted->native;
    const auto spacing = fitted && !native ? snapshot.cellSize() : mesh.cellSize;
    const double legacyLengthScale = std::max(spacing.x, spacing.y);
    const auto linkGeometry =
        native ? P3::LinkGeometry::ProjectedNormalCenters : P3::LinkGeometry::Orthogonal;
    const Sampling::Lookup sampling = fitted
                                          ? (native ? Sampling::Lookup(mesh, fitted->snapshot,Sampling::NativeBoundaryWindow::OneSided)
                                                    : Sampling::Lookup(mesh, fitted->snapshot, fitted->scope))
                                          : Sampling::Lookup(mesh, Sampling::Policy::GlobalFaceSymmetric);
    Stencil result;
    result.epsilon.resize(segments.size());
    result.correction.resize(segments.size());
    result.physicalOwner.resize(segments.size());
    result.boundaryRhs.assign(mesh.cellIndexNum(), 0.0);
    if (conductivity)
        result.electricFlux.resize(segments.size());
    const auto nativeGeometry = prepareNativeFaces(mesh, graph, fitted, config,
                                                   surfaceCharge != nullptr, conductivity != nullptr);
    const NativePartition noNativePartition;
    std::map<std::pair<bool, int>, std::size_t> physical;
    for (std::size_t s = 0; s < segments.size(); ++s)
    {
        const auto &face = segments[s];
        const bool interior = isInterior(face, s, fitted);
        const double faceMeasure = CartesianMeasure::face(face.Af, face.faceCenter.x, config.axisymmetric);
        const auto key = std::make_pair(face.isHorizontal, face.edgeLabel);
        const auto inserted = physical.emplace(key, s);
        const std::size_t first = inserted.first->second;
        result.physicalOwner[s] = first;
        if (!inserted.second)
        {
            const auto &mate = segments[first];
            require(interior && isInterior(mate, first, fitted) && face.K == mate.L && face.L == mate.K &&
                        face.Af == mate.Af && face.df == mate.df && face.n.x == -mate.n.x &&
                        face.n.y == -mate.n.y,
                    "internal physical face does not have exactly reversed owners and metrics");
            result.epsilon[s] = result.epsilon[first];
            result.correction[s] = result.correction[first];
            result.correction[s].constant = -result.correction[first].constant;
            for (auto &coefficient : result.correction[s].coefficients)
                coefficient.value = -coefficient.value;
            continue;
        }

        const auto endpoints = faceEndpoints(mesh, face);
        const auto owner = mesh.cells(face.K).centerCoord();
        if (!interior)
        {
            require(fitted || face.L == -1, "domain boundary unexpectedly has a neighboring cell");
            if (symmetryAxis && face.boundary == Boundary::west)
            {
                // The actual r=0 surface has zero revolved area. Its regular
                // homogeneous-Neumann flux is zero even when the interface
                // meets the axis. No material phase or wall-contact certificate
                // is invented, and no interface trace is evaluated at r=0.
                require(!xAxis && face.faceCenter.x == 0.0 && endpoints.first.x == 0.0 &&
                            endpoints.second.x == 0.0 && face.n.x == -1.0 && face.n.y == 0.0 &&
                            faceMeasure == 0.0 && std::isfinite(face.Af) && face.Af > 0.0 &&
                            std::isfinite(face.df) && face.df > 0.0,
                        "symmetry-axis face requires exact zero area and homogeneous radial Neumann data");
                // A positive placeholder is required by the generic coefficient
                // port; multiplied by zero area it is not a phase assignment.
                result.epsilon[s] = config.permittivityPositive;
                if (conductivity)
                    result.electricFlux[s].zeroAreaAxis = true;
                ++report.zeroAreaAxisFaces;
                ++report.zeroNeumannFaces;
                continue;
            }
            if(fitted && fitted->native && symmetryAxis
                && fitted->boundaries[s].type==BCType::HomogeneousNeumann
                && face.isHorizontal && face.n.x==0 && std::abs(face.n.y)==1) {
                // An explicitly retained axial reflection plane has even
                // scalar geometry. Homogeneous electric Neumann data impose
                // E.n=0 in EACH phase, including a mixed boundary face. Its
                // D and J fluxes are therefore exactly zero without assigning
                // the entire face to either phase or adding a contact law.
                bool reflected=false;
                for(const auto* surface:*fitted->surfaces) for(const auto& plane:surface->nativeChain().symmetryPlanes())
                    reflected|=plane.normalAxis==1 && plane.coordinate==face.faceCenter.y
                        && endpoints.first.y==plane.coordinate && endpoints.second.y==plane.coordinate
                        && (plane.coordinate==snapshot.minCoord().y || plane.coordinate==snapshot.maxCoord().y);
                if(reflected) {
                    require(fitted->boundaries[s].value==0,"reflection electric boundary must have zero normal data");
                    result.epsilon[s]=config.permittivityPositive; // Unused by the exact zero-Neumann row.
                    if(conductivity)result.electricFlux[s].zeroNormalBoundary=true;
                    ++report.certifiedBoundaryFaces;++report.zeroNeumannFaces;
                    continue;
                }
            }
            const bool negative = fitted ? certifyFittedBoundary(mesh, face, endpoints, *fitted)
                                         : certifyBoundary(mesh, snapshot, face, endpoints);
            if (conductivity)
                result.electricFlux[s].negative = negative;
            result.epsilon[s] = negative ? config.permittivityNegative : config.permittivityPositive;
            ++report.certifiedBoundaryFaces;
            // Dirichlet values are supplied through the original BC object;
            // no manufactured boundary constant or wall trace is added here.
            // Both the constant and every coefficient on zero-N walls stay 0.
            if (fitted ? fitted->boundaries[s].type == BCType::HomogeneousNeumann
                       : !lowerElectrode(face.boundary, xAxis) && !upperElectrode(face.boundary, xAxis))
                ++report.zeroNeumannFaces;
        }
        else
        {
            const auto neighbor = mesh.cells(face.L).centerCoord();
            const auto* nativeData = native ? &nativeGeometry.at(s) : nullptr;
            const auto& nativeLink = nativeData ? nativeData->link : noNativePartition;
            const auto& nativeFace = nativeData ? nativeData->face : noNativePartition;
            const double linkFraction =
                native ? nativeLink.negativeFraction : snapshot.negativeFraction(owner, neighbor);
            const double faceFraction = native ? nativeFace.negativeFraction
                                               : snapshot.negativeFraction(endpoints.first, endpoints.second);
            const bool coarseFine = nativeData && nativeData->coarseFine;
            const bool supportNeedsTrace = nativeData && nativeData->supportNeedsTrace;
            const bool cut = (config.permittivityNegative != config.permittivityPositive || surfaceCharge) &&
                             (native ? nativeLink.mixed || nativeFace.mixed || supportNeedsTrace
                                     : ((linkFraction > 0.0 && linkFraction < 1.0) ||
                                        (faceFraction > 0.0 && faceFraction < 1.0)));
            const bool splitCurrent =
                conductivity && (native ? nativeFace.mixed : faceFraction > 0.0 && faceFraction < 1.0);
            if (conductivity && !splitCurrent)
            {
                require(faceFraction == 0.0 || faceFraction == 1.0,
                        "single-phase electric flux has a mixed face");
                result.electricFlux[s].negative = native ? nativeFace.negative : faceFraction == 1.0;
            }
            // The inactive normal is only an argument for the link-mixture
            // helper; no fabricated interface frame is used by a trace.
            LevelSetHermiteGeometry::Frame frame;
            frame.normal = Vector2D(1.0, 0.0);
            if (cut || splitCurrent)
                frame = native   ? nativeData->frame
                        : fitted ? snapshot.projectedCellLocalFrame(face.faceCenter)
                                 : snapshot.projectedFrame(face.faceCenter);
            const auto mixture = Flux::evaluateGeometry(face.n, face.Af, face.df, linkFraction, faceFraction,
                                                        config.permittivityNegative,
                                                        config.permittivityPositive, frame.normal);
            result.epsilon[s] = mixture.kappa;
            if (coarseFine && !cut)
            {
                // -epsilon*A*dV/dn minus the ACTUAL-center TPFA term. This
                // correction and its exact opposite enter both FV rows.
                const double sign = face.isHorizontal ? face.n.y : face.n.x;
                std::map<int, long double> correction;
                for (const auto &[cell, coefficient] : nativeData->gradient)
                    correction[cell] =
                        -static_cast<long double>(mixture.kappa) * faceMeasure * sign * coefficient;
                const long double tpfa = static_cast<long double>(mixture.kappa) * faceMeasure / face.df;
                correction[mesh.cells(face.K).index()] -= tpfa;
                correction[mesh.cells(face.L).index()] += tpfa;
                for (const auto &[cell, value] : correction)
                    if (value != 0)
                        result.correction[s].coefficients.push_back(
                            {cell, checked(value, "native electric gradient coefficient overflow")});
                ++report.correctedPhysicalFaces;
            }
            if (cut || splitCurrent)
            {
                const auto window = native   ? sampling.faceWindow(mesh, face)
                                    : fitted ? sampling.faceWindow(face.faceCenter, face.isHorizontal)
                                             : sampling.faceWindow(face.faceCenter, face.isHorizontal,
                                                                   [&](const Vector2D &point) {
                                                                       return snapshot.negative(point);
                                                                   });
                const double curvatureDerivative =
                    fitted ? retainedCurvatureDerivative(snapshot, frame)
                           : snapshot.reconstruction() == LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
                                 ? snapshot.curvatureTangentialDerivative(frame.point)
                                 : snapshot.reconstructedCurvatureTangentialDerivative(frame.point);
                const auto &fine =
                    mesh.cells(mesh.cells(face.K).level() >= mesh.cells(face.L).level() ? face.K : face.L);
                const double lengthScale =
                    native ? std::max(fine.size().x, fine.size().y) : legacyLengthScale;
                const NativeTrace trace = [&]() -> NativeTrace {
                    try
                    {
                        if(native && axisNeckWindow(mesh,snapshot,*fitted->surfaces,frame,window,config))
                            return AxisTransmission::Reconstruction(frame.point,frame.normal,frame.curvature,curvatureDerivative,
                                config.permittivityNegative,config.permittivityPositive,window.samples,lengthScale);
                        return P3::Reconstruction(frame.point, frame.normal, frame.curvature,
                                                  curvatureDerivative, config.permittivityNegative,
                                                  config.permittivityPositive, window.samples, lengthScale,
                                                  config.axisymmetric);
                    }
                    catch (const std::exception &error)
                    {
                        throw std::runtime_error("P3 physical face (" + std::to_string(face.faceCenter.x) +
                                                 "," + std::to_string(face.faceCenter.y) +
                                                 "), K=" + std::to_string(face.K) +
                                                 ", L=" + std::to_string(face.L) + ", native support jump=" +
                                                 std::to_string(supportNeedsTrace) + ": " + error.what());
                    }
                }();
                report.maxTraceCondition = std::max(report.maxTraceCondition, trace.conditionEstimate());
                report.axisNeckHarmonicFaces+=trace.axis()!=nullptr;
                report.maxTraceReproductionError =
                    std::max(report.maxTraceReproductionError, trace.normalizedReproductionError());
                std::vector<P3::SegmentPiece> pieces = nativeFace.pieces;
                if (!native)
                    for (const auto &piece : snapshot.partitionSegment(endpoints.first, endpoints.second))
                        pieces.push_back({piece.begin, piece.end,
                                          piece.negative ? P3::Phase::Negative : P3::Phase::Positive});
                const auto ownerPhase = phase(snapshot, owner), neighborPhase = phase(snapshot, neighbor);
                P3::Jump jump;
                if (surfaceCharge)
                    jump = knownJump(*surfaceCharge, frame);
                if (splitCurrent)
                    result.electricFlux[s] = phaseElectricFlux(trace, window, pieces, face.n, config,
                                                               surfaceCharge ? &jump : nullptr);
                // Equal epsilon and no prescribed charge may need a phase
                // current split but must retain the original TPFA operator.
                if (!cut)
                {
                    // This is an interior face, so it has no electrode RHS;
                    // retain the ordinary conductance validity gate as well.
                    const double conductance = result.epsilon[s] * faceMeasure / face.df;
                    require(std::isfinite(conductance) && conductance > 0.0,
                            "non-finite/nonpositive sharp face conductance");
                    continue;
                }
                std::vector<double> weights;
                if (surfaceCharge)
                {
                    auto affine =
                        trace.affineCorrection(owner, ownerPhase, neighbor, neighborPhase, face.n, face.Af,
                                               face.df, mixture.kappa, pieces, jump, linkGeometry);
                    weights = std::move(affine.weights);
                    result.correction[s].constant = affine.constant;
                }
                else
                    weights = trace.correctionWeights(owner, ownerPhase, neighbor, neighborPhase, face.n,
                                                      face.Af, face.df, mixture.kappa, pieces, linkGeometry);
                require(weights.size() == window.indices.size(), "global P3 correction/sample size mismatch");
                auto &coefficients = result.correction[s].coefficients;
                for (std::size_t k = 0; k < weights.size(); ++k)
                    if (weights[k] != 0.0)
                        coefficients.push_back({window.indices[k], weights[k]});
                // A sampling window is Cartesian-ordered, not necessarily
                // solver-index ordered. Sort once before exact mate reversal.
                std::sort(coefficients.begin(), coefficients.end(),
                          [](const auto &a, const auto &b) { return a.cellIndex < b.cellIndex; });
                if (diagnostics)
                {
                    // Observe already computed values only. No extra geometry
                    // query, QR, alternate weights or symmetry enforcement.
                    DielectricSharpFaceDiagnostic record;
                    record.graphSegmentIndex = s;
                    record.horizontal = face.isHorizontal;
                    record.edgeIndex = face.isHorizontal ? mesh.hEdges(face.edgeLabel).index()
                                                         : mesh.vEdges(face.edgeLabel).index();
                    record.ownerCellIndex = mesh.cells(face.K).index();
                    record.neighborCellIndex = mesh.cells(face.L).index();
                    record.ownerNegative = ownerPhase == P3::Phase::Negative;
                    record.neighborNegative = neighborPhase == P3::Phase::Negative;
                    record.faceCenter = face.faceCenter;
                    record.outwardNormal = face.n;
                    record.faceArea = face.Af;
                    record.centerDistance = face.df;
                    record.epsilonFace = mixture.kappa;
                    record.physicalFaceMeasure = faceMeasure;
                    record.linkNegativeFraction = linkFraction;
                    record.faceNegativeFraction = faceFraction;
                    record.interfacePoint = frame.point;
                    record.interfaceNormal = frame.normal;
                    record.curvature = frame.curvature;
                    record.curvatureTangentialDerivative = curvatureDerivative;
                    record.condition = trace.conditionEstimate();
                    record.reproduction = trace.normalizedReproductionError();
                    record.correctionConstant = result.correction[s].constant;
                    record.sampleCellIndices = window.indices;
                    for (const auto &sample : window.samples)
                        record.sampleNegative.push_back(sample.phase == P3::Phase::Negative ? 1 : 0);
                    for (const auto &piece : pieces)
                        record.facePieces.push_back(
                            {piece.begin, piece.end, piece.phase == P3::Phase::Negative});
                    for (const auto &coefficient : coefficients)
                        record.correction.push_back({coefficient.cellIndex, coefficient.value});
                    diagnostics->correctedFaces.push_back(std::move(record));
                }
                ++report.correctedPhysicalFaces;
                if (prepared)
                {
                    require(fitted && surfaceCharge, "prepared face fits require fitted prescribed charge");
                    prepared->fits.push_back({s, frame, window, trace, pieces, owner, neighbor, ownerPhase,
                                              neighborPhase, splitCurrent, linkGeometry});
                }
            }
        }
        const double distance = fitted ? electricDistance(mesh, face, interior) : face.df;
        const double conductance = result.epsilon[s] * faceMeasure / distance;
        const bool axis = config.axisymmetric && face.faceCenter.x == 0.0 && face.boundary == Boundary::west;
        require(std::isfinite(conductance) && (axis ? conductance == 0.0 : conductance > 0.0),
                "invalid sharp face conductance (zero allowed only on the symmetry axis)");
        if (fitted ? fitted->boundaries[s].type == BCType::Dirichlet : upperElectrode(face.boundary, xAxis))
        {
            auto &rhs = result.boundaryRhs[mesh.cells(face.K).index()];
            rhs += conductance * (fitted ? fitted->boundaries[s].value : voltage);
            require(std::isfinite(rhs), "shifted electrode RHS overflow");
        }
    }
    require(fitted || report.certifiedBoundaryFaces + report.zeroAreaAxisFaces ==
                          static_cast<std::size_t>(2 * mesh.Nx + 2 * mesh.Ny),
            report.zeroAreaAxisFaces
                ? "not every domain boundary face received a no-contact or zero-area symmetry-axis check"
                : "not every domain boundary face received a strict no-contact certificate");
    require(fitted || report.zeroNeumannFaces == static_cast<std::size_t>(2 * (xAxis ? mesh.Nx : mesh.Ny)),
            "zero-Neumann boundary face count mismatch");
    return result;
}

Stencil updatePreparedCharge(const FaceSegmentGraph &graph, const DielectricConfig &config,
                             const SurfaceChargeSource &charge, const PreparedFaces &prepared)
{
    Stencil result = prepared.stencil;
    for (const auto &fit : prepared.fits)
    {
        const auto s = fit.segment;
        const auto &face = graph.segments()[s];
        const auto jump = knownJump(charge, fit.frame);
        if (fit.splitCurrent)
            result.electricFlux[s] =
                phaseElectricFlux(fit.trace, fit.window, fit.pieces, face.n, config, &jump);
        // Keep the original arithmetic and validity gates. QR, phase queries,
        // projection and face partitions are not repeated; no basis-response
        // interpolation or approximate matrix-freezing replaces this formula.
        const auto affine = fit.trace.affineCorrection(fit.owner, fit.ownerPhase, fit.neighbor,
                                                       fit.neighborPhase, face.n, face.Af, face.df,
                                                       result.epsilon[s], fit.pieces, jump, fit.linkGeometry);
        result.correction[s].constant = affine.constant;
    }
    for (std::size_t s = 0; s < result.physicalOwner.size(); ++s)
        if (result.physicalOwner[s] != s)
            result.correction[s].constant = -result.correction[result.physicalOwner[s]].constant;
    return result;
}

// Exact values, not a hash or sparse-pattern identity. The retained Snapshot
// must remain immutable through all aliases (the existing Snapshot contract).
// Actual mesh identities/coordinates and fitted BCs also belong to the key;
// changing charge, conductivity magnitudes or acceptance tolerances does not
// change the geometric fit. All public input gates still execute before use.
std::vector<double> fittedAssemblyKey(const QuadtreeMesh &mesh, const std::optional<UniformLeafGrid> &layout,
                                      const FaceSegmentGraph &graph, const FittedContext &fitted,
                                      const DielectricConfig &config, bool requestsCurrent)
{
    std::vector<double> key;
    if (fitted.native)
    {
        key.reserve(50u * mesh.cellIndexNum() + 17u * graph.segments().size());
        const auto append = [&](auto value) { key.push_back(static_cast<double>(value)); };
        append(-1);
        append(config.axisymmetric);
        append(config.nativeAxisNeckHarmonics);
        append(requestsCurrent);
        append(config.permittivityNegative);
        append(config.permittivityPositive);
        append(config.electrodeAxis == "x");
        append(config.lowerPotential);
        append(config.upperPotential);
        append(mesh.cellIndexNum());
        for (int index = 0; index < mesh.cellIndexNum(); ++index)
        {
            const int label = mesh.cellIndexToLabel[index];
            const auto &cell = mesh.cells(label);
            append(label);
            append(cell.index());
            append(cell.level());
            append(cell.area());
            append(cell.centerCoord().x);
            append(cell.centerCoord().y);
            append(cell.size().x);
            append(cell.size().y);
            append(fitted.snapshot->nativeCellActive(index));
            for (int vertex : cell.cellToVert)
            {
                append(vertex);
                append(mesh.verts(vertex).coord().x);
                append(mesh.verts(vertex).coord().y);
            }
            for (const auto &side : cell.neiCell)
                for (int neighbor : side)
                    append(neighbor);
        }
        append(fitted.completeNativeCollection);
        append(fitted.surfaces->size());
        for(const auto* surface:*fitted.surfaces) {
        const auto& incidence=surface->nativeChain().cellToGraphs();
        append(incidence.size());
        for(const auto& pieces:incidence) {
            append(pieces.size());
            for(int panel:pieces)append(panel);
        }
        const auto& nativeChain=surface->nativeChain();
        append(nativeChain.graphs().size());
        for(const auto& piece:nativeChain.graphs()) {
            append(piece.nativeCellIndex);append(piece.coordinateIsX);append(piece.controlledIntegration);
            append(piece.localLower.x);append(piece.localLower.y);append(piece.localUpper.x);append(piece.localUpper.y);
            append(piece.begin.side);append(piece.begin.rootOrdinal);append(piece.end.side);append(piece.end.rootOrdinal);
        }
        append(nativeChain.closed());
        append(static_cast<int>(nativeChain.beginEndpointKind()));
        append(static_cast<int>(nativeChain.endEndpointKind()));
        append(nativeChain.symmetryPlanes().size());
        for(const auto& plane:nativeChain.symmetryPlanes()) {append(plane.normalAxis);append(plane.coordinate);}
        append(surface->boundaryContacts().size());
        for (auto pin : surface->boundaryContacts())
        {
            append(pin.x);
            append(pin.y);
        }
        }
        append(graph.segments().size());
        for (std::size_t s = 0; s < graph.segments().size(); ++s)
        {
            const auto &face = graph.segments()[s];
            const auto &edge = face.isHorizontal ? mesh.hEdges(face.edgeLabel) : mesh.vEdges(face.edgeLabel);
            append(face.K);
            append(face.L);
            append(face.boundary);
            append(face.side);
            append(face.isHorizontal);
            append(face.edgeLabel);
            append(face.Af);
            append(face.df);
            append(face.n.x); append(face.n.y); append(face.faceCenter.x); append(face.faceCenter.y);
            append(edge.index()); append(edge.masterCell); append(edge.slaveCell);
            append(fitted.boundaries[s].type); append(fitted.boundaries[s].value);
        }
        return key;
    }
    require(layout.has_value(),"uniform fitted cache key has no actual lattice");
    const auto& grid=*layout;
    // Exact entry count of the serialized data below. Avoid repeated growth
    // and copies of this full identity key during a charge iteration.
    key.reserve(16u + grid.cells().size() + grid.vertices().size()
        + grid.vertical().size() + grid.horizontal().size()
        + static_cast<std::size_t>(mesh.cellIndexNum()) + 8u * grid.cells().size()
        + 2u * grid.vertices().size() + fitted.panels.cellToPanel().size()
        + 2u * fitted.panels.boundaryContacts().size() + 14u * graph.segments().size());
    const auto append = [&](auto value) { key.push_back(static_cast<double>(value)); };
    append(fitted.scope); append(config.axisymmetric); append(requestsCurrent);
    append(config.permittivityNegative); append(config.permittivityPositive);
    append(config.electrodeAxis=="x"); append(config.lowerPotential); append(config.upperPotential);
    for (const auto* labels : {&grid.cells(), &grid.vertices(), &grid.vertical(), &grid.horizontal()})
    {
        append(labels->size());
        for (int label : *labels) append(label);
    }
    append(mesh.cellIndexNum());
    for (int k=0; k<mesh.cellIndexNum(); ++k) append(mesh.cellIndexToLabel[k]);
    for (int label : grid.cells())
    {
        const auto& cell=mesh.cells(label);
        append(cell.index()); append(cell.centerCoord().x); append(cell.centerCoord().y); append(cell.area());
        for (int vertex : cell.cellToVert) append(vertex);
    }
    for (int label : grid.vertices())
    {
        append(mesh.verts(label).coord().x); append(mesh.verts(label).coord().y);
    }
    append(fitted.panels.cellToPanel().size());
    for (int panel : fitted.panels.cellToPanel()) append(panel);
    append(fitted.panels.boundaryContacts().size());
    for (auto pin : fitted.panels.boundaryContacts()) { append(pin.x); append(pin.y); }
    append(graph.segments().size());
    for (std::size_t s=0; s<graph.segments().size(); ++s)
    {
        const auto& face=graph.segments()[s];
        append(face.K); append(face.L); append(face.boundary); append(face.side);
        append(face.isHorizontal); append(face.edgeLabel); append(face.Af); append(face.df);
        append(face.n.x); append(face.n.y); append(face.faceCenter.x); append(face.faceCenter.y);
        append(fitted.boundaries[s].type); append(fitted.boundaries[s].value);
    }
    return key;
}

struct Fields
{
    std::vector<double> epsilon, potential, ex, ey;
    std::vector<DielectricInterfaceTraceReport> traces;
};

std::vector<double> diagnosticPermittivity(const InterfaceGeometry& geometry,
    const DielectricConfig& config, int count, DielectricReport& report)
{
    const auto& fraction = geometry.positiveFraction();
    require(fraction.size() == static_cast<std::size_t>(count), "shared PL area-fraction size mismatch");
    std::vector<double> epsilon(count);
    report.minPermittivity = std::numeric_limits<double>::infinity();
    report.maxPermittivity = 0.0;
    for (int k = 0; k < count; ++k)
    {
        const double h = fraction[k];
        require(std::isfinite(h) && h >= 0.0 && h <= 1.0, "invalid shared PL positive area fraction");
        // This arithmetic PL-H cell mixture is for cell output and D/epsilon
        // field diagnostics only. The FV operator uses actual Hermite segment
        // fractions; no diffuse/tanh or cell-mixture permittivity enters it.
        epsilon[k] = h <= 0.5
            ? config.permittivityNegative + h * (config.permittivityPositive - config.permittivityNegative)
            : config.permittivityPositive + (1.0 - h) * (config.permittivityNegative - config.permittivityPositive);
        require(std::isfinite(epsilon[k]) && epsilon[k] > 0.0, "invalid diagnostic PL cell permittivity");
        report.minPermittivity = std::min(report.minPermittivity, epsilon[k]);
        report.maxPermittivity = std::max(report.maxPermittivity, epsilon[k]);
    }
    return epsilon;
}

PicDomainBC electrodes(const DielectricConfig& config)
{
    PicDomainBC bc;
    bc.left.type = bc.right.type = bc.bottom.type = bc.top.type = PicBCType::Neumann;
    const double voltage = config.upperPotential - config.lowerPotential;
    if (config.electrodeAxis == "x")
    {
        bc.left.type = bc.right.type = PicBCType::Dirichlet;
        bc.left.value = 0.0; bc.right.value = voltage;
    }
    else
    {
        bc.bottom.type = bc.top.type = PicBCType::Dirichlet;
        bc.bottom.value = 0.0; bc.top.value = voltage;
    }
    return bc;
}

void reconstructFluxFields(const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const Stencil& stencil, const DielectricConfig& config, const Eigen::VectorXd& values,
    Fields& fields, DielectricReport& report, const SurfaceChargeSource* surfaceCharge,
    std::vector<double>* evaluatedFlux, const FittedContext* fitted = nullptr)
{
    const int count = mesh.cellIndexNum();
    const auto& segments = graph.segments();
    const bool xAxis = config.electrodeAxis == "x";
    const double voltage = config.upperPotential - config.lowerPotential;
    std::vector<double> flux(segments.size(), 0.0), weightX(count, 0.0), weightY(count, 0.0);
    std::vector<long double> balance(count, 0.0L), displacementX(count, 0.0L), displacementY(count, 0.0L);
    long double boundarySum = 0.0L, lowerSum = 0.0L, upperSum = 0.0L, boundaryAbs = 0.0L;
    double rhsNorm = 0.0;
    if (surfaceCharge)
    {
        // Normalize the independent physical residual by the actual affine
        // algebraic RHS: electrode + integrated charge - face constants.
        std::vector<long double> rhs(stencil.boundaryRhs.begin(), stencil.boundaryRhs.end());
        for (int k = 0; k < count; ++k) rhs[k] += surfaceCharge->cellIntegratedCharge[k];
        for (std::size_t s = 0; s < segments.size(); ++s)
            rhs[mesh.cells(segments[s].K).index()] -= stencil.correction[s].constant;
        for (long double value : rhs) rhsNorm = std::hypot(rhsNorm, checked(value, "charged affine RHS overflow"));
    }
    else
        for (double value : stencil.boundaryRhs) rhsNorm = std::hypot(rhsNorm, value);
    require(std::isfinite(rhsNorm), "independent electrode RHS norm overflow");

    for (std::size_t s = 0; s < segments.size(); ++s)
    {
        const auto& face = segments[s];
        const int owner = mesh.cells(face.K).index();
        const double faceMeasure = CartesianMeasure::face(face.Af, face.faceCenter.x, config.axisymmetric);
        if (stencil.physicalOwner[s] != s)
            flux[s] = -flux[stencil.physicalOwner[s]];
        else if (fitted ? fitted->boundaries[s].type != BCType::HomogeneousNeumann
            : face.boundary == Boundary::non || lowerElectrode(face.boundary, xAxis)
                 || upperElectrode(face.boundary, xAxis))
        {
            // Independently evaluate the full integrated FV outward D flux,
            // including every nonlocal linear correction. Not TPFA-only E,
            // not an analytic trace, and not a solver success flag as residual.
            const bool interior = isInterior(face,s,fitted);
            const double distance = fitted ? electricDistance(mesh,face,interior) : face.df;
            const double conductance = stencil.epsilon[s] * faceMeasure / distance;
            long double q = static_cast<long double>(conductance) * values[owner];
            if (interior)
                q -= static_cast<long double>(conductance) * values[mesh.cells(face.L).index()];
            else if (fitted || upperElectrode(face.boundary, xAxis))
                q -= static_cast<long double>(conductance) * (fitted ? fitted->boundaries[s].value : voltage);
            q += stencil.correction[s].constant;
            for (const auto& coefficient : stencil.correction[s].coefficients)
                q += static_cast<long double>(coefficient.value) * values[coefficient.cellIndex];
            flux[s] = checked(q, "non-finite full corrected displacement flux");
        }
        // A certified homogeneous-Neumann segment has exactly zero flux and
        // correction; it does not receive the Dirichlet owner diagonal term.
        balance[owner] += flux[s];
        if (config.axisymmetric)
        {
            // Recover a cell-center diagnostic from opposite face densities,
            // not an r-weighted average that doubles a linear radial field in
            // the first cell. Regularity gives D_r=0 on the actual axis.
            const long double density = faceMeasure == 0.0 ? 0.0L
                : static_cast<long double>(flux[s]) / faceMeasure;
            displacementX[owner] += face.n.x * density;
            displacementY[owner] += face.n.y * density;
            weightX[owner] += face.n.x * face.n.x;
            weightY[owner] += face.n.y * face.n.y;
        }
        else
        {
            displacementX[owner] += static_cast<long double>(face.n.x) * flux[s];
            displacementY[owner] += static_cast<long double>(face.n.y) * flux[s];
            weightX[owner] += face.Af * face.n.x * face.n.x;
            weightY[owner] += face.Af * face.n.y * face.n.y;
        }
        if (!isInterior(face,s,fitted))
        {
            boundarySum += flux[s]; boundaryAbs += std::abs(static_cast<long double>(flux[s]));
            if (!fitted && lowerElectrode(face.boundary, xAxis)) lowerSum += flux[s];
            if (!fitted && upperElectrode(face.boundary, xAxis)) upperSum += flux[s];
        }
    }
    double residualNorm = 0.0;
    fields.ex.resize(count); fields.ey.resize(count);
    for (int k = 0; k < count; ++k)
    {
        if (fitted && graph.cellOffsets()[k] == graph.cellOffsets()[k+1]) continue;
        const double residual = checked(surfaceCharge ? balance[k] - surfaceCharge->cellIntegratedCharge[k] : balance[k],
            "non-finite independent integrated cell residual");
        report.cellIntegratedFluxResidual = std::max(report.cellIntegratedFluxResidual, std::abs(residual));
        residualNorm = std::hypot(residualNorm, residual);
        require(std::isfinite(weightX[k]) && weightX[k] > 0.0 && std::isfinite(weightY[k]) && weightY[k] > 0.0,
            "invalid cell displacement reconstruction weights");
        fields.ex[k] = checked((displacementX[k] / weightX[k]) / fields.epsilon[k], "non-finite reconstructed cell Ex");
        fields.ey[k] = checked((displacementY[k] / weightY[k]) / fields.epsilon[k], "non-finite reconstructed cell Ey");
    }
    report.relativeResidual = rhsNorm > 0.0 ? residualNorm / rhsNorm : residualNorm;
    const auto fluxChecks = config.sharpFluxChecks.value_or(DielectricFluxChecks());
    fluxChecks.validate();
    require(std::isfinite(report.relativeResidual) && report.relativeResidual <= fluxChecks.reconstructedPoisson,
        "independently evaluated corrected Poisson relative residual exceeds reconstructedPoisson tolerance");
    report.netBoundaryFlux = checked(boundarySum, "net boundary displacement flux overflow");
    report.lowerBoundaryFlux = checked(lowerSum, "lower electrode displacement flux overflow");
    report.upperBoundaryFlux = checked(upperSum, "upper electrode displacement flux overflow");
    const double absoluteBoundaryFlux = checked(boundaryAbs, "absolute boundary displacement flux overflow");
    if (surfaceCharge)
    {
        long double charge = 0.0L, chargeAbs = 0.0L;
        for (double value : surfaceCharge->cellIntegratedCharge) { charge += value; chargeAbs += std::abs(value); }
        report.prescribedSurfaceCharge = checked(charge, "total prescribed surface charge overflow");
        report.surfaceChargeClosureResidual = checked(boundarySum - charge, "charged Gauss-law closure overflow");
        const double scale = checked(boundaryAbs + chargeAbs, "charged Gauss-law closure scale overflow");
        report.boundaryFluxImbalance = scale > 0.0
            ? std::abs(report.surfaceChargeClosureResidual) / scale : std::abs(report.surfaceChargeClosureResidual);
    }
    else
        report.boundaryFluxImbalance = absoluteBoundaryFlux > 0.0
            ? std::abs(report.netBoundaryFlux) / absoluteBoundaryFlux : std::abs(report.netBoundaryFlux);
    const double fluxTolerance = std::max(fluxChecks.gaussClosure,
        128.0 * std::numeric_limits<double>::epsilon() * count);
    require(std::isfinite(report.boundaryFluxImbalance) && report.boundaryFluxImbalance <= fluxTolerance,
        surfaceCharge ? "integrated boundary displacement flux fails prescribed-charge closure"
                      : "integrated boundary displacement flux fails charge-free closure");
    // cellIntegratedFluxResidual is an absolute, code-unit diagnostic, not a
    // new dimensionless acceptance gate borrowed from a normalized fixture.
    if (evaluatedFlux) *evaluatedFlux = std::move(flux);
}

std::shared_ptr<const DielectricConductionReport> conductionCurrent(const QuadtreeMesh& mesh,
    const FaceSegmentGraph& graph, const Stencil& stencil, const DielectricConfig& config,
    const DielectricConductivity& conductivity, const Eigen::VectorXd& values,
    const std::vector<double>& displacementFlux, const FittedContext* fitted = nullptr,
    std::vector<double>* evaluatedCurrent = nullptr)
{
    const auto& segments = graph.segments();
    require(stencil.electricFlux.size() == segments.size() && displacementFlux.size() == segments.size(),
        "incomplete actual FV electric-flux decomposition");
    auto result = std::make_shared<DielectricConductionReport>();
    const int count = mesh.cellIndexNum();
    std::vector<long double> cell(count, 0.0L);
    std::vector<double> current(segments.size(), 0.0);
    long double boundary = 0.0L, lower = 0.0L, upper = 0.0L, absoluteFace = 0.0L;
    const bool xAxis = config.electrodeAxis == "x";
    for (std::size_t s = 0; s < segments.size(); ++s)
    {
        const auto& face = segments[s];
        const std::size_t first = stencil.physicalOwner[s];
        if (first != s) current[s] = -current[first];
        else
        {
            const auto& split = stencil.electricFlux[s];
            std::array<double, 2> electric{{0.0, 0.0}};
            if (split.zeroAreaAxis)
                require(displacementFlux[s] == 0.0, "zero-area symmetry-axis displacement flux is not exactly zero");
            else if(split.zeroNormalBoundary)
                require(displacementFlux[s]==0.0,"reflection-plane displacement flux is not exactly zero");
            else if (split.mixed)
            {
                require(!split.indices.empty(), "mixed phase electric-flux stencil has no samples");
                for (int side = 0; side < 2; ++side)
                {
                    require(split.weights[side].size() == split.indices.size(), "phase electric-flux map size mismatch");
                    long double value = split.constant[side];
                    for (std::size_t i = 1; i < split.indices.size(); ++i)
                        value += split.weights[side][i]
                            * (static_cast<long double>(values[split.indices[i]]) - values[split.indices[0]]);
                    electric[side] = checked(value, "actual phase electric-flux evaluation overflow");
                }
                const long double remainder = static_cast<long double>(displacementFlux[s])
                    - static_cast<long double>(config.permittivityNegative) * electric[0]
                    - static_cast<long double>(config.permittivityPositive) * electric[1];
                const long double denominator = static_cast<long double>(config.permittivityNegative) * split.negativeMeasure
                    + static_cast<long double>(config.permittivityPositive) * split.positiveMeasure;
                require(std::isfinite(denominator) && denominator > 0.0L, "phase flux remainder denominator overflow or underflow");
                electric[0] = checked(static_cast<long double>(electric[0]) + (split.negativeMeasure / denominator) * remainder,
                    "negative electric-flux remainder overflow");
                electric[1] = checked(static_cast<long double>(electric[1]) + (split.positiveMeasure / denominator) * remainder,
                    "positive electric-flux remainder overflow");
            }
            else
            {
                const int side = split.negative ? 0 : 1;
                electric[side] = checked(static_cast<long double>(displacementFlux[s])
                    / (split.negative ? config.permittivityNegative : config.permittivityPositive),
                    "single-phase electric-flux evaluation overflow");
            }
            const long double dm = static_cast<long double>(config.permittivityNegative) * electric[0];
            const long double dp = static_cast<long double>(config.permittivityPositive) * electric[1];
            const double residual = checked(dm + dp - displacementFlux[s], "phase displacement recombination overflow");
            result->maxDisplacementSplitResidual = std::max(result->maxDisplacementSplitResidual, std::abs(residual));
            const long double scale = std::abs(dm) + std::abs(dp) + std::abs(static_cast<long double>(displacementFlux[s]));
            require(std::isfinite(scale) && std::abs(static_cast<long double>(residual))
                <= 4096.0L * std::numeric_limits<double>::epsilon() * scale,
                "phase electric fluxes do not reproduce the actual corrected displacement flux");
            const double jm = checked(static_cast<long double>(conductivity.negative) * electric[0], "negative phase current overflow");
            const double jp = checked(static_cast<long double>(conductivity.positive) * electric[1], "positive phase current overflow");
            require((conductivity.negative == 0.0 || electric[0] == 0.0 || jm != 0.0)
                && (conductivity.positive == 0.0 || electric[1] == 0.0 || jp != 0.0), "nonzero phase current underflow");
            current[s] = checked(static_cast<long double>(jm) + jp, "total face conduction current overflow");
            absoluteFace += std::abs(static_cast<long double>(current[s]));
        }
        cell[mesh.cells(face.K).index()] -= current[s];
        if (!isInterior(face,s,fitted))
        {
            boundary += current[s];
            if (!fitted && lowerElectrode(face.boundary, xAxis)) lower += current[s];
            if (!fitted && upperElectrode(face.boundary, xAxis)) upper += current[s];
        }
    }
    result->cellIntegratedCurrent.resize(count);
    long double total = 0.0L, absoluteCell = 0.0L;
    for (int k = 0; k < count; ++k)
    {
        result->cellIntegratedCurrent[k] = checked(cell[k], "integrated cell conduction current overflow");
        total += result->cellIntegratedCurrent[k];
        absoluteCell += std::abs(static_cast<long double>(result->cellIntegratedCurrent[k]));
    }
    result->lowerBoundaryCurrent = checked(lower, "lower electrode current overflow");
    result->upperBoundaryCurrent = checked(upper, "upper electrode current overflow");
    result->netBoundaryCurrent = checked(boundary, "net boundary current overflow");
    result->absoluteFaceCurrent = checked(absoluteFace, "absolute physical-face current overflow");
    result->totalCellCurrent = checked(total, "total cell conduction current overflow");
    result->closureResidual = checked(total + boundary, "charge/electrode current closure overflow");
    const double closureScale = checked(std::max(absoluteFace, absoluteCell), "current closure scale overflow");
    result->relativeClosure = closureScale > 0.0 ? std::abs(result->closureResidual) / closureScale : std::abs(result->closureResidual);
    require(std::isfinite(result->relativeClosure) && result->relativeClosure <= 4096.0 * std::numeric_limits<double>::epsilon(),
        "integrated cell and electrode currents fail conservation at roundoff");
    if (evaluatedCurrent) *evaluatedCurrent = std::move(current);
    return result;
}

Vector2D maxwellAction(Vector2D electric, double epsilon, Vector2D normal)
{
    const long double x = electric.x, y = electric.y;
    const long double squared = x * x + y * y;
    const double xx = checked(epsilon * (x * x - 0.5L * squared), "non-finite Cartesian Maxwell Txx");
    const double xy = checked(epsilon * x * y, "non-finite Cartesian Maxwell Txy");
    const double yy = checked(epsilon * (y * y - 0.5L * squared), "non-finite Cartesian Maxwell Tyy");
    return {checked(static_cast<long double>(xx) * normal.x + static_cast<long double>(xy) * normal.y,
                    "non-finite Cartesian Maxwell traction x"),
            checked(static_cast<long double>(xy) * normal.x + static_cast<long double>(yy) * normal.y,
                    "non-finite Cartesian Maxwell traction y")};
}

double normalTraction(const P3::State& state, Vector2D normal, const DielectricConfig& config)
{
    const Vector2D tangent(-normal.y, normal.x);
    const double en = state.normalDisplacement / config.permittivityNegative;
    const double ep = state.normalDisplacement / config.permittivityPositive;
    require(std::isfinite(en) && std::isfinite(ep), "non-finite phase-normal electric field");
    const Vector2D negative = tangent * state.tangentialElectricField + normal * en;
    const Vector2D positive = tangent * state.tangentialElectricField + normal * ep;
    const auto jump = maxwellAction(positive, config.permittivityPositive, normal)
        - maxwellAction(negative, config.permittivityNegative, normal);
    const double tau = checked(static_cast<long double>(jump.x) * normal.x
        + static_cast<long double>(jump.y) * normal.y, "non-finite Maxwell normal traction jump");
    const double tangentJump = checked(static_cast<long double>(jump.x) * tangent.x
        + static_cast<long double>(jump.y) * tangent.y, "non-finite Maxwell tangent traction jump");
    // Independent common-Dn/Et scalar identity checks the Cartesian tensor
    // action and sign. Tangential continuity here is a local transmission
    // constraint of the fit, not an independent global-accuracy claim.
    const long double et = state.tangentialElectricField, dn = state.normalDisplacement;
    const double common = checked(0.5L * (dn * ep - dn * en
        - et * et * (config.permittivityPositive - config.permittivityNegative)),
        "non-finite common-trace Maxwell identity");
    const long double magnitude = static_cast<long double>(config.permittivityNegative)
        * (static_cast<long double>(negative.x) * negative.x + static_cast<long double>(negative.y) * negative.y)
        + static_cast<long double>(config.permittivityPositive)
        * (static_cast<long double>(positive.x) * positive.x + static_cast<long double>(positive.y) * positive.y);
    const double scale = std::max(1e-30, checked(magnitude, "non-finite Maxwell identity scale"));
    const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon() * scale;
    require(std::abs(tau - common) <= tolerance && std::abs(tangentJump) <= tolerance,
        "Cartesian Maxwell jump violates the common-Dn/Et transmission identity");
    return tau;
}

struct ChargedTraction
{
    double normal = 0.0, tangent = 0.0, positiveDisplacement = 0.0;
};

ChargedTraction chargedTraction(const P3::State& state, const P3::Jump& charge,
    Vector2D normal, const DielectricConfig& config)
{
    const Vector2D tangent(-normal.y, normal.x);
    const double dn = state.normalDisplacement;
    const double dp = checked(static_cast<long double>(dn) + charge.value, "positive normal displacement overflow");
    const double en = dn / config.permittivityNegative, ep = dp / config.permittivityPositive;
    require(std::isfinite(en) && std::isfinite(ep), "non-finite charged phase-normal electric field");
    const Vector2D negative = tangent * state.tangentialElectricField + normal * en;
    const Vector2D positive = tangent * state.tangentialElectricField + normal * ep;
    const auto jump = maxwellAction(positive, config.permittivityPositive, normal)
        - maxwellAction(negative, config.permittivityNegative, normal);
    const double normalJump = checked(static_cast<long double>(jump.x) * normal.x
        + static_cast<long double>(jump.y) * normal.y, "non-finite charged Maxwell normal traction");
    const double tangentJump = checked(static_cast<long double>(jump.x) * tangent.x
        + static_cast<long double>(jump.y) * tangent.y, "non-finite charged Maxwell tangent traction");
    const long double et = state.tangentialElectricField;
    const double normalIdentity = checked(0.5L * (static_cast<long double>(dp) * ep
        - static_cast<long double>(dn) * en - et * et * (config.permittivityPositive - config.permittivityNegative)),
        "non-finite charged Maxwell normal identity");
    const double tangentIdentity = checked(static_cast<long double>(charge.value) * et,
        "non-finite surface-charge tangential traction");
    require(charge.value == 0.0 || et == 0.0L || tangentIdentity != 0.0,
        "nonzero surface-charge tangential traction underflow");
    const long double magnitude = static_cast<long double>(config.permittivityNegative)
        * (static_cast<long double>(negative.x) * negative.x + static_cast<long double>(negative.y) * negative.y)
        + static_cast<long double>(config.permittivityPositive)
        * (static_cast<long double>(positive.x) * positive.x + static_cast<long double>(positive.y) * positive.y);
    const double scale = std::max(1e-30, checked(magnitude, "non-finite charged Maxwell identity scale"));
    const double tolerance = 4096.0 * std::numeric_limits<double>::epsilon() * scale;
    require(std::abs(normalJump - normalIdentity) <= tolerance && std::abs(tangentJump - tangentIdentity) <= tolerance,
        "Cartesian Maxwell jump violates the prescribed-charge normal/tangential identities");
    return {normalJump, tangentIdentity, dp};
}

struct AxisCapFit
{
    Eigen::Matrix<double,5,1> coefficients;
    double h=0,offset=0,condition=0,reproduction=0,relativeResidual=0;
    std::vector<int> sampleCellIndices;
};
AxisCapFit prepareAxisCapFit(const DielectricConfig& config,
    const LevelSetHermiteGeometry::Snapshot& snapshot,const InterfaceTraceGeometry::Point& query,
    const Sampling::Window& window,const Eigen::VectorXd& values)
{
    const auto& frame=query.frame;
    const auto bounds=snapshot.nativeCellBounds(frame.nativeCellIndex);
    const double h=std::max(bounds[1].x-bounds[0].x,bounds[1].y-bounds[0].y);
    require(config.axisymmetric && snapshot.radialAxisSymmetry() && bounds[0].x==0
        && frame.point.x>=0 && frame.point.x<=bounds[1].x,"axis-cap traction escaped its declared radial cell");
    std::vector<int> indices;
    for(std::size_t i=0;i<window.indices.size();++i)
        if(window.samples[i].phase==P3::Phase::Positive)indices.push_back(static_cast<int>(i));
    require(indices.size()>=5,"axis-cap exterior potential requires at least five actual gas samples");
    Eigen::MatrixXd design(indices.size(),5);Eigen::VectorXd rhs(indices.size());
    const double offset=values[window.indices[indices.front()]];
    for(std::size_t i=0;i<indices.size();++i) {
        const auto& point=window.samples[indices[i]].point;
        const double r=point.x/h,z=(point.y-frame.point.y)/h;
        // Exact regular RZ harmonic basis: 1, z, z^2-r^2/2,
        // z^3-3*z*r^2/2, z^4-3*z^2*r^2+3*r^4/8. All solve the
        // cylindrical Laplace equation, including the regular r=0 limit.
        // The extra quartic term resolves the first omitted even radial
        // variation in a curved exterior field without widening the stencil.
        // The origin in z is the retained query height.
        design.row(i)<<1,z,z*z-.5*r*r,z*z*z-1.5*z*r*r,
            z*z*z*z-3*z*z*r*r+.375*r*r*r*r;
        rhs[i]=values[window.indices[indices[i]]]-offset;
    }
    require(design.allFinite() && rhs.allFinite(),"nonfinite axis-cap exterior fit data");
    const double tolerance=128*std::numeric_limits<double>::epsilon();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);qr.setThreshold(tolerance);
    require(qr.rank()==5,"rank-deficient axis-cap exterior harmonic fit");
    const Eigen::MatrixXd inverse=qr.solve(Eigen::MatrixXd::Identity(indices.size(),indices.size()));
    const double condition=design.cwiseAbs().rowwise().sum().maxCoeff()*inverse.cwiseAbs().rowwise().sum().maxCoeff();
    const double reproduction=(inverse*design-Eigen::Matrix<double,5,5>::Identity()).cwiseAbs().maxCoeff();
    require(inverse.allFinite() && std::isfinite(condition) && condition>=1 && condition<=1e8
        && std::isfinite(reproduction) && reproduction<=256*std::numeric_limits<double>::epsilon()*condition,
        "axis-cap exterior fit failed its conditioning or harmonic reproduction gate");
    const Eigen::Matrix<double,5,1> a=qr.solve(rhs);
    AxisCapFit result;result.coefficients=a;result.h=h;result.offset=offset;
    result.condition=condition;result.reproduction=reproduction;
    const double residual=(design*a-rhs).norm(),variation=rhs.norm();
    result.relativeResidual=variation>0?residual/variation:residual;
    require(std::isfinite(result.relativeResidual),"axis-cap model residual overflow");
    for(int i:indices)result.sampleCellIndices.push_back(window.indices[i]);
    return result;
}
DielectricInterfaceTraceReport reconstructAxisCapTrace(const DielectricConfig& config,
    const InterfaceTraceGeometry::Point& query,const AxisCapFit& fit,const SurfaceChargeSource& chargeSource)
{
    const auto& frame=query.frame;const auto& a=fit.coefficients;
    const double h=fit.h;
    const double radius=frame.point.x/h;
    // Analytic gradient, Hessian and third derivatives at local z=0.
    // Retain all quartic contributions in the transmission jet, not only E.
    const Vector2D e{(a[2]*radius-1.5*a[4]*radius*radius*radius)/h,
        (-a[1]+1.5*a[3]*radius*radius)/h};
    const auto applyJacobian=[&](Vector2D v) {
        return Vector2D{((a[2]-4.5*a[4]*radius*radius)*v.x+3*a[3]*radius*v.y)/(h*h),
            (3*a[3]*radius*v.x+(-2*a[2]+6*a[4]*radius*radius)*v.y)/(h*h)};
    };
    const auto dot=[](Vector2D a,Vector2D b){return a.x*b.x+a.y*b.y;};
    const auto n=frame.normal;const Vector2D t{-n.y,n.x};
    const auto jt=applyJacobian(t),jn=applyJacobian(n);
    const Vector2D second{(-9*a[4]*radius*t.x*t.x+6*a[3]*t.x*t.y+12*a[4]*radius*t.y*t.y)/(h*h*h),
        (3*a[3]*t.x*t.x+24*a[4]*radius*t.x*t.y-6*a[3]*t.y*t.y)/(h*h*h)};
    const double et=dot(e,t),en=dot(e,n),k=frame.curvature,ks=query.curvatureTangentialDerivative;
    const auto charge=knownJump(chargeSource,frame);
    P3::State state;
    state.interfacePotential=fit.offset+a[0]-.5*a[2]*radius*radius+.375*a[4]*radius*radius*radius*radius;
    state.tangentialElectricField=et;state.normalDisplacement=config.permittivityPositive*en-charge.value;
    state.tangentialFieldDerivative=dot(jt,t)-k*en;
    state.normalDisplacementDerivative=config.permittivityPositive*(dot(jt,n)+k*et)-charge.firstDerivative;
    state.tangentialFieldSecondDerivative=dot(second,t)-3*k*dot(jt,n)-ks*en-k*k*et;
    state.normalDisplacementSecondDerivative=config.permittivityPositive*(dot(second,n)-k*dot(jn,n)
        +2*k*dot(jt,t)+ks*et-k*k*en)-charge.secondDerivative;
    for(double v:{state.interfacePotential,state.tangentialElectricField,state.normalDisplacement,
        state.tangentialFieldDerivative,state.normalDisplacementDerivative,
        state.tangentialFieldSecondDerivative,state.normalDisplacementSecondDerivative})
        require(std::isfinite(v),"axis-cap transmission jet overflow");
    const auto traction=chargedTraction(state,charge,n,config);
    DielectricInterfaceTraceReport result;
    result.interfacePotential=state.interfacePotential;result.tangentialElectricField=et;
    result.normalDisplacement=state.normalDisplacement;result.positiveNormalDisplacement=traction.positiveDisplacement;
    result.tangentialFieldDerivative=state.tangentialFieldDerivative;
    result.normalDisplacementDerivative=state.normalDisplacementDerivative;
    result.tangentialFieldSecondDerivative=state.tangentialFieldSecondDerivative;
    result.normalDisplacementSecondDerivative=state.normalDisplacementSecondDerivative;
    result.normalTraction=traction.normal;result.tangentialTraction=traction.tangent;
    result.surfaceChargeDensity=charge.value;result.condition=fit.condition;result.reproduction=fit.reproduction;
    result.axisCapExteriorModel=true;
    result.modelRelativeResidual=fit.relativeResidual;result.sampleCellIndices=fit.sampleCellIndices;
    return result;
}

struct PreparedLoadTrace
{
    Sampling::Window window;
    std::unique_ptr<NativeTrace> fit;
    std::optional<AxisCapFit> axisCap;
    std::exception_ptr failure;
};

DielectricInterfaceTraceReport reconstructTrace(const DielectricConfig& config,
    const LevelSetHermiteGeometry::Snapshot& snapshot, const LevelSetHermiteGeometry::Frame& frame,
    double curvatureTangentialDerivative,
    const Sampling::Lookup& sampling, double lengthScale, double potentialOffset, bool cellLocal,
    const Eigen::VectorXd& values, DielectricReport& report,
    const SurfaceChargeSource* surfaceCharge,const PreparedLoadTrace* prepared=nullptr)
{
        // Exactly the same retained frame/kappa/kappa_s as Hydro capillarity.
        // Preserve the old global trace window unless the caller explicitly
        // selects FittedActive, whose source-local window cannot read solids.
        const auto localWindow = prepared ? Sampling::Window{} : cellLocal ? sampling.window(frame)
            : sampling.window(frame.point, [&](const Vector2D& point) { return snapshot.negative(point); });
        const auto& window=prepared?prepared->window:localWindow;
        double fitLength = lengthScale;
        if (snapshot.usesNativeLeaves())
        {
            const auto bounds = snapshot.nativeCellBounds(frame.nativeCellIndex);
            fitLength = std::max(bounds[1].x-bounds[0].x, bounds[1].y-bounds[0].y);
        }
        const auto localFit=prepared?nullptr:std::make_unique<NativeTrace>(frame.point, frame.normal, frame.curvature,
            curvatureTangentialDerivative, config.permittivityNegative,
            config.permittivityPositive, window.samples, fitLength, config.axisymmetric);
        const auto& fit=prepared?*prepared->fit:*localFit;
        P3::Jump charge;
        P3::State state;
        if (surfaceCharge)
        {
            charge = knownJump(*surfaceCharge, frame);
            state = fit.evaluate(Sampling::gather(window, values), charge);
        }
        else state = fit.evaluate(Sampling::gather(window, values));
        ++report.localTraceFits;
        report.maxTraceCondition = std::max(report.maxTraceCondition, fit.conditionEstimate());
        report.maxTraceReproductionError = std::max(report.maxTraceReproductionError, fit.normalizedReproductionError());
        double tau = 0.0;
        ChargedTraction fullTraction;
        if (surfaceCharge)
        {
            fullTraction = chargedTraction(state, charge, frame.normal, config);
            tau = fullTraction.normal;
        }
        else tau = normalTraction(state, frame.normal, config);

        DielectricInterfaceTraceReport trace;
        trace.interfacePotential = checked(static_cast<long double>(state.interfacePotential) + potentialOffset,
            "restored local interface potential overflow");
        trace.tangentialElectricField = state.tangentialElectricField;
        trace.normalDisplacement = state.normalDisplacement;
        trace.tangentialFieldDerivative = state.tangentialFieldDerivative;
        trace.normalDisplacementDerivative = state.normalDisplacementDerivative;
        trace.tangentialFieldSecondDerivative = state.tangentialFieldSecondDerivative;
        trace.normalDisplacementSecondDerivative = state.normalDisplacementSecondDerivative;
        trace.normalTraction = tau; trace.condition = fit.conditionEstimate();
        trace.reproduction = fit.normalizedReproductionError();
        trace.sampleCellIndices = window.indices;
        if(const auto* axis=fit.axis()) {
            trace.axisNeckHarmonicModel=true;
            trace.modelRelativeResidual=axis->relativeResidual(Sampling::gather(window,values),charge);
        }
        if (surfaceCharge)
        {
            trace.positiveNormalDisplacement = fullTraction.positiveDisplacement;
            trace.surfaceChargeDensity = charge.value;
            trace.tangentialTraction = fullTraction.tangent;
        }
        return trace;
}

void reconstructTraces(const QuadtreeMesh& mesh, const DielectricConfig& config,
    const InterfaceTraceGeometry& geometry, const Eigen::VectorXd& values, Fields& fields,
    DielectricReport& report, const SurfaceChargeSource* surfaceCharge)
{
    const auto& snapshot = *geometry.source();
    const bool symmetryAxis = config.axisymmetric && mesh.minCoord.x == 0.0 && snapshot.radialAxisSymmetry();
    const Sampling::Lookup sampling(mesh, symmetryAxis ? Sampling::Policy::RZAxisOneSided : Sampling::Policy::TraceStrict);
    const double lengthScale = std::max(mesh.cellSize.x, mesh.cellSize.y);
    fields.traces.reserve(geometry.points().size());
    for (const auto& point : geometry.points())
    {
        const auto& frame = point.frame;
        require(!(frame.hasCellLocalPoint && frame.hasCellChartPoint), "ambiguous trace chart provenance");
        const auto actual = frame.hasCellChartPoint
            ? snapshot.jetAtCellChart(frame.cellI, frame.cellJ, frame.cellChartPoint)
            : frame.hasCellLocalPoint ? snapshot.jetAtCellLocal(frame.cellI, frame.cellJ, frame.cellLocalPoint)
            : snapshot.jet(frame.point);
        const double curvatureDerivative = snapshot.reconstruction() == LevelSetHermiteGeometry::Reconstruction::EnrichedQuartic
            ? snapshot.curvatureTangentialDerivative(frame.point)
            : snapshot.reconstructedCurvatureTangentialDerivative(frame.point);
        require(actual.point.x == frame.point.x && actual.point.y == frame.point.y
            && actual.normal.x == frame.normal.x && actual.normal.y == frame.normal.y
            && actual.curvature == frame.curvature && actual.levelSetResidual == frame.levelSetResidual
            && std::isfinite(point.curvatureTangentialDerivative) && point.curvatureTangentialDerivative == curvatureDerivative
            && std::isfinite(frame.stationarityResidual) && frame.stationarityResidual >= 0.0
            && std::isfinite(point.endpointKktResidual) && point.endpointKktResidual >= 0.0,
            "trace query does not match the retained scalar jet");
        fields.traces.push_back(reconstructTrace(config, snapshot, frame, point.curvatureTangentialDerivative,
            sampling, lengthScale, config.lowerPotential, false, values, report, surfaceCharge));
    }
    require(report.localTraceFits == geometry.points().size(), "not every query received one actual-potential trace");
}

void finishDiagnostics(const QuadtreeMesh& mesh,
    const DielectricConfig& config, const Eigen::VectorXd& values, Fields& fields, DielectricReport& report)
{
    const int count = mesh.cellIndexNum();
    fields.potential.resize(count);
    report.minPotential = std::numeric_limits<double>::infinity();
    report.maxPotential = -std::numeric_limits<double>::infinity();
    long double energy = 0.0L;
    for (int k = 0; k < count; ++k)
    {
        fields.potential[k] = checked(static_cast<long double>(values[k]) + config.lowerPotential,
            "restored cell potential overflow");
        const double squared = checked(static_cast<long double>(fields.ex[k]) * fields.ex[k]
            + static_cast<long double>(fields.ey[k]) * fields.ey[k], "diagnostic electric field squared overflow");
        const auto& cell = mesh.cells(mesh.cellIndexToLabel[k]);
        const double area = CartesianMeasure::cell(cell.area(), cell.centerCoord().x, config.axisymmetric);
        energy += 0.5L * fields.epsilon[k] * squared * area;
        report.maxField = std::max(report.maxField, std::sqrt(squared));
        report.minPotential = std::min(report.minPotential, fields.potential[k]);
        report.maxPotential = std::max(report.maxPotential, fields.potential[k]);
    }
    report.fieldEnergy = checked(energy, "diagnostic cell electric energy overflow");
    require(std::isfinite(report.maxField)
        && std::isfinite(report.minPotential) && std::isfinite(report.maxPotential), "non-finite sharp electric diagnostics");
}
} // namespace

void DielectricElectrostatics::solveSharp(const QuadtreeMesh& mesh, const VertField<double>& levelSet,
    const DielectricConfig& config, const InterfaceGeometry& sharedGeometry, DielectricReport& report,
    DielectricSharpAssemblyDiagnostics* assemblyDiagnostics, const SurfaceChargeSource* surfaceCharge,
    const DielectricConductivity* conductivity, const std::shared_ptr<const InterfaceTraceGeometry>& traceQueries,
    std::shared_ptr<const DielectricFields>& publication)
{
    report.usedSharpInterface = true;
    if (config.sharpFluxChecks) config.sharpFluxChecks->validate();
    report.stage = "sharp-geometry";
    sharedGeometry.validateSource(mesh, levelSet);
    require(!traceQueries || traceQueries->source() == sharedGeometry.source(),
        "trace queries must retain the solved scalar source");
    const int count = mesh.cellIndexNum();
    if (conductivity)
        require(std::isfinite(conductivity->negative) && conductivity->negative >= 0.0
            && std::isfinite(conductivity->positive) && conductivity->positive >= 0.0,
            "phase conductivities must be finite and nonnegative");
    if (surfaceCharge)
    {
        report.stage = "sharp-charge-input";
        require(surfaceCharge->axisymmetric == config.axisymmetric,
            "surface-charge integral measure must match the electric coordinate system");
        require(surfaceCharge->cellIntegratedCharge.size() == static_cast<std::size_t>(count),
            "integrated surface-charge size must match solver cells");
        require(static_cast<bool>(surfaceCharge->evaluate), "prescribed surface-charge evaluator is missing");
        for (double value : surfaceCharge->cellIntegratedCharge)
            require(std::isfinite(value), "integrated surface charge is not finite");
        report.usedSurfaceCharge = true;
    }
    Fields next;
    report.stage = "sharp-material";
    next.epsilon = diagnosticPermittivity(sharedGeometry, config, count, report);
    report.stage = "sharp-assembly";
    DielectricSharpAssemblyDiagnostics nextAssembly;
    const auto stencil = assemble(mesh, graph_, config, sharedGeometry.snapshot(), report,
        assemblyDiagnostics ? &nextAssembly : nullptr, surfaceCharge, conductivity);
    require(!assemblyDiagnostics || nextAssembly.correctedFaces.size() == report.correctedPhysicalFaces,
        "incomplete actual global assembly diagnostics");

    CellField<double> rhs, shiftedPotential;
    rhs.initialize(mesh, 0.0);
    if (surfaceCharge)
        for (int k = 0; k < count; ++k)
        {
            const double charge = surfaceCharge->cellIntegratedCharge[k];
            const auto& cell = mesh.cells(mesh.cellIndexToLabel[k]);
            const double volume = CartesianMeasure::cell(cell.area(), cell.centerCoord().x, config.axisymmetric);
            rhs.value[k] = checked(static_cast<long double>(charge) / volume,
                "prescribed cell surface-charge density overflow");
            require(charge == 0.0 || rhs.value[k] != 0.0, "nonzero prescribed cell charge density underflow");
        }
    Solver::Options options;
    options.useSuperLU = true; options.requireSuperLU = true;
    options.useCG = false; options.useIC = false;
    options.reuseDirectPattern = false; options.reuseIterativePattern = false;
    options.freezeDirectMatrixUntilPatternChange = false;
    options.useTemporalExtrapolatedGuess = false;
    options.relTol = config.linearResidualTolerance;
    options.axisymmetric = config.axisymmetric;
    const auto bc = electrodes(config);
    Solver::SolveReport linear;
    report.stage = "sharp-solve";
    const bool solved = solver_.solve(mesh, graph_, bc, rhs, shiftedPotential, &stencil.epsilon, nullptr,
        &linear, options, nullptr, nullptr, nullptr, &stencil.correction);
    report.usedSuperLU = linear.usedSuperLU;
    // Retain the actual algebraic failure evidence if the solve is rejected;
    // a successful step replaces this with the independent flux residual below.
    report.relativeResidual = linear.finalResidual;
    require(solved && linear.success, "corrected Poisson solve failed: " + linear.message);
    require(linear.usedSuperLU && !linear.usedCG && !linear.directMatrixFrozen,
        "sharp electrostatics requires actual unfrozen SuperLU execution");
    require(std::isfinite(linear.finalResidual) && linear.finalResidual <= config.linearResidualTolerance,
        "actual corrected Poisson algebraic residual exceeds linearResidualTolerance");
    require(shiftedPotential.value && shiftedPotential.size() == count, "corrected Poisson returned invalid field storage");
    Eigen::VectorXd values(count);
    for (int k = 0; k < count; ++k)
    {
        require(std::isfinite(shiftedPotential.value[k]), "non-finite actual corrected Poisson potential");
        values[k] = shiftedPotential.value[k];
    }
    report.stage = "sharp-flux";
    std::vector<double> actualDisplacementFlux;
    reconstructFluxFields(mesh, graph_, stencil, config, values, next, report, surfaceCharge,
        conductivity ? &actualDisplacementFlux : nullptr);
    std::shared_ptr<const DielectricConductionReport> nextCurrent;
    if (conductivity)
    {
        report.stage = "sharp-conduction-current";
        nextCurrent = conductionCurrent(mesh, graph_, stencil, config, *conductivity, values, actualDisplacementFlux);
    }
    report.stage = "sharp-traces";
    if (traceQueries) reconstructTraces(mesh, config, *traceQueries, values, next, report, surfaceCharge);
    report.stage = "sharp-diagnostics";
    finishDiagnostics(mesh, config, values, next, report);

    auto published = std::make_shared<DielectricFields>(mesh);
    published->levelSet.assign(levelSet.value.get(), levelSet.value.get() + levelSet.size());
    published->axisymmetric = config.axisymmetric;
    published->potential = std::move(next.potential);
    published->permittivity = std::move(next.epsilon);
    published->ex = std::move(next.ex); published->ey = std::move(next.ey);
    published->interfaceTraces = std::move(next.traces);
    published->traceGeometry = traceQueries;
    // Complete every allocation before publishing fields or observations.
    report.stage = "complete";
    report.reason = surfaceCharge
        ? "prescribed-charge sharp Hermite/P3 electrostatics and full Maxwell traction evaluated"
        : "charge-free sharp Hermite/P3 electrostatics and actual-potential Maxwell traction evaluated";
    if (assemblyDiagnostics) assemblyDiagnostics->correctedFaces.swap(nextAssembly.correctedFaces);
    publication = std::move(published);
    report.conductionCurrent = std::move(nextCurrent);
    report.success = true;
}

std::vector<Vector2D> DielectricFields::interfaceTractions() const
{
    require(traceGeometry && interfaceTraces.size() == traceGeometry->points().size(), "sharp tractions have not been published");
    std::vector<Vector2D> result;
    result.reserve(interfaceTraces.size());
    for (const auto& trace : interfaceTraces) result.emplace_back(trace.normalTraction, trace.tangentialTraction);
    return result;
}

struct FittedSharpDielectricSolver::PreparedAssembly
{
    const QuadtreeMesh* mesh = nullptr;
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot;
    std::vector<double> key;
    PreparedFaces faces;
    // Same immutable scalar, actual cell centers, activity and permittivity
    // parameters as the full key. Not the geometric/FV face permittivity.
    // Filled after successful flux reconstruction; clear() owns invalidation.
    std::vector<double> diagnosticCellPermittivity;
};

FittedSharpDielectricSolver::FittedSharpDielectricSolver() = default;
FittedSharpDielectricSolver::~FittedSharpDielectricSolver() = default;

void FittedSharpDielectricSolver::clear()
{
    prepared_.reset();
    solver_.clearDirectCache();
}

FittedSharpDielectricResult solveFittedSharpDielectricElectrostatics(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const std::vector<Solver::FittedBoundary>& boundaries,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, CartesianGridScope scope,
    const DielectricConfig& config, const FittedPotentialInput& previousPotential,
    const SurfaceChargeSource& surfaceCharge, const DielectricConductivity* conductivity)
{
    // Compatibility one-shot entry uses the same assembler/solver and owns
    // no state beyond this call. A charge iteration instead retains an instance.
    FittedSharpDielectricSolver solver;
    return solver.solve(mesh,graph,boundaries,std::move(snapshot),panels,scope,
        config,previousPotential,surfaceCharge,conductivity);
}

FittedSharpDielectricResult FittedSharpDielectricSolver::solve(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const std::vector<Solver::FittedBoundary>& boundaries,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, CartesianGridScope scope,
    const DielectricConfig& config, const FittedPotentialInput& previousPotential,
    const SurfaceChargeSource& surfaceCharge, const DielectricConductivity* conductivity)
{
    return solveImpl(mesh,graph,boundaries,std::move(snapshot),{&panels},scope,config,
        previousPotential,surfaceCharge,conductivity);
}

FittedSharpDielectricResult FittedSharpDielectricSolver::solveNative(
    const QuadtreeMesh& mesh, const FaceSegmentGraph& graph,
    const std::vector<Solver::FittedBoundary>& boundaries,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, const DielectricConfig& config,
    const FittedPotentialInput& previousPotential, const SurfaceChargeSource& surfaceCharge,
    const DielectricConductivity* conductivity)
{
    return solveImpl(mesh,graph,boundaries,std::move(snapshot),{&panels},std::nullopt,config,
        previousPotential,surfaceCharge,conductivity);
}

FittedSharpDielectricResult FittedSharpDielectricSolver::solveNative(
    const QuadtreeMesh& mesh,const FaceSegmentGraph& graph,const std::vector<Solver::FittedBoundary>& boundaries,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const std::vector<const LevelSetInterfacePanels*>& surfaces,const DielectricConfig& config,
    const FittedPotentialInput& previousPotential,const SurfaceChargeSource& surfaceCharge,
    const DielectricConductivity* conductivity)
{
    return solveImpl(mesh,graph,boundaries,std::move(snapshot),surfaces,std::nullopt,config,
        previousPotential,surfaceCharge,conductivity,true);
}

FittedSharpDielectricResult FittedSharpDielectricSolver::solveImpl(
    const QuadtreeMesh &mesh, const FaceSegmentGraph &graph,
    const std::vector<Solver::FittedBoundary> &boundaries,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const std::vector<const LevelSetInterfacePanels*>& surfaces,
    std::optional<CartesianGridScope> scope, const DielectricConfig &config,
    const FittedPotentialInput &previousPotential, const SurfaceChargeSource &surfaceCharge,
    const DielectricConductivity *conductivity,bool completeNativeCollection)
{
    FittedSharpDielectricResult result;
    auto timedStage = std::chrono::steady_clock::now();
    const auto elapsedStage = [&]() {
        const auto now = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(now - timedStage).count();
        timedStage = now;
        return seconds;
    };
    auto &report = result.report;
    report.stage = "fitted-sharp-input";
    try
    {
        require(!surfaces.empty() && surfaces.front(),"fitted P3 requires a nonempty surface collection");
        require(!config.nativeAxisNeckHarmonics || !scope.has_value(),"axis-neck harmonics require native retained geometry");
        const auto& panels=*surfaces.front();
        require(snapshot && panels.usesSource(*snapshot),
                "fitted P3 requires the retained charge-panel scalar");
        require(!scope.has_value() || surfaces.size()==1,"multiple fitted interfaces require native geometry");
        for(std::size_t s=0;s<surfaces.size();++s) {
            const auto* surface=surfaces[s];
            require(surface && surface->usesSource(*snapshot)
                && surface->axisymmetric()==panels.axisymmetric()
                && surface->usesNativeLeaves()==snapshot->usesNativeLeaves(),
                "fitted interface collection mixes source, coordinates or leaf layouts");
            for(std::size_t k=0;k<s;++k)require(surface!=surfaces[k],"duplicate fitted interface component");
        }
        require(config.interfaceDiscretization == "sharpHermiteP3",
                "fitted charged solve requires explicit sharpHermiteP3");
        require(!config.boundaryFluxTolerance && !config.chargeClosureTolerance,
                "fitted electricity does not consume ordinary boundaryFlux/chargeClosure controls");
        if (config.sharpFluxChecks)
            config.sharpFluxChecks->validate();
        require(config.model == "leakyDielectric" || config.model == "perfectDielectric",
                "unsupported fitted dielectric model");
        require(std::isfinite(config.permittivityNegative) && config.permittivityNegative > 0.0 &&
                    std::isfinite(config.permittivityPositive) && config.permittivityPositive > 0.0,
                "phase permittivities must be finite and positive");
        require(std::isfinite(config.linearResidualTolerance) && config.linearResidualTolerance > 0.0 &&
                    config.linearResidualTolerance < 1.0,
                "invalid fitted linearResidualTolerance");
        require(
            config.axisymmetric == panels.axisymmetric() &&
                config.axisymmetric == surfaceCharge.axisymmetric &&
                (!config.axisymmetric || (snapshot->radialAxisSymmetry() && snapshot->minCoord().x == 0.0)),
            "fitted charge, panels and electric measures/symmetry must agree");
        const bool native = !scope.has_value();
        require(!config.nativeAxisNeckHarmonics || native,"axis-neck harmonics require native fitted geometry");
        require(native == snapshot->usesNativeLeaves(), "fitted entry and scalar leaf addressing disagree");
        std::optional<UniformLeafGrid> layout;
        const int count = mesh.cellIndexNum();
        if (native)
        {
            validateCartesianLeafMesh(mesh);
            // Own/check the actual source bounds and indices before consuming
            // activity or a cached matrix; the source pointer alone is not a
            // certificate for mutable mesh geometry.
            const Sampling::Lookup checkedSource(mesh, snapshot);
        }
        else
        {
        layout.emplace(mesh,*scope);
        const auto& grid=*layout;
        require(grid.nx()==snapshot->nx() && grid.ny()==snapshot->ny()
            && grid.spacing().x==snapshot->cellSize().x && grid.spacing().y==snapshot->cellSize().y
            && grid.minCoord().x==snapshot->minCoord().x && grid.minCoord().y==snapshot->minCoord().y
            && grid.maxCoord().x==snapshot->maxCoord().x && grid.maxCoord().y==snapshot->maxCoord().y,
            "fitted scalar does not match actual uniform leaf lattice");
        }
        require(previousPotential.size()==count
            && surfaceCharge.cellIntegratedCharge.size()==static_cast<std::size_t>(count)
            && static_cast<bool>(surfaceCharge.evaluate), "missing fitted potential or prescribed charge data");
        std::vector<unsigned char> active(count,0);
        if (native)
            for (int k=0;k<count;++k) active[k]=snapshot->nativeCellActive(k)?1:0;
        else
            for(int j=0;j<layout->ny();++j) for(int i=0;i<layout->nx();++i)
                active[mesh.cells(layout->cells()[i+layout->nx()*j]).index()]=snapshot->cellActive(i,j)?1:0;
        // Rebuild only the topology view; no solve, geometry replacement or
        // alternate material mask. This validates identities before indexing
        // caller graph entries, including solid-neighbor boundary segments.
        FaceSegmentGraph expected;
        expected.build(mesh,[&](const QuadtreeMesh&,const Cell& cell){return active[cell.index()]!=0;});
        require(graph.cellOffsets()==expected.cellOffsets()
            && graph.segments().size()==expected.segments().size()
            && boundaries.size()==graph.segments().size(), "fitted graph does not match the scalar activity");
        std::size_t boundaryCount=0;
        for(std::size_t s=0;s<boundaries.size();++s)
        {
            const auto& a=graph.segments()[s]; const auto& b=expected.segments()[s];
            require(a.K==b.K && a.L==b.L && a.boundary==b.boundary && a.side==b.side
                && a.isHorizontal==b.isHorizontal && a.edgeLabel==b.edgeLabel && a.Af==b.Af && a.df==b.df
                && a.n.x==b.n.x && a.n.y==b.n.y && a.faceCenter.x==b.faceCenter.x && a.faceCenter.y==b.faceCenter.y,
                "fitted graph identities or physical metrics changed");
            const bool interior=a.L>=0 && active[mesh.cells(a.L).index()];
            const auto& bc=boundaries[s];
            require(std::isfinite(bc.value) && (bc.type==BCType::Interior || bc.type==BCType::Dirichlet
                || bc.type==BCType::HomogeneousNeumann), "invalid fitted electric boundary");
            require(interior ? bc.type==BCType::Interior && bc.value==0.0
                : bc.type!=BCType::Interior && (bc.type==BCType::Dirichlet || bc.value==0.0),
                "fitted boundary classification disagrees with actual activity");
            if(!interior) ++boundaryCount;
            if(config.axisymmetric && a.faceCenter.x==0.0)
                require(!interior && a.boundary==Boundary::west && bc.type==BCType::HomogeneousNeumann && bc.value==0.0,
                    "actual zero-area axis requires homogeneous Neumann data");
        }
        for(int k=0;k<count;++k)
        {
            require(previousPotential.label(k)==mesh.cellIndexToLabel[k], "fitted initial potential has stale native cell labels");
            require(std::isfinite(surfaceCharge.cellIntegratedCharge[k]), "non-finite prescribed cell charge");
            require(active[k] || surfaceCharge.cellIntegratedCharge[k]==0.0, "prescribed charge lies in an inactive cell");
            require(std::isfinite(previousPotential[k]), "fitted initial potential, including retained inactive payload, must be finite");
        }
        if(conductivity) require(std::isfinite(conductivity->negative) && conductivity->negative>=0.0
            && std::isfinite(conductivity->positive) && conductivity->positive>=0.0,
            "phase conductivities must be finite and nonnegative");
        const FittedContext fitted{snapshot,panels,boundaries,scope.value_or(CartesianGridScope::UniformLeaves),native,
            &surfaces,completeNativeCollection};
        report.usedSharpInterface=true; report.usedSurfaceCharge=true;
        report.stage="fitted-sharp-assembly";
        auto key=fittedAssemblyKey(mesh,layout,graph,fitted,config,conductivity!=nullptr);
        const bool reuse=prepared_ && prepared_->mesh==&mesh && prepared_->snapshot==snapshot
            && prepared_->key.size()==key.size()
            && std::memcmp(prepared_->key.data(),key.data(),key.size()*sizeof(double))==0;
        result.timings.input_s = elapsedStage();
        Stencil stencil;
        if (reuse)
        {
            report=prepared_->faces.assemblyReport;
            stencil=updatePreparedCharge(graph,config,surfaceCharge,prepared_->faces);
            result.geometryAssemblyReused=true;
        }
        else
        {
            if(completeNativeCollection)requireCompleteNativeCollection(mesh,fitted);
            auto prepared=std::make_unique<PreparedAssembly>();
            prepared->mesh=&mesh; prepared->snapshot=snapshot; prepared->key=std::move(key);
            stencil=assemble(mesh,graph,config,*snapshot,report,nullptr,&surfaceCharge,conductivity,&fitted,
                &prepared->faces);
            prepared->faces.stencil=stencil;
            prepared->faces.assemblyReport=report;
            prepared_=std::move(prepared);
        }
        require(report.certifiedBoundaryFaces+report.zeroAreaAxisFaces==boundaryCount,
            "not every actual fitted boundary received a phase/contact or zero-axis check");
        result.timings.geometry_s = elapsedStage();
        CellField<double> rhs(mesh,0.0);
        auto potential=std::make_shared<CellField<double>>(*previousPotential.field_);
        for(int k=0;k<count;++k) if(active[k])
        {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[k]);
            const double volume=CartesianMeasure::cell(cell.area(),cell.centerCoord().x,config.axisymmetric);
            const double charge=surfaceCharge.cellIntegratedCharge[k];
            rhs.value[k]=checked(static_cast<long double>(charge)/volume,"fitted cell charge density overflow");
            require(charge==0.0 || rhs.value[k]!=0.0,"nonzero fitted cell charge density underflow");
        }
        Solver::Options options;
        options.useSuperLU=true; options.requireSuperLU=true; options.useCG=false; options.useIC=false;
        options.reuseDirectPattern=false; options.reuseIterativePattern=false;
        options.freezeDirectMatrixUntilPatternChange=false; options.useTemporalExtrapolatedGuess=false;
        options.relTol=config.linearResidualTolerance; options.axisymmetric=config.axisymmetric;
        options.reuseFittedFactorization=true;
        auto& linear=result.linearSolve;
        report.stage="fitted-sharp-solve";
        bool solved = false;
        if (reuse)
        {
            // The full live-geometry/material/boundary checks above still ran.
            // updatePreparedCharge copies the privately owned stencil and ONLY
            // changes its affine constants; no mutable external coefficient
            // vector can enter this owned-operator RHS path.
            std::vector<double> constants(stencil.correction.size());
            for (std::size_t s = 0; s < constants.size(); ++s)
                constants[s] = stencil.correction[s].constant;
            solved = solver_.solvePreparedFittedRhs(rhs, *potential, &constants,
                options.relTol, &linear);
        }
        else
            solved = solver_.solveFittedBoundary(mesh,graph,rhs,*potential,stencil.epsilon,
                boundaries,stencil.correction,options,&linear);
        report.usedSuperLU=linear.usedSuperLU; report.relativeResidual=linear.finalResidual;
        require(solved && linear.success && linear.usedSuperLU && !linear.usedCG && !linear.directMatrixFrozen,
            "fitted corrected Poisson solve failed: "+linear.message);
        require(std::isfinite(linear.finalResidual) && linear.finalResidual<=config.linearResidualTolerance,
            "fitted corrected algebraic residual exceeds linearResidualTolerance");
        result.timings.linear_s = elapsedStage();
        Eigen::VectorXd values=Eigen::VectorXd::Zero(count);
        Fields diagnostic;
        const bool reuseCellPhase = !prepared_->diagnosticCellPermittivity.empty();
        if (reuseCellPhase)
        {
            require(prepared_->diagnosticCellPermittivity.size() == static_cast<std::size_t>(count),
                "prepared diagnostic cell permittivity indexing mismatch");
            diagnostic.epsilon = prepared_->diagnosticCellPermittivity;
        }
        else
            diagnostic.epsilon.assign(count,std::numeric_limits<double>::quiet_NaN());
        report.minPotential=std::numeric_limits<double>::infinity();
        report.maxPotential=-report.minPotential;
        for(int k=0;k<count;++k) if(active[k])
        {
            values[k]=potential->value[k];
            require(std::isfinite(values[k]),"non-finite fitted solved potential");
            report.minPotential=std::min(report.minPotential,values[k]);
            report.maxPotential=std::max(report.maxPotential,values[k]);
            // Cell E is not a phase trace. This center-phase coefficient is
            // used ONLY by the existing diagnostic D/epsilon recovery below.
            if (!reuseCellPhase)
            {
                diagnostic.epsilon[k]=snapshot->negative(mesh.cells(mesh.cellIndexToLabel[k]).centerCoord())
                    ? config.permittivityNegative : config.permittivityPositive;
                ++result.diagnosticCellPhaseQueries;
            }
        }
        report.stage="fitted-sharp-flux";
        reconstructFluxFields(mesh,graph,stencil,config,values,diagnostic,report,&surfaceCharge,
            &result.displacementFlux,&fitted);
        if (!reuseCellPhase)
            prepared_->diagnosticCellPermittivity = std::move(diagnostic.epsilon);
        result.timings.flux_s = elapsedStage();
        if(conductivity)
        {
            report.stage="fitted-sharp-conduction-current";
            report.conductionCurrent=conductionCurrent(mesh,graph,stencil,config,*conductivity,values,
                result.displacementFlux,&fitted,&result.conductionCurrentBySegment);
        }
        result.timings.current_s = elapsedStage();
        result.epsilonBySegment=std::move(stencil.epsilon);
        result.corrections=std::move(stencil.correction);
        report.stage="complete";
        report.reason="fitted prescribed-charge P3 potential, corrected FV flux and requested phase currents; no momentum loads";
        result.potential=FittedElectricPotential(std::move(potential),snapshot);
        report.success=true;
        return result;
    }
    catch(const std::exception& error)
    {
        clear();
        throw std::runtime_error("[FittedSharpDielectricElectrostatics] stage="+report.stage+": "+error.what());
    }
    catch(...)
    {
        clear();
        throw;
    }
}

std::vector<Vector2D> FittedSharpDielectricTractions::tractionComponents() const
{
    std::vector<Vector2D> result;
    result.reserve(interfaceTraces.size());
    for (const auto& trace : interfaceTraces)
        result.push_back({trace.normalTraction, trace.tangentialTraction});
    return result;
}

static FittedSharpDielectricTractions reconstructFittedTractions(
    const QuadtreeMesh &mesh, std::shared_ptr<const InterfaceTraceGeometry> traceGeometry,
    const LevelSetInterfacePanels &panels, std::optional<CartesianGridScope> scope,
    const DielectricConfig &config, const FittedElectricPotential &absolutePotential,
    const SurfaceChargeSource &surfaceCharge,const FittedAxisCapTractionModel* axisCapModel=nullptr)
{
    const auto snapshot = traceGeometry ? traceGeometry->source() : nullptr;
    DielectricReport report;
    report.stage = "fitted-Maxwell-input";
    try
    {
        require(snapshot && panels.usesSource(*snapshot),
                "Maxwell loads require the retained charge-panel scalar");
        require(config.interfaceDiscretization == "sharpHermiteP3" &&
                    (config.model == "leakyDielectric" || config.model == "perfectDielectric"),
                "fitted Maxwell loads require an explicit sharpHermiteP3 dielectric model");
        require(std::isfinite(config.permittivityNegative) && config.permittivityNegative > 0.0 &&
                    std::isfinite(config.permittivityPositive) && config.permittivityPositive > 0.0,
                "Maxwell phase permittivities must be finite and positive");
        require(
            config.axisymmetric == panels.axisymmetric() &&
                config.axisymmetric == surfaceCharge.axisymmetric &&
                (!config.axisymmetric || (snapshot->radialAxisSymmetry() && snapshot->minCoord().x == 0.0)),
            "Maxwell charge, panels and physical measures/symmetry must agree");
        const bool native = !scope.has_value();
        require(snapshot->usesNativeLeaves() == native && panels.usesNativeLeaves() == native,
                "Maxwell entry point and retained scalar/panel layout disagree");
        std::vector<unsigned char> modelCell(mesh.cellIndexNum(),0);
        require(!config.nativeAxisNeckHarmonics || native,"axis-neck Maxwell reconstruction requires native geometry");
        if(axisCapModel) {
            require(native && config.axisymmetric && axisCapModel->source==snapshot
                && !axisCapModel->cells.empty(),"axis-cap traction model requires its current native RZ scalar");
            for(int cell:axisCapModel->cells) {
                require(cell>=0 && cell<mesh.cellIndexNum() && snapshot->nativeCellActive(cell)
                    && snapshot->nativeCellBounds(cell)[0].x==0 && !modelCell[cell],
                    "axis-cap traction cells must be distinct active leaves touching r=0");
                modelCell[cell]=1;
            }
        }
        const auto grid = native ? nullptr : std::make_unique<UniformLeafGrid>(mesh, *scope);
        // This constructor validates the identical actual lattice and supplies
        // only active native samples. It does not extrapolate solid potential.
        const Sampling::Lookup sampling =
            native ? Sampling::Lookup(mesh, snapshot,Sampling::NativeBoundaryWindow::OneSided)
                   : Sampling::Lookup(mesh, snapshot, *scope);
        const int count = mesh.cellIndexNum();
        require(absolutePotential.source() == snapshot && absolutePotential.size() == count &&
                    surfaceCharge.cellIntegratedCharge.size() == static_cast<std::size_t>(count) &&
                    static_cast<bool>(surfaceCharge.evaluate) &&
                    (native?panels.nativeChain().cellToGraphs().size():panels.cellToPanel().size())
                        == static_cast<std::size_t>(count),
                "missing or incompatible fitted Maxwell potential/charge/panel storage");
        std::vector<unsigned char> active(count, 0);
        std::vector<std::array<int, 2>> logicalCell(count);
        Eigen::VectorXd values = Eigen::VectorXd::Zero(count);
        if (native)
        {
            require(snapshot->nativeCellCount() == static_cast<std::size_t>(count),
                    "Maxwell native scalar has an incompatible leaf count");
            for (int cell = 0; cell < count; ++cell)
            {
                active[cell] = snapshot->nativeCellActive(cell) ? 1 : 0;
                require(std::isfinite(surfaceCharge.cellIntegratedCharge[cell]) &&
                            (active[cell] || surfaceCharge.cellIntegratedCharge[cell] == 0.0),
                        "non-finite prescribed charge or charge in an inactive cell");
                if (!active[cell])
                    continue;
                require(std::isfinite(absolutePotential[cell]), "non-finite active absolute potential");
                values[cell] = absolutePotential[cell];
            }
        }
        else
            for (int j = 0; j < grid->ny(); ++j)
                for (int i = 0; i < grid->nx(); ++i)
                {
                    const int label = grid->cells()[static_cast<std::size_t>(j) * grid->nx() + i];
                    const int cell = mesh.cells(label).index();
                    require(cell >= 0 && cell < count && mesh.cellIndexToLabel[cell] == label,
                            "stale native Maxwell cell index/label");
                    logicalCell[cell] = {{i, j}};
                    active[cell] = snapshot->cellActive(i, j) ? 1 : 0;
                    require(std::isfinite(surfaceCharge.cellIntegratedCharge[cell]) &&
                                (active[cell] || surfaceCharge.cellIntegratedCharge[cell] == 0.0),
                            "non-finite prescribed charge or charge in an inactive cell");
                    if (!active[cell])
                        continue; // Inactive potential can be NaN; it is not electric data.
                    require(std::isfinite(absolutePotential[cell]),
                            "non-finite active absolute potential");
                    values[cell] = absolutePotential[cell];
                }
        for (std::size_t p = 0; p < panels.panels().size(); ++p)
        {
            const int cell = panels.panels()[p].cellIndex;
            if (native)
            {
                require(cell >= 0 && cell < count && active[cell] &&
                            panels.nativeChain().graphs()[p].nativeCellIndex == cell,
                        "retained Maxwell graph does not match its active native leaf");
                continue;
            }
            const auto &graph = panels.panelGraph(static_cast<int>(p));
            require(cell >= 0 && cell < count && active[cell] &&
                        panels.cellToPanel()[cell] == static_cast<int>(p) &&
                        graph.cellI == logicalCell[cell][0] && graph.cellJ == logicalCell[cell][1],
                    "retained Maxwell panel does not match the actual active native cell");
        }
        const auto h = native ? panels.nativeChain().minimumActiveSpacing() : grid->spacing();
        std::size_t sample = 0;
        for (const auto &point : traceGeometry->points())
        {
            report.stage = "fitted-Maxwell-geometry-sample-" + std::to_string(sample++);
            const auto &frame = point.frame;
            require(frame.hasCellLocalPoint != frame.hasCellChartPoint,
                    "Maxwell load requires exactly one local or boundary-chart provenance");
            if (native)
            {
                require(frame.hasCellChartPoint && !frame.hasCellLocalPoint && frame.cellI == -1 &&
                            frame.cellJ == -1 && frame.nativeCellIndex >= 0 &&
                            frame.nativeCellIndex < count && active[frame.nativeCellIndex],
                        "native Maxwell frame lacks actual active leaf/chart provenance");
                const int panel = panels.nativeChain().graphAt(frame.nativeCellIndex,frame.cellChartPoint);
                require(panel >= 0 && panel < static_cast<int>(panels.panels().size()),
                        "native Maxwell frame has no retained free-interface incidence");
                panels.nativeChain().graphFraction(panel, frame.cellChartPoint);
            }
            else
            {
                require(frame.cellI >= 0 && frame.cellI < grid->nx() && frame.cellJ >= 0 &&
                            frame.cellJ < grid->ny() && snapshot->cellActive(frame.cellI, frame.cellJ),
                        "Maxwell load frame lacks actual active cell provenance");
                const int frameCell =
                    mesh.cells(
                            grid->cells()[static_cast<std::size_t>(frame.cellJ) * grid->nx() + frame.cellI])
                        .index();
                const int panel = panels.cellToPanel()[frameCell];
                require(panel >= 0 && panel < static_cast<int>(panels.panels().size()),
                        "Maxwell load frame has no retained free-interface panel incidence");
                const auto &graph = panels.panelGraph(panel);
                const double first = graph.coordinateIsX ? graph.localBegin.x : graph.localBegin.y;
                const double last = graph.coordinateIsX ? graph.localEnd.x : graph.localEnd.y;
                bool originalRoot = false, interior = false;
                if (frame.hasCellChartPoint)
                {
                    namespace Chart = Quantum::Geometry;
                    const auto point = LevelSetHermiteGeometry::cellChartPoint(frame);
                    const auto q = graph.coordinateIsX ? point.x : point.y;
                    originalRoot = Chart::cellChartEqual(point, graph.localBegin) ||
                                   Chart::cellChartEqual(point, graph.localEnd);
                    interior = Chart::cellChartInterior(point) &&
                               Chart::compareCellChartCoordinate(q, {std::min(first, last), false}) > 0 &&
                               Chart::compareCellChartCoordinate(q, {std::max(first, last), false}) < 0;
                }
                else
                {
                    const auto local = frame.cellLocalPoint;
                    const double q = graph.coordinateIsX ? local.x : local.y;
                    originalRoot = (local.x == graph.localBegin.x && local.y == graph.localBegin.y) ||
                                   (local.x == graph.localEnd.x && local.y == graph.localEnd.y);
                    interior = std::isfinite(local.x) && std::isfinite(local.y) && local.x > 0.0 &&
                               local.x < 1.0 && local.y > 0.0 && local.y < 1.0 && q > std::min(first, last) &&
                               q < std::max(first, last);
                }
                require(originalRoot || interior,
                        "Maxwell local/chart frame lies outside its original retained graph interval");
            }
            require(std::isfinite(frame.levelSetResidual) && frame.levelSetResidual >= 0.0 &&
                        std::isfinite(frame.stationarityResidual) && frame.stationarityResidual >= 0.0 &&
                        std::isfinite(point.endpointKktResidual) && point.endpointKktResidual >= 0.0 &&
                        std::isfinite(point.curvatureTangentialDerivative) &&
                        point.curvatureTangentialDerivative == retainedCurvatureDerivative(*snapshot, frame),
                    "invalid Maxwell local projection diagnostics or stale raw local curvature derivative");
            // Sampling::window below rechecks point, normal, curvature and
            // scalar residual against the matching local/chart jet. Projection selection and
            // its endpoint KKT condition remain the shared builder's contract;
            // do not reproject the rounded display point in this consumer.
        }
        FittedSharpDielectricTractions result;
        result.geometry = traceGeometry;
        result.interfaceTraces.reserve(traceGeometry->points().size());
        // Prepare independent windows and small matrix factorizations on this
        // immutable source. Charge callbacks stay serial: callers need not
        // make their own query implementation concurrently callable. Publish
        // results/errors only in the original requested-point order.
        std::vector<PreparedLoadTrace> prepared(traceGeometry->points().size());
        const auto prepare=[&](int i) {
            auto& item=prepared[i];const auto& point=traceGeometry->points()[i];
            try {
                item.window=sampling.window(point.frame);
                if(native && axisNeckWindow(mesh,*snapshot,{&panels},point.frame,item.window,config)) {
                    const auto bounds=snapshot->nativeCellBounds(point.frame.nativeCellIndex);
                    const double scale=std::max(bounds[1].x-bounds[0].x,bounds[1].y-bounds[0].y);
                    item.fit=std::make_unique<NativeTrace>(AxisTransmission::Reconstruction(point.frame.point,point.frame.normal,
                        point.frame.curvature,point.curvatureTangentialDerivative,config.permittivityNegative,
                        config.permittivityPositive,item.window.samples,scale));
                }
                else if(native && modelCell[point.frame.nativeCellIndex])
                    item.axisCap=prepareAxisCapFit(config,*snapshot,point,item.window,values);
                else {
                    double scale=std::max(h.x,h.y);
                    if(native) {const auto bounds=snapshot->nativeCellBounds(point.frame.nativeCellIndex);
                        scale=std::max(bounds[1].x-bounds[0].x,bounds[1].y-bounds[0].y);}
                    item.fit=std::make_unique<NativeTrace>(point.frame.point,point.frame.normal,
                        point.frame.curvature,point.curvatureTangentialDerivative,config.permittivityNegative,
                        config.permittivityPositive,item.window.samples,scale,config.axisymmetric);
                }
            } catch(...) {item.failure=std::current_exception();}
        };
        // Observe the caller's floating-point environment before worker tasks,
        // as in the native electric face preparation above.
        if(!prepared.empty())prepare(0);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 8) if(prepared.size()>32)
#endif
        for(int i=1;i<static_cast<int>(prepared.size());++i)prepare(i);
        sample = 0;
        for (const auto& point : traceGeometry->points())
        {
            report.stage = "fitted-Maxwell-trace-sample-" + std::to_string(sample);
            const auto& item=prepared[sample++];
            if(item.failure)std::rethrow_exception(item.failure);
            if(item.axisCap) {
                const auto trace=reconstructAxisCapTrace(config,point,*item.axisCap,surfaceCharge);
                ++report.localTraceFits;
                report.maxTraceCondition=std::max(report.maxTraceCondition,trace.condition);
                report.maxTraceReproductionError=std::max(report.maxTraceReproductionError,trace.reproduction);
                result.interfaceTraces.push_back(trace);
            } else result.interfaceTraces.push_back(reconstructTrace(config, *snapshot, point.frame,
                    point.curvatureTangentialDerivative, sampling, std::max(h.x, h.y), 0.0, true,
                    values, report, &surfaceCharge,&item));
        }
        auto& diagnostics = result.diagnostics;
        diagnostics.samplePoints = traceGeometry->points().size();
        diagnostics.localTraceFits = report.localTraceFits;
        diagnostics.maxTraceCondition = report.maxTraceCondition;
        diagnostics.maxTraceReproductionError = report.maxTraceReproductionError;
        for (const auto& trace : result.interfaceTraces)
        {
            if(trace.axisNeckHarmonicModel) {
                ++diagnostics.axisNeckHarmonicSamples;
                diagnostics.maxAxisNeckRelativeResidual=std::max(diagnostics.maxAxisNeckRelativeResidual,trace.modelRelativeResidual);
            }
            diagnostics.maxNormalTraction = std::max(diagnostics.maxNormalTraction, std::abs(trace.normalTraction));
            diagnostics.maxTangentialTraction = std::max(diagnostics.maxTangentialTraction, std::abs(trace.tangentialTraction));
            if(trace.axisCapExteriorModel) {
                ++diagnostics.axisCapExteriorSamples;
                diagnostics.maxAxisCapRelativeResidual=std::max(diagnostics.maxAxisCapRelativeResidual,trace.modelRelativeResidual);
            }
        }
        require(diagnostics.localTraceFits == diagnostics.samplePoints,
            "not every requested point received one actual-potential trace");
        return result;
    }
    catch (const std::exception &error)
    {
        throw std::runtime_error("[FittedSharpDielectricTractions] stage=" + report.stage + ": " + error.what());
    }
}

FittedSharpDielectricTractions reconstructFittedSharpDielectricTractions(
    const QuadtreeMesh &mesh, std::shared_ptr<const InterfaceTraceGeometry> geometry,
    const LevelSetInterfacePanels &panels, CartesianGridScope scope, const DielectricConfig &config,
    const FittedElectricPotential &potential, const SurfaceChargeSource &charge)
{
    return reconstructFittedTractions(mesh, std::move(geometry), panels, scope, config, potential, charge);
}

FittedSharpDielectricTractions reconstructNativeFittedSharpDielectricTractions(
    const QuadtreeMesh &mesh, std::shared_ptr<const InterfaceTraceGeometry> geometry,
    const LevelSetInterfacePanels &panels, const DielectricConfig &config, const FittedElectricPotential &potential,
    const SurfaceChargeSource &charge,const FittedAxisCapTractionModel* axisCapModel)
{
    return reconstructFittedTractions(mesh, std::move(geometry), panels, std::nullopt, config, potential, charge,axisCapModel);
}
