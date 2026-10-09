#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"

#include <limits>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace
{
template<class Field, class Labels>
void checkNative(const Field& field, int count, Labels labels, bool finite = true)
{
    if (field.size() != count || !field.value)
        throw std::invalid_argument("Hydro state field does not match the native mesh");
    for (int i = 0; i < count; ++i)
        if (field.label(i) != labels[i] || (finite && !std::isfinite(field.value[i])))
            throw std::invalid_argument("Hydro state field has stale labels or non-finite values");
}

void checkMaterial(double density, double viscosity)
{
    if (!std::isfinite(density) || density <= 0.0 || !std::isfinite(viscosity) || viscosity < 0.0)
        throw std::invalid_argument("Hydro requires finite density > 0 and viscosity >= 0");
}
}

void Hydrodynamics::validateStateGeometry(const QuadtreeMesh& mesh, const HydroStepGeometry& geometry) const
{
    const bool fitted = geometry.kind == HydroDomainKind::FittedRZ;
    if ((!fitted && geometry.kind != HydroDomainKind::CompleteXY && geometry.kind != HydroDomainKind::CompleteRZ)
        || fitted != (geometry.fitted != nullptr) || mesh.cellIndexNum() <= 0)
        throw std::invalid_argument("invalid Hydro state geometry capability");
    if (stepProgress_.geometry && *stepProgress_.geometry != geometry.kind)
        throw std::invalid_argument("cannot change Hydro state coordinate/domain capability");
    if (stepMeshIdentity_) stepMeshIdentity_->validate(mesh);
    if (fitted)
    {
        if (!geometry.fitted->axisymmetric()) throw std::invalid_argument("fitted Hydro state requires RZ");
        geometry.fitted->validateMeshIdentity(mesh);
        if (fittedDomain_ && !fittedDomain_->sameDefinition(*geometry.fitted))
            throw std::invalid_argument("Hydro state fluid domain/boundaries changed");
    }
    checkNative(p, mesh.cellIndexNum(), mesh.cellIndexToLabel, false);
    checkNative(ux, mesh.vEdgeIndexNum(), mesh.vEdgeIndexToLabel, false);
    checkNative(uy, mesh.hEdgeIndexNum(), mesh.hEdgeIndexToLabel, false);
}

void Hydrodynamics::initializeCompleteVelocity(QuadtreeMesh& mesh,
    const HydroStepGeometry& geometry, HydroInitialVelocity profile, double amplitude,
    const CellField<double>* initialPressure)
{
    if (geometry.kind != HydroDomainKind::CompleteXY && geometry.kind != HydroDomainKind::CompleteRZ)
        throw std::invalid_argument("prescribed Hydro initial velocity requires a complete domain");
    if (profile != HydroInitialVelocity::Zero && profile != HydroInitialVelocity::Smooth && profile!=HydroInitialVelocity::Vortex)
        throw std::invalid_argument("unknown Hydro initial velocity profile");
    if (!std::isfinite(amplitude) || amplitude < 0.0)
        throw std::invalid_argument("Hydro initial velocity amplitude must be finite and nonnegative");
    if (profile == HydroInitialVelocity::Zero)
    {
        initializeState(mesh, geometry, {initialPressure, nullptr, nullptr});
        return;
    }
    VEdgeField<double> initialX(mesh, 0.0);
    HEdgeField<double> initialY(mesh, 0.0);
    const double pi = std::acos(-1.0);
    if(profile==HydroInitialVelocity::Vortex) {
        const bool rz=geometry.kind==HydroDomainKind::CompleteRZ;
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
        const double lx=mesh.maxCoord.x-mesh.minCoord.x,ly=mesh.maxCoord.y-mesh.minCoord.y;
        const auto stream=[&](Vector2D p){
            if(p.x==mesh.minCoord.x || p.x==mesh.maxCoord.x || p.y==mesh.minCoord.y || p.y==mesh.maxCoord.y)return 0.;
            const double x=(p.x-mesh.minCoord.x)/lx,y=(p.y-mesh.minCoord.y)/ly;
            return rz?2*pi*amplitude*lx*lx*x*x*std::pow(1-x*x,3)*std::sin(pi*y)
                :amplitude*lx/pi*(std::sin(pi*x)*std::sin(pi*y)+.17*std::sin(2*pi*x)*std::sin(3*pi*y));};
        for(const auto face:domain.faces(mesh)){
            const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            (face.horizontal?initialY.value:initialX.value)[edge.index()]=face.outwardSign?0
                :(face.horizontal?1:-1)*(stream(b)-stream(a))/face.area;
        }
        initializeState(mesh,geometry,{initialPressure,&initialX,&initialY});return;
    }
    const auto sample = [&](const Vector2D& point, bool radial) {
        const double x = (point.x - mesh.minCoord.x) / (mesh.maxCoord.x - mesh.minCoord.x);
        const double y = (point.y - mesh.minCoord.y) / (mesh.maxCoord.y - mesh.minCoord.y);
        if (radial) return amplitude * std::sin(pi * x) * std::sin(pi * y);
        // Even at r=0 for the RZ axial component; a non-solenoidal smooth seed
        // deliberately exercises pressure projection rather than a zero RHS.
        return 0.5 * amplitude * std::cos(pi * x) * std::sin(2.0 * pi * y);
    };
    for (int i = 0; i < mesh.vEdgeIndexNum(); ++i)
    {
        const Edge& edge = mesh.vEdges(mesh.vEdgeIndexToLabel[i]);
        initialX.value[i] = edge.boundary() == Boundary::non ? sample(edge.midpointCoord(), true) : 0.0;
    }
    for (int i = 0; i < mesh.hEdgeIndexNum(); ++i)
    {
        const Edge& edge = mesh.hEdges(mesh.hEdgeIndexToLabel[i]);
        initialY.value[i] = edge.boundary() == Boundary::non ? sample(edge.midpointCoord(), false) : 0.0;
    }
    initializeState(mesh, geometry, {initialPressure, &initialX, &initialY});
}

void Hydrodynamics::initializeState(QuadtreeMesh& mesh, const HydroStepGeometry& geometry,
    const HydroInitialFields& initial)
{
    if (initialStateSet_ || stepProgress_.stopped || stepEntry_ == StepEntry::Legacy
        || stepProgress_.completedSteps != 0 || fixedDt_ != 0.0 || pressureMode_ != PressureMode::Unset)
        throw std::logic_error("Hydro initial fields must be supplied once before any time-step history");
    validateStateGeometry(mesh, geometry);
    if ((initial.velocityX == nullptr) != (initial.velocityY == nullptr))
        throw std::invalid_argument("Hydro initial velocity requires both native components");
    if (initial.pressure) checkNative(*initial.pressure, mesh.cellIndexNum(), mesh.cellIndexToLabel);
    if (initial.velocityX)
    {
        checkNative(*initial.velocityX, mesh.vEdgeIndexNum(), mesh.vEdgeIndexToLabel);
        checkNative(*initial.velocityY, mesh.hEdgeIndexNum(), mesh.hEdgeIndexToLabel);
    }
    for (const auto* field : {&uxPre, &uxMid}) checkNative(*field, mesh.vEdgeIndexNum(), mesh.vEdgeIndexToLabel, false);
    for (const auto* field : {&uyPre, &uyMid}) checkNative(*field, mesh.hEdgeIndexNum(), mesh.hEdgeIndexToLabel, false);
    for (const auto* field : {&divU, &uMag}) checkNative(*field, mesh.cellIndexNum(), mesh.cellIndexToLabel, false);
    const bool fitted = geometry.fitted != nullptr;
    auto domain = fittedDomain_;
    auto identity = stepMeshIdentity_;
    if (fitted && !domain) domain = std::make_shared<const HydroFittedDomain>(*geometry.fitted);
    if (!fitted && !identity) identity = std::make_shared<const NativeMeshIdentity>(mesh);

    VEdgeField<double> nextX(mesh, 0.0);
    HEdgeField<double> nextY(mesh, 0.0);
    CellField<double> nextDivergence(mesh, 0.0), nextMagnitude(mesh, 0.0);
    if (initial.velocityX)
    {
        for (int i = 0; i < nextX.size(); ++i) nextX.value[i] = initial.velocityX->value[i];
        for (int i = 0; i < nextY.size(); ++i) nextY.value[i] = initial.velocityY->value[i];
        if (fitted)
        {
            Eigen::VectorXd velocity(domain->nativeByVelocityRow().size());
            for (const auto& face : domain->faces())
            {
                const double value = face.horizontal ? nextY.value[face.nativeRow - nextX.size()] : nextX.value[face.nativeRow];
                if (face.velocityRow >= 0) velocity[face.velocityRow] = value;
                else if (value != 0.0) throw std::invalid_argument("constrained fitted initial velocity must be zero");
            }
            const Eigen::VectorXd flux = domain->divergence() * velocity;
            for (int row = 0; row < flux.size(); ++row)
            {
                const int i = domain->cellIndexByRow()[row];
                const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
                nextDivergence.value[i] = flux[row] / domain->cellVolumes()[row];
                const double ur = .5 * (nextX.value[mesh.vEdges(cell.cellToVEdge[0][0]).index()]
                    + nextX.value[mesh.vEdges(cell.cellToVEdge[1][0]).index()]);
                const double uz = .5 * (nextY.value[mesh.hEdges(cell.cellToHEdge[0][0]).index()]
                    + nextY.value[mesh.hEdges(cell.cellToHEdge[1][0]).index()]);
                nextMagnitude.value[i] = std::hypot(ur, uz);
            }
        }
        else
        {
            bool hangingFaces = false;
            for (int i = 0; i < mesh.cellIndexNum(); ++i)
            {
                const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
                for (int side = 0; side < 2; ++side)
                    hangingFaces |= cell.cellToVEdge[side][1] != -1 || cell.cellToHEdge[side][1] != -1;
            }
            const bool rz = geometry.kind == HydroDomainKind::CompleteRZ;
            if (hangingFaces)
                velocityDiagnostics(mesh, nextX, nextY, rz, true, &nextDivergence, &nextMagnitude);
            else
                // Exact original uniform initial-field diagnostic arithmetic.
                // The general native diagnostic remains the existing AMR path.
                for (int i = 0; i < mesh.cellIndexNum(); ++i)
                {
                    const auto& cell = mesh.cells(mesh.cellIndexToLabel[i]);
                    const auto& west = mesh.vEdges(cell.cellToVEdge[0][0]);
                    const auto& east = mesh.vEdges(cell.cellToVEdge[1][0]);
                    const auto& south = mesh.hEdges(cell.cellToHEdge[0][0]);
                    const auto& north = mesh.hEdges(cell.cellToHEdge[1][0]);
                    const double uw = nextX.value[west.index()], ue = nextX.value[east.index()];
                    const double vs = nextY.value[south.index()], vn = nextY.value[north.index()];
                    const double radialFlux = rz
                        ? (east.midpointCoord().x * ue - west.midpointCoord().x * uw) / cell.centerCoord().x
                        : ue - uw;
                    nextDivergence.value[i] = radialFlux / cell.size().x + (vn - vs) / cell.size().y;
                    nextMagnitude.value[i] = std::hypot(0.5 * (uw + ue), 0.5 * (vs + vn));
                }
        }
    }
    checkNative(nextDivergence, mesh.cellIndexNum(), mesh.cellIndexToLabel);
    checkNative(nextMagnitude, mesh.cellIndexNum(), mesh.cellIndexToLabel);
    // Publish only after all checks/allocations. Keep p's actual initial value;
    // neither an analytic pressure nor a projection multiplier is substituted.
    for (int i = 0; i < p.size(); ++i)
    {
        p.value[i] = initial.pressure ? initial.pressure->value[i] : 0.0;
        divU.value[i] = nextDivergence.value[i]; uMag.value[i] = nextMagnitude.value[i];
    }
    for (int i = 0; i < ux.size(); ++i) ux.value[i] = uxPre.value[i] = uxMid.value[i] = nextX.value[i];
    for (int i = 0; i < uy.size(); ++i) uy.value[i] = uyPre.value[i] = uyMid.value[i] = nextY.value[i];
    fittedDomain_ = std::move(domain); stepMeshIdentity_ = std::move(identity);
    stepProgress_.geometry = geometry.kind;
    initialStateSet_ = true; stepEntry_ = StepEntry::Managed;
}

void Hydrodynamics::updateMaterials(QuadtreeMesh& mesh, const HydroMaterialUpdate& input)
{
    if (stepProgress_.stopped || stepEntry_ == StepEntry::Legacy)
        throw std::logic_error("cannot update managed materials after a failed or legacy Hydro step");
    validateStateGeometry(mesh, input.geometry);
    const double allowance = 64.0 * std::numeric_limits<double>::epsilon()
        * std::max({std::abs(stepProgress_.endTime), fixedDt_, std::abs(input.sourceTime)});
    if (!std::isfinite(input.sourceTime) || input.sourceTime < 0.0
        || std::abs(input.sourceTime - stepProgress_.endTime) > allowance)
        throw std::invalid_argument("Hydro material source time must match the current state time");
    if (input.evolution != HydroMaterialEvolution::Frozen && input.evolution != HydroMaterialEvolution::TimeDependent)
        throw std::invalid_argument("invalid Hydro material evolution contract");
    const auto revision = materialState_ ? materialState_->revision : 0;
    if (revision == std::numeric_limits<std::uint64_t>::max()) throw std::overflow_error("Hydro material revision overflow");
    std::visit([&](const auto& material)
    {
        using Material = std::decay_t<decltype(material)>;
        if constexpr (std::is_same_v<Material, HydroConstantMaterial>)
            checkMaterial(material.density, material.viscosity);
        else if constexpr (std::is_same_v<Material, HydroNativeMaterial>)
        {
            if (!material.density || !material.viscosity || !material.vertexDensity || !material.vertexViscosity)
                throw std::invalid_argument("native Hydro materials require all four coefficient fields");
            checkNative(*material.density, mesh.cellIndexNum(), mesh.cellIndexToLabel);
            checkNative(*material.viscosity, mesh.cellIndexNum(), mesh.cellIndexToLabel);
            checkNative(*material.vertexDensity, mesh.vertIndexNum(), mesh.vertIndexToLabel);
            checkNative(*material.vertexViscosity, mesh.vertIndexNum(), mesh.vertIndexToLabel);
            for (int i = 0; i < mesh.cellIndexNum(); ++i)
                if (!input.geometry.fitted || input.geometry.fitted->cellRowByIndex()[i] >= 0)
                    checkMaterial(material.density->value[i], material.viscosity->value[i]);
            for (int i = 0; i < mesh.vertIndexNum(); ++i)
                checkMaterial(material.vertexDensity->value[i], material.vertexViscosity->value[i]);
        }
        else if constexpr (std::is_same_v<Material, HydroLevelSetMaterial>)
        {
            if (input.geometry.fitted || !material.levelSet || !std::isfinite(material.thickness) || material.thickness <= 0.0)
                throw std::invalid_argument("vertex-mixture materials require a complete domain and a finite positive thickness");
            checkNative(*material.levelSet, mesh.vertIndexNum(), mesh.vertIndexToLabel);
            checkMaterial(material.densityNegative, material.viscosityNegative);
            checkMaterial(material.densityPositive, material.viscosityPositive);
        }
        else if constexpr (std::is_same_v<Material, HydroConservativeMaterial> || std::is_same_v<Material, HydroParametricMaterial>)
        {
            if (input.geometry.fitted || !material.domain || !material.inventory)
                throw std::invalid_argument("conservative Hydro materials require a complete native domain and inventory");
            if constexpr (std::is_same_v<Material, HydroParametricMaterial>)
                InterfaceParametricMotion::validateDomain(mesh,*material.domain,*material.inventory);
            else material.inventory->validate(mesh,*material.domain);
            if (material.domain->cellIndexByRow().size()!=static_cast<std::size_t>(mesh.cellIndexNum())
                || material.domain->axisymmetric()!=(input.geometry.kind==HydroDomainKind::CompleteRZ))
                throw std::invalid_argument("conservative Hydro material domain differs from the complete XY/RZ solver");
            const auto source=material.inventory->progress();
            if (source.completedSteps!=stepProgress_.completedSteps || std::abs(source.endTime-input.sourceTime)>allowance)
                throw std::invalid_argument("conservative Hydro material inventory has a different physical step/time");
            checkMaterial(material.densityNegative, material.viscosityNegative);
            checkMaterial(material.densityPositive, material.viscosityPositive);
        }
        else
        {
            if (!input.geometry.fitted || !material.phaseByCell
                || material.phaseByCell->size() != static_cast<std::size_t>(mesh.cellIndexNum())
                || material.first.phaseId == material.second.phaseId)
                throw std::invalid_argument("classified Hydro materials require a fitted domain and two bound phases");
            checkMaterial(material.first.density, material.first.viscosity);
            checkMaterial(material.second.density, material.second.viscosity);
            for (int i : input.geometry.fitted->cellIndexByRow())
                if ((*material.phaseByCell)[i] != material.first.phaseId && (*material.phaseByCell)[i] != material.second.phaseId)
                    throw std::invalid_argument("active fluid cell has an unbound material phase");
        }
    }, input.model);
    auto domain = fittedDomain_;
    auto identity = stepMeshIdentity_;
    if (input.geometry.fitted && !domain) domain = std::make_shared<const HydroFittedDomain>(*input.geometry.fitted);
    if (!input.geometry.fitted && !identity) identity = std::make_shared<const NativeMeshIdentity>(mesh);
    stepEntry_ = StepEntry::Managed;
    try
    {
        std::visit([&](const auto& material)
        {
            using Material = std::decay_t<decltype(material)>;
            if constexpr (std::is_same_v<Material, HydroConstantMaterial>)
                assignConstantMaterialFields(mesh, material.density, material.viscosity);
            else if constexpr (std::is_same_v<Material, HydroNativeMaterial>)
            {
                // Snapshot first: callers may supply views of existing fields,
                // including swapped buffers. No source may change mid-copy.
                CellField<double> nextDensity(*material.density), nextViscosity(*material.viscosity);
                VertField<double> nextVertexDensity(*material.vertexDensity), nextVertexViscosity(*material.vertexViscosity);
                rho = nextDensity;
                mu = nextViscosity;
                vertRho = nextVertexDensity;
                vertMu = nextVertexViscosity;
            }
            else if constexpr (std::is_same_v<Material, HydroLevelSetMaterial>)
                assignLevelSetMaterialFields(mesh, *material.levelSet, material.densityNegative, material.densityPositive,
                    material.viscosityNegative, material.viscosityPositive, material.thickness);
            else if constexpr (std::is_same_v<Material, HydroConservativeMaterial> || std::is_same_v<Material, HydroParametricMaterial>)
            {
                const auto& phase=*material.inventory;
                rho.initialize(mesh,0);mu.initialize(mesh,0);
                vertRho.initialize(mesh,0);vertMu.initialize(mesh,0);
                std::vector<int> slotByCell;
                if constexpr (std::is_same_v<Material, HydroParametricMaterial>) {
                    slotByCell.resize(mesh.cellIndexNum());const auto& cells=InterfaceParametricMotion::cellIndices(phase);
                    for(std::size_t slot=0;slot<cells.size();++slot)slotByCell[cells[slot]]=static_cast<int>(slot);
                }
                for (int c:material.domain->cellIndexByRow()) {
                    const long double capacity=material.domain->cellVolumes()[material.domain->cellRowByIndex()[c]];
                    const long double liquid=[&]() -> long double {
                        if constexpr (std::is_same_v<Material, HydroParametricMaterial>)return phase.parametricGeometry().inventory()[slotByCell[c]];
                        else return phase.liquidVolumes()[c];
                    }();
                    const long double fraction=liquid/capacity;
                    // Convex sums preserve the pure phase limits and avoid a
                    // cancellation-prone density difference at large ratios.
                    rho.value[c]=static_cast<double>(fraction*material.densityNegative+(1-fraction)*material.densityPositive);
                    mu.value[c]=static_cast<double>(fraction*material.viscosityNegative+(1-fraction)*material.viscosityPositive);
                    checkMaterial(rho.value[c],mu.value[c]);
                }
                for (int v=0;v<mesh.vertIndexNum();++v) {
                    long double density=0,viscosity=0;int count=0;
                    for (int label:mesh.verts(mesh.vertIndexToLabel[v]).vertToCell) if(label>=0) {
                        const int c=mesh.cells(label).index();density+=rho.value[c];viscosity+=mu.value[c];++count;
                    }
                    if(!count) throw std::invalid_argument("conservative Hydro material vertex has no incident cell");
                    vertRho.value[v]=static_cast<double>(density/count);vertMu.value[v]=static_cast<double>(viscosity/count);
                }
            }
            else
                assignHydroFittedMaterials(mesh, *input.geometry.fitted, {rho, mu, vertRho, vertMu},
                    *material.phaseByCell, material.first, material.second);
        }, input.model);
        checkNative(rho, mesh.cellIndexNum(), mesh.cellIndexToLabel);
        checkNative(mu, mesh.cellIndexNum(), mesh.cellIndexToLabel);
        checkNative(vertRho, mesh.vertIndexNum(), mesh.vertIndexToLabel);
        checkNative(vertMu, mesh.vertIndexNum(), mesh.vertIndexToLabel);
        for (int i = 0; i < rho.size(); ++i)
            if (!domain || domain->cellRowByIndex()[i] >= 0) checkMaterial(rho.value[i], mu.value[i]);
        const bool classified = std::holds_alternative<HydroClassifiedMaterial>(input.model);
        for (int i = 0; i < vertRho.size(); ++i)
            if (classified && vertRho.value[i] == 0.0 && vertMu.value[i] >= 0.0) continue;
            else checkMaterial(vertRho.value[i], vertMu.value[i]);
    }
    catch (...)
    {
        stepProgress_.stopped = true;
        throw;
    }
    fittedDomain_ = std::move(domain); stepMeshIdentity_ = std::move(identity);
    stepProgress_.geometry = input.geometry.kind;
    materialState_ = HydroMaterialState{revision + 1, input.sourceTime, input.evolution};
    stepEntry_ = StepEntry::Managed;
}
