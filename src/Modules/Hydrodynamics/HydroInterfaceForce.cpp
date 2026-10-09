#include "Modules/Hydrodynamics/HydroInterfaceForce.h"
#include "Modules/Hydrodynamics/HydroBalancedForce.h"
#include "Core/Field/FieldChecks.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Modules/Hydrodynamics/HydroInterfaceLoadGeometry.h"
#include "Modules/Interface/Compatibility/UniformInterfaceContour.h"
#include "Modules/Hydrodynamics/HydroNodalTraction.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <vector>

namespace HydroInterfaceForce
{
using FieldChecks::checkFinite;

// Pre-advection capillarity/traction kernel. Closed-node quadrature remains
// reference-only; active sigma is applied on faces once by the owning assembly.
static void build(QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>& geometry,
    int interfaceIndex, bool rz, bool axisContact,
    const Options& options, CellField<double>& forceX, CellField<double>& forceY,
    VEdgeField<double>& faceForceX, HEdgeField<double>& faceForceY,
    Report& diagnostic, std::string& stage, const HydroInterfaceLoadGeometry* loadGeometry);

Assembly::Assembly(QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>* geometry,
    int interfaceIndex, bool rz, bool axisContact, const Options& options,
    Report& diagnostic, std::string& stage, const HydroInterfaceLoadGeometry* loadGeometry)
    : mesh_(mesh), interfaceEnabled_(options.surfaceTension != 0.0 || options.normalTraction != 0.0
        || options.tangentialTraction != 0.0), hasFaces_(options.surfaceTension != 0.0 || loadGeometry)
{
    // Preserve the explicit Hermite zero-load path and its original clearing
    // order. A later coupling result must not select a different inlet.
    faceX_.initialize(mesh, 0.0);
    faceY_.initialize(mesh, 0.0);
    if (loadGeometry && !interfaceEnabled_)
    {
        cellX_.initialize(mesh, 0.0);
        cellY_.initialize(mesh, 0.0);
    }
    if (interfaceEnabled_ || axisContact || options.boundaryContact)
    {
        if (!geometry || interfaceIndex < 0 || static_cast<std::size_t>(interfaceIndex) >= geometry->levelSets.size())
            throw std::invalid_argument("interface load assembly requires its bound geometry");
        build(mesh, *geometry, interfaceIndex, rz, axisContact, options,
            cellX_, cellY_, faceX_, faceY_, diagnostic, stage, loadGeometry);
    }
}

void Assembly::requireOpen() const
{
    if (failed_ || published_) throw std::logic_error("Hydro load assembly is failed or already published");
}

void Assembly::addCouplingCells(const CellField<double>& x, const CellField<double>& y)
{
    requireOpen();
    if (cellsAdded_) throw std::logic_error("Hydro coupling cells may be added only once");
    failed_ = true;
    checkFinite(x, "coupling force density x", mesh_.cellIndexNum());
    checkFinite(y, "coupling force density y", mesh_.cellIndexNum());
    if (!interfaceEnabled_)
    {
        cellX_.initialize(mesh_, 0.0);
        cellY_.initialize(mesh_, 0.0);
    }
    for (int i = 0; i < mesh_.cellIndexNum(); ++i)
    {
        cellX_.value[i] += x.value[i];
        cellY_.value[i] += y.value[i];
    }
    checkFinite(cellX_, "combined force density x", mesh_.cellIndexNum());
    checkFinite(cellY_, "combined force density y", mesh_.cellIndexNum());
    cellsAdded_ = true;
    failed_ = false;
}

void Assembly::addCouplingFaces(const VEdgeField<double>& x, const HEdgeField<double>& y, bool supplied)
{
    requireOpen();
    if (!cellsAdded_ || facesAdded_) throw std::logic_error("Hydro coupling faces require one preceding cell result");
    failed_ = true;
    checkFinite(x, "coupling face force density x", mesh_.vEdgeIndexNum());
    checkFinite(y, "coupling face force density y", mesh_.hEdgeIndexNum());
    for (int i = 0; i < mesh_.vEdgeIndexNum(); ++i)
    {
        if (!supplied && x.value[i] != 0.0)
            throw std::runtime_error("coupling returned false with nonzero x face force");
        faceX_.value[i] += x.value[i];
    }
    for (int i = 0; i < mesh_.hEdgeIndexNum(); ++i)
    {
        if (!supplied && y.value[i] != 0.0)
            throw std::runtime_error("coupling returned false with nonzero y face force");
        faceY_.value[i] += y.value[i];
    }
    checkFinite(faceX_, "combined face force density x", mesh_.vEdgeIndexNum());
    checkFinite(faceY_, "combined face force density y", mesh_.hEdgeIndexNum());
    hasFaces_ = hasFaces_ || supplied;
    facesAdded_ = true;
    failed_ = false;
}

const VEdgeField<double>& Assembly::faceX() const
{
    if (!published_) throw std::logic_error("Hydro load diagnostics require publication");
    return faceX_;
}

const HEdgeField<double>& Assembly::faceY() const
{
    if (!published_) throw std::logic_error("Hydro load diagnostics require publication");
    return faceY_;
}

HydroForceDensity Assembly::publish(double evaluationTime)
{
    requireOpen();
    if (cellsAdded_ && !facesAdded_) throw std::logic_error("Hydro coupling face stage is incomplete");
    if (!std::isfinite(evaluationTime) || evaluationTime < 0.0)
        throw std::invalid_argument("Hydro load evaluation time must be finite and nonnegative");
    HydroForceDensity load;
    load.evaluationTime = evaluationTime;
    if (hasFaces_ || interfaceEnabled_ || cellsAdded_) { load.cellX = &cellX_; load.cellY = &cellY_; }
    if (hasFaces_) { load.faceX = &faceX_; load.faceY = &faceY_; }
    published_ = true;
    return load;
}

static void build(QuadtreeMesh& mesh, const ImmersedGeometry<QuadtreeMesh>& geometry,
    int interfaceIndex, bool rz, bool axisContact,
    const Options& options, CellField<double>& forceX, CellField<double>& forceY,
    VEdgeField<double>& faceForceX, HEdgeField<double>& faceForceY,
    Report& diagnostic, std::string& stage,
    const HydroInterfaceLoadGeometry* loadGeometry)
{
    const bool boundaryContact = options.boundaryContact;
    if (axisContact || boundaryContact)
    {
        stage = boundaryContact ? "interface-boundary-force" : "interface-axis-force";
        diagnostic.legacyReferenceAvailable = false;
        if (!loadGeometry || (!boundaryContact && (!rz || !options.axisContactSupported)))
            throw std::runtime_error("axis-contact capillarity requires RZ regularized shared Hermite geometry without reinitialization");
        if (options.normalTraction != 0.0 || options.tangentialTraction != 0.0)
            throw std::runtime_error("prescribed normal/tangential traction on an axis-contact interface is unsupported");
        // The legacy closed-node/Peskin reference is undefined on these open
        // axis/wall chains. Do not manufacture reference geometry or forces.
        forceX.initialize(mesh, 0.0);
        forceY.initialize(mesh, 0.0);
        stage = "interface-balanced-force";
        const auto active = HydroBalancedForce::buildCapillaryForce(mesh, geometry.levelSets[interfaceIndex],
            *loadGeometry, options.surfaceTension, faceForceX, faceForceY);
        diagnostic.activeCapillaryFaceX = active.forceX;
        diagnostic.activeCapillaryFaceY = active.forceY;
        diagnostic.maxActiveCapillaryFaceDensity = active.maxAbsComponent;
        return;
    }
    stage = "interface-geometry";
    const auto contour = UniformInterfaceContour::build(mesh, geometry, interfaceIndex, rz);
    const auto& geometryReport = contour.report();
    Field<double> tractionX(geometryReport.nodeCount, 0.0), tractionY(geometryReport.nodeCount, 0.0);
    diagnostic.nodeCount = geometryReport.nodeCount;
    diagnostic.surfaceMeasure = geometryReport.surfaceMeasure;
    diagnostic.maxAbsCurvature = geometryReport.maxAbsCurvature;
    if (diagnostic.nodeCount <= 0 || !std::isfinite(diagnostic.surfaceMeasure) ||
        !(diagnostic.surfaceMeasure > 0.0) || !std::isfinite(diagnostic.maxAbsCurvature))
        throw std::runtime_error("interface geometry returned invalid node/measure/curvature diagnostics");

    stage = "interface-traction";
    const auto& nodes = contour.nodes();
    if (nodes.size() != static_cast<std::size_t>(diagnostic.nodeCount))
        throw std::runtime_error("interface geometry node count does not match its label map");
    checkFinite(tractionX, "interface traction x", diagnostic.nodeCount);
    checkFinite(tractionY, "interface traction y", diagnostic.nodeCount);
    for (const auto& node : nodes)
    {
        const int index = node.index();
        const auto normal = node.normal();
        const double normalTraction = options.surfaceTension * node.curvature() + options.normalTraction;
        if (index < 0 || index >= diagnostic.nodeCount || !std::isfinite(normal.x) ||
            !std::isfinite(normal.y) || !std::isfinite(normalTraction))
            throw std::runtime_error("invalid interface node index, normal or normal traction");
        // The signed curvature follows the geometry normal. Tangent is (-ny,nx).
        // Prescribed traction is independent of surface tension and may later be
        // supplied by an electric module through the same force-density port.
        tractionX.value[index] =
            normalTraction * normal.x - options.tangentialTraction * normal.y;
        tractionY.value[index] =
            normalTraction * normal.y + options.tangentialTraction * normal.x;
    }
    checkFinite(tractionX, "interface traction x", diagnostic.nodeCount);
    checkFinite(tractionY, "interface traction y", diagnostic.nodeCount);

    stage = "interface-force-spreading";
    // Preserve the original node and cell/Peskin quadrature definitions in the
    // existing CSV columns. These are explicitly LEGACY REFERENCE diagnostics,
    // not the active sigma load after switching capillarity to staggered faces.
    CellField<double> legacyForceX(mesh, 0.0), legacyForceY(mesh, 0.0);
    const auto forceReport = HydroNodalTraction::spread(mesh, contour, tractionX, tractionY, legacyForceX, legacyForceY);
    checkFinite(legacyForceX, "legacy reference interface force density x", mesh.cellIndexNum());
    checkFinite(legacyForceY, "legacy reference interface force density y", mesh.cellIndexNum());
    diagnostic.forceX = forceReport.forceX;
    diagnostic.forceY = forceReport.forceY;
    diagnostic.spreadForceX = forceReport.spreadForceX;
    diagnostic.spreadForceY = forceReport.spreadForceY;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
        diagnostic.maxForceDensity = std::max(diagnostic.maxForceDensity,
            std::hypot(legacyForceX.value[i], legacyForceY.value[i]));
    if (!std::isfinite(diagnostic.forceX) || !std::isfinite(diagnostic.forceY) ||
        !std::isfinite(diagnostic.spreadForceX) || !std::isfinite(diagnostic.spreadForceY) ||
        !std::isfinite(diagnostic.maxForceDensity))
        throw std::runtime_error("interface force integrals or maximum force density are non-finite");

    // Only prescribed non-sigma normal/tangential traction remains in the
    // active cell-force port. Never apply the old cell sigma force as well.
    stage = "interface-prescribed-traction";
    forceX.initialize(mesh, 0.0);
    forceY.initialize(mesh, 0.0);
    for (const auto& node : nodes)
    {
        const auto normal = node.normal();
        tractionX.value[node.index()] =
            options.normalTraction * normal.x - options.tangentialTraction * normal.y;
        tractionY.value[node.index()] =
            options.normalTraction * normal.y + options.tangentialTraction * normal.x;
    }
    if (options.normalTraction != 0.0 || options.tangentialTraction != 0.0)
        HydroNodalTraction::spread(mesh, contour, tractionX, tractionY, forceX, forceY);
    checkFinite(forceX, "active prescribed traction density x", mesh.cellIndexNum());
    checkFinite(forceY, "active prescribed traction density y", mesh.cellIndexNum());

    stage = "interface-balanced-force";
    if (loadGeometry)
    {
        const auto active = HydroBalancedForce::buildCapillaryForce(mesh, geometry.levelSets[interfaceIndex],
            *loadGeometry, options.surfaceTension, faceForceX, faceForceY);
        diagnostic.activeCapillaryFaceX = active.forceX;
        diagnostic.activeCapillaryFaceY = active.forceY;
        diagnostic.maxActiveCapillaryFaceDensity = active.maxAbsComponent;
        return;
    }
    if (options.surfaceTension != 0.0)
    {
        std::vector<HydroBalancedForce::Node> balanceNodes;
        balanceNodes.reserve(nodes.size());
        std::set<int> visited;
        int current = nodes.front().label();
        for (std::size_t i = 0; i < nodes.size(); ++i)
        {
            if (!visited.insert(current).second)
                throw std::runtime_error("capillary contour repeats before all nodes are visited");
            const auto& node = nodes[current];
            balanceNodes.push_back({node.coord(), node.curvature()});
            current = node.upNodeLabel();
        }
        if (current != nodes.front().label()) throw std::runtime_error("capillary contour does not close");
        const auto active = HydroBalancedForce::buildCapillaryForce(mesh, geometry.levelSets[interfaceIndex],
            balanceNodes, options.surfaceTension, rz, faceForceX, faceForceY);
        diagnostic.activeCapillaryFaceX = active.forceX;
        diagnostic.activeCapillaryFaceY = active.forceY;
        diagnostic.maxActiveCapillaryFaceDensity = active.maxAbsComponent;
    }
}
}
