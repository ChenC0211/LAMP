#include "Modules/Hydrodynamics/HydroFittedCapillary.h"
#include "Modules/Hydrodynamics/HydroFittedDomain.h"
#include "Modules/Hydrodynamics/HydroFittedInterfaceLoadGeometry.h"
#include "Modules/Hydrodynamics/HydroBalancedForce.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Geometry/LevelSetPhaseFraction.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[HydroFittedCapillary] ") + message);
}

template<class Query>
auto checkedQuery(const char* kind, int index, Vector2D point, Query query) -> decltype(query())
{
    try { return query(); }
    catch (const std::exception& error)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "[HydroFittedCapillary] " << kind << " " << index
            << " at (" << point.x << ',' << point.y << "): " << error.what();
        throw std::runtime_error(message.str());
    }
}
}

HydroFittedCapillaryState::HydroFittedCapillaryState(const QuadtreeMesh& mesh)
    : forceX(mesh, 0.0), forceY(mesh, 0.0) {}

void assembleFittedCapillaryForce(const QuadtreeMesh& mesh, const HydroFittedDomain& domain,
    const Quantum::Hydrodynamics::HydroFittedInterfaceLoadGeometry& geometry,
    double sigma, HydroFittedCapillaryState& state)
{
    state.source.reset();
    require(std::isfinite(sigma) && sigma >= 0.0, "nonnegative finite surface tension required");
    geometry.validate(mesh, domain, geometry.traceGeometry());
    state.halfWidth = geometry.halfWidth; state.sigma = sigma;
    state.indicator = geometry.indicator;
    state.faces = geometry.faces();
    state.maxGradientIdentityError = geometry.maxGradientIdentityError;
    state.endpointProjections = geometry.endpointProjections;
    const auto report = HydroBalancedForce::buildCapillaryFaceForce(mesh,
        geometry.faces(), sigma, state.forceX, state.forceY);
    state.maxAbsComponent = report.maxAbsComponent;
    state.integralX = report.forceX; state.integralY = report.forceY;
}

HydroFittedCapillaryState buildNativeFittedCapillaryForce(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> snapshot,
    const LevelSetInterfacePanels& panels, double halfWidth, double sigma)
{
    require(std::isfinite(sigma) && sigma >= 0.0, "nonnegative finite surface tension required");
    const auto geometry = Quantum::Hydrodynamics::buildNativeFittedInterfaceLoadGeometry(
        mesh, domain, std::move(snapshot), panels, halfWidth,
        Quantum::Hydrodynamics::FittedTransverseGradientClosure::ActiveHermiteWallGradient);
    HydroFittedCapillaryState state(mesh);
    assembleFittedCapillaryForce(mesh, domain, geometry, sigma, state);
    return state;
}

HydroFittedCapillaryState buildFittedCapillaryForce(const QuadtreeMesh& mesh,
    const HydroFittedDomain& domain, const LevelSetHermiteGeometry::Snapshot& snapshot,
    const std::vector<Vector2D>& boundaryContacts, double halfWidth, double sigma)
{
    require(std::isfinite(halfWidth) && halfWidth > 0.0 && std::isfinite(sigma) && sigma >= 0.0,
        "positive finite half-width and nonnegative finite surface tension required");
    require(domain.axisymmetric() && snapshot.radialAxisSymmetry() && snapshot.restrictedDomain(),
        "fitted capillary port currently requires an explicit active RZ scalar");
    domain.validateMeshIdentity(mesh);
    const UniformLeafGrid grid(mesh, CartesianGridScope::UniformLeaves);
    require(snapshot.nx() == grid.nx() && snapshot.ny() == grid.ny()
        && snapshot.minCoord().x == grid.minCoord().x && snapshot.minCoord().y == grid.minCoord().y
        && snapshot.maxCoord().x == grid.maxCoord().x && snapshot.maxCoord().y == grid.maxCoord().y
        && snapshot.cellSize().x == grid.spacing().x && snapshot.cellSize().y == grid.spacing().y,
        "snapshot and fitted mesh layouts disagree");
    std::vector<std::array<int, 2>> logicalCell(mesh.cellIndexNum());
    for (int j = 0; j < grid.ny(); ++j) for (int i = 0; i < grid.nx(); ++i)
    {
        const int k = mesh.cells(grid.cells()[static_cast<std::size_t>(j) * grid.nx() + i]).index();
        logicalCell[k] = {{i, j}};
        require(snapshot.cellActive(i, j) == (domain.cellRowByIndex()[k] >= 0),
            "scalar and fitted pressure domains must have identical activity");
    }
    // Reconstruct through the same public panel kernel so a stale or unrelated
    // caller-owned panel object cannot silently select another zero curve.
    const LevelSetInterfacePanels chain(mesh, snapshot, true, boundaryContacts, CartesianGridScope::UniformLeaves);
    require(!chain.closed() && !chain.panels().empty(), "fitted capillarity requires the declared open interface");
    std::vector<LevelSetHermiteGeometry::PatchEndpoint> endpoints;
    for (int end = 0; end < 2; ++end)
    {
        const auto& panel = end ? chain.panels().back() : chain.panels().front();
        const auto ij = logicalCell[panel.cellIndex];
        endpoints.push_back({end ? panel.end : panel.begin, ij[0], ij[1], end ? -1 : 1});
    }
    HydroFittedCapillaryState state(mesh);
    state.halfWidth = halfWidth; state.sigma = sigma;
    state.indicator.assign(mesh.cellIndexNum(), std::numeric_limits<double>::quiet_NaN());
    Eigen::VectorXd h(domain.cellIndexByRow().size());
    for (std::size_t row = 0; row < domain.cellIndexByRow().size(); ++row)
    {
        const int k = domain.cellIndexByRow()[row];
        const auto point = mesh.cells(mesh.cellIndexToLabel[k]).centerCoord();
        h[static_cast<Eigen::Index>(row)] = state.indicator[k]
            = LevelSetPhaseFraction::regularizedPositiveIndicator(checkedQuery("cell indicator", k, point,
                [&] { return snapshot.value(point); }), halfWidth);
    }
    const Eigen::VectorXd pressureLoad = domain.pressureGradientLoad() * h;
    require(pressureLoad.allFinite(), "non-finite fitted C(H)");
    double identityScale = 0.0;
    for (const auto& native : domain.faces())
    {
        if (native.velocityRow < 0) continue; // Exact solid/axis constraints.
        const int low = native.negativeCell, high = native.positiveCell;
        if (low < 0 || high < 0)
        {
            const int owner = low >= 0 ? low : high;
            require(owner >= 0 && (native.boundary.type == HydroFittedDomain::BoundaryType::PressureReservoir
                || native.boundary.type == HydroFittedDomain::BoundaryType::PressureOutlet),
                "unsupported free boundary force row");
            const auto center = mesh.cells(mesh.cellIndexToLabel[owner]).centerCoord();
            require(checkedQuery("pressure port owner", owner, center,
                [&] { return snapshot.closedCellOutsideBand(center, halfWidth); })
                && (state.indicator[owner] == 0.0 || state.indicator[owner] == 1.0),
                "capillary band reaches a pressure port; boundary force/phase closure is required");
            continue; // No arbitrary capillary pressure is imposed at ports.
        }
        const double integratedGradient = pressureLoad[native.velocityRow];
        const double independent = native.area * (state.indicator[high] - state.indicator[low]);
        state.maxGradientIdentityError = std::max(state.maxGradientIdentityError, std::abs(integratedGradient - independent));
        identityScale = std::max({identityScale, std::abs(integratedGradient), std::abs(independent)});
        if (integratedGradient == 0.0) continue;
        HydroInterfaceLoadGeometry::Face face;
        face.horizontal = native.horizontal;
        face.edgeIndex = native.horizontal ? native.nativeRow - mesh.vEdgeIndexNum() : native.nativeRow;
        face.lowCellIndex = low; face.highCellIndex = high; face.point = native.center;
        face.gradient = integratedGradient / native.dualVolume;
        face.dualMeasure = native.dualVolume;
        face.dualArea = native.normalSpacing * native.tangentSpacing;
        require(std::isfinite(face.gradient) && face.gradient != 0.0, "unrepresentable fitted pressure gradient density");
        const auto projection = checkedQuery("native face projection", native.nativeRow, face.point,
            [&] { return snapshot.patchBoundaryProjectedFrame(face.point, endpoints); });
        face.frame = projection.frame;
        face.boundaryEndpointProjection = projection.endpoint;
        face.endpointKktResidual = projection.kktResidual;
        state.endpointProjections += projection.endpoint ? 1 : 0;
        require(face.frame.point.x >= 0.0 && std::isfinite(face.frame.curvature), "invalid RZ projected capillary frame");
        const bool axis = face.frame.point.x == 0.0;
        require(!axis || face.frame.normal.x == 0.0, "axis curvature requires the actual even-scalar pole");
        face.surfaceCurvature = face.frame.curvature + (axis ? face.frame.curvature : face.frame.normal.x / face.frame.point.x);
        require(std::isfinite(face.surfaceCurvature), "non-finite RZ surface curvature");
        // No reconstructed kappa_s query: normal capillary force does not
        // consume D3, and restricted-domain cross-cell D3 is unsupported.
        state.faces.push_back(face);
    }
    require(state.maxGradientIdentityError <= 512.0 * std::numeric_limits<double>::epsilon() * identityScale,
        "actual C(H) disagrees with independently assembled native pressure-face load");
    const auto report = HydroBalancedForce::buildCapillaryFaceForce(mesh, state.faces, sigma, state.forceX, state.forceY);
    state.maxAbsComponent = report.maxAbsComponent;
    state.integralX = report.forceX; state.integralY = report.forceY;
    return state;
}
