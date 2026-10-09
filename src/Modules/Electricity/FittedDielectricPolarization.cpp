#include "Modules/Electricity/FittedDielectricPolarization.h"
#include "Modules/Electricity/DielectricPolarizationForce.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
void require(bool condition, const char* reason)
{
    if (!condition) throw std::invalid_argument(std::string("[FittedDielectricPolarization] ") + reason);
}
}

FittedDielectricPolarization::FittedDielectricPolarization(const QuadtreeMesh& mesh)
    : fieldX(mesh.cellIndexNum(), 0.0), fieldY(mesh.cellIndexNum(), 0.0) {}

FittedDielectricPolarization buildFittedDielectricPolarization(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain, const FaceSegmentGraph& graph,
    const FittedPotentialInput& potential, const std::vector<double>& epsilonCell,
    const std::vector<double>& epsilonSegment,
    const std::vector<PolarizationPoissonFVSolverCellCentered::FittedBoundary>& boundary,
    double relativeTolerance)
{
    using BC = PolarizationPoissonFVSolverCellCentered::FittedBoundaryType;
    domain.validateMeshIdentity(mesh);
    const int count = mesh.cellIndexNum(), nv = mesh.vEdgeIndexNum();
    const auto& segments = graph.segments();
    const auto& offsets = graph.cellOffsets();
    require(potential.size() == count && epsilonCell.size() == static_cast<std::size_t>(count)
        && epsilonSegment.size() == segments.size() && boundary.size() == segments.size(), "electric input dimensions disagree");
    require(std::isfinite(relativeTolerance) && relativeTolerance > 0.0 && relativeTolerance < 1.0,
        "positive finite relative tolerance required");
    require(offsets.size() == static_cast<std::size_t>(count) + 1 && offsets.front() == 0
        && offsets.back() == segments.size(), "invalid electric graph offsets");
    for (int k = 0; k < count; ++k)
    {
        require(offsets[k] <= offsets[k+1], "unordered electric cell offsets");
        const int row = domain.cellRowByIndex()[k];
        const auto faceCount = offsets[k+1] - offsets[k];
        // A 2:1 coarse leaf has split sides: count real segments, not four
        // logical sides. Exact owner/edge identity and complete oriented
        // incidence are still checked below before returning any result.
        require(row < 0 ? faceCount == 0 : faceCount >= 4 && faceCount <= 8,
            "electric/domain active cell graphs disagree");
        if (row < 0) continue; // No solid potential/epsilon enters a field reconstruction.
        require(potential.label(k) == mesh.cellIndexToLabel[k] && std::isfinite(potential[k])
            && std::isfinite(epsilonCell[k]) && epsilonCell[k] > 0.0,
            "invalid active electric material or potential");
    }
    FittedDielectricPolarization result(mesh);
    result.geometry = domain;
    result.permittivityByCell = epsilonCell;
    result.tangentByNative.assign(domain.nativeFaceCount(), 0.0);
    result.coefficientByNative.assign(domain.nativeFaceCount(), 0.0);
    std::vector<double> dx(count, 0.0), dy(count, 0.0), wx(count, 0.0), wy(count, 0.0);
    std::vector<double> residual(count, 0.0), rhs(count, 0.0), dn(domain.faces(mesh).size(), 0.0);
    std::vector<unsigned> incidence(domain.faces(mesh).size(), 0);
    double boundaryFlux = 0.0, absoluteBoundaryFlux = 0.0;
    for (std::size_t s = 0; s < segments.size(); ++s)
    {
        const auto& segment = segments[s];
        require(segment.K >= 0 && segment.K < mesh.cellLabelCapacity(), "invalid electric owner label");
        const int k = mesh.cells(segment.K).index();
        require(k >= 0 && k < count && domain.cellRowByIndex()[k] >= 0
            && mesh.cellIndexToLabel[k] == segment.K && s >= offsets[k] && s < offsets[k+1], "electric owner is not the actual active cell");
        const auto n = faceSideToNormal(segment.side);
        const bool horizontal = segment.side == FaceSide::South || segment.side == FaceSide::North;
        require(segment.isHorizontal == horizontal && segment.n.x == n.x && segment.n.y == n.y
            && (std::abs(n.x) + std::abs(n.y) == 1.0), "invalid oriented electric face normal");
        require(segment.edgeLabel >= 0 && segment.edgeLabel < (horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity()),
            "invalid electric edge label");
        const auto& edge = horizontal ? mesh.hEdges(segment.edgeLabel) : mesh.vEdges(segment.edgeLabel);
        const int native = edge.index() + (horizontal ? nv : 0);
        require(edge.index() >= 0 && native < static_cast<int>(domain.faces(mesh).size()), "electric edge is not active native storage");
        const auto& face = domain.faces(mesh)[native];
        const bool positive = n.x + n.y > 0.0;
        require(segment.K == (positive ? edge.slaveCell : edge.masterCell)
            && segment.L == (positive ? edge.masterCell : edge.slaveCell)
            && face.edgeLabel == segment.edgeLabel && face.center.x == segment.faceCenter.x
            && face.center.y == segment.faceCenter.y && edge.length() == segment.Af,
            "electric graph is stale or does not match fitted native geometry");
        const unsigned bit = positive ? 1u : 2u;
        require((incidence[native] & bit) == 0, "duplicate oriented electric face");
        const int other = positive ? face.positiveCell : face.negativeCell;
        require((positive ? face.negativeCell : face.positiveCell) == k, "electric owner/domain incidence disagreement");
        const bool interior = other >= 0;
        const auto& bc = boundary[s];
        require((interior && bc.type == BC::Interior) || (!interior && (bc.type == BC::Dirichlet || bc.type == BC::HomogeneousNeumann)),
            "electric BC type disagrees with actual fitted fluid face");
        require(std::isfinite(bc.value) && (bc.type == BC::Dirichlet || bc.value == 0.0), "invalid electric BC value");
        require(face.area != 0.0 || (domain.axisymmetric() && !horizontal && face.center.x == 0.0
            && mesh.minCoord.x == 0.0 && segment.side == FaceSide::West && segment.L == -1
            && bc.type == BC::HomogeneousNeumann && bc.value == 0.0),
            "zero-area axis requires the actual homogeneous-Neumann symmetry condition");
        const auto center = mesh.cells(segment.K).centerCoord();
        const double ownerDistance = std::abs(horizontal ? face.center.y - center.y : face.center.x - center.x);
        const auto opposite = interior ? mesh.cells(mesh.cellIndexToLabel[other]).centerCoord() : face.center;
        const double distance = interior ? std::abs(horizontal ? opposite.y - center.y : opposite.x - center.x) : ownerDistance;
        require(std::isfinite(distance) && distance > 0.0 && std::isfinite(ownerDistance) && ownerDistance > 0.0,
            "invalid actual electric half/full link");
        require(!interior || segment.df == distance, "electric interior link differs from solved FV link");
        const double expectedEpsilon = interior
            ? (ownerDistance + std::abs(horizontal ? opposite.y - face.center.y : opposite.x - face.center.x))
                / (ownerDistance / epsilonCell[k] + std::abs(horizontal ? opposite.y - face.center.y : opposite.x - face.center.x) / epsilonCell[other])
            : epsilonCell[k];
        require(std::isfinite(epsilonSegment[s]) && epsilonSegment[s] > 0.0 && epsilonSegment[s] == expectedEpsilon,
            "force must consume the same cell-harmonic coefficients as the electric solve");
        const double drop = interior ? potential[k] - potential[other]
            : bc.type == BC::Dirichlet ? potential[k] - bc.value : 0.0;
        const double displacement = epsilonSegment[s] * drop / distance;
        const double flux = face.area * displacement;
        require(std::isfinite(displacement) && std::isfinite(flux), "non-finite solved electric displacement flux");
        const double cartesianD = (positive ? 1.0 : -1.0) * displacement;
        require(incidence[native] == 0 || dn[native] == cartesianD, "opposite electric fluxes disagree on one native face");
        dn[native] = cartesianD; incidence[native] |= bit;
        // Actual RZ face area already appears here, including exact axis A=0.
        // Never divide a zero-area axis flux by its area.
        dx[k] += n.x * flux; dy[k] += n.y * flux;
        wx[k] += face.area * n.x * n.x; wy[k] += face.area * n.y * n.y;
        residual[k] += flux;
        if (!interior)
        {
            boundaryFlux += flux; absoluteBoundaryFlux += std::abs(flux);
            if (bc.type == BC::Dirichlet) rhs[k] += epsilonSegment[s] * face.area * bc.value / distance;
        }
    }
    double residualNorm = 0.0, rhsNorm = 0.0;
    for (int k : domain.cellIndexByRow())
    {
        require(std::isfinite(wx[k]) && wx[k] > 0.0 && std::isfinite(wy[k]) && wy[k] > 0.0,
            "invalid active electric reconstruction measures");
        result.fieldX[k] = (dx[k] / wx[k]) / epsilonCell[k];
        result.fieldY[k] = (dy[k] / wy[k]) / epsilonCell[k];
        require(std::isfinite(result.fieldX[k]) && std::isfinite(result.fieldY[k])
            && std::isfinite(residual[k]) && std::isfinite(rhs[k]), "non-finite reconstructed electric field or Gauss residual");
        residualNorm = std::hypot(residualNorm, residual[k]); rhsNorm = std::hypot(rhsNorm, rhs[k]);
    }
    result.chargeFreeResidual = rhsNorm > 0.0 ? residualNorm / rhsNorm : residualNorm;
    result.boundaryFluxRelativeClosure = absoluteBoundaryFlux > 0.0 ? std::abs(boundaryFlux) / absoluteBoundaryFlux : std::abs(boundaryFlux);
    require(std::isfinite(result.chargeFreeResidual) && result.chargeFreeResidual <= relativeTolerance,
        "actual FV flux divergence is not a solved charge-free electric state");
    require(std::isfinite(result.boundaryFluxRelativeClosure) && result.boundaryFluxRelativeClosure
        <= std::max(relativeTolerance, 128.0 * std::numeric_limits<double>::epsilon() * domain.cellIndexByRow().size()),
        "actual charge-free boundary displacement fails closure");
    for (const auto& face : domain.faces(mesh))
    {
        require(incidence[face.nativeRow] == ((face.negativeCell >= 0 ? 1u : 0u) | (face.positiveCell >= 0 ? 2u : 0u)),
            "missing actual oriented electric faces");
        if (face.negativeCell < 0 || face.positiveCell < 0) continue;
        const int low = face.negativeCell, high = face.positiveCell;
        if (epsilonCell[low] == epsilonCell[high]) continue;
        const double tangent = face.horizontal ? .5 * result.fieldX[low] + .5 * result.fieldX[high]
            : .5 * result.fieldY[low] + .5 * result.fieldY[high];
        result.tangentByNative[face.nativeRow] = tangent;
        result.coefficientByNative[face.nativeRow] = DielectricPolarizationForce::gradientCoefficient(
            dn[face.nativeRow], tangent, epsilonCell[low], epsilonCell[high]);
        require(std::isfinite(result.coefficientByNative[face.nativeRow]), "non-finite polarization coefficient");
    }
    result.displacementByNative = std::move(dn);
    return result;
}
