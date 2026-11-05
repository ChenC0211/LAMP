#pragma once

#include "Modules/Interface/HermiteGalsAmrClosure.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include <cstddef>
#include <memory>
#include <vector>

namespace LAMP::Interface
{
struct HermiteGalsBubblePrediction
{
    HermiteGalsIdentity primaryIdentity;
    // Same order as candidate.primary().geometryCells. Temporary prediction,
    // never an extra transported state. Gauge: zero-center cardinal extension
    // phi_boundary, hence c equals the Q4 center value in THIS API only.
    std::vector<double> coefficients;
};
enum class HermiteGalsBubbleVolumeConstraint { InventoryRoundoff, GeometryAccuracy };
struct HermiteGalsBubbleAccuracyBudget
{
    int cell = -1;
    double maximumHessianComponentChange = 0; // physical xx/xy/yy; zero freezes the predicted center
};
struct HermiteGalsBubbleConfig
{
    double relativeVolumeTolerance = 1e-10; // geometric residual / capacity; legacy mode also imposes inventory roundoff
    double relativeMomentEstimateTolerance = 1e-11; // per-cell GL8/subdivision estimate / capacity
    double relativeLedgerTolerance = 1e-10;
    double maxDisplacementCells = 0.125;
    double retainedDerivativeFraction = 0.25;
    int maxIterations = 64;
    std::vector<Vector2D> boundaryContacts;
    std::vector<LevelSetInterfacePanels::SymmetryPlane> symmetryPlanes;
    // GeometryAccuracy covers cut AND certified pure geometry; pure centers
    // stay fixed and M-V is transported by the bounded inventory source.
    // A full driver must use that shared flux, not the legacy geometry flux.
    HermiteGalsBubbleVolumeConstraint volumeConstraint = HermiteGalsBubbleVolumeConstraint::InventoryRoundoff;
    std::vector<HermiteGalsBubbleAccuracyBudget> accuracyBudgets; // complete candidate geometry-cell order in GeometryAccuracy
};
struct HermiteGalsBubbleCellReport
{
    int cell = -1, phase = 0;
    bool corrected = false;
    double predictedCoefficient = 0, coefficient = 0, lowerCoefficient = 0, upperCoefficient = 0;
    int partialAxis = -1, partialSign = 0;
    double minimumPhysicalDerivative = 0, maximumDisplacementPerCoefficient = 0;
    long double initialVolume = 0, geometricVolume = 0, rawVolume = 0, targetVolume = 0;
    long double minimumVolume = 0, maximumVolume = 0, finalDerivative = 0;
    bool directionalExpansion = false;
    double symmetricLowerCoefficient = 0, symmetricUpperCoefficient = 0;
    double maximumHessianComponentChange = 0, hessianComponentChangeBound = 0;
};
struct HermiteGalsBubbleReport
{
    std::size_t correctedCells = 0, geometryEvaluations = 0, newtonSteps = 0, bisectionSteps = 0, endpointSteps = 0;
    std::size_t directionallyExpandedCells = 0;
    int iterations = 0;
    long double initialMaxRelativeResidual = 0, finalMaxRelativeResidual = 0;
    long double relativeWorkingVolumeTolerance = 0; // cut-cell scalar root
    long double relativeAcceptedVolumeTolerance = 0; // cut acceptance / root activation; pure tolerance is separate
    long double relativePureVolumeTolerance = 0; // no volume-changing pure-cell degree of freedom
    long double maxPureRelativeResidual = 0, maxRelativeMomentEstimate = 0;
    long double signedResidual = 0, l1Residual = 0;
    double maxCoefficientCorrection = 0, maxDisplacementBoundCells = 0;
    HermiteGalsBubbleVolumeConstraint volumeConstraint = HermiteGalsBubbleVolumeConstraint::InventoryRoundoff;
    double maxHessianComponentChangeBound = 0;
};
class HermiteGalsBubbleResult
{
public:
    HermiteGalsIdentity primaryIdentity() const;
    HermiteGalsTimeLayer timeLayer() const;
    const HermiteGalsClosure& closure() const;
    const std::vector<HermiteGalsBubbleCellReport>& cells() const;
    const HermiteGalsBubbleReport& report() const;
    void validateCandidate(const QuadtreeMesh& mesh, const HermiteGalsCandidate& candidate) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsBubbleResult(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsBubbleResult closeHermiteGalsBubble(const QuadtreeMesh&,
        const HermiteGalsCandidate&, const HermiteGalsAmrClosureResult&,
        const HermiteGalsBubblePrediction&, const HermiteGalsBubbleConfig&);
};
// M is read from the unchanged candidate. Only derived bubble coefficients
// are returned. Trial geometry remains private; publication belongs to I11.
HermiteGalsBubbleResult closeHermiteGalsBubble(const QuadtreeMesh& mesh,
    const HermiteGalsCandidate& candidate, const HermiteGalsAmrClosureResult& mixed,
    const HermiteGalsBubblePrediction& prediction, const HermiteGalsBubbleConfig& config = {});
}
