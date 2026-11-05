#include "Modules/Interface/HermiteGalsBubble.h"

#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Modules/Interface/Detail/HermiteGalsVolumeRoundoff.h"
#include "Core/Geometry/ActiveQuarticHermitePatch.h"
#include "Core/Geometry/LevelSetHermiteGeometry.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Interface
{
namespace
{
using Patch = LAMP::Geometry::ActiveQuarticHermitePatch;
using Snapshot = LevelSetHermiteGeometry::Snapshot;
using Interval = Snapshot::QuarticCenterInterval;
void require(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
bool tolerance(double value) { return std::isfinite(value) && value > 0 && value < 1; }
struct Evaluation
{
    std::shared_ptr<const Snapshot> geometry;
    InterfaceCellMeasure::DiscoveredBubbleVolumes measured;
    std::vector<long double> volume; // compact candidate order, certified pure values
    long double maxMomentEstimate = 0;
};
struct Root
{
    bool active = false;
    Interval certificate;
    bool expanded = false;
    double symmetricLower = 0, symmetricUpper = 0;
    double lo = 0, hi = 0, previousStep = 0;
    long double vlo = 0, vhi = 0;
};
}
struct HermiteGalsBubbleResult::Data
{
    explicit Data(const HermiteGalsAmrClosureResult& input) : mixed(input) {}
    HermiteGalsAmrClosureResult mixed;
    HermiteGalsClosure closure;
    std::vector<HermiteGalsBubbleCellReport> cells;
    HermiteGalsBubbleReport report;
};
HermiteGalsBubbleResult::HermiteGalsBubbleResult(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
const HermiteGalsBubbleResult::Data& HermiteGalsBubbleResult::checked() const
{
    if (!data_) throw std::invalid_argument("[Interface/HermiteGalsBubble] moved-from result");
    return *data_;
}
HermiteGalsIdentity HermiteGalsBubbleResult::primaryIdentity() const { return checked().mixed.primaryIdentity(); }
HermiteGalsTimeLayer HermiteGalsBubbleResult::timeLayer() const { return checked().mixed.timeLayer(); }
const HermiteGalsClosure& HermiteGalsBubbleResult::closure() const { return checked().closure; }
const std::vector<HermiteGalsBubbleCellReport>& HermiteGalsBubbleResult::cells() const { return checked().cells; }
const HermiteGalsBubbleReport& HermiteGalsBubbleResult::report() const { return checked().report; }
void HermiteGalsBubbleResult::validateCandidate(const QuadtreeMesh& mesh, const HermiteGalsCandidate& candidate) const
{
    const auto& data = checked();
    require(candidate.identity() == data.mixed.primaryIdentity(), "[Interface/HermiteGalsBubble] result revision mismatch");
    candidate.validateMesh(mesh);
    data.mixed.skeleton().validateMesh(mesh);
    require(!candidate.hasGeometry(), "[Interface/HermiteGalsBubble] candidate is already published");
}

HermiteGalsBubbleResult closeHermiteGalsBubble(const QuadtreeMesh& mesh,
    const HermiteGalsCandidate& candidate, const HermiteGalsAmrClosureResult& mixed,
    const HermiteGalsBubblePrediction& prediction, const HermiteGalsBubbleConfig& config)
{
    HermiteGalsIdentity id;
    int diagnosticCell = -1;
    const char* stage = "input";
    try {
        id = candidate.identity(); candidate.validateMesh(mesh);
        require(!candidate.hasGeometry(), "requires an unpublished candidate");
        require(mixed.primaryIdentity() == id && prediction.primaryIdentity == id, "primary revision mismatch");
        mixed.skeleton().validateMesh(mesh);
        const auto time = candidate.timeLayer(), mt = mixed.timeLayer();
        require(time.step == mt.step && time.begin == mt.begin && time.end == mt.end && time.position == mt.position,
            "mixed closure time mismatch");
        require(tolerance(config.relativeVolumeTolerance) && tolerance(config.relativeMomentEstimateTolerance)
            && tolerance(config.relativeLedgerTolerance) && tolerance(config.retainedDerivativeFraction)
            && std::isfinite(config.maxDisplacementCells) && config.maxDisplacementCells > 0
            && config.maxIterations > 0, "invalid bubble configuration");
        require(config.volumeConstraint == HermiteGalsBubbleVolumeConstraint::InventoryRoundoff
            || config.volumeConstraint == HermiteGalsBubbleVolumeConstraint::GeometryAccuracy,"invalid bubble volume constraint");
        const bool accuracy=config.volumeConstraint == HermiteGalsBubbleVolumeConstraint::GeometryAccuracy;
        // The explicit legacy contract solves below its roundoff-only flux
        // consumer's budget. GeometryAccuracy instead honors the requested
        // geometry tolerance; its full driver needs bounded shared fluxes.
        const long double workingVolumeTolerance=accuracy?config.relativeVolumeTolerance:std::min<long double>(
            config.relativeVolumeTolerance,Detail::hermiteInventoryRoundoffBudget/4);
        // A certified pure scalar stays unchanged. GeometryAccuracy retains
        // a bounded inventory discrepancy using the same geometry budget;
        // its consumer must use the bounded inventory transport representation.
        // The legacy consumer still requires its original roundoff budget.
        const long double pureVolumeTolerance=accuracy?config.relativeVolumeTolerance:std::min<long double>(
            config.relativeVolumeTolerance,Detail::hermiteInventoryRoundoffBudget);
        const long double acceptedVolumeTolerance=accuracy?config.relativeVolumeTolerance:pureVolumeTolerance;
        const auto& p = candidate.primary(); const auto& cells = p.geometryCells;
        require(!cells.empty() && cells.size() == p.phaseVolumes.size(), "inventory requires complete scalar support");
        require(prediction.coefficients.size() == cells.size(), "prediction size mismatch");
        require(accuracy?config.accuracyBudgets.size()==cells.size():config.accuracyBudgets.empty(),
            "accuracy budget layout does not match volume constraint");
        require(mixed.vertexMixedDerivatives().size() == p.vertices.size(), "mixed closure size mismatch");
        Patch patch;
        patch.activeCells.assign(mesh.cellIndexNum(), 0);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        patch.nodes.assign(mesh.vertIndexNum(), {nan,nan,nan,nan});
        patch.cellCenters.assign(mesh.cellIndexNum(), nan);
        for (std::size_t k = 0; k < cells.size(); ++k) {
            diagnosticCell = cells[k];
            require(p.phaseVolumes[k].cell == cells[k], "inventory/scalar support mismatch");
            require(std::isfinite(prediction.coefficients[k]), "nonfinite prediction");
            patch.activeCells[cells[k]] = 1; patch.cellCenters[cells[k]] = prediction.coefficients[k];
        }
        for (std::size_t k = 0; k < p.vertices.size(); ++k) {
            const auto& v = p.vertices[k];
            patch.nodes[v.vertex] = {v.phi,v.dx,v.dy,mixed.vertexMixedDerivatives()[k]};
        }
        for (const auto& e : p.edges) patch.edgeMidpoints.push_back({e.vertices,e.phi,e.dn});
        const bool radial = p.measure == HermiteGalsMeasure::AxisymmetricVolume;
        FaceSegmentGraph graph;
        graph.build(mesh, [&](const QuadtreeMesh&, const Cell& cell) { return patch.activeCells[cell.index()] != 0; });
        const auto domain = NativeLeafDomain::build(mesh, graph, radial);
        const auto& capacity = domain.cellVolumes();
        require(domain.cellIndexByRow() == cells, "native domain order mismatch");
        long double totalCapacity = 0;
        for (std::size_t k = 0; k < cells.size(); ++k) {
            diagnosticCell = cells[k]; const auto target = p.phaseVolumes[k].liquid;
            require(std::isfinite(target) && target >= 0 && target <= capacity[k], "inventory outside strict capacity bounds");
            totalCapacity += capacity[k];
        }
        require(std::isfinite(totalCapacity) && totalCapacity > 0, "unrepresentable total capacity");
        std::size_t evaluations = 0, newtonSteps = 0, bisectionSteps = 0, endpointSteps = 0;
        long double maxEstimate = 0;
        const auto evaluate = [&](const Patch& trial) {
            Evaluation e;
            e.geometry = std::make_shared<const Snapshot>(Snapshot::fromNativeLeaves(mesh, trial, radial));
            e.measured = InterfaceCellMeasure::discoverBubbleVolumes(mesh, domain, e.geometry,
                config.boundaryContacts, config.symmetryPlanes);
            ++evaluations;
            const auto& measure = e.measured.measure;
            require(std::abs(measure.fixedBoundary.assemblyClosure) <= config.relativeLedgerTolerance * totalCapacity,
                "global boundary ledger defect");
            for (std::size_t k = 0; k < cells.size(); ++k) {
                const int c = cells[k]; diagnosticCell = c;
                const auto raw = measure.negativeVolumeByCell[c];
                const auto error = measure.momentErrorEstimateByCell[c] / capacity[k];
                require(std::isfinite(raw) && std::isfinite(error) && error >= 0
                    && error <= config.relativeMomentEstimateTolerance, "cell volume quadrature estimate exceeds tolerance");
                e.maxMomentEstimate = std::max(e.maxMomentEstimate, error);
                require(std::abs(measure.fixedBoundary.capacityByCell[c] - capacity[k])
                    <= config.relativeLedgerTolerance * capacity[k], "cell boundary capacity defect");
                const int phase = e.measured.cellPhase[c];
                long double volume = raw;
                if (phase) {
                    require(phase == -1 || phase == 1, "invalid pure certificate");
                    volume = phase < 0 ? capacity[k] : 0;
                    require(std::abs(volume - raw) <= config.relativeLedgerTolerance * capacity[k], "pure certificate/raw volume mismatch");
                } else require(raw >= 0 && raw <= capacity[k], "cut volume outside strict capacity bounds");
                const auto derivative = e.measured.centerDerivative[c];
                require(std::isfinite(derivative) && derivative <= 0, "bubble volume derivative is positive/nonfinite");
                e.volume.push_back(volume);
            }
            maxEstimate = std::max(maxEstimate, e.maxMomentEstimate);
            return e;
        };
        stage = "predictor measure";
        auto current = evaluate(patch);
        const auto initialVolume = current.volume;
        const auto initialPhase = current.measured.cellPhase;
        std::vector<double> hessianFactors(cells.size(),0),accuracyLower=prediction.coefficients,accuracyUpper=prediction.coefficients;
        if(accuracy) for(std::size_t k=0;k<cells.size();++k){
            const int c=cells[k];diagnosticCell=c;const auto& budget=config.accuracyBudgets[k];
            require(budget.cell==c && std::isfinite(budget.maximumHessianComponentChange)
                && budget.maximumHessianComponentChange>=0,"invalid cell Hessian accuracy budget");
            const auto bounds=current.geometry->nativeCellBounds(c);const auto h=bounds[1]-bounds[0];
            // B=f(xi)f(eta), f=16*t^2*(1-t)^2: |f|<=1,
            // |f'|<=4, |f''|<=32 on [0,1]. Physical component bound.
            const double factor=std::max({32/(h.x*h.x),32/(h.y*h.y),16/(h.x*h.y)});
            const double radius=budget.maximumHessianComponentChange/factor,predicted=prediction.coefficients[k];
            require(std::isfinite(factor) && factor>0 && std::isfinite(radius),"unrepresentable Hessian accuracy interval");
            hessianFactors[k]=factor;
            double lo=predicted-radius,hi=predicted+radius;
            require(std::isfinite(lo) && std::isfinite(hi),"unrepresentable Hessian accuracy interval");
            // Inward rounding of the stored endpoints keeps the requested
            // derivative budget, including cancellation near large centers.
            if(predicted-lo>radius || (predicted-lo)*factor>budget.maximumHessianComponentChange)lo=std::nextafter(lo,predicted);
            if(hi-predicted>radius || (hi-predicted)*factor>budget.maximumHessianComponentChange)hi=std::nextafter(hi,predicted);
            require(lo<=predicted && predicted<=hi && (predicted-lo)*factor<=budget.maximumHessianComponentChange
                && (hi-predicted)*factor<=budget.maximumHessianComponentChange,"invalid Hessian accuracy interval");
            accuracyLower[k]=lo;accuracyUpper[k]=hi;
        }
        const auto restrictAccuracy=[&](std::size_t k,Interval interval){
            if(accuracy){interval.lower=std::max(interval.lower,accuracyLower[k]);interval.upper=std::min(interval.upper,accuracyUpper[k]);}
            return interval;
        };
        const auto components = current.measured.components;
        const auto checkTopology = [&](const Evaluation& e) {
            require(e.measured.components == components, "component topology changed inside certified interval");
            for (int c : cells) require(e.measured.cellPhase[c] == initialPhase[c], "cell phase topology changed inside certified interval");
        };
        std::vector<Root> roots(cells.size());
        auto lowerPatch = patch, upperPatch = patch;
        std::size_t active = 0;
        stage = "topology/amplitude guard";
        for (std::size_t k = 0; k < cells.size(); ++k) {
            const int c = cells[k]; diagnosticCell = c; auto& root = roots[k];
            const auto error = std::abs(current.volume[k] - p.phaseVolumes[k].liquid);
            if (initialPhase[c]) { require(error <= pureVolumeTolerance * capacity[k],
                accuracy?"pure-cell mismatch exceeds geometry volume budget":"pure-cell mismatch would create an unmarked component"); continue; }
            // Keep predictors already within the selected geometry budget.
            // Chasing smaller residuals can amplify tiny dV/dc into large
            // coefficient changes, without improving geometric accuracy.
            if (error <= acceptedVolumeTolerance * capacity[k]) continue;
            if(accuracy)require(accuracyLower[k]<prediction.coefficients[k] && prediction.coefficients[k]<accuracyUpper[k],
                "no representable nontrivial Hessian accuracy interval");
            root.certificate = restrictAccuracy(k,current.geometry->nativeQuarticCenterInterval(c, config.maxDisplacementCells, config.retainedDerivativeFraction));
            root.lo = root.certificate.lower; root.hi = root.certificate.upper;
            root.symmetricLower = root.lo; root.symmetricUpper = root.hi;
            root.previousStep = root.hi - root.lo;
            require(std::isfinite(root.previousStep) && root.lo < patch.cellCenters[c] && patch.cellCenters[c] < root.hi,
                "no representable nontrivial certified interval");
            root.active = true; ++active;
            lowerPatch.cellCenters[c] = root.lo; upperPatch.cellCenters[c] = root.hi;
        }
        std::vector<long double> minimum = initialVolume, maximum = initialVolume;
        if (active) {
            stage = "volume bracket";
            auto lower = evaluate(lowerPatch), upper = evaluate(upperPatch);
            checkTopology(lower); checkTopology(upper);
            bool expanded = false;
            for (std::size_t k = 0; k < cells.size(); ++k) if (roots[k].active) {
                const int c = cells[k]; diagnosticCell = c; auto& root = roots[k];
                const long double tol = workingVolumeTolerance * capacity[k], target = p.phaseVolumes[k].liquid;
                require(upper.volume[k] <= initialVolume[k] + tol && initialVolume[k] <= lower.volume[k] + tol,
                    "volume bracket violated monotonicity");
                const bool decrease = target > lower.volume[k] + tol, increase = target < upper.volume[k] - tol;
                if (!decrease && !increase) continue;
                // M(c) is nonincreasing. Prove the entire requested one-sided
                // interval from the ORIGINAL predictor; never accumulate
                // several displacement budgets or recenter the certificate.
                const auto direction = decrease ? Snapshot::QuarticCenterDirection::Decrease : Snapshot::QuarticCenterDirection::Increase;
                const auto directed = restrictAccuracy(k,current.geometry->nativeQuarticCenterInterval(c,
                    config.maxDisplacementCells, config.retainedDerivativeFraction, direction));
                const auto predicted = patch.cellCenters[c];
                require(directed.lower <= predicted && predicted <= directed.upper && directed.lower < directed.upper,
                    "invalid directed center certificate");
                require(decrease ? directed.lower < root.lo : directed.upper > root.hi,
                    accuracy?"target outside Hessian accuracy/geometry volume budgets":"target outside certified attainable volume interval: directed search cannot extend bracket");
                root.certificate = directed; root.lo = directed.lower; root.hi = directed.upper;
                root.previousStep = root.hi - root.lo;
                require(std::isfinite(root.previousStep) && root.previousStep > 0, "unrepresentable directed interval");
                root.expanded = true; expanded = true;
                lowerPatch.cellCenters[c] = root.lo; upperPatch.cellCenters[c] = root.hi;
            }
            if (expanded) {
                lower = evaluate(lowerPatch); upper = evaluate(upperPatch);
                checkTopology(lower); checkTopology(upper);
            }
            for (std::size_t k = 0; k < cells.size(); ++k) if (roots[k].active) {
                diagnosticCell = cells[k]; auto& root = roots[k];
                root.vlo = lower.volume[k]; root.vhi = upper.volume[k];
                minimum[k] = root.vhi; maximum[k] = root.vlo;
                const long double tol = workingVolumeTolerance * capacity[k], target = p.phaseVolumes[k].liquid;
                require(root.vhi <= initialVolume[k] + tol && initialVolume[k] <= root.vlo + tol,
                    "volume bracket violated monotonicity");
                require(target >= root.vhi - tol && target <= root.vlo + tol, "target outside certified attainable volume interval");
            }
        }
        stage = "safeguarded root";
        int iterations = 0;
        for (;;) {
            bool changed = false;
            for (std::size_t k = 0; k < cells.size(); ++k) if (roots[k].active) {
                const int c = cells[k]; diagnosticCell = c; auto& root = roots[k];
                const long double target = p.phaseVolumes[k].liquid, error = current.volume[k] - target;
                const long double tol = workingVolumeTolerance * capacity[k];
                require(current.volume[k] >= root.vhi - tol && current.volume[k] <= root.vlo + tol,
                    "root evaluation violated retained volume bracket");
                if (std::abs(error) <= tol) continue;
                require(iterations < config.maxIterations, "scalar root iteration budget exhausted");
                const double center = patch.cellCenters[c];
                if (error > 0) { root.lo = center; root.vlo = current.volume[k]; }
                else { root.hi = center; root.vhi = current.volume[k]; }
                double next;
                if (std::abs(root.vlo - target) <= tol) { next = root.lo; ++endpointSteps; }
                else if (std::abs(root.vhi - target) <= tol) { next = root.hi; ++endpointSteps; }
                else {
                    const auto derivative = current.measured.centerDerivative[c];
                    next = derivative < 0 ? static_cast<double>(center - error / derivative) : nan;
                    if (!std::isfinite(next) || next <= root.lo || next >= root.hi
                        || std::abs(next - center) > .5 * root.previousStep) {
                        next = .5 * root.lo + .5 * root.hi; ++bisectionSteps;
                    } else ++newtonSteps;
                }
                require(std::isfinite(next) && next != center && next >= root.lo && next <= root.hi,
                    "scalar root stagnated at storage resolution");
                root.previousStep = std::abs(next - center);
                patch.cellCenters[c] = next; changed = true;
            }
            if (!changed) break;
            current = evaluate(patch); checkTopology(current); ++iterations;
        }
        stage = "final closure";
        // All decisions use local solve data above, never a Report field.
        for (std::size_t k = 0; k < cells.size(); ++k) {
            diagnosticCell = cells[k];
            const auto finalTolerance=initialPhase[cells[k]]?pureVolumeTolerance:(roots[k].active?workingVolumeTolerance:acceptedVolumeTolerance);
            require(std::abs(current.volume[k] - p.phaseVolumes[k].liquid) <= finalTolerance * capacity[k],
                "final cell volume residual exceeds tolerance");
            if(accuracy)require(std::abs(patch.cellCenters[cells[k]]-prediction.coefficients[k])*hessianFactors[k]
                <=config.accuracyBudgets[k].maximumHessianComponentChange,"final Hessian correction exceeds accuracy budget");
        }
        auto result = std::make_shared<HermiteGalsBubbleResult::Data>(mixed);
        result->closure.vertexMixedDerivatives = mixed.vertexMixedDerivatives();
        auto& report = result->report;
        report.relativeWorkingVolumeTolerance = workingVolumeTolerance;
        report.relativeAcceptedVolumeTolerance = acceptedVolumeTolerance;
        report.relativePureVolumeTolerance = pureVolumeTolerance;
        report.volumeConstraint = config.volumeConstraint;
        report.correctedCells = active; report.geometryEvaluations = evaluations; report.iterations = iterations;
        report.newtonSteps = newtonSteps; report.bisectionSteps = bisectionSteps; report.endpointSteps = endpointSteps;
        report.maxRelativeMomentEstimate = maxEstimate;
        for (std::size_t k = 0; k < cells.size(); ++k) {
            const int c = cells[k]; const auto& root = roots[k]; const auto& cert = root.certificate;
            const double coefficient = patch.cellCenters[c], predicted = prediction.coefficients[k];
            result->closure.cellBubbleCoefficients.push_back(coefficient);
            const long double target = p.phaseVolumes[k].liquid, residual = current.volume[k] - target;
            result->cells.push_back({c,initialPhase[c],root.active,predicted,coefficient,
                root.active ? cert.lower : predicted,root.active ? cert.upper : predicted,
                cert.partialAxis,cert.partialSign,cert.minimumPhysicalDerivative,cert.maximumDisplacementPerCenter,
                initialVolume[k],current.volume[k],current.measured.measure.negativeVolumeByCell[c],target,
                minimum[k],maximum[k],current.measured.centerDerivative[c],root.expanded,
                root.active ? root.symmetricLower : predicted,root.active ? root.symmetricUpper : predicted,
                accuracy?config.accuracyBudgets[k].maximumHessianComponentChange:0,
                std::abs(coefficient-predicted)*hessianFactors[k]});
            report.maxHessianComponentChangeBound=std::max(report.maxHessianComponentChangeBound,
                std::abs(coefficient-predicted)*hessianFactors[k]);
            report.directionallyExpandedCells += root.expanded;
            report.initialMaxRelativeResidual = std::max(report.initialMaxRelativeResidual, std::abs(initialVolume[k] - target) / capacity[k]);
            report.finalMaxRelativeResidual = std::max(report.finalMaxRelativeResidual, std::abs(residual) / capacity[k]);
            if (initialPhase[c]) report.maxPureRelativeResidual = std::max(report.maxPureRelativeResidual, std::abs(residual) / capacity[k]);
            report.signedResidual += residual; report.l1Residual += std::abs(residual);
            report.maxCoefficientCorrection = std::max(report.maxCoefficientCorrection, std::abs(coefficient - predicted));
            if (root.active) {
                const auto bounds = current.geometry->nativeCellBounds(c); const auto h = bounds[1] - bounds[0];
                report.maxDisplacementBoundCells = std::max(report.maxDisplacementBoundCells,
                    std::abs(coefficient - predicted) * cert.maximumDisplacementPerCenter / std::min(h.x,h.y));
            }
        }
        return HermiteGalsBubbleResult(std::move(result));
    } catch (const std::exception& e) {
        throw std::invalid_argument("[Interface/HermiteGalsBubble primary=" + std::to_string(id.value())
            + " cell=" + std::to_string(diagnosticCell) + " stage=" + stage + "] " + e.what());
    }
}
}
