#pragma once

#include "Modules/Hydrodynamics/CchPressureLoad.h"
#include "Modules/Hydrodynamics/CchRegularizedColor.h"
#include <memory>
#include <array>

namespace Quantum::Hydrodynamics
{
struct CchCsfCurvatureConfig
{
    std::size_t maxCandidateQueriesPerFace=262144;
    std::size_t maxVisibilityCellQueriesPerFace=262144;
    double relativeGradientRoundoff=256*std::numeric_limits<double>::epsilon(); // Audit only; never filters S.
};
enum class CchCsfSupportStatus {Projected,ZeroGradient,PureBoundaryReaction,SigmaZero};
struct CchCsfFaceSupportReport
{
    Vector2D center{0,0};double radius=0;
    int pressureComponent=-1;
    std::size_t candidateQueries=0,supportCount=0,visibilityRejected=0,visibilityCellQueries=0;
    long double totalWeight=0;
    double curvatureMin=0,curvatureMax=0,maxSupportDistance=0;
    double gradientRoundoffScale=0,colorGradientErrorEstimate=0;
    bool curvatureDefined=false,boundary=false,nearBoundary=false,nearAxis=false,crossesAmr=false,subResolutionGradient=false;
    CchCsfSupportStatus status=CchCsfSupportStatus::ZeroGradient;
};
struct CchCsfCurvatureReport
{
    std::size_t sourceNodes=0,selectedNodes=0,candidateQueries=0,visibilityRejected=0;
    std::size_t projectedFaces=0,pureBoundaryRows=0,nearBoundaryFaces=0,nearAxisFaces=0,amrFaces=0;
    bool projectionSkippedSigmaZero=false;
};
// Read-only balanced CSF load on the retained authoritative C=-D_I^T.
// The output is the ACTUAL momentum RHS of A*u+C*p=rhs: add it directly.
// Native boundary rows are retained; the existing saddle handles wall reactions.
class CchBalancedCsfLoad
{
public:
    static CchBalancedCsfLoad constantCurvature(const QuadtreeMesh&,CchPressureLoad,
        CchRegularizedColor,double surfaceTension,double curvature);
    static CchBalancedCsfLoad capillary(const QuadtreeMesh&,CchPressureLoad,
        CchRegularizedColor,double surfaceTension,CchCsfCurvatureConfig={});
    // Pure positive-average arithmetic, also used by capillary(). This is NOT
    // a geometry/source override port; capillary reads only its publication.
    static double positiveCurvatureAverage(const std::vector<double>& weights,const std::vector<double>& curvatures);
    // Exact closed-segment/rectangle-union predicate. The factory supplies
    // only its actual selected leaves; this is not a support override port.
    static bool selectedSegmentVisible(Vector2D begin,Vector2D end,
        const std::vector<std::array<Vector2D,2>>& selectedBoxes);
    const CchPressureLoad& pressureLoad() const noexcept;
    const CchRegularizedColor& color() const noexcept;
    const CchDofLayout& layout() const noexcept;
    double surfaceTension() const noexcept;
    double constantCurvature() const; // rejects projected loads, including skipped sigma=0
    bool isCapillary() const noexcept;
    const CchCsfCurvatureConfig& curvatureConfig() const;
    const std::vector<double>& curvatureAverages() const; // Only defined where support.curvatureDefined is true.
    const std::vector<CchCsfFaceSupportReport>& supportReports() const;
    const CchCsfCurvatureReport& curvatureReport() const;
    const std::vector<double>& colorGradient() const noexcept; // C*c, integral native face load
    const std::vector<double>& momentumRhsLoad() const noexcept; // +sigma*kappa*C*c
    Interface::HermiteGalsIdentity geometryIdentity() const;
    Interface::HermiteGalsTimeLayer timeLayer() const;
    void validateFor(const QuadtreeMesh&,const CchPressureLoad& actualPressure,
        Interface::HermiteGalsIdentity actualGeometry,Interface::HermiteGalsTimeLayer actualClock) const;
    // Also requires the SAME immutable color owner, including its epsilon.
    void validateFor(const QuadtreeMesh&,const CchPressureLoad&,const CchRegularizedColor& actualColor) const;
private:
    struct Data;
    explicit CchBalancedCsfLoad(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
}
