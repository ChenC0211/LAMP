#pragma once

#include "Modules/Interface/HermiteGalsPublication.h"
#include "Modules/Interface/HermiteGalsTrace.h"
#include <cstddef>
#include <functional>
#include <memory>

namespace Quantum::Interface
{
struct HermiteGalsBoundaryFirstJet { double phi = 0; Vector2D gradient; };
// Neutral consumer of an already resolved map and its SAME derivative.
// Owns its immutable callable and explicit provenance; no Hydro dependency.
class HermiteGalsResolvedMotion
{
public:
    using Query=std::function<LevelSetBoundaryTrace::JetResult(Vector2D)>;
    HermiteGalsResolvedMotion(HermiteGalsIdentity source,HermiteGalsTimeLayer sourceTime,
        HermiteGalsIdentity velocity,HermiteGalsTimeLayer targetTime,Query);
    HermiteGalsIdentity sourceIdentity() const noexcept {return source_;}
    HermiteGalsIdentity velocityIdentity() const noexcept {return velocity_;}
    HermiteGalsTimeLayer sourceTime() const noexcept {return sourceTime_;}
    HermiteGalsTimeLayer targetTime() const noexcept {return targetTime_;}
    LevelSetBoundaryTrace::JetResult trace(Vector2D) const;
private:
    HermiteGalsIdentity source_,velocity_;
    HermiteGalsTimeLayer sourceTime_,targetTime_;
    Query query_;
};
// Actual prescribed data at the OUTER arrival point/time, not a source-field
// extension or an interior reset. Captures must be immutable and repeatable.
using HermiteGalsBoundaryFirstJetQuery = std::function<HermiteGalsBoundaryFirstJet(
    Vector2D, HermiteGalsTimeLayer, Boundary)>;
enum class HermiteGalsBoundaryTreatment { InflowHistory, PrescribedOuterFirstJet };
// Passive steady-clock diagnostics. Never inputs to numerical acceptance.
// The seven phases partition preparation; total also includes timer/report overhead.
struct HermiteGalsFirstJetStepTimings
{
    double input=0,boundary=0,vertices=0,edges=0,centers=0,mixed=0,publication=0,total=0;
};
struct HermiteGalsFirstJetStepReport
{
    std::size_t transportedVertices = 0, transportedEdges = 0, tracedCenters = 0;
    std::size_t prescribedVertices = 0, prescribedEdges = 0;
    std::size_t sourceCells = 0, sourceNativeFaces = 0, sourceCutCells = 0, sourceGraphs = 0;
    HermiteGalsFirstJetStepTimings timings;
};
class HermiteGalsFirstJetStep
{
public:
    HermiteGalsIdentity sourceIdentity() const;
    HermiteGalsIdentity arrivalVelocityIdentity() const;
    const HermiteGalsPublication& publication() const;
    const HermiteGalsFirstJetStepReport& report() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsFirstJetStep(std::shared_ptr<const Data>);
    const Data& checked() const;
    friend HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh&,
        const HermiteGalsPublication&, const LevelSetBoundaryTrace::Domain&,
        const HermiteGalsVelocityView&, const HermiteGalsVelocityView&,
        const HermiteGalsBoundaryFirstJetQuery&, const HermiteGalsPublicationConfig&);
    friend HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh&,
        const HermiteGalsPublication&,const HermiteGalsResolvedMotion&,
        const HermiteGalsBoundaryFirstJetQuery&,const HermiteGalsPublicationConfig&,HermiteGalsBoundaryTreatment);
};
// Prepare one Midpoint/End candidate from an accepted Initial/End publication.
// Arrival velocity time is the explicit target clock; its geometry identity
// identifies the velocity provenance, not the newly generated publication.
// No accepted commit, regrid, inventory correction, output or fallback.
HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source, const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity, const HermiteGalsVelocityView& arrivalVelocity,
    const HermiteGalsBoundaryFirstJetQuery& boundary = {}, const HermiteGalsPublicationConfig& config = {});
HermiteGalsFirstJetStep prepareHermiteGalsFirstJetStep(const QuadtreeMesh& mesh,
    const HermiteGalsPublication& source,const HermiteGalsResolvedMotion&,
    const HermiteGalsBoundaryFirstJetQuery& boundary = {},const HermiteGalsPublicationConfig& config = {},
    HermiteGalsBoundaryTreatment treatment=HermiteGalsBoundaryTreatment::InflowHistory);
}
