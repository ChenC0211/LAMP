#pragma once

#include "Modules/Interface/HermiteGalsState.h"
#include <cstddef>
#include <vector>

namespace Quantum::Interface
{
struct HermiteGalsMixedClosureReport
{
    std::size_t vertices=0;
    // max |Dx(phi_y)-Dy(phi_x)| / max(1,|Dx(phi_y)|,|Dy(phi_x)|).
    // Diagnostic only; does not select weights or change acceptance.
    double maxRelativeDirectionalDisagreement=0;
};
class HermiteGalsMixedClosureResult;

// Complete uniform leaf skeleton, >=2 cells in each direction. Candidate must
// have no published geometry. Reads only this primary revision's first jet;
// no old closure, scalar reconstruction, velocity or Hessian transport.
HermiteGalsMixedClosureResult closeHermiteGalsMixedUniform(
    const QuadtreeMesh& mesh,const HermiteGalsCandidate& candidate);

class HermiteGalsMixedClosureResult
{
public:
    // Physical phi_xy, in candidate.primary().vertices order.
    const std::vector<double>& vertexMixedDerivatives() const noexcept { return mixed_; }
    HermiteGalsIdentity primaryIdentity() const noexcept { return primary_; }
    HermiteGalsTimeLayer timeLayer() const noexcept { return time_; }
    const HermiteGalsMixedClosureReport& report() const noexcept { return report_; }
private:
    HermiteGalsMixedClosureResult(HermiteGalsIdentity primary,HermiteGalsTimeLayer time,
        std::vector<double> mixed,HermiteGalsMixedClosureReport report);
    HermiteGalsIdentity primary_;
    HermiteGalsTimeLayer time_;
    std::vector<double> mixed_;
    HermiteGalsMixedClosureReport report_;
    friend HermiteGalsMixedClosureResult closeHermiteGalsMixedUniform(const QuadtreeMesh&,const HermiteGalsCandidate&);
};
} // namespace Quantum::Interface
