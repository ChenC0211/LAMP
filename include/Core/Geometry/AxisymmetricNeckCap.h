#pragma once
#include <array>
#include <optional>

namespace Quantum::Geometry
{
// Local terminal model, not a pinch-time detector. Two equal axial cells in
// 0<=r<=radialWidth retain their SEPARATE revolved liquid volumes. The radial
// graph in reconstruct() is r^2=radialWidth^2*q(z/axialWidth), represented by
// two C1 cubics. reconstructMatched() replaces the first piece by a bicubic
// scalar bridge that retains the complete supplied source trace and C.
// It matches the supplied left radius/slope, closes halfway through the last
// cell, and has an even scalar extension at the right reflection section.
// Positive q must remain strictly inside the first radial cell. No volume is
// clipped, transferred to another cell, or inferred from a reference shape.
class AxisymmetricNeckCap
{
public:
    struct Input {
        double radialWidth=0,axialWidth=0;
        double leftRadiusSquared=0,leftRadiusSquaredDerivative=0;
        std::array<long double,2> liquidVolumes{};
    };
    struct AxialJet {double value=0,derivative=0;}; // q and dq/dt, t=z/h.
    // Normalized source phi and h_z*dphi/dz at r=0,h_r. The other two
    // quantities are h_r*dphi/dr and h_r*h_z*d2phi/drdz at r=h_r. Axis
    // radial/mixed derivatives are exactly zero. All data are scalar data,
    // not a reconstructed radius or a signed-distance assumption.
    struct LeftTrace {AxialJet axis,outer;double outerRadial=0,outerRadialAxial=0;};
    struct ScalarJet {double value=0,dr=0,dz=0,drr=0,drz=0,dzz=0;};
    struct Report {
        std::array<long double,2> representedVolumes{},volumeResiduals{},volumeTolerances{};
        double tipLocalZ=0,tipLowerZ=0,tipUpperZ=0;
        double maxRadiusUpperBound=0,axisCurvature=0;
        double scalarSecondDerivativeJump=0; // At r=0; matched bridge may vary with r. C1, not C2.
        int matchingIterations=0;
        long double matchingQuadratureError=0; // Numerical estimate, not a bound.
    };
    static constexpr double closureFraction=.5;
    // Incompatible inventories/jets, extra roots, radial-wall contact,
    // degenerate tips and unresolved arithmetic throw. Inputs are untouched.
    static AxisymmetricNeckCap reconstruct(const Input&);
    // Preserve the complete existing cubic radial trace and its axial
    // derivative at the left seam. Last cell keeps the squared-radius cap;
    // first-cell volume determines one middle-section radius by a monotone,
    // bracketed solve. Strict radial monotonicity is required in the bridge.
    // An unmatched source trace is never silently replaced by a distance jet.
    static AxisymmetricNeckCap reconstructMatched(const Input&,const LeftTrace&);
    // Interior-neck building block: prescribe phi(r,2h)=(r/hr)^2+gasValue
    // and dphi/dz=0, instead of prescribing the half-cell tip location.
    // The two cell inventories determine the middle radius and slope; the
    // tip must remain simple and interior. Oppositely oriented caps with
    // the SAME positive scalar scale and gasValue then share a complete C1
    // gas trace. This alone does not create topology, IDs or a physical event.
    static AxisymmetricNeckCap reconstructMatchedGap(const Input&,const LeftTrace&,double gasValue);
    const Input& input() const noexcept {return input_;}
    const Report& report() const noexcept {return report_;}
    const std::array<AxialJet,3>& axialJets() const noexcept {return jets_;}
    bool matchesScalarTrace() const noexcept {return leftTrace_.has_value();}
    // phi=(r/radialWidth)^2-q. Query is restricted to this two-cell band.
    // Both sides of the internal seam share exactly the same stored C1 jet;
    // dzz at that seam is the right-sided value.
    ScalarJet scalar(double radius,double localZ) const;
private:
    AxisymmetricNeckCap()=default;
    Input input_;
    std::array<AxialJet,3> jets_{};
    std::optional<LeftTrace> leftTrace_;
    bool prescribedGap_=false;
    Report report_;
    void completeReport();
};
}
