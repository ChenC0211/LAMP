#pragma once

#include "Core/Base/Vector2D.h"
#include "Core/Geometry/UnitPolynomialRoots.h"

class LevelSetMappedSegment;

// Retained normal-flux trace of ONE actual native face, including a 2:1 child
// face. Owns its chart and polynomial; no uniform-grid/owner-cell addressing.
// The polynomial is the existing fitted Hydro Hermite face primitive, not a
// new momentum discretization. Positive normal is +x/+y or +r/+z. XY uses
// unit depth; RZ includes 2*pi*r. No extra metric is applied by consumers.
class HydroFittedFaceFlux
{
  public:
    struct Chart
    {
        Vector2D begin, midpoint, end;
    };

    // q is the full native integrated normal flux. dBegin/dEnd are the
    // existing streamfunction tangential jets; on an axis horizontal face
    // they are derivatives in xi^2, otherwise in the whole-face parameter.
    HydroFittedFaceFlux(Chart chart, bool horizontal, bool radialPolynomial, long double q,
                        long double dBegin, long double dEnd, bool axisymmetric = true);

    // Prepare the SAME face polynomial from physical endpoint velocity traces
    // and its independent FV face flux. Endpoints are supplied once/shared by
    // the caller; this does not infer wall/axis/contact laws or construct a
    // volume velocity field. Actual half-edge metrics are retained. The
    // default preserves fitted RZ; axisymmetric=false selects XY, including
    // negative x and nonzero normal transport through x=0.
    static HydroFittedFaceFlux fromEndpointVelocities(Chart chart, bool horizontal, long double q,
                                                      Vector2D beginVelocity, Vector2D endVelocity,
                                                      bool axisymmetric = true);

    const Chart &chart() const noexcept
    {
        return chart_;
    }
    long double nativeFlux() const noexcept
    {
        return q;
    }

    // half 0: begin -> midpoint; half 1: midpoint -> end.
    // t=anchor0+(anchor1-anchor0)*s. Reversal changes the integral sign.
    // Keep interval factors separate even when composed endpoints round equal.
    long double integrate(int half, double anchor0, double anchor1, double sBegin, double sEnd) const;
    long double integrate(int half, double begin, double end) const
    {
        return integrate(half, 0.0, 1.0, begin, end);
    }

    // Segment anchors must exactly equal the corresponding positive half.
    // Integrates only negative-phase ranges, including endpoint-relative
    // sub-ULP ranges. This is floating integration, not an interval enclosure.
    long double integrate(int half, const LevelSetMappedSegment &segment) const;

    // Normal velocity corresponding to this same integral polynomial, on a
    // retained half-face. Includes its physical RZ metric; the horizontal-axis
    // limit is evaluated analytically. This does not supply tangential speed.
    double normalVelocity(int half, double parameter) const;

    // P'(t) of the SAME half-face flux primitive, including its RZ metric.
    // The near-end offset stays separate; fromEnd preserves the +r/+z sign.
    // Unlike normalVelocity, this is flux per unit half-face parameter.
    long double density(int half, bool fromEnd, long double offset) const;
    // Derivative in the positive half-face parameter, including when evaluated
    // through the reversed endpoint chart. Uses the same retained primitive.
    long double densityDerivative(int half, bool fromEnd, long double offset) const;
    // Roots of the SAME retained primitive's density polynomial, in the
    // positive half-face chart. Includes endpoint roots and their enclosures;
    // a zero trace is explicit. Exact basis combinations preserve endpoint
    // cancellations; unresolved multiplicity throws in Core.
    Quantum::Geometry::UnitPolynomialRoots densityRoots(int half) const;
    std::vector<Quantum::Geometry::UnitPolynomialTerm> densityPolynomial(int half) const;

  private:
    Chart chart_;
    bool radialPolynomial;
    bool axisymmetric_;
    long double q, d0, d1, sign;
};
