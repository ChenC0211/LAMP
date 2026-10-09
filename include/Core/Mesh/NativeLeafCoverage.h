#pragma once

#include "Core/Mesh/CartesianLeafLocator.h"

class NativeLeafDomain;

// Closed union of selected native rectangles. Geometry only: no extrapolated
// field, ghost gas, contact angle or boundary condition. Borrowed inputs must
// outlive this view and remain unchanged. Segment coverage uses actual cell
// bounds and exact rational intervals, including zero-length corner contacts.
class NativeLeafCoverage
{
public:
    NativeLeafCoverage(const QuadtreeMesh& mesh, const NativeLeafDomain& domain);
    const CartesianLeafLocator& locator() const noexcept { return locator_; }
    bool covers(Vector2D a, Vector2D b) const;
    bool contains(Vector2D p) const { return covers(p,p); }
private:
    const QuadtreeMesh& mesh_;
    const NativeLeafDomain& domain_;
    CartesianLeafLocator locator_;
    bool complete_;
};
