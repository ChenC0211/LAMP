#pragma once
#include "Modules/Interface/GeometricVofPublication.h"
#include <functional>

namespace LAMP::Interface
{
using GeometricVofInitialValueQuery=std::function<double(Vector2D)>;
struct GeometricVofInitializationConfig
{
    unsigned subcellsPerAxis=8; // Strictly 1..32, once at t=0 only.
    GeometricVofConfig geometry;
};
// Finite-sampling piecewise-linear INITIAL inventory approximation on actual
// native rectangles. Each fixed subrectangle has a four-triangle centre fan.
// Both stocks are independently measured by Core cuts of phi and -phi; no
// high-order field or global interface topology is constructed.
// Smooth resolved regular zero sets normally have second-order geometric
// accuracy in the subcell spacing. This is not an error certificate: thin
// phases/closed islands between query sites may be missed. These stocks are
// NOT promised to equal any historical quadratic/Hermite initialization.
// The query is consumed here, never retained or resampled during transport.
class GeometricVofInitialization
{
public:
    // Whole-leaf selected fluid domains are supplied by Core. Only selected
    // leaves are queried; excluded solid leaves carry no invented gas stock.
    static GeometricVofPublication prepare(const QuadtreeMesh&,NativeLeafDomain,
        const GeometricVofInitialValueQuery&,const GeometricVofInitializationConfig& = {});
    // Complete-domain convenience consumer of the same initializer above.
    static GeometricVofPublication prepare(const QuadtreeMesh&,bool axisymmetric,
        const GeometricVofInitialValueQuery&,const GeometricVofInitializationConfig& = {});
};
}
