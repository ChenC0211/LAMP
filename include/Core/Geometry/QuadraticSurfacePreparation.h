#pragma once

#include "Core/Geometry/QuadraticSurfaceMeasure.h"

namespace LAMP::Geometry::QuadraticSurfaceMeasure
{
using MaterialCurves=std::vector<std::shared_ptr<const Detail::QuadraticParametricData>>;
// Same immutable physical graph and ALL real clipping arcs, prepared by Core.
// A later physical measure may change XY/RZ or coincident-side ownership,
// while retaining exactly these arc events. No caller-supplied area is accepted.
class SurfaceArcPreparation
{
public:
    static SurfaceArcPreparation prepare(const Rectangle& source,const Graph&,const Rectangle& region,bool axisymmetric);
    const Result& surface() const noexcept;
    Result evaluate(const Rectangle& source,const Graph&,const Rectangle& region,
        bool axisymmetric,unsigned includedCoincidentSides=All) const;
private:
    struct Data;
    explicit SurfaceArcPreparation(std::shared_ptr<const Data>);
    std::shared_ptr<const Data> data_;
};
struct SurfacePreparationTiming
{
    double totalSeconds=0,eventSeconds=0,areaIntegralSeconds=0;
    double derivativeSeconds=0,materialCurveSeconds=0;
};
struct PreparedSurface
{
    Sensitivity sensitivity;
    SurfacePreparationTiming timing; // Passive elapsed, with nested phases; not CPU or a certificate.
    MaterialCurves materialCurves; // Same actual patch order and original exact side-root identities.
};
// Compute the original physical surface (and optionally its original moving-
// boundary sensitivity), retaining the SAME clipping events for material curves.
// No mesh/state/phase owner, across-query cache, or accepted-stock substitution.
// fixedCornerEndpoints selects only the existing constrained derivative chart;
// it is ignored when withDerivative is false. Legacy inputs retain their bridge.
PreparedSurface prepareWithMaterialCurves(const Rectangle& source,const Graph&,const Rectangle& region,
    bool axisymmetric,bool withDerivative,bool fixedCornerEndpoints=false,unsigned includedCoincidentSides=All);
// Original physical operation, with optional curve retention and passive phases.
// False retention performs no extra material-curve construction.
PreparedSurface prepareWithTiming(const Rectangle& source,const Graph&,const Rectangle& region,
    bool axisymmetric,bool withDerivative,bool fixedCornerEndpoints=false,
    bool retainMaterialCurves=false,unsigned includedCoincidentSides=All);
}
