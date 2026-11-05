#pragma once
#include <vector>
class QuadtreeMesh;
namespace LAMP::Hydrodynamics
{
class HydroTwoPhaseState;
// Actual accepted regular reconstruction samples in native-cell order.
// omega=du_y/dx-du_x/dy, including the same r,z component convention.
// h=min(actual dx,dy); h|omega| has code_velocity units. This diagnostic
// neither redefines LocalPLIC transport nor uses its curl-free RT0 interior.
struct HydroVofRefinementIndicators
{
    std::vector<double> vorticity,hAbsVorticity;
    double maximumHAbsVorticity=0;
};
class HydroVofRefinement
{
public:
    static HydroVofRefinementIndicators sample(const QuadtreeMesh&,const HydroTwoPhaseState&);
};
}
