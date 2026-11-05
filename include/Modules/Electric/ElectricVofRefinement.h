#pragma once
#include "Modules/Electric/ElectricVofField.h"
#include <cstddef>
#include <vector>

namespace LAMP::Electric
{
// Passive, same-End observation of the independent signed cell charge.
// Entries are native-cell indexed; excluded-domain entries remain zero.
// eta_i=max_f h_normal,i*|rho_j-rho_i|/|centroid_normal,j-centroid_normal,i|,
// rho=Q/actual cell volume. RZ uses full 2*pi*r volumes and volume centroids.
// The indicator has charge-density units. No floor, normalization, field solve
// or geometry iteration is applied; boundary faces imply no exterior charge.
class ElectricVofRefinementIndicators
{
public:
    const std::vector<long double>& densityJump() const noexcept { return densityJump_; }
    long double maximumDensityJump() const noexcept { return maximumDensityJump_; }
    std::size_t interiorFaces() const noexcept { return interiorFaces_; }
    void validateFor(const QuadtreeMesh&,const Interface::GeometricVofPublication&,
        const ElectricVofCharge&) const;
private:
    friend class ElectricVofRefinement;
    ElectricVofRefinementIndicators(const Interface::GeometricVofPublication&,
        const ElectricVofCharge&,std::vector<long double>,long double,std::size_t);
    Interface::GeometricVofPublication geometry_;
    ElectricVofCharge charge_;
    std::vector<long double> densityJump_;
    long double maximumDensityJump_=0;
    std::size_t interiorFaces_=0;
};
class ElectricVofRefinement
{
public:
    // Requires the genuine geometry/charge owners at one complete End clock.
    // One cell pass and one actual native-subface pass; O(N+F) work/storage.
    static ElectricVofRefinementIndicators sample(const QuadtreeMesh&,
        const Interface::GeometricVofPublication&,const ElectricVofCharge&);
};
}
