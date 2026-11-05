#pragma once
#include "Core/Mesh/DyadicGridCoordinates.h"
#include <array>
#include <map>
#include <vector>

class QuadtreeMesh;
// Sparse index of the actual leaves by their existing dyadic ancestry.
// Missing same-level cells return -1: no virtual children, restriction,
// prolongation, material values or solid extension are manufactured here.
// Like UniformLeafGrid, owns its maps and must be rebuilt after adaptation.
class CartesianLeafLevelIndex
{
public:
    struct Address { int i=0,j=0; unsigned level=0; };
    explicit CartesianLeafLevelIndex(const QuadtreeMesh&);
    const Address& address(int nativeCell) const { return addresses_.at(nativeCell); }
    int cell(unsigned level,int i,int j) const;
    std::array<int,2> extent(unsigned level) const;
    Vector2D spacing(unsigned level) const { return coordinates_.cellSize(level); }
private:
    LAMP::Mesh::DyadicGridCoordinates coordinates_;
    int nx_,ny_;
    std::vector<Address> addresses_;
    std::map<std::array<int,3>,int> cells_;
};
