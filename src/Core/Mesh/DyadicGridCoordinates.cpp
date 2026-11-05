#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/Mesh/BackgroundMesh.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace LAMP::Mesh
{
namespace
{
constexpr std::uint64_t exactIntegerLimit=std::uint64_t{1}<<53;
void require(bool condition,const char* reason)
{
    if(!condition) throw std::runtime_error(std::string("[DyadicGridCoordinates] ")+reason);
}
std::uint64_t countAtLevel(std::uint64_t count,unsigned level)
{
    require(level<=53 && count<= (exactIntegerLimit>>level),
        "refinement address exceeds exactly representable integer range");
    return count<<level;
}
double refinedSpacing(double base,unsigned level)
{
    // countAtLevel has already bounded level. A positive subnormal result
    // can still have lost bits while halving: multiplying its address by
    // two would then select a different point on the next-level chart.
    const double spacing=std::ldexp(base,-static_cast<int>(level));
    require(std::isfinite(spacing) && spacing>0.0,"unrepresentable refined spacing");
    require(std::ldexp(spacing,static_cast<int>(level))==base,
        "refined spacing loses bits through dyadic underflow");
    return spacing;
}
}

DyadicGridCoordinates::DyadicGridCoordinates(Vector2D minimum,Vector2D maximum,int nx,int ny)
    : minimum_(minimum),maximum_(maximum),nx_(nx),ny_(ny)
{
    require(nx>0 && ny>0 && std::isfinite(minimum.x) && std::isfinite(minimum.y)
        && std::isfinite(maximum.x) && std::isfinite(maximum.y)
        && maximum.x>minimum.x && maximum.y>minimum.y,"invalid Cartesian domain");
    spacing_={(maximum.x-minimum.x)/nx,(maximum.y-minimum.y)/ny};
    require(std::isfinite(spacing_.x) && std::isfinite(spacing_.y)
        && spacing_.x>0.0 && spacing_.y>0.0,"unrepresentable base spacing");
}

double DyadicGridCoordinates::axis(std::uint64_t index,unsigned level,int direction) const
{
    const auto count=countAtLevel(direction==0?nx_:ny_,level);
    require(index<=count,"address is outside the Cartesian domain");
    const double lower=direction==0?minimum_.x:minimum_.y;
    const double upper=direction==0?maximum_.x:maximum_.y;
    const double spacing=refinedSpacing(direction==0?spacing_.x:spacing_.y,level);
    const auto value=[&](std::uint64_t k) {
        if(k==0) return lower;
        if(k==count) return upper;
        const double offset=static_cast<double>(k)*spacing;
        return lower+offset;
    };
    const double coordinate=value(index);
    require(std::isfinite(coordinate) && coordinate>=lower && coordinate<=upper,
        "computed coordinate left the finite domain");
    require((index==0 || value(index-1)<coordinate) && (index==count || coordinate<value(index+1)),
        "adjacent logical coordinates collapse in binary64");
    return coordinate;
}

Vector2D DyadicGridCoordinates::point(std::uint64_t i,std::uint64_t j,unsigned level) const
{
    return {axis(i,level,0),axis(j,level,1)};
}

Vector2D DyadicGridCoordinates::cellSize(unsigned level) const
{
    countAtLevel(nx_,level); countAtLevel(ny_,level);
    return {refinedSpacing(spacing_.x,level),refinedSpacing(spacing_.y,level)};
}

DyadicCellAddress cellAddress(const BackgroundMesh& mesh,const Cell& cell)
{
    const DyadicGridCoordinates coordinates(mesh.minCoord,mesh.maxCoord,mesh.Nx,mesh.Ny);
    const Cell* current=&cell;
    std::vector<unsigned> path;
    while(current->parent>=0)
    {
        require(path.size()<53 && current->parent<mesh.cellLabelCapacity(),"invalid or cyclic parent ancestry");
        const auto& parent=mesh.cells(current->parent);
        require(parent.level()+1==current->level(),"parent and child levels disagree");
        unsigned found=4;
        for(unsigned k=0;k<4;++k) if(parent.son[k]==current->label())
        {
            require(found==4,"child appears twice in its parent"); found=k;
        }
        require(found<4,"parent does not reference its child");
        path.push_back(found); current=&parent;
    }
    require(current->level()==0,"ancestry does not terminate at a base cell");
    const auto center=current->centerCoord();
    const auto spacing=coordinates.cellSize(0);
    const double i=std::round((center.x-mesh.minCoord.x)/spacing.x-0.5);
    const double j=std::round((center.y-mesh.minCoord.y)/spacing.y-0.5);
    require(std::isfinite(i) && std::isfinite(j) && i>=0 && j>=0 && i<mesh.Nx && j<mesh.Ny,
        "base cell has no valid Cartesian address");
    DyadicCellAddress result{static_cast<std::uint64_t>(i),static_cast<std::uint64_t>(j),0};
    const auto expected=coordinates.point(2*result.i+1,2*result.j+1,1);
    require(center.x==expected.x && center.y==expected.y,
        "base cell uses legacy recursive coordinates; rebuild the mesh from case geometry");
    constexpr unsigned dx[4]={0,1,1,0},dy[4]={0,0,1,1};
    for(auto entry=path.rbegin();entry!=path.rend();++entry)
    {
        require(result.i<=exactIntegerLimit/2 && result.j<=exactIntegerLimit/2,"cell address overflow");
        result.i=2*result.i+dx[*entry]; result.j=2*result.j+dy[*entry]; ++result.level;
    }
    coordinates.cellSize(result.level);
    return result;
}
}
