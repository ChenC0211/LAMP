#include "Modules/Hydrodynamics/Detail/HydroNativeVelocityTrace.h"
#include "Core/Base/Interpolation.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace HydroFittedDetail
{
namespace {
void require(bool ok,const char* why) { if(!ok)throw std::invalid_argument(why); }
const Cell& cellAt(const QuadtreeMesh& mesh,int index) { return mesh.cells(mesh.cellIndexToLabel[index]); }
int nativeAt(const QuadtreeMesh& mesh,int component,int label)
{return label<0?-1:component==0?mesh.vEdges(label).index():mesh.vEdgeIndexNum()+mesh.hEdges(label).index();}
}
std::vector<int> nativeTraceSupport(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const std::vector<int>& seed)
{
    std::vector<int> cells;
    for(int index:seed)if(index>=0 && domain.cellRowByIndex()[index]>=0)cells.push_back(index);
    std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    const auto original=cells;
    for(int index:original)for(const auto& side:cellAt(mesh,index).neiCell)for(int label:side)
        if(label>=0 && domain.cellRowByIndex()[mesh.cells(label).index()]>=0)cells.push_back(mesh.cells(label).index());
    std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    return cells;
}
NativeVelocityTrace nativeVelocityTrace(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::vector<int>& cells,int component,Vector2D point)
{
    std::vector<int> rows;Vector2D scale{0,0};
    for(int index:cells) {
        const auto& cell=cellAt(mesh,index);
        scale.x=std::max(scale.x,cell.size().x);scale.y=std::max(scale.y,cell.size().y);
        for(int side=0;side<2;++side)for(int part=0;part<2;++part) {
            const int label=component==0?cell.cellToVEdge[side][part]:cell.cellToHEdge[side][part];
            const int native=nativeAt(mesh,component,label);if(native>=0)rows.push_back(native);
        }
    }
    std::sort(rows.begin(),rows.end());rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    require(rows.size()>=3 && scale.x>0 && scale.y>0,"native linear trace lacks active support");
    std::vector<Vector2D> points;
    for(int row:rows) {
        const auto face=domain.face(mesh,row);
        require(face.negativeCell>=0 || face.positiveCell>=0,"linear trace references solid-only storage");
        points.push_back({(face.center.x-point.x)/scale.x,(face.center.y-point.y)/scale.y});
    }
    long double rank=0;
    for(std::size_t i=1;i<points.size();++i)for(std::size_t j=i+1;j<points.size();++j)
        rank=std::max(rank,std::abs((static_cast<long double>(points[i].x)-points[0].x)*(points[j].y-points[0].y)
            -(static_cast<long double>(points[j].x)-points[0].x)*(points[i].y-points[0].y)));
    require(rank>0,"singular native linear trace");
    const auto value=Interpolation::leastSquareInterpolationCoeff(points,{0,0});
    const auto dx=Interpolation::leastSquareInterpolationCoeff(points,{1,0});
    const auto dy=Interpolation::leastSquareInterpolationCoeff(points,{0,1});
    NativeVelocityTrace result;std::array<std::array<long double,3>,3> moments{};
    for(std::size_t i=0;i<rows.size();++i) {
        const std::array<double,3> weight{value[i],dx[i]-value[i],dy[i]-value[i]};
        for(int k=0;k<3;++k) {
            require(std::isfinite(weight[k]),"nonfinite native linear trace coefficient");
            moments[k][0]+=weight[k];moments[k][1]+=weight[k]*points[i].x;moments[k][2]+=weight[k]*points[i].y;
        }
        result.value.emplace_back(rows[i],value[i]);
        result.gradient[0].emplace_back(rows[i],weight[1]/scale.x);
        result.gradient[1].emplace_back(rows[i],weight[2]/scale.y);
    }
    for(int k=0;k<3;++k)for(int p=0;p<3;++p)
        require(std::abs(moments[k][p]-(k==p?1.L:0.L))<1e-10L,"native linear trace lost constant/affine consistency");
    return result;
}
}
