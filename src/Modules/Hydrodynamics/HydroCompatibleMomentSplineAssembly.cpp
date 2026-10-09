#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineAssembly.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

namespace Quantum::Hydrodynamics {
namespace {
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
}
std::shared_ptr<HydroCompatibleMomentSpline::NativeLayout::CoordinateAssembly>
HydroCompatibleMomentSpline::NativeLayout::coordinateAssembly(const HydroCompatibleMomentComplement* complement) const{
    std::lock_guard<std::mutex> lock(mutex);
    if(!assemblyCache||bool(assemblyCache->complement)!=bool(complement)
        ||(complement && &assemblyCache->complement->basisMap()!=&complement->basisMap()))
        assemblyCache=std::make_shared<CoordinateAssembly>(*this,complement);
    return assemblyCache;
}
HydroCompatibleMomentSpline::NativeLayout::CoordinateAssembly::CoordinateAssembly(
    const NativeLayout& layout,const HydroCompatibleMomentComplement* extra){
    if(extra)complement=*extra;
    std::vector<int> nativeOffsets;auto nativeTransform=layout.localCoordinateMap(extra,nativeOffsets);
    initialize(std::move(nativeTransform),std::move(nativeOffsets));geometricGram.resize(cells.size());
}
MomentSplineCoordinateAssembly::MomentSplineCoordinateAssembly(Matrix map,std::vector<int> boundaries){
    initialize(std::move(map),std::move(boundaries));
}
void MomentSplineCoordinateAssembly::initialize(Matrix map,std::vector<int> boundaries){
    require(boundaries.size()>1&&boundaries.front()==0&&boundaries.back()==map.rows(),"invalid complete local assembly offsets");
    for(std::size_t k=1;k<boundaries.size();k++)require(boundaries[k]>boundaries[k-1],"empty or overlapping local assembly block");
    transform=std::move(map);offsets=std::move(boundaries);rowTransform=transform;
    const int size=int(transform.cols());cells.resize(offsets.size()-1);
    std::vector<std::vector<int>> active(cells.size()),columns(size);
    for(int slot=0;slot<int(cells.size());slot++){
        auto& ids=active[slot];
        for(int row=offsets[slot];row<offsets[slot+1];row++)
            for(decltype(rowTransform)::InnerIterator it(rowTransform,row);it;++it)if(it.value()!=0)ids.push_back(int(it.col()));
        std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
        std::vector<Eigen::Triplet<double>> entries;
        for(int row=offsets[slot];row<offsets[slot+1];row++)
            for(decltype(rowTransform)::InnerIterator it(rowTransform,row);it;++it)if(it.value()!=0)
                entries.emplace_back(row-offsets[slot],int(std::lower_bound(ids.begin(),ids.end(),it.col())-ids.begin()),it.value());
        auto& cell=cells[slot];cell.map.resize(offsets[slot+1]-offsets[slot],ids.size());cell.map.setFromTriplets(entries.begin(),entries.end());
        for(int col:ids)columns[col].insert(columns[col].end(),ids.begin(),ids.end());
    }
    Eigen::Index count=0;
    for(auto& rows:columns){std::sort(rows.begin(),rows.end());rows.erase(std::unique(rows.begin(),rows.end()),rows.end());count+=rows.size();}
    require(count<=std::numeric_limits<int>::max(),"complete assembly pattern exceeds sparse index range");
    pattern.resize(size,size);pattern.reserve(count);
    for(int col=0;col<size;col++){pattern.startVec(col);for(int row:columns[col])pattern.insertBackByOuterInner(col,row)=0;}
    pattern.finalize();
    for(int slot=0;slot<int(cells.size());slot++){
        const auto& ids=active[slot];auto& positions=cells[slot].positions;positions.resize(ids.size()*ids.size());
        for(int j=0;j<int(ids.size());j++){
            const int col=ids[j];const auto& rows=columns[col];
            for(int i=0;i<int(ids.size());i++)positions[i+ids.size()*j]=pattern.outerIndexPtr()[col]
                +int(std::lower_bound(rows.begin(),rows.end(),ids[i])-rows.begin());
        }
    }
}
MomentSplineCoordinateAssembly::Matrix MomentSplineCoordinateAssembly::assemble(
    const std::function<const Eigen::MatrixXd&(int)>& block,const std::vector<double>* scales) const{
    require(!scales||scales->size()==cells.size(),"complete assembly weights changed the native layout");
    std::vector<Eigen::MatrixXd> mapped(cells.size());std::vector<std::exception_ptr> failures(cells.size());
#pragma omp parallel for schedule(dynamic,1) if(cells.size()>1)
    for(int slot=0;slot<int(cells.size());slot++){try{
        require(!scales||std::isfinite((*scales)[slot]),"nonfinite complete assembly weight");
        if(scales&&(*scales)[slot]==0)continue;
        const auto& b=block(slot);const auto& map=cells[slot].map;
        require(b.rows()==map.rows()&&b.cols()==map.rows(),"complete local matrix changed its coordinate space");
        if(b.isZero(0))continue;
        Eigen::MatrixXd scaled;if(scales)scaled=(*scales)[slot]*b;
        const Eigen::MatrixXd first=map.transpose()*(scales?scaled:b);
        mapped[slot]=first*map;
        require(mapped[slot].allFinite(),"nonfinite assembled local form");
    }catch(...){failures[slot]=std::current_exception();}}
    bool nonzero=false;
    for(int slot=0;slot<int(cells.size());slot++){if(failures[slot])std::rethrow_exception(failures[slot]);nonzero|=mapped[slot].size()!=0;}
    if(!nonzero)return Matrix(pattern.rows(),pattern.cols());
    Matrix result=pattern.cast<double>();
    require(result.nonZeros()==pattern.nonZeros(),"complete assembly lost its cached structure");
    // Preserve native cell order for shared global coefficients. No atomic
    // additions or scheduling-dependent physical reductions are introduced.
    for(int slot=0;slot<int(cells.size());slot++)if(mapped[slot].size()){
        const auto& values=mapped[slot];const auto& positions=cells[slot].positions;
        for(Eigen::Index k=0;k<values.size();k++)result.valuePtr()[positions[k]]+=values.data()[k];
    }
    result.prune(0.);return result;
}
}
