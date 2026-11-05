#include "Core/Mesh/NativeLeafOverlap.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>

namespace LAMP::Mesh
{
namespace
{
namespace P=Geometry::PlanarCellCut;
void require(bool ok,const std::string& message)
{if(!ok)throw std::invalid_argument("[NativeLeafOverlap] "+message);}
P::Rectangle rectangle(const QuadtreeMesh& mesh,int native)
{
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
    // The native vertices are the shared dyadic sites used by phase geometry.
    // Rebuilding a boundary from a rounded center and nominal half size can
    // move it by an ULP, creating non-nested boxes on an unchanged tree.
    return {mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
}
bool contains(const P::Rectangle& a,const P::Rectangle& b)
{return a.lower.x<=b.lower.x && a.lower.y<=b.lower.y && a.upper.x>=b.upper.x && a.upper.y>=b.upper.y;}
bool equal(const P::Rectangle& a,const P::Rectangle& b){return contains(a,b)&&contains(b,a);}
std::string nestingContext(const QuadtreeMesh& source,int si,const P::Rectangle& sourceBox,
    const QuadtreeMesh& target,int ti,const P::Rectangle& targetBox)
{
    std::ostringstream text;text<<std::setprecision(std::numeric_limits<double>::max_digits10);
    const auto leaf=[&](const char* side,const QuadtreeMesh& mesh,int native,const P::Rectangle& box)
    {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);
        const auto center=cell.centerCoord(),size=cell.size();
        text<<" ["<<side<<" native="<<native<<" label="<<cell.label()<<" parent="<<cell.parent
            <<" level="<<cell.level()<<" center=("<<center.x<<','<<center.y<<") size=("<<size.x<<','<<size.y
            <<") rect=(("<<box.lower.x<<','<<box.lower.y<<"),("<<box.upper.x<<','<<box.upper.y<<")) vertices=";
        for(int k=0;k<4;++k){const auto p=mesh.verts(cell.cellToVert[k]).coord();
            text<<'('<<cell.cellToVert[k]<<':'<<p.x<<','<<p.y<<')';}
        try {
            const auto address=cellAddress(mesh,cell);
            const DyadicGridCoordinates coordinates(mesh.minCoord,mesh.maxCoord,mesh.Nx,mesh.Ny);
            const auto lo=coordinates.point(address.i,address.j,address.level);
            const auto hi=coordinates.point(address.i+1,address.j+1,address.level);
            text<<" address=("<<address.i<<','<<address.j<<','<<address.level<<") dyadic=(("
                <<lo.x<<','<<lo.y<<"),("<<hi.x<<','<<hi.y<<"))";
        } catch(const std::exception& error){text<<" address_error="<<error.what();}
        text<<']';
    };
    leaf("source",source,si,sourceBox);leaf("target",target,ti,targetBox);
    const auto comparisons=[&](const char* side,const P::Rectangle& a,const P::Rectangle& b){
        text<<" ["<<side<<" contains4="<<(a.lower.x<=b.lower.x)<<','<<(a.lower.y<=b.lower.y)<<','
            <<(a.upper.x>=b.upper.x)<<','<<(a.upper.y>=b.upper.y)<<']';};
    comparisons("source",sourceBox,targetBox);comparisons("target",targetBox,sourceBox);
    return text.str();
}
}
NativeLeafOverlap NativeLeafOverlap::build(const QuadtreeMesh& source,const NativeLeafDomain& sourceDomain,
    const QuadtreeMesh& target,const NativeLeafDomain& targetDomain)
{
    sourceDomain.validateMeshIdentity(source);targetDomain.validateMeshIdentity(target);
    require(sourceDomain.axisymmetric()==targetDomain.axisymmetric(),"physical measure mismatch");
    require(source.Nx==target.Nx && source.Ny==target.Ny && source.minCoord.x==target.minCoord.x &&
        source.minCoord.y==target.minCoord.y && source.maxCoord.x==target.maxCoord.x && source.maxCoord.y==target.maxCoord.y &&
        source.cellSize.x==target.cellSize.x && source.cellSize.y==target.cellSize.y,"different background lattice/domain");
    const auto& sourceActive=sourceDomain.cellRowByIndex();
    const auto& targetActive=targetDomain.cellRowByIndex();
    require(sourceActive.size()==static_cast<std::size_t>(source.cellIndexNum()) &&
        targetActive.size()==static_cast<std::size_t>(target.cellIndexNum()),"selected maps must cover native indices");
    NativeLeafOverlap out;out.sourceIdentity_=sourceDomain.meshIdentity();out.targetIdentity_=targetDomain.meshIdentity();
    out.bySource_.resize(source.cellIndexNum());out.byTarget_.resize(target.cellIndexNum());
    for(int ci=0;ci<source.cellIndexNum();++ci)
        out.source_.push_back(sourceActive[ci]>=0?rectangle(source,ci):P::Rectangle{});
    for(int ci=0;ci<target.cellIndexNum();++ci)
        out.target_.push_back(targetActive[ci]>=0?rectangle(target,ci):P::Rectangle{});
    const bool rz=sourceDomain.axisymmetric();
    const auto add=[&](int si,int ti,NativeLeafOverlapPiece::Kind kind)
    {
        require(sourceActive[si]>=0 && targetActive[ti]>=0,"selected fluid domains differ across actual overlap");
        const auto box=kind==NativeLeafOverlapPiece::Kind::Refined?out.target_[ti]:out.source_[si];
        const auto index=out.pieces_.size();
        out.pieces_.push_back({kind,si,ti,box,P::capacity(box,rz)});
        out.bySource_[si].push_back(index);out.byTarget_[ti].push_back(index);
    };
    for(int ti=0;ti<target.cellIndexNum();++ti)
    {
        if(targetActive[ti]<0)continue;
        const auto point=target.cells(target.cellIndexToLabel[ti]).centerCoord();
        const auto label=source.findCellLabel(point,0);
        require(label>=0 && label<source.cellLabelCapacity(),"target leaf has no source coverage");
        const auto si=source.cells(label).index();require(si>=0 && si<source.cellIndexNum(),"source lookup is not a leaf");
        const auto sourceBox=sourceActive[si]>=0?out.source_[si]:rectangle(source,si);
        if(contains(sourceBox,out.target_[ti]))
            add(si,ti,equal(sourceBox,out.target_[ti])?NativeLeafOverlapPiece::Kind::Unchanged:NativeLeafOverlapPiece::Kind::Refined);
        else if(!contains(out.target_[ti],sourceBox)) require(false,"source/target leaves are not nested"
            +nestingContext(source,si,sourceBox,target,ti,out.target_[ti]));
    }
    for(int si=0;si<source.cellIndexNum();++si)
    {
        if(sourceActive[si]<0)continue;
        const auto point=source.cells(source.cellIndexToLabel[si]).centerCoord();
        const auto label=target.findCellLabel(point,0);
        require(label>=0 && label<target.cellLabelCapacity(),"source leaf has no target coverage");
        const auto ti=target.cells(label).index();require(ti>=0 && ti<target.cellIndexNum(),"target lookup is not a leaf");
        const auto targetBox=targetActive[ti]>=0?out.target_[ti]:rectangle(target,ti);
        if(contains(targetBox,out.source_[si]) && !equal(out.source_[si],targetBox))
            add(si,ti,NativeLeafOverlapPiece::Kind::Coarsened);
        else if(!contains(out.source_[si],targetBox)) require(false,"target/source leaves are not nested"
            +nestingContext(source,si,out.source_[si],target,ti,targetBox));
    }
    const auto check=[&](const auto& groups,const auto& boxes,const auto& active)
    {
        for(std::size_t cell=0;cell<groups.size();++cell)
        {
            if(active[cell]<0){require(groups[cell].empty(),"inactive leaf has physical overlap");continue;}
            require(!groups[cell].empty(),"uncovered selected leaf");
            long double sum=0,correction=0;
            for(const auto index:groups[cell])
            {const auto y=out.pieces_[index].capacity-correction,t=sum+y;correction=(t-sum)-y;sum=t;}
            const auto capacity=P::capacity(boxes[cell],rz),residual=std::abs(sum-capacity)/capacity;
            out.coverage_=std::max(out.coverage_,residual);
            require(residual<=512*std::numeric_limits<long double>::epsilon(),"incomplete/duplicate physical leaf coverage");
        }
    };
    check(out.bySource_,out.source_,sourceActive);check(out.byTarget_,out.target_,targetActive);return out;
}
void NativeLeafOverlap::validateFor(const QuadtreeMesh& source,const QuadtreeMesh& target)const
{require(sourceIdentity_ && targetIdentity_,"unbound overlap");sourceIdentity_->validate(source);targetIdentity_->validate(target);}
}
