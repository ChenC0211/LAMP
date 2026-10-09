#include "Core/Geometry/ActiveHermiteTransfer.h"
#include "Core/Geometry/Detail/LevelSetLeafHermiteScalar.h"
#include "Core/Geometry/Detail/NativeHermiteTransferSampling.h"
#include "Core/Mesh/CartesianLeafMesh.h"
#include <map>

namespace ActiveHermiteTransfer
{
namespace
{
namespace Scalar = Quantum::Geometry::Detail::LevelSetHermiteScalar;
namespace Exact = Quantum::Geometry::Detail::ExactExpansion;
using Source = Scalar::LeafSnapshot;
using Sampler = Quantum::Geometry::Detail::NativeHermiteTransfer::DatumSampler<Source>;
using Box = std::array<Vector2D,2>;
using Point = std::array<double,2>;
void require(bool valid, const char* reason)
{ if (!valid) throw std::invalid_argument(std::string("[ActiveHermiteTransfer] ")+reason); }
Point key(Vector2D p) { return {{p.x,p.y}}; }
Box bounds(const QuadtreeMesh& mesh,int c)
{
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
    return {{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()}};
}
Box intersection(const Box& a,const Box& b)
{ return {{{std::max(a[0].x,b[0].x),std::max(a[0].y,b[0].y)},
           {std::min(a[1].x,b[1].x),std::min(a[1].y,b[1].y)}}}; }
bool nonempty(const Box& b) { return b[0].x<b[1].x && b[0].y<b[1].y; }
bool contains(const Box& a,const Box& b)
{ return a[0].x<=b[0].x && a[0].y<=b[0].y && a[1].x>=b[1].x && a[1].y>=b[1].y; }
double up(double x)
{
    const double result=std::nextafter(x,std::numeric_limits<double>::infinity());
    require(std::isfinite(result),"nonfinite scalar error bound");
    return result;
}
double addUpper(double a,double b) { return a==0?b:b==0?a:up(a+b); }
double multiplyUpper(double a,double b) { return a==0 || b==0?0:up(a*b); }
double defect(const Source& a,int ca,const Source& b,int cb,Vector2D p,unsigned dx,unsigned dy)
{
    const auto av=a.restrictPhysical(ca,p,p,dx,dy),bv=b.restrictPhysical(cb,p,p,dx,dy);
    const auto numerator=Exact::subtract(Exact::multiply(av.numerator[0],bv.denominator),
        Exact::multiply(bv.numerator[0],av.denominator));
    if (Exact::sign(numerator)==0) return 0;
    const auto n=Exact::enclosure(numerator);
    const auto d=Exact::enclosure(Exact::multiply(av.denominator,bv.denominator));
    require(d.lo>0,"error bound denominator has no positive enclosure");
    return up(std::max(std::abs(n.lo),std::abs(n.hi))/d.lo);
}
double differenceBound(const Source& a,int ca,const Source& b,int cb,const Box& box)
{
    std::array<double,4> jets{};
    for(double x:{box[0].x,box[1].x}) for(double y:{box[0].y,box[1].y})
        for(unsigned j=0;j<4;++j)
            jets[j]=std::max(jets[j],defect(a,ca,b,cb,{x,y},j%2,j/2));
    // For cubic Hermite basis on [0,1], h00+h01=1 and
    // |h10|+|h11|=t(1-t)<=1/4. Tensor products give this
    // whole-rectangle bound, not just a sampled error estimate.
    const double hx=up(box[1].x-box[0].x),hy=up(box[1].y-box[0].y);
    const double wx=multiplyUpper(hx,.25),wy=multiplyUpper(hy,.25);
    return addUpper(addUpper(jets[0],multiplyUpper(wx,jets[1])),
        addUpper(multiplyUpper(wy,jets[2]),multiplyUpper(multiplyUpper(wx,wy),jets[3])));
}
}

Result prepare(const QuadtreeMesh& sourceMesh,const Quantum::Geometry::ActiveHermitePatch& patch,
    const QuadtreeMesh& target,const std::vector<std::uint8_t>& active,bool radial)
{
    Exact::validateEnvironment();
    require(key(sourceMesh.minCoord)==key(target.minCoord) && key(sourceMesh.maxCoord)==key(target.maxCoord)
        && sourceMesh.Nx==target.Nx && sourceMesh.Ny==target.Ny,"different background lattice/domain");
    require(active.size()==static_cast<std::size_t>(target.cellIndexNum()),"wrong target activity size");
    for(auto bit:active) require(bit<=1,"target activity is not binary");
    validateCartesianLeafMesh(target);
    const Source source(sourceMesh,patch,radial);
    std::vector<int> oldSector(sourceMesh.cellIndexNum(),-1);
    for(std::size_t i=0;i<patch.cornerSectors.size();++i)oldSector[patch.cornerSectors[i].cellIndex]=static_cast<int>(i);
    Result result;
    result.patch.activeCells=active;
    result.sourceCellsByTarget.resize(target.cellIndexNum());
    std::vector<std::uint8_t> used(target.vertIndexNum(),0),oldUsed(sourceMesh.vertIndexNum(),0);
    for(int c=0;c<sourceMesh.cellIndexNum();++c) if(source.cellActive(c))
        for(int label:sourceMesh.cells(sourceMesh.cellIndexToLabel[c]).cellToVert)
            oldUsed[sourceMesh.verts(label).index()]=1;
    for(int c=0;c<target.cellIndexNum();++c) {
        const auto box=bounds(target,c);
        auto& overlapping=result.sourceCellsByTarget[c];
        for(int old:source.leafCandidates(box[0],box[1])) {
            const auto previous=source.leafBounds(old);
            if(!nonempty(intersection(box,previous))) continue;
            require(contains(box,previous) || contains(previous,box),"leaves are not nested parent/child boxes");
            require(source.cellActive(old)==(active[c]!=0),"remap changes fixed active/solid coverage");
            overlapping.push_back(old);
        }
        require(!overlapping.empty(),"target leaf has no source coverage");
        std::sort(overlapping.begin(),overlapping.end());
        if(active[c]) for(int label:target.cells(target.cellIndexToLabel[c]).cellToVert)
            used[target.verts(label).index()]=1;
    }
    std::map<Point,int> oldVertices,nextVertices;
    for(int v=0;v<sourceMesh.vertIndexNum();++v)
        if(oldUsed[v]) oldVertices.emplace(key(sourceMesh.verts(sourceMesh.vertIndexToLabel[v]).coord()),v);
    for(int v=0;v<target.vertIndexNum();++v)
        nextVertices.emplace(key(target.verts(target.vertIndexToLabel[v]).coord()),v);
    result.sourceToTargetVertex.assign(sourceMesh.vertIndexNum(),-1);
    for(int v=0;v<sourceMesh.vertIndexNum();++v) {
        const auto found=nextVertices.find(key(sourceMesh.verts(sourceMesh.vertIndexToLabel[v]).coord()));
        if(found!=nextVertices.end()) result.sourceToTargetVertex[v]=found->second;
    }
    const double nan=std::numeric_limits<double>::quiet_NaN();
    result.patch.nodes.assign(target.vertIndexNum(),{nan,nan,nan,nan});
    const Sampler sampler(source);
    for(int v=0;v<target.vertIndexNum();++v) {
        const auto& vertex=target.verts(target.vertIndexToLabel[v]);
        if(!used[v] || vertex.supCell>=0) continue;
        const auto found=oldVertices.find(key(vertex.coord()));
        if(found!=oldVertices.end() && sourceMesh.verts(sourceMesh.vertIndexToLabel[found->second]).supCell<0) {
            result.patch.nodes[v]=patch.nodes[found->second];++result.copiedVertices;
        } else {
            std::array<double,4> jet;
            if(found!=oldVertices.end()) jet=source.nodeJet(found->second);
            else {
                auto point=sampler.locate(vertex.coord(),vertex.coord());
                if(oldSector[point.cell]>=0) {
                    // At a C0 sector seam the cap/ordinary chart owns the
                    // shared jet. Only the pure sector descendant may change.
                    for(int c:source.leafCandidates(vertex.coord(),vertex.coord())) {
                        if(!source.cellActive(c) || oldSector[c]>=0)continue;
                        const auto b=source.leafBounds(c);const auto p=vertex.coord();
                        if(p.x<b[0].x || p.x>b[1].x || p.y<b[0].y || p.y>b[1].y)continue;
                        const auto datum=[&](unsigned dx,unsigned dy){const auto value=source.restrictPhysical(c,p,p,dx,dy);
                            const double sampled=Exact::estimate(value.numerator[0])/Exact::estimate(value.denominator);
                            require(std::isfinite(sampled) && (sampled!=0 || Exact::sign(value.numerator[0])==0),
                                "nonfinite or underflowed non-sector source datum");return sampled;};
                        jet={{datum(0,0),datum(1,0),datum(0,1),datum(1,1)}};point.cell=-1;break;
                    }
                }
                if(point.cell>=0)jet={{sampler.value(point),sampler.value(point,1,0),sampler.value(point,0,1),sampler.value(point,1,1)}};
            }
            result.patch.nodes[v]={jet[0],jet[1],jet[2],jet[3]};++result.sampledVertices;
        }
    }
    for(const auto& sector:patch.cornerSectors) {
        const int vertex=result.sourceToTargetVertex.at(sector.vertexIndex);
        require(vertex>=0 && target.verts(target.vertIndexToLabel[vertex]).supCell<0,"remap lost a regular declared corner");
        const auto p=target.verts(target.vertIndexToLabel[vertex]).coord();int owner=-1;
        for(int c=0;c<target.cellIndexNum();++c)if(active[c]) {
            const auto& overlap=result.sourceCellsByTarget[c];
            if(std::find(overlap.begin(),overlap.end(),sector.cellIndex)==overlap.end())continue;
            const auto b=bounds(target,c);
            if((p.x==b[0].x || p.x==b[1].x) && (p.y==b[0].y || p.y==b[1].y)) {
                require(owner<0,"remap has ambiguous corner-sector ownership");owner=c;
            }
        }
        require(owner>=0,"remap lost its corner-sector leaf");
        auto next=sector;next.cellIndex=owner;next.vertexIndex=vertex;
        result.patch.cornerSectors.push_back(next);
    }
    const Source candidate(target,result.patch,radial);
    result.scalarDifferenceBound.assign(target.cellIndexNum(),nan);
    result.pureSectorExtension.assign(target.cellIndexNum(),0);
    for(int c=0;c<target.cellIndexNum();++c) if(active[c]) {
        double bound=0;
        for(int old:result.sourceCellsByTarget[c]) bound=std::max(bound,
            differenceBound(source,old,candidate,c,intersection(source.leafBounds(old),candidate.leafBounds(c))));
        result.scalarDifferenceBound[c]=bound;
        result.maxScalarDifferenceBound=std::max(result.maxScalarDifferenceBound,bound);
        const auto& overlap=result.sourceCellsByTarget[c];
        if(overlap.size()==1 && oldSector[overlap[0]]>=0) {
            const auto& sector=patch.cornerSectors[oldSector[overlap[0]]];
            const auto point=sourceMesh.verts(sourceMesh.vertIndexToLabel[sector.vertexIndex]).coord();
            const auto b=bounds(target,c);unsigned mask=0;
            if(point.x==b[0].x && point.y==b[0].y)mask=1;
            if(point.x==b[1].x && point.y==b[0].y)mask=2;
            if(point.x==b[1].x && point.y==b[1].y)mask=4;
            if(point.x==b[0].x && point.y==b[1].y)mask=8;
            require(candidate.cellSignExceptVertices(c,mask)==sector.phaseSign,
                "corner-sector reextension changed phase or introduced an unbound zero");
            result.pureSectorExtension[c]=sector.phaseSign;
        } else result.maxNonExtensionDifferenceBound=std::max(result.maxNonExtensionDifferenceBound,bound);
    }
    return result;
}
}
