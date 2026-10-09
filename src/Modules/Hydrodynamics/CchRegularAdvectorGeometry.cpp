#include "Modules/Hydrodynamics/Detail/CchRegularAdvectorGeometry.h"
#include "Modules/Hydrodynamics/CchSpatialAdvector.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <Eigen/SVD>
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

namespace Quantum::Hydrodynamics::Detail
{
namespace
{
void require(bool ok,const char* why)
{if(!ok)throw std::invalid_argument(std::string("[CchRegularAdvectorGeometry] ")+why);}
void fp()
{
    require(std::fegetround()==FE_TONEAREST,"round-to-nearest required; mode not changed");
#if defined(_M_X64) || defined(__SSE2__)
    require((_mm_getcsr() & ((1U<<15)|(1U<<6)))==0,"gradual underflow required; mode not changed");
#endif
}
double coord(Vector2D p,int d){return d==0?p.x:p.y;}
constexpr int exponents[6][2]={{0,0},{1,0},{0,1},{2,0},{1,1},{0,2}};
constexpr double basisPower[5][5]={{1,0,-11,18,-8},{0,0,-5,14,-8},{0,1,-4,5,-2},{0,0,1,-3,2},{0,0,16,-32,16}};
std::array<double,5> basis(double s,bool derivative=false)
{
    std::array<double,5> out{};
    for(int i=0;i<5;++i)for(int k=4;k>=(derivative?1:0);--k)
        out[i]=out[i]*s+basisPower[i][k]*(derivative?k:1);
    return out;
}
constexpr double gauss[5]={.0469100770306680036,.230765344947158454,.5,.769234655052841546,.953089922969331996};
constexpr double gaussWeight[5]={.118463442528094544,.239314335249683234,.284444444444444444,.239314335249683234,.118463442528094544};
std::array<double,6> mean(const CchRegularFace& f,Vector2D p,Vector2D h,bool rz)
{
    std::array<double,6> out{};const double radius=(f.lower.x+f.upper.x)/2;
    for(int q=0;q<5;++q){const auto x=f.lower+(f.upper-f.lower)*gauss[q];
        const double a=(x.x-p.x)/h.x,b=(x.y-p.y)/h.y;
        const std::array<double,6> v{1,a,b,a*a,a*b,b*b};
        const double w=gaussWeight[q]*(rz?x.x/radius:1);
        for(int k=0;k<6;++k)out[k]+=w*v[k];}
    return out;
}
std::vector<int> columns(Vector2D p,int component,bool rz,bool wall)
{
    std::vector<int> out;
    for(int k=0;k<6;++k){if(wall && exponents[k][component]==0)continue;
        if(rz && p.x==0 && exponents[k][0]%2!=(component==0?1:0))continue;
        out.push_back(k);}
    return out;
}
// Orthonormal coefficient directions for the declared lower physical space,
// intersected with the same axis/wall conditions as the full P2 columns.
Eigen::MatrixXd prioritySpace(const CchRegularRow& row,bool rz,int degree=2)
{
    std::vector<Eigen::VectorXd> vectors;
    const auto add=[&](int a,double x,int b=-1,double y=0.) {
        Eigen::VectorXd v=Eigen::VectorXd::Zero(row.columns.size());
        for(auto term:{std::pair<int,double>{a,x},{b,y}})if(term.second!=0){
            const auto it=std::find(row.columns.begin(),row.columns.end(),term.first);
            if(it==row.columns.end())return; // This entire physical direction is incompatible.
            v[it-row.columns.begin()]=term.second;
        }
        if(v.norm()>0)vectors.push_back(v.normalized());
    };
    if(!rz){add(0,1);if(degree>=1){add(1,1);add(2,1);}}
    else if(row.component==0){if(degree>=1)add(0,row.point.x,1,row.scale.x);} // physical r.
    else {
        add(0,1);if(degree>=1)add(2,1);
        // Remove the constant component of r^2; this spans the same space
        // with 1 and avoids a nearly parallel pair far from the axis.
        if(degree>=2&&!row.zeroNormal)add(1,2*row.point.x,3,row.scale.x);
    }
    Eigen::MatrixXd out(row.columns.size(),vectors.size());
    for(std::size_t k=0;k<vectors.size();++k)out.col(k)=vectors[k];return out;
}
// Resolve observable lower-degree directions first. A missing transverse
// slope must not contaminate a constant whose observations are available.
// The threshold belongs to the original lower observation matrix: scaling a
// roundoff-only projected residual by its own largest value invents rank.
Eigen::MatrixXd identifiableInverse(const CchRegularRow& row,bool rz,
    const Eigen::MatrixXd& observations,double threshold,std::size_t& svdCount)
{
    const int count=static_cast<int>(observations.rows());
    Eigen::MatrixXd inverse=Eigen::MatrixXd::Zero(observations.cols(),count);
    Eigen::MatrixXd seen(count,0);int previous=0;
    for(int degree=0;degree<=2;++degree){
        const Eigen::MatrixXd space=prioritySpace(row,rz,degree);
        const int added=static_cast<int>(space.cols())-previous;
        if(added==0)continue;
        const Eigen::MatrixXd directions=space.rightCols(added),values=observations*directions;
        previous=static_cast<int>(space.cols());
        const Eigen::MatrixXd residual=values-seen*(seen.transpose()*values);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(residual,Eigen::ComputeThinU|Eigen::ComputeThinV);++svdCount;
        int rank=0;while(rank<svd.singularValues().size()&&svd.singularValues()[rank]>threshold)++rank;
        const Eigen::MatrixXd coefficients=directions-inverse*values;
        for(int k=0;k<rank;++k)inverse.noalias()+=coefficients*svd.matrixV().col(k)*(svd.matrixU().col(k).transpose()/svd.singularValues()[k]);
        const int old=static_cast<int>(seen.cols());seen.conservativeResize(count,old+rank);
        seen.rightCols(rank)=svd.matrixU().leftCols(rank);
    }
    require(seen.cols()==row.priorityRank,"inconsistent identifiable lower-space rank");
    return inverse;
}
struct Datum
{
    CchRegularRow::Kind kind;
    int index=0;
    Vector2D point{0,0};
    std::set<int> seeds;
    std::array<std::vector<int>,2> wallFaces;
    std::array<int,2> base{},wall{};
    std::array<int,2> normalPort{{-1,-1}};
};
void sortUnique(std::vector<double>& v)
{std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());}
}
struct CchRegularAdvectorGeometry::Data
{
    explicit Data(CchMomentumGeometry g):geometry(std::move(g)){}
    CchMomentumGeometry geometry;
    std::vector<Vector2D> vertices;
    std::vector<int> vertexParent,usedVertices;
    std::vector<CchRegularCell> cells;
    std::vector<CchRegularEdge> edges;
    std::vector<CchRegularFace> faces;
    std::vector<CchRegularRow> rows;
    std::vector<CchRegularNormalPort> normalPorts;
    std::vector<Datum> datums;
    std::map<std::pair<CchRegularRow::Kind,int>,std::size_t> datumIndex;
    CchRegularPreparationReport report;
};
CchRegularAdvectorGeometry::CchRegularAdvectorGeometry(std::shared_ptr<const Data> d):data_(std::move(d)){}
const CchMomentumGeometry& CchRegularAdvectorGeometry::momentumGeometry() const noexcept{return data_->geometry;}
const std::vector<CchRegularCell>& CchRegularAdvectorGeometry::cells() const noexcept{return data_->cells;}
const std::vector<CchRegularEdge>& CchRegularAdvectorGeometry::edges() const noexcept{return data_->edges;}
const std::vector<CchRegularFace>& CchRegularAdvectorGeometry::faces() const noexcept{return data_->faces;}
const std::vector<CchRegularRow>& CchRegularAdvectorGeometry::rows() const noexcept{return data_->rows;}
const std::vector<CchRegularNormalPort>& CchRegularAdvectorGeometry::normalPorts() const noexcept{return data_->normalPorts;}
const CchRegularPreparationReport& CchRegularAdvectorGeometry::report() const noexcept{return data_->report;}
bool CchRegularAdvectorGeometry::sameIdentity(const CchRegularAdvectorGeometry& other) const noexcept{return data_==other.data_;}
void CchRegularAdvectorGeometry::validateFor(const QuadtreeMesh& mesh,const CchMomentumGeometry& expected) const
{fp();momentumGeometry().validateMesh(mesh);require(momentumGeometry().sameIdentity(expected),"different momentum geometry instance");}

CchRegularAdvectorGeometry CchRegularAdvectorGeometry::prepare(const QuadtreeMesh& mesh,CchMomentumGeometry geometry)
{
    fp();geometry.validateMesh(mesh);auto d=std::make_shared<Data>(std::move(geometry));
    const auto& l=d->geometry.layout();const bool rz=l.domain().axisymmetric();
    d->vertices.resize(mesh.vertIndexNum());d->vertexParent.assign(mesh.vertIndexNum(),-1);
    d->cells.resize(l.pressureCount());d->faces.resize(l.velocityCount());
    std::map<int,std::set<int>> vertexCells;
    std::map<std::array<int,2>,std::set<int>> edgeCells;
    constexpr int corners[4][2]={{0,1},{3,2},{0,3},{1,2}};
    for(int p=0;p<l.pressureCount();++p){auto& c=d->cells[p];c.nativeCell=l.nativeCell(p);
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c.nativeCell]);
        for(int k=0;k<4;++k){const auto& v=mesh.verts(cell.cellToVert[k]);c.vertices[k]=v.index();d->vertices[v.index()]=v.coord();vertexCells[v.index()].insert(p);}
        c.lower=d->vertices[c.vertices[0]];c.upper=d->vertices[c.vertices[2]];
        for(const auto& k:corners)edgeCells[{c.vertices[k[0]],c.vertices[k[1]]}].insert(p);
    }
    std::vector<std::set<int>> seeds;
    std::map<std::array<int,2>,int> edgeId;
    std::map<std::pair<int,double>,std::vector<int>> groups;
    for(const auto& kv:edgeCells){CchRegularEdge e;e.vertices=kv.first;e.lower=d->vertices[e.vertices[0]];e.upper=d->vertices[e.vertices[1]];
        e.normalAxis=e.lower.y==e.upper.y?1:0;const int i=static_cast<int>(d->edges.size());edgeId[kv.first]=i;
        d->edges.push_back(e);seeds.push_back(kv.second);groups[{e.normalAxis,coord(e.lower,e.normalAxis)}].push_back(i);}
    for(auto& c:d->cells)for(int k=0;k<4;++k)c.edges[k]=edgeId.at({c.vertices[corners[k][0]],c.vertices[corners[k][1]]});
    for(const auto& group:groups)for(int i:group.second){auto& e=d->edges[i];const int t=1-e.normalAxis;double longest=coord(e.upper-e.lower,t);
        for(int j:group.second){const auto& other=d->edges[j];const double length=coord(other.upper-other.lower,t);
            if(length>longest && coord(other.lower,t)<=coord(e.lower,t) && coord(other.upper,t)>=coord(e.upper,t)){e.parent=j;longest=length;}}}
    for(std::size_t i=0;i<d->edges.size();++i){const auto& e=d->edges[i];if(e.parent<0)continue;const auto& master=d->edges[e.parent];
        seeds[e.parent].insert(seeds[i].begin(),seeds[i].end());const int t=1-e.normalAxis;
        for(int v:e.vertices)if(coord(master.lower,t)<coord(d->vertices[v],t) && coord(d->vertices[v],t)<coord(master.upper,t)){
            require(d->vertexParent[v]<0 || d->vertexParent[v]==e.parent,"incompatible active master sides");d->vertexParent[v]=e.parent;}}
    std::vector<std::set<int>> neighbors(d->cells.size());
    for(int u=0;u<l.velocityCount();++u){auto& f=d->faces[u];const auto actual=l.domain().face(mesh,l.nativeFace(u));
        const auto& edge=actual.horizontal?mesh.hEdges(actual.edgeLabel):mesh.vEdges(actual.edgeLabel);
        f.normalAxis=actual.horizontal?1:0;f.lower=mesh.verts(edge.edgeToVert[0]).coord();f.upper=mesh.verts(edge.edgeToVert[1]).coord();
        if(coord(f.lower,1-f.normalAxis)>coord(f.upper,1-f.normalAxis))std::swap(f.lower,f.upper);
        f.negativePressure=actual.negativeCell<0?-1:l.pressureDof(actual.negativeCell);f.positivePressure=actual.positiveCell<0?-1:l.pressureDof(actual.positiveCell);
        f.area=d->geometry.metrics().faceMeasures()[u];f.external=actual.outwardSign!=0;
        const int ids[2]={f.negativePressure,f.positivePressure};
        for(int side=0;side<2;++side)if(ids[side]>=0){auto& c=d->cells[ids[side]];c.faces.push_back(u);f.supportWidth[side]=coord(c.upper-c.lower,f.normalAxis)/2;
            c.supportKnots[f.normalAxis].push_back(coord(c.lower+c.upper,f.normalAxis)/2);
            for(int k=0;k<2;++k){c.supportKnots[k].push_back(coord(f.lower,k));c.supportKnots[k].push_back(coord(f.upper,k));}}
        if(ids[0]>=0 && ids[1]>=0){neighbors[ids[0]].insert(ids[1]);neighbors[ids[1]].insert(ids[0]);}
        if(f.area>0){const auto& cell=d->cells[ids[0]>=0?ids[0]:ids[1]];const auto h=cell.upper-cell.lower;
            for(int q=0;q<5;++q){const auto x=f.lower+(f.upper-f.lower)*gauss[q];const auto bx=basis((x.x-cell.lower.x)/h.x),by=basis((x.y-cell.lower.y)/h.y);
                const double w=gaussWeight[q]*(rz?x.x/((f.lower.x+f.upper.x)/2):1);
                for(int y=0;y<5;++y)for(int xindex=0;xindex<5;++xindex)f.meanWeights[5*y+xindex]+=w*bx[xindex]*by[y];}}
    }
    for(auto& cell:d->cells)for(auto& knots:cell.supportKnots)sortUnique(knots);
    // Coplanar physical subfaces, partitioned at gaps in the selected domain.
    // Preparation alone owns the moment solves; apply only multiplies weights.
    using LineKey=std::pair<int,double>;
    std::map<LineKey,std::vector<int>> lines;
    std::map<LineKey,std::vector<std::pair<int,int>>> segments;
    for(std::size_t u=0;u<d->faces.size();++u){const auto& f=d->faces[u];
        if(f.area>0)lines[{f.normalAxis,coord(f.lower,f.normalAxis)}].push_back(static_cast<int>(u));}
    for(auto& item:lines){auto& ids=item.second;const int t=1-item.first.first;
        std::sort(ids.begin(),ids.end(),[&](int a,int b){return coord(d->faces[a].lower,t)<coord(d->faces[b].lower,t);});
        int begin=0;for(int i=1;i<=static_cast<int>(ids.size());++i)
            if(i==static_cast<int>(ids.size())||coord(d->faces[ids[i-1]].upper,t)!=coord(d->faces[ids[i]].lower,t)){
                segments[item.first].emplace_back(begin,i);begin=i;}}
    for(const auto& kv:vertexCells){d->usedVertices.push_back(kv.first);if(d->vertexParent[kv.first]<0)d->datums.push_back({CchRegularRow::Kind::Vertex,kv.first,d->vertices[kv.first],kv.second});}
    for(std::size_t e=0;e<d->edges.size();++e)if(d->edges[e].parent<0)d->datums.push_back({CchRegularRow::Kind::Edge,static_cast<int>(e),(d->edges[e].lower+d->edges[e].upper)*.5,seeds[e]});
    for(std::size_t p=0;p<d->cells.size();++p)d->datums.push_back({CchRegularRow::Kind::Center,static_cast<int>(p),(d->cells[p].lower+d->cells[p].upper)*.5,{static_cast<int>(p)}});
    // Late boundary binding selects a precomputed finite mode, never a fit.
    std::vector<std::pair<int,std::vector<int>>> wallEdges;
    for(std::size_t e=0;e<d->edges.size();++e){const auto& edge=d->edges[e];if(edge.parent>=0)continue;std::vector<int> walls;const int a=edge.normalAxis,t=1-a;
        for(std::size_t u=0;u<d->faces.size();++u){const auto& f=d->faces[u];if(f.external && f.area>0 && f.normalAxis==a && coord(f.lower,a)==coord(edge.lower,a) && std::min(coord(f.upper,t),coord(edge.upper,t))>std::max(coord(f.lower,t),coord(edge.lower,t)))walls.push_back(static_cast<int>(u));}
        if(!walls.empty())wallEdges.emplace_back(static_cast<int>(e),std::move(walls));}
    for(auto& datum:d->datums){Vector2D h{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
        for(int p:datum.seeds){const auto w=d->cells[p].upper-d->cells[p].lower;h.x=std::min(h.x,w.x);h.y=std::min(h.y,w.y);}
        for(const auto& item:wallEdges){const auto& e=d->edges[item.first];const int a=e.normalAxis,t=1-a;
            if(coord(datum.point,a)==coord(e.lower,a) && coord(datum.point,t)>=coord(e.lower,t) && coord(datum.point,t)<=coord(e.upper,t))
                datum.wallFaces[a].insert(datum.wallFaces[a].end(),item.second.begin(),item.second.end());}
        for(int component=0;component<2;++component)for(int mode=0;mode<2;++mode){const auto cols=columns(datum.point,component,rz,mode==1);
            if(mode==1 && (datum.wallFaces[component].empty() || cols==columns(datum.point,component,rz,false))){datum.wall[component]=datum.base[component];continue;}
            CchRegularRow row;row.kind=datum.kind;row.index=datum.index;row.component=component;row.zeroNormal=mode==1;row.point=datum.point;row.scale=h;row.columns=cols;
            std::set<int> support=datum.seeds;Eigen::MatrixXd m,inverse;Eigen::JacobiSVD<Eigen::MatrixXd> svd;
            for(int ring=1;ring<=4;++ring){const auto before=support;for(int p:before)support.insert(neighbors[p].begin(),neighbors[p].end());
                std::set<int> faces;for(int p:support)for(int u:d->cells[p].faces)if(d->faces[u].normalAxis==component && d->faces[u].area>0)faces.insert(u);
                require(!faces.empty(),"no native component observations");row.support.assign(faces.begin(),faces.end());m.resize(faces.size(),cols.size());
                for(std::size_t i=0;i<row.support.size();++i){const auto moments=mean(d->faces[row.support[i]],datum.point,h,rz);for(std::size_t j=0;j<cols.size();++j)m(i,j)=moments[cols[j]];}
                svd.compute(m,Eigen::ComputeThinU|Eigen::ComputeFullV);++d->report.svdCount;
                svd.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max<std::size_t>(row.support.size(),6));row.rank=static_cast<int>(svd.rank());row.rings=ring;
                bool exhausted=true;for(int p:support)for(int q:neighbors[p])if(!support.count(q))exhausted=false;
                if(row.rank==static_cast<int>(cols.size()) || ring==4 || exhausted)break;
            }
            inverse=Eigen::MatrixXd::Zero(cols.size(),row.support.size());
            for(int k=0;k<row.rank;++k)inverse.noalias()+=svd.matrixV().col(k)*(svd.matrixU().col(k).transpose()/svd.singularValues()[k]);
            const Eigen::MatrixXd lower=prioritySpace(row,rz),observed=m*lower;
            row.priorityDimension=static_cast<int>(lower.cols());row.priorityRank=row.priorityDimension;
            if(row.rank<static_cast<int>(cols.size()) && lower.cols()>0){
                Eigen::JacobiSVD<Eigen::MatrixXd> low(observed,Eigen::ComputeThinU|Eigen::ComputeThinV);++d->report.svdCount;
                low.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max<std::size_t>(row.support.size(),6));row.priorityRank=static_cast<int>(low.rank());
                if(row.priorityRank==row.priorityDimension){
                    Eigen::MatrixXd right=Eigen::MatrixXd::Zero(lower.cols(),row.support.size());
                    for(int k=0;k<row.priorityRank;++k)right.noalias()+=low.matrixV().col(k)*(low.matrixU().col(k).transpose()/low.singularValues()[k]);
                    // Preserve the already accepted fully identifiable weights.
                    const Eigen::MatrixXd correction=(lower-inverse*observed)*right;
                    inverse+=correction;
                }else{
                    const double threshold=low.threshold()*low.singularValues()[0];
                    const Eigen::MatrixXd hierarchy=identifiableInverse(row,rz,m,threshold,d->report.svdCount);
                    // Observation-null correction preserves the full least-squares
                    // fit while reproducing the identifiable degree hierarchy.
                    const Eigen::MatrixXd correction=hierarchy-inverse*(m*hierarchy);
                    inverse+=correction;
                }
            }
            row.priorityDeficit=(inverse*observed-lower).norm();
            row.coefficientNorm=inverse.norm();row.deficit=(Eigen::MatrixXd::Identity(cols.size(),cols.size())-inverse*m).norm();
            if(row.rank<static_cast<int>(cols.size())){++d->report.rankDeficient;const auto null=svd.matrixV().col(cols.size()-1);for(std::size_t j=0;j<cols.size();++j)row.nullPolynomial[cols[j]]=null[j];row.nullResidual=(m*null).cwiseAbs().maxCoeff();}
            constexpr int outputs[4]={0,1,2,4};const double scales[4]={1,h.x,h.y,h.x*h.y};
            for(int k=0;k<4;++k){auto& w=row.weights[k];w.assign(row.support.size(),0);const auto it=std::find(cols.begin(),cols.end(),outputs[k]);
                if(it!=cols.end())for(std::size_t i=0;i<w.size();++i){w[i]=inverse(it-cols.begin(),i)/scales[k];require(std::isfinite(w[i]),"nonfinite prepared weight");}}
            if(mode==0&&datum.kind==CchRegularRow::Kind::Vertex){
                CchRegularNormalPort port;port.vertex=datum.index;port.component=component;
                const int tangent=1-component;const double target=coord(datum.point,tangent),scale=coord(h,tangent);
                // Axis radial normal velocity and its axial derivative vanish.
                if(!(rz&&component==0&&datum.point.x==0)){
                    const LineKey key{component,coord(datum.point,component)};const auto found=lines.find(key);
                    require(found!=lines.end(),"missing normal-port coordinate line");const auto& ids=found->second;
                    const auto& spans=segments.at(key);
                    const auto span=std::lower_bound(spans.begin(),spans.end(),target,[&](const auto& s,double x){return coord(d->faces[ids[s.second-1]].upper,tangent)<x;});
                    require(span!=spans.end()&&coord(d->faces[ids[span->first]].lower,tangent)<=target,"normal port outside contiguous native span");
                    const auto position=std::lower_bound(ids.begin()+span->first,ids.begin()+span->second,target,[&](int u,double x){return coord(d->faces[u].lower+d->faces[u].upper,tangent)/2<x;});
                    int left=static_cast<int>(position-ids.begin())-1,right=left+1;
                    std::vector<int> observations;
                    while(observations.size()<4&&(left>=span->first||right<span->second)){
                        if(right==span->second||(left>=span->first&&target-coord(d->faces[ids[left]].lower+d->faces[ids[left]].upper,tangent)/2<=coord(d->faces[ids[right]].lower+d->faces[ids[right]].upper,tangent)/2-target))observations.push_back(ids[left--]);
                        else observations.push_back(ids[right++]);}
                    const bool evenAxis=rz&&component==1&&datum.point.x==0;
                    const std::vector<int> powers=evenAxis?std::vector<int>{0,2}:std::vector<int>{0,1,2,3};
                    if(observations.size()>powers.size())observations.resize(powers.size());
                    const int count=static_cast<int>(observations.size()),degree=static_cast<int>(powers.size());
                    Eigen::MatrixXd moments=Eigen::MatrixXd::Zero(count,degree);
                    for(int i=0;i<count;++i){const auto& f=d->faces[observations[i]];
                        for(int q=0;q<5;++q){const auto x=f.lower+(f.upper-f.lower)*gauss[q];const double z=(coord(x,tangent)-target)/scale;
                            const double weight=gaussWeight[q]*(rz?x.x/((f.lower.x+f.upper.x)/2):1);
                            for(int k=0;k<degree;++k)moments(i,k)+=weight*std::pow(z,powers[k]);}}
                    Eigen::JacobiSVD<Eigen::MatrixXd> portSvd(moments.leftCols(count),Eigen::ComputeThinU|Eigen::ComputeThinV);++d->report.svdCount;
                    portSvd.setThreshold(128*std::numeric_limits<double>::epsilon()*std::max(count,6));
                    require(portSvd.rank()==count,"deficient native line moments");port.observedModes=count;
                    const Eigen::MatrixXd rightInverse=portSvd.solve(Eigen::MatrixXd::Identity(count,count));
                    // Missing higher line modes are completed by the existing
                    // identifiable P2 seed. Observable line modes always obey
                    // their physical moments. A cubic seed is exactly zero.
                    std::map<int,std::array<double,2>> sparse;
                    for(int i=0;i<count;++i){auto& w=sparse[observations[i]];w[0]+=rightInverse(0,i);
                        if(!evenAxis&&count>1)w[1]+=rightInverse(1,i)/scale;}
                    for(int k=count;k<degree;++k){const int power=powers[k];if(power>2)continue;
                        const int polynomial=power==0?0:(power==1?1+tangent:(tangent==0?3:5));
                        const auto col=std::find(cols.begin(),cols.end(),polynomial);if(col==cols.end())continue;
                        const Eigen::VectorXd correction=-rightInverse*moments.col(k);
                        for(std::size_t i=0;i<row.support.size();++i){const double seed=inverse(col-cols.begin(),i);auto& w=sparse[row.support[i]];
                            w[0]+=correction[0]*seed;
                            if(!evenAxis)w[1]+=(count>1?correction[1]:(power==1?1.:0.))*seed/scale;}}
                    for(const auto& item:sparse){port.support.push_back(item.first);for(int k=0;k<2;++k){require(std::isfinite(item.second[k]),"nonfinite native port weight");port.weights[k].push_back(item.second[k]);}}
                }
                datum.normalPort[component]=static_cast<int>(d->normalPorts.size());
                ++d->report.normalPorts;d->report.normalPortWeights+=2*port.support.size();
                d->report.weights+=2*port.support.size();d->normalPorts.push_back(std::move(port));
            }
            d->report.weights+=4*row.support.size();d->report.maxSupport=std::max(d->report.maxSupport,row.support.size());d->report.maxCoefficientNorm=std::max(d->report.maxCoefficientNorm,row.coefficientNorm);
            (mode==0?datum.base:datum.wall)[component]=static_cast<int>(d->rows.size());d->rows.push_back(std::move(row));
        }
    }
    for(std::size_t i=0;i<d->datums.size();++i)d->datumIndex.emplace(std::make_pair(d->datums[i].kind,d->datums[i].index),i);
    return CchRegularAdvectorGeometry(std::move(d));
}

CchRegularCoefficients CchRegularAdvectorGeometry::apply(const QuadtreeMesh& mesh,const CchNativeValues& values,const CchAdvectorBoundary* boundary) const
{
    fp();const auto& d=*data_;values.validateFor(d.geometry.layout(),mesh);if(boundary)boundary->validateFor(mesh,d.geometry);
    const auto& velocity=values.normalVelocity();std::vector<bool> wall(d.faces.size(),false);
    if(boundary)for(int u:boundary->velocityDofs()){require(velocity[u]==0,"nonzero prescribed wall velocity");wall[u]=true;}
    for(std::size_t u=0;u<d.faces.size();++u)if(d.faces[u].area==0)require(velocity[u]==0,"nonzero axis velocity");
    using Jet=std::array<std::array<double,4>,2>;
    std::vector<Jet> nodes(d.vertices.size());
    std::vector<std::array<std::array<double,2>,2>> mid(d.edges.size());
    std::vector<std::array<double,2>> centers(d.cells.size());
    for(const auto& datum:d.datums)for(int c=0;c<2;++c){bool prescribed=false;for(int f:datum.wallFaces[c])prescribed=prescribed||wall[f];
        const auto& row=d.rows[prescribed?datum.wall[c]:datum.base[c]];std::array<double,4> value{};
        for(int k=0;k<4;++k)for(std::size_t i=0;i<row.support.size();++i)value[k]+=row.weights[k][i]*velocity[row.support[i]];
        if(datum.kind==CchRegularRow::Kind::Vertex){
            const auto& port=d.normalPorts[datum.normalPort[c]];value[0]=0;value[2-c]=0;
            if(!prescribed)for(std::size_t i=0;i<port.support.size();++i){value[0]+=port.weights[0][i]*velocity[port.support[i]];value[2-c]+=port.weights[1][i]*velocity[port.support[i]];}
            nodes[datum.index][c]=value;}
        else if(datum.kind==CchRegularRow::Kind::Edge)mid[datum.index][c]={value[0],value[1+d.edges[datum.index].normalAxis]};
        else centers[datum.index][c]=value[0];}
    std::vector<unsigned char> nodeState(d.vertices.size()),edgeState(d.edges.size());
    std::function<void(int)> node,edge;
    std::function<Jet(int,double)> edgeJet=[&](int e,double s){edge(e);const auto& ge=d.edges[e];const int a=ge.vertices[0],b=ge.vertices[1],n=ge.normalAxis,t=1-n;
        node(a);node(b);const double h=coord(ge.upper-ge.lower,t);const auto bv=basis(s),bd=basis(s,true);Jet result{};
        for(int c=0;c<2;++c){const std::array<double,5> v{nodes[a][c][0],nodes[b][c][0],h*nodes[a][c][1+t],h*nodes[b][c][1+t],mid[e][c][0]};
            const std::array<double,5> dn{nodes[a][c][1+n],nodes[b][c][1+n],h*nodes[a][c][3],h*nodes[b][c][3],mid[e][c][1]};
            for(int k=0;k<5;++k){result[c][0]+=bv[k]*v[k];result[c][1+t]+=bd[k]*v[k]/h;result[c][1+n]+=bv[k]*dn[k];result[c][3]+=bd[k]*dn[k]/h;}}
        return result;};
    node=[&](int v){if(nodeState[v]==2)return;require(nodeState[v]==0,"cyclic active node constraint");nodeState[v]=1;const int e=d.vertexParent[v];
        if(e>=0){const auto& ge=d.edges[e];const int t=1-ge.normalAxis;nodes[v]=edgeJet(e,(coord(d.vertices[v],t)-coord(ge.lower,t))/coord(ge.upper-ge.lower,t));}nodeState[v]=2;};
    edge=[&](int e){if(edgeState[e]==2)return;require(edgeState[e]==0,"cyclic active edge constraint");edgeState[e]=1;const auto& ge=d.edges[e];
        if(ge.parent>=0){const auto& master=d.edges[ge.parent];const int t=1-ge.normalAxis;const auto jet=edgeJet(ge.parent,(coord((ge.lower+ge.upper)*.5,t)-coord(master.lower,t))/coord(master.upper-master.lower,t));
            for(int c=0;c<2;++c)mid[e][c]={jet[c][0],jet[c][1+ge.normalAxis]};}edgeState[e]=2;};
    for(int v:d.usedVertices)node(v);for(std::size_t e=0;e<d.edges.size();++e)edge(static_cast<int>(e));
    CchRegularCoefficients out;out.geometry_=data_;out.cellDofs.resize(d.cells.size());out.faceCorrections.resize(d.faces.size());
    constexpr int corners[4][2]={{0,0},{1,0},{1,1},{0,1}},mids[4][2]={{4,0},{4,1},{0,4},{1,4}};
    for(std::size_t p=0;p<d.cells.size();++p){const auto& cell=d.cells[p];const auto h=cell.upper-cell.lower;
        for(int c=0;c<2;++c){auto& raw=out.cellDofs[p][c];for(int k=0;k<4;++k){const int x=corners[k][0],y=corners[k][1];const auto& j=nodes[cell.vertices[k]][c];
                raw[5*y+x]=j[0];raw[5*y+x+2]=h.x*j[1];raw[5*(y+2)+x]=h.y*j[2];raw[5*(y+2)+x+2]=h.x*h.y*j[3];}
            for(int k=0;k<4;++k){const int x=mids[k][0],y=mids[k][1];const auto& z=mid[cell.edges[k]][c];raw[5*y+x]=z[0];if(k<2)raw[5*(y+2)+x]=z[1]*h.y;else raw[5*y+x+2]=z[1]*h.x;}
            raw[24]=centers[p][c];for(double v:raw)require(std::isfinite(v),"nonfinite Q4 datum");}}
    for(std::size_t u=0;u<d.faces.size();++u){const auto& f=d.faces[u];if(f.area==0)continue;const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure;
        double meanValue=0;for(int k=0;k<25;++k)meanValue+=f.meanWeights[k]*out.cellDofs[p][f.normalAxis][k];out.faceCorrections[u]=15./8*(velocity[u]-meanValue);
        require(std::isfinite(out.faceCorrections[u]),"nonfinite native mean correction");}
    return out;
}

namespace
{
void inside(const CchRegularCell& cell,Vector2D point)
{
    require(std::isfinite(point.x)&&std::isfinite(point.y)&&point.x>=cell.lower.x&&point.x<=cell.upper.x
        &&point.y>=cell.lower.y&&point.y<=cell.upper.y,"point outside regular cell");
}
double supportWidth(const CchRegularFace& face,int pressure)
{
    require(pressure==face.negativePressure||pressure==face.positivePressure,"missing regular face support");
    const auto width=face.supportWidth[pressure==face.positivePressure?1:0];
    require(std::isfinite(width)&&width>0,"invalid regular face support width");return width;
}
// Exact algebraic restriction q(origin+scale*t), never a sampled fit.
std::array<long double,5> compose(const double* power,long double origin,long double scale)
{
    std::array<long double,5> out{};
    for(int k=4;k>=0;--k){for(int j=4;j>0;--j)out[j]=origin*out[j]+scale*out[j-1];out[0]=origin*out[0]+power[k];}
    return out;
}
}
CchRegularCellStencil CchRegularAdvectorGeometry::cellStencil(const QuadtreeMesh& mesh,int pressure,const CchAdvectorBoundary* boundary) const
{
    fp();const auto& d=*data_;d.geometry.validateMesh(mesh);if(boundary)boundary->validateFor(mesh,d.geometry);
    require(pressure>=0&&static_cast<std::size_t>(pressure)<d.cells.size(),"stencil pressure index out of range");
    using Sparse=std::map<int,double>;
    using Jet=std::array<std::array<Sparse,4>,2>;
    using Mid=std::array<std::array<Sparse,2>,2>;
    using Raw=std::array<std::array<Sparse,25>,2>;
    const auto add=[](Sparse& into,const Sparse& from,double factor){if(factor!=0)for(const auto& w:from)into[w.first]+=factor*w.second;};
    std::set<int> wall;if(boundary)wall.insert(boundary->velocityDofs().begin(),boundary->velocityDofs().end());
    const auto allowed=[&](int u){return d.faces[u].area>0&&wall.count(u)==0;};
    const auto datumJet=[&](CchRegularRow::Kind kind,int index){Jet jet{};const auto& datum=d.datums[d.datumIndex.at({kind,index})];
        for(int c=0;c<2;++c){bool prescribed=false;for(int f:datum.wallFaces[c])prescribed=prescribed||wall.count(f)!=0;
            const auto& row=d.rows[prescribed?datum.wall[c]:datum.base[c]];
            for(int k=0;k<4;++k)for(std::size_t i=0;i<row.support.size();++i)if(allowed(row.support[i]))jet[c][k][row.support[i]]+=row.weights[k][i];
            if(kind==CchRegularRow::Kind::Vertex){const auto& port=d.normalPorts[datum.normalPort[c]];jet[c][0].clear();jet[c][2-c].clear();
                if(!prescribed)for(std::size_t i=0;i<port.support.size();++i)if(allowed(port.support[i])){
                    jet[c][0][port.support[i]]+=port.weights[0][i];jet[c][2-c][port.support[i]]+=port.weights[1][i];}}
        }return jet;};
    std::map<int,Jet> nodes;std::map<int,Mid> mids;std::map<int,Raw> rawCells;
    std::set<int> pendingNodes,pendingEdges;
    std::function<const Jet&(int)> node;std::function<const Mid&(int)> edge;
    const auto edgeJet=[&](int e,double s){const auto& mid=edge(e);const auto& ge=d.edges[e];
        const auto& a=node(ge.vertices[0]);const auto& b=node(ge.vertices[1]);const int n=ge.normalAxis,t=1-n;
        const double h=coord(ge.upper-ge.lower,t);const auto bv=basis(s),bd=basis(s,true);Jet result{};
        for(int c=0;c<2;++c){const Sparse* values[5]={&a[c][0],&b[c][0],&a[c][1+t],&b[c][1+t],&mid[c][0]};
            const Sparse* normal[5]={&a[c][1+n],&b[c][1+n],&a[c][3],&b[c][3],&mid[c][1]};
            for(int k=0;k<5;++k){const double scale=(k==2||k==3)?h:1;
                add(result[c][0],*values[k],bv[k]*scale);add(result[c][1+t],*values[k],bd[k]*scale/h);
                add(result[c][1+n],*normal[k],bv[k]*scale);add(result[c][3],*normal[k],bd[k]*scale/h);}}
        return result;};
    node=[&](int v)->const Jet&{const auto found=nodes.find(v);if(found!=nodes.end())return found->second;
        require(pendingNodes.insert(v).second,"cyclic local node constraint");Jet jet{};const int e=d.vertexParent[v];
        if(e<0)jet=datumJet(CchRegularRow::Kind::Vertex,v);
        else {const auto& ge=d.edges[e];const int t=1-ge.normalAxis;jet=edgeJet(e,(coord(d.vertices[v],t)-coord(ge.lower,t))/coord(ge.upper-ge.lower,t));}
        pendingNodes.erase(v);return nodes.emplace(v,std::move(jet)).first->second;};
    edge=[&](int e)->const Mid&{const auto found=mids.find(e);if(found!=mids.end())return found->second;
        require(pendingEdges.insert(e).second,"cyclic local edge constraint");const auto& ge=d.edges[e];Mid mid{};
        if(ge.parent<0){const auto jet=datumJet(CchRegularRow::Kind::Edge,e);for(int c=0;c<2;++c)mid[c]={jet[c][0],jet[c][1+ge.normalAxis]};}
        else {const auto& master=d.edges[ge.parent];const int t=1-ge.normalAxis;
            const auto jet=edgeJet(ge.parent,(coord((ge.lower+ge.upper)*.5,t)-coord(master.lower,t))/coord(master.upper-master.lower,t));
            for(int c=0;c<2;++c)mid[c]={jet[c][0],jet[c][1+ge.normalAxis]};}
        pendingEdges.erase(e);return mids.emplace(e,std::move(mid)).first->second;};
    const auto rawCell=[&](int p)->const Raw&{const auto found=rawCells.find(p);if(found!=rawCells.end())return found->second;
        const auto& cell=d.cells[p];const auto h=cell.upper-cell.lower;Raw raw{};
        constexpr int corners[4][2]={{0,0},{1,0},{1,1},{0,1}},positions[4][2]={{4,0},{4,1},{0,4},{1,4}};
        for(int c=0;c<2;++c){for(int k=0;k<4;++k){const int x=corners[k][0],y=corners[k][1];const auto& j=node(cell.vertices[k])[c];
                add(raw[c][5*y+x],j[0],1);add(raw[c][5*y+x+2],j[1],h.x);add(raw[c][5*(y+2)+x],j[2],h.y);add(raw[c][5*(y+2)+x+2],j[3],h.x*h.y);}
            for(int k=0;k<4;++k){const int x=positions[k][0],y=positions[k][1];const auto& z=edge(cell.edges[k])[c];
                add(raw[c][5*y+x],z[0],1);if(k<2)add(raw[c][5*(y+2)+x],z[1],h.y);else add(raw[c][5*y+x+2],z[1],h.x);}
            raw[c][24]=datumJet(CchRegularRow::Kind::Center,p)[c][0];}
        return rawCells.emplace(p,std::move(raw)).first->second;};
    const auto pack=[](const Sparse& values){std::vector<CchRegularScalarWeight> result;for(const auto& w:values)if(w.second!=0){require(std::isfinite(w.second),"nonfinite local source stencil");result.push_back({w.first,w.second});}return result;};
    CchRegularCellStencil out;out.pressure=pressure;out.geometry_=data_;const auto& raw=rawCell(pressure);
    for(int c=0;c<2;++c)for(int k=0;k<25;++k)out.cellDofs[c][k]=pack(raw[c][k]);
    for(int u:d.cells[pressure].faces){const auto& f=d.faces[u];Sparse correction;
        if(f.area>0){const int p=f.negativePressure>=0?f.negativePressure:f.positivePressure;const auto& owner=rawCell(p);
            if(allowed(u))correction[u]=15./8;for(int k=0;k<25;++k)add(correction,owner[f.normalAxis][k],-15./8*f.meanWeights[k]);}
        out.faceCorrections.push_back(pack(correction));}
    return out;
}
std::vector<CchRegularJetWeight> CchRegularAdvectorGeometry::cellJetWeights(const CchRegularCellStencil& stencil,Vector2D point) const
{
    fp();require(stencil.geometry_==data_,"different local stencil geometry identity");
    require(stencil.pressure>=0&&static_cast<std::size_t>(stencil.pressure)<cells().size(),"stencil pressure index out of range");
    const auto& cell=cells()[stencil.pressure];inside(cell,point);require(stencil.faceCorrections.size()==cell.faces.size(),"different local stencil face inventory");
    const auto h=cell.upper-cell.lower;const auto bx=basis((point.x-cell.lower.x)/h.x),by=basis((point.y-cell.lower.y)/h.y);
    const auto dx=basis((point.x-cell.lower.x)/h.x,true),dy=basis((point.y-cell.lower.y)/h.y,true);
    std::map<int,std::array<long double,6>> weights;
    const auto add=[&](const auto& list,int c,double value,double x,double y){for(const auto& w:list){auto& jet=weights[w.velocityDof];
        jet[c]+=static_cast<long double>(w.weight)*value;jet[2+2*c]+=static_cast<long double>(w.weight)*x;jet[3+2*c]+=static_cast<long double>(w.weight)*y;}};
    for(int c=0;c<2;++c)for(int y=0;y<5;++y)for(int x=0;x<5;++x)add(stencil.cellDofs[c][5*y+x],c,bx[x]*by[y],dx[x]*by[y]/h.x,bx[x]*dy[y]/h.y);
    for(std::size_t i=0;i<cell.faces.size();++i){const auto& face=faces()[cell.faces[i]];const int n=face.normalAxis,t=1-n;
        if(face.area==0)continue;const auto width=supportWidth(face,stencil.pressure),length=coord(face.upper-face.lower,t);
        const auto s=(coord(point,t)-coord(face.lower,t))/length,xi=(coord(point,n)-coord(face.lower,n))/width;if(s<0||s>1||xi<-1||xi>1)continue;
        const auto g=16*s*s*(1-s)*(1-s),b=(1-xi*xi)*(1-xi*xi),gt=32*s*(1-s)*(1-2*s)/length,bn=-4*xi*(1-xi*xi)/width;
        add(stencil.faceCorrections[i],n,b*g,n==0?bn*g:b*gt,n==0?b*gt:bn*g);}
    std::vector<CchRegularJetWeight> out;for(const auto& w:weights){CchRegularJetWeight row;row.velocityDof=w.first;bool nonzero=false;
        for(int k=0;k<6;++k){row.jet[k]=static_cast<double>(w.second[k]);require(std::isfinite(row.jet[k]),"nonfinite local jet weight");nonzero=nonzero||row.jet[k]!=0;}
        if(nonzero)out.push_back(row);}return out;
}
void CchRegularAdvectorGeometry::validateCoefficients(const CchRegularCoefficients& coefficients,int pressure) const
{
    fp();require(coefficients.geometry_==data_,"different regular coefficient geometry identity");
    require(coefficients.cellDofs.size()==cells().size()&&coefficients.faceCorrections.size()==faces().size(),"different regular coefficient inventory");
    require(pressure>=0&&static_cast<std::size_t>(pressure)<cells().size(),"regular pressure index out of range");
}
CchAdvectorCellJet CchRegularAdvectorGeometry::cellJet(const CchRegularCoefficients& coefficients,int pressure,Vector2D point) const
{
    validateCoefficients(coefficients,pressure);const auto& cell=cells()[pressure];inside(cell,point);
    const auto h=cell.upper-cell.lower;const auto bx=basis((point.x-cell.lower.x)/h.x),by=basis((point.y-cell.lower.y)/h.y);
    const auto dx=basis((point.x-cell.lower.x)/h.x,true),dy=basis((point.y-cell.lower.y)/h.y,true);
    std::array<long double,2> value{},jx{},jy{};
    for(int component=0;component<2;++component)for(int y=0;y<5;++y)for(int x=0;x<5;++x){const auto datum=coefficients.cellDofs[pressure][component][5*y+x];
        require(std::isfinite(datum),"nonfinite regular Q4 datum");value[component]+=datum*bx[x]*by[y];
        jx[component]+=datum*dx[x]*by[y]/h.x;jy[component]+=datum*bx[x]*dy[y]/h.y;}
    for(int u:cell.faces){const auto& face=faces()[u];const int n=face.normalAxis,t=1-n;const auto delta=coefficients.faceCorrections[u];
        require(std::isfinite(delta),"nonfinite regular face correction");if(delta==0)continue;
        const auto width=supportWidth(face,pressure),length=coord(face.upper-face.lower,t);
        const auto s=(coord(point,t)-coord(face.lower,t))/length,xi=(coord(point,n)-coord(face.lower,n))/width;
        if(s<0||s>1||xi<-1||xi>1)continue;
        const auto g=16*s*s*(1-s)*(1-s),b=(1-xi*xi)*(1-xi*xi);
        const auto gt=32*s*(1-s)*(1-2*s)/length,bn=-4*xi*(1-xi*xi)/width;
        value[n]+=delta*b*g;(n==0?jx:jy)[n]+=delta*bn*g;(t==0?jx:jy)[n]+=delta*b*gt;
    }
    Vector2D lower,upper;
    for(int axis=0;axis<2;++axis){const auto& knots=cell.supportKnots[axis];const auto coordinate=coord(point,axis);
        require(knots.size()>=2&&knots.front()==coord(cell.lower,axis)&&knots.back()==coord(cell.upper,axis),"incomplete regular support partition");
        auto index=static_cast<std::size_t>(std::upper_bound(knots.begin(),knots.end(),coordinate)-knots.begin());if(index==knots.size())--index;
        (axis==0?lower.x:lower.y)=knots[index-1];(axis==0?upper.x:upper.y)=knots[index];
    }
    for(int c=0;c<2;++c)require(std::isfinite(static_cast<double>(value[c]))&&std::isfinite(static_cast<double>(jx[c]))&&std::isfinite(static_cast<double>(jy[c])),"nonfinite regular jet");
    return {{{static_cast<double>(value[0]),static_cast<double>(value[1])},
        {static_cast<double>(jx[0]),static_cast<double>(jy[0]),static_cast<double>(jx[1]),static_cast<double>(jy[1])}},lower,upper};
}
CchNormalPolynomial CchRegularAdvectorGeometry::normalPolynomial(const CchRegularCoefficients& coefficients,int pressure,
    int normal,Vector2D begin,Vector2D end) const
{
    validateCoefficients(coefficients,pressure);require(normal==0||normal==1,"invalid regular normal axis");
    const auto& cell=cells()[pressure];inside(cell,begin);inside(cell,end);const int tangent=1-normal;
    require(coord(begin,normal)==coord(end,normal)&&coord(end,tangent)>coord(begin,tangent),"invalid canonical regular normal interval");
    for(double knot:cell.supportKnots[tangent])require(!(coord(begin,tangent)<knot&&knot<coord(end,tangent)),"normal interval crosses regular support knot");
    CchNormalPolynomial out;out.length=static_cast<long double>(coord(end,tangent))-coord(begin,tangent);
    out.radialMeasure=momentumGeometry().layout().domain().axisymmetric();out.radialTangent=out.radialMeasure&&normal==1;out.radiusBegin=begin.x;out.radiusEnd=end.x;
    const auto h=cell.upper-cell.lower;const auto normalBasis=basis((coord(begin,normal)-coord(cell.lower,normal))/coord(h,normal));
    for(int reverse=0;reverse<2;++reverse){const auto origin=reverse?end:begin;auto& result=reverse?out.endCoefficients:out.beginCoefficients;
        const long double sign=reverse?-1:1,z=(static_cast<long double>(coord(origin,tangent))-coord(cell.lower,tangent))/coord(h,tangent);
        for(int j=0;j<5;++j){const auto restriction=compose(basisPower[j],z,sign/coord(h,tangent));long double weight=0;
            for(int i=0;i<5;++i){const auto datum=coefficients.cellDofs[pressure][normal][normal==0?5*j+i:5*i+j];
                require(std::isfinite(datum),"nonfinite regular Q4 datum");weight+=datum*normalBasis[i];}
            for(int k=0;k<5;++k)result[k]+=weight*restriction[k];}
        for(int u:cell.faces){const auto& face=faces()[u];if(face.normalAxis!=normal)continue;
            const auto delta=coefficients.faceCorrections[u];require(std::isfinite(delta),"nonfinite regular face correction");if(delta==0)continue;
            const auto width=supportWidth(face,pressure),length=coord(face.upper-face.lower,tangent);
            const auto xi=(coord(origin,normal)-coord(face.lower,normal))/width;
            const auto mid=(static_cast<long double>(coord(begin,tangent))+coord(end,tangent))/2;
            if(xi<-1||xi>1||mid<coord(face.lower,tangent)||mid>coord(face.upper,tangent))continue;
            const auto restriction=compose(basisPower[4],(static_cast<long double>(coord(origin,tangent))-coord(face.lower,tangent))/length,sign/length);
            const long double weight=delta*(1-xi*xi)*(1-xi*xi);for(int k=0;k<5;++k)result[k]+=weight*restriction[k];
        }
        for(auto value:result)require(std::isfinite(value)&&std::isfinite(static_cast<double>(value)),"nonfinite regular normal polynomial");
    }
    return out;
}
}
