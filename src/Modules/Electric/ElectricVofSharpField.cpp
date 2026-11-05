#include "Modules/Electric/ElectricVofSharpField.h"
#include "Modules/Electric/ElectricSurfaceChargeTransport.h"
#include "Solvers/Poisson/ElectricLocalOperator.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <Eigen/Cholesky>
#include <Eigen/SparseCore>
#ifdef EIGEN_SUPERLU_SUPPORT
#include <Eigen/SuperLUSupport>
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace LAMP::Electric
{
namespace
{
namespace E=LAMP::ElectricDiscretization;
namespace FV=LAMP::Geometry::Detail::FiniteVolume;
using G=Interface::VofSmoothGeometry;
using Clock=LAMP::TimeLayer;
using Timer=std::chrono::steady_clock;
using Triple=std::array<double,3>;
using ExactTriple=std::array<FV::Dyadic,3>;
constexpr std::array<double,16> gx{{
-.98940093499164993260,-.94457502307323257608,-.86563120238783174388,-.75540440835500303390,
-.61787624440264374845,-.45801677765722738634,-.28160355077925891323,-.09501250983763744019,
 .09501250983763744019, .28160355077925891323, .45801677765722738634, .61787624440264374845,
 .75540440835500303390, .86563120238783174388, .94457502307323257608, .98940093499164993260}};
constexpr std::array<double,16> gw{{
.02715245941175409485,.06225352393864789286,.09515851168249278481,.12462897125553387205,
.14959598881657673208,.16915651939500253819,.18260341504492358887,.18945061045506849629,
.18945061045506849629,.18260341504492358887,.16915651939500253819,.14959598881657673208,
.12462897125553387205,.09515851168249278481,.06225352393864789286,.02715245941175409485}};
void require(bool ok,const char* reason)
{if(!ok)throw std::invalid_argument(std::string("[ElectricVofSharpField] ")+reason);}
bool sameClock(Clock a,Clock b)
{return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;}
double elapsed(Timer::time_point t){return std::chrono::duration<double>(Timer::now()-t).count();}
Vector2D asVector(Geometry::PlanarCellCut::Point p)
{return {static_cast<double>(p.x),static_cast<double>(p.y)};}
Vector2D relativePoint(Geometry::PlanarCellCut::Point p,Vector2D origin)
{
    if(origin.x==0&&origin.y==0)return asVector(p);
    return {static_cast<double>(p.x-static_cast<long double>(origin.x)),
        static_cast<double>(p.y-static_cast<long double>(origin.y))};
}
double relative(double error,double scale)
{require(std::isfinite(error)&&std::isfinite(scale)&&scale>=0,"nonfinite residual");return error==0?0:(scale>0?error/scale:std::numeric_limits<double>::infinity());}
Triple triple(const Eigen::Vector3d& x){return {{x[0],x[1],x[2]}};}
Eigen::Vector3d vector(const Triple& x){return {x[0],x[1],x[2]};}
double exposed(const FV::Dyadic& value)
{const double rounded=FV::roundDouble(value).value;require(std::isfinite(rounded),"unrepresentable physical flux or charge");return rounded;}
FV::Dyadic exactRowDot(const Eigen::MatrixXd& matrix,Eigen::Index row,const Eigen::VectorXd& values)
{FV::Dyadic result;for(Eigen::Index j=0;j<values.size();++j)result.add(FV::Dyadic::product({matrix(row,j),values[j]}));return result;}
struct Face
{
    int nativeRow=-1,arc=-1,boundary=-1,offset=-1;
    unsigned phase=0;
    std::vector<E::FaceQuadraturePoint> quadrature; // canonical n, ordinary ds.
    double length=0;
    bool prescribed=false;
    Eigen::Matrix3d mass=Eigen::Matrix3d::Zero();
    Eigen::Vector3d trace=Eigen::Vector3d::Zero();
    std::array<Eigen::Vector3d,2> normalElectricCoefficient{{Eigen::Vector3d::Zero(),Eigen::Vector3d::Zero()}};
    Eigen::Vector3d chargeDensityCoefficient=Eigen::Vector3d::Zero();
    std::array<int,2> elements{{-1,-1}}; // Gamma liquid/gas or native neg/pos.
};
std::vector<E::FaceQuadraturePoint> curveQuadrature(const Interface::VofSmoothArc& arc)
{
    std::vector<E::FaceQuadraturePoint> q;q.reserve(gx.size());
    for(std::size_t j=0;j<gx.size();++j){const double t=(1+gx[j])/2;
        const auto d=arc.derivative(t);const double speed=static_cast<double>(std::hypot(d.x,d.y));
        require(std::isfinite(speed)&&speed>0,"degenerate curve derivative");
        q.push_back({asVector(arc.sample(t)),asVector(arc.normal(t)),gw[j]*speed/2,gx[j]});}
    return q;
}
void finishFace(Face& f,bool rz)
{
    for(const auto& q:f.quadrature){const auto p=E::faceBasis(q.coordinate);
        f.length+=q.weight;f.mass.noalias()+=q.weight*E::physicalWeight(q.point,rz)*p*p.transpose();}
    Eigen::LLT<Eigen::Matrix3d> factor(f.mass);
    require(std::isfinite(f.length)&&f.length>0&&factor.info()==Eigen::Success,
        "face P2 mass lacks positive rank; unsupported small cut");
}
unsigned sourcePhase(const G& g,int cell,Vector2D p)
{
    return g.phaseAt(cell,{static_cast<long double>(p.x),static_cast<long double>(p.y)});
}
struct Element
{
    int cell=-1;unsigned phase=0;
    Vector2D center,coordinateOrigin{0,0};double scale=0;
    E::LocalGeometry geometry;
    std::vector<int> faces;
    std::shared_ptr<const E::LocalOperator> op;
    Eigen::MatrixXd combined,cellFromFaces,condensed;
    E::CellVector cellPolynomial=E::CellVector::Zero();
    E::PotentialVector potential=E::PotentialVector::Zero();
    Eigen::VectorXd localValues,unitResidual;
};
}
struct ElectricVofSurfaceCharge::Data
{
    G geometry;Clock clock;std::vector<Triple> moments;std::vector<double> integrated;
    Data(const G& g,Clock t,std::vector<Triple> q):geometry(g),clock(t),moments(std::move(q))
    {for(const auto& m:moments)integrated.push_back(m[0]);}
};
ElectricVofSurfaceCharge::ElectricVofSurfaceCharge(std::shared_ptr<const Data> d):data_(std::move(d)){}
ElectricVofSurfaceCharge ElectricVofSurfaceCharge::fromPrepared(const QuadtreeMesh& mesh,
    const G& g,Clock t,std::vector<Triple> q)
{
    g.validateFor(mesh,g.source());t.time();require(q.size()==g.arcs().size(),"surface charge arc count mismatch");
    for(const auto& m:q)for(double x:m)require(std::isfinite(x),"nonfinite surface moment");
    return ElectricVofSurfaceCharge(std::make_shared<Data>(g,t,std::move(q)));
}
ElectricVofSurfaceCharge ElectricVofSurfaceCharge::fromMoments(const QuadtreeMesh& mesh,
    const G& g,Clock t,const std::vector<Triple>& q)
{
    require(sameClock(t,g.sourceTime())&&t.position!=Clock::Position::Midpoint,
        "initial surface moments require the genuine accepted geometry clock");
    return fromPrepared(mesh,g,t,q);
}
ElectricVofSurfaceCharge ElectricVofSurfaceCharge::fromDensity(const QuadtreeMesh& mesh,
    const G& g,const std::function<double(std::size_t,double)>& density)
{
    require(bool(density),"missing surface density provider");std::vector<Triple> q(g.arcs().size());
    for(std::size_t a=0;a<q.size();++a){Eigen::Vector3d m=Eigen::Vector3d::Zero();
        for(const auto& node:curveQuadrature(g.arcs()[a])){const double value=density(a,(1+node.coordinate)/2);
            require(std::isfinite(value),"nonfinite surface density");
            m.noalias()+=node.weight*E::physicalWeight(node.point,g.source().domain().axisymmetric())*value*E::faceBasis(node.coordinate);}
        q[a]=triple(m);}
    return fromMoments(mesh,g,g.sourceTime(),q);
}
const G& ElectricVofSurfaceCharge::geometry() const{require(bool(data_),"consumed surface owner");return data_->geometry;}
Clock ElectricVofSurfaceCharge::timeLayer() const{require(bool(data_),"consumed surface owner");return data_->clock;}
const std::vector<Triple>& ElectricVofSurfaceCharge::moments() const{require(bool(data_),"consumed surface owner");return data_->moments;}
const std::vector<double>& ElectricVofSurfaceCharge::integratedCharge() const{require(bool(data_),"consumed surface owner");return data_->integrated;}
bool ElectricVofSurfaceCharge::sameInstance(const ElectricVofSurfaceCharge& q) const noexcept{return data_&&data_==q.data_;}
void ElectricVofSurfaceCharge::validateFor(const QuadtreeMesh& mesh,const G& g,Clock t) const
{require(bool(data_)&&data_->geometry.sameInstance(g)&&sameClock(t,data_->clock),"foreign surface geometry/clock");g.validateFor(mesh,g.source());t.time();}
struct ElectricVofSharpField::Data
{
    G geometry;ElectricVofSurfaceCharge source,candidate,input;Clock clock;double dt=0,potentialReference=0;
    std::shared_ptr<const ElectricSurfaceChargeTransport> advection;
    ElectricVofSharpProblem problem;ElectricVofSharpReport report;
    std::vector<Face> faces;std::vector<Element> elements;
    std::vector<std::array<int,2>> elementByCell;
    std::vector<int> interfaceFace;
    Data(const G& g,const ElectricVofSurfaceCharge& q,Clock t,const ElectricVofSharpProblem& p)
        :geometry(g),source(q),candidate(q),input(q),clock(t),problem(p){}
};
ElectricVofSharpField::ElectricVofSharpField(std::shared_ptr<const Data> d):data_(std::move(d)){}
ElectricVofSharpField ElectricVofSharpField::prepareElectrostatic(const QuadtreeMesh& m,const G& g,
    const ElectricVofSurfaceCharge& q,const ElectricVofSharpProblem& p,const ElectricVofSharpControls& c)
{return prepare(m,g,q,q.timeLayer(),p,c,true);}
ElectricVofSharpField ElectricVofSharpField::prepareOhmic(const QuadtreeMesh& m,const G& g,
    const ElectricVofSurfaceCharge& q,Clock t,const ElectricVofSharpProblem& p,const ElectricVofSharpControls& c)
{return prepare(m,g,q,t,p,c,false);}
ElectricVofSharpField ElectricVofSharpField::prepareOhmic(const QuadtreeMesh& m,
    const ElectricSurfaceChargeTransport& transport,const ElectricVofSharpProblem& p,const ElectricVofSharpControls& c)
{
    require(transport.interfaceTransport()!=nullptr,"Ohmic advection requires a genuine advancing Interface receipt");
    const auto& input=transport.chargeCandidate();
    return prepare(m,input.geometry(),input,input.timeLayer(),p,c,false,&transport);
}
ElectricVofSharpField ElectricVofSharpField::prepare(const QuadtreeMesh& mesh,const G& g,
    const ElectricVofSurfaceCharge& q,Clock clock,const ElectricVofSharpProblem& p,
    const ElectricVofSharpControls& controls,bool snapshot,const ElectricSurfaceChargeTransport* advection)
{
    const auto started=Timer::now();g.validateFor(mesh,g.source());q.validateFor(mesh,g,q.timeLayer());clock.time();
    const auto& domain=g.source().domain();const bool rz=domain.axisymmetric();
    require(domain.uniformLeaves()&&g.cells().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "E1 requires actual complete same-level smooth geometry");
    require(std::isfinite(controls.relativeLinearTolerance)&&controls.relativeLinearTolerance>0&&controls.relativeLinearTolerance<=1e-11
        &&std::isfinite(controls.relativeConservationTolerance)&&controls.relativeConservationTolerance>0
        &&controls.relativeConservationTolerance<=1e-10,"invalid sharp field tolerance");
    for(unsigned phase=0;phase<2;++phase)require(std::isfinite(p.epsilon[phase])&&p.epsilon[phase]>0
        &&std::isfinite(p.conductivity[phase])&&p.conductivity[phase]>=0,"invalid sharp phase material");
    for(const auto& b:p.boundaries)require((b.type==ElectricVofSharpBoundary::Type::Dirichlet&&bool(b.value))
        ||b.type==ElectricVofSharpBoundary::Type::HomogeneousNeumann,"invalid explicit boundary provider");
    auto d=std::make_shared<Data>(g,q,clock,p);
    if(advection){require(q.sameInstance(advection->chargeCandidate())&&sameClock(clock,g.sourceTime()),
            "advected implicit input lost its actual target geometry/clock");
        d->advection=std::make_shared<ElectricSurfaceChargeTransport>(*advection);d->source=advection->sourceCharge();
        d->source.validateFor(mesh,d->source.geometry(),d->source.timeLayer());}
    if(snapshot)require(sameClock(clock,q.timeLayer()),"snapshot clock mismatch");
    else{require(clock.position==Clock::Position::Midpoint||clock.position==Clock::Position::End,"Ohmic target must be an advancing stage");
        require(d->source.timeLayer().step<std::numeric_limits<std::int64_t>::max()
            &&clock.step==d->source.timeLayer().step+1
            &&d->source.timeLayer().time()==clock.begin&&d->source.timeLayer().position!=Clock::Position::Midpoint,
            "frozen Ohmic source must be accepted at interval begin");d->dt=clock.time()-clock.begin;}
    d->elementByCell.assign(g.cells().size(),{{-1,-1}});
    for(const auto& cell:g.cells())for(unsigned phase=0;phase<2;++phase)if(cell.inventory[phase]>0){
        Element e;e.cell=cell.nativeCell;e.phase=phase;const auto& r=cell.rectangle;
        if(rz)e.geometry.center={r.lower.x+(r.upper.x-r.lower.x)/2,r.lower.y+(r.upper.y-r.lower.y)/2};
        else{
            e.coordinateOrigin=r.lower;e.geometry.coordinateOrigin=e.coordinateOrigin;
            e.geometry.center={static_cast<double>((static_cast<long double>(r.upper.x)-r.lower.x)/2),
                static_cast<double>((static_cast<long double>(r.upper.y)-r.lower.y)/2)};
        }
        e.geometry.scale=std::min(r.upper.x-r.lower.x,r.upper.y-r.lower.y);e.geometry.axisymmetric=rz;
        const auto nativeCenter=e.geometry.center;long double area=0,firstX=0,firstY=0;
        for(const auto& v:cell.phaseQuadrature[phase]){const auto point=relativePoint(v.point,e.coordinateOrigin);const double metric=E::physicalWeight(point,rz);
            require(metric>0&&std::isfinite(metric),"volume node lies on axis; ordinary area unavailable");
            const double weight=static_cast<double>(v.weight/metric);
            require(std::isfinite(weight)&&weight>0,"unsupported phase volume weight");e.geometry.volume.push_back({point,weight});
            area+=weight;firstX+=weight*(static_cast<long double>(point.x)-nativeCenter.x);
            firstY+=weight*(static_cast<long double>(point.y)-nativeCenter.y);}
        require(!e.geometry.volume.empty(),"positive phase inventory has no actual quadrature");
        require(std::isfinite(area)&&area>0,"nonfinite phase ordinary area");
        // The centroid stays in the declared element frame. In XY, no world
        // addition/RN precedes the basis, extent or phase conditioning.
        e.geometry.center={static_cast<double>(nativeCenter.x+firstX/area),static_cast<double>(nativeCenter.y+firstY/area)};
        e.center=e.geometry.center;e.scale=e.geometry.scale;
        d->elementByCell.at(static_cast<std::size_t>(e.cell))[phase]=static_cast<int>(d->elements.size());d->elements.push_back(std::move(e));}
    std::map<int,std::size_t> nodeByFace;
    for(std::size_t n=0;n<g.nodes().size();++n)require(nodeByFace.emplace(g.nodes()[n].nativeRow,n).second,
        "multiple smooth crossings per native face are unsupported");
    auto append=[&](Face f){finishFace(f,rz);const int index=static_cast<int>(d->faces.size());
        for(unsigned side=0;side<2;++side)if(f.elements[side]>=0){auto& e=d->elements[static_cast<std::size_t>(f.elements[side])];
            E::FaceQuadrature quadrature;quadrature.length=f.length;quadrature.nodes=f.quadrature;
            // Native canonical negative->positive; Gamma canonical liquid->gas.
            const double sign=side==0?1:-1;for(auto& n:quadrature.nodes){
                n.point=relativePoint({n.point.x,n.point.y},e.coordinateOrigin);
                n.outwardNormal.x*=sign;n.outwardNormal.y*=sign;}
            e.faces.push_back(index);e.geometry.faces.push_back(std::move(quadrature));}
        d->faces.push_back(std::move(f));};
    for(int row=0;row<domain.nativeFaceCount();++row){const auto nf=domain.face(mesh,row);
        if(nf.area==0){require(rz&&nf.center.x==0,"zero measure nonaxis face");
            require(p.boundaries[0].type==ElectricVofSharpBoundary::Type::HomogeneousNeumann,"RZ axis requires homogeneous Neumann");continue;}
        const double low=(nf.horizontal?nf.center.x:nf.center.y)-nf.tangentSpacing/2;
        const double high=low+nf.tangentSpacing;std::vector<double> cuts{low,high};
        const auto found=nodeByFace.find(row);
        if(found!=nodeByFace.end()){const auto point=g.nodes()[found->second].point;
            const double cut=static_cast<double>(nf.horizontal?point.x:point.y);
            require(cut>low&&cut<high,"smooth crossing must be strictly inside native face");cuts.insert(cuts.begin()+1,cut);}
        for(std::size_t range=1;range<cuts.size();++range){const double a=cuts[range-1],b=cuts[range];
            const Vector2D mid=nf.horizontal?Vector2D{(a+b)/2,nf.center.y}:Vector2D{nf.center.x,(a+b)/2};
            unsigned phase=0;
            if(found!=nodeByFace.end()){const auto& node=g.nodes()[found->second];
                const long double value=node.tangent.y*(static_cast<long double>(mid.x)-node.point.x)
                    -node.tangent.x*(static_cast<long double>(mid.y)-node.point.y);require(value!=0,"tangent phase split is unresolved");phase=value<0?0:1;}
            else{const int cell=nf.negativeCell>=0?nf.negativeCell:nf.positiveCell;phase=sourcePhase(g,cell,mid);
                if(nf.negativeCell>=0&&nf.positiveCell>=0)require(phase==sourcePhase(g,nf.positiveCell,mid),
                    "native face phase classification disagrees without a smooth crossing");}
            Face f;f.nativeRow=row;f.phase=phase;
            if(nf.negativeCell<0||nf.positiveCell<0)f.boundary=nf.horizontal?(nf.negativeCell<0?1:3):(nf.negativeCell<0?0:2);
            for(unsigned side=0;side<2;++side){const int cell=side==0?nf.negativeCell:nf.positiveCell;
                if(cell>=0){f.elements[side]=d->elementByCell.at(static_cast<std::size_t>(cell))[phase];
                    require(f.elements[side]>=0,"native smooth phase face lacks its phase element");}}
            const Vector2D normal=nf.horizontal?Vector2D{0,1}:Vector2D{1,0};
            for(std::size_t j=0;j<gx.size();++j){const double x=a+(b-a)*(1+gx[j])/2;
                f.quadrature.push_back({nf.horizontal?Vector2D{x,nf.center.y}:Vector2D{nf.center.x,x},normal,gw[j]*(b-a)/2,gx[j]});}
            append(std::move(f));}}
    d->interfaceFace.resize(g.arcs().size());
    for(std::size_t a=0;a<g.arcs().size();++a){const auto& arc=g.arcs()[a];Face f;f.arc=static_cast<int>(a);
        f.elements=d->elementByCell.at(static_cast<std::size_t>(arc.nativeCell));
        require(f.elements[0]>=0&&f.elements[1]>=0,"smooth interface lacks two positive phases");f.quadrature=curveQuadrature(arc);
        d->interfaceFace[a]=static_cast<int>(d->faces.size());append(std::move(f));}
    int free=0;bool anchored=false;
    for(auto& f:d->faces){f.prescribed=f.boundary>=0&&p.boundaries[static_cast<std::size_t>(f.boundary)].type==ElectricVofSharpBoundary::Type::Dirichlet;
        if(f.prescribed){if(!anchored){d->potentialReference=p.boundaries[static_cast<std::size_t>(f.boundary)].value(f.quadrature.front().point);
                require(std::isfinite(d->potentialReference),"nonfinite potential reference");}
            anchored=true;Eigen::Vector3d rhs=Eigen::Vector3d::Zero();
            // A fixed additive voltage gauge is removed before assembly.
            // Electric fields and physical conservation scales are gauge invariant.
            for(const auto& n:f.quadrature){const double value=p.boundaries[static_cast<std::size_t>(f.boundary)].value(n.point)-d->potentialReference;
                require(std::isfinite(value),"nonfinite prescribed potential");rhs.noalias()+=n.weight*E::physicalWeight(n.point,rz)*value*E::faceBasis(n.coordinate);}
            f.trace=f.mass.llt().solve(rhs);require(f.trace.allFinite(),"invalid Dirichlet projection");}
        else{require(free<=std::numeric_limits<int>::max()-3,"skeleton size overflow");f.offset=free;free+=3;}}
    require(anchored,"sharp field requires an explicit physical Dirichlet anchor");
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(free);std::vector<Eigen::Triplet<double>> entries;
    for(const auto& f:d->faces)if(f.arc>=0)rhs.segment<3>(f.offset)+=vector(q.moments()[static_cast<std::size_t>(f.arc)]);
    using CellFactor=Eigen::LLT<Eigen::Matrix<double,6,6>>;
    std::vector<CellFactor> cellFactors;cellFactors.reserve(d->elements.size());
    for(auto& e:d->elements){
        // One equivalent P2/P3 basis scale for every nonempty phase, from
        // its actual supplied positive volume and boundary quadrature nodes.
        double xmin=std::numeric_limits<double>::infinity(),xmax=-xmin,ymin=xmin,ymax=-xmin;
        const auto include=[&](Vector2D point){xmin=std::min(xmin,point.x);xmax=std::max(xmax,point.x);
            ymin=std::min(ymin,point.y);ymax=std::max(ymax,point.y);};
        for(const auto& node:e.geometry.volume)if(node.weight>0)include(node.point);
        for(const auto& face:e.geometry.faces)for(const auto& node:face.nodes)if(node.weight>0)include(node.point);
        const double extent=std::hypot(xmax-xmin,ymax-ymin);
        require(std::isfinite(extent)&&extent>0,"phase supplied-quadrature node extent is not positive finite");
        e.geometry.scale=extent;e.scale=extent;
        try{e.op=std::make_shared<E::LocalOperator>(std::move(e.geometry),1.0);}
        catch(const std::exception& failure){throw std::invalid_argument("[ElectricVofSharpField nativeCell="+std::to_string(e.cell)
            +" phase="+std::to_string(e.phase)+"] "+failure.what());}
        const auto& unit=e.op->stiffnessMatrix();
        // Sum in the common UNCONDENSED cell+face space; never add Schur complements.
        e.combined=p.epsilon[e.phase]*unit+d->dt*p.conductivity[e.phase]*unit;
        const Eigen::Index nf=e.combined.rows()-6;cellFactors.emplace_back(e.combined.topLeftCorner<6,6>());const auto& factor=cellFactors.back();
        require(factor.info()==Eigen::Success,"phase cell block lacks positive rank; unsupported small cut");
        e.cellFromFaces=factor.solve(-e.combined.topRightCorner(6,nf));
        e.condensed=e.combined.bottomRightCorner(nf,nf)+e.combined.bottomLeftCorner(nf,6)*e.cellFromFaces;
        require(e.cellFromFaces.allFinite()&&e.condensed.allFinite(),"nonfinite phase cell elimination");
        for(std::size_t i=0;i<e.faces.size();++i){const auto& fi=d->faces[static_cast<std::size_t>(e.faces[i])];if(fi.prescribed)continue;
            for(std::size_t j=0;j<e.faces.size();++j){const auto& fj=d->faces[static_cast<std::size_t>(e.faces[j])];
                const auto block=e.condensed.block<3,3>(3*static_cast<Eigen::Index>(i),3*static_cast<Eigen::Index>(j));
                if(fj.prescribed)rhs.segment<3>(fi.offset).noalias()-=block*fj.trace;
                else for(int k=0;k<3;++k)for(int l=0;l<3;++l)entries.emplace_back(fi.offset+k,fj.offset+l,block(k,l));}}}
    Eigen::SparseMatrix<double> matrix(free,free);matrix.setFromTriplets(entries.begin(),entries.end());matrix.makeCompressed();
    d->report.elements=d->elements.size();d->report.sharedFaces=d->faces.size();d->report.interfaceFaces=g.arcs().size();
    d->report.skeletonDofs=static_cast<std::size_t>(free);d->report.nonzeros=static_cast<std::size_t>(matrix.nonZeros());
    Eigen::SparseMatrix<double> transpose=matrix.transpose();
    const double matrixNorm=matrix.norm();require(free==0||(std::isfinite(matrixNorm)&&matrixNorm>0),"empty or nonfinite sharp skeleton operator");
    d->report.relativeSymmetryDefect=free==0?0:(matrix-transpose).norm()/matrixNorm;
    require(d->report.relativeSymmetryDefect<1e-12,"sharp skeleton symmetry defect exceeds original gate");
    d->report.assembleSeconds=elapsed(started);Eigen::VectorXd solution=Eigen::VectorXd::Zero(free);
#ifdef EIGEN_SUPERLU_SUPPORT
    Eigen::SuperLU<Eigen::SparseMatrix<double>> solver;
    if(free>0){auto timer=Timer::now();solver.compute(matrix);
        require(solver.info()==Eigen::Success,"sharp skeleton SuperLU factorization failed; no fallback");d->report.factorSeconds=elapsed(timer);
        timer=Timer::now();solution=solver.solve(rhs);d->report.solveSeconds=elapsed(timer);
        require(solver.info()==Eigen::Success&&solution.allFinite(),"sharp skeleton SuperLU solve failed; no fallback");d->report.usedSuperLU=true;}
#else
    require(free==0,"sharp skeleton requires compiled SuperLU support; no fallback");
#endif
    d->report.relativeLinearResidual=relative((matrix*solution-rhs).norm(),matrix.norm()*solution.norm()+rhs.norm());
    require(d->report.relativeLinearResidual<=controls.relativeLinearTolerance,"sharp skeleton linear residual failed");
    const auto recovery=Timer::now();double correctionSolveSeconds=0;
    for(auto& f:d->faces)if(!f.prescribed)f.trace=solution.segment<3>(f.offset);
    for(auto& e:d->elements){Eigen::VectorXd faces(3*static_cast<Eigen::Index>(e.faces.size()));
        for(std::size_t i=0;i<e.faces.size();++i)faces.segment<3>(3*static_cast<Eigen::Index>(i))=d->faces[static_cast<std::size_t>(e.faces[i])].trace;
        e.cellPolynomial=e.cellFromFaces*faces;}
    std::vector<Eigen::Vector3d> balance(d->faces.size(),Eigen::Vector3d::Zero()),scale(d->faces.size(),Eigen::Vector3d::Zero());
    std::vector<Triple> gauss(g.arcs().size()),conductive(g.arcs().size());
    std::vector<ExactTriple> exactGauss(g.arcs().size()),exactConductive(g.arcs().size());
    std::vector<double> gammaFluxScale(g.arcs().size(),0);
    int worstLocalCell=-1;unsigned worstLocalPhase=0;
    std::size_t worstNative=d->faces.size(),worstGamma=g.arcs().size(),worstUpdate=g.arcs().size();
    FV::Dyadic exactBoundaryDisplacement,exactBoundaryCurrent,exactBoundaryMagnitude,exactBoundaryCurrentMagnitude;
    FV::Dyadic globalGaussBudget,globalUpdateBudget;std::vector<Triple> updated;
    const auto audit=[&](){
    std::fill(balance.begin(),balance.end(),Eigen::Vector3d::Zero());std::fill(scale.begin(),scale.end(),Eigen::Vector3d::Zero());
    exactGauss.assign(g.arcs().size(),ExactTriple{});exactConductive.assign(g.arcs().size(),ExactTriple{});
    std::fill(gammaFluxScale.begin(),gammaFluxScale.end(),0);worstLocalCell=-1;worstLocalPhase=0;
    worstNative=d->faces.size();worstGamma=g.arcs().size();worstUpdate=g.arcs().size();
    exactBoundaryDisplacement={};exactBoundaryCurrent={};exactBoundaryMagnitude={};exactBoundaryCurrentMagnitude={};
    d->report.maximumRelativeLocalResidual=0;d->report.maximumRelativeNativeBalance=0;
    d->report.maximumRelativeGaussResidual=0;d->report.maximumRelativeChargeUpdateResidual=0;
    for(auto& e:d->elements){const Eigen::Index nf=3*static_cast<Eigen::Index>(e.faces.size());Eigen::VectorXd faces(nf);
        for(std::size_t i=0;i<e.faces.size();++i)faces.segment<3>(3*static_cast<Eigen::Index>(i))=d->faces[static_cast<std::size_t>(e.faces[i])].trace;
        e.potential=e.op->reconstruct(e.cellPolynomial,faces);
        e.localValues.resize(6+nf);e.localValues.head<6>()=e.cellPolynomial;e.localValues.tail(nf)=faces;
        const Eigen::VectorXd residual=e.combined*e.localValues;
        e.unitResidual=Eigen::VectorXd::Zero(e.localValues.size());
        const Eigen::VectorXd budget=e.combined.cwiseAbs()*e.localValues.cwiseAbs();
        const double localRelative=relative(residual.head<6>().norm(),budget.head<6>().norm());
        if(localRelative>d->report.maximumRelativeLocalResidual){worstLocalCell=e.cell;worstLocalPhase=e.phase;}
        d->report.maximumRelativeLocalResidual=std::max(d->report.maximumRelativeLocalResidual,localRelative);
        for(std::size_t i=0;i<e.faces.size();++i){const auto index=static_cast<std::size_t>(e.faces[i]);const auto& f=d->faces[index];
            const Eigen::Vector3d r=residual.segment<3>(6+3*static_cast<Eigen::Index>(i));
            balance[index]+=r;scale[index]+=budget.segment<3>(6+3*static_cast<Eigen::Index>(i));
            if(f.arc>=0||f.boundary>=0){Eigen::Vector3d physical=Eigen::Vector3d::Zero();
                for(int k=0;k<(f.arc>=0?3:1);++k){const auto row=6+3*static_cast<Eigen::Index>(i)+k;
                    const auto unit=exactRowDot(e.op->stiffnessMatrix(),row,e.localValues);
                    e.unitResidual[row]=exposed(unit);
                    const auto displacement=unit.multiplied(p.epsilon[e.phase]),current=unit.multiplied(p.conductivity[e.phase]);
                    if(f.arc>=0){const auto a=static_cast<std::size_t>(f.arc);physical[k]=exposed(displacement);
                        exactGauss[a][static_cast<std::size_t>(k)].add(displacement);exactConductive[a][static_cast<std::size_t>(k)].add(current);}
                    if(f.boundary>=0){exactBoundaryDisplacement.add(displacement,-1);exactBoundaryMagnitude.add(displacement.absolute());
                        exactBoundaryCurrent.add(current,-1);exactBoundaryCurrentMagnitude.add(current.absolute());}}
                if(f.arc>=0)gammaFluxScale[static_cast<std::size_t>(f.arc)]+=physical.norm();}}
    }
    updated=q.moments();
    for(std::size_t f=0;f<d->faces.size();++f)if(!d->faces[f].prescribed){const int a=d->faces[f].arc;
        const Eigen::Vector3d load=a<0?Eigen::Vector3d::Zero():vector(q.moments()[static_cast<std::size_t>(a)]);
        const double error=relative((balance[f]-load).norm(),scale[f].norm()+load.norm());
        if(error>d->report.maximumRelativeNativeBalance)worstNative=f;
        d->report.maximumRelativeNativeBalance=std::max(d->report.maximumRelativeNativeBalance,error);}
    for(std::size_t a=0;a<updated.size();++a){Eigen::Vector3d gaussDefect,updateDefect;
        for(unsigned k=0;k<3;++k){gauss[a][k]=exposed(exactGauss[a][k]);conductive[a][k]=exposed(exactConductive[a][k]);
            const auto transfer=exactConductive[a][k].multiplied(d->dt);
            if(transfer.compare(FV::Dyadic{})!=0){auto next=FV::Dyadic(q.moments()[a][k]);next.add(transfer,-1);updated[a][k]=exposed(next);}
            auto defect=FV::Dyadic(updated[a][k]);defect.add(exactGauss[a][k],-1);gaussDefect[k]=exposed(defect);
            auto update=FV::Dyadic(updated[a][k]);update.add(q.moments()[a][k],-1);update.add(transfer);updateDefect[k]=exposed(update);}
        const auto error=gaussDefect.norm();
        const auto budget=vector(updated[a]).norm()+gammaFluxScale[a];
        const double rel=relative(error,budget);
        if(rel>d->report.maximumRelativeGaussResidual)worstGamma=a;
        d->report.maximumRelativeGaussResidual=std::max(d->report.maximumRelativeGaussResidual,rel);
        const double updateError=updateDefect.norm();
        const double updateScale=vector(updated[a]).norm()+vector(q.moments()[a]).norm()+d->dt*vector(conductive[a]).norm();
        const double updateRelative=relative(updateError,updateScale);
        if(updateRelative>d->report.maximumRelativeChargeUpdateResidual)worstUpdate=a;
        d->report.maximumRelativeChargeUpdateResidual=std::max(d->report.maximumRelativeChargeUpdateResidual,updateRelative);}
    FV::Dyadic totalCharge,chargeMagnitude,inputCharge,inputMagnitude;
    for(const auto& m:updated){totalCharge.add(m[0]);chargeMagnitude.add(std::abs(m[0]));}
    for(const auto& m:q.moments()){inputCharge.add(m[0]);inputMagnitude.add(std::abs(m[0]));}
    d->report.boundaryDisplacement=exposed(exactBoundaryDisplacement);d->report.boundaryCurrent=exposed(exactBoundaryCurrent);
    auto globalGauss=exactBoundaryDisplacement;globalGauss.add(totalCharge,-1);
    globalGaussBudget=exactBoundaryMagnitude;globalGaussBudget.add(chargeMagnitude);
    auto globalUpdate=totalCharge;globalUpdate.add(inputCharge,-1);globalUpdate.add(exactBoundaryCurrent.multiplied(d->dt));
    globalUpdateBudget=chargeMagnitude;globalUpdateBudget.add(inputMagnitude);globalUpdateBudget.add(exactBoundaryCurrentMagnitude.multiplied(d->dt));
    d->report.globalGaussResidual=exposed(globalGauss);
    d->report.relativeGlobalGaussResidual=relative(std::abs(d->report.globalGaussResidual),
        exposed(globalGaussBudget));
    d->report.globalChargeUpdateResidual=exposed(globalUpdate);
    d->report.relativeGlobalChargeUpdateResidual=relative(std::abs(d->report.globalChargeUpdateResidual),
        exposed(globalUpdateBudget));
    return d->report.maximumRelativeLocalResidual<=controls.relativeConservationTolerance
        &&d->report.maximumRelativeNativeBalance<=controls.relativeConservationTolerance
        &&d->report.maximumRelativeGaussResidual<=controls.relativeConservationTolerance
        &&d->report.maximumRelativeChargeUpdateResidual<=controls.relativeConservationTolerance
        &&d->report.relativeGlobalGaussResidual<=controls.relativeConservationTolerance
        &&d->report.relativeGlobalChargeUpdateResidual<=controls.relativeConservationTolerance;
    };
    bool conserved=audit();d->report.initialMaximumRelativeGaussResidual=d->report.maximumRelativeGaussResidual;
    if(!conserved){
        const auto defectStarted=Timer::now();
        std::vector<E::CellVector> cellResidual(d->elements.size()),cellResponse(d->elements.size());
        std::vector<FV::Dyadic> freeResidual(static_cast<std::size_t>(free));
        const auto measureResidual=[&](){
            freeResidual.assign(static_cast<std::size_t>(free),FV::Dyadic{});double norm=0;
            for(const auto& face:d->faces)if(face.arc>=0)for(unsigned k=0;k<3;++k)
                freeResidual[static_cast<std::size_t>(face.offset)+k].add(q.moments()[static_cast<std::size_t>(face.arc)][k],-1);
            for(std::size_t ei=0;ei<d->elements.size();++ei){const auto& e=d->elements[ei];
                for(unsigned k=0;k<6;++k){cellResidual[ei][k]=exposed(exactRowDot(e.combined,k,e.localValues));norm=std::hypot(norm,cellResidual[ei][k]);}
                for(std::size_t i=0;i<e.faces.size();++i){const auto& face=d->faces[static_cast<std::size_t>(e.faces[i])];if(face.prescribed)continue;
                    for(unsigned k=0;k<3;++k)freeResidual[static_cast<std::size_t>(face.offset)+k].add(exactRowDot(e.combined,6+3*static_cast<Eigen::Index>(i)+k,e.localValues));}}
            for(const auto& value:freeResidual)norm=std::hypot(norm,exposed(value));return norm;
        };
        d->report.defectResidualMeasured=true;d->report.uncondensedResidualBefore=measureResidual();
        d->report.uncondensedResidualAfter=d->report.uncondensedResidualBefore;
        for(unsigned pass=0;pass<2&&!conserved;++pass){
            std::vector<FV::Dyadic> condensedResidual=freeResidual;
            for(std::size_t ei=0;ei<d->elements.size();++ei){const auto& e=d->elements[ei];
                cellResponse[ei]=cellFactors[ei].solve(cellResidual[ei]);require(cellResponse[ei].allFinite(),"nonfinite same-C cell defect response");
                for(std::size_t i=0;i<e.faces.size();++i){const auto& face=d->faces[static_cast<std::size_t>(e.faces[i])];if(face.prescribed)continue;
                    for(unsigned k=0;k<3;++k)for(unsigned j=0;j<6;++j)
                        condensedResidual[static_cast<std::size_t>(face.offset)+k].add(FV::Dyadic::product({e.combined(6+3*static_cast<Eigen::Index>(i)+k,j),cellResponse[ei][j]}),-1);}}
            Eigen::VectorXd correction=Eigen::VectorXd::Zero(free),correctionRhs(free);
            for(int k=0;k<free;++k)correctionRhs[k]=exposed(condensedResidual[static_cast<std::size_t>(k)].multiplied(-1));
            if(free>0){
#ifdef EIGEN_SUPERLU_SUPPORT
                const auto timer=Timer::now();correction=solver.solve(correctionRhs);const double seconds=elapsed(timer);
                d->report.solveSeconds+=seconds;correctionSolveSeconds+=seconds;
                require(solver.info()==Eigen::Success&&correction.allFinite(),"same-C Schur defect backsolve failed; no fallback");
#else
                require(false,"same-C Schur defect requires compiled SuperLU support");
#endif
            }
            for(std::size_t ei=0;ei<d->elements.size();++ei){auto& e=d->elements[ei];Eigen::VectorXd local=Eigen::VectorXd::Zero(3*static_cast<Eigen::Index>(e.faces.size()));
                for(std::size_t i=0;i<e.faces.size();++i){const auto& face=d->faces[static_cast<std::size_t>(e.faces[i])];
                    if(!face.prescribed)local.segment<3>(3*static_cast<Eigen::Index>(i))=correction.segment<3>(face.offset);}
                for(unsigned k=0;k<6;++k){FV::Dyadic value(e.cellPolynomial[k]);value.add(cellResponse[ei][k],-1);
                    value.add(exactRowDot(e.cellFromFaces,k,local));e.cellPolynomial[k]=exposed(value);}}
            for(int k=0;k<free;++k){FV::Dyadic value(solution[k]);value.add(correction[k]);solution[k]=exposed(value);}
            for(auto& face:d->faces)if(!face.prescribed)face.trace=solution.segment<3>(face.offset);
            ++d->report.defectCorrections;
            conserved=audit();d->report.uncondensedResidualAfter=measureResidual();
            d->report.relativeLinearResidual=relative((matrix*solution-rhs).norm(),matrix.norm()*solution.norm()+rhs.norm());
            require(d->report.relativeLinearResidual<=controls.relativeLinearTolerance,"same-C corrected skeleton linear residual failed");
        }
        d->report.defectCorrectionSeconds=elapsed(defectStarted);
    }
    if(!conserved){
        const auto sourceClock=d->source.timeLayer(),geometryClock=g.sourceTime();
        std::ostringstream diagnostic;diagnostic<<std::setprecision(std::numeric_limits<long double>::max_digits10)
            <<"sharp cell/native/Gamma conservation gate failed limit="<<controls.relativeConservationTolerance
            <<" local="<<d->report.maximumRelativeLocalResidual<<" native="<<d->report.maximumRelativeNativeBalance
            <<" Gamma="<<d->report.maximumRelativeGaussResidual<<" update="<<d->report.maximumRelativeChargeUpdateResidual
            <<" globalGauss="<<d->report.relativeGlobalGaussResidual<<" globalUpdate="<<d->report.relativeGlobalChargeUpdateResidual
            <<" linear="<<d->report.relativeLinearResidual<<" symmetry="<<d->report.relativeSymmetryDefect<<" snapshot="<<snapshot<<" dt="<<d->dt
            <<" defectCorrections="<<d->report.defectCorrections<<" uncondensedBefore="<<d->report.uncondensedResidualBefore
            <<" uncondensedAfter="<<d->report.uncondensedResidualAfter
            <<" geometryClock=("<<geometryClock.step<<','<<geometryClock.begin<<','<<geometryClock.end<<','<<static_cast<int>(geometryClock.position)<<')'
            <<" sourceClock=("<<sourceClock.step<<','<<sourceClock.begin<<','<<sourceClock.end<<','<<static_cast<int>(sourceClock.position)<<')'
            <<" targetClock=("<<clock.step<<','<<clock.begin<<','<<clock.end<<','<<static_cast<int>(clock.position)<<')'
            <<" worstLocal=("<<worstLocalCell<<','<<worstLocalPhase<<") worstNative="<<worstNative<<" worstUpdate="<<worstUpdate;
        if(worstNative<d->faces.size())diagnostic<<" nativeRow="<<d->faces[worstNative].nativeRow;
        if(worstGamma<g.arcs().size()){
            const auto arc=worstGamma;const auto& cell=g.cells()[static_cast<std::size_t>(g.arcs()[arc].nativeCell)];
            Triple physicalDefect{};ExactTriple combinedRows,legacyMaterial;
            for(unsigned k=0;k<3;++k){auto defect=FV::Dyadic(updated[arc][k]);defect.add(exactGauss[arc][k],-1);physicalDefect[k]=exposed(defect);}
            diagnostic<<" worstGamma=(arc="<<arc<<",cell="<<cell.nativeCell<<",error="<<vector(physicalDefect).norm()
                <<",physicalBudget="<<vector(updated[arc]).norm()+gammaFluxScale[arc]<<",fluxBudget="<<gammaFluxScale[arc]<<')'
                <<" Qadv=("<<vector(q.moments()[arc]).transpose()<<") Qnew=("<<vector(updated[arc]).transpose()
                <<") Gauss=("<<vector(gauss[arc]).transpose()<<") J=("<<vector(conductive[arc]).transpose()<<')';
            const auto& face=d->faces[static_cast<std::size_t>(d->interfaceFace[arc])];
            for(const int index:face.elements){const auto& e=d->elements[static_cast<std::size_t>(index)];const auto& actual=e.op->geometry();
                const auto local=std::distance(e.faces.begin(),std::find(e.faces.begin(),e.faces.end(),d->interfaceFace[arc]));
                const Eigen::VectorXd legacy=e.op->stiffnessMatrix()*e.localValues;
                for(unsigned k=0;k<3;++k){const auto row=6+3*local+k;combinedRows[k].add(exactRowDot(e.combined,row,e.localValues));
                    legacyMaterial[k].add(FV::Dyadic::product({p.epsilon[e.phase],legacy[row]}));
                    legacyMaterial[k].add(FV::Dyadic::product({d->dt,p.conductivity[e.phase],legacy[row]}));}
                double xmin=std::numeric_limits<double>::infinity(),xmax=-xmin,ymin=xmin,ymax=-xmin;long double area=0;
                const auto include=[&](Vector2D point){xmin=std::min(xmin,point.x);xmax=std::max(xmax,point.x);
                    ymin=std::min(ymin,point.y);ymax=std::max(ymax,point.y);};
                for(const auto& node:actual.volume)if(node.weight>0){area+=node.weight;include(node.point);}
                for(const auto& boundary:actual.faces)for(const auto& node:boundary.nodes)if(node.weight>0)include(node.point);
                diagnostic<<" phase=("<<e.phase<<",inventory="<<cell.inventory[e.phase]<<",ordinaryArea="<<area<<",basisScale="<<e.scale
                    <<",nodeExtent="<<std::hypot(xmax-xmin,ymax-ymin)<<",bbox="<<xmin<<','<<xmax<<','<<ymin<<','<<ymax<<')';}
            Triple solveResidual{},materialDifference{},dotDifference{},chargeRounding{},splitClosure{};
            for(unsigned k=0;k<3;++k){auto implicit=exactGauss[arc][k];implicit.add(exactConductive[arc][k].multiplied(d->dt));
                auto r=combinedRows[k];r.add(q.moments()[arc][k],-1);solveResidual[k]=exposed(r);
                auto m=combinedRows[k];m.add(implicit,-1);materialDifference[k]=exposed(m);
                auto old=legacyMaterial[k];old.add(implicit,-1);dotDifference[k]=exposed(old);
                auto qround=FV::Dyadic(updated[arc][k]);qround.add(q.moments()[arc][k],-1);
                qround.add(exactConductive[arc][k].multiplied(d->dt));chargeRounding[k]=exposed(qround);
                auto closure=m;closure.add(r,-1);closure.add(qround);closure.add(updated[arc][k],-1);closure.add(exactGauss[arc][k]);splitClosure[k]=exposed(closure);}
            diagnostic<<" R_combined=("<<vector(solveResidual).transpose()<<") M_material=("<<vector(materialDifference).transpose()
                <<") oldDotCombinedMinusExact=("<<vector(dotDifference).transpose()<<") QnewRounding=("<<vector(chargeRounding).transpose()
                <<") MR_closure=("<<vector(splitClosure).transpose()<<')';
        }
        diagnostic<<" globalGaussAbs="<<d->report.globalGaussResidual<<" globalGaussBudget="<<exposed(globalGaussBudget)
            <<" globalUpdateAbs="<<d->report.globalChargeUpdateResidual<<" globalUpdateBudget="<<exposed(globalUpdateBudget);
        const auto reason=diagnostic.str();require(false,reason.c_str());
    }
    if(!snapshot)d->candidate=ElectricVofSurfaceCharge::fromPrepared(mesh,g,clock,std::move(updated));
    for(std::size_t a=0;a<g.arcs().size();++a){auto& f=d->faces[static_cast<std::size_t>(d->interfaceFace[a])];
        const Eigen::LLT<Eigen::Matrix3d> factor(f.mass);
        for(unsigned phase=0;phase<2;++phase){const auto& element=d->elements[static_cast<std::size_t>(f.elements[phase])];
            const auto found=std::find(element.faces.begin(),element.faces.end(),d->interfaceFace[a]);
            require(found!=element.faces.end(),"Gamma absent from phase element");const auto local=std::distance(element.faces.begin(),found);
            const Eigen::Vector3d residual=element.unitResidual.segment<3>(6+3*local);
            f.normalElectricCoefficient[phase]=(phase==0?-1:1)*factor.solve(residual);}
        f.chargeDensityCoefficient=factor.solve(vector(d->candidate.moments()[a]));}
    // Result queries require recovered polynomials and shared Gamma traces,
    // never the local factors, matrix blocks or quadrature work arrays.
    for(auto& e:d->elements){e.op.reset();e.combined.resize(0,0);e.cellFromFaces.resize(0,0);e.condensed.resize(0,0);
        e.localValues.resize(0);e.unitResidual.resize(0);std::vector<int>().swap(e.faces);}
    std::vector<Face> gamma;gamma.reserve(g.arcs().size());
    for(std::size_t a=0;a<g.arcs().size();++a){auto face=std::move(d->faces[static_cast<std::size_t>(d->interfaceFace[a])]);
        std::vector<E::FaceQuadraturePoint>().swap(face.quadrature);gamma.push_back(std::move(face));d->interfaceFace[a]=static_cast<int>(a);}
    d->faces=std::move(gamma);
    d->report.recoverySeconds=elapsed(recovery)-correctionSolveSeconds;return ElectricVofSharpField(std::move(d));
}
const G& ElectricVofSharpField::geometry() const{require(bool(data_),"consumed sharp field");return data_->geometry;}
const ElectricVofSurfaceCharge& ElectricVofSharpField::sourceCharge() const{require(bool(data_),"consumed sharp field");return data_->source;}
const ElectricVofSurfaceCharge& ElectricVofSharpField::chargeCandidate() const{require(bool(data_),"consumed sharp field");return data_->candidate;}
const ElectricVofSurfaceCharge& ElectricVofSharpField::implicitInputCharge() const{require(bool(data_),"consumed sharp field");return data_->input;}
const ElectricSurfaceChargeTransport* ElectricVofSharpField::advection() const{require(bool(data_),"consumed sharp field");return data_->advection.get();}
Clock ElectricVofSharpField::candidateTime() const{require(bool(data_),"consumed sharp field");return data_->clock;}
double ElectricVofSharpField::effectiveDt() const{require(bool(data_),"consumed sharp field");return data_->dt;}
const ElectricVofSharpProblem& ElectricVofSharpField::materials() const{require(bool(data_),"consumed sharp field");return data_->problem;}
const ElectricVofSharpReport& ElectricVofSharpField::report() const{require(bool(data_),"consumed sharp field");return data_->report;}
bool ElectricVofSharpField::sameInstance(const ElectricVofSharpField& f) const noexcept{return data_&&data_==f.data_;}
Vector2D ElectricVofSharpField::sampleElectricField(int cell,unsigned phase,Vector2D point) const
{return sampleElectricFieldAt(cell,phase,{{point.x,point.y},{0,0}});}
Vector2D ElectricVofSharpField::sampleElectricFieldAt(int cell,unsigned phase,
    Geometry::PlanarCellCut::AnchoredPoint point) const
{
    require(bool(data_)&&cell>=0&&static_cast<std::size_t>(cell)<data_->elementByCell.size()&&phase<2,"invalid sharp element query");
    require(std::isfinite(point.anchor.x)&&std::isfinite(point.anchor.y)
        &&std::isfinite(point.offset.x)&&std::isfinite(point.offset.y),"nonfinite sharp query site");
    const auto& rectangle=data_->geometry.cells()[static_cast<std::size_t>(cell)].rectangle;
    // Four independent signed distances. In particular, a tiny upper-face
    // excursion must not disappear in a rounded lower-distance/width test.
    const auto fromLower=Geometry::PlanarCellCut::anchoredDifference(point,
        {{rectangle.lower.x,rectangle.lower.y},{0,0}});
    const auto toUpper=Geometry::PlanarCellCut::anchoredDifference(
        {{rectangle.upper.x,rectangle.upper.y},{0,0}},point);
    require(fromLower.x>=0&&fromLower.y>=0&&toUpper.x>=0&&toUpper.y>=0,
        "sharp field point must remain in selected native cell");
    const int index=data_->elementByCell[static_cast<std::size_t>(cell)][phase];require(index>=0,"requested sharp phase is empty");
    const auto& e=data_->elements[static_cast<std::size_t>(index)];
    // Contract the complete actual site minus origin and local center before
    // its one binary64 rounding; never round a partial long-double frame.
    const auto centered=Geometry::PlanarCellCut::anchoredDifferenceDouble(point,
        {{e.coordinateOrigin.x,e.coordinateOrigin.y},{0,0}},{e.center.x,e.center.y});
    const Eigen::Vector2d gradient=E::p3Gradient(centered,{0,0},e.scale)*e.potential;
    return {-gradient[0],-gradient[1]};
}
ElectricVofSharpInterfaceField ElectricVofSharpField::interfaceField(std::size_t arc,double t) const
{
    require(bool(data_)&&arc<data_->interfaceFace.size()&&std::isfinite(t)&&t>=0&&t<=1,"invalid sharp interface query");
    const auto& f=data_->faces[static_cast<std::size_t>(data_->interfaceFace[arc])];
    const auto& curve=data_->geometry.arcs()[arc];const auto derivative=curve.derivative(t);
    const double speed=static_cast<double>(std::hypot(derivative.x,derivative.y));require(speed>0,"degenerate sharp tangent");
    ElectricVofSharpInterfaceField result;result.normal=asVector(curve.normal(t));
    result.tangent={static_cast<double>(derivative.x)/speed,static_cast<double>(derivative.y)/speed};
    const double coordinate=2*t-1;const auto basis=E::faceBasis(coordinate);
    result.tangentialElectricField=-(2*f.trace[1]+6*coordinate*f.trace[2])/speed;
    for(unsigned phase=0;phase<2;++phase){
        result.normalElectricField[phase]=basis.dot(f.normalElectricCoefficient[phase]);
        result.electricField[phase]={result.normalElectricField[phase]*result.normal.x+result.tangentialElectricField*result.tangent.x,
            result.normalElectricField[phase]*result.normal.y+result.tangentialElectricField*result.tangent.y};}
    result.surfaceChargeDensity=basis.dot(f.chargeDensityCoefficient);return result;
}
}

