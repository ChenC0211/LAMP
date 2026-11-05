#include "Modules/Hydrodynamics/HydroMaxwellLoad.h"
#include "Modules/Electric/ElectricSurfaceChargeTransport.h"
#include "Modules/Interface/GeometricVofInitialization.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/Detail/FiniteVolumeClosure.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
using namespace LAMP::Hydrodynamics;
using namespace LAMP::Interface;
using namespace LAMP::Electric;
using Clock=LAMP::TimeLayer;
unsigned assertions=0;
void check(bool yes,const char* message){++assertions;if(!yes)throw std::runtime_error(message);}
template<class F>void rejects(F&& f,const char* message)
{bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}check(rejected,message);}
void close(long double a,long double b,long double scale,const char* message,long double tolerance=1e-12L)
{check(std::isfinite(a)&&std::isfinite(b)&&std::isfinite(scale)&&scale>=0,message);check(scale==0?a==b:std::abs(a-b)<=tolerance*scale,message);}
struct ReferenceMaxwellWork {long double momentumWork=0,volumeWork=0,boundaryWork=0,residual=0,relativeResidual=0,absoluteScale=0;};
long double nativePower(const HydroMaxwellLoad& load,const HydroNativeValues& values)
{
    LAMP::Geometry::Detail::FiniteVolume::Dyadic sum;const auto& u=values.normalVelocity();const auto& rhs=load.momentumRhsLoad();
    check(u.size()==rhs.size(),"sharp actual native power size");
    for(std::size_t i=0;i<u.size();++i)sum.add(LAMP::Geometry::Detail::FiniteVolume::Dyadic::product({u[i],rhs[i]}));
    return sum.rounded();
}
struct Tensor {long double xx,xy,yy;};
Tensor exactTensor(unsigned phase)
{
    const long double x=-.37L,y=phase==0?-.8L:-1.4L,eps=phase==0?2.3L:.9L;
    return {eps*(x*x-y*y)/2,eps*x*y,eps*(y*y-x*x)/2};
}
// A has a named variable-RT0 Green pairing. Keep this public literal
// reference separate so the original B1 strain reference below is retained.
ReferenceMaxwellWork greenWork(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,const HydroNativeValues& values)
{
    namespace FV=LAMP::Geometry::Detail::FiniteVolume;
    check(!load.isSharp()&&load.report().pairing==HydroMaxwellPairing::Rt0Green,"A reference requires the declared Green pairing");
    load.validateFor(mesh,load.pressureLoad(),load.phaseGeometry(),load.stress(),load.phaseGeometry().timeLayer());
    const auto& op=load.geometry();const auto& layout=op.layout();values.validateFor(layout,mesh);
    const auto& u=values.normalVelocity();const auto& v=op.stressVolumeVelocityMatrix();const auto& b=op.stressFaceVelocityMatrix();
    FV::Dyadic momentum,volume,boundary,scale,budget;
    const auto add=[&](FV::Dyadic& sum,std::initializer_list<long double> factors){const auto term=FV::Dyadic::product(factors);
        sum.add(term);scale.add(term.absolute());budget.add(term.absolute().multiplied(1e-12L));};
    for(std::size_t i=0;i<u.size();++i)add(momentum,{u[i],load.momentumRhsLoad()[i]});
    for(std::size_t q=0;q<op.points().size();++q){const auto& point=op.points()[q];
        const auto sample=load.stress().sample(layout.nativeCell(point.cellDof),point.point);
        const std::array<double,2> div{{sample.tensorDivergence.x,sample.tensorDivergence.y}};
        for(int axis=0;axis<2;++axis)for(HydroStrainOperator::Matrix::InnerIterator it(v,static_cast<Eigen::Index>(q)*2+axis);it;++it)
            add(volume,{point.physicalWeight,div[axis],it.value(),u[it.col()]});}
    std::size_t cellSidePoints=0;
    for(std::size_t q=0;q<op.stressFacePoints().size();++q){const auto& face=op.stressFacePoints()[q];if(face.physicalWeight==0)continue;
        const bool vertical=layout.normalAxis(face.velocityDof)==HydroDofLayout::NormalAxis::First;
        for(int side=0;side<2;++side){const int owner=side==0?face.negativeCellDof:face.positiveCellDof;if(owner<0)continue;++cellSidePoints;
            const auto t=load.stress().sample(layout.nativeCell(owner),face.point).tensor;
            const std::array<double,2> canonical{{vertical?t.xx:t.xy,vertical?t.xy:t.yy}};
            for(int axis=0;axis<2;++axis)for(HydroStrainOperator::Matrix::InnerIterator it(b,static_cast<Eigen::Index>(q)*2+axis);it;++it){
                add(volume,{(side==0?-1.L:1.L)*face.physicalWeight,canonical[axis],it.value(),u[it.col()]});
                if(face.negativeCellDof<0||face.positiveCellDof<0)
                    add(boundary,{face.physicalWeight,(side==0?1.L:-1.L)*canonical[axis],it.value(),u[it.col()]});}}}
    check(cellSidePoints==load.report().cellSidePoints,"A Green reference missed an actual child-face side");
    auto residual=momentum;residual.add(volume,-1);residual.add(boundary,-1);
    check(residual.absolute().compare(budget)<=0,"original per-quadrature Maxwell work gate failed");
    ReferenceMaxwellWork result;result.momentumWork=momentum.rounded();result.volumeWork=volume.rounded();result.boundaryWork=boundary.rounded();
    result.residual=residual.rounded();result.absoluteScale=scale.rounded();
    result.relativeResidual=result.absoluteScale==0?0:residual.absolute().rounded()/result.absoluteScale;return result;
}
ReferenceMaxwellWork originalWork(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,const HydroNativeValues& values)
{
    if(!load.isSharp())return greenWork(mesh,load,values);
    check(load.report().pairing==HydroMaxwellPairing::CellMeanStrain,"B1 reference requires its original strain pairing");
    namespace FV=LAMP::Geometry::Detail::FiniteVolume;
    const auto& op=load.geometry();const auto& layout=op.layout();const auto& u=values.normalVelocity();
    const auto& points=op.points();const auto& faces=op.stressFacePoints();const auto& e=op.matrix();
    const auto& b=op.stressFaceVelocityMatrix();const bool rz=layout.domain().axisymmetric();
    const auto tensorAt=[&](int cell,Vector2D point){
        if(load.isSharp()){const auto& geometry=load.sharpStress().field().geometry();
            return load.sharpStress().sample(cell,geometry.phaseAt(cell,{point.x,point.y}),point).tensor;}
        const auto t=load.stress().sample(cell,point).tensor;
        return ElectricVofSharpMaxwellTensor{t.xx,t.xy,t.yy,t.hoop};};
    FV::Dyadic momentum,volume,boundary,scale,budget;
    const auto add=[&](FV::Dyadic& sum,std::initializer_list<long double> factors){const auto term=FV::Dyadic::product(factors);
        sum.add(term);scale.add(term.absolute());budget.add(term.absolute().multiplied(1e-12L));};
    for(std::size_t i=0;i<u.size();++i)add(momentum,{u[i],load.momentumRhsLoad()[i]});
    for(std::size_t q=0;q<points.size();++q){const auto& point=points[q];const auto t=tensorAt(layout.nativeCell(point.cellDof),point.point);
        const double shear=FV::roundDouble(FV::Dyadic::product({std::sqrt(2.),t.xy})).value;
        const std::array<double,4> packed{{t.xx,t.yy,rz?t.hoop:shear,shear}};
        for(int k=0;k<op.components();++k)for(HydroStrainOperator::Matrix::InnerIterator it(e,static_cast<Eigen::Index>(q)*op.components()+k);it;++it)
            add(volume,{-point.physicalWeight,packed[k],it.value(),u[it.col()]});}
    for(std::size_t q=0;q<faces.size();++q){const auto& face=faces[q];
        if((face.negativeCellDof>=0&&face.positiveCellDof>=0)||face.physicalWeight==0)continue;
        const int owner=face.negativeCellDof>=0?face.negativeCellDof:face.positiveCellDof;
        const auto t=tensorAt(layout.nativeCell(owner),face.point);const double sign=face.negativeCellDof>=0?1.:-1.;
        const bool vertical=layout.normalAxis(face.velocityDof)==HydroDofLayout::NormalAxis::First;
        const std::array<double,2> traction{{sign*(vertical?t.xx:t.xy),sign*(vertical?t.xy:t.yy)}};
        for(unsigned axis=0;axis<2;++axis)for(HydroStrainOperator::Matrix::InnerIterator it(b,static_cast<Eigen::Index>(q)*2+axis);it;++it)
            add(boundary,{face.physicalWeight,traction[axis],it.value(),u[it.col()]});}
    auto residual=momentum;residual.add(volume,-1);residual.add(boundary,-1);
    check(residual.absolute().compare(budget)<=0,"original per-quadrature Maxwell work gate failed");
    ReferenceMaxwellWork result;result.momentumWork=momentum.rounded();result.volumeWork=volume.rounded();result.boundaryWork=boundary.rounded();
    result.residual=residual.rounded();result.absoluteScale=scale.rounded();
    result.relativeResidual=result.absoluteScale==0?0:residual.absolute().rounded()/result.absoluteScale;return result;
}
void exactWorkEquivalence(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,const HydroNativeValues& values)
{
    const auto a=nativePower(load,values);const auto b=originalWork(mesh,load,values);
    check(a==b.momentumWork&&std::signbit(a)==std::signbit(b.momentumWork)&&b.relativeResidual<=1e-12L,
        "actual RHS exact-dot power or independent Green/B1 virtual-work proof changed");
}
void verifyWork(const QuadtreeMesh& mesh,const HydroMaxwellLoad& load,bool flat)
{
    const auto& op=load.geometry();const auto& layout=op.layout();
    std::vector<double> u(static_cast<std::size_t>(layout.velocityCount())),p(static_cast<std::size_t>(layout.pressureCount()));
    for(std::size_t i=0;i<u.size();++i)u[i]=std::sin(.37*static_cast<double>(i+1))+.03*std::cos(.11*static_cast<double>(i+1));
    for(std::size_t i=0;i<p.size();++i)p[i]=.2+.01*static_cast<double>(i);
    const HydroNativeValues values(layout,u,p);const auto work=originalWork(mesh,load,values);
    exactWorkEquivalence(mesh,load,values);
    check(work.relativeResidual<=1e-12L,"independent native virtual-work identity failed");
    const auto e=op.apply(mesh,values);const auto& v=op.stressFaceVelocityMatrix();
    const Eigen::Map<const Eigen::VectorXd> native(u.data(),static_cast<Eigen::Index>(u.size()));const Eigen::VectorXd faceVelocity=v*native;
    long double volume=0,boundary=0,momentum=0,scale=0;
    std::vector<std::array<long double,2>> phaseWeight(static_cast<std::size_t>(layout.pressureCount()));
    std::vector<long double> integratedDivergence(phaseWeight.size()),divergenceScale(phaseWeight.size());
    const auto& smooth=load.sharpStress().field().geometry();
    for(std::size_t q=0;q<op.points().size();++q){const auto& point=op.points()[q];const int cell=layout.nativeCell(point.cellDof);
        const auto s=flat?(point.point.y<.423?0u:1u):smooth.phaseAt(cell,{point.point.x,point.point.y});
        phaseWeight[static_cast<std::size_t>(point.cellDof)][s]+=point.physicalWeight;
        const auto actual=load.sharpStress().sample(cell,s,point.point).tensor;
        const Tensor t=flat?exactTensor(s):Tensor{actual.xx,actual.xy,actual.yy};
        const auto base=q*3;const auto term=-point.physicalWeight*(t.xx*e[base]+t.yy*e[base+1]+std::sqrt(2.L)*t.xy*e[base+2]);
        volume+=term;scale+=std::abs(term);
        integratedDivergence[static_cast<std::size_t>(point.cellDof)]+=point.physicalWeight*(static_cast<long double>(e[base])+e[base+1]);
        divergenceScale[static_cast<std::size_t>(point.cellDof)]+=point.physicalWeight*(std::abs(static_cast<long double>(e[base]))+std::abs(static_cast<long double>(e[base+1])));}
    const auto actualDivergence=op.metrics().divergence().apply(mesh,values);
    for(int row=0;row<layout.pressureCount();++row){const auto& inventory=smooth.cells().at(static_cast<std::size_t>(layout.nativeCell(row))).inventory;
        for(unsigned s=0;s<2;++s)close(phaseWeight[static_cast<std::size_t>(row)][s],inventory[s],inventory[s],"partitioned strain lost actual own-phase volume");
        close(integratedDivergence[static_cast<std::size_t>(row)],actualDivergence[static_cast<std::size_t>(row)],divergenceScale[static_cast<std::size_t>(row)]+std::abs(actualDivergence[static_cast<std::size_t>(row)]),"Ehat trace is not the original integral divergence");}
    for(std::size_t q=0;q<op.stressFacePoints().size();++q){const auto& face=op.stressFacePoints()[q];
        if(face.negativeCellDof>=0&&face.positiveCellDof>=0)continue;
        const int owner=face.negativeCellDof>=0?face.negativeCellDof:face.positiveCellDof;
        const int cell=layout.nativeCell(owner);const auto s=flat?(face.point.y<.423?0u:1u):smooth.phaseAt(cell,{face.point.x,face.point.y});
        const auto actual=load.sharpStress().sample(cell,s,face.point).tensor;
        const Tensor t=flat?exactTensor(s):Tensor{actual.xx,actual.xy,actual.yy};
        const long double sign=face.negativeCellDof>=0?1:-1;const bool vertical=layout.normalAxis(face.velocityDof)==HydroDofLayout::NormalAxis::First;
        const auto tx=sign*(vertical?t.xx:t.xy),ty=sign*(vertical?t.xy:t.yy);
        const auto term=face.physicalWeight*(tx*faceVelocity[static_cast<Eigen::Index>(2*q)]+ty*faceVelocity[static_cast<Eigen::Index>(2*q+1)]);
        boundary+=term;scale+=std::abs(term);}
    for(std::size_t i=0;i<u.size();++i){const auto term=static_cast<long double>(u[i])*load.momentumRhsLoad()[i];momentum+=term;scale+=std::abs(term);}
    close(momentum,work.momentumWork,scale,"independent momentum work differs from actual native RHS");
    close(work.volumeWork,volume,scale,"phase-resolved volume work differs from independent tensor contraction",flat?3e-8L:1e-12L);
    close(work.boundaryWork,boundary,scale,"physical boundary work differs from independent traction contraction",flat?3e-8L:1e-12L);
    close(momentum,volume+boundary,scale,"sharp load failed independent native virtual work",flat?3e-8L:1e-12L);
    const auto& d=load.pressureLoad().divergence().matrix();const auto& c=load.pressureLoad().matrix();
    for(Eigen::Index row=0;row<d.outerSize();++row)for(HydroIntegralDivergence::Matrix::InnerIterator it(d,row);it;++it)
        check(c.coeff(it.col(),it.row())==-it.value(),"original C=-D_I^T pressure pair changed");
    const auto zero=HydroEnergyViscosity::zero(mesh,op.metrics());
    check(zero.matrix().nonZeros()==0&&zero.strain().points().empty()&&zero.dissipation(mesh,values)==0,"mu=0 is not the exact original zero path");
    std::cout<<"sharp_work flat="<<flat<<" volume_points="<<op.points().size()<<" face_points="<<op.stressFacePoints().size()
        <<" virtual_work_residual="<<static_cast<double>(work.relativeResidual)<<'\n';
}
void run(QuadtreeMesh& mesh,bool curved)
{
    const auto source=GeometricVofInitialization::prepare(mesh,false,[=](Vector2D p){return p.y-(curved?.421+.031*std::sin(2*std::acos(-1.)*p.x):.423);});
    const auto stocks=source.plicCells();VofSmoothGeometryConfig gc;gc.allowBoundaryEndpoints=true;
    const auto oldSmooth=VofSmoothGeometry::prepare(mesh,source,gc);
    const double density=curved?0:2.3*.8-.9*1.4;
    const auto accepted=ElectricVofSurfaceCharge::fromDensity(mesh,oldSmooth,[=](std::size_t,double){return density;});const auto moments=accepted.moments();
    Clock end{1,0,.02,Clock::Position::End};std::vector<double> zeros(static_cast<std::size_t>(source.domain().nativeFaceCount()));
    const auto slab=GeometricVofAdvectorSlab::prepare(mesh,source,end,zeros,zeros);
    const auto transport=GeometricVofTransport::prepareLocalPlicSweep(mesh,source,slab,1,1);
    const auto layout=HydroDofLayout::build(mesh,source.domain());const auto divergence=HydroIntegralDivergence::build(mesh,layout);
    const auto metrics=HydroPhysicalDivergence::build(mesh,divergence);const auto pressure=HydroPressureLoad::build(divergence);
    const auto phase=HydroPhaseGeometry::build(mesh,metrics,transport.midpoint());
    const auto midpoint=VofSmoothGeometry::prepare(mesh,transport.midpoint(),gc);
    const auto moved=ElectricSurfaceChargeTransport::prepareStage(mesh,accepted,transport,midpoint);
    ElectricVofSharpProblem problem;problem.epsilon=curved?std::array<double,2>{{1,1}}:std::array<double,2>{{2.3,.9}};
    for(auto& b:problem.boundaries)b.value=[=](Vector2D p){return .37*p.x+(curved?.8:(p.y<=.423?.8:1.4))*(p.y-.423);};
    const auto field=ElectricVofSharpField::prepareOhmic(mesh,moved,problem);const auto stress=ElectricVofSharpMaxwellStress::prepare(mesh,field);
    const auto load=HydroMaxwellLoad::prepareSharp(mesh,pressure,phase,stress);
    check(load.isSharp()&&load.sharpStress().sameInstance(stress),"sharp load lost actual stress owner");
    check(field.advection()!=nullptr&&field.sourceCharge().sameInstance(accepted)&&field.effectiveDt()==.01,"sharp Midpoint used a forged charge clock");
    load.validateSharpFor(mesh,pressure,phase,stress,transport.midpoint().timeLayer());
    rejects([&]{load.stress();},"sharp load masqueraded as A stress");
    rejects([&]{load.validateSharpFor(mesh,pressure,phase,stress,end);},"End field was accepted as applied Midpoint stress");
    const auto foreignPressure=HydroPressureLoad::build(divergence);
    rejects([&]{load.validateSharpFor(mesh,foreignPressure,phase,stress,phase.timeLayer());},"equivalent pressure owner replaced actual C");
    const auto foreignStress=ElectricVofSharpMaxwellStress::prepare(mesh,field);
    rejects([&]{load.validateSharpFor(mesh,pressure,phase,foreignStress,phase.timeLayer());},"equivalent sharp stress replaced actual owner");
    const auto initialField=ElectricVofSharpField::prepareElectrostatic(mesh,oldSmooth,accepted,problem);
    rejects([&]{HydroMaxwellLoad::prepareSharp(mesh,pressure,phase,ElectricVofSharpMaxwellStress::prepare(mesh,initialField));},"Initial stress was relabeled as Midpoint");
    verifyWork(mesh,load,!curved);
    if(!curved){
        std::vector<double> emptyCharge(source.alpha().size());const auto oldCharge=ElectricVofCharge::fromIntegrated(mesh,source,emptyCharge);
        ElectricVofMaterials material;std::vector<ElectricVofBoundary> boundaries(zeros.size());
        for(int row=0;row<source.domain().nativeFaceCount();++row){const auto f=source.domain().face(mesh,row);
            if(f.negativeCell<0||f.positiveCell<0){boundaries[row].type=ElectricVofBoundary::Type::Dirichlet;boundaries[row].value=f.center.x+2*f.center.y;}}
        const auto aField=ElectricVofField::prepareOhmic(mesh,transport.midpoint(),oldCharge,phase.timeLayer(),material,boundaries);
        const auto aStress=ElectricVofMaxwellStress::prepare(mesh,aField);const auto aLoad=HydroMaxwellLoad::prepare(mesh,pressure,phase,aStress);
        check(!aLoad.isSharp() && &aLoad.stress().cells()==&aStress.cells(),"A load stress semantics changed");
        rejects([&]{aLoad.sharpStress();},"A load masqueraded as a sharp stress");
        const HydroNativeValues values(layout,std::vector<double>(static_cast<std::size_t>(layout.velocityCount()),.1),std::vector<double>(static_cast<std::size_t>(layout.pressureCount())));
        check(greenWork(mesh,aLoad,values).relativeResidual<=1e-12L,"independent A Green virtual-work proof failed");
        exactWorkEquivalence(mesh,aLoad,values);
    }
    check(accepted.moments()==moments,"load/rejection changed accepted surface charge");
    for(std::size_t i=0;i<stocks.size();++i)check(source.plicCells()[i].inventory==stocks[i].inventory,"load changed accepted VOF inventory");
}
}
int main()
{
    try{DictionaryMesh d;d.Nx=d.Ny=8;d.maxCellNum=512;d.maxInitialIterations=0;QuadtreeMesh mesh(d);
        run(mesh,false);run(mesh,true);std::cout<<"hydro_sharp_maxwell_load_contract PASS assertions="<<assertions<<'\n';return 0;}
    catch(const std::exception& error){std::cerr<<"hydro_sharp_maxwell_load_contract FAIL "<<error.what()<<'\n';return 1;}
}

