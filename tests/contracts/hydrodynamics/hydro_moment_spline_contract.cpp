#include "Modules/Hydrodynamics/HydroCompatibleMomentSpline.h"
#include "Modules/Hydrodynamics/Detail/HydroCompatibleMomentBlockSolve.h"
#include "Modules/Hydrodynamics/Detail/HydroCompatibleMomentSplineBasis.h"
#include "Modules/Hydrodynamics/Detail/HydroCompatibleMomentSplineAssembly.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentFlow.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentBoundary.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTraction.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMassCurrent.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include "hydro_complete_field_observation.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/Geometry/ParametricHermiteCurve.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Core/Field/HEdgeField.h"
#include "Core/Field/VEdgeField.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <atomic>
#include <chrono>
#include <numeric>
#include <Eigen/QR>
#include <Eigen/SuperLUSupport>

namespace {
using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
using Spline=Quantum::Hydrodynamics::HydroCompatibleMomentSpline;
using Flow=Quantum::Hydrodynamics::HydroCompatibleMomentFlow;
void require(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
DictionaryMesh input(int n){DictionaryMesh p;p.Nx=p.Ny=n;p.minCoord={0,0};p.maxCoord={1,1};p.maxCellNum=10000;p.maxInitialIterations=0;return p;}
void kineticQuadrature(double* rx,double* rw){
    const double pi=std::acos(-1.);
    for(int k=0;k<12;k++){double x=std::cos(pi*(k+.75)/12.5),derivative=0;
        for(int iteration=0;iteration<30;iteration++){double a=1,b=x;for(int degree=2;degree<=12;degree++){const double c=((2*degree-1)*x*b-(degree-1)*a)/degree;a=b;b=c;}
            derivative=12*(x*b-a)/(x*x-1);const double next=x-b/derivative;if(next==x)break;x=next;}rx[k]=x;rw[k]=2/((1-x*x)*derivative*derivative);}
}
Spline::FirstJet exact(Vector2D p){
    const double x=p.x,y=p.y,f=x*x*(1-x)*(1-x),fx=2*x-6*x*x+4*x*x*x,fxx=2-12*x+12*x*x;
    const double g=y*y*(1-y)*(1-y),gy=2*y-6*y*y+4*y*y*y,gyy=2-12*y+12*y*y;
    return {{32*f*gy,-32*fx*g},{32*fx*gy,-32*fxx*g},{32*f*gyy,-32*fx*gy}};
}
void algebra(){
    double maxValue=0,maxJet=0,maxMoments=0,maxTrace=0,maxJoin=0,maxWork=0;
    for(int n:{8,16}){
        QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
        const auto u=space.interpolateVelocity(mesh,domain,[](Vector2D p){const auto v=exact(p).value;return Eigen::Vector2d(v.x,v.y);});
        const auto spline=Spline::build(mesh,domain,space,u,0);const auto recovered=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=spline.firstJet(p).value;return Eigen::Vector2d(v.x,v.y);});
        maxMoments=std::max(maxMoments,(recovered-u).norm()/u.norm());const auto flux=spline.normalFlux();const auto& rule=Quantum::Geometry::Parametric::gaussEight();
        for(const auto face:domain.faces(mesh)){
            const auto& info=space.faces()[face.nativeRow];const auto chart=flux.faceChart(face.nativeRow);double moments[2]{};
            for(int half=0;half<2;half++)for(int k=0;k<8;k++){
                const double t=.5*(1+rule.x[k]);const auto a=half?chart.midpoint:chart.begin,b=half?chart.end:chart.midpoint;const Vector2D p=a+(b-a)*t;
                const double density=double(flux.normalFluxDensity(face.nativeRow,half,false,t));const auto sample=spline.firstJet(p);const double normal=face.horizontal?sample.value.y:sample.value.x;
                maxTrace=std::max(maxTrace,std::abs(2*density/face.area-normal));
                const double theta=((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;
                moments[0]+=.5*rule.w[k]*density/face.area;moments[1]+=.5*rule.w[k]*density*theta/face.area;
            }
            maxTrace=std::max({maxTrace,std::abs(moments[0]-u[info.first]),std::abs(moments[1]-u[info.first+1])});
            if(!face.outwardSign){const auto p=chart.midpoint;const double delta=1e-8/n;const Vector2D a=face.horizontal?Vector2D{p.x,p.y-delta}:Vector2D{p.x-delta,p.y},b=face.horizontal?Vector2D{p.x,p.y+delta}:Vector2D{p.x+delta,p.y};
                const auto va=spline.firstJet(a),vb=spline.firstJet(b),ea=exact(a),eb=exact(b);
                maxJoin=std::max({maxJoin,std::hypot((va.dx.x-vb.dx.x)-(ea.dx.x-eb.dx.x),(va.dx.y-vb.dx.y)-(ea.dx.y-eb.dx.y)),std::hypot((va.dy.x-vb.dy.x)-(ea.dy.x-eb.dy.x),(va.dy.y-vb.dy.y)-(ea.dy.y-eb.dy.y))});}
        }
        std::vector<Spline::PointForce> forces;double work=0;
        for(int j=0;j<31;j++)for(int i=0;i<29;i++){
            const Vector2D p{(i+.37)/29,(j+.29)/31},force{std::sin(i+.2*j),std::cos(j-.3*i)};const auto got=spline.firstJet(p),ref=exact(p);
            maxValue=std::max(maxValue,std::hypot(got.value.x-ref.value.x,got.value.y-ref.value.y));
            maxJet=std::max({maxJet,std::hypot(got.dx.x-ref.dx.x,got.dx.y-ref.dx.y),std::hypot(got.dy.x-ref.dy.x,got.dy.y-ref.dy.y)});
            require(got.dx.x+got.dy.y==0,"streamfunction velocity is not identically divergence free");forces.push_back({p,force});work+=force.x*got.value.x+force.y*got.value.y;
        }
        maxWork=std::max(maxWork,std::abs(u.dot(spline.transposePointForces(forces))-work)/std::max(1.,std::abs(work)));
        // A non-polynomial no-slip streamfunction checks moment preservation
        // beyond the exactly reproduced tensor-quartic flow.
        const auto nonlinear=space.interpolateVelocity(mesh,domain,[](Vector2D p){const auto one=[](double x){const double f=x*x*(1-x)*(1-x),d=2*x-6*x*x+4*x*x*x;return std::array<double,2>{f*(1+.2*std::sin(5*x)),d*(1+.2*std::sin(5*x))+f*std::cos(5*x)};};const auto a=one(p.x),b=one(p.y);return Eigen::Vector2d(32*a[0]*b[1],-32*a[1]*b[0]);});
        const auto other=spline.withVelocity(mesh,domain,space,nonlinear,.1);const auto again=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=other.firstJet(p).value;return Eigen::Vector2d(v.x,v.y);});
        maxMoments=std::max(maxMoments,(again-nonlinear).norm()/nonlinear.norm());
    }
    std::cout<<std::setprecision(17)<<"MOMENT_C1_SPLINE value="<<maxValue<<" jet="<<maxJet<<" moments="<<maxMoments<<" normal_trace="<<maxTrace<<" join="<<maxJoin<<" virtual_work="<<maxWork<<'\n';
    require(maxValue<2e-11&&maxJet<2e-9&&maxMoments<2e-11&&maxTrace<2e-11&&maxJoin<2e-9&&maxWork<2e-11,"moment-preserving C1 velocity or transpose failed");
}
void trajectory(){
    constexpr int n=8,steps=10;constexpr double dt=.04;
    QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const auto initial=space.interpolateVelocity(mesh,domain,[](Vector2D p){const auto v=exact(p).value;return Eigen::Vector2d(v.x,v.y);});
    auto stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),.02));
    auto hydro=Flow::create(mesh,domain,space,std::move(stress),1.3,initial);
    auto spline=*hydro.state().physicalVelocity;
    auto interface=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.43,.57},{.13,.13},24));
    const auto first=interface.parametricGeometry().geometry();const double initialVolume=first.area();double maxVolume=0,maxFraction=0,maxMotion=0,maxMassClosure=0;
    std::ofstream output,history;if(const char* prefix=std::getenv("QUANTUM_MOMENT_C1_OUTPUT")){
        output.open(std::string(prefix)+"_curves.csv");output<<std::setprecision(17)<<"step,time,parameter,x,y\n";
        history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)<<"step,time,volume_drift,fraction_error,kinetic_energy\n";}
    const auto save=[&](int step){if(output)for(int k=0;k<=512;k++){const double s=double(k)/512;const auto p=interface.parametricGeometry().geometry().sample(s);output<<step<<','<<step*dt<<','<<s<<','<<p.x<<','<<p.y<<'\n';}};save(0);
    for(int step=1;step<=steps;step++){
        HydroStepTime time;time.index=step;time.begin=hydro.state().time;time.dt=dt;
        auto load=hydro.makeLoad(mesh,domain,Eigen::VectorXd::Zero(initial.size()),time.begin+.5*dt);auto candidate=hydro.prepareNavierStokesStep(mesh,domain,time,load);
        auto nextSpline=*candidate.state().physicalVelocity;
        const auto map=HydroInterfaceTransport::timeLinearParametricFlow(spline,nextSpline,dt);
        const auto begin=spline.normalFlux(),end=nextSpline.normalFlux();
        const InterfaceStepTime interval{step,time.begin,time.begin+dt,dt};InterfaceParametricProjection::Options options;
        auto moved=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,interface,interval,map,begin,end,options);
        const auto mass=Quantum::Hydrodynamics::HydroCompatibleMomentTransport::phaseMassFlux(mesh,domain,space,interface,moved.state,interval,1000,1);
        maxMassClosure=std::max(maxMassClosure,mass.maximumRelativeClosure);
        hydro.accept(mesh,domain,std::move(candidate));interface=std::move(moved.state);spline=std::move(nextSpline);
        const auto& state=interface.parametricGeometry();double volume=0,fraction=0;const auto geometric=state.geometricVolumes();
        for(int k=0;k<int(state.inventory().size());k++){volume+=state.inventory()[k];fraction=std::max(fraction,std::abs(geometric[k]-state.inventory()[k])*n*n);}
        maxVolume=std::max(maxVolume,std::abs(volume-initialVolume));maxFraction=std::max(maxFraction,fraction);
        const double energy=hydro.kineticEnergy();
        if(history)history<<step<<','<<hydro.state().time<<','<<volume-initialVolume<<','<<fraction<<','<<energy<<'\n';
        if(step==5||step==steps)save(step);
        std::cout<<std::setprecision(17)<<"MOMENT_C1_INTERFACE_STEP step="<<step<<" time="<<hydro.state().time<<" volume_drift="<<volume-initialVolume<<" fraction_error="<<fraction<<" energy="<<energy<<" moment_mass_closure="<<mass.maximumRelativeClosure<<" moment_quadrature_error="<<moved.flux.momentQuadratureError<<std::endl;
    }
    for(int k=0;k<512;k++){const double s=double(k)/512;const auto a=first.sample(s),b=interface.parametricGeometry().geometry().sample(s);maxMotion=std::max(maxMotion,std::hypot(a.x-b.x,a.y-b.y));}
    std::cout<<std::setprecision(17)<<"MOMENT_C1_INTERFACE_COMPLETE steps="<<steps<<" time="<<hydro.state().time<<" volume_drift="<<maxVolume<<" fraction_error="<<maxFraction<<" maximum_material_displacement="<<maxMotion<<" moment_mass_closure="<<maxMassClosure<<" passive_interface=1 variable_mass=0 capillary_feedback=0\n";
    require(maxVolume<2e-13&&maxFraction<1e-10&&maxMotion>.005,"moment C1 interface trajectory incomplete or inconsistent");
}
// Replay the actual flow snapshots without repeating interface space-time
// quadrature. A failure prints the retained polynomial, not rounded samples.
void fluxRoots(bool integralMode=false){
    QuadtreeMesh mesh(input(8));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const auto initial=space.interpolateVelocity(mesh,domain,[](Vector2D p){const auto v=exact(p).value;return Eigen::Vector2d(v.x,v.y);});
    auto stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),.02));
    auto hydro=Flow::create(mesh,domain,space,std::move(stress),1.3,initial);auto spline=*hydro.state().physicalVelocity;
    for(int step=1;step<=2;step++){
        HydroStepTime time;time.index=step;time.begin=hydro.state().time;time.dt=.04;
        const auto load=hydro.makeLoad(mesh,domain,Eigen::VectorXd::Zero(initial.size()),time.begin+.5*time.dt);
        auto candidate=hydro.prepareNavierStokesStep(mesh,domain,time,load);auto next=*candidate.state().physicalVelocity;
        const auto begin=spline.normalFlux(),end=next.normalFlux();
        for(const auto face:domain.faces(mesh))for(int half=0;half<2;half++){
            std::vector<double> times;for(int k=0;k<=8;k++)times.push_back(k/8.);
            for(double s:{0.,1.}){
                const auto a=Quantum::Geometry::evaluateUnitPolynomial(begin.normalFluxPolynomial(face.nativeRow,half),s);
                const auto b=Quantum::Geometry::evaluateUnitPolynomial(end.normalFluxPolynomial(face.nativeRow,half),s);
                if((a<0&&b>0)||(a>0&&b<0)){const double root=double(a/(a-b));if(root>0&&root<1){times.push_back(root);times.push_back(std::nextafter(root,0.));times.push_back(std::nextafter(root,1.));}}
            }
            for(const double theta:times){auto terms=begin.normalFluxPolynomial(face.nativeRow,half);for(auto& term:terms)term.weight*=1-theta;
            for(auto term:end.normalFluxPolynomial(face.nativeRow,half)){term.weight*=theta;terms.push_back(term);}
            if(integralMode){
                double scale=0;for(const auto& term:terms)for(double c:term.coefficients)scale+=std::abs(double(term.factor)*term.weight*c);
                const double tolerance=32*std::numeric_limits<double>::epsilon()*scale;
                const auto positive=Quantum::Geometry::integratePositiveUnitPolynomial(terms,0,1,tolerance);
                auto negative=terms;for(auto& term:negative)term.weight=-term.weight;
                const auto reverse=Quantum::Geometry::integratePositiveUnitPolynomial(negative,0,1,tolerance);
                const auto signedValue=Quantum::Geometry::integrateUnitPolynomial(terms,0,1,0,1);
                require(std::abs(positive.value-reverse.value-signedValue)<=positive.absoluteError+reverse.absoluteError+8*std::numeric_limits<double>::epsilon()*std::abs(signedValue),
                    "bounded positive flux does not retain the signed physical trace");
                continue;
            }
            try{Quantum::Geometry::isolateUnitPolynomial(terms);}catch(const std::exception&){
                std::cout<<std::setprecision(17)<<"MOMENT_C1_FLUX_ROOT_FAILURE step="<<step<<" face="<<face.nativeRow<<" horizontal="<<face.horizontal<<" center_x="<<face.center.x<<" center_y="<<face.center.y<<" half="<<half<<" theta="<<theta<<std::endl;
                auto derivative=terms;for(auto& term:derivative){for(std::size_t k=0;k+1<term.coefficients.size();k++)term.coefficients[k]=(k+1)*term.coefficients[k+1];term.coefficients.back()=0;}
                const auto atEnd=Quantum::Geometry::evaluateUnitPolynomial(terms,1),slope=Quantum::Geometry::evaluateUnitPolynomial(derivative,1);
                std::cout<<"MOMENT_C1_FLUX_ENDPOINT value="<<atEnd<<" derivative="<<slope<<" linearized_offset="<<(slope==0?0:atEnd/slope)<<std::endl;
                for(const auto& term:terms){std::cout<<"TERM factor="<<term.factor<<" weight="<<term.weight;for(double c:term.coefficients)std::cout<<" "<<c;std::cout<<std::endl;}throw;
            }
            }
        }
        hydro.accept(mesh,domain,std::move(candidate));spline=std::move(next);
    }
    std::cout<<"MOMENT_C1_FLUX_REPLAY steps=2 temporal_samples=9 endpoint_crossing_neighbors=1 integral_mode="<<integralMode<<" passed=1\n";
}
void positiveIntegrals(){
    using namespace Quantum::Geometry;UnitPolynomialTerm p;p.factor=1;p.coefficients[0]=-.25;p.coefficients[1]=1;
    // A binary-exact full-face flux must not spend an artificial relative
    // roundoff allowance when AMR subdivides the shared plane budget.
    UnitPolynomialTerm constant;constant.factor=std::ldexp(1.,-7);constant.coefficients[0]=1;
    const auto dyadic=integratePositiveUnitPolynomial({constant},0,1,0);
    require(dyadic.value==std::ldexp(1.L,-7)&&dyadic.absoluteError==0,"exact constant flux has a fictitious integration error");
    UnitPolynomialTerm linear;linear.factor=1;linear.coefficients[1]=1;
    const auto affine=integratePositiveTimeLinearUnitPolynomial({constant},{linear},0,1,0);
    require(affine.value==.25390625L&&affine.absoluteError==0,"exact time-linear positive flux has a fictitious integration error");
    UnitPolynomialTerm quadratic;quadratic.factor=1;quadratic.coefficients[2]=1;
    const auto third=integratePositiveUnitPolynomial({quadratic},0,1,1e-16L);
    const long double exactResidual=std::ldexp(1.L,-std::numeric_limits<long double>::digits-1)/3;
    require(third.value==1.L/3&&third.absoluteError>=exactResidual&&third.absoluteError<=1.01L*exactResidual,
        "nonrepresentable integral rounding is not enclosed tightly");
    const auto positive=integratePositiveUnitPolynomial({p},0,1,1e-14L);
    auto q=p;q.weight=-1;const auto reversal=integratePositiveTimeLinearUnitPolynomial({p},{q},0,1,1e-14L);
    UnitPolynomialTerm eighth;eighth.factor=1;eighth.coefficients[8]=1;
    const auto high=integratePositiveUnitPolynomial({eighth},0,1,1e-14L);
    require(std::abs(positive.value-.28125L)<=positive.absoluteError+1e-17L
        &&std::abs(reversal.value-.078125L)<=reversal.absoluteError+1e-17L
        &&std::abs(high.value-1.L/9)<=high.absoluteError,"positive/time-reversing/eighth-degree integral differs from analytic values");
    // A narrow positive island invisible to a fixed Gauss rule. Keep its
    // small constant as a separate term, so the input does not round it away.
    UnitPolynomialTerm square;square.factor=1;square.coefficients[0]=-.140625;square.coefficients[1]=.75;square.coefficients[2]=-1;
    UnitPolynomialTerm small;small.factor=std::ldexp(1.,-40);small.coefficients[0]=1;
    const auto island=integratePositiveUnitPolynomial({square,small},0,1,1e-21L);
    const long double expected=std::ldexp(4.L/3,-60);
    require(island.value>0&&std::abs(island.value-expected)<=island.absoluteError+1e-30L,"thin sign-changing flux island was omitted");
    // The formerly unrepresentable near-endpoint sign interval. Its true
    // positive volume is tiny; the INTEGRAL uncertainty, not root precision,
    // decides acceptance. No input coefficient is removed.
    UnitPolynomialTerm endpoint;endpoint.factor=1;endpoint.coefficients[0]=-1;endpoint.coefficients[1]=1;
    small.factor=1e-29;const auto terminal=integratePositiveUnitPolynomial({endpoint,small},0,1,1e-18L);
    require(std::abs(terminal.value-.5e-58L)<=terminal.absoluteError+1e-70L,"terminal positive-volume bound excludes the actual thin interval");
    // A full degree-31 product with an analytic beta-function integral.
    // Exact cancellation of large extra terms exercises the shared trace
    // evaluation; partition/weight sums alone could hide swapped entries.
    std::array<long double,21> products{},bounds{};int shared=0;
    for(int half=0;half<2;half++){
        UnitPolynomialTerm power;power.factor=std::ldexp(1.,-11);double choose=1;
        for(int j=0;j<=11;j++){power.coefficients[j]=half?choose:(j==11?1:0);choose=choose*(11-j)/(j+1);}
        auto large=power;large.factor=std::ldexp(1.,30);auto negative=large;negative.weight=-1;
        const std::vector<UnitPolynomialTerm> terms{large,power,negative};
        for(const auto& interval:std::array<std::array<double,2>,3>{{{0,1},{.125,.875},{.5,std::nextafter(.5,1.)}}}){
            const auto all=integrateUnitPolynomialBernsteinBasis(terms,20,half,interval[0],interval[1]);
            for(int k=0;k<=20;k++){
                const auto one=integrateUnitPolynomialBernstein(terms,20,k,half,interval[0],interval[1]);
                require(all[k].value==one.value&&all[k].absoluteError==one.absoluteError,"shared Bernstein rule changed a scalar value or error estimate");++shared;
                if(interval[0]==0){products[k]+=.5L*all[k].value;bounds[k]+=.5L*all[k].absoluteError;}
            }
        }
    }
    long double productError=0;
    for(int k=0;k<=20;k++){long double exact=1.L/21;for(int j=1;j<=11;j++)exact*=static_cast<long double>(k+j)/(21+j);
        productError=std::max(productError,std::abs(products[k]-exact));
        require(std::abs(products[k]-exact)<=bounds[k]+2e-17L,"shared complete products differ from analytic beta integrals");}
    std::cout<<std::setprecision(17)<<"MOMENT_SHARED_PRODUCTS scalar_comparisons="<<shared<<" bitwise_equal=1 maximum_analytic_error="<<double(productError)<<'\n';
    std::cout<<std::setprecision(17)<<"MOMENT_POSITIVE_FLUX exact_constant="<<double(dyadic.value)<<" constant_error="<<double(dyadic.absoluteError)
        <<" third_error="<<double(third.absoluteError)<<" analytic="<<double(positive.value)<<" reversal="<<double(reversal.value)<<" eighth="<<double(high.value)
        <<" thin_island="<<double(island.value)<<" thin_error_bound="<<double(island.absoluteError)<<" endpoint_value="<<double(terminal.value)<<" endpoint_error_bound="<<double(terminal.absoluteError)<<'\n';
}
}
void hydroMomentSplineContract(){algebra();}
void hydroMomentPortSplineContract(){
    using Boundary=Quantum::Hydrodynamics::HydroCompatibleMomentBoundary;
    using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
    std::ofstream samples;
    if(const char* output=std::getenv("QUANTUM_MOMENT_PORT_SPLINE_OUTPUT")){
        samples.open(output);require(bool(samples),"cannot write port C1 samples");
        samples<<std::setprecision(17)<<"rz,mixed,step,time,x,y,u,v\n";
    }
    for(bool rz:{false,true})for(bool mixed:{false,true}){
        QuadtreeMesh mesh(input(4));if(mixed){DictionaryMesh::RefineBox box;box.min={0,.5};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});}
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto space=Space::build(mesh,domain);
        const std::vector<HydroBoundaryPatch> ports{{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,1,-.1},
            {"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
        auto flow=Flow::create(mesh,domain,space,Stress::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),.02)),1.,Eigen::VectorXd::Zero(space.mass().rows()),0,ports);
        auto spline=*flow.state().physicalVelocity;
        double moment=0,trace=0,wall=0,join=0,div=0,adjoint=0,flowError=0,energyError=0,minimumOutlet=1,maximumInlet=-1;
        const auto& g=Quantum::Geometry::Parametric::gaussEight();
        for(int step=0;step<=3;step++){
            if(step){const HydroStepTime time{step,flow.state().time,.01,HydroTimePolicy::Fixed};
                const auto load=flow.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),time.begin+.5*time.dt);
                auto candidate=flow.prepareNavierStokesStep(mesh,domain,time,load);flow.accept(mesh,domain,std::move(candidate));
                const auto next=*flow.state().physicalVelocity;require(next.sameLayout(spline),"port C1 space changed with velocity");spline=next;
            }
            const auto& u=flow.state().velocity;const auto recovered=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=spline.firstJet(p).value;return Eigen::Vector2d{v.x,v.y};});
            moment=std::max(moment,(recovered-u).norm()/u.norm());const auto flux=spline.normalFlux();std::vector<Spline::PointForce> forces;double work=0,inflow=0,outflow=0,kinetic=0;
            for(const auto f:domain.faces(mesh)){
                const auto chart=flux.faceChart(f.nativeRow);const double length=std::hypot(chart.end.x-chart.begin.x,chart.end.y-chart.begin.y);
                for(int half=0;half<2;half++)for(int k=0;k<8;k++){
                    const double t=.5+.5*g.x[k];const auto p=chart.begin+(chart.end-chart.begin)*(.5*(half+t));const auto v=spline.firstJet(p);
                    const double density=double(flux.normalFluxDensity(f.nativeRow,half,false,t)),expected=.5*length*(rz?2*std::acos(-1.)*p.x:1)*(f.horizontal?v.value.y:v.value.x);
                    if(f.outwardSign){const double normal=f.outwardSign*(f.horizontal?v.value.y:v.value.x);const auto kind=flow.boundary().faces()[f.nativeRow].kind;
                        if(kind==Boundary::Kind::PressureOutlet)minimumOutlet=std::min(minimumOutlet,normal);
                        if(kind==Boundary::Kind::NormalVelocity)maximumInlet=std::max(maximumInlet,normal);}
                    trace=std::max(trace,std::abs(density-expected));if(f.outwardSign&&f.area>0&&flow.boundary().faces()[f.nativeRow].kind==Boundary::Kind::Wall)wall=std::max(wall,std::hypot(v.value.x,v.value.y));
                    if(rz&&p.x==0)require(v.value.x==0&&v.dx.y==0,"port C1 axis lost its even regular limit");
                    if(f.outwardSign){const Vector2D force{.13+.07*p.y,.09-.03*p.x},dx{.03,-.07},dy{-.02,.05};forces.push_back({p,force,dx,dy});
                        work+=force.x*v.value.x+force.y*v.value.y+dx.x*v.dx.x+dx.y*v.dx.y+dy.x*v.dy.x+dy.y*v.dy.y;}
                }
                const double q=f.outwardSign*double(flux.integratedNormalFlux(f.nativeRow,0,0,1)+flux.integratedNormalFlux(f.nativeRow,1,0,1));
                inflow+=std::max(0.,-q);outflow+=std::max(0.,q);
                if(!f.outwardSign){const auto p=chart.begin+(chart.end-chart.begin)*.371;const double eps=1e-11;
                    const auto a=spline.firstJet(p+(f.horizontal?Vector2D{0,-eps}:Vector2D{-eps,0})),b=spline.firstJet(p+(f.horizontal?Vector2D{0,eps}:Vector2D{eps,0}));
                    join=std::max({join,std::hypot(a.value.x-b.value.x,a.value.y-b.value.y),std::hypot(a.dx.x-b.dx.x,a.dx.y-b.dx.y),std::hypot(a.dy.x-b.dy.x,a.dy.y-b.dy.y)});
                }
            }
            HydroCompleteObservation::volume(mesh,domain,[&](Vector2D p,double dv){const auto v=spline.firstJet(p);
                div=std::max(div,std::abs(v.dx.x+v.dy.y+(rz?v.value.x/p.x:0)));
                kinetic+=dv*(v.value.x*v.value.x+v.value.y*v.value.y);});
            adjoint=std::max(adjoint,std::abs(spline.coordinates().dot(spline.transposePointForces(forces))-work)/std::max(1.,std::abs(work)));
            const double expected=.1*(rz?std::acos(-1.):1);flowError=std::max({flowError,std::abs(inflow-expected)/expected,std::abs(outflow-expected)/expected});
            energyError=std::max(energyError,std::abs(kinetic/(2*flow.kineticEnergy())-1));
            if(samples&&step==3)for(int j=0;j<=48;j++)for(int i=0;i<=48;i++){const Vector2D p{i/48.,j/48.};const auto v=spline.firstJet(p).value;samples<<rz<<','<<mixed<<','<<step<<','<<flow.state().time<<','<<p.x<<','<<p.y<<','<<v.x<<','<<v.y<<'\n';}
        }
        // Wall/port identity belongs to the reconstruction, even for a zero
        // state. A foreign boundary owner cannot silently become a closed one.
        const Eigen::VectorXd zero=Eigen::VectorXd::Zero(space.mass().rows());const auto closed=Spline::build(mesh,domain,space,zero,0).nativeVelocity();
        require(!closed.sameLayout(spline)&&spline.withCoordinates(Eigen::VectorXd::Zero(spline.coordinates().size()),0).sameLayout(spline),"port partition lost at zero flow");
        Boundary unbound;bool rejected=false;try{Spline::build(mesh,domain,space,zero,0,&unbound);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"port C1 accepted foreign boundary");
        std::cout<<std::setprecision(17)<<"MOMENT_PORT_C1 metric="<<(rz?"RZ":"XY")<<" mixed="<<mixed<<" cells="<<mesh.cellIndexNum()<<" steps=3 time="<<flow.state().time
            <<" moments="<<moment<<" trace="<<trace<<" wall="<<wall<<" join="<<join<<" divergence="<<div<<" virtual_work="<<adjoint<<" flow_relative="<<flowError<<" kinetic_relative="<<energyError<<" sampled_min_outlet="<<minimumOutlet<<" sampled_max_inlet="<<maximumInlet<<std::endl;
        require(moment<2e-9&&trace<2e-10&&wall<2e-10&&join<2e-7&&div<2e-11&&adjoint<2e-9&&flowError<2e-9&&energyError<2e-9,"open complete C1 field lost physical moments, trace, energy or work");
    }
    // Partial inlet and wall junctions, with all solenoidal moment directions
    // on a small mixed mesh. No prescribed smooth profile hides missing modes.
    for(bool rz:{false,true}){
        QuadtreeMesh mesh(input(2));DictionaryMesh::RefineBox box;box.min={0,.5};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto space=Space::build(mesh,domain);
        const auto stress=Stress::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),.02));
        const auto boundary=Boundary::build(mesh,domain,space,stress,{{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,.5,-.1},{"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,0}});
        const int nv=int(space.mass().rows()),nd=int(space.divergence().rows());std::vector<int> walls;
        for(const auto f:domain.faces(mesh))if(f.outwardSign&&boundary.faces()[f.nativeRow].kind==Boundary::Kind::Wall&&space.faces()[f.nativeRow].first>=0){walls.push_back(space.faces()[f.nativeRow].first);walls.push_back(walls.back()+1);}
        Eigen::MatrixXd constraints=Eigen::MatrixXd::Zero(nd+walls.size(),nv);constraints.topRows(nd)=space.divergence();
        for(int k=0;k<nd;k++)constraints.row(k)/=constraints.row(k).norm();for(int k=0;k<int(walls.size());k++)constraints(nd+k,walls[k])=1;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(constraints,Eigen::ComputeFullV);int rank=0;for(int k=0;k<svd.singularValues().size();k++)if(svd.singularValues()[k]>1e-12*svd.singularValues()[0])rank++;
        Eigen::MatrixXd basis=svd.matrixV().rightCols(nv-rank);for(int row:walls)basis.row(row).setZero();require((constraints*basis).norm()<2e-11,"partial-port basis is not solenoidal");
        const auto spline=Spline::build(mesh,domain,space,Eigen::VectorXd::Zero(nv),0,&boundary);double moment=0,rx[12],rw[12];kineticQuadrature(rx,rw);
        Eigen::MatrixXd values(2*mesh.cellIndexNum()*12*12,basis.cols());
        for(int col=0;col<basis.cols();col++){
            const auto field=spline.withVelocity(mesh,domain,space,basis.col(col),0);int row=0;
            const auto recovered=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=field.firstJet(p).value;return Eigen::Vector2d{v.x,v.y};});moment=std::max(moment,(recovered-basis.col(col)).norm());
            for(int index:domain.cellIndexByRow()){
                const auto& c=mesh.cells(mesh.cellIndexToLabel[index]);const auto center=c.centerCoord(),h=c.size();
                for(int a=0;a<12;a++)for(int b=0;b<12;b++){const Vector2D p{center.x+.5*h.x*rx[a],center.y+.5*h.y*rx[b]};const auto v=field.firstJet(p).value;
                    const double root=std::sqrt(.25*h.x*h.y*rw[a]*rw[b]*(rz?2*std::acos(-1.)*p.x:1));values(row++,col)=root*v.x;values(row++,col)=root*v.y;}
            }
        }
        const Eigen::MatrixXd actual=values.transpose()*values,owned=basis.transpose()*space.mass()*basis;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen(actual,owned,Eigen::EigenvaluesOnly);
        require(eigen.info()==Eigen::Success&&eigen.eigenvalues().allFinite(),"partial-port kinetic metric failed");
        std::cout<<"MOMENT_PORT_C1_SPACE metric="<<(rz?"RZ":"XY")<<" cells="<<mesh.cellIndexNum()<<" solenoidal_dofs="<<basis.cols()<<" moments="<<moment<<" minimum="<<eigen.eigenvalues()[0]<<" maximum="<<eigen.eigenvalues().tail(1)[0]<<std::endl;
        require(moment<2e-9&&eigen.eigenvalues()[0]>0&&eigen.eigenvalues().tail(1)[0]<32,"partial-port reconstruction lacks bounded moment lifting");
    }
}
void hydroMomentSupplyContract(){
    namespace P=Quantum::Geometry::Parametric;using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;using Traction=Quantum::Hydrodynamics::HydroCompatibleMomentTraction;
    std::ofstream history,curves;if(const char* prefix=std::getenv("QUANTUM_MOMENT_SUPPLY_OUTPUT")){
        history.open(std::string(prefix)+"_history.csv");curves.open(std::string(prefix)+"_curves.csv");require(history&&curves,"cannot write coupled supplied interface");
        history<<std::setprecision(17)<<"sigma,step,time,volume,volume_error,mass_error,fraction,kinetic,surface,energy_defect,surface_work_defect,boundary_traction,boundary_reaction,boundary_kinetic,inflow_dissipation,material_work,equation,iterations\n";
        curves<<std::setprecision(17)<<"sigma,step,time,s,r,z\n";
    }
    for(double sigma:{0.,.1}){
        QuadtreeMesh mesh(input(4));DictionaryMesh::RefineBox box;box.min={0,.5};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
        auto interface=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricSupplyCap(.5,.3,8));
        std::vector<HydroBoundaryPatch> ports{{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,.5,-.1},{"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
        ports[0].incomingPhase=-1;ports[1].incomingPhase=1;ports[1].incomingVelocity=Vector2D{0,0};
        auto flow=Flow::createTwoPhase(mesh,domain,space,std::move(interface),Material::Properties{5,1,.02,.01},Eigen::VectorXd::Zero(space.mass().rows()),sigma,ports);
        const auto m0=flow.kineticMass();const double k0=flow.kineticEnergy(),v0=.5*std::acos(-1.)*.25*.3;
        // Legacy moment initialization norm, not complete-field kinetic
        // stationarity. Its boundary initialization is tested separately.
        std::vector<int> free;for(int i=0;i<m0.rows();i++)if(!flow.boundary().fixed()[i])free.push_back(i);
        const Eigen::VectorXd weighted=m0*flow.state().velocity;Eigen::VectorXd g(free.size());Eigen::MatrixXd d(space.divergence().rows(),free.size());
        for(int i=0;i<int(free.size());i++){g[i]=weighted[free[i]];d.col(i)=Eigen::VectorXd(space.divergence().col(free[i]));}
        const Eigen::VectorXd multiplier=d.transpose().completeOrthogonalDecomposition().solve(g);
        const double initialStationarity=(d.transpose()*multiplier-g).norm()/weighted.norm();require(initialStationarity<2e-10,"two-phase supplied initialization did not use cut-material mass");
        const auto spline=*flow.state().physicalVelocity;
        const auto complement=Quantum::Hydrodynamics::HydroCompatibleMomentComplement::build(mesh,domain,space,flow.boundary().fixed());
        const auto force=Traction::build(mesh,domain,space,spline,complement,flow.material(),[](const P::SurfaceSample& q){return Traction::Components{.1*q.curvature,0};},.1);
        long double rate=0;const auto& curve=flow.interfaceState().parametricGeometry().geometry();double nodes[12],weights[12];kineticQuadrature(nodes,weights);
        for(const auto& arc:curve.arcs){std::vector<double> cuts{0,1};
            for(const auto face:domain.faces(mesh)){const int axis=face.horizontal?1:0;auto equation=arc.coordinate(axis);equation[0]-=axis?face.center.y:face.center.x;
                for(double t:P::roots(equation,0,1))cuts.push_back(t);}
            std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
            for(std::size_t part=1;part<cuts.size();part++)for(int panel=0;panel<16;panel++)for(int q=0;q<12;q++){
            const double width=(cuts[part]-cuts[part-1])/16,t=cuts[part-1]+width*(panel+.5*(1+nodes[q])),w=.5*weights[q]*width;const auto surface=P::sampleSurfaceRZ(arc,t);
            const auto v=spline.firstJet({surface.point.x,surface.point.y});const auto tangent=surface.derivative;
            const P::Point change{v.dx.x*tangent.x+v.dy.x*tangent.y,v.dx.y*tangent.x+v.dy.y*tangent.y};
            rate+=w*2*std::acos(-1.)*(v.value.x*surface.speed+surface.point.x*P::dot(tangent,change)/surface.speed);
        }}
        const double work=spline.coordinates().dot(force.variation.head(spline.coordinates().size()))
            +flow.state().velocity.dot(space.pressureLoad()*force.pressurePotential),workDefect=std::abs(work+.1*double(rate));
        std::cout<<std::setprecision(17)<<"MOMENT_SUPPLY_WORK sigma="<<sigma<<" initial_stationarity="<<initialStationarity<<" force_work="<<work<<" area_rate="<<double(rate)<<" defect="<<workDefect<<std::endl;
        require(workDefect<2e-10,"attached capillary work omitted the wet boundary term");
        const auto save=[&](int step){if(curves)for(int k=0;k<=128;k++){const auto p=flow.interfaceState().parametricGeometry().geometry().sample(k/128.);curves<<sigma<<','<<step<<','<<flow.state().time<<','<<k/128.<<','<<p.x<<','<<p.y<<'\n';}};save(0);
        Flow::Options options;options.interfaceProjection.materialSegments=8;
        double maxVolume=0,maxMass=0,maxFraction=0,maxEquation=0,maxEnergy=0,maxMetric=0,maxPortWork=0,totalSurfaceDefect=0;
        for(int step=1;step<=3;step++){
            constexpr double dt=.002;const double begin=flow.state().time;const auto old=flow.interfaceState().parametricGeometry().geometryHandle();
            const auto load=flow.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),begin+.5*dt);
            auto candidate=flow.prepareNavierStokesStep(mesh,domain,{step,begin,dt,HydroTimePolicy::Fixed},load,options);const auto report=candidate.report();
            require(flow.state().time==begin&&flow.interfaceState().parametricGeometry().geometryHandle()==old,"supplied midpoint preparation changed accepted state");flow.accept(mesh,domain,std::move(candidate));
            const auto& state=flow.interfaceState();const auto& inventory=state.parametricGeometry().inventory();const double volume=std::accumulate(inventory.begin(),inventory.end(),0.),expected=v0+flow.state().time*std::acos(-1.)*.25*.1;
            double mass=0;for(const auto& cell:flow.material().cells())mass+=cell.density(0,0);const double massError=std::abs(mass-(std::acos(-1.)+4*expected));
            maxVolume=std::max(maxVolume,std::abs(volume-expected));maxMass=std::max(maxMass,massError);maxFraction=std::max(maxFraction,state.parametricGeometry().report().maximumInventoryFractionError);
            maxEquation=std::max(maxEquation,report.equationResidual);maxMetric=std::max(maxMetric,report.massMetricDefect);maxEnergy=std::max(maxEnergy,std::abs(report.energyDefect)/k0);
            maxPortWork=std::max(maxPortWork,std::abs(report.materialMetricWork-report.boundaryKineticFlux-report.boundaryInflowDissipation)/k0);totalSurfaceDefect+=report.surfaceWorkDefect;
            if(history)history<<sigma<<','<<step<<','<<flow.state().time<<','<<volume<<','<<std::abs(volume-expected)<<','<<massError<<','<<state.parametricGeometry().report().maximumInventoryFractionError<<','<<report.kineticAfter<<','<<report.surfaceEnergyAfter<<','<<report.energyDefect<<','<<report.surfaceWorkDefect<<','<<report.boundaryTractionWork<<','<<report.boundaryReactionWork<<','<<report.boundaryKineticFlux<<','<<report.boundaryInflowDissipation<<','<<report.materialMetricWork<<','<<report.equationResidual<<','<<report.nonlinearIterations<<'\n';save(step);
            std::cout<<"MOMENT_SUPPLY_STEP sigma="<<sigma<<" step="<<step<<" volume_error="<<maxVolume<<" mass_error="<<maxMass<<" energy="<<maxEnergy<<" metric="<<maxMetric<<" boundary_work="<<maxPortWork<<" equation="<<maxEquation<<std::endl;
        }
        const double massChange=(flow.kineticMass()-m0).norm()/m0.norm();
        std::cout<<"MOMENT_SUPPLY sigma="<<sigma<<" time="<<flow.state().time<<" volume_error="<<maxVolume<<" mass_error="<<maxMass<<" fraction="<<maxFraction<<" equation="<<maxEquation<<" energy="<<maxEnergy<<" metric="<<maxMetric<<" boundary_work="<<maxPortWork<<" mass_change="<<massChange<<" surface_work_defect="<<totalSurfaceDefect<<std::endl;
        require(maxVolume<2e-12&&maxMass<2e-11&&maxFraction<3.1e-12&&maxEquation<1e-10&&maxEnergy<2e-9&&maxMetric<2e-10&&maxPortWork<2e-9&&massChange>1e-5,"active two-phase supplied midpoint failed");
    }
}
void hydroParametricSupplyContract(){
    namespace P=Quantum::Geometry::Parametric;using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
    constexpr double radius=.5,height=.3,speed=.3;
    std::ofstream history,curves;if(const char* prefix=std::getenv("QUANTUM_PARAMETRIC_SUPPLY_OUTPUT")){
        history.open(std::string(prefix)+"_history.csv");curves.open(std::string(prefix)+"_curves.csv");require(history&&curves,"cannot write supplied interface result");
        history<<std::setprecision(17)<<"native,step,time,volume,expected_volume,fraction,balance\n";curves<<std::setprecision(17)<<"native,step,time,s,r,z\n";
    }
    for(bool native:{false,true}){
        QuadtreeMesh mesh(input(4));DictionaryMesh::RefineBox box;box.min={0,.5};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
        auto initial=InterfacePrescribedData::parametricSupplyCap(radius,height,8);
        require(std::abs(initial.area()-2*radius*height/3)<2e-15,"phase closure is missing from supplied cap area");
        const auto surface=Quantum::Hydrodynamics::HydroCompatibleMomentTraction::measure(initial,true);
        const double expectedSurface=std::acos(-1.)*std::pow(radius,4)/(6*height*height)*(std::pow(1+4*height*height/(radius*radius),1.5)-1);
        require(std::abs(surface.length-expectedSurface)<2e-12,"inlet disk was included in free capillary area");
        auto state=InterfaceParametricMotion::initialize(mesh,domain,std::move(initial));const double initialVolume=.5*std::acos(-1.)*radius*radius*height;
        const auto liquid=[&]{const auto& v=state.parametricGeometry().inventory();return std::accumulate(v.begin(),v.end(),0.);};
        require(std::abs(liquid()-initialVolume)<2e-14,"supplied cap initial RZ measure failed");
        const std::vector<HydroBoundaryPatch> ports{{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,radius,-.1},{"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
        auto fluid=Flow::create(mesh,domain,space,Stress::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),.02)),1.,Eigen::VectorXd::Zero(space.mass().rows()),0,ports);
        auto spline=*fluid.state().physicalVelocity;
        using Terms=std::vector<Quantum::Geometry::UnitPolynomialTerm>;
        auto charts=std::make_shared<std::vector<InterfaceNormalFluxView::FaceChart>>(domain.nativeFaceCount());
        auto terms=std::make_shared<std::vector<std::array<Terms,2>>>(domain.nativeFaceCount());
        for(const auto f:domain.faces(mesh)){
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const Vector2D a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            (*charts)[f.nativeRow]={a,edge.midpointCoord(),b};if(!f.horizontal||a.x>=radius)continue;require(b.x<=radius,"prescribed supply trace crosses an unsplit rim");
            for(int half=0;half<2;half++){Quantum::Geometry::UnitPolynomialTerm p;const double lo=a.x+.5*half*(b.x-a.x),length=.5*(b.x-a.x);p.factor=2*std::acos(-1.)*speed*length;
                for(int degree:{1,3,5}){const double scale=degree==1?1:degree==3?-2/(radius*radius):1/std::pow(radius,4);
                    for(int k=0;k<=degree;k++){double choose=1;for(int j=1;j<=k;j++)choose*=double(degree+1-j)/j;p.coefficients[k]+=scale*choose*std::pow(lo,degree-k)*std::pow(length,k);}}
                (*terms)[f.nativeRow][half].push_back(p);
            }
        }
        InterfaceNormalFluxView::Queries queries;queries.retainedMomentDegree=2;queries.chart=[charts](int row){return charts->at(row);};
        queries.polynomial=[terms](int row,int half){return terms->at(row).at(half);};
        queries.integral=[terms](int row,int half,double a,double b,double c,double d){return Quantum::Geometry::integrateUnitPolynomial(terms->at(row).at(half),a,b,c,d);};
        queries.density=[terms](int row,int half,bool end,long double x){return Quantum::Geometry::evaluateUnitPolynomial(terms->at(row).at(half),double(end?1-x:x));};
        queries.phaseIntegral=[](int,int,const LevelSetMappedSegment&)->long double{throw std::invalid_argument("supplied cap requires parametric geometry");};
        const InterfaceNormalFluxView analyticFlux(std::move(queries));
        const InterfaceParametricTransport::Flow analytic=[=](Vector2D p,double t){const double q=std::max(0.,1-p.x*p.x/(radius*radius));
            const double u=speed*q*q,du=p.x<radius?-4*speed*p.x*q/(radius*radius):0;
            return InterfaceParametricTransport::MapJet{{p.x,p.y+t*u},{1,t*du},{0,1},{0,u},{0,t*u},0};};
        const auto save=[&](int step){if(curves)for(int k=0;k<=128;k++){const double s=k/128.;const auto p=state.parametricGeometry().geometry().sample(s);curves<<native<<','<<step<<','<<state.progress().endTime<<','<<s<<','<<p.x<<','<<p.y<<'\n';}};save(0);
        double maxVolume=0,maxFraction=0,maxPosition=0,maxCurvature=0;const int steps=native?3:5;const double dt=native?.002:.04;
        InterfaceParametricProjection::Options options;options.materialSegments=8;
        for(int step=1;step<=steps;step++){
            const double begin=state.progress().endTime,end=begin+dt;
            if(native){const auto load=fluid.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),begin+.5*dt);
                auto candidate=fluid.prepareNavierStokesStep(mesh,domain,{step,begin,dt,HydroTimePolicy::Fixed},load);fluid.accept(mesh,domain,std::move(candidate));
                spline=*fluid.state().midpointPhysicalVelocity;}
            const auto velocity=native?HydroInterfaceTransport::frozenParametricFlow(spline):analytic;const auto flux=native?spline.normalFlux():analyticFlux;
            const auto oldGeometry=state.parametricGeometry().geometryHandle();
            auto next=InterfaceParametricMotion::prepare(mesh,domain,state,{step,begin,end,dt},dt,velocity,flux,options);
            require(state.parametricGeometry().geometryHandle()==oldGeometry&&state.progress().endTime==begin,"supply preparation mutated accepted state");
            const double balance=next.inventoryBalance;state=std::move(next.state);const auto& geometry=state.parametricGeometry().geometry();
            require(geometry.boundary.fixed(0)&&geometry.boundary.axis(1)&&geometry.sample(0).x==radius&&geometry.sample(0).y==0&&geometry.sample(1).x==0,"supplied endpoint ownership lost");
            const double volume=liquid(),expected=initialVolume+end*std::acos(-1.)*radius*radius*(native?.1:speed/3);
            maxVolume=std::max(maxVolume,std::abs(volume-expected));maxFraction=std::max(maxFraction,state.parametricGeometry().report().maximumInventoryFractionError);
            if(!native)for(int k=0;k<=128;k++){
                const double s=k/128.,r=radius*(1-s),q=1-r*r/(radius*radius),z=height*q+speed*end*q*q;
                maxPosition=std::max(maxPosition,P::norm(geometry.sample(s)-P::Point{r,z}));
                const double dz=-2*height*r/(radius*radius)-4*speed*end*r*q/(radius*radius),ddz=-2*height/(radius*radius)-4*speed*end/(radius*radius)*(1-3*r*r/(radius*radius));
                const double curvature=r==0?-2*ddz:-ddz/std::pow(1+dz*dz,1.5)-dz/(r*std::sqrt(1+dz*dz));
                const auto& arc=geometry.arcs[geometry.index(s)];const auto got=P::sampleSurfaceRZ(arc,(s-arc.a)/(arc.b-arc.a),s==1);maxCurvature=std::max(maxCurvature,std::abs(got.curvature-curvature));
            }
            if(history)history<<native<<','<<step<<','<<end<<','<<volume<<','<<expected<<','<<state.parametricGeometry().report().maximumInventoryFractionError<<','<<balance<<'\n';save(step);
            std::cout<<std::setprecision(17)<<"PARAMETRIC_SUPPLY_STEP native="<<native<<" step="<<step<<" time="<<end<<" volume="<<volume<<" expected="<<expected<<" fraction="<<maxFraction<<" balance="<<balance<<std::endl;
        }
        std::cout<<"PARAMETRIC_SUPPLY native="<<native<<" steps="<<steps<<" time="<<state.progress().endTime<<" volume_error="<<maxVolume<<" fraction="<<maxFraction<<" free_area="<<surface.length;
        if(!native)std::cout<<" position="<<maxPosition<<" curvature="<<maxCurvature;std::cout<<std::endl;
        require(maxVolume<2e-12&&maxFraction<3.1e-12&&maxPosition<2e-10&&maxCurvature<2e-7,"supplied Hermite geometry/inventory or curvature failed");
        if(!native){
            const double begin=state.progress().endTime;
            auto back=InterfaceParametricMotion::prepare(mesh,domain,state,{steps+1,begin,2*begin,begin},-begin,analytic,analyticFlux,options);
            const auto reference=InterfacePrescribedData::parametricSupplyCap(radius,height,8);double shape=0;
            for(int k=0;k<=128;k++)shape=std::max(shape,P::norm(back.state.parametricGeometry().geometry().sample(k/128.)-reference.sample(k/128.)));
            const auto& inventory=back.state.parametricGeometry().inventory();const double volume=std::accumulate(inventory.begin(),inventory.end(),0.);
            std::cout<<"PARAMETRIC_SUPPLY_RETURN volume_error="<<std::abs(volume-initialVolume)<<" position="<<shape<<" fraction="<<back.state.parametricGeometry().report().maximumInventoryFractionError<<std::endl;
            require(std::abs(volume-initialVolume)<2e-12&&shape<2e-10,"reversed inlet phase flux did not recover the supplied cap");
            const auto oldGeometry=state.parametricGeometry().geometryHandle();bool rejected=false;
            const InterfaceParametricTransport::Flow movingRim=[&](Vector2D p,double t){auto q=analytic(p,t);q.point.y+=.01*t;q.displacement.y+=.01*t;q.timeDerivative.y+=.01;return q;};
            try{auto bad=InterfaceParametricMotion::prepare(mesh,domain,state,{steps+1,begin,begin+dt,dt},dt,movingRim,analyticFlux,options);}
            catch(const std::runtime_error& e){rejected=std::string(e.what()).find("fixed supply endpoint")!=std::string::npos;}
            require(rejected&&state.parametricGeometry().geometryHandle()==oldGeometry&&state.progress().endTime==begin,"moving supplied endpoint was not rejected without mutation");
        }
    }
}
void hydroMomentAmrSplineContract(){
    bool energyWithinBudget=true;
    const auto analytic=[](Vector2D p,bool rz){
        if(!rz)return exact(p);
        const double r=p.x,z=p.y,F=4*std::pow(1-r*r,2),Fr=-16*r*(1-r*r),Frr=-16+48*r*r;
        const double g=z*z*(1-z)*(1-z),gz=2*z-6*z*z+4*z*z*z,gzz=2-12*z+12*z*z;
        return Spline::FirstJet{{r*F*gz,-(2*F+r*Fr)*g},{(F+r*Fr)*gz,-(3*Fr+r*Frr)*g},{r*F*gzz,-(2*F+r*Fr)*gz}};
    };
    for(bool rz:{false,true})for(int n:{2,4}){
        QuadtreeMesh mesh(input(n));DictionaryMesh::RefineBox box;box.min={0,.5};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto space=Space::build(mesh,domain);
        require(!domain.uniformLeaves(),"native C1 control must cross coarse/fine faces");
        const auto u=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=analytic(p,rz).value;return Eigen::Vector2d(v.x,v.y);});
        const auto spline=Spline::build(mesh,domain,space,u,0);const auto again=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=spline.firstJet(p).value;return Eigen::Vector2d(v.x,v.y);});
        const double moment=(again-u).norm()/u.norm();double value=0,gradient=0,join=0,trace=0,div=0,work=0;int coarseFine=0;
        std::vector<Spline::PointForce> forces;
        for(int j=0;j<17;j++)for(int i=0;i<17;i++){
            const Vector2D p{(i+.31)/17,(j+.43)/17};const auto v=spline.firstJet(p),a=analytic(p,rz);
            value=std::max(value,std::hypot(v.value.x-a.value.x,v.value.y-a.value.y));
            gradient=std::max({gradient,std::hypot(v.dx.x-a.dx.x,v.dx.y-a.dx.y),std::hypot(v.dy.x-a.dy.x,v.dy.y-a.dy.y)});
            div=std::max(div,std::abs(v.dx.x+v.dy.y+(rz?v.value.x/p.x:0)));
            forces.push_back({p,{std::sin(i+.2*j),std::cos(j-.3*i)},{.1*std::cos(i),-.2*std::sin(j)},{.15*std::sin(i),.12*std::cos(j)}});
        }
        for(Vector2D p:std::vector<Vector2D>{{0,.43},{1,.43},{.31,0},{.31,1},{0,0},{0,1}})forces.push_back({p,{.3,-.2},{.1,.7},{-.2,.4}});
        for(const auto& f:forces){const auto v=spline.firstJet(f.point);work+=f.measure.x*v.value.x+f.measure.y*v.value.y+f.dx.x*v.dx.x+f.dx.y*v.dx.y+f.dy.x*v.dy.x+f.dy.y*v.dy.y;}
        const double adjoint=std::abs(u.dot(spline.transposePointForces(forces))-work)/std::max(1.,std::abs(work));
        const auto flux=spline.normalFlux();const auto& g=Quantum::Geometry::Parametric::gaussEight();
        for(const auto face:domain.faces(mesh)){
            const auto chart=flux.faceChart(face.nativeRow);const double length=face.horizontal?chart.end.x-chart.begin.x:chart.end.y-chart.begin.y;
            if(!face.outwardSign){
                const auto ca=mesh.cells(mesh.cellIndexToLabel[face.negativeCell]).size(),cb=mesh.cells(mesh.cellIndexToLabel[face.positiveCell]).size();
                if(ca.x!=cb.x||ca.y!=cb.y)coarseFine++;
                for(double t:{.13,.5,.87}){const auto p=chart.begin+(chart.end-chart.begin)*t;auto left=p,right=p;
                    if(face.horizontal){left.y=std::nextafter(p.y,-1.);right.y=std::nextafter(p.y,2.);}else{left.x=std::nextafter(p.x,-1.);right.x=std::nextafter(p.x,2.);}
                    const auto a=spline.firstJet(left),b=spline.firstJet(right);
                    join=std::max({join,std::hypot(a.value.x-b.value.x,a.value.y-b.value.y),std::hypot(a.dx.x-b.dx.x,a.dx.y-b.dx.y),std::hypot(a.dy.x-b.dy.x,a.dy.y-b.dy.y)});
                }
            }
            for(int half=0;half<2;half++)for(int k=0;k<8;k++){
                const double t=.5+.5*g.x[k];const auto a=half?chart.midpoint:chart.begin,b=half?chart.end:chart.midpoint,p=a+(b-a)*t;
                const auto v=spline.firstJet(p).value;const double direct=.5*length*(rz?2*std::acos(-1.)*p.x:1)*(face.horizontal?v.y:v.x);
                trace=std::max(trace,std::abs(double(flux.normalFluxDensity(face.nativeRow,half,false,t))-direct));
            }
        }
        const auto axis=spline.firstJet({0,.43});if(rz)require(axis.value.x==0&&axis.dx.y==0&&axis.dy.x==0&&std::abs(axis.value.y)>.01,"native C1 lost regular nonzero axis velocity");
        std::cout<<std::setprecision(17)<<"MOMENT_AMR_C1 metric="<<(rz?"RZ":"XY")<<" base="<<n<<" cells="<<mesh.cellIndexNum()<<" coarse_fine="<<coarseFine
            <<" moments="<<moment<<" value="<<value<<" gradient="<<gradient<<" join="<<join<<" trace="<<trace<<" divergence="<<div<<" virtual_work="<<adjoint<<std::endl;
        require(coarseFine>0&&moment<2e-10&&join<2e-8&&trace<2e-10&&div<2e-12&&adjoint<2e-10,"native C1 moment/trace/continuity/work contract failed");
        // Entire finite-dimensional solenoidal space, not just the smooth
        // velocity above: compare continuous kinetic energy with native Hodge.
        const int nv=int(space.mass().rows()),nd=int(space.divergence().rows());std::vector<int> walls;
        for(const auto f:domain.faces(mesh))if(f.outwardSign&&space.faces()[f.nativeRow].first>=0){walls.push_back(space.faces()[f.nativeRow].first);walls.push_back(walls.back()+1);}
        Eigen::MatrixXd constraints=Eigen::MatrixXd::Zero(nd+walls.size(),nv);constraints.topRows(nd)=space.divergence();
        for(int k=0;k<nd;k++)constraints.row(k)/=constraints.row(k).norm();for(int k=0;k<int(walls.size());k++)constraints(nd+k,walls[k])=1;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(constraints,Eigen::ComputeFullV);int rank=0;
        for(int k=0;k<svd.singularValues().size();k++)if(svd.singularValues()[k]>1e-12*svd.singularValues()[0])rank++;
        Eigen::MatrixXd basis=svd.matrixV().rightCols(nv-rank);for(int row:walls)basis.row(row).setZero();
        require((constraints*basis).norm()<2e-11,"native C1 energy basis violates solenoidality");
        double rx[12],rw[12];kineticQuadrature(rx,rw);Eigen::MatrixXd values(2*mesh.cellIndexNum()*12*12,basis.cols());
        for(int column=0;column<basis.cols();column++){
            const auto field=spline.withVelocity(mesh,domain,space,basis.col(column),0);int row=0;
            for(int index:domain.cellIndexByRow()){
                const auto& c=mesh.cells(mesh.cellIndexToLabel[index]);const auto center=c.centerCoord(),h=c.size();
                for(int a=0;a<12;a++)for(int b=0;b<12;b++){
                    const Vector2D p{center.x+.5*h.x*rx[a],center.y+.5*h.y*rx[b]};const auto v=field.firstJet(p).value;
                    const double root=std::sqrt(.25*h.x*h.y*rw[a]*rw[b]*(rz?2*std::acos(-1.)*p.x:1));values(row++,column)=root*v.x;values(row++,column)=root*v.y;
                }
            }
        }
        const Eigen::MatrixXd actual=values.transpose()*values,owned=basis.transpose()*space.mass()*basis;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen(actual,owned,Eigen::EigenvaluesOnly);
        require(eigen.info()==Eigen::Success&&eigen.eigenvalues().allFinite()&&eigen.eigenvalues()[0]>0,"native C1 kinetic norm lost injectivity");
        std::cout<<"MOMENT_AMR_KINETIC metric="<<(rz?"RZ":"XY")<<" base="<<n<<" cells="<<mesh.cellIndexNum()<<" solenoidal_dofs="<<basis.cols()
            <<" minimum="<<eigen.eigenvalues()[0]<<" maximum="<<eigen.eigenvalues().tail(1)[0]<<std::endl;
        energyWithinBudget=energyWithinBudget&&eigen.eigenvalues().tail(1)[0]<32;
    }
    require(energyWithinBudget,"native C1 lifting energy exceeds the finite-mesh control bound");
}
void hydroMomentSplineInterfaceContract(){trajectory();}
// One fixed-material Stokes step is a component check of the physical mixed
// forms, NOT another production flow route. No advection, moving geometry or
// capillary load is asserted by this test. The actual inlet and outlet stay.
static void physicalMidpointControl(const Flow& flow,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Spline::NativeVelocity& old){
    using Complement=Quantum::Hydrodynamics::HydroCompatibleMomentComplement;
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    using Matrix=Space::Matrix;
    const auto& space=flow.space();const auto& boundary=flow.boundary();
    const auto complement=Complement::build(mesh,domain,space,boundary.fixed());
    const auto forms=old.mixedMaterialForms(mesh,domain,space,flow.material(),complement);
    const int nq=forms.solenoidalCoordinates,ne=forms.complementCoordinates,nv=nq+ne,np=int(space.divergence().rows());
    const int ng=boundary.hasPressureOutlet()?0:1;
    require(ne+ng==np,"fixed-material pressure control has inconsistent gauge");
    // Keep every original wall/axis and native inlet condition for checking;
    // only exact row dependencies are removed from the multiplier system.
    const auto& h=old.kinematicConstraints();const Eigen::MatrixXd j=old.momentMap();
    int nb=int(h.rows());for(int row=0;row<j.rows();row++)nb+=boundary.fixed()[row];
    Eigen::MatrixXd all=Eigen::MatrixXd::Zero(nb,nq);Eigen::VectorXd allValues=Eigen::VectorXd::Zero(nb);
    all.topRows(h.rows())=h;int row=int(h.rows());
    for(int k=0;k<j.rows();k++)if(boundary.fixed()[k]){all.row(row)=j.row(k);allValues[row++]=boundary.values()[k];}
    for(int k=0;k<nb;k++){const double norm=all.row(k).norm();if(norm){all.row(k)/=norm;allValues[k]/=norm;}}
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> independent(all.transpose());independent.setThreshold(1e-10);
    const int rank=int(independent.rank());Eigen::MatrixXd b(rank,nq);Eigen::VectorXd g(rank);
    for(int k=0;k<rank;k++){const int source=independent.colsPermutation().indices()[k];b.row(k)=all.row(source);g[k]=allValues[source];}
    const double initialBoundary=(all*old.coordinates()-allValues).norm();require(initialBoundary<2e-10,"Stokes control changed its actual inlet initial data");
    constexpr double dt=.002;Matrix a=2/dt*forms.mass+forms.viscosity;
    const int total=nv+np+rank+ng;std::vector<Eigen::Triplet<double>> entries;
    for(int col=0;col<a.outerSize();col++)for(Matrix::InnerIterator it(a,col);it;++it)entries.emplace_back(it.row(),it.col(),it.value());
    const auto& c=complement.pressureLoad();
    for(int col=0;col<c.outerSize();col++)for(Matrix::InnerIterator it(c,col);it;++it){entries.emplace_back(nq+it.row(),nv+col,it.value());entries.emplace_back(nv+col,nq+it.row(),it.value());}
    for(int k=0;k<rank;k++)for(int col=0;col<nq;col++)if(b(k,col)){entries.emplace_back(nv+np+k,col,b(k,col));entries.emplace_back(col,nv+np+k,b(k,col));}
    if(ng)for(int cell=0;cell<mesh.cellIndexNum();cell++){
        const double volume=domain.cellVolumes()[domain.cellRowByIndex()[cell]];
        entries.emplace_back(nv+3*cell,total-1,volume);entries.emplace_back(total-1,nv+3*cell,volume);
        if(domain.axisymmetric()){const auto& leaf=mesh.cells(mesh.cellIndexToLabel[cell]);const double linear=volume*leaf.size().x/(12*leaf.centerCoord().x);
            entries.emplace_back(nv+3*cell+1,total-1,linear);entries.emplace_back(total-1,nv+3*cell+1,linear);}
    }
    Matrix system(total,total);system.setFromTriplets(entries.begin(),entries.end());
    Eigen::VectorXd initial=Eigen::VectorXd::Zero(nv),force=initial;initial.head(nq)=old.coordinates();
    const Eigen::VectorXd nativeForce=boundary.traction().head(space.mass().rows());
    force.head(nq)=old.momentMap().transpose()*nativeForce;force.tail(ne)=complement.nativeMoments().transpose()*nativeForce;
    Eigen::VectorXd rhs=Eigen::VectorXd::Zero(total);rhs.head(nv)=2/dt*forms.mass*initial+force;rhs.segment(nv+np,rank)=g;
    Eigen::SuperLU<Matrix> solver;solver.compute(system);require(solver.info()==Eigen::Success,"physical mixed Stokes factorization failed");
    const Eigen::VectorXd solution=solver.solve(rhs);require(solver.info()==Eigen::Success&&solution.allFinite(),"physical mixed Stokes solve failed");
    Quantum::Hydrodynamics::HydroCompatibleMomentBlockSolve block;
    block.compute(system,nq,ne,np,rank,ng,old.cellInteriorCoordinates());const Eigen::VectorXd blocked=block.solve(rhs);
    require(block.statistics().retainedVelocityUnknowns==nq+rank-9*int(old.cellInteriorCoordinates().size()),
        "complete block solver did not eliminate exactly the owned cell interiors");
    const double blockEquation=(system*blocked-rhs).norm()/rhs.norm(),blockDifference=(blocked-solution).norm()/solution.norm();
    std::cout<<std::setprecision(17)<<"MOMENT_PHYSICAL_BLOCK_SOLVE metric="<<(domain.axisymmetric()?"RZ":"XY")
        <<" closed="<<bool(ng)<<" equation="<<blockEquation<<" full_solution_difference="<<blockDifference<<" eta="<<blocked.segment(nq,ne).norm()<<std::endl;
    require(blockEquation<2e-10&&blockDifference<2e-9&&blocked.segment(nq,ne).isZero(0),
        "exact block elimination changed the complete mixed Stokes solution");
    // Change the physical mass/time coefficient and its right hand side,
    // then compare retained-factor solves with fresh FULL saddle solves.
    // This catches both stale momentum matrices and stale pressure coupling.
    for(double ratio:{1.01,1.03}){
        auto changedEntries=entries;const double delta=2/dt*(ratio-1);
        for(int col=0;col<forms.mass.outerSize();col++)for(Matrix::InnerIterator it(forms.mass,col);it;++it)
            changedEntries.emplace_back(it.row(),it.col(),delta*it.value());
        Matrix changed(total,total);changed.setFromTriplets(changedEntries.begin(),changedEntries.end());
        Eigen::VectorXd changedRhs=rhs;changedRhs.head(nv)+=delta*(forms.mass*initial);
        solver.compute(changed);require(solver.info()==Eigen::Success,"updated full saddle factorization failed");
        const Eigen::VectorXd reference=solver.solve(changedRhs);
        require(solver.info()==Eigen::Success&&reference.allFinite(),"updated full saddle reference failed");
        block.compute(changed,nq,ne,np,rank,ng,old.cellInteriorCoordinates());
        const Eigen::VectorXd reused=block.solve(changedRhs,1e-12*changedRhs.norm());
        const double equation=(changed*reused-changedRhs).norm()/changedRhs.norm();
        const double difference=(reused-reference).norm()/reference.norm();const auto& stats=block.statistics();
        std::cout<<"MOMENT_RETAINED_BLOCK_SOLVE metric="<<(domain.axisymmetric()?"RZ":"XY")<<" closed="<<bool(ng)
            <<" mass_ratio="<<ratio<<" equation="<<equation<<" full_solution_difference="<<difference
            <<" velocity_factors="<<stats.velocityFactorizations<<" pressure_factors="<<stats.pressureFactorizations
            <<" krylov="<<stats.krylovIterations<<std::endl;
        require(equation<2e-10&&difference<2e-9&&reused.segment(nq,ne).isZero(0),
            "retained factor changed the current complete mixed solution");
        require(stats.velocityFactorizations==1&&stats.pressureFactorizations==1&&stats.krylovIterations>0,
            "nearby physical systems did not reuse their block factors");
    }
    // An explicit direct solve refreshes only the changed velocity block.
    block.compute(system,nq,ne,np,rank,ng,old.cellInteriorCoordinates());
    const Eigen::VectorXd refreshed=block.solve(rhs,0);const auto& refreshedStats=block.statistics();
    require(refreshed.allFinite()&&(system*refreshed-rhs).norm()/rhs.norm()<2e-10
        &&(refreshed-solution).norm()/solution.norm()<2e-9
        &&refreshedStats.velocityFactorizations==2&&refreshedStats.pressureFactorizations==1,
        "direct refresh did not retain the unchanged physical pressure block");
    const Eigen::VectorXd mid=solution.head(nv),next=2*mid-initial;
    const auto accepted=old.withCoordinates(next.head(nq),old.time()+dt);
    const double equation=(system*solution-rhs).norm()/rhs.norm(),constraints=(all*accepted.coordinates()-allValues).norm();
    const double eta=mid.tail(ne).norm(),energy0=.5*initial.dot(forms.mass*initial),energy1=.5*next.dot(forms.mass*next);
    const double dissipation=dt*mid.dot(forms.viscosity*mid),tractionWork=dt*mid.dot(force),reactionWork=-dt*solution.segment(nv+np,rank).dot(g);
    const double energy=std::abs(energy1-energy0+dissipation-tractionWork-reactionWork)/(energy0+energy1+dissipation+std::abs(tractionWork)+std::abs(reactionWork));
    const double change=(accepted.coordinates()-old.coordinates()).norm()/old.coordinates().norm();
    const double divergence=(space.divergence()*accepted.nativeMoments()).norm();
    // The inviscid limit is the SAME quadrature path, with no floor for mu.
    auto inviscidProperties=flow.material().properties();inviscidProperties.liquidViscosity=0;inviscidProperties.gasViscosity=0;
    const auto inviscid=Material::buildFrozen(mesh,domain,space,flow.material().boundaryHandle(),old.time(),inviscidProperties);
    const auto zero=old.mixedMaterialForms(mesh,domain,space,inviscid,complement);
    const double massChange=(zero.mass-forms.mass).norm()/forms.mass.norm();
    std::cout<<std::setprecision(17)<<"MOMENT_PHYSICAL_MIDPOINT metric="<<(domain.axisymmetric()?"RZ":"XY")<<" cells="<<mesh.cellIndexNum()
        <<" velocity_coordinates="<<nv<<" pressure_rows="<<np<<" boundary_rank="<<rank<<" equation="<<equation<<" boundary="<<constraints<<" eta="<<eta
        <<" divergence="<<divergence<<" energy_before="<<energy0<<" energy_after="<<energy1<<" dissipated="<<dissipation<<" outlet_work="<<tractionWork
        <<" inlet_work="<<reactionWork<<" energy_balance="<<energy<<" velocity_change="<<change<<" pressure_norm="<<solution.segment(nv,np).norm()
        <<" zero_mu_matrix="<<zero.viscosity.norm()<<" zero_mu_mass_change="<<massChange<<" formal_ns=0"<<std::endl;
    require(equation<2e-10&&constraints<2e-9&&eta<2e-10&&divergence<2e-10&&energy<2e-9&&change>1e-5&&zero.viscosity.nonZeros()==0&&massChange<2e-10,
        "physical mixed Stokes step lost conservation, energy, boundary data or exact zero viscosity");
}
static void physicalComplementControl(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Space& space,
    const Quantum::Hydrodynamics::HydroCompatibleMomentBoundary& boundary,const Quantum::Hydrodynamics::HydroCompatibleMomentMaterial& material,
    const Spline::NativeVelocity& state,const Spline::NativeMaterialForms& original){
    using Complement=Quantum::Hydrodynamics::HydroCompatibleMomentComplement;
    for(bool closed:{false,true}){
        auto fixed=boundary.fixed();if(closed)for(const auto f:domain.faces(mesh))if(f.outwardSign&&f.area>0)for(int k=0;k<2;k++)fixed[space.faces()[f.nativeRow].first+k]=true;
        const auto complement=Complement::build(mesh,domain,space,fixed);const int ne=int(complement.divergence().cols());
        Eigen::VectorXd eta(ne);for(int k=0;k<ne;k++)eta[k]=std::sin(.31*k+.19);
        const Eigen::VectorXd coefficients=complement.basisMap()*eta;
        const auto local=[&](int cell,Vector2D p)->Eigen::VectorXd{const auto& rows=complement.cells()[cell].basis;Eigen::VectorXd v(rows.size());
            for(int k=0;k<int(rows.size());k++)v[k]=coefficients[rows[k]];return complement.cellJets(cell,p)*v;};
        double valueJump=0,jetJump=0,axis=0;Space::Partition partition;
        for(const auto& cell:complement.cells())for(int d=0;d<2;d++)for(double x:cell.partition[d])partition[d].push_back(x);
        for(auto& axisCuts:partition){std::sort(axisCuts.begin(),axisCuts.end());axisCuts.erase(std::unique(axisCuts.begin(),axisCuts.end()),axisCuts.end());}
        for(const auto f:domain.faces(mesh)){
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            for(double t:{0.,.13,.5,.87,1.}){
                const auto p=a+(b-a)*t;if(f.negativeCell>=0&&f.positiveCell>=0){const auto left=local(f.negativeCell,p),right=local(f.positiveCell,p);
                    valueJump=std::max(valueJump,(left.head(2)-right.head(2)).norm());jetJump=std::max(jetJump,(left.tail(4)-right.tail(4)).norm());}
                if(domain.axisymmetric()&&p.x==0){const auto v=complement.firstJet(p,eta);axis=std::max({axis,std::abs(v.value.x),std::abs(v.dx.y)});}
            }
        }
        const auto actual=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=complement.firstJet(p,eta).value;return Eigen::Vector2d(v.x,v.y);},partition);
        const Eigen::VectorXd expected=complement.nativeMoments()*eta;
        const double moments=(actual-expected).norm()/expected.norm();
        double massBlock=0,viscosityBlock=0,energyDifference=0,fullEnergy=0,fullDissipation=0,symmetry=0;
        if(!closed){
            const auto forms=state.mixedMaterialForms(mesh,domain,space,material,complement);const int nq=forms.solenoidalCoordinates;
            massBlock=(Space::Matrix(forms.mass.topLeftCorner(nq,nq))-original.mass).norm()/original.mass.norm();
            viscosityBlock=(Space::Matrix(forms.viscosity.topLeftCorner(nq,nq))-original.viscosity).norm()/original.viscosity.norm();
            Eigen::VectorXd combined(nq+ne);combined.head(nq)=state.coordinates();combined.tail(ne)=eta;
            fullEnergy=.5*combined.dot(forms.mass*combined);fullDissipation=combined.dot(forms.viscosity*combined);
            double localEnergy=0;
            // Independently assemble the complementary diagonal from its
            // local coefficients; qq and cross work use the global forms.
            for(const auto& cell:forms.cells){const int count=int(complement.cells()[cell.index].basis.size());
                const auto& rows=complement.cells()[cell.index].basis;
                Eigen::VectorXd e(count);for(int k=0;k<count;k++)e[k]=coefficients[rows[k]];
                localEnergy+=.5*e.dot(cell.mass.bottomRightCorner(count,count)*e);
            }
            const Eigen::VectorXd q=state.coordinates();
            localEnergy+=.5*q.dot(forms.mass.topLeftCorner(nq,nq)*q)+q.dot(forms.mass.topRightCorner(nq,ne)*eta);
            energyDifference=std::abs(localEnergy-fullEnergy)/fullEnergy;
            symmetry=std::max((forms.mass-Space::Matrix(forms.mass.transpose())).norm()/forms.mass.norm(),
                (forms.viscosity-Space::Matrix(forms.viscosity.transpose())).norm()/forms.viscosity.norm());
        }
        std::cout<<std::setprecision(17)<<"MOMENT_PHYSICAL_COMPLEMENT metric="<<(domain.axisymmetric()?"RZ":"XY")<<" closed="<<closed<<" cells="<<mesh.cellIndexNum()
            <<" pressure_rows="<<complement.divergence().rows()<<" complement_coordinates="<<ne<<" commuting="<<complement.commutingResidual()
            <<" physical_moments="<<moments<<" value_jump="<<valueJump<<" jet_jump="<<jetJump<<" axis="<<axis<<" mass_block="<<massBlock
            <<" viscosity_block="<<viscosityBlock<<" mixed_energy="<<fullEnergy<<" mixed_dissipation="<<fullDissipation<<" assembly_energy="<<energyDifference<<" symmetry="<<symmetry<<std::endl;
        require(complement.commutingResidual()<2e-12&&moments<2e-10&&valueJump<2e-10&&jetJump<2e-8&&axis<2e-12&&massBlock<2e-10&&viscosityBlock<2e-10
            &&energyDifference<2e-10&&symmetry<2e-12&&(closed||(fullEnergy>0&&fullDissipation>0)),"physical pressure complement lost continuity, moments or full material cross forms");
    }
}
static void nativeMassCurrentControl(const Flow& flow,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const Spline::NativeVelocity& velocity,const InterfaceState& next,const InterfaceStepTime& time,bool completeMomentum=false,
    const Quantum::Hydrodynamics::HydroCompatibleMomentBoundary* suppliedBoundary=nullptr,int direction=1){
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;
    using Current=Quantum::Hydrodynamics::HydroCompatibleMomentMassCurrent;
    const auto& previous=flow.interfaceState();const auto& space=flow.space();
    const auto& boundary=suppliedBoundary?*suppliedBoundary:flow.boundary();
    const auto wind=velocity.withCoordinates(velocity.coordinates(),time.begin+.5*time.dt);
    Flow::Options options;options.interfaceProjection.materialSegments=16;
    const InterfaceStepTime half{time.index,time.begin,wind.time(),.5*time.dt};
    const auto middle=InterfaceParametricMotion::prepare(mesh,domain,previous,half,half.dt,
        HydroInterfaceTransport::frozenParametricFlow(velocity),velocity.normalFlux(),options.interfaceProjection);
    Material::Properties properties;properties.liquidDensity=1000;properties.gasDensity=1;
    const auto midpoint=Material::buildFrozen(mesh,domain,space,middle.state.parametricGeometry().boundaryHandle(),half.end,properties);
    const auto complement=Quantum::Hydrodynamics::HydroCompatibleMomentComplement::build(mesh,domain,space,boundary.fixed());
    const auto current=Current::build(mesh,domain,space,wind,previous,next,time,midpoint,2e-12,&complement);
    std::optional<Spline::NativeMomentumTransport> momentum;Spline::NativeVolumeTransport assembled;
    if(completeMomentum){momentum=wind.nativeMomentumTransport(mesh,domain,space,current,boundary,&complement);assembled=*momentum;}
    else assembled=wind.nativeVolumeTransport(mesh,domain,space,current,&complement);
    Eigen::VectorXd coordinates=Eigen::VectorXd::Zero(assembled.volume.rows());coordinates.head(wind.coordinates().size())=wind.coordinates();
    const Eigen::VectorXd assembledLoad=assembled.volume*coordinates;
    Eigen::VectorXd eta(assembled.complementCoordinates);for(int k=0;k<eta.size();k++)eta[k]=std::sin(.371*(k+1));
    const Eigen::VectorXd etaCoefficients=complement.basisMap()*eta;
    const auto mass=Transport::phaseMassFlux(mesh,domain,space,previous,next,time,1000,1);
    std::vector<Eigen::Vector3d> face(mesh.cellIndexNum(),Eigen::Vector3d::Zero());
    for(const auto f:domain.faces(mesh)){const auto products=velocity.faceProducts(f.nativeRow);Eigen::Vector3d value=Eigen::Vector3d::Zero();
        for(int k=0;k<=InterfaceParametricTransport::FaceVolumes::productDegree;k++){
            value[0]+=mass.products[f.nativeRow][k]*products.x[k];value[1]+=mass.products[f.nativeRow][k]*products.y[k];value[2]+=mass.products[f.nativeRow][k]*products.speedSquared[k];}
        if(f.negativeCell>=0)face[f.negativeCell]+=value;if(f.positiveCell>=0)face[f.positiveCell]-=value;
    }
    double constant=0,gradient=0,reference=0,work=0,energy=0,orthogonal=0,correctionSquared=0,referenceSquared=0,quad=0;
    Eigen::Vector3d firstGlobal=Eigen::Vector3d::Zero(),lastGlobal=firstGlobal;std::size_t nodes=0;int modes=0;
    long double directVolumeWork=0,directComplementWork=0;double pureCorrectionSquared=0,pureReferenceSquared=0;
    for(const auto& cell:current.cells()){
        const auto integral=[&](const auto& points){Eigen::Vector3d result=Eigen::Vector3d::Zero();
            for(const auto& p:points){const auto u=wind.firstJet(p.point).value;result+=p.massWeight*Eigen::Vector3d(u.x,u.y,u.x*u.x+u.y*u.y);}return result;};
        const auto first=integral(cell.before),last=integral(cell.after);firstGlobal+=first;lastGlobal+=last;
        Eigen::Vector3d bulk=Eigen::Vector3d::Zero();
        const auto& complementRows=complement.cells()[cell.index].basis;Eigen::VectorXd e(complementRows.size());
        for(int k=0;k<e.size();k++)e[k]=etaCoefficients[complementRows[k]];
        for(const auto& p:cell.current){const auto u=wind.firstJet(p.point);bulk[0]+=p.weight.x*u.dx.x+p.weight.y*u.dy.x;bulk[1]+=p.weight.x*u.dx.y+p.weight.y*u.dy.y;
            bulk[2]+=2*p.weight.x*(u.value.x*u.dx.x+u.value.y*u.dx.y)+2*p.weight.y*(u.value.x*u.dy.x+u.value.y*u.dy.y);
            const Eigen::VectorXd v=complement.cellJets(cell.index,p.point)*e;
            directComplementWork-=p.weight.x*(v[2]*u.value.x+v[3]*u.value.y)+p.weight.y*(v[4]*u.value.x+v[5]*u.value.y);
        }
        directVolumeWork-=.5*bulk[2];
        // Independent zero-normal, weighted-divergence-free polynomial
        // bubble. A gradient correction must leave its mass-current work.
        long double curlWork=0;
        for(const auto& p:cell.correction){double x=(p.point.x-cell.lower.x)/cell.size.x,y=(p.point.y-cell.lower.y)/cell.size.y;
            const double bx=x*x*(1-x)*(1-x),by=y*y*(1-y)*(1-y),dx=2*x*(1-x)*(1-2*x)/cell.size.x,dy=2*y*(1-y)*(1-2*y)/cell.size.y;
            const double metric=domain.axisymmetric()?p.point.x:1.;
            curlWork+=(static_cast<long double>(p.weight.x)*bx*dy-static_cast<long double>(p.weight.y)*dx*by)/metric;}
        const double scale=1000*domain.cellVolumes()[domain.cellRowByIndex()[cell.index]];
        const Eigen::Vector3d residual=last-first+face[cell.index]-bulk;
        work=std::max(work,residual.head<2>().norm()/scale);energy=std::max(energy,std::abs(residual[2])/scale);
        orthogonal=std::max(orthogonal,double(std::abs(curlWork))/scale);
        constant=std::max(constant,cell.constantDefect);gradient=std::max(gradient,cell.gradientDefect);reference=std::max(reference,cell.referenceDefect);
        correctionSquared+=cell.correctionNorm*cell.correctionNorm;referenceSquared+=cell.referenceNorm*cell.referenceNorm;
        const double density=cell.before.front().density;
        const auto pure=[&](const auto& points){return std::all_of(points.begin(),points.end(),[&](const auto& p){return p.density==density;});};
        if(pure(cell.before)&&pure(cell.midpoint)&&pure(cell.after)){pureCorrectionSquared+=cell.correctionNorm*cell.correctionNorm;pureReferenceSquared+=cell.referenceNorm*cell.referenceNorm;}
        quad=std::max(quad,cell.quadratureEstimate);nodes+=cell.current.size();modes=std::max(modes,int(cell.potential.size()));
    }
    const auto before=Material::buildFrozen(mesh,domain,space,previous.parametricGeometry().boundaryHandle(),time.begin,properties);
    const auto after=Material::buildFrozen(mesh,domain,space,next.parametricGeometry().boundaryHandle(),time.end,properties);
    const auto firstState=velocity.withCoordinates(velocity.coordinates(),time.begin),lastState=velocity.withCoordinates(velocity.coordinates(),time.end);
    const auto firstForms=completeMomentum?firstState.mixedMaterialForms(mesh,domain,space,before,complement):firstState.nativeMaterialForms(mesh,domain,space,before);
    const auto lastForms=completeMomentum?lastState.mixedMaterialForms(mesh,domain,space,after,complement):lastState.nativeMaterialForms(mesh,domain,space,after);
    const Eigen::Vector3d firstReference(firstForms.physicalMomentum.x,firstForms.physicalMomentum.y,2*firstForms.kinetic);
    const Eigen::Vector3d lastReference(lastForms.physicalMomentum.x,lastForms.physicalMomentum.y,2*lastForms.kinetic);
    const double material=std::max((firstGlobal-firstReference).norm()/firstReference.norm(),(lastGlobal-lastReference).norm()/lastReference.norm());
    const double volumeWork=double(std::abs(coordinates.dot(assembledLoad)-directVolumeWork)/std::max(std::abs(directVolumeWork),1e-30L));
    const double complementWork=double(std::abs(eta.dot(assembledLoad.tail(eta.size()))-directComplementWork)/std::max(std::abs(directComplementWork),1e-30L));
    std::cout<<std::setprecision(17)<<"MOMENT_MASS_CURRENT metric="<<(domain.axisymmetric()?"RZ":"XY")<<" direction="<<direction<<" cells="<<mesh.cellIndexNum()<<" modes="<<modes
        <<" constant="<<constant<<" gradient="<<gradient<<" reference="<<reference<<" momentum_work="<<work<<" energy_work="<<energy
        <<" material="<<material<<" curl_work="<<orthogonal<<" relative_correction="<<std::sqrt(correctionSquared/referenceSquared)
        <<" quadrature="<<quad<<" nodes="<<nodes<<" volume_work="<<volumeWork<<" complement_work="<<complementWork
        <<" complement_load="<<assembledLoad.tail(eta.size()).norm()<<" pure_correction="<<std::sqrt(pureCorrectionSquared/pureReferenceSquared)<<" formal_ns=0"<<std::endl;
    require(constant<2e-11&&gradient<2e-12&&work<2e-10&&energy<2e-10&&material<2e-10&&orthogonal<2e-12&&volumeWork<2e-9&&complementWork<2e-9&&assembledLoad.tail(eta.size()).norm()>1e-10,
        "complete volume mass current lost material continuity, physical work or rotational complement");
    if(momentum){
        const auto& m=*momentum;const Space::Matrix delta=lastForms.mass-firstForms.mass;
        const Space::Matrix defect=m.central+Space::Matrix(m.central.transpose())+delta-m.boundaryMetric;
        const double fullMetric=defect.norm()/(firstForms.mass.norm()+lastForms.mass.norm());
        const int nq=int(wind.coordinates().size());
        const double physicalMetric=Space::Matrix(defect.topLeftCorner(nq,nq)).norm()
            /(Space::Matrix(firstForms.mass.topLeftCorner(nq,nq)).norm()+Space::Matrix(lastForms.mass.topLeftCorner(nq,nq)).norm());
        const double transportWork=coordinates.dot(m.central*coordinates-m.boundaryLoad);
        const double materialWork=.5*coordinates.dot(delta*coordinates);
        const double energyScale=std::abs(transportWork)+std::abs(materialWork)+std::abs(m.boundaryKineticFlux)+std::abs(m.boundaryInflowDissipation);
        const double balance=std::abs(transportWork+materialWork-m.boundaryKineticFlux-m.boundaryInflowDissipation)/energyScale;
        const double geometricWork=.5*coordinates.dot(m.geometricMetric*coordinates);
        const double accountedBalance=std::abs(transportWork+materialWork-m.boundaryKineticFlux-m.boundaryInflowDissipation-geometricWork)/energyScale;
        const double accountedMetric=(defect-m.geometricMetric).norm()/(firstForms.mass.norm()+lastForms.mass.norm());
        long double directKinetic=0,directDissipation=0;
        for(const auto f:domain.faces(mesh))if(f.outwardSign&&f.area>0){
            const auto& port=boundary.faces()[f.nativeRow];if(port.kind==Quantum::Hydrodynamics::HydroCompatibleMomentBoundary::Kind::Wall)continue;
            const auto trace=velocity.faceProducts(f.nativeRow);const double lower=(f.horizontal?f.center.x:f.center.y)-.5*f.tangentSpacing;
            const auto one=InterfaceParametricTransport::FaceVolumes::constantProduct(f.horizontal,domain.axisymmetric(),lower,lower+f.tangentSpacing);
            const Eigen::VectorXd out=f.outwardSign>0?mass.exteriorPositiveProducts[f.nativeRow]:Eigen::VectorXd(mass.exteriorPositiveProducts[f.nativeRow]-mass.products[f.nativeRow]);
            const Eigen::VectorXd in=out-f.outwardSign*mass.products[f.nativeRow];
            Vector2D exterior{0,0};if(port.incomingVelocity)exterior=*port.incomingVelocity;
            for(int k=0;k<=20;k++){const long double squared=(exterior.x*exterior.x+exterior.y*exterior.y)*one[k];
                directKinetic+=.5L*(out[k]*trace.speedSquared[k]-in[k]*squared);
                directDissipation+=.5L*in[k]*(trace.speedSquared[k]-2*exterior.x*trace.x[k]-2*exterior.y*trace.y[k]+squared);}
        }
        const double boundaryWork=double((std::abs(directKinetic-m.boundaryKineticFlux)+std::abs(directDissipation-m.boundaryInflowDissipation))/energyScale);
        std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_MOMENTUM metric="<<(domain.axisymmetric()?"RZ":"XY")
            <<" direction="<<direction<<" cells="<<mesh.cellIndexNum()<<" coordinates="<<m.central.rows()<<" trace="<<m.traceDefect<<" full_metric="<<fullMetric<<" physical_metric="<<physicalMetric
            <<" energy_balance="<<balance<<" boundary_work="<<boundaryWork<<" kinetic_flux="<<m.boundaryKineticFlux<<" inflow_dissipation="<<m.boundaryInflowDissipation
            <<" boundary_volume="<<m.boundaryVolume<<" geometric_work="<<geometricWork<<" accounted_energy="<<accountedBalance<<" accounted_metric="<<accountedMetric
            <<" raw_energy_gate="<<(balance<2e-9)<<" pressure_backflow="<<m.pressureBackflowFaces<<" bidirectional="<<m.bidirectionalFaces<<" formal_ns=0"<<std::endl;
        require(m.traceDefect<2e-10&&fullMetric<2e-10&&physicalMetric<2e-10&&accountedBalance<2e-9&&accountedMetric<2e-10&&boundaryWork<2e-9&&m.boundaryInflowDissipation>=0
            &&(direction>0||m.pressureBackflowFaces>0),
            "complete native central momentum lost moving-material or open-boundary work balance");
    }
}
static void nativePhaseFluxControl(const Flow& flow,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Spline::NativeVelocity& velocity,bool volumeCurrent=false,bool completeMomentum=false){
    using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;using Volumes=InterfaceParametricTransport::FaceVolumes;
    const auto flux=velocity.normalFlux();
    const auto motion=HydroInterfaceTransport::frozenParametricFlow(velocity);
    const double dt=.002;const auto& old=flow.interfaceState();const double begin=flow.state().time;
    const InterfaceStepTime time{old.progress().completedSteps+1,begin,begin+dt,dt};
    Flow::Options options;options.interfaceProjection.materialSegments=16;
    const auto next=InterfaceParametricMotion::prepare(mesh,domain,old,time,dt,motion,flux,options.interfaceProjection);
    const auto mass=Transport::phaseMassFlux(mesh,domain,flow.space(),old,next.state,time,1000,1);
    require(next.state.parametricFluxes().hasVelocityProducts&&mass.products.size()==std::size_t(domain.nativeFaceCount()),"complete velocity did not publish its actual shared phase functionals");
    double trace=0,energyTrace=0,totalMomentum=0,totalEnergy=0;const auto& rule=Quantum::Geometry::Parametric::gaussEight();
    for(const auto f:domain.faces(mesh)){
        if(f.area==0)continue;const auto chart=flux.faceChart(f.nativeRow);const auto product=velocity.faceProducts(f.nativeRow);const double lower=f.horizontal?chart.begin.x:chart.begin.y,upper=f.horizontal?chart.end.x:chart.end.y;
        long double directX=0,directY=0,directEnergy=0;
        // Independent physical point values and flux density, with composite
        // Gauss integration, check the algebraic Bernstein trace conversion.
        for(int panel=0;panel<4;panel++)for(int k=0;k<8;k++){
            const double t=(panel+.5+.5*rule.x[k])/4;const auto p=chart.begin+(chart.end-chart.begin)*t;const auto u=velocity.firstJet(p).value;
            double vx=0,vy=0,speed=0;for(int j=0;j<=Volumes::productDegree;j++){const double b=Volumes::productWeight(j,t,f.horizontal,domain.axisymmetric(),lower,upper);
                vx+=product.x[j]*b;vy+=product.y[j]*b;speed+=product.speedSquared[j]*b;}
            trace=std::max(trace,std::hypot(vx-u.x,vy-u.y));energyTrace=std::max(energyTrace,std::abs(speed-u.x*u.x-u.y*u.y));
            const int half=t<.5?0:1;const long double w=dt*.25*rule.w[k]*flux.normalFluxDensity(f.nativeRow,half,false,2*t-half);
            directX+=w*u.x;directY+=w*u.y;directEnergy+=w*(u.x*u.x+u.y*u.y);
        }
        const auto& retained=InterfaceParametricMotion::faceVolumes(next.state,f.nativeRow);long double x=0,y=0,energy=0;
        for(int k=0;k<=Volumes::productDegree;k++){x+=product.x[k]*static_cast<long double>(retained.totalProducts[k]);y+=product.y[k]*static_cast<long double>(retained.totalProducts[k]);energy+=product.speedSquared[k]*static_cast<long double>(retained.totalProducts[k]);}
        totalMomentum=std::max(totalMomentum,std::hypot(double(x-directX),double(y-directY)));totalEnergy=std::max(totalEnergy,double(std::abs(energy-directEnergy)));
    }
    double movement=0;for(int i=0;i<64;i++){const auto a=old.parametricGeometry().geometry().sample(i/64.),b=next.state.parametricGeometry().geometry().sample(i/64.);movement=std::max(movement,Quantum::Geometry::Parametric::norm(b-a));}
    std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_FLUX metric="<<(domain.axisymmetric()?"RZ":"XY")<<" cells="<<mesh.cellIndexNum()<<" functionals="<<Volumes::productDegree+1
        <<" trace="<<trace<<" energy_trace="<<energyTrace<<" total_momentum="<<totalMomentum<<" total_energy="<<totalEnergy
        <<" shared_constant="<<mass.productConstantDefect<<" mass_closure="<<mass.maximumRelativeClosure<<" inventory="<<next.inventoryBalance<<" movement="<<movement
        <<" fraction="<<next.state.parametricGeometry().report().maximumInventoryFractionError<<" weighted_quadrature="<<next.flux.momentQuadratureError<<" formal_ns=0"<<std::endl;
    require(trace<2e-10&&energyTrace<2e-10&&totalMomentum<2e-12&&totalEnergy<2e-12&&mass.productConstantDefect<2e-12&&mass.maximumRelativeClosure<2e-12&&std::abs(next.inventoryBalance)<2e-13&&movement>1e-7,
        "complete velocity phase motion or shared mass/work face integrals failed");
    if(volumeCurrent)nativeMassCurrentControl(flow,mesh,domain,velocity,next.state,time,completeMomentum);
    if(completeMomentum){
        using Boundary=Quantum::Hydrodynamics::HydroCompatibleMomentBoundary;
        // A separate prescribed-flow control: reverse the SAME full field
        // and the inlet normal data, retaining the actual material geometry
        // and the declared gas/zero-momentum reservoir at pressure outlets.
        // This is not a second accepted NS trajectory.
        std::vector<HydroBoundaryPatch> patches;
        for(const auto f:domain.faces(mesh))if(f.outwardSign&&f.area>0){
            const auto& b=flow.boundary().faces()[f.nativeRow];if(b.kind==Boundary::Kind::Wall)continue;
            HydroBoundaryPatch p;p.name="reversed_"+std::to_string(f.nativeRow);
            p.side=f.horizontal?(f.outwardSign>0?HydroBoundarySide::North:HydroBoundarySide::South):(f.outwardSign>0?HydroBoundarySide::East:HydroBoundarySide::West);
            p.type=b.kind==Boundary::Kind::NormalVelocity?HydroBoundaryType::NormalVelocity:HydroBoundaryType::PressureOutlet;
            p.start=(f.horizontal?f.center.x:f.center.y)-.5*f.tangentSpacing;p.end=p.start+f.tangentSpacing;
            p.value=b.kind==Boundary::Kind::NormalVelocity?-b.outwardVelocity:b.pressure;p.incomingPhase=b.incomingPhase;
            if(b.kind==Boundary::Kind::PressureOutlet)p.incomingVelocity=b.incomingVelocity;
            patches.push_back(p);
        }
        const auto reversed=Boundary::build(mesh,domain,flow.space(),flow.viscosity(),patches);
        const auto backward=velocity.withCoordinates(-velocity.coordinates(),velocity.time());
        const auto reverseStep=InterfaceParametricMotion::prepare(mesh,domain,old,time,dt,
            HydroInterfaceTransport::frozenParametricFlow(backward),backward.normalFlux(),options.interfaceProjection);
        nativeMassCurrentControl(flow,mesh,domain,backward,reverseStep.state,time,true,&reversed,-1);
    }
}
enum class NativeTransferCheck {None,Restriction,Material,State,PhysicalComplement,PhaseFlux,VolumeCurrent,FaceMomentum};
static void nativeRefinementControl(const Flow& flow,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const QuadtreeMesh& target,const NativeLeafDomain& targetDomain,const std::vector<HydroBoundaryPatch>& ports,NativeTransferCheck check){
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
    using Boundary=Quantum::Hydrodynamics::HydroCompatibleMomentBoundary;
    const bool rz=domain.axisymmetric();const auto space=Space::build(target,targetDomain);
    const auto material=Material::buildFrozen(target,targetDomain,space,flow.interfaceState().parametricGeometry().boundaryHandle(),flow.state().time,flow.material().properties());
    const auto stress=Stress::build(target,targetDomain,space,material.viscosityMoments());
    const auto boundary=Boundary::build(target,targetDomain,space,stress,ports);
    const auto original=*flow.state().physicalVelocity;
    if(check==NativeTransferCheck::PhaseFlux||check==NativeTransferCheck::VolumeCurrent||check==NativeTransferCheck::FaceMomentum)
        nativePhaseFluxControl(flow,mesh,domain,original,check!=NativeTransferCheck::PhaseFlux,check==NativeTransferCheck::FaceMomentum);
    if(check==NativeTransferCheck::FaceMomentum)return;
    const auto retained=original.refine(target,targetDomain,space,&boundary);
    // This comparison deliberately does not accept a Hydro state. Native
    // moments alone do not carry all the retained polynomial information.
    const auto rebuilt=Spline::build(target,targetDomain,space,retained.nativeMoments(),retained.time(),&boundary);
    double value=0,jet=0,lostValue=0,lostJet=0,divergence=0,trace=0,inletMean=0,inletLinear=0;
    const auto compare=[&](Vector2D p){const auto a=original.firstJet(p),b=retained.firstJet(p),c=rebuilt.firstJet(p);
        value=std::max(value,std::hypot(a.value.x-b.value.x,a.value.y-b.value.y));
        jet=std::max({jet,std::hypot(a.dx.x-b.dx.x,a.dx.y-b.dx.y),std::hypot(a.dy.x-b.dy.x,a.dy.y-b.dy.y)});
        lostValue=std::max(lostValue,std::hypot(c.value.x-b.value.x,c.value.y-b.value.y));
        lostJet=std::max({lostJet,std::hypot(c.dx.x-b.dx.x,c.dx.y-b.dx.y),std::hypot(c.dy.x-b.dy.x,c.dy.y-b.dy.y)});
        divergence=std::max(divergence,std::abs(b.dx.x+b.dy.y+(rz?(p.x?b.value.x/p.x:b.dx.x):0)));
    };
    double nodes[12],weights[12];kineticQuadrature(nodes,weights);long double oldNorm=0,newNorm=0;
    for(int index:targetDomain.cellIndexByRow()){
        const auto& c=target.cells(target.cellIndexToLabel[index]);const auto p=c.centerCoord(),h=c.size();
        for(int a=0;a<12;a++)for(int b=0;b<12;b++){
            const Vector2D q{p.x+.5*h.x*nodes[a],p.y+.5*h.y*nodes[b]};compare(q);
            const auto u=original.firstJet(q).value,v=retained.firstJet(q).value;
            const long double w=.25*h.x*h.y*weights[a]*weights[b]*(rz?2*std::acos(-1.)*q.x:1);
            oldNorm+=w*(u.x*u.x+u.y*u.y);newNorm+=w*(v.x*v.x+v.y*v.y);
        }
    }
    const auto flux=retained.normalFlux();const auto& rule=Quantum::Geometry::Parametric::gaussEight();
    for(const auto f:targetDomain.faces(target)){
        const auto chart=flux.faceChart(f.nativeRow);const double length=std::hypot(chart.end.x-chart.begin.x,chart.end.y-chart.begin.y);
        for(const auto p:{chart.begin,chart.midpoint,chart.end})compare(p);
        for(int half=0;half<2;half++)for(int k=0;k<8;k++){
            const double t=.5+.5*rule.x[k];const auto p=chart.begin+(chart.end-chart.begin)*(.5*(half+t));compare(p);
            const auto u=retained.firstJet(p).value;const double expected=.5*length*(rz?2*std::acos(-1.)*p.x:1)*(f.horizontal?u.y:u.x);
            trace=std::max(trace,std::abs(double(flux.normalFluxDensity(f.nativeRow,half,false,t))-expected));
        }
        if(f.outwardSign&&f.area>0&&boundary.faces()[f.nativeRow].kind==Boundary::Kind::NormalVelocity){
            const int row=space.faces()[f.nativeRow].first;
            inletMean=std::max(inletMean,std::abs(retained.nativeMoments()[row]-boundary.values()[row]));
            inletLinear=std::max(inletLinear,std::abs(retained.nativeMoments()[row+1]-boundary.values()[row+1]));
        }
    }
    const auto coarseMoments=flow.space().interpolateVelocity(mesh,domain,[&](Vector2D p){const auto u=retained.firstJet(p).value;return Eigen::Vector2d(u.x,u.y);});
    const double moments=(coarseMoments-flow.state().velocity).norm()/flow.state().velocity.norm();
    const double nativeDiv=(space.divergence()*retained.nativeMoments()).norm();
    const double normDifference=double(std::abs(newNorm-oldNorm)/oldNorm);
    const auto rebuiltMoments=space.interpolateVelocity(target,targetDomain,[&](Vector2D p){const auto u=rebuilt.firstJet(p).value;return Eigen::Vector2d(u.x,u.y);});
    const double rebuiltMomentError=(rebuiltMoments-retained.nativeMoments()).norm()/retained.nativeMoments().norm();
    if(const char* prefix=std::getenv("QUANTUM_MOMENT_REFINEMENT_OUTPUT")){
        std::ofstream out(std::string(prefix)+(rz?"_rz_inlet.csv":"_xy_inlet.csv"));require(bool(out),"cannot write native refinement profile");
        out<<std::setprecision(17)<<"x,old_normal,retained_normal,rebuilt_normal,prescribed_mean\n";
        for(int k=0;k<=400;k++){const Vector2D p{.5*k/400,0};out<<p.x<<','<<original.firstJet(p).value.y<<','<<retained.firstJet(p).value.y<<','<<rebuilt.firstJet(p).value.y<<",0.1\n";}
    }
    std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_REFINEMENT metric="<<(rz?"RZ":"XY")<<" cells="<<mesh.cellIndexNum()<<" target_cells="<<target.cellIndexNum()
        <<" value="<<value<<" jet="<<jet<<" coarse_moments="<<moments<<" point_divergence="<<divergence<<" native_divergence="<<nativeDiv
        <<" normal_trace="<<trace<<" same_quadrature_norm_difference="<<normDifference<<" rebuilt_moment_error="<<rebuiltMomentError
        <<" lost_value="<<lostValue<<" lost_jet="<<lostJet<<" child_inlet_mean_mismatch="<<inletMean<<" child_inlet_linear_mismatch="<<inletLinear
        <<" old_native_energy="<<.5*flow.state().velocity.dot(flow.material().kineticMass()*flow.state().velocity)
        <<" refined_native_energy="<<.5*retained.nativeMoments().dot(material.kineticMass()*retained.nativeMoments())<<" accepted=0"<<std::endl;
    require(value<2e-10&&jet<2e-8&&moments<2e-10&&divergence<2e-9&&nativeDiv<2e-11&&trace<2e-11&&normDifference<2e-11&&rebuiltMomentError<2e-9
        &&retained.time()==flow.state().time,"native polynomial restriction lost the original physical field or its moments");
    if(check==NativeTransferCheck::Material||check==NativeTransferCheck::State||check==NativeTransferCheck::PhysicalComplement){
        const auto coarse=original.nativeMaterialForms(mesh,domain,flow.space(),flow.material());
        const auto fine=retained.nativeMaterialForms(target,targetDomain,space,material);
        const double energy=std::abs(coarse.kinetic-fine.kinetic)/coarse.kinetic;
        const double dissipation=std::abs(coarse.dissipation-fine.dissipation)/coarse.dissipation;
        const double momentum=std::hypot(coarse.physicalMomentum.x-fine.physicalMomentum.x,coarse.physicalMomentum.y-fine.physicalMomentum.y)
            /std::hypot(coarse.physicalMomentum.x,coarse.physicalMomentum.y);
        double affineEnergy=0,affineDissipation=0,constantEnergy=0,constantDissipation=0,symmetry=0;
        for(const auto& forms:fine.cells){
            const auto& c=target.cells(target.cellIndexToLabel[forms.index]);const auto h=c.size(),center=c.centerCoord(),lo=center-h*.5;
            Eigen::VectorXd affine=Eigen::VectorXd::Zero(81),constant=affine;
            // Independent exact streamfunctions: XY psi=x*y, psi=-x;
            // RZ Psi=r^2*z, Psi=-r^2/2. Their velocities are affine and
            // have exact kinetic/dissipation integrals from phase moments.
            for(int x=0;x<2;x++)for(int y=0;y<2;y++){
                const double r=lo.x+x*h.x,z=lo.y+y*h.y;
                for(int a=0;a<3;a++)for(int b=0;b<3;b++){
                    const int row=9*(3*x+a)+3*y+b;
                    if(rz){const double radial=a==0?r*r:(a==1?2*r*h.x:2*h.x*h.x);
                        affine[row]=radial*(b==0?z:(b==1?h.y:0));constant[row]=b==0?-.5*radial:0;
                    }else{affine[row]=(a==0?r:(a==1?h.x:0))*(b==0?z:(b==1?h.y:0));constant[row]=b==0?(a==0?-r:(a==1?-h.x:0)):0;}
                }
            }
            const auto& owned=material.cells()[forms.index];
            Eigen::Vector3d vx,vy;vx<<center.x,h.x,0;vy<<(rz?-2.:-1.)*center.y,0,(rz?-2.:-1.)*h.y;
            const double expected=.5*(vx.dot(owned.density*vx)+vy.dot(owned.density*vy));
            affineEnergy=std::max(affineEnergy,std::abs(.5*affine.dot(forms.mass*affine)-expected)/expected);
            const double expectedStrain=(rz?12.:4.)*owned.viscosity(0,0);
            affineDissipation=std::max(affineDissipation,std::abs(affine.dot(forms.viscosity*affine)-expectedStrain)/expectedStrain);
            constantEnergy=std::max(constantEnergy,std::abs(.5*constant.dot(forms.mass*constant)-.5*owned.density(0,0))/(.5*owned.density(0,0)));
            constantDissipation=std::max(constantDissipation,std::abs(constant.dot(forms.viscosity*constant))/expectedStrain);
            symmetry=std::max({symmetry,(forms.mass-forms.mass.transpose()).norm()/forms.mass.norm(),(forms.viscosity-forms.viscosity.transpose()).norm()/forms.viscosity.norm()});
        }
        std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_MATERIAL metric="<<(rz?"RZ":"XY")<<" coarse_energy="<<coarse.kinetic<<" fine_energy="<<fine.kinetic
            <<" energy_difference="<<energy<<" dissipation_difference="<<dissipation<<" momentum_difference="<<momentum<<" phase_moments="<<std::max(coarse.maximumPhaseMomentError,fine.maximumPhaseMomentError)
            <<" quadrature_estimate="<<std::max(coarse.quadratureEstimate,fine.quadratureEstimate)<<" coarse_points="<<coarse.quadraturePoints<<" fine_points="<<fine.quadraturePoints
            <<" affine_energy="<<affineEnergy<<" affine_dissipation="<<affineDissipation<<" constant_energy="<<constantEnergy<<" constant_dissipation="<<constantDissipation<<" symmetry="<<symmetry
            <<" source_coordinates="<<coarse.mass.rows()<<" target_coordinates="<<fine.mass.rows()<<" accepted=0"<<std::endl;
        require(energy<2e-10&&dissipation<2e-10&&momentum<2e-10&&affineEnergy<2e-10&&affineDissipation<2e-9&&constantEnergy<2e-10&&constantDissipation<2e-9&&symmetry<2e-14,
            "complete native material forms lost physical refinement or exact affine integration");
        if(check==NativeTransferCheck::State||check==NativeTransferCheck::PhysicalComplement){
            const auto sourceState=original,fineState=retained,minimalState=rebuilt.nativeVelocity();
            if(check==NativeTransferCheck::PhysicalComplement){physicalComplementControl(target,targetDomain,space,boundary,material,fineState,fine);physicalMidpointControl(flow,mesh,domain,sourceState);}
            const auto detail=fineState.withCoordinates(fineState.coordinates()-minimalState.coordinates(),fineState.time());
            const double momentDifference=(fineState.nativeMoments()-retained.nativeMoments()).norm()/retained.nativeMoments().norm();
            const double zeroMoments=detail.nativeMoments().norm()/retained.nativeMoments().norm();
            const auto& map=fineState.momentMap();
            const double commute=(space.divergence()*map).norm()/(space.divergence().norm()*map.norm());
            double stateValue=0,stateJet=0,work=0;std::vector<Spline::PointForce> forces;
            for(int index:targetDomain.cellIndexByRow()){
                const auto& cell=target.cells(target.cellIndexToLabel[index]);const auto center=cell.centerCoord(),h=cell.size();
                for(double x:{-.37,.03,.31})for(double y:{-.29,.11,.38}){
                    const Vector2D p{center.x+x*h.x,center.y+y*h.y};const auto expected=retained.firstJet(p),actual=fineState.firstJet(p),v=detail.firstJet(p);
                    stateValue=std::max(stateValue,std::hypot(expected.value.x-actual.value.x,expected.value.y-actual.value.y));
                    stateJet=std::max({stateJet,std::hypot(expected.dx.x-actual.dx.x,expected.dx.y-actual.dx.y),std::hypot(expected.dy.x-actual.dy.x,expected.dy.y-actual.dy.y)});
                    const double w=targetDomain.cellVolumes()[targetDomain.cellRowByIndex()[index]]/9,hh=h.x*h.y;
                    forces.push_back({p,v.value*w,v.dx*(w*hh),v.dy*(w*hh)});
                    work+=w*(v.value.x*v.value.x+v.value.y*v.value.y+hh*(v.dx.x*v.dx.x+v.dx.y*v.dx.y+v.dy.x*v.dy.x+v.dy.y*v.dy.y));
                }
            }
            const auto load=detail.transposePointForces(forces);
            const double virtualWork=std::abs(detail.coordinates().dot(load)-work)/work;
            const double oldMomentWork=std::abs(detail.nativeMoments().dot(rebuilt.transposePointForces(forces)))/work;
            const double coarseEnergy=.5*sourceState.coordinates().dot(coarse.mass*sourceState.coordinates());
            const double fineEnergy=.5*fineState.coordinates().dot(fine.mass*fineState.coordinates());
            const double coordinateEnergy=std::max(std::abs(coarseEnergy-coarse.kinetic)/coarse.kinetic,std::abs(fineEnergy-fine.kinetic)/fine.kinetic);
            const double detailEnergy=.5*detail.coordinates().dot(fine.mass*detail.coordinates());
            bool wallDenied=false;const auto& constraints=fineState.kinematicConstraints();
            require(constraints.rows()>0,"complete state control has no physical kinematic conditions");
            Eigen::VectorXd bad=fineState.coordinates();bad+=Eigen::MatrixXd(constraints).row(0).transpose();
            try{auto invalid=fineState.withCoordinates(bad,fineState.time());}catch(const std::invalid_argument&){wallDenied=true;}
            std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_STATE metric="<<(rz?"RZ":"XY")<<" coordinate_value="<<stateValue<<" coordinate_jet="<<stateJet
                <<" retained_moments="<<momentDifference<<" detail_moments="<<zeroMoments<<" divergence_commutator="<<commute<<" complete_work="<<work
                <<" virtual_work="<<virtualWork<<" old_moment_work_fraction="<<oldMomentWork<<" coordinate_energy="<<coordinateEnergy<<" detail_energy="<<detailEnergy
                <<" kinematic_violation_rejected="<<wallDenied<<" accepted=0"<<std::endl;
            require(stateValue<2e-10&&stateJet<2e-8&&momentDifference<2e-10&&zeroMoments<2e-10&&commute<2e-12&&virtualWork<2e-10
                &&oldMomentWork<2e-10&&coordinateEnergy<2e-10&&detailEnergy>0&&work>0&&wallDenied,
                "complete coordinate state lost detail, material energy, kinematic conditions or work");
            for(bool closed:{false,true}){
                auto fixed=boundary.fixed();if(closed)for(const auto f:targetDomain.faces(target))if(f.outwardSign&&f.area>0)for(int k=0;k<2;k++)fixed[space.faces()[f.nativeRow].first+k]=true;
                const auto lift=space.divergenceComplement(target,targetDomain,fixed);
                Eigen::VectorXd eta(lift.moments.cols()),pressure(lift.divergence.rows()),u=Eigen::VectorXd::Zero(space.mass().rows());
                for(int k=0;k<eta.size();k++)eta[k]=std::sin(.31*k+.19);for(int k=0;k<pressure.size();k++)pressure[k]=std::cos(.27*k+.37);
                for(int k=0;k<u.size();k++)if(!fixed[k])u[k]=std::sin(.73*k+.13);
                const Eigen::VectorXd d=space.divergence()*u;Eigen::VectorXd independent(eta.size());
                independent=Eigen::MatrixXd(lift.divergence).colPivHouseholderQr().solve(d);
                const Eigen::VectorXd solenoidal=u-lift.moments*independent;
                const double decomposition=(space.divergence()*solenoidal).norm()/d.norm();
                const double power=std::abs(pressure.dot(space.divergence()*(lift.moments*eta))+eta.dot(lift.pressureLoad*pressure))
                    /(pressure.norm()*(lift.divergence*eta).norm());
                std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_PRESSURE_COMPLEMENT metric="<<(rz?"RZ":"XY")<<" closed="<<closed<<" pressure_rows="<<pressure.size()
                    <<" complement_coordinates="<<eta.size()<<" commuting="<<lift.commutingResidual<<" decomposition="<<decomposition<<" adjoint_work="<<power<<std::endl;
                require(eta.size()==pressure.size()-(closed?1:0)&&lift.commutingResidual<2e-12&&decomposition<2e-12&&power<2e-12,
                    "complete velocity pressure complement changed native conservation or its adjoint");
            }
        }
    }
}
static void jointMomentMeshTransfer(bool open,bool through=false,NativeTransferCheck check=NativeTransferCheck::None){
    const bool exactRefinement=check!=NativeTransferCheck::None;
    constexpr double dt=.002;std::ofstream history;
    if(const char* file=std::getenv("QUANTUM_MOMENT_TRANSFER_OUTPUT")){history.open(std::string(file)+(through?".through.csv":""));require(bool(history),"cannot write mesh-transfer history");
        history<<std::setprecision(17)<<"rz,operation,step,time,cells,volume,kinetic,transfer_energy,equation,fraction\n";}
    const char* metric=std::getenv("QUANTUM_NATIVE_REMESH_METRIC");
    require(!metric||std::string(metric)=="XY"||std::string(metric)=="RZ","unknown complete remesh metric selection");
    for(bool rz:{false,true}){
        if(metric&&std::string(metric)!=(rz?"RZ":"XY"))continue;
        if(through&&!rz)continue;
        const bool wet=open&&rz&&!through;
        auto dictionary=input(4);dictionary.minLevel_neg=dictionary.minLevel_pos={0};dictionary.maxLevel_local={1};
        auto mesh=std::make_unique<QuadtreeMesh>(dictionary);FaceSegmentGraph graph;graph.build(*mesh);
        auto domain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*mesh,graph,rz));auto space=Space::build(*mesh,*domain);
        auto interface=InterfaceParametricMotion::initialize(*mesh,*domain,rz?(wet?InterfacePrescribedData::parametricSupplyCap(.5,.3,8):InterfacePrescribedData::parametricSpheroid(.5,{.22,.17},16)):InterfacePrescribedData::parametricEllipse({.48,.52},{.18,.14},16));
        std::vector<HydroBoundaryPatch> ports;
        if(open){ports={{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,through?1.:.5,-.1},{"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
            ports[0].incomingPhase=wet?-1:1;ports[1].incomingPhase=1;ports[1].incomingVelocity=Vector2D{0,0};}
        Eigen::VectorXd initial=Eigen::VectorXd::Zero(space.mass().rows());
        if(through)initial=space.interpolateVelocity(*mesh,*domain,[](Vector2D){return Eigen::Vector2d(0,.1);});
        auto flow=Flow::createTwoPhase(*mesh,*domain,std::move(space),std::move(interface),{1,.1,through?0.:.01,through?0.:.001},initial,through?0:(open?.1:1),ports);
        Flow::Options options;const auto liquid=[&]{const auto& c=flow.interfaceState().parametricGeometry();return std::accumulate(c.inventory().begin(),c.inventory().end(),0.);};
        if(wet)options.interfaceProjection.materialSegments=8;
        const double supplyRate=wet?std::acos(-1.)*.25*.1:0;
        const auto volumeError=[&](double initial){return std::abs(liquid()-initial-supplyRate*flow.state().time);};
        const auto pressureIntegral=[&]{double integral=0;for(int c=0;c<mesh->cellIndexNum();c++){
            const auto& cell=mesh->cells(mesh->cellIndexToLabel[c]);const double v=domain->cellVolumes()[domain->cellRowByIndex()[c]];
            integral+=v*(flow.state().pressure[3*c]+(rz?cell.size().x/(12*cell.centerCoord().x)*flow.state().pressure[3*c+1]:0));}return integral;};
        const auto verifyPorts=[&]{if(!open)return;double maximumNormalError=0;
            require(flow.boundary().open(),"joint remap lost open boundary ownership");int inlet=0,outlet=0;
            for(const auto f:domain->faces(*mesh))if(f.outwardSign&&f.area>0){const auto& b=flow.boundary().faces()[f.nativeRow];const auto row=flow.space().faces()[f.nativeRow].first;
                if(b.kind==Quantum::Hydrodynamics::HydroCompatibleMomentBoundary::Kind::NormalVelocity){++inlet;require(b.incomingPhase==(wet?-1:1)&&b.incomingVelocity&&b.incomingVelocity->y==.1,"remap lost inlet material/momentum reservoir");}
                if(b.kind==Quantum::Hydrodynamics::HydroCompatibleMomentBoundary::Kind::PressureOutlet){++outlet;require(b.pressure==.2&&b.incomingPhase==1&&b.incomingVelocity&&b.incomingVelocity->x==0&&b.incomingVelocity->y==0&&!flow.boundary().fixed()[row],"remap lost physical pressure outlet or backflow reservoir");}
                for(const auto* u:{&flow.state().velocity,&flow.state().previousVelocity,&flow.state().midpointVelocity})if(u->size())for(int k=0;k<2;k++)if(flow.boundary().fixed()[row+k]){
                    require(std::isfinite((*u)[row+k]),"nonfinite prescribed normal history");
                    maximumNormalError=std::max(maximumNormalError,std::abs((*u)[row+k]-flow.boundary().values()[row+k]));}}
            // U is the integrated moment of the evolved complete field.
            // Check the physical trace, not bitwise equality produced by
            // overwriting its derived cache with the input boundary value.
            std::cout<<std::setprecision(17)<<"MOMENT_BOUNDARY_HISTORY metric="<<(rz?"RZ":"XY")<<" step="<<flow.state().acceptedSteps
                <<" maximum_normal_error="<<maximumNormalError<<std::endl;
            require(maximumNormalError<2e-12,"remap violated a prescribed normal history");
            require(inlet>0&&outlet>0,"remap removed a physical port");
            if(wet){const auto p=flow.interfaceState().parametricGeometry().geometry().sample(0);require(p.x==.5&&p.y==0,"joint remap moved the fixed wet rim");}
            if(through){const auto exact=flow.space().interpolateVelocity(*mesh,*domain,[](Vector2D){return Eigen::Vector2d(0,.1);});
                require((flow.state().velocity-exact).norm()<1e-10*exact.norm(),"joint RZ transfer lost uniform two-phase throughflow");}
        };
        const double initialVolume=liquid();double maximumVolume=0,maximumEquation=0,maximumMomentum=0,maximumFraction=0;
        const auto advance=[&]{const HydroStepTime time{flow.state().acceptedSteps+1,flow.state().time,dt,HydroTimePolicy::Fixed};
            auto load=flow.makeLoad(*mesh,*domain,Eigen::VectorXd::Zero(flow.space().mass().rows()),time.begin+.5*dt);
            auto candidate=flow.prepareNavierStokesStep(*mesh,*domain,time,load,options);const auto report=candidate.report();flow.accept(*mesh,*domain,std::move(candidate));
            maximumVolume=std::max(maximumVolume,volumeError(initialVolume));maximumEquation=std::max(maximumEquation,report.equationResidual);verifyPorts();
            maximumFraction=std::max(maximumFraction,flow.interfaceState().parametricGeometry().report().maximumInventoryFractionError);
            if(history){history<<rz<<",step,"<<flow.state().acceptedSteps<<','<<flow.state().time<<','<<mesh->cellIndexNum()<<','<<liquid()<<','<<report.kineticAfter<<",0,"<<report.equationResidual<<','<<maximumFraction<<'\n';history.flush();}
            std::cout<<std::setprecision(17)<<"MOMENT_MESH_TRANSFER_STEP metric="<<(rz?"RZ":"XY")<<" step="<<flow.state().acceptedSteps<<" time="<<flow.state().time<<" cells="<<mesh->cellIndexNum()<<" kinetic="<<report.kineticAfter<<" equation="<<report.equationResidual<<std::endl;
        };
        advance();
        for(bool refine:{true,false}){
            const auto oldGeometry=flow.interfaceState().parametricGeometry().geometryHandle();const auto oldState=flow.state();const auto oldInventory=flow.interfaceState().parametricGeometry().inventory();
            const double oldPressureIntegral=pressureIntegral(),oldInitialCorrection=flow.initialBoundaryCorrection();
            auto target=mesh->clone();CellField<double> gradient(*target,0.),levels(*target,0.);
            for(int c=0;c<target->cellIndexNum();c++){const auto p=target->cells(target->cellIndexToLabel[c]).centerCoord();if(refine&&p.x<.5&&(open||p.y>.5))levels.value[c]=1;}
            std::vector<VertField<double>> scalars;std::vector<VertField<double>*> vertices;std::vector<HEdgeField<double>*> horizontal;
            std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> fields{&levels};target->adapt(scalars,gradient,vertices,horizontal,vertical,fields,&levels);
            require(target->cellIndexNum()==(refine?(open?40:28):16),"joint remap did not execute the actual refiner/coarsener");
            FaceSegmentGraph targetGraph;targetGraph.build(*target);auto targetDomain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*target,targetGraph,rz));
            if(exactRefinement){
                nativeRefinementControl(flow,*mesh,*domain,*target,*targetDomain,ports,check);
                require(flow.interfaceState().parametricGeometry().geometryHandle()==oldGeometry&&flow.interfaceState().parametricGeometry().inventory()==oldInventory
                    &&flow.state().velocity==oldState.velocity&&flow.state().previousVelocity==oldState.previousVelocity&&flow.state().midpointVelocity==oldState.midpointVelocity
                    &&flow.state().pressure==oldState.pressure&&flow.state().viscousTrace==oldState.viscousTrace&&flow.state().time==oldState.time&&flow.state().acceptedSteps==oldState.acceptedSteps,
                    "native refinement source changed the accepted physical state");
                break;
            }
            if(wet){
                // The original supplied-cap attempt is retained as a failed
                // physical target, not tuned into a positive control. Its
                // constrained minimum energy exceeds the old kinetic bound.
                bool denied=false;std::string reason;
                try{auto candidate=flow.prepareMeshTransfer(*mesh,*domain,*target,*targetDomain,options);}
                catch(const std::invalid_argument& e){reason=e.what();denied=reason.find("moment remap energy/conservation infeasible:")==0;}
                require(denied&&flow.interfaceState().parametricGeometry().geometryHandle()==oldGeometry&&flow.interfaceState().parametricGeometry().inventory()==oldInventory
                    &&flow.state().velocity==oldState.velocity&&flow.state().previousVelocity==oldState.previousVelocity&&flow.state().midpointVelocity==oldState.midpointVelocity
                    &&flow.state().pressure==oldState.pressure&&flow.state().viscousTrace==oldState.viscousTrace&&flow.state().time==oldState.time&&flow.state().acceptedSteps==oldState.acceptedSteps,
                    "infeasible supplied RZ mesh target was accepted or mutated the state");
                std::cout<<"MOMENT_PORT_TRANSFER_REJECTED metric=RZ cells="<<mesh->cellIndexNum()<<" target_cells="<<target->cellIndexNum()<<" unchanged=1 reason="<<reason<<std::endl;
                advance();break;
            }
            auto candidate=flow.prepareMeshTransfer(*mesh,*domain,*target,*targetDomain,options);const auto report=candidate.report();
            require(flow.interfaceState().parametricGeometry().geometryHandle()==oldGeometry&&flow.interfaceState().parametricGeometry().inventory()==oldInventory&&flow.state().velocity==oldState.velocity,"joint remap preparation changed accepted state");
            bool foreign=false;try{flow.acceptMeshTransfer(*mesh,*domain,*target,*domain,std::move(candidate));}catch(const std::exception&){foreign=true;}
            require(foreign&&flow.interfaceState().parametricGeometry().geometryHandle()==oldGeometry,"foreign target remap changed accepted state");
            flow.acceptMeshTransfer(*mesh,*domain,*target,*targetDomain,std::move(candidate));
            mesh=std::move(target);domain=std::move(targetDomain);
            verifyPorts();require(std::abs(pressureIntegral()-oldPressureIntegral)<5e-12*std::max(1.,std::abs(oldPressureIntegral))&&flow.initialBoundaryCorrection()==oldInitialCorrection,"joint remap changed physical pressure mean or initialization provenance");
            require(flow.interfaceState().parametricGeometry().geometryHandle()==oldGeometry&&flow.state().time==oldState.time&&flow.state().previousTime==oldState.previousTime
                &&flow.state().acceptedSteps==oldState.acceptedSteps&&flow.state().lastDt==oldState.lastDt&&flow.state().stageTime==oldState.stageTime
                &&flow.state().continuousHistory==oldState.continuousHistory&&flow.state().remappedStage,"joint remap changed physical clock, curve or stage provenance");
            require(report.levels.size()==3&&flow.interfaceState().progress().completedSteps==oldState.acceptedSteps&&flow.interfaceState().progress().endTime==oldState.time,"joint remap lost Interface/Hydro history");
            bool consumed=false;try{flow.acceptMeshTransfer(*mesh,*domain,*mesh,*domain,std::move(candidate));}catch(const std::invalid_argument&){consumed=true;}require(consumed,"joint remap accepted a consumed transaction");
            maximumVolume=std::max(maximumVolume,volumeError(initialVolume));maximumFraction=std::max(maximumFraction,report.interface.maximumFractionError);
            for(const auto& level:report.levels){maximumMomentum=std::max(maximumMomentum,level.momentumResidual);
                require(level.kineticAfter<=level.kineticBefore*(1+1e-12),"joint remap created kinetic energy");
                std::cout<<"MOMENT_MESH_TRANSFER_LEVEL metric="<<(rz?"RZ":"XY")<<" refine="<<refine<<" time="<<level.time<<" patches="<<level.patches<<" momentum="<<level.momentumResidual<<" divergence="<<level.divergenceResidual<<" kinetic_before="<<level.kineticBefore<<" kinetic_after="<<level.kineticAfter<<" minimum_kinetic="<<level.minimumKinetic<<" unconstrained_kinetic="<<level.unconstrainedKinetic<<" retained_nullspace="<<level.retainedNullspace<<" energy_limited="<<level.energyLimited<<" correction="<<level.correctionNorm<<std::endl;}
            const double energy=report.levels.front().kineticAfter,change=energy-report.levels.front().kineticBefore;
            if(history){history<<rz<<','<<(refine?"refine":"coarsen")<<','<<flow.state().acceptedSteps<<','<<flow.state().time<<','<<mesh->cellIndexNum()<<','<<liquid()<<','<<energy<<','<<change<<",0,"<<report.interface.maximumFractionError<<'\n';history.flush();}
            std::cout<<"MOMENT_MESH_TRANSFER metric="<<(rz?"RZ":"XY")<<" refine="<<refine<<" cells="<<mesh->cellIndexNum()<<" inventory_balance="<<report.interface.inventoryBalance<<" fraction="<<report.interface.maximumFractionError<<" geometry_identical=1 clock_identical=1 consumed_rejected="<<consumed<<std::endl;
            advance();
        }
        if(exactRefinement)continue;
        require(maximumVolume<2e-13&&maximumEquation<=1e-10&&maximumMomentum<2e-10&&maximumFraction<3.1e-12&&flow.state().acceptedSteps==(wet?2:3),"joint mesh transfer or subsequent physical step exceeded original budgets");
        std::cout<<"MOMENT_MESH_TRANSFER_COMPLETE metric="<<(rz?"RZ":"XY")<<" open="<<open<<" through="<<through<<" steps="<<flow.state().acceptedSteps<<" time="<<flow.state().time<<" refined_and_coarsened="<<!wet<<" volume_drift="<<maximumVolume<<" fraction="<<maximumFraction<<" momentum="<<maximumMomentum<<" equation="<<maximumEquation<<std::endl;
    }
}
void hydroMomentMeshTransferContract(){jointMomentMeshTransfer(false);}
void hydroMomentPortMeshTransferContract(){jointMomentMeshTransfer(true);jointMomentMeshTransfer(true,true);}
void hydroMomentNativeRefinementContract(){jointMomentMeshTransfer(true,false,NativeTransferCheck::Restriction);}
void hydroMomentCartesianMaterialGramContract(){
    // Independent dense T^T B T reference: shared/hanging-style coordinates,
    // an unmapped local row, empty global columns, nonsymmetric forms and
    // signed constant-mode weights. Reuse the same plan with changed B.
    {
        const std::vector<int> offsets{0,5,11,15};std::vector<Eigen::Triplet<double>> entries;
        for(int row=0;row<15;row++)if(row!=4){
            entries.emplace_back(row,row%7,.7+.01*row);
            if(row%3==0)entries.emplace_back(row,(row+2)%7,-.3);
        }
        Space::Matrix map(15,9);map.setFromTriplets(entries.begin(),entries.end());
        const Quantum::Hydrodynamics::MomentSplineCoordinateAssembly plan(map,offsets);
        std::vector<Eigen::MatrixXd> blocks;
        for(int slot=0;slot<3;slot++){
            const int n=offsets[slot+1]-offsets[slot];Eigen::MatrixXd b(n,n);
            for(int i=0;i<n;i++)for(int j=0;j<n;j++)b(i,j)=std::sin(.3*(i+1)+.17*(j+1)+slot)+(i==j);
            blocks.push_back(std::move(b));
        }
        const std::vector<double> weights{1.,-.2,0.};
        for(int pass=0;pass<3;pass++){
            if(pass==1)for(auto& b:blocks)b.array()+=.13;
            if(pass==2)for(auto& b:blocks)b.setZero();
            const auto* scales=pass==1?&weights:nullptr;
            Eigen::MatrixXd local=Eigen::MatrixXd::Zero(15,15);
            for(int slot=0;slot<3;slot++)local.block(offsets[slot],offsets[slot],blocks[slot].rows(),blocks[slot].cols())
                =(scales?weights[slot]:1.)*blocks[slot];
            const Eigen::MatrixXd expected=Eigen::MatrixXd(map.transpose())*local*Eigen::MatrixXd(map);
            const auto actual=plan.assemble([&](int slot)->const Eigen::MatrixXd&{return blocks[slot];},scales);
            const double error=(Eigen::MatrixXd(actual)-expected).norm()/(1+expected.norm());
            require(error<3e-14&&actual.rows()==9&&actual.cols()==9,"cached coordinate assembly differs from full dense contraction");
            require(pass!=2||actual.nonZeros()==0,"zero local forms acquired nonzero global storage");
            std::cout<<std::setprecision(17)<<"MOMENT_COORDINATE_ASSEMBLY pass="<<pass<<" relative_error="<<error<<" nonzeros="<<actual.nonZeros()<<std::endl;
        }
    }
    using Gram=Quantum::Hydrodynamics::MomentSplineBasis::CartesianGram;
    const Quantum::Hydrodynamics::MomentSplineBasis::GaussTwelve rule;
    // Independent derivative samples of a complete tensor polynomial. The
    // reference forms dense physical velocity/strain rows at EVERY node.
    // Include a phase jump, spatially varying viscosity, anisotropic cells,
    // and a batch boundary inside one horizontal slice.
    const auto axis=[](double x,double h){Gram::Axes a;
        for(int k=0;k<9;k++)for(int d=0;d<3;d++){
            double value=k<d?0:std::pow(x,k-d);for(int j=0;j<d;j++)value*=double(k-j)/h;a(d,k)=value;
        }return a;};
    for(double hx:{1.,1./64})for(bool viscous:{false,true}){
        const double hy=.3*hx;Eigen::MatrixXd mass=Eigen::MatrixXd::Zero(81,81),viscosity=mass;
        Eigen::MatrixXd velocityRows(2*12*24,81),strainRows=Eigen::MatrixXd::Zero(3*12*24,81);
        int node=0;std::optional<Gram> gram;
        for(int j=0;j<12;j++){
            const double y=.5*(1+rule.x[j]);const auto by=axis(y,hy);gram.emplace(by);
            for(int phase=0;phase<2;phase++)for(int i=0;i<12;i++){
                const double left=phase?.37:0,right=phase?1:.37,x=left+.5*(right-left)*(1+rule.x[i]);
                const auto ax=axis(x,hx);const double w=.25*hx*hy*(right-left)*rule.w[i]*rule.w[j];
                const double rho=phase?.001:1.,mu=viscous?(phase?.01:1.)*(1+x*x+.3*x*y+.2*y*y):0;
                gram->add(ax,w*rho,w*mu);
                for(int a=0;a<9;a++)for(int b=0;b<9;b++){
                    const int k=9*a+b;
                    velocityRows(2*node,k)=std::sqrt(w*rho)*ax(0,a)*by(1,b);
                    velocityRows(2*node+1,k)=-std::sqrt(w*rho)*ax(1,a)*by(0,b);
                    strainRows(3*node,k)=std::sqrt(2*w*mu)*ax(1,a)*by(1,b);
                    strainRows(3*node+1,k)=-std::sqrt(2*w*mu)*ax(1,a)*by(1,b);
                    strainRows(3*node+2,k)=std::sqrt(w*mu)*(ax(0,a)*by(2,b)-ax(2,a)*by(0,b));
                }
                if(i==4){gram->addTo(mass,viscosity);gram.emplace(by);}
                ++node;
            }
            gram->addTo(mass,viscosity);
        }
        const Eigen::MatrixXd referenceMass=velocityRows.transpose()*velocityRows,referenceViscosity=strainRows.transpose()*strainRows;
        const auto scaled=[](const Eigen::MatrixXd& value,const Eigen::MatrixXd& reference){double result=0;
            for(int i=0;i<81;i++)for(int j=0;j<81;j++){
                const double scale=std::sqrt(reference(i,i)*reference(j,j)),difference=std::abs(value(i,j)-reference(i,j));
                if(scale>0)result=std::max(result,difference/scale);else require(difference==0,"tensor Gram changed a zero polynomial direction");
            }return result;};
        const double massError=scaled(mass,referenceMass),viscosityError=scaled(viscosity,referenceViscosity);
        const double symmetry=std::max((mass-mass.transpose()).norm()/mass.norm(),viscous?(viscosity-viscosity.transpose()).norm()/viscosity.norm():0.);
        double energyError=0;
        for(int direction=0;direction<8;direction++){
            Eigen::VectorXd q(81);for(int k=0;k<81;k++)q[k]=std::sin((k+1)*(.31+.07*direction));
            const double kinetic=(velocityRows*q).squaredNorm(),dissipation=(strainRows*q).squaredNorm();
            energyError=std::max(energyError,std::abs(q.dot(mass*q)-kinetic)/kinetic);
            if(viscous)energyError=std::max(energyError,std::abs(q.dot(viscosity*q)-dissipation)/dissipation);
        }
        std::cout<<std::setprecision(17)<<"MOMENT_CARTESIAN_GRAM hx="<<hx<<" hy="<<hy<<" viscous="<<viscous
            <<" points="<<node<<" mass="<<massError<<" viscosity="<<viscosityError<<" symmetry="<<symmetry<<" energy="<<energyError<<std::endl;
        require(massError<5e-13&&viscosityError<5e-13&&symmetry<5e-15&&energyError<5e-12,"tensor contraction differs from positive physical quadrature");
    }
}
void hydroMomentNativeMaterialContract(){jointMomentMeshTransfer(true,false,NativeTransferCheck::Material);}
void hydroMomentNativeStateContract(){jointMomentMeshTransfer(true,false,NativeTransferCheck::State);}
void hydroMomentPhysicalComplementContract(){jointMomentMeshTransfer(true,false,NativeTransferCheck::PhysicalComplement);}
void hydroMomentBlockSolveContract(){
    for(bool open:{false,true})for(bool rz:{false,true}){
        QuadtreeMesh mesh(input(4));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);auto space=Space::build(mesh,domain);
        auto interface=InterfaceParametricMotion::initialize(mesh,domain,rz?InterfacePrescribedData::parametricSpheroid(.5,{.22,.17},16)
            :InterfacePrescribedData::parametricEllipse({.48,.52},{.18,.14},16));
        std::vector<HydroBoundaryPatch> ports;
        if(open){ports={{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,.5,-.1},
            {"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
            ports[0].incomingPhase=ports[1].incomingPhase=1;ports[1].incomingVelocity=Vector2D{0,0};}
        auto initial=space.interpolateVelocity(mesh,domain,[&](Vector2D p)->Eigen::Vector2d{
            if(open)return Eigen::Vector2d::Zero();if(!rz){const auto u=exact(p).value;return {u.x,u.y};}
            const double r=p.x,z=p.y,F=(1-r*r)*(1-r*r),Fr=-4*r*(1-r*r),g=z*z*(1-z)*(1-z),gz=2*z*(1-z)*(1-2*z);
            return {r*F*gz,-(2*F+r*Fr)*g};
        });
        const auto flow=Flow::createTwoPhase(mesh,domain,std::move(space),std::move(interface),{1,.1,.01,.001},std::move(initial),0,ports);
        physicalMidpointControl(flow,mesh,domain,*flow.state().physicalVelocity);
    }
}
void hydroMomentNativeFluxContract(){jointMomentMeshTransfer(true,false,NativeTransferCheck::PhaseFlux);}
void hydroMomentMassCurrentContract(){jointMomentMeshTransfer(true,false,NativeTransferCheck::VolumeCurrent);}
void hydroMomentNativeMomentumContract(){jointMomentMeshTransfer(true,false,NativeTransferCheck::FaceMomentum);}
void hydroMomentCapillaryMetricContract(){
    namespace P=Quantum::Geometry::Parametric;
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    using Traction=Quantum::Hydrodynamics::HydroCompatibleMomentTraction;
    using Stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
    using Boundary=Quantum::Hydrodynamics::HydroCompatibleMomentBoundary;
    constexpr double sigma=.7;
    for(bool rz:{false,true})for(bool mixed:{false,true}){
        auto dictionary=input(4);dictionary.minLevel_neg=dictionary.minLevel_pos={0};dictionary.maxLevel_local={1};QuadtreeMesh mesh(dictionary);
        if(mixed){CellField<double> gradient(mesh,0.),levels(mesh,0.);
            for(int c=0;c<mesh.cellIndexNum();c++){const auto p=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();if(p.x<.5&&p.y>.5)levels.value[c]=1;}
            std::vector<VertField<double>> scalars;std::vector<VertField<double>*> vertices;std::vector<HEdgeField<double>*> horizontal;
            std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> fields{&levels};mesh.adapt(scalars,gradient,vertices,horizontal,vertical,fields,&levels);
            require(mesh.cellIndexNum()==28,"capillary metric did not use mixed leaves");}
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto space=Space::build(mesh,domain);
        auto interface=InterfaceParametricMotion::initialize(mesh,domain,rz?InterfacePrescribedData::parametricSpheroid(.5,{.22,.17},64):InterfacePrescribedData::parametricEllipse({.48,.52},{.18,.14},64));
        const auto moments=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const double x=p.x,y=p.y;
            if(rz)return Eigen::Vector2d(2*x*std::pow(1-x*x,2)*y*(1-y)*(1-2*y),-2*(1-x*x)*(1-3*x*x)*y*y*(1-y)*(1-y));
            return Eigen::Vector2d(2*x*x*(1-x)*(1-x)*y*(1-y)*(1-2*y),-2*x*(1-x)*(1-2*x)*y*y*(1-y)*(1-y));});
        const auto velocity=Spline::initializeNative(mesh,domain,space,moments,0);
        const auto material=Material::buildFrozen(mesh,domain,space,interface.parametricGeometry().boundaryHandle(),0,{1,1,0,0});
        const auto metric=Traction::iterationMetric(mesh,domain,space,velocity,material,sigma);
        const Space::Matrix difference=metric-Space::Matrix(metric.transpose());const double symmetry=difference.norm()/metric.norm();
        const double quadratic=velocity.coordinates().dot(metric*velocity.coordinates());
        // Independent centered differences of the actual scalar normal speed,
        // integrated without the metric's support cuts or Gauss quadrature.
        // This does not reuse its analytic curvature/velocity-jet derivative.
        const auto integrate=[&](int samples){long double value=0;
            for(const auto& c:material.boundaryHandle()->components())for(const auto& arc:c.curve->arcs)for(int i=0;i<samples;i++){
                const double t=(i+.5)/samples;const auto at=rz?P::sampleSurfaceRZ(arc,t,c.curve->boundary):P::sampleSurface(arc,t);
                const auto normalSpeed=[&](double s){const auto q=rz?P::sampleSurfaceRZ(arc,s,c.curve->boundary):P::sampleSurface(arc,s);
                    const auto u=velocity.firstJet({q.point.x,q.point.y}).value;return u.x*q.normal.x+u.y*q.normal.y;};
                constexpr double increment=1e-5;
                const double derivative=(normalSpeed(t+increment)-normalSpeed(t-increment))/(2*increment*at.speed),gradient=derivative*derivative;
                value+=sigma*at.speed*(rz?2*std::acos(-1.)*at.point.x:1.)*gradient/samples;
            }return double(value);};
        const double coarse=integrate(128),fine=integrate(256),quadrature=std::abs(fine-coarse)/fine,error=std::abs(quadratic-fine)/fine;
        std::cout<<std::setprecision(17)<<"MOMENT_CAPILLARY_METRIC metric="<<(rz?"RZ":"XY")<<" mixed="<<mixed<<" cells="<<mesh.cellIndexNum()
            <<" quadratic="<<quadratic<<" reference="<<fine<<" relative_error="<<error<<" quadrature_difference="<<quadrature<<" symmetry="<<symmetry<<std::endl;
        require(quadratic>0&&symmetry<1e-12&&quadrature<2e-5&&error<std::max(8*quadrature,2e-9),"capillary metric changed normal shape work or symmetry");
    }
    // Rotation only relabels an exact circle. Its residual shape work must
    // vanish as the SAME Hermite circle is refined, rather than paying the
    // finite vector-gradient penalty of the earlier iteration approximation.
    QuadtreeMesh mesh(input(4));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const auto stress=Stress::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),0));
    const auto boundary=Boundary::build(mesh,domain,space,stress,{
        {"west",HydroBoundarySide::West,HydroBoundaryType::PressureOutlet,0,1,0},
        {"east",HydroBoundarySide::East,HydroBoundaryType::PressureOutlet,0,1,0},
        {"south",HydroBoundarySide::South,HydroBoundaryType::PressureOutlet,0,1,0},
        {"north",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,0}});
    const auto moments=space.interpolateVelocity(mesh,domain,[](Vector2D p){return Eigen::Vector2d(.5-p.y,p.x-.5);});
    const auto velocity=Spline::initializeNative(mesh,domain,space,moments,0,&boundary);
    double previous=0;
    for(int segments:{32,64}){
        auto interface=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.5,.5},{.18,.18},segments));
        const auto material=Material::buildFrozen(mesh,domain,space,interface.parametricGeometry().boundaryHandle(),0,{1,1,0,0});
        const auto metric=Traction::iterationMetric(mesh,domain,space,velocity,material,sigma);
        const double work=std::abs(velocity.coordinates().dot(metric*velocity.coordinates()))/(sigma*2*std::acos(-1.)*.18);
        long double direct=0;
        for(const auto& arc:material.boundaryHandle()->components().front().curve->arcs)for(int i=0;i<256;i++){
            const double t=(i+.5)/256;const auto at=P::sampleSurface(arc,t);
            const auto normalSpeed=[&](double s){const auto q=P::sampleSurface(arc,s);const auto u=velocity.firstJet({q.point.x,q.point.y}).value;
                return u.x*q.normal.x+u.y*q.normal.y;};
            constexpr double increment=1e-5;
            const double derivative=(normalSpeed(t+increment)-normalSpeed(t-increment))/(2*increment*at.speed);
            direct+=at.speed*derivative*derivative/256;
        }
        const double shapeWork=double(direct)/(2*std::acos(-1.)*.18);
        std::cout<<std::setprecision(17)<<"MOMENT_CAPILLARY_REPARAMETRIZATION segments="<<segments<<" matrix_relative_work="<<work<<" direct_relative_work="<<shapeWork<<std::endl;
        // Test decay in the directly integrated nonnegative integrand. Forming
        // q^T K q for this near-null mode cancels large sparse matrix entries.
        require(work<1e-8&&shapeWork<1e-8&&(segments==32||shapeWork<previous/16),"capillary iteration penalized tangential circle rotation");previous=shapeWork;
    }
}
void hydroMomentNativeCoupledContract(){
    constexpr double dt=.002;
    for(bool open:{true,false})for(bool rz:{false,true}){
        QuadtreeMesh mesh(input(4));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);auto space=Space::build(mesh,domain);
        const bool wet=open&&rz;
        auto interface=InterfaceParametricMotion::initialize(mesh,domain,rz?(wet?InterfacePrescribedData::parametricSupplyCap(.5,.3,8):InterfacePrescribedData::parametricSpheroid(.5,{.22,.17},16)):InterfacePrescribedData::parametricEllipse({.48,.52},{.18,.14},16));
        std::vector<HydroBoundaryPatch> ports;
        if(open){ports={{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,.5,-.1},{"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
            ports[0].incomingPhase=wet?-1:1;ports[1].incomingPhase=1;ports[1].incomingVelocity=Vector2D{0,0};}
        const Eigen::VectorXd zero=Eigen::VectorXd::Zero(space.mass().rows());
        auto flow=Flow::createTwoPhase(mesh,domain,std::move(space),std::move(interface),{1,.1,.01,.001},zero,.1,ports);
        Flow::Options options;if(wet)options.interfaceProjection.materialSegments=8;
        const auto liquid=[&]{const auto& inventory=flow.interfaceState().parametricGeometry().inventory();return std::accumulate(inventory.begin(),inventory.end(),0.);};
        const double initialVolume=liquid(),supply=wet?std::acos(-1.)*.025:0;
        require(flow.state().physicalVelocity&&flow.state().previousPhysicalVelocity&&!flow.state().midpointPhysicalVelocity,"complete initial ownership missing");
        for(int step=1;step<=(open?2:1);step++){
            const auto old=flow.state();const auto oldGeometry=flow.interfaceState().parametricGeometry().boundaryHandle();
            const HydroStepTime time{step,old.time,dt,HydroTimePolicy::Fixed};const auto load=flow.makeLoad(mesh,domain,zero,old.time+.5*dt);
            auto candidate=flow.prepareNavierStokesStep(mesh,domain,time,load,options);const auto report=candidate.report();
            require(flow.state().time==old.time&&flow.state().physicalVelocity->coordinates()==old.physicalVelocity->coordinates()
                &&flow.interfaceState().parametricGeometry().boundaryHandle()==oldGeometry,"complete preparation mutated accepted state");
            const auto& next=candidate.state();const auto& q=next.physicalVelocity->coordinates();
            const double midpoint=(next.midpointPhysicalVelocity->coordinates()-.5*(q+old.physicalVelocity->coordinates())).norm();
            const double momentCache=(next.velocity-next.physicalVelocity->nativeMoments()).norm();
            const double workScale=report.kineticBefore+report.kineticAfter+std::abs(report.viscousDissipation)+std::abs(report.capillaryWork)+std::abs(report.boundaryReactionWork)+std::abs(report.boundaryTractionWork);
            const double energy=std::abs(report.energyDefect)/std::max(workScale,1e-300);
            const double spatial=report.materialMetricWork-report.boundaryKineticFlux-report.boundaryInflowDissipation-report.geometricMassWork;
            flow.accept(mesh,domain,std::move(candidate));bool consumed=false;try{flow.accept(mesh,domain,std::move(candidate));}catch(const std::invalid_argument&){consumed=true;}
            const double volume=std::abs(liquid()-initialVolume-supply*flow.state().time);
            const double kineticCache=std::abs(flow.kineticEnergy()-report.kineticAfter)/std::max(report.kineticAfter,1e-300);
            std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_COUPLED metric="<<(rz?"RZ":"XY")<<" open="<<open<<" step="<<step<<" time="<<flow.state().time
                <<" coordinates="<<q.size()<<" iterations="<<report.nonlinearIterations<<" equation="<<report.equationResidual<<" divergence="<<report.divergenceResidual
                <<" boundary="<<report.boundaryConstraintResidual<<" complement="<<report.complementResidual<<" energy="<<energy<<" spatial_work="<<spatial
                <<" kinetic="<<report.kineticAfter<<" capillary="<<report.capillaryWork<<" surface_work="<<report.surfaceWorkDefect
                <<" temporal="<<report.temporalMassWork<<" geometric="<<report.geometricMassWork<<" material_metric="<<report.massMetricDefect
                <<" volume="<<volume<<" midpoint="<<midpoint<<" cache="<<momentCache<<" kinetic_cache="<<kineticCache<<" consumed="<<consumed<<std::endl;
            require(midpoint<2e-12&&momentCache==0&&kineticCache<2e-12&&consumed&&volume<2e-12&&energy<2e-8&&report.equationResidual<=1e-10,
                "complete coupled midpoint lost ownership, material work, volume or equation closure");
        }
    }
}
namespace {
struct NativeManufactured {Spline::FirstJet velocity;Vector2D laplacian,gradient;double pressure;};
NativeManufactured nativeManufactured(Vector2D p,bool rz,bool polynomial){
    const double pi=std::acos(-1.),x=p.x,y=p.y;
    const auto shape=[&](double x){return polynomial?std::array<double,4>{x*x*(1-x)*(1-x),2*x-6*x*x+4*x*x*x,2-12*x+12*x*x,-12+24*x}
        :std::array<double,4>{std::pow(std::sin(pi*x),2),pi*std::sin(2*pi*x),2*pi*pi*std::cos(2*pi*x),-4*pi*pi*pi*std::sin(2*pi*x)};};
    const auto f=shape(x),g=shape(y);const double s=polynomial?8.:.05;
    NativeManufactured out;
    if(!rz){out.velocity={{s*f[0]*g[1],-s*f[1]*g[0]},{s*f[1]*g[1],-s*f[2]*g[0]},{s*f[0]*g[2],-s*f[1]*g[1]}};
        out.laplacian={s*(f[2]*g[1]+f[0]*g[3]),-s*(f[3]*g[0]+f[1]*g[2])};}
    else {const double a=x-2*x*x*x+std::pow(x,5),ar=1-6*x*x+5*std::pow(x,4),b=-2+8*x*x-6*std::pow(x,4),br=16*x-24*x*x*x;
        out.velocity={{s*a*g[1],s*b*g[0]},{s*ar*g[1],s*br*g[0]},{s*a*g[2],s*b*g[1]}};
        out.laplacian={s*((-16*x+24*x*x*x)*g[1]+a*g[3]),s*((32-96*x*x)*g[0]+b*g[2])};}
    if(polynomial){out.pressure=rz?y-.5:x+y-1;out.gradient={rz?0.:1.,1};}
    else {out.pressure=std::cos(pi*x)*std::cos(pi*y);out.gradient={-pi*std::sin(pi*x)*std::cos(pi*y),-pi*std::cos(pi*x)*std::sin(pi*y)};}
    return out;
}
struct NativeAccuracyResult {double velocity=0,pressure=0,gradient=0;};
NativeAccuracyResult nativeAccuracy(int n,bool rz,bool temporal,double dt,bool mixed){
    QuadtreeMesh mesh(input(n));
    if(mixed){DictionaryMesh::RefineBox box;box.min={0,0};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});
        require(mesh.cellIndexNum()==5*n*n/2,"manufactured mixed grid did not refine the declared left half");}
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);auto space=Space::build(mesh,domain);
    require(!mixed||!domain.uniformLeaves(),"manufactured mixed grid has no coarse-fine interface");
    auto interface=InterfaceParametricMotion::initialize(mesh,domain,rz?InterfacePrescribedData::parametricSpheroid(.5,{.18,.16},16):InterfacePrescribedData::parametricEllipse({.46,.53},{.16,.13},16));
    const Eigen::VectorXd zero=Eigen::VectorXd::Zero(space.mass().rows());
    // Equal materials isolate Hydro accuracy but retain the real complete
    // moving-interface owner, phase fluxes and common midpoint code path.
    constexpr double mu=.02;auto flow=Flow::createTwoPhase(mesh,domain,std::move(space),std::move(interface),{1,1,mu,mu},zero,0);
    const double end=temporal?.08:.02;const int steps=int(std::llround(end/dt));Flow::Options options;
    double equation=0,divergence=0,energy=0;int iterations=0;
    for(int step=1;step<=steps;step++){
        const double begin=flow.state().time,mid=begin+.5*dt,a=temporal?std::sin(4*mid):mid,ad=temporal?4*std::cos(4*mid):1.;
        const auto load=flow.makeBodyForceLoad(mesh,domain,[&](Vector2D at){const auto e=nativeManufactured(at,rz,temporal);const auto v=e.velocity;
            return v.value*ad+(v.dx*v.value.x+v.dy*v.value.y)*(a*a)-e.laplacian*(mu*a)+e.gradient;},mid);
        auto next=flow.prepareNavierStokesStep(mesh,domain,{step,begin,dt,HydroTimePolicy::Fixed},load,options);const auto report=next.report();
        equation=std::max(equation,report.equationResidual);divergence=std::max(divergence,report.divergenceResidual);iterations+=report.nonlinearIterations;
        const double scale=report.kineticBefore+report.kineticAfter+std::abs(report.externalWork)+std::abs(report.viscousDissipation);
        energy=std::max(energy,std::abs(report.energyDefect)/std::max(scale,1e-300));flow.accept(mesh,domain,std::move(next));
        std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_ACCURACY_STEP metric="<<(rz?"RZ":"XY")<<" temporal="<<temporal<<(mixed?" mixed=1":"")<<" n="<<n<<" dt="<<dt<<" step="<<step<<" equation="<<report.equationResidual<<std::endl;
    }
    long double ue=0,un=0,pe=0,pn=0,ge=0,gn=0;const double amplitude=temporal?std::sin(4*end):end;const auto& rule=Quantum::Geometry::Parametric::gaussEight();
    std::ofstream output;if(const char* prefix=std::getenv("QUANTUM_NATIVE_ACCURACY_OUTPUT")){
        output.open(std::string(prefix)+(rz?"_rz":"_xy")+"_n"+std::to_string(n)+"_dt"+std::to_string(dt)+".csv");
        output<<std::setprecision(17)<<"x,y,ux,uy,exact_ux,exact_uy,pressure,exact_pressure\n";}
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
        for(int i=0;i<8;i++)for(int j=0;j<8;j++){
            const Vector2D at=center+Vector2D{.5*h.x*rule.x[i],.5*h.y*rule.x[j]};const auto e=nativeManufactured(at,rz,temporal);const auto v=flow.state().physicalVelocity->firstJet(at);
            const double w=.25*h.x*h.y*rule.w[i]*rule.w[j]*(rz?2*std::acos(-1.)*at.x:1.);
            const Eigen::Vector3d basis{1,(at.x-center.x)/h.x,(at.y-center.y)/h.y};const double p=flow.state().pressure.segment<3>(3*c).dot(basis);
            const auto du=v.value-e.velocity.value*amplitude,dx=v.dx-e.velocity.dx*amplitude,dy=v.dy-e.velocity.dy*amplitude;
            const auto square=[](Vector2D a){return a.x*a.x+a.y*a.y;};
            ue+=w*square(du);un+=w*amplitude*amplitude*square(e.velocity.value);pe+=w*(p-e.pressure)*(p-e.pressure);pn+=w*e.pressure*e.pressure;
            ge+=w*(square(dx)+square(dy));gn+=w*amplitude*amplitude*(square(e.velocity.dx)+square(e.velocity.dy));
            if(output)output<<at.x<<','<<at.y<<','<<v.value.x<<','<<v.value.y<<','<<amplitude*e.velocity.value.x<<','<<amplitude*e.velocity.value.y<<','<<p<<','<<e.pressure<<'\n';
        }
    }
    NativeAccuracyResult result{std::sqrt(double(ue/un)),std::sqrt(double(pe/pn)),std::sqrt(double(ge/gn))};
    std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_ACCURACY metric="<<(rz?"RZ":"XY")<<" temporal="<<temporal<<(mixed?" mixed=1":"")<<" n="<<n<<" cells="<<mesh.cellIndexNum()<<" dt="<<dt<<" steps="<<steps<<" time="<<flow.state().time
        <<" velocity="<<result.velocity<<" pressure="<<result.pressure<<" gradient="<<result.gradient<<" equation="<<equation<<" divergence="<<divergence<<" energy="<<energy<<" iterations="<<iterations<<std::endl;
    require(equation<=1e-10&&divergence<=1e-10&&energy<2e-8,"manufactured complete trajectory failed equation or work balance");return result;
}
}
void hydroMomentNativeAccuracyContract(){
    const char* family=std::getenv("QUANTUM_NATIVE_ACCURACY_FAMILY");const bool temporal=family&&std::string(family)=="temporal";
    const bool mixed=family&&std::string(family)=="mixed";
    require(!family||temporal||mixed||std::string(family)=="spatial","unknown complete accuracy family");
    const char* metric=std::getenv("QUANTUM_NATIVE_ACCURACY_METRIC");require(!metric||std::string(metric)=="XY"||std::string(metric)=="RZ","unknown complete accuracy metric");
    bool pass=true;
    for(bool rz:{false,true}){if(metric&&std::string(metric)!=(rz?"RZ":"XY"))continue;
        std::vector<NativeAccuracyResult> results;for(int level=0;level<3;level++)results.push_back(nativeAccuracy(temporal?4:2*(1<<level),rz,temporal,temporal?.02/(1<<level):.01,mixed));
        for(int level=1;level<3;level++){
            const auto& a=results[level-1];const auto& b=results[level];const double u=std::log2(a.velocity/b.velocity),p=std::log2(a.pressure/b.pressure),g=std::log2(a.gradient/b.gradient);
            std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_ACCURACY_ORDER metric="<<(rz?"RZ":"XY")<<" temporal="<<temporal<<(mixed?" mixed=1":"")<<" level="<<level<<" velocity="<<u<<" pressure="<<p<<" gradient="<<g<<std::endl;
            pass=pass&&u>1.7&&p>1.7&&g>1.7;
        }
    }
    require(pass,"complete physical velocity/pressure did not reach the declared second-order accuracy floor");
}
void hydroMomentNativeRemeshContract(){
    constexpr double dt=.002;
    const char* metric=std::getenv("QUANTUM_NATIVE_REMESH_METRIC");
    require(!metric||std::string(metric)=="XY"||std::string(metric)=="RZ","unknown complete remesh metric selection");
    for(bool rz:{false,true}){
        if(metric&&std::string(metric)!=(rz?"RZ":"XY"))continue;
        auto dictionary=input(4);dictionary.minLevel_neg=dictionary.minLevel_pos={0};dictionary.maxLevel_local={1};
        auto mesh=std::make_unique<QuadtreeMesh>(dictionary);FaceSegmentGraph graph;graph.build(*mesh);
        auto domain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*mesh,graph,rz));auto space=Space::build(*mesh,*domain);
        auto interface=InterfaceParametricMotion::initialize(*mesh,*domain,rz?InterfacePrescribedData::parametricSupplyCap(.5,.3,8):InterfacePrescribedData::parametricEllipse({.48,.52},{.18,.14},16));
        std::vector<HydroBoundaryPatch> ports{{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,.5,-.1},{"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
        ports[0].incomingPhase=rz?-1:1;ports[1].incomingPhase=1;ports[1].incomingVelocity=Vector2D{0,0};
        const Eigen::VectorXd zero=Eigen::VectorXd::Zero(space.mass().rows());
        auto flow=Flow::createTwoPhase(*mesh,*domain,std::move(space),std::move(interface),{1,.1,.01,.001},zero,.1,ports);
        Flow::Options options;if(rz)options.interfaceProjection.materialSegments=8;
        const auto liquid=[&]{const auto& inventory=flow.interfaceState().parametricGeometry().inventory();return std::accumulate(inventory.begin(),inventory.end(),0.);};
        const double volume0=liquid(),supply=rz?std::acos(-1.)*.025:0;
        const auto advance=[&]{const auto& state=flow.state();const HydroStepTime time{state.acceptedSteps+1,state.time,dt,HydroTimePolicy::Fixed};
            const auto load=flow.makeLoad(*mesh,*domain,Eigen::VectorXd::Zero(flow.space().mass().rows()),time.begin+.5*dt);
            auto candidate=flow.prepareNavierStokesStep(*mesh,*domain,time,load,options);const auto r=candidate.report();flow.accept(*mesh,*domain,std::move(candidate));
            const auto& next=flow.state();const double boundary=(flow.boundary().normalConstraints()*next.velocity-flow.boundary().normalConstraintValues()).norm();
            const double volume=std::abs(liquid()-volume0-supply*next.time),energyScale=r.kineticBefore+r.kineticAfter+std::abs(r.capillaryWork)+std::abs(r.boundaryReactionWork)+std::abs(r.boundaryTractionWork)+std::abs(r.viscousDissipation);
            const double energy=std::abs(r.energyDefect)/energyScale,spatial=r.materialMetricWork-r.boundaryKineticFlux-r.boundaryInflowDissipation-r.geometricMassWork;
            std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_REMESH_STEP metric="<<(rz?"RZ":"XY")<<" cells="<<mesh->cellIndexNum()<<" step="<<next.acceptedSteps<<" time="<<next.time
                <<" equation="<<r.equationResidual<<" divergence="<<r.divergenceResidual<<" iterations="<<r.nonlinearIterations<<" energy="<<energy
                <<" kinetic="<<r.kineticAfter<<" boundary="<<boundary<<" spatial_work="<<spatial<<" surface_work="<<r.surfaceWorkDefect<<" volume="<<volume<<std::endl;
            require(r.equationResidual<=1e-10&&boundary<2e-10&&volume<2e-12&&energy<2e-8,"complete remesh continuation lost physical step or original inlet constraints");
        };
        advance();const auto sourceState=flow.state();const auto sourceGeometry=flow.interfaceState().parametricGeometry().boundaryHandle();
        const auto sourceInventory=flow.interfaceState().parametricGeometry().inventory();const double correction=flow.initialBoundaryCorrection();
        const int boundaryRows=int(flow.boundary().normalConstraints().rows());
        auto target=mesh->clone();CellField<double> gradient(*target,0.),levels(*target,0.);
        for(int cell=0;cell<target->cellIndexNum();cell++)if(target->cells(target->cellIndexToLabel[cell]).centerCoord().x<.5)levels.value[cell]=1;
        std::vector<VertField<double>> scalars;std::vector<VertField<double>*> vertices;std::vector<HEdgeField<double>*> horizontal;
        std::vector<VEdgeField<double>*> vertical;std::vector<CellField<double>*> fields{&levels};target->adapt(scalars,gradient,vertices,horizontal,vertical,fields,&levels);
        require(target->cellIndexNum()==40,"complete remesh control did not execute the actual refiner");FaceSegmentGraph targetGraph;targetGraph.build(*target);
        auto targetDomain=std::make_unique<NativeLeafDomain>(NativeLeafDomain::build(*target,targetGraph,rz));
        auto candidate=flow.prepareMeshTransfer(*mesh,*domain,*target,*targetDomain,options);const auto report=candidate.report();
        require(flow.state().physicalVelocity->coordinates()==sourceState.physicalVelocity->coordinates()&&flow.state().time==sourceState.time
            &&flow.interfaceState().parametricGeometry().boundaryHandle()==sourceGeometry&&flow.interfaceState().parametricGeometry().inventory()==sourceInventory,"complete remesh preparation mutated source");
        bool foreign=false;try{flow.acceptMeshTransfer(*mesh,*domain,*target,*domain,std::move(candidate));}catch(const std::exception&){foreign=true;}
        require(foreign,"complete remesh accepted a foreign target identity");
        flow.acceptMeshTransfer(*mesh,*domain,*target,*targetDomain,std::move(candidate));
        auto oldMesh=std::move(mesh);auto oldDomain=std::move(domain);mesh=std::move(target);domain=std::move(targetDomain);
        const auto& next=flow.state();require(next.time==sourceState.time&&next.previousTime==sourceState.previousTime&&next.stageTime==sourceState.stageTime
            &&next.acceptedSteps==sourceState.acceptedSteps&&next.remappedStage&&flow.interfaceState().parametricGeometry().boundaryHandle()==sourceGeometry
            &&flow.initialBoundaryCorrection()==correction&&flow.boundary().normalConstraints().rows()==boundaryRows,"complete remesh changed clocks, geometry, initialization or boundary data count");
        double valueError=0,jetError=0,pressureError=0,childMeanDifference=0;
        const Spline::NativeVelocity* beforeFields[]={&*sourceState.physicalVelocity,&*sourceState.previousPhysicalVelocity,&*sourceState.midpointPhysicalVelocity};
        const Spline::NativeVelocity* afterFields[]={&*next.physicalVelocity,&*next.previousPhysicalVelocity,&*next.midpointPhysicalVelocity};
        for(int layer=0;layer<3;layer++){
            for(int cell=0;cell<mesh->cellIndexNum();cell++){const auto& leaf=mesh->cells(mesh->cellIndexToLabel[cell]);const auto center=leaf.centerCoord(),h=leaf.size();
                for(double x:{-.39,0.,.37})for(double y:{-.41,.05,.42}){const auto a=beforeFields[layer]->firstJet({center.x+x*h.x,center.y+y*h.y}),b=afterFields[layer]->firstJet({center.x+x*h.x,center.y+y*h.y});
                    valueError=std::max(valueError,std::hypot(a.value.x-b.value.x,a.value.y-b.value.y));jetError=std::max({jetError,std::hypot(a.dx.x-b.dx.x,a.dx.y-b.dx.y),std::hypot(a.dy.x-b.dy.x,a.dy.y-b.dy.y)});}
                const int parent=beforeFields[layer]->firstJet(center).nativeCell;const auto& old=oldMesh->cells(oldMesh->cellIndexToLabel[parent]);const auto c=old.centerCoord(),s=old.size();
                pressureError=std::max(pressureError,std::abs(next.pressure[3*cell]-sourceState.pressure.segment<3>(3*parent).dot(Eigen::Vector3d(1,(center.x-c.x)/s.x,(center.y-c.y)/s.y))));
            }
            const auto& r=report.levels.at(layer);
            std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_REMESH_LEVEL metric="<<(rz?"RZ":"XY")<<" layer="<<layer<<" time="<<r.time<<" momentum="<<r.momentumResidual
                <<" divergence="<<r.divergenceResidual<<" boundary="<<r.constraintResidual<<" kinetic_before="<<r.kineticBefore<<" kinetic_after="<<r.kineticAfter<<" exact="<<r.exactFieldRestriction<<std::endl;
            require(r.exactFieldRestriction&&!r.energyLimited&&r.momentumResidual<=3.2e-11,"complete remesh replaced exact restriction by a moment correction");
        }
        for(const auto face:domain->faces(*mesh))if(face.outwardSign&&face.area>0&&flow.boundary().faces()[face.nativeRow].kind==Quantum::Hydrodynamics::HydroCompatibleMomentBoundary::Kind::NormalVelocity){
            const int row=flow.space().faces()[face.nativeRow].first;childMeanDifference=std::max(childMeanDifference,std::abs(next.velocity[row]-flow.boundary().values()[row]));}
        bool consumed=false;try{flow.acceptMeshTransfer(*mesh,*domain,*mesh,*domain,std::move(candidate));}catch(const std::invalid_argument&){consumed=true;}
        std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_REMESH metric="<<(rz?"RZ":"XY")<<" source_cells=16 target_cells=40 value="<<valueError<<" jet="<<jetError<<" pressure="<<pressureError
            <<" reference_boundary_rows="<<boundaryRows<<" child_nominal_difference="<<childMeanDifference<<" consumed="<<consumed<<" geometry_identical=1 clock_identical=1"<<std::endl;
        require(valueError<2e-11&&jetError<2e-9&&pressureError<2e-11&&consumed&&childMeanDifference>1e-5,"complete remesh lost physical field/history or replaced inherited inlet conditions");
        advance();
    }
}
void hydroMomentNativeDensityTransportContract(){
    using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;
    constexpr double speed=.25,dt=.125,rhoLiquid=1000,rhoGas=1;constexpr int steps=8;
    const char* metric=std::getenv("QUANTUM_NATIVE_DENSITY_METRIC");require(!metric||std::string(metric)=="XY"||std::string(metric)=="RZ","unknown density-transport metric");
    for(bool rz:{false,true}){if(metric&&std::string(metric)!=(rz?"RZ":"XY"))continue;
    for(bool viscous:{false,true}){
        QuadtreeMesh mesh(input(4));DictionaryMesh::RefineBox box;box.min={0,0};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});
        require(mesh.cellIndexNum()==40,"density-transport grid is not the declared 2:1 grid");FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);auto space=Space::build(mesh,domain);
        auto interface=InterfaceParametricMotion::initialize(mesh,domain,rz?InterfacePrescribedData::parametricSpheroid(.30,{.18,.12},16):InterfacePrescribedData::parametricEllipse({.46,.30},{.16,.12},16));
        std::vector<HydroBoundaryPatch> ports{{"inlet",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,1,-speed},
            {"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,0},
            {"side_east",HydroBoundarySide::East,HydroBoundaryType::NormalVelocity,0,1,0}};
        if(!rz)ports.push_back({"side_west",HydroBoundarySide::West,HydroBoundaryType::NormalVelocity,0,1,0});
        for(auto& p:ports){p.incomingPhase=1;p.incomingVelocity=Vector2D{0,speed};}
        const auto initial=space.interpolateVelocity(mesh,domain,[&](Vector2D){return Eigen::Vector2d{0,speed};});
        const double muLiquid=viscous?1e-4:0,muGas=viscous?1e-6:0;
        auto flow=Flow::createTwoPhase(mesh,domain,std::move(space),std::move(interface),{rhoLiquid,rhoGas,muLiquid,muGas},initial,0,ports);
        const auto curve0=flow.interfaceState().parametricGeometry().geometry();const double energy0=flow.kineticEnergy();
        const auto liquid=[&](){const auto& v=flow.interfaceState().parametricGeometry().inventory();return std::accumulate(v.begin(),v.end(),0.);};const double volume0=liquid();
        const auto& rule=Quantum::Geometry::Parametric::gaussEight();
        double maximumVelocity=0,maximumGradient=0,maximumPressure=0,maximumLocalMomentum=0,maximumMass=0,maximumVolume=0,maximumShape=0,maximumEnergy=0,maximumFraction=0;
        double cumulativeRaw=0,cumulativeTemporal=0,cumulativeGeometric=0;
        std::ofstream curves,history;if(const char* prefix=std::getenv("QUANTUM_NATIVE_DENSITY_OUTPUT")){
            const std::string file=std::string(prefix)+(rz?"_rz":"_xy")+(viscous?"_low_mu":"_zero_mu");curves.open(file+"_curves.csv");history.open(file+"_history.csv");
            curves<<std::setprecision(17)<<"step,time,parameter,x,y,exact_x,exact_y\n";
            history<<std::setprecision(17)<<"step,time,velocity,gradient,pressure,local_momentum,mass,volume,shape,kinetic,raw_energy,temporal,geometric,viscous,fraction,equation\n";
        }
        const auto save=[&](int step){if(curves)for(int k=0;k<=512;k++){const double s=double(k)/512;const auto a=flow.interfaceState().parametricGeometry().geometry().sample(s),b=curve0.sample(s);
            curves<<step<<','<<flow.state().time<<','<<s<<','<<a.x<<','<<a.y<<','<<b.x<<','<<b.y+speed*flow.state().time<<'\n';}if(curves)curves.flush();};save(0);
        for(int step=1;step<=steps;step++){
            const double begin=flow.state().time;const auto before=flow.physicalCellBalance();const auto load=flow.makeLoad(mesh,domain,Eigen::VectorXd::Zero(initial.size()),begin+.5*dt);
            auto candidate=flow.prepareNavierStokesStep(mesh,domain,{step,begin,dt,HydroTimePolicy::Fixed},load);const auto r=candidate.report();
            const auto mass=Transport::phaseMassFlux(mesh,domain,flow.space(),flow.interfaceState(),candidate.interfaceState(),{step,begin,begin+dt,dt},rhoLiquid,rhoGas);
            flow.accept(mesh,domain,std::move(candidate));const auto after=flow.physicalCellBalance();
            std::vector<double> net(mesh.cellIndexNum(),0.);for(const auto f:domain.faces(mesh)){
                if(f.negativeCell>=0)net[f.negativeCell]+=mass.face[f.nativeRow][0];if(f.positiveCell>=0)net[f.positiveCell]-=mass.face[f.nativeRow][0];}
            double velocity=0,gradient=0,pressure=0,localMomentum=0;
            for(int c=0;c<mesh.cellIndexNum();c++){
                Eigen::Vector2d residual=after.momentum.col(c)-before.momentum.col(c);residual[1]+=speed*net[c];
                localMomentum=std::max(localMomentum,residual.norm()/(speed*(mass.massBefore[c]+mass.massAfter[c])));
                const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),h=cell.size();
                for(int i=0;i<8;i++)for(int j=0;j<8;j++){
                    const Vector2D at=center+Vector2D{.5*h.x*rule.x[i],.5*h.y*rule.x[j]};const auto v=flow.state().physicalVelocity->firstJet(at);
                    velocity=std::max(velocity,std::hypot(v.value.x,v.value.y-speed)/speed);
                    gradient=std::max(gradient,std::min(h.x,h.y)*std::sqrt(v.dx.x*v.dx.x+v.dx.y*v.dx.y+v.dy.x*v.dy.x+v.dy.y*v.dy.y)/speed);
                    const Eigen::Vector3d basis{1,(at.x-center.x)/h.x,(at.y-center.y)/h.y};pressure=std::max(pressure,std::abs(flow.state().pressure.segment<3>(3*c).dot(basis)));
                }
            }
            double shape=0;for(int k=0;k<=512;k++){const double s=double(k)/512;const auto a=flow.interfaceState().parametricGeometry().geometry().sample(s),b=curve0.sample(s);shape=std::max(shape,std::hypot(a.x-b.x,a.y-b.y-speed*flow.state().time));}
            const double volume=std::abs(liquid()-volume0),fraction=flow.interfaceState().parametricGeometry().report().maximumInventoryFractionError;
            // Raw physical balance: do not subtract the diagnosed geometric
            // or temporal mass work to manufacture a smaller energy error.
            const double raw=r.kineticAfter-r.kineticBefore+r.viscousDissipation+r.boundaryKineticFlux+r.boundaryInflowDissipation
                -r.externalWork-r.pressureWork-r.boundaryTractionWork-r.boundaryReactionWork-r.capillaryWork;
            cumulativeRaw+=raw;cumulativeTemporal+=r.temporalMassWork;cumulativeGeometric+=r.geometricMassWork;
            maximumVelocity=std::max(maximumVelocity,velocity);maximumGradient=std::max(maximumGradient,gradient);maximumPressure=std::max(maximumPressure,pressure);
            maximumLocalMomentum=std::max(maximumLocalMomentum,localMomentum);maximumMass=std::max(maximumMass,mass.maximumRelativeClosure);maximumVolume=std::max(maximumVolume,volume);
            maximumShape=std::max(maximumShape,shape);maximumEnergy=std::max(maximumEnergy,std::abs(r.kineticAfter-energy0)/energy0);maximumFraction=std::max(maximumFraction,fraction);
            std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_DENSITY_STEP metric="<<(rz?"RZ":"XY")<<" viscous="<<viscous<<" step="<<step<<" time="<<flow.state().time
                <<" velocity="<<velocity<<" gradient="<<gradient<<" pressure="<<pressure<<" local_momentum="<<localMomentum<<" mass="<<mass.maximumRelativeClosure<<" volume="<<volume<<" shape="<<shape
                <<" kinetic="<<r.kineticAfter<<" raw_energy="<<raw<<" cumulative_raw="<<cumulativeRaw/energy0<<" temporal="<<r.temporalMassWork<<" geometric="<<r.geometricMassWork
                <<" viscous_work="<<r.viscousDissipation<<" fraction="<<fraction<<" equation="<<r.equationResidual<<" iterations="<<r.nonlinearIterations<<std::endl;
            if(history){history<<step<<','<<flow.state().time<<','<<velocity<<','<<gradient<<','<<pressure<<','<<localMomentum<<','<<mass.maximumRelativeClosure<<','<<volume<<','<<shape<<','<<r.kineticAfter<<','<<raw<<','<<r.temporalMassWork<<','<<r.geometricMassWork<<','<<r.viscousDissipation<<','<<fraction<<','<<r.equationResidual<<'\n';history.flush();}
            if(step==4||step==steps)save(step);
            require(std::abs(after.kinetic.sum()-flow.kineticEnergy())/energy0<2e-12,"owned physical cell energies disagree with complete kinetic energy");
            require(viscous||r.viscousDissipation==0,"zero viscosity acquired artificial strain dissipation");
            require(velocity<2e-8&&gradient<2e-7&&pressure/(rhoLiquid*speed*speed)<2e-8&&localMomentum<2e-8&&mass.maximumRelativeClosure<2e-11
                &&volume<2e-12&&shape<2e-8&&fraction<3.1e-12&&maximumEnergy<2e-8&&std::abs(cumulativeRaw)/energy0<2e-8,
                "high-density actual transport did not preserve its exact uniform state or uncorrected physical budgets");
        }
        std::cout<<std::setprecision(17)<<"MOMENT_NATIVE_DENSITY_COMPLETE metric="<<(rz?"RZ":"XY")<<" viscous="<<viscous<<" rho_l="<<rhoLiquid<<" rho_g="<<rhoGas<<" mu_l="<<muLiquid<<" mu_g="<<muGas
            <<" cells="<<mesh.cellIndexNum()<<" steps="<<flow.state().acceptedSteps<<" time="<<flow.state().time<<" displacement="<<speed*flow.state().time
            <<" velocity="<<maximumVelocity<<" gradient="<<maximumGradient<<" pressure="<<maximumPressure<<" local_momentum="<<maximumLocalMomentum<<" mass="<<maximumMass
            <<" volume="<<maximumVolume<<" shape="<<maximumShape<<" energy_drift="<<maximumEnergy<<" raw_energy="<<cumulativeRaw/energy0<<" temporal="<<cumulativeTemporal/energy0<<" geometric="<<cumulativeGeometric/energy0<<" fraction="<<maximumFraction<<std::endl;
    }}
}
void hydroMomentNeckReadinessContract(){
    namespace IP=InterfaceParametricProjection;
    {
        // Exact vertical seam: z(r) is singular although the curve and its
        // curvature are regular. Both seam orientations must close directly.
        namespace P=Quantum::Geometry::Parametric;constexpr double h=.0625,r=.025,neck=.0005,a=r-neck;
        P::Curve c;c.boundary=P::CurveBoundary::axisMeridian();
        // Retained curvature is deliberately half the removed-side curvature
        // at the C1 knot. The cap must match the retained one-sided value.
        P::Arc lower{0,.25,{{{0,.5-1.5*h},{2*r-2*a/3,0},{0,.5*h},{2*a-2*r,0},{r-4*a/3,0}}}};
        P::Arc middle{.25,.75,{{{r,.5-h},{0,2*h},{-16*a,0},{32*a,0},{-16*a,0}}}};
        auto upper=P::restrictedArc(lower,.25,0);upper.a=.75;upper.b=1;for(auto& p:upper.p)p.y=-p.y;upper.p[0].y+=1;
        c.arcs={lower,middle,upper};IP::Grid grid;grid.nx=grid.ny=16;grid.axisymmetric=true;
        auto original=IP::initialize(std::move(c),grid);const auto boundary=original.boundaryHandle();
        auto assessed=IP::assessNeckSplit(original,{1,.5,.1,1e-10});const auto& g=assessed.report;
        require(assessed.readiness==IP::NeckReadiness::Ready&&assessed.split&&g.lowerPole<g.upperPole,"regular vertical seam did not admit its parametric cap");
        require(g.seamPosition<2e-13&&g.seamTangent<2e-11&&g.seamCurvature<2e-9&&g.maximumInventoryFractionError<1e-10
            &&original.boundaryHandle()==boundary&&assessed.split->candidate.inventory()==original.inventory(),"parametric cap lost seam geometry or independent inventory");
        const double exact=2*std::acos(-1.)*h*(neck*neck+2*neck*a*(2./3-1./5)+a*a*(4./5-4./7+1./9));
        require(std::abs(g.modelRegionVolume-exact)<2e-16&&std::abs(g.seamRadialTangent[0])<1e-13&&std::abs(g.seamRadialTangent[1])<1e-13
            &&std::abs(g.seamSignedCurvature[0]-4*a/(h*h))<2e-10&&std::abs(g.seamSignedCurvature[1]-4*a/(h*h))<2e-10,"cap did not consume the retained one-sided seam curvature");
        std::cout<<std::setprecision(17)<<"PARAMETRIC_NECK_VERTICAL lower_pole="<<g.lowerPole<<" upper_pole="<<g.upperPole<<" volume="<<g.modelRegionVolume<<" seam_curvature_error="<<g.seamCurvature<<" fraction="<<g.maximumInventoryFractionError<<std::endl;
    }
    QuadtreeMesh mesh(input(16));DictionaryMesh::RefineBox box;box.min={0,.375};box.max={.125,.625};box.level=1;mesh.applyRefineBoxes({box});
    FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
    const auto initial=[&](double radius){return InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricNeckedSpheroid(.5,{.35,.3},radius/.35,32));};
    // A nearly cylindrical neck has too much local liquid to fit two regular
    // axis caps in this support. Pending must remain a physical, valid state.
    auto interface=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricNeckedSpheroid(.5,{.005,.3},.4,32));
    const auto geometry=interface.parametricGeometry().boundaryHandle();const auto inventory=interface.parametricGeometry().inventory();
    const auto sites=IP::inspectNecks(interface.parametricGeometry(),.1);require(sites.size()==1&&sites[0].localScale==.03125,"readiness control did not reproduce the fine neck support");
    const auto cap=InterfaceParametricMotion::assessNeckSplit(mesh,domain,interface,sites[0].options);
    require(cap.readiness==IP::NeckReadiness::CapPending&&!cap.state,"infeasible volume/C2 cap was published");
    require(interface.parametricGeometry().boundaryHandle()==geometry&&interface.parametricGeometry().inventory()==inventory&&interface.parametricNeckEvents().empty(),"pending cap changed accepted geometry or inventory");
    bool rejected=false;try{auto denied=InterfaceParametricMotion::prepareNeckSplit(mesh,domain,interface,sites[0].options);}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"strict event preparation accepted an infeasible cap");
    for(double radius:{.0005,.004}){
        auto state=initial(radius);const auto observed=IP::inspectNecks(state.parametricGeometry(),.1);require(observed.size()==1,"control neck was lost");
        auto assessment=InterfaceParametricMotion::assessNeckSplit(mesh,domain,state,observed[0].options);
        require(assessment.readiness==(radius==.0005?IP::NeckReadiness::Ready:IP::NeckReadiness::ScalePending),"joint readiness control returned the wrong geometric state");
        require(bool(assessment.state)==(radius==.0005),"pending geometry owns an event candidate");
    }
    auto flow=Flow::createTwoPhase(mesh,domain,space,std::move(interface),{5,1,.02,.01},Eigen::VectorXd::Zero(space.mass().rows()),.1);
    auto assessment=flow.assessNeckSplit(mesh,domain,sites[0].options);
    require(assessment.readiness==Flow::NeckReadiness::CapPending&&!assessment.candidate&&flow.state().time==0&&flow.state().velocity.isZero(),"coupled pending event changed accepted fluid state");
    constexpr double dt=1e-5;Flow::Options options;options.interfaceProjection.materialSegments=32;
    auto load=flow.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),.5*dt);
    auto next=flow.prepareNavierStokesStep(mesh,domain,{1,0,dt},load,options);const auto report=next.report();flow.accept(mesh,domain,std::move(next));
    const auto& evolved=flow.interfaceState().parametricGeometry();
    const double volume=std::accumulate(evolved.inventory().begin(),evolved.inventory().end(),0.)-std::accumulate(inventory.begin(),inventory.end(),0.);
    require(flow.state().time==dt&&flow.state().acceptedSteps==1&&flow.state().velocity.norm()>0&&evolved.componentCount()==1&&flow.interfaceState().parametricNeckEvents().empty()
        &&std::abs(volume)<2e-13&&report.equationResidual<1e-10,"finite pending neck did not continue through the same physical NS step");
    std::cout<<std::setprecision(17)<<"MOMENT_NECK_READINESS cells="<<mesh.cellIndexNum()<<" h="<<sites[0].localScale<<" radius="<<sites[0].radius<<" lower_pole="<<cap.report.lowerPole<<" upper_pole="<<cap.report.upperPole
        <<" cap_pending=1 scale_pending=1 ready_control=1 physical_steps=1 time="<<flow.state().time<<" volume_error="<<volume<<" equation="<<report.equationResidual<<" velocity_norm="<<flow.state().velocity.norm()<<std::endl;
}
void hydroMomentNeckContract(){
    namespace P=Quantum::Geometry::Parametric;namespace IP=InterfaceParametricProjection;
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    std::ofstream history,curves;if(const char* prefix=std::getenv("QUANTUM_MOMENT_NECK_OUTPUT")){
        history.open(std::string(prefix)+"_history.csv");curves.open(std::string(prefix)+"_curves.csv");require(history&&curves,"cannot open coupled neck outputs");
        history<<std::setprecision(17)<<"supplied,stage,step,time,components,liquid,kinetic,surface\n";
        curves<<std::setprecision(17)<<"supplied,stage,step,time,component,s,r,z,curvature_model\n";
    }
    for(bool supplied:{false,true}){
        constexpr double dt=.0001,sigma=.1;QuadtreeMesh mesh(input(8));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
        auto curve=InterfacePrescribedData::parametricNeckedSpheroid(supplied?.375:.5,{.35,supplied?.5:.3},(supplied?.002:.008)/.35,64);
        std::vector<HydroBoundaryPatch> ports;
        if(supplied){
            double begin=-1;for(const auto& a:curve.arcs){auto p=a.coordinate(1);for(double t:P::roots(p,0,1))begin=a.a+(a.b-a.a)*t;}
            require(begin>0&&begin<1,"supplied neck initial seam missing");P::Curve cropped;cropped.boundary=P::CurveBoundary::suppliedMeridian({.125,0},{0,0});
            for(const auto& a:curve.arcs)if(a.b>begin){auto part=P::restrictedArc(a,std::max(begin,a.a),a.b);part.a=(part.a-begin)/(1-begin);part.b=(part.b-begin)/(1-begin);cropped.arcs.push_back(part);}
            const double radial=.125/cropped.arcs.front().p[0].x;for(auto& a:cropped.arcs)for(auto& p:a.p)p.x*=radial;cropped.arcs.front().p[0]={.125,0};cropped.arcs.front().a=0;cropped.arcs.back().b=1;curve=std::move(cropped);
            ports={{"supply",HydroBoundarySide::South,HydroBoundaryType::NormalVelocity,0,.125,-.04},{"outlet",HydroBoundarySide::North,HydroBoundaryType::PressureOutlet,0,1,.2}};
            ports[0].incomingPhase=-1;ports[1].incomingPhase=1;ports[1].incomingVelocity=Vector2D{0,0};
        }
        const auto initial=space.interpolateVelocity(mesh,domain,[&](Vector2D p)->Eigen::Vector2d{
            if(supplied)return Eigen::Vector2d::Zero();const double r=p.x,z=p.y,q=r*r,g=z*z*(1-z)*(1-z),gy=2*z*(1-z)*(1-2*z);
            return {.24*r*(1-q)*(1-q)*gy,.24*(-2+8*q-6*q*q)*g};});
        auto interface=InterfaceParametricMotion::initialize(mesh,domain,std::move(curve));
        auto flow=Flow::createTwoPhase(mesh,domain,space,std::move(interface),{5,1,.02,.01},initial,sigma,ports);
        const auto volume=[](const IP::Candidate& c,std::size_t i){const auto& v=c.componentInventory(i);return std::accumulate(v.begin(),v.end(),0.);};
        const double v0=volume(flow.interfaceState().parametricGeometry(),0),k0=flow.kineticEnergy();
        double viscous=0,capillary=0,temporal=0,metric=0,boundary=0,maxEquation=0,maxVolume=0,maxFraction=0,eventK=0;
        Flow::NeckReport event;const auto save=[&](const char* stage){const auto& state=flow.state();const auto& shape=flow.interfaceState().parametricGeometry();
            if(history){history<<supplied<<','<<stage<<','<<state.acceptedSteps<<','<<state.time<<','<<shape.componentCount()<<','<<std::accumulate(shape.inventory().begin(),shape.inventory().end(),0.)<<','<<flow.kineticEnergy()<<','<<sigma*Quantum::Hydrodynamics::HydroCompatibleMomentTraction::measure(*shape.boundaryHandle(),true).length<<'\n';history.flush();}
            if(curves)for(const auto& c:shape.boundaryHandle()->components())for(int k=0;k<=256;k++){const auto p=c.curve->sample(k/256.);curves<<supplied<<','<<stage<<','<<state.acceptedSteps<<','<<state.time<<','<<c.id<<','<<k/256.<<','<<p.x<<','<<p.y<<','<<flow.interfaceState().parametricModelCurvature(c.id,k/256.)<<'\n';}
        };save("initial");
        Flow::Options options;options.interfaceProjection.materialSegments=32;
        for(int step=1;step<=4;step++){
            const double begin=flow.state().time;const auto load=flow.makeLoad(mesh,domain,Eigen::VectorXd::Zero(initial.size()),begin+.5*dt);
            auto next=flow.prepareNavierStokesStep(mesh,domain,{step,begin,dt},load,options);const auto r=next.report();flow.accept(mesh,domain,std::move(next));
            viscous+=r.viscousDissipation;capillary+=r.capillaryWork;temporal+=r.temporalMassWork;metric+=r.materialMetricWork;boundary+=r.boundaryTractionWork+r.boundaryReactionWork;
            maxEquation=std::max(maxEquation,r.equationResidual);maxFraction=std::max(maxFraction,flow.interfaceState().parametricGeometry().report().maximumInventoryFractionError);
            const auto& shape=flow.interfaceState().parametricGeometry();const double rate=supplied?std::acos(-1.)*.125*.125*.04:0;
            maxVolume=std::max(maxVolume,std::abs(std::accumulate(shape.inventory().begin(),shape.inventory().end(),0.)-v0-rate*flow.state().time));
            if(step>1)for(int i=0;i<2;i++)maxVolume=std::max(maxVolume,std::abs(volume(shape,i)-event.geometry.childVolumes[i]-(i==0?rate*(step-1)*dt:0)));
            save("physical");
            if(step==1){
                const auto detected=IP::detectNecks(shape,.1);require(detected.size()==1,"actual neck detector did not find the one unresolved neck");
                const auto oldGeometry=shape.boundaryHandle();const auto oldInventory=shape.inventory();const auto oldState=flow.state();const auto oldMass=flow.kineticMass();
                auto denied=detected[0];denied.maximumRadiusCells=1e-5;bool rejected=false;
                try{auto fail=flow.prepareNeckSplit(mesh,domain,denied);}catch(const std::invalid_argument&){rejected=true;}
                require(rejected&&flow.interfaceState().parametricGeometry().boundaryHandle()==oldGeometry&&flow.state().velocity==oldState.velocity,"failed coupled event changed accepted state");
                auto candidate=flow.prepareNeckSplit(mesh,domain,detected[0]);event=candidate.report();
                require(candidate.state().previousVelocity==oldState.previousVelocity&&candidate.state().midpointVelocity==oldState.midpointVelocity&&candidate.state().pressure==oldState.pressure
                    &&candidate.state().stageTime==oldState.stageTime&&candidate.state().time==oldState.time&&candidate.state().acceptedSteps==oldState.acceptedSteps,"event rewrote historical time layers");
                require(oldState.continuousHistory&&!candidate.state().continuousHistory,"instantaneous neck event retained a smooth predictor history");
                require(candidate.interfaceState().parametricGeometry().inventory()==oldInventory&&flow.interfaceState().parametricGeometry().boundaryHandle()==oldGeometry,"coupled event changed inventory or committed before acceptance");
                flow.acceptNeckSplit(mesh,domain,std::move(candidate));eventK=event.kineticChange;
                bool consumed=false;try{flow.acceptNeckSplit(mesh,domain,std::move(candidate));}catch(const std::invalid_argument&){consumed=true;}require(consumed,"coupled event accepted twice");
                for(int i=0;i<2;i++){const auto range=event.geometry.capIntervals[i];require(flow.interfaceState().parametricModelCurvature(event.geometry.children[i],.5*(range[0]+range[1])),"event cap curvature lost model provenance");}
                require(event.momentum.momentumResidual<2e-10&&event.momentum.divergenceResidual<1e-11&&event.relativeMassMatrixChange>1e-10
                    &&event.totalEnergyChange<=1e-12&&event.momentum.kineticBefore>0,"coupled event failed momentum, energy or actual material change");
                save("event");std::cout<<std::setprecision(17)<<"MOMENT_NECK_EVENT supplied="<<supplied<<" time="<<flow.state().time<<" section="<<event.geometry.section<<" mass_change="<<event.relativeMassMatrixChange
                    <<" momentum="<<event.momentum.momentumResidual<<" divergence="<<event.momentum.divergenceResidual<<" kinetic_before="<<event.momentum.kineticBefore<<" kinetic_after="<<event.momentum.kineticAfter
                    <<" surface_before="<<event.surfaceBefore<<" surface_after="<<event.surfaceAfter<<" model_dissipation="<<event.modelDissipation<<" limited="<<event.momentum.energyLimited<<std::endl;
            }
        }
        const double k=flow.kineticEnergy(),balance=k-k0-eventK+viscous-capillary+metric-temporal-boundary;
        require(maxEquation<1e-10&&maxVolume<2e-13&&maxFraction<1e-10&&std::abs(balance)<2e-10&&flow.interfaceState().parametricGeometry().componentCount()==2,"coupled post-neck physical trajectory failed");
        std::cout<<"MOMENT_NECK_COMPLETE supplied="<<supplied<<" steps=4 post_event_steps=3 time="<<flow.state().time<<" volume_error="<<maxVolume<<" fraction="<<maxFraction<<" equation="<<maxEquation<<" kinetic_balance="<<balance<<std::endl;
    }
}
void hydroParametricNeckContract(){
    namespace P=Quantum::Geometry::Parametric;namespace IP=InterfaceParametricProjection;namespace IT=InterfaceParametricTransport;
    std::ofstream curves,history;
    if(const char* prefix=std::getenv("QUANTUM_PARAMETRIC_NECK_OUTPUT")){
        curves.open(std::string(prefix)+"_curves.csv");history.open(std::string(prefix)+"_history.csv");require(curves&&history,"cannot open neck outputs");
        curves<<std::setprecision(17)<<"n,stage,step,time,component,parameter,x,y\n";
        history<<std::setprecision(17)<<"n,stage,step,time,components,total,lower,upper,fraction\n";
    }
    for(int n:{8,16}){
        constexpr int segments=64,steps=14;constexpr double dt=.01,speed=.1;
        QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);
        P::Curve initial;initial.boundary=P::CurveBoundary::axisMeridian();const double pi=std::acos(-1.),radius=.35,amplitude=1-.002/radius;
        const auto position=[&](double s){const double q=std::sin(pi*s);return P::Point{s==0||s==1?0:radius*q*(1-amplitude*q*q),.498-.3*std::cos(pi*s)};};
        const auto derivative=[&](double s){const double q=std::sin(pi*s);return P::Point{pi*radius*std::cos(pi*s)*(1-3*amplitude*q*q),s==0||s==1?0:.3*pi*q};};
        for(int k=0;k<segments;k++){
            const double a=double(k)/segments,b=double(k+1)/segments,h=b-a;const auto p=position(a),delta=position(b)-p,d0=derivative(a)*h,d1=derivative(b)*h;
            P::Arc arc{a,b,{{p,d0,delta*3-d0*2-d1,delta*(-2)+d0+d1,{}}}};
            const auto bubble=(position(.5*(a+b))-arc.sampleFraction(.5))*16;arc.p[2]+=bubble;arc.p[3]+=bubble*(-2);arc.p[4]=bubble;initial.arcs.push_back(arc);
        }
        auto state=InterfaceParametricMotion::initialize(mesh,domain,std::move(initial));
        const InterfacePrescribedData::AxialColumn velocity(1,speed);const auto flux=velocity.normalFlux(mesh,domain);
        const IT::Flow map=[=](Vector2D p,double duration){const double q=1-p.x*p.x,u=speed*std::pow(q,5),dr=-10*speed*p.x*std::pow(q,4);
            return IT::MapJet{{p.x,p.y+duration*u},{1,duration*dr},{0,1},{0,u},{0,duration*u}};};
        const auto volume=[](const IP::Candidate& c,std::size_t i){const auto& v=c.componentInventory(i);return std::accumulate(v.begin(),v.end(),0.);};
        const double initialVolume=volume(state.parametricGeometry(),0);double maxDrift=0,maxFraction=0;
        const auto save=[&](const char* stage,bool shape){const auto& c=state.parametricGeometry();
            if(history){history<<n<<','<<stage<<','<<state.progress().completedSteps<<','<<state.progress().endTime<<','<<c.componentCount()<<','<<std::accumulate(c.inventory().begin(),c.inventory().end(),0.)<<','<<volume(c,0)<<','<<(c.componentCount()>1?volume(c,1):0)<<','<<c.report().maximumInventoryFractionError<<'\n';history.flush();}
            if(curves&&shape)for(std::size_t i=0;i<c.componentCount();i++)for(int k=0;k<=512;k++){const auto p=c.boundaryHandle()->components()[i].curve->sample(k/512.);
                curves<<n<<','<<stage<<','<<state.progress().completedSteps<<','<<state.progress().endTime<<','<<c.boundaryHandle()->components()[i].id<<','<<k/512.<<','<<p.x<<','<<p.y<<'\n';}
        };
        save("initial",true);IP::NeckSplitReport event;std::array<double,2> childVolume{};
        for(int step=1;step<=steps;step++){
            const InterfaceStepTime time{step,state.progress().endTime,step*dt,dt};
            auto prepared=InterfaceParametricMotion::prepare(mesh,domain,state,time,dt,map,flux,{});
            require(prepared.state.parametricFollows(state,time),"neck trajectory lost ordinary-step lineage");state=std::move(prepared.state);
            const auto& c=state.parametricGeometry();maxFraction=std::max(maxFraction,c.report().maximumInventoryFractionError);
            maxDrift=std::max(maxDrift,std::abs(std::accumulate(c.inventory().begin(),c.inventory().end(),0.)-initialVolume));
            if(step>2){require(c.componentCount()==2&&state.parametricNeckEvents().size()==1,"post-event state lost component ancestry");
                for(int i=0;i<2;i++){require(c.boundaryHandle()->components()[i].id==event.children[i],"post-event identity changed");maxDrift=std::max(maxDrift,std::abs(volume(c,i)-childVolume[i]));}}
            save("transport",step==2||step==steps);
            if(step==2){
                const auto before=c.boundaryHandle();const auto inventory=c.inventory();IP::NeckSplitOptions options;options.parent=1;options.section=.5;
                auto denied=options;denied.maximumRadiusCells=1e-4;bool rejected=false;
                try{auto fail=InterfaceParametricMotion::prepareNeckSplit(mesh,domain,state,denied);}catch(const std::invalid_argument&){rejected=true;}
                require(rejected&&state.parametricGeometry().boundaryHandle()==before&&state.parametricGeometry().inventory()==inventory,"failed neck model changed accepted state");
                auto candidate=InterfaceParametricMotion::prepareNeckSplit(mesh,domain,state,options);event=candidate.parametricNeckEvents().back().report;
                require(state.parametricGeometry().boundaryHandle()==before&&candidate.progress().completedSteps==2&&candidate.progress().endTime==time.end
                    &&candidate.parametricGeometry().inventory()==inventory,"zero-time neck event changed clock or cell inventory");
                bool noFlux=false;try{candidate.parametricFluxes();}catch(const std::logic_error&){noFlux=true;}require(noFlux,"zero-time neck event manufactured time-integrated flux");
                double remotePosition=0,remoteCurvature=0;
                for(int child=0;child<2;child++){
                    const auto& actual=*candidate.parametricGeometry().boundaryHandle()->components()[child].curve;
                    const auto parentRange=event.retainedParentIntervals[child];const double childBegin=child?event.capIntervals[1][1]:0,childEnd=child?1:event.capIntervals[0][0];
                    for(int k=0;k<=128;k++){const double f=k/128.,s=parentRange[0]+f*(parentRange[1]-parentRange[0]),t=childBegin+f*(childEnd-childBegin);
                        const auto a=before->single()->sample(s),b=actual.sample(t);remotePosition=std::max(remotePosition,P::norm(a-b));
                        const auto oldArc=P::restrictedArc(before->single()->arcs[before->single()->index(s)],before->single()->arcs[before->single()->index(s)].a,before->single()->arcs[before->single()->index(s)].b);
                        const auto& newArc=actual.arcs[actual.index(t)];
                        const auto qa=P::sampleSurfaceRZ(oldArc,(s-oldArc.a)/(oldArc.b-oldArc.a),s==0||s==1),qb=P::sampleSurfaceRZ(newArc,(t-newArc.a)/(newArc.b-newArc.a),t==0||t==1);
                        remoteCurvature=std::max(remoteCurvature,std::abs(qa.curvature-qb.curvature));
                    }
                }
                using Surface=Quantum::Hydrodynamics::HydroCompatibleMomentTraction;
                const double areaBefore=Surface::measure(*before,true).length,areaAfter=Surface::measure(*candidate.parametricGeometry().boundaryHandle(),true).length;
                require(remotePosition<2e-12&&remoteCurvature<2e-8&&std::isfinite(areaBefore)&&std::isfinite(areaAfter),"neck model changed geometry outside its declared region");
                state=std::move(candidate);for(int i=0;i<2;i++)childVolume[i]=volume(state.parametricGeometry(),i);
                bool retired=false;try{auto fail=InterfaceParametricMotion::prepareNeckSplit(mesh,domain,state,options);}catch(const std::invalid_argument&){retired=true;}require(retired,"neck event reused a retired parent");
                require(event.maximumCellBalance==0&&event.maximumInventoryFractionError<1e-10&&event.seamPosition<2e-11&&event.seamTangent<2e-8&&event.seamCurvature<2e-8
                    &&event.lowerPole<event.upperPole,"local cap conservation or seam geometry failed");
                save("split",true);
                std::cout<<std::setprecision(17)<<"PARAMETRIC_NECK_EVENT n="<<n<<" time="<<time.end<<" parent="<<event.parent<<" lower_id="<<event.children[0]<<" upper_id="<<event.children[1]
                    <<" radius="<<event.neckRadius<<" section_radius="<<event.sectionRadius<<" h="<<event.localScale<<" lower_pole="<<event.lowerPole<<" upper_pole="<<event.upperPole
                    <<" model_volume="<<event.modelRegionVolume<<" total_volume="<<event.parentVolume<<" lower_volume="<<event.childVolumes[0]<<" upper_volume="<<event.childVolumes[1]
                    <<" cell_balance="<<event.maximumCellBalance<<" fraction="<<event.maximumInventoryFractionError<<" seam_position="<<event.seamPosition<<" seam_tangent="<<event.seamTangent<<" seam_curvature="<<event.seamCurvature
                    <<" remote_position="<<remotePosition<<" remote_curvature="<<remoteCurvature<<" surface_before="<<areaBefore<<" surface_after="<<areaAfter<<" surface_change="<<areaAfter-areaBefore<<std::endl;
            }
        }
        require(state.progress().completedSteps==steps&&state.progress().endTime==steps*dt&&maxDrift<2e-13&&maxFraction<1e-10,"complete neck/post-event trajectory failed");
        QuadtreeMesh target(input(n==8?16:8));FaceSegmentGraph targetGraph;targetGraph.build(target);const auto targetDomain=NativeLeafDomain::build(target,targetGraph,true);
        const auto transferred=InterfaceParametricMotion::prepareMeshTransfer(mesh,domain,state,target,targetDomain,{});
        require(transferred.state.progress().endTime==state.progress().endTime&&transferred.state.parametricGeometry().boundaryHandle()==state.parametricGeometry().boundaryHandle()
            &&transferred.state.parametricNeckEvents().size()==1&&transferred.state.parametricNeckEvents()[0].time==.02
            &&transferred.state.parametricNeckEvents()[0].report.children==event.children&&transferred.state.parametricNeckEvents()[0].report.capIntervals==event.capIntervals,
            "post-event remap lost geometry, clock or material ancestry");
        std::cout<<"PARAMETRIC_NECK_TRANSFER n="<<n<<" target="<<(n==8?16:8)<<" fraction="<<transferred.report.maximumFractionError<<" balance="<<transferred.report.inventoryBalance<<" ancestry=1\n";
        std::cout<<"PARAMETRIC_NECK_COMPLETE n="<<n<<" steps="<<steps<<" time="<<state.progress().endTime<<" post_event_steps=12 volume_error="<<maxDrift<<" fraction="<<maxFraction
            <<" components="<<state.parametricGeometry().componentCount()<<" events="<<state.parametricNeckEvents().size()<<std::endl;
    }
}
void hydroMomentComponentsContract(){
    namespace P=Quantum::Geometry::Parametric;namespace IP=InterfaceParametricProjection;namespace IT=InterfaceParametricTransport;
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;using Traction=Quantum::Hydrodynamics::HydroCompatibleMomentTraction;
    std::ofstream curves,history;
    if(const char* prefix=std::getenv("QUANTUM_MOMENT_COMPONENTS_OUTPUT")){
        curves.open(std::string(prefix)+"_curves.csv");history.open(std::string(prefix)+"_history.csv");require(curves&&history,"cannot open component outputs");
        curves<<std::setprecision(17)<<"rz,step,time,component,parameter,x,y\n";
        history<<std::setprecision(17)<<"rz,step,time,volume_1,volume_2,fraction,equation,kinetic,surface,energy_defect,surface_work_defect\n";
    }
    for(bool rz:{false,true}){
        std::vector<P::Curve> initial;
        initial.push_back(rz?InterfacePrescribedData::parametricSpheroid(.42,{.18,.12},16):InterfacePrescribedData::parametricEllipse({.47,.42},{.14,.12},16));
        initial.push_back(rz?InterfacePrescribedData::parametricSpheroid(.67,{.12,.09},16):InterfacePrescribedData::parametricEllipse({.47,.67},{.10,.09},16));
        IP::Grid grid;grid.nx=grid.ny=8;grid.axisymmetric=rz;
        auto translated=IP::initialize(initial,grid);const auto original=translated;const auto volume=[](const IP::Candidate& c,std::size_t i){const auto& v=c.componentInventory(i);return std::accumulate(v.begin(),v.end(),0.);};
        bool rejected=false;try{IP::initialize(std::vector<P::Curve>{initial[0],initial[0]},grid);}catch(const std::exception&){rejected=true;}
        require(rejected&&original.componentCount()==2,"overlapping components accepted");
        IT::NativeFlux flux;flux.grid=grid;
        flux.integral=[=](bool horizontal,int,int,double lo,double hi){return horizontal?(rz?std::acos(-1.)*(hi*hi-lo*lo):hi-lo)*.5:0.;};
        flux.density=[=](bool horizontal,int,int,double p){return horizontal?.5*(rz?2*std::acos(-1.)*p:1.):0.;};
        flux.densityMagnitudeBound=[=](bool horizontal,int,int){return horizontal?.5*(rz?2*std::acos(-1.):1.):0.;};
        const IT::Flow map=[](Vector2D p,double t){return IT::MapJet{{p.x,p.y+.5*t},{1,0},{0,1},{0,.5},{0,.5*t}};};
        IP::Options projection;double translationDrift=0,translationFraction=0;
        for(int step=0;step<6;step++){
            auto next=IT::prepare(translated,step<3?.04:-.04,map,flux,projection);translated=std::move(next.candidate);
            for(std::size_t i=0;i<2;i++){require(translated.boundaryHandle()->components()[i].id==i+1,"translation changed component identity");translationDrift=std::max(translationDrift,std::abs(volume(translated,i)-volume(original,i)));}
            translationFraction=std::max(translationFraction,translated.report().maximumInventoryFractionError);
        }
        double returned=0;for(std::size_t i=0;i<2;i++)for(int k=0;k<128;k++)returned=std::max(returned,P::norm(translated.boundaryHandle()->components()[i].curve->sample(k/128.)-initial[i].sample(k/128.)));
        require(translationDrift<2e-13&&translationFraction<3.1e-12&&returned<2e-8,"multicomponent translation/return exceeded original budgets");
        std::cout<<std::setprecision(17)<<"MOMENT_COMPONENT_TRANSLATION rz="<<rz<<" steps=6 volume_error="<<translationDrift<<" fraction="<<translationFraction<<" return_position="<<returned<<" overlap_rejected="<<rejected<<std::endl;
        QuadtreeMesh mesh(input(4));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto space=Space::build(mesh,domain);
        auto state=InterfaceParametricMotion::initialize(mesh,domain,initial);const auto boundary=state.parametricGeometry().boundaryHandle();
        const double v0[2]={volume(state.parametricGeometry(),0),volume(state.parametricGeometry(),1)};
        const Eigen::VectorXd zero=Eigen::VectorXd::Zero(space.mass().rows());const auto spline=Spline::build(mesh,domain,space,zero,0);
        const auto material=Material::buildFrozen(mesh,domain,space,boundary,0,{5,1,.02,.01});
        const auto pressure=Traction::build(mesh,domain,space,spline,material,[](const P::SurfaceSample& q){return Traction::Components{q.point.y<.56?2.:3.,0};});
        Eigen::VectorXd expected=Eigen::VectorXd::Zero(space.divergence().rows());double sum=0,capacity=0;int shared=0;
        for(int c=0;c<mesh.cellIndexNum();c++){const auto& m=material.cells()[c];const Eigen::Vector3d load=2*m.componentLiquid[0].col(0)+3*m.componentLiquid[1].col(0);
            expected.segment<3>(3*c)=m.volume.ldlt().solve(load);sum+=load[0];capacity+=m.volume(0,0);if(m.componentLiquid[0](0,0)>1e-10&&m.componentLiquid[1](0,0)>1e-10)shared++;}
        for(int c=0;c<mesh.cellIndexNum();c++)expected[3*c]-=sum/capacity;
        const double pressureError=(expected-pressure.pressurePotential).norm(),variation=pressure.variation.norm();
        require(shared>0&&pressureError<2e-12&&variation<2e-11,"component constant stresses do not occupy their distinct compatible pressure ranges");
        auto flow=Flow::createTwoPhase(mesh,domain,space,std::move(state),{5,1,.02,.01},zero,.1);const Space::Matrix mass0=flow.kineticMass();
        Flow::Options options;options.interfaceProjection.materialSegments=16;double drift=0,fraction=0,equation=0,energy=0,surfaceError=0;
        const auto save=[&](int step){if(curves)for(const auto& c:flow.interfaceState().parametricGeometry().boundaryHandle()->components())for(int k=0;k<=128;k++){const auto p=c.curve->sample(k/128.);curves<<rz<<','<<step<<','<<flow.state().time<<','<<c.id<<','<<k/128.<<','<<p.x<<','<<p.y<<'\n';}};save(0);
        for(int step=1;step<=3;step++){
            constexpr double dt=.001;const HydroStepTime time{step,flow.state().time,dt};const auto old=flow.interfaceState().parametricGeometry().boundaryHandle();
            const auto load=flow.makeLoad(mesh,domain,zero,time.begin+.5*dt);auto candidate=flow.prepareNavierStokesStep(mesh,domain,time,load,options);const auto r=candidate.report();
            require(flow.interfaceState().parametricGeometry().boundaryHandle()==old,"multicomponent preparation changed accepted geometry");flow.accept(mesh,domain,std::move(candidate));
            const auto& current=flow.interfaceState().parametricGeometry();for(std::size_t i=0;i<2;i++){require(current.boundaryHandle()->components()[i].id==i+1,"midpoint changed component identity");drift=std::max(drift,std::abs(volume(current,i)-v0[i]));}
            fraction=std::max(fraction,current.report().maximumInventoryFractionError);equation=std::max(equation,r.equationResidual);energy=std::max(energy,std::abs(r.energyDefect));surfaceError=std::max(surfaceError,std::abs(r.surfaceWorkDefect));
            if(history){history<<rz<<','<<step<<','<<flow.state().time<<','<<volume(current,0)<<','<<volume(current,1)<<','<<current.report().maximumInventoryFractionError<<','<<r.equationResidual<<','<<r.kineticAfter<<','<<r.surfaceEnergyAfter<<','<<r.energyDefect<<','<<r.surfaceWorkDefect<<'\n';history.flush();}save(step);
            std::cout<<"MOMENT_COMPONENT_STEP rz="<<rz<<" step="<<step<<" iterations="<<r.nonlinearIterations<<" volume_error="<<drift<<" fraction="<<fraction<<" equation="<<r.equationResidual<<std::endl;
        }
        double motion=0;for(std::size_t i=0;i<2;i++)for(int k=0;k<128;k++)motion=std::max(motion,P::norm(flow.interfaceState().parametricGeometry().boundaryHandle()->components()[i].curve->sample(k/128.)-initial[i].sample(k/128.)));
        const double massChange=(flow.kineticMass()-mass0).norm()/mass0.norm();
        std::cout<<"MOMENT_COMPONENTS rz="<<rz<<" steps=3 shared_cells="<<shared<<" pressure_error="<<pressureError<<" constant_variation="<<variation<<" volume_error="<<drift<<" fraction="<<fraction<<" equation="<<equation<<" kinetic_budget="<<energy<<" surface_work_defect="<<surfaceError<<" mass_change="<<massChange<<" motion="<<motion<<std::endl;
        require(drift<2e-13&&fraction<3.1e-12&&equation<1e-10&&energy<2e-10&&massChange>1e-9&&motion>1e-8,"multicomponent active midpoint failed");
    }
}
void hydroParametricComponentsTransferContract(){
    namespace P=Quantum::Geometry::Parametric;
    for(bool rz:{false,true}){
        QuadtreeMesh source(input(4)),target(input(4));DictionaryMesh::RefineBox box;box.min={0,.5};box.max={.5,1};box.level=1;target.applyRefineBoxes({box});
        FaceSegmentGraph sg,tg;sg.build(source);tg.build(target);const auto sd=NativeLeafDomain::build(source,sg,rz),td=NativeLeafDomain::build(target,tg,rz);
        std::vector<P::Curve> curves;
        curves.push_back(rz?InterfacePrescribedData::parametricSpheroid(.42,{.18,.12},16):InterfacePrescribedData::parametricEllipse({.47,.42},{.14,.12},16));
        curves.push_back(rz?InterfacePrescribedData::parametricSpheroid(.67,{.12,.09},16):InterfacePrescribedData::parametricEllipse({.47,.67},{.10,.09},16));
        const auto original=InterfaceParametricMotion::initialize(source,sd,std::move(curves));const auto identity=original.parametricGeometry().boundaryHandle();
        auto refined=InterfaceParametricMotion::prepareMeshTransfer(source,sd,original,target,td,{});
        auto returned=InterfaceParametricMotion::prepareMeshTransfer(target,td,refined.state,source,sd,{});
        require(original.parametricGeometry().boundaryHandle()==identity&&refined.state.parametricGeometry().boundaryHandle()==identity&&returned.state.parametricGeometry().boundaryHandle()==identity,"component remap changed authoritative geometry");
        double drift=0,local=0;for(std::size_t i=0;i<2;i++){
            const auto& a=original.parametricGeometry().componentInventory(i);const auto& b=refined.state.parametricGeometry().componentInventory(i);const auto& c=returned.state.parametricGeometry().componentInventory(i);
            const double v=std::accumulate(a.begin(),a.end(),0.);drift=std::max({drift,std::abs(std::accumulate(b.begin(),b.end(),0.)-v),std::abs(std::accumulate(c.begin(),c.end(),0.)-v)});
            for(std::size_t k=0;k<a.size();k++)local=std::max(local,std::abs(a[k]-c[k])/original.parametricGeometry().grid().capacity(k));
        }
        require(drift<2e-13&&local<3.1e-12&&refined.report.overlaps.size()==28&&refined.state.progress().endTime==0,"component remap lost per-material inventory or duplicated geometric overlaps");
        std::cout<<std::setprecision(17)<<"PARAMETRIC_COMPONENT_TRANSFER rz="<<rz<<" leaves=16,28,16 volume_error="<<drift<<" local_return_fraction="<<local<<" fraction="<<refined.report.maximumFractionError<<" overlaps="<<refined.report.overlaps.size()<<std::endl;
    }
}
void hydroMomentAmrTwoPhaseContract(){
    namespace P=Quantum::Geometry::Parametric;constexpr int steps=6;constexpr double dt=.002;
    std::ofstream curves,history,cells;
    if(const char* prefix=std::getenv("QUANTUM_MOMENT_AMR_OUTPUT")){
        curves.open(std::string(prefix)+"_curves.csv");history.open(std::string(prefix)+"_history.csv");cells.open(std::string(prefix)+"_cells.csv");
        require(curves&&history&&cells,"cannot open native C1 two-phase output");curves<<std::setprecision(17)<<"rz,step,time,parameter,x,y\n";
        history<<std::setprecision(17)<<"rz,step,time,kinetic,surface,viscous,volume_drift,fraction_error,equation,mass_closure,metric_defect,kinetic_budget,total_energy_budget,segments\n";
        cells<<std::setprecision(17)<<"rz,cell,left,right,bottom,top\n";
    }
    for(bool rz:{false,true}){
        QuadtreeMesh mesh(input(4));DictionaryMesh::RefineBox box;box.min={0,.5};box.max={.5,1};box.level=1;mesh.applyRefineBoxes({box});
        FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);const auto space=Space::build(mesh,domain);
        auto initial=InterfaceParametricMotion::initialize(mesh,domain,rz?InterfacePrescribedData::parametricSpheroid(.5,{.22,.17},16):InterfacePrescribedData::parametricEllipse({.48,.52},{.18,.14},16));
        const auto first=initial.parametricGeometry().geometry();const double volume0=std::accumulate(initial.parametricGeometry().inventory().begin(),initial.parametricGeometry().inventory().end(),0.);
        auto hydro=Flow::createTwoPhase(mesh,domain,space,std::move(initial),{1,.1,.01,.001},Eigen::VectorXd::Zero(space.mass().rows()),1);
        const Space::Matrix firstMass=hydro.kineticMass();Flow::Options options;options.interfaceProjection.materialSegments=16;
        if(cells)for(int index:domain.cellIndexByRow()){const auto& c=mesh.cells(mesh.cellIndexToLabel[index]);const auto p=c.centerCoord(),h=c.size();cells<<rz<<','<<index<<','<<p.x-.5*h.x<<','<<p.x+.5*h.x<<','<<p.y-.5*h.y<<','<<p.y+.5*h.y<<'\n';}
        const auto save=[&](int step){if(curves){for(int k=0;k<=256;k++){const double s=k/256.;const auto p=hydro.interfaceState().parametricGeometry().geometry().sample(s);curves<<rz<<','<<step<<','<<hydro.state().time<<','<<s<<','<<p.x<<','<<p.y<<'\n';}curves.flush();}};save(0);
        double drift=0,fraction=0,equation=0,mass=0,metric=0,kineticBudget=0,surfaceBudget=0,viscous=0,temporal=0,material=0,energy0=0,total=0;
        for(int step=1;step<=steps;step++){
            const HydroStepTime time{step,hydro.state().time,dt,HydroTimePolicy::Fixed};const auto load=hydro.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),time.begin+.5*dt);
            const auto old=hydro.interfaceState().parametricGeometry().geometryHandle();auto candidate=hydro.prepareNavierStokesStep(mesh,domain,time,load,options);const auto r=candidate.report();
            require(hydro.interfaceState().parametricGeometry().geometryHandle()==old,"native C1 candidate changed accepted geometry");hydro.accept(mesh,domain,std::move(candidate));
            const auto& phase=hydro.interfaceState().parametricGeometry();const double volume=std::accumulate(phase.inventory().begin(),phase.inventory().end(),0.);
            drift=std::max(drift,std::abs(volume-volume0));fraction=std::max(fraction,phase.report().maximumInventoryFractionError);equation=std::max(equation,r.equationResidual);
            mass=std::max(mass,r.phaseMassClosure);metric=std::max(metric,r.massMetricDefect);kineticBudget=std::max(kineticBudget,std::abs(r.energyDefect));surfaceBudget=std::max(surfaceBudget,std::abs(r.surfaceWorkDefect));
            viscous+=r.viscousDissipation;temporal+=r.temporalMassWork;material+=r.materialMetricWork;if(step==1)energy0=r.surfaceEnergyBefore;
            const double budget=r.kineticAfter+r.surfaceEnergyAfter+viscous-temporal+material-energy0;total=std::max(total,std::abs(budget));
            require(r.transportDissipation==0&&r.viscousDissipation>=0,"native C1 midpoint added hidden dissipation");
            if(history){history<<rz<<','<<step<<','<<hydro.state().time<<','<<r.kineticAfter<<','<<r.surfaceEnergyAfter<<','<<viscous<<','<<volume-volume0<<','<<phase.report().maximumInventoryFractionError<<','<<r.equationResidual<<','<<r.phaseMassClosure<<','<<r.massMetricDefect<<','<<r.energyDefect<<','<<budget<<','<<phase.geometry().arcs.size()<<'\n';history.flush();}save(step);
            std::cout<<std::setprecision(17)<<"MOMENT_AMR_TWO_PHASE_STEP metric="<<(rz?"RZ":"XY")<<" step="<<step<<" time="<<hydro.state().time<<" iterations="<<r.nonlinearIterations<<" kinetic="<<r.kineticAfter<<" volume_drift="<<volume-volume0<<" fraction="<<phase.report().maximumInventoryFractionError<<" equation="<<r.equationResidual<<" energy_budget="<<budget<<std::endl;
        }
        double motion=0;for(int k=0;k<=256;k++)motion=std::max(motion,P::norm(hydro.interfaceState().parametricGeometry().geometry().sample(k/256.)-first.sample(k/256.)));
        const double massChange=(hydro.kineticMass()-firstMass).norm()/firstMass.norm();
        std::cout<<"MOMENT_AMR_TWO_PHASE metric="<<(rz?"RZ":"XY")<<" steps="<<steps<<" time="<<hydro.state().time<<" cells="<<mesh.cellIndexNum()<<" volume_drift="<<drift<<" fraction="<<fraction<<" equation="<<equation<<" mass_closure="<<mass<<" metric="<<metric<<" kinetic_budget="<<kineticBudget<<" surface_work_defect="<<surfaceBudget<<" total_energy_budget="<<total<<" mass_change="<<massChange<<" motion="<<motion<<" capillary_feedback=1 dynamic_mesh=0\n";
        require(drift<2e-13&&fraction<3.1e-12&&equation<1e-10&&mass<1e-12&&metric<2e-10&&kineticBudget<2e-10&&massChange>1e-6&&motion>1e-4,"native C1 active AMR midpoint exceeds physical budgets");
    }
}
void hydroMomentSplineFluxRootsContract(){fluxRoots();}
void hydroMomentPositiveFluxContract(){positiveIntegrals();fluxRoots(true);}
void hydroMomentRZSplineContract(){
    const auto analytic=[](Vector2D p,bool oscillatory){
        const double r=p.x,z=p.y,f=std::pow(1-r*r,2),fr=-4*r*(1-r*r),frr=-4+12*r*r;
        const double g=z*z*(1-z)*(1-z),gz=2*z-6*z*z+4*z*z*z,gzz=2-12*z+12*z*z;
        const double k=oscillatory?1+.2*std::cos(3*r):1,kr=oscillatory?-.6*std::sin(3*r):0,krr=oscillatory?-1.8*std::cos(3*r):0;
        const double F=4*f*k,Fr=4*(fr*k+f*kr),Frr=4*(frr*k+2*fr*kr+f*krr);
        return Spline::FirstJet{{r*F*gz,-(2*F+r*Fr)*g},{(F+r*Fr)*gz,-(3*Fr+r*Frr)*g},{r*F*gzz,-(2*F+r*Fr)*gz}};
    };
    double maximumValue=0,maximumJet=0,maximumMoments=0,maximumTrace=0,maximumJoin=0,maximumDivergence=0,maximumWork=0,previousValue=0,previousJet=0;
    for(int n:{4,8,16}){
        double gridValue=0,gridJet=0;
        QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
        for(bool nonpolynomial:{false,true}){
            const auto u=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=analytic(p,nonpolynomial).value;return Eigen::Vector2d(v.x,v.y);});
            const auto spline=Spline::build(mesh,domain,space,u,.125);
            const auto recovered=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=spline.firstJet(p).value;return Eigen::Vector2d(v.x,v.y);});
            maximumMoments=std::max(maximumMoments,(recovered-u).norm()/u.norm());
            std::vector<Spline::PointForce> forces;double work=0;
            for(int j=0;j<19;j++)for(int i=0;i<17;i++){
                const Vector2D p{(i+.37)/17,(j+.29)/19};const auto v=spline.firstJet(p),a=analytic(p,nonpolynomial);
                gridValue=std::max(gridValue,std::hypot(v.value.x-a.value.x,v.value.y-a.value.y));
                gridJet=std::max({gridJet,std::hypot(v.dx.x-a.dx.x,v.dx.y-a.dx.y),std::hypot(v.dy.x-a.dy.x,v.dy.y-a.dy.y)});
                maximumDivergence=std::max(maximumDivergence,std::abs(v.dx.x+v.dy.y+v.value.x/p.x));
                forces.push_back({p,{std::sin(i+.2*j),std::cos(j-.3*i)},{.1*std::cos(i),-.2*std::sin(j)},{.15*std::sin(i),.12*std::cos(j)}});
            }
            for(Vector2D p:std::vector<Vector2D>{{0,.43},{1,.43},{.31,0},{.31,1},{0,0},{0,1}})forces.push_back({p,{.3,-.2},{.1,.7},{-.2,.4}});
            for(const auto& f:forces){const auto v=spline.firstJet(f.point);work+=f.measure.x*v.value.x+f.measure.y*v.value.y+f.dx.x*v.dx.x+f.dx.y*v.dx.y+f.dy.x*v.dy.x+f.dy.y*v.dy.y;}
            maximumWork=std::max(maximumWork,std::abs(u.dot(spline.transposePointForces(forces))-work)/std::max(1.,std::abs(work)));
            const auto atAxis=spline.firstJet({0,.43});require(atAxis.value.x==0&&atAxis.dx.y==0&&atAxis.dy.x==0&&std::abs(atAxis.value.y)>.05,"RZ C1 reconstruction imposed no slip or a singular axis");
            maximumDivergence=std::max(maximumDivergence,std::abs(2*atAxis.dx.x+atAxis.dy.y));
            for(int k=1;k<n;k++)for(int j=0;j<7;j++)for(int axis=0;axis<2;axis++){
                Vector2D left{double(k)/n,(j+.41)/7},right=left;
                left.x=std::nextafter(left.x,0.);right.x=std::nextafter(right.x,1.);if(axis){std::swap(left.x,left.y);std::swap(right.x,right.y);}
                const auto a=spline.firstJet(left),b=spline.firstJet(right);maximumJoin=std::max({maximumJoin,std::hypot(a.value.x-b.value.x,a.value.y-b.value.y),std::hypot(a.dx.x-b.dx.x,a.dx.y-b.dx.y),std::hypot(a.dy.x-b.dy.x,a.dy.y-b.dy.y)});
            }
            const auto flux=spline.normalFlux();const auto& rule=Quantum::Geometry::Parametric::gaussEight();const double twoPi=2*std::acos(-1.);
            for(const auto face:domain.faces(mesh)){
                const auto& info=space.faces()[face.nativeRow];const auto chart=flux.faceChart(face.nativeRow);const double length=face.horizontal?chart.end.x-chart.begin.x:chart.end.y-chart.begin.y;double moments[2]{};
                for(int half=0;half<2;half++)for(int k=0;k<8;k++){
                    const double t=.5*(1+rule.x[k]);const auto a=half?chart.midpoint:chart.begin,b=half?chart.end:chart.midpoint;const auto p=a+(b-a)*t;
                    const double density=double(flux.normalFluxDensity(face.nativeRow,half,false,t));const auto v=spline.firstJet(p).value;
                    maximumTrace=std::max(maximumTrace,std::abs(2*density/length-twoPi*p.x*(face.horizontal?v.y:v.x)));
                    if(info.first>=0){moments[0]+=.5*rule.w[k]*density/info.area;moments[1]+=.5*rule.w[k]*density/info.area*((face.horizontal?p.x:p.y)-info.tangentCenter)/info.tangentScale;}
                }
                if(info.first>=0)maximumTrace=std::max({maximumTrace,std::abs(moments[0]-u[info.first]),std::abs(moments[1]-u[info.first+1])});
            }
        }
        // Minimum physical energy does not interpolate every quartic chi.
        // Test approximation to independent smooth velocities and gradients;
        // retain the exact owned-moment, C1, axis and adjoint checks separately.
        const double valueOrder=previousValue?std::log2(previousValue/gridValue):0,jetOrder=previousJet?std::log2(previousJet/gridJet):0;
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_C1_ACCURACY n="<<n<<" value="<<gridValue<<" first_jet="<<gridJet<<" value_order="<<valueOrder<<" first_jet_order="<<jetOrder<<std::endl;
        if(previousValue)require(valueOrder>1.8&&jetOrder>.8,"RZ minimum-energy common velocity lacks second-order value/first-order gradient accuracy");
        previousValue=gridValue;previousJet=gridJet;maximumValue=std::max(maximumValue,gridValue);maximumJet=std::max(maximumJet,gridJet);
    }
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_C1_SPLINE value="<<maximumValue<<" first_jet="<<maximumJet<<" all_moments="<<maximumMoments
        <<" weighted_trace="<<maximumTrace<<" C1_join="<<maximumJoin<<" weighted_divergence="<<maximumDivergence<<" value_gradient_work="<<maximumWork<<std::endl;
    require(maximumMoments<2e-10&&maximumTrace<2e-10&&maximumJoin<2e-8&&maximumDivergence<2e-12&&maximumWork<2e-10,"RZ moment C1 geometry velocity or adjoint work failed");
    // Consume actual RZ Navier--Stokes midpoints, not just interpolated data.
    // This checks the motion port; RZ phase inventory is still a later layer.
    QuadtreeMesh mesh(input(8));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
    const auto initial=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=analytic(p,false).value;return Eigen::Vector2d(v.x,v.y);});
    auto stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),.02));
    auto owner=Flow::create(mesh,domain,space,std::move(stress),1.3,initial);const auto base=Spline::build(mesh,domain,space,initial,0);
    Vector2D point{.413,.47},axisPoint{0,.47};const auto first=point;double momentError=0,metricError=0,equation=0;
    for(int step=1;step<=3;step++){
        HydroStepTime time{step,owner.state().time,.01,HydroTimePolicy::Fixed};
        const auto load=owner.makeLoad(mesh,domain,Eigen::VectorXd::Zero(initial.size()),time.begin+.5*time.dt);
        auto candidate=owner.prepareNavierStokesStep(mesh,domain,time,load);equation=std::max(equation,candidate.report().equationResidual);owner.accept(mesh,domain,std::move(candidate));
        const auto middle=*owner.state().midpointPhysicalVelocity;
        const auto recovered=space.interpolateVelocity(mesh,domain,[&](Vector2D p){const auto v=middle.firstJet(p).value;return Eigen::Vector2d(v.x,v.y);});
        momentError=std::max(momentError,(recovered-owner.state().midpointVelocity).norm()/owner.state().midpointVelocity.norm());
        const auto map=HydroInterfaceTransport::frozenParametricFlow(middle);const auto moved=map(point,time.dt),axisMoved=map(axisPoint,time.dt);
        const double determinant=moved.dx.x*moved.dy.y-moved.dx.y*moved.dy.x;
        metricError=std::max({metricError,std::abs(moved.point.x/point.x*determinant-1),std::abs(axisMoved.dx.x*axisMoved.dx.x*axisMoved.dy.y-1)});
        require(axisMoved.point.x==0&&axisMoved.dx.y==0&&axisMoved.dy.x==0,"actual midpoint motion lost its regular axis");point=moved.point;axisPoint=axisMoved.point;
    }
    const double displacement=std::hypot(point.x-first.x,point.y-first.y);
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_C1_FLOW steps=3 time="<<owner.state().time<<" midpoint_moments="<<momentError<<" weighted_map_jacobian="<<metricError
        <<" equation="<<equation<<" particle_displacement="<<displacement<<" axis_position="<<axisPoint.x<<" phase_inventory=0\n";
    require(momentError<2e-10&&metricError<2e-8&&displacement>1e-4,"actual RZ midpoint geometry-velocity port failed");
}
void hydroMomentRZInterfaceContract(bool motion){
    namespace P=Quantum::Geometry::Parametric;
    using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    using Traction=Quantum::Hydrodynamics::HydroCompatibleMomentTraction;
    using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;
    constexpr int n=8,steps=6;constexpr double dt=.02;const double pi=std::acos(-1.);
    QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
    double previousCurvature=0;
    for(int segments:{16,32}){
        auto sphere=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricSpheroid(.5,{.2,.2},segments));
        const auto& candidate=sphere.parametricGeometry();const auto surface=Traction::measure(candidate.geometry(),true);
        const double volume=std::accumulate(candidate.inventory().begin(),candidate.inventory().end(),0.);
        double curvature=0,weight=0;const auto& rule=P::gaussEight();
        for(const auto& arc:candidate.geometry().arcs)for(int k=0;k<8;k++){const auto q=P::sampleSurfaceRZ(arc,.5*(1+rule.x[k]));
            curvature+=rule.w[k]*q.point.x*q.speed*std::pow(q.curvature-10,2);weight+=rule.w[k]*q.point.x*q.speed;}
        curvature=std::sqrt(curvature/weight)/10;
        const double pole=std::max(std::abs(P::sampleSurfaceRZ(candidate.geometry().arcs.front(),0,true).curvature-10),std::abs(P::sampleSurfaceRZ(candidate.geometry().arcs.back(),1,true).curvature-10))/10;
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_SPHERE segments="<<segments<<" volume_relative="<<std::abs(volume/(4*pi*.008/3)-1)
            <<" area_relative="<<std::abs(surface.length/(4*pi*.04)-1)<<" curvature_relative_rms="<<curvature<<" pole_curvature_relative="<<pole<<std::endl;
        require(std::abs(volume/(4*pi*.008/3)-1)<2e-6&&std::abs(surface.length/(4*pi*.04)-1)<2e-6&&curvature<.001&&pole<.001,"RZ sphere measure or total curvature failed");
        if(previousCurvature>0)require(curvature<.3*previousCurvature,"RZ sphere curvature did not improve with material refinement");previousCurvature=curvature;
    }
    const auto initial=space.interpolateVelocity(mesh,domain,[](Vector2D p){const double r=p.x,z=p.y,f=1-r*r,g=z*z*(1-z)*(1-z),gz=2*z*(1-z)*(1-2*z);
        return Eigen::Vector2d(4*r*f*f*gz,-8*f*f*g+16*r*r*f*g);});
    const auto base=Spline::build(mesh,domain,space,initial,0);
    auto interface=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricSpheroid(.5,{.22,.17},16));
    const auto first=interface.parametricGeometry().geometry();const double firstVolume=std::accumulate(interface.parametricGeometry().inventory().begin(),interface.parametricGeometry().inventory().end(),0.);
    const auto material=Material::buildFrozen(mesh,domain,space,interface.parametricGeometry().geometryHandle(),0,{5,.1,.02,.01});
    const auto constant=Traction::build(mesh,domain,space,base,material,[](const P::SurfaceSample&){return Traction::Components{3,0};});
    const auto force=Traction::build(mesh,domain,space,base,material,[](const P::SurfaceSample& q){return Traction::Components{q.curvature,0};},1);
    long double independentWork=0;const auto& rule=P::gaussEight();
    for(const auto& arc:first.arcs)for(int panel=0;panel<128;panel++)for(int k=0;k<8;k++){
        const double t=(panel+.5+.5*rule.x[k])/128;const auto q=P::sampleSurfaceRZ(arc,t);const auto u=base.firstJet({q.point.x,q.point.y}).value;
        independentWork+=(.5L*rule.w[k]/128)*2*pi*q.point.x*(q.curvature-force.report.meanNormal)*(-q.derivative.y*u.x+q.derivative.x*u.y);
    }
    const double virtualWork=std::abs(initial.dot(force.variation)-double(independentWork));
    require(virtualWork<2e-9&&constant.variation.norm()<2e-11,"RZ sharp traction or hoop surface work is inconsistent");
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_SURFACE_WORK defect="<<virtualWork<<" constant_variation="<<constant.variation.norm()<<" surface_area="<<force.report.length<<std::endl;
    if(!motion)return;
    auto stress=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity::buildCellwise(mesh,domain,space,std::vector<double>(mesh.cellIndexNum(),.02));
    auto hydro=Flow::create(mesh,domain,space,std::move(stress),1.3,initial);
    InterfaceParametricProjection::Options options;options.materialSegments=16;
    std::ofstream curves,history;if(const char* prefix=std::getenv("QUANTUM_MOMENT_RZ_INTERFACE_OUTPUT")){
        curves.open(std::string(prefix)+"_curves.csv");curves<<std::setprecision(17)<<"step,time,parameter,r,z\n";
        history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)<<"step,time,volume_drift,fraction_error,mass_flux_closure,axis_radius,axis_tangent,material_segments\n";}
    const auto save=[&](int step){if(curves)for(int k=0;k<=256;k++){const double s=double(k)/256;const auto p=interface.parametricGeometry().geometry().sample(s);curves<<step<<','<<hydro.state().time<<','<<s<<','<<p.x<<','<<p.y<<'\n';}if(curves)curves.flush();};save(0);
    double drift=0,fraction=0,mass=0,axisRadius=0,axisTangent=0;
    for(int step=1;step<=steps;step++){
        const HydroStepTime time{step,hydro.state().time,dt,HydroTimePolicy::Fixed};const auto zero=hydro.makeLoad(mesh,domain,Eigen::VectorXd::Zero(initial.size()),time.begin+.5*dt);
        auto advanced=hydro.prepareNavierStokesStep(mesh,domain,time,zero);hydro.accept(mesh,domain,std::move(advanced));
        const auto middle=*hydro.state().midpointPhysicalVelocity;
        const auto normal=middle.normalFlux();const auto map=HydroInterfaceTransport::frozenParametricFlow(middle);
        const InterfaceStepTime clock{step,time.begin,time.begin+dt,dt};const auto oldIdentity=interface.parametricGeometry().geometryHandle();
        auto next=InterfaceParametricMotion::prepare(mesh,domain,interface,clock,dt,map,normal,options);
        require(interface.parametricGeometry().geometryHandle()==oldIdentity&&next.state.parametricFollows(interface,clock),"RZ motion lost transactional geometry/time identity");
        const auto flux=Transport::phaseMassFlux(mesh,domain,space,interface,next.state,clock,5,.1);mass=std::max(mass,flux.maximumRelativeClosure);
        interface=std::move(next.state);const auto& candidate=interface.parametricGeometry();const auto& curve=candidate.geometry();
        require(curve.boundary.axisPair()&&candidate.grid().axisymmetric,"RZ motion discarded its measure or physical meridian topology");
        const double volume=std::accumulate(candidate.inventory().begin(),candidate.inventory().end(),0.);
        drift=std::max(drift,std::abs(volume-firstVolume));fraction=std::max(fraction,candidate.report().maximumInventoryFractionError);
        for(const auto& pole:std::vector<std::pair<const P::Arc*,double>>{{&curve.arcs.front(),0},{&curve.arcs.back(),1}}){const auto p=pole.first->sampleFraction(pole.second),d=pole.first->sampleFraction(pole.second,1);
            axisRadius=std::max(axisRadius,std::abs(p.x));axisTangent=std::max(axisTangent,std::abs(d.y));P::sampleSurfaceRZ(*pole.first,pole.second,true);}
        if(history){history<<step<<','<<hydro.state().time<<','<<volume-firstVolume<<','<<candidate.report().maximumInventoryFractionError<<','<<flux.maximumRelativeClosure<<','<<axisRadius<<','<<axisTangent<<','<<curve.arcs.size()<<'\n';history.flush();}save(step);
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_INTERFACE_STEP step="<<step<<" time="<<hydro.state().time<<" volume_drift="<<volume-firstVolume<<" fraction_error="<<candidate.report().maximumInventoryFractionError<<" arcs="<<curve.arcs.size()<<std::endl;
    }
    double displacement=0;for(int k=0;k<=256;k++)displacement=std::max(displacement,P::norm(interface.parametricGeometry().geometry().sample(k/256.)-first.sample(k/256.)));
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_INTERFACE steps="<<steps<<" time="<<hydro.state().time<<" volume_drift="<<drift<<" fraction_error="<<fraction<<" mass_flux_closure="<<mass
        <<" axis_radius="<<axisRadius<<" axis_tangent="<<axisTangent<<" displacement="<<displacement<<" active_two_phase=0\n";
    require(drift<2e-13&&fraction<3.1e-12&&mass<1e-12&&axisRadius<2e-12&&axisTangent<2e-8&&displacement>.01,"RZ physical meridian inventory trajectory failed");
}
void hydroMomentRZTwoPhaseContract(){
    namespace P=Quantum::Geometry::Parametric;
    constexpr int n=8,steps=6;constexpr double dt=.002;
    QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
    auto initial=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricSpheroid(.5,{.22,.17},16));
    const auto first=initial.parametricGeometry().geometry();const double firstVolume=std::accumulate(initial.parametricGeometry().inventory().begin(),initial.parametricGeometry().inventory().end(),0.);
    auto hydro=Flow::createTwoPhase(mesh,domain,space,std::move(initial),{1,.1,.01,.001},Eigen::VectorXd::Zero(space.mass().rows()),1);
    const Space::Matrix firstMass=hydro.kineticMass();Flow::Options options;options.interfaceProjection.materialSegments=16;
    std::ofstream curves,history;if(const char* prefix=std::getenv("QUANTUM_MOMENT_RZ_TWO_PHASE_OUTPUT")){
        curves.open(std::string(prefix)+"_curves.csv");curves<<std::setprecision(17)<<"step,time,parameter,r,z\n";
        history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)<<"step,time,kinetic_energy,surface_energy,viscous_dissipation,volume_drift,fraction_error,equation,mass_closure,mass_metric_defect,kinetic_budget,total_energy_budget,material_segments\n";}
    const auto save=[&](int step){if(curves){for(int k=0;k<=256;k++){const double s=k/256.;const auto p=hydro.interfaceState().parametricGeometry().geometry().sample(s);curves<<step<<','<<hydro.state().time<<','<<s<<','<<p.x<<','<<p.y<<'\n';}curves.flush();}};save(0);
    double maxVolume=0,maxFraction=0,maxEquation=0,maxMass=0,maxMetric=0,maxKineticBudget=0,maxSurfaceBudget=0,viscous=0,temporal=0,materialWork=0,totalEnergy=0;
    for(int step=1;step<=steps;step++){
        const HydroStepTime time{step,hydro.state().time,dt,HydroTimePolicy::Fixed};const auto zero=hydro.makeLoad(mesh,domain,Eigen::VectorXd::Zero(space.mass().rows()),time.begin+.5*dt);
        const auto old=hydro.interfaceState().parametricGeometry().geometryHandle();auto candidate=hydro.prepareNavierStokesStep(mesh,domain,time,zero,options);const auto report=candidate.report();
        require(hydro.interfaceState().parametricGeometry().geometryHandle()==old,"RZ coupled preparation changed accepted geometry");hydro.accept(mesh,domain,std::move(candidate));
        const auto& phase=hydro.interfaceState().parametricGeometry();const double volume=std::accumulate(phase.inventory().begin(),phase.inventory().end(),0.);
        maxVolume=std::max(maxVolume,std::abs(volume-firstVolume));maxFraction=std::max(maxFraction,phase.report().maximumInventoryFractionError);
        maxEquation=std::max(maxEquation,report.equationResidual);maxMass=std::max(maxMass,report.phaseMassClosure);maxMetric=std::max(maxMetric,report.massMetricDefect);
        maxKineticBudget=std::max(maxKineticBudget,std::abs(report.energyDefect));maxSurfaceBudget=std::max(maxSurfaceBudget,std::abs(report.surfaceWorkDefect));
        // The flow report already contains time-integrated viscous work.
        viscous+=report.viscousDissipation;temporal+=report.temporalMassWork;materialWork+=report.materialMetricWork;
        if(step==1)totalEnergy=report.surfaceEnergyBefore;const double budget=report.kineticAfter+report.surfaceEnergyAfter+viscous-temporal+materialWork-totalEnergy;
        if(history){history<<step<<','<<hydro.state().time<<','<<report.kineticAfter<<','<<report.surfaceEnergyAfter<<','<<viscous<<','<<volume-firstVolume<<','<<phase.report().maximumInventoryFractionError
            <<','<<report.equationResidual<<','<<report.phaseMassClosure<<','<<report.massMetricDefect<<','<<report.energyDefect<<','<<budget<<','<<phase.geometry().arcs.size()<<'\n';history.flush();}save(step);
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_TWO_PHASE_STEP step="<<step<<" time="<<hydro.state().time<<" iterations="<<report.nonlinearIterations<<" kinetic="<<report.kineticAfter
            <<" volume_drift="<<volume-firstVolume<<" fraction_error="<<phase.report().maximumInventoryFractionError<<" equation="<<report.equationResidual<<" metric="<<report.massMetricDefect<<" total_energy_budget="<<budget<<" arcs="<<phase.geometry().arcs.size()<<std::endl;
    }
    double motion=0;for(int k=0;k<=256;k++)motion=std::max(motion,P::norm(hydro.interfaceState().parametricGeometry().geometry().sample(k/256.)-first.sample(k/256.)));
    const double massChange=(hydro.kineticMass()-firstMass).norm()/firstMass.norm();
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_TWO_PHASE steps="<<steps<<" time="<<hydro.state().time<<" volume_drift="<<maxVolume<<" fraction_error="<<maxFraction<<" equation="<<maxEquation
        <<" phase_mass_closure="<<maxMass<<" mass_metric_defect="<<maxMetric<<" kinetic_budget="<<maxKineticBudget<<" surface_work_defect="<<maxSurfaceBudget<<" mass_change="<<massChange<<" motion="<<motion<<" variable_mass=1 capillary_feedback=1\n";
    require(maxVolume<2e-13&&maxFraction<3.1e-12&&maxEquation<1e-10&&maxMass<1e-12&&maxMetric<2e-10&&maxKineticBudget<2e-10&&massChange>1e-6&&motion>1e-4,"RZ actual moving-mass/capillary midpoint failed");
}
// Read the complete candidate q, never reconstruct it from its lower moments.
// This diagnostic explicitly supplies uniform XY/free-slip boundary data;
// these diagnostic files do not describe arbitrary boundary conditions.
void hydroMomentFreeSlipPhaseProbe(const char* path,const std::string& task){
    const bool resolutionOnly=task=="resolution",materialOnly=task=="material",initialLayers=task=="initialLayers";
    require(task=="velocity"||resolutionOnly||materialOnly||initialLayers,"unknown saved phase diagnostic");
    namespace P=Quantum::Geometry::Parametric;
    using Boundary=Quantum::Hydrodynamics::HydroCompatibleMomentBoundary;
    using Viscosity=Quantum::Hydrodynamics::HydroCompatibleMomentViscosity;
    std::ifstream in(path);std::string magic;in>>magic;
    const bool trialMaterial=magic=="MOMENT_MATERIAL_COMPLETE_SNAPSHOT_V1";
    require(magic=="MOMENT_PHASE_COMPLETE_CANDIDATE_V2"||(materialOnly&&trialMaterial),"phase probe requires complete V2 coordinates or a material-only snapshot");
    int rz,uniform,cells,count,step;double begin,dt,l,b,r,t,hx,hy;
    in>>rz>>uniform>>cells>>count>>step>>begin>>dt>>l>>b>>r>>t>>hx>>hy;
    require(bool(in)&&!rz&&uniform&&cells>0&&count>0&&r>l&&t>b&&hx>0&&hy>0&&dt>0,"phase probe requires a uniform XY candidate");
    const int nx=int(std::lround((r-l)/hx)),ny=int(std::lround((t-b)/hy));
    require(nx>0&&ny>0&&nx*ny==cells,"phase probe grid dimensions disagree");
    Eigen::VectorXd wind(count);for(int i=0;i<count;i++)in>>wind[i];
    P::Curve curve;std::size_t arcs;bool axisEnds;in>>axisEnds>>arcs;require(!axisEnds&&arcs>0,"phase probe requires a closed XY curve");
    for(std::size_t i=0;i<arcs;i++){P::Arc arc;in>>arc.a>>arc.b;for(auto& p:arc.p)in>>p.x>>p.y;curve.arcs.push_back(arc);}
    std::size_t knots;in>>knots;double previous=0,minimum=1;std::vector<double> retained;
    for(std::size_t i=0;i<knots;i++){double value;in>>value;require(value>=previous&&value<=1,"invalid saved material knots");if(value>previous)minimum=std::min(minimum,value-previous);previous=value;retained.push_back(value);}
    require(bool(in)&&wind.allFinite(),"incomplete complete-velocity checkpoint");
    auto dictionary=input(nx);dictionary.Ny=ny;dictionary.minCoord={l,b};dictionary.maxCoord={r,t};dictionary.maxCellNum=std::max(10000,cells);
    QuadtreeMesh mesh(dictionary);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const std::vector<HydroBoundaryPatch> walls{{"west",HydroBoundarySide::West,HydroBoundaryType::FreeSlip,b,t},
        {"east",HydroBoundarySide::East,HydroBoundaryType::FreeSlip,b,t},{"south",HydroBoundarySide::South,HydroBoundaryType::FreeSlip,l,r},
        {"north",HydroBoundarySide::North,HydroBoundaryType::FreeSlip,l,r}};
    const auto viscosity=Viscosity::buildCellwise(mesh,domain,space,std::vector<double>(cells,0));
    const auto boundary=Boundary::build(mesh,domain,space,viscosity,walls);
    const auto zero=Spline::initializeNative(mesh,domain,space,Eigen::VectorXd::Zero(space.mass().rows()),begin+.5*dt,&boundary);
    require(zero.coordinates().size()==wind.size(),"phase probe boundary/layout differs from saved complete field");
    const auto field=zero.withCoordinates(wind,begin+.5*dt);
    if(initialLayers){
        require(!trialMaterial&&step==1&&begin==0&&arcs==64&&nx==64&&ny==64&&dt==.001&&l==-.5&&b==-.5&&r==.5&&t==.5,
            "initial-layer replay requires original B3 step one; later independent inventories are absent from V2");
        const char* prefix=std::getenv("QUANTUM_MOMENT_LAYER_OUTPUT");require(prefix,"initial-layer replay requires an explicit output prefix");
        const auto initial=InterfaceParametricMotion::initialize(mesh,domain,curve);
        const auto& initialized=initial.parametricGeometry().geometry();require(initialized.arcs.size()==curve.arcs.size(),"initial replay changed the saved geometry");
        for(std::size_t a=0;a<arcs;a++)for(std::size_t k=0;k<curve.arcs[a].p.size();k++)
            require(initialized.arcs[a].p[k].x==curve.arcs[a].p[k].x&&initialized.arcs[a].p[k].y==curve.arcs[a].p[k].y,"initial replay changed a saved curve coefficient");
        InterfaceParametricProjection::Options options;options.materialSegments=64;options.relativeH2=.001*std::pow(16./64,2);options.minimumMaterialKnots=retained;
        const auto map=HydroInterfaceTransport::frozenParametricFlow(field);const auto flux=field.normalFlux();
        const auto started=std::chrono::steady_clock::now();
        auto half=InterfaceParametricMotion::prepare(mesh,domain,initial,{step,0,.5*dt,.5*dt},.5*dt,map,flux,options);
        auto full=InterfaceParametricMotion::prepare(mesh,domain,initial,{step,0,dt,dt},dt,map,flux,options);
        // Same nested half/end space as the production coupled owner. This
        // replay seeds inventory only at the actual initial time, never from
        // a later transported curve, and does not accept either candidate.
        for(;;){auto common=options.minimumMaterialKnots;
            for(const auto* state:{&half.state,&full.state})for(const auto& arc:state->parametricGeometry().geometry().arcs)common.push_back(arc.a);
            common.push_back(1);std::sort(common.begin(),common.end());common.erase(std::unique(common.begin(),common.end()),common.end());
            if(common==options.minimumMaterialKnots)break;options.minimumMaterialKnots=std::move(common);
            InterfaceParametricMotion::refineMaterialSpace(half,options.minimumMaterialKnots);InterfaceParametricMotion::refineMaterialSpace(full,options.minimumMaterialKnots);
        }
        const auto write=[&](const InterfaceState& state,const char* layer,double at){
            const auto& geometry=state.parametricGeometry().geometry();std::ofstream out(std::string(prefix)+"."+layer+".dat");require(bool(out),"cannot write initial-layer replay snapshot");
            out<<std::setprecision(17)<<"MOMENT_MATERIAL_COMPLETE_SNAPSHOT_V1\n"<<rz<<' '<<uniform<<' '<<cells<<' '<<count<<' '<<step<<' '<<at<<' '<<dt<<'\n'
                <<l<<' '<<b<<' '<<r<<' '<<t<<' '<<hx<<' '<<hy<<'\n';for(int i=0;i<count;i++)out<<wind[i]<<'\n';
            out<<0<<' '<<geometry.arcs.size()<<'\n';for(const auto& arc:geometry.arcs){out<<arc.a<<' '<<arc.b;for(const auto& point:arc.p)out<<' '<<point.x<<' '<<point.y;out<<'\n';}
            out<<options.minimumMaterialKnots.size()<<'\n';for(double s:options.minimumMaterialKnots)out<<s<<'\n';require(bool(out),"initial-layer replay snapshot write failed");
        };
        write(half.state,"half",.5*dt);write(full.state,"full",dt);
        std::cout<<std::setprecision(17)<<"MOMENT_INITIAL_LAYER_REPLAY source="<<path<<" half_arcs="<<half.state.parametricGeometry().geometry().arcs.size()
            <<" full_arcs="<<full.state.parametricGeometry().geometry().arcs.size()<<" half_balance="<<half.inventoryBalance<<" full_balance="<<full.inventoryBalance
            <<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<" phase_flux_evaluated=1 accepted=0"<<std::endl;
        return;
    }
    if(materialOnly){
        using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
        using Complement=Quantum::Hydrodynamics::HydroCompatibleMomentComplement;
        require(nx==64&&ny==64&&dt==.001&&l==-.5&&b==-.5&&r==.5&&t==.5,"material replay declares original B3 properties and grid");
        const auto material=Material::buildFrozen(mesh,domain,space,std::make_shared<const P::Curve>(curve),begin,{1,.001,0,0});
        const auto complement=Complement::build(mesh,domain,space,boundary.fixed());
        // Evaluate the saved candidate coordinate vector on the recorded
        // geometry's time slice. This is a frozen algebraic probe, not an
        // accepted velocity or a newly prepared midpoint material state.
        const auto sampledField=field.withCoordinates(wind,begin);
        const char* writePath=std::getenv("QUANTUM_MATERIAL_REFERENCE_WRITE");
        const char* readPath=std::getenv("QUANTUM_MATERIAL_REFERENCE_READ");
        require(!(writePath&&readPath),"select only one material reference operation");
        for(int pass=0;pass<2;pass++){
            const auto started=std::chrono::steady_clock::now();
            const auto forms=sampledField.mixedMaterialForms(mesh,domain,space,material,complement);
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            double maximum=0,relative=0;
            if(writePath&&pass==0){
                std::ofstream out(writePath,std::ios::binary);require(bool(out),"cannot create material matrix reference");
                for(const auto& cell:forms.cells){const int n=int(cell.mass.rows());
                    out.write(reinterpret_cast<const char*>(&cell.index),sizeof(int));out.write(reinterpret_cast<const char*>(&n),sizeof(int));
                    out.write(reinterpret_cast<const char*>(cell.mass.data()),sizeof(double)*n*n);}
                require(bool(out),"incomplete material matrix reference write");
            }
            if(readPath){
                std::ifstream in(readPath,std::ios::binary);require(bool(in),"cannot read material matrix reference");
                long double differenceSquared=0,referenceSquared=0;
                for(const auto& cell:forms.cells){int index=-1,n=0;in.read(reinterpret_cast<char*>(&index),sizeof(int));in.read(reinterpret_cast<char*>(&n),sizeof(int));
                    require(index==cell.index&&n==cell.mass.rows(),"material reference cell layout differs");
                    Eigen::MatrixXd reference(n,n);in.read(reinterpret_cast<char*>(reference.data()),sizeof(double)*n*n);require(bool(in)&&reference.allFinite(),"incomplete material matrix reference");
                    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
                        const double difference=cell.mass(i,j)-reference(i,j),scale=std::sqrt(reference(i,i)*reference(j,j));
                        if(scale>0)maximum=std::max(maximum,std::abs(difference)/scale);else require(difference==0,"reference has a changed zero mass direction");
                        differenceSquared+=static_cast<long double>(difference)*difference;
                        referenceSquared+=static_cast<long double>(reference(i,j))*reference(i,j);
                    }
                }
                require(in.peek()==std::char_traits<char>::eof(),"unexpected material reference tail");
                relative=std::sqrt(double(differenceSquared/referenceSquared));
                require(maximum<2e-11&&relative<2e-12,"material contraction differs from the saved full cell matrices");
            }
            std::cout<<std::setprecision(17)<<"MOMENT_MATERIAL_REPLAY source_step="<<step<<" geometry_time="<<begin<<" trial_geometry="<<trialMaterial<<" pass="<<pass
                <<" cells="<<cells<<" arcs="<<arcs<<" seconds="<<seconds<<" points="<<forms.quadraturePoints
                <<" phase_error="<<forms.maximumPhaseMomentError<<" estimate="<<forms.quadratureEstimate<<" kinetic="<<forms.kinetic
                <<" mass_norm="<<forms.mass.norm()<<" reference_compared="<<bool(readPath)<<" entry_error="<<maximum<<" matrix_error="<<relative
                <<" accepted=0 phase_flux_evaluated=0"<<std::endl;
        }
        return;
    }
    if(resolutionOnly){
        // Geometry-only diagnosis of the ORIGINAL B3 configuration, also
        // after accepted steps. The coarse space/tolerance belong to the
        // case, not to the growing number of arcs in a saved source curve.
        // No FV inventory is constructed, replaced or accepted here.
        require(step>0&&begin>=0&&nx==64&&ny==64&&dt==.001&&l==-.5&&b==-.5&&r==.5&&t==.5,
            "free-space replay declares the original B3 projection settings only");
        InterfaceParametricProjection::Grid grid{nx,ny};grid.left=l;grid.bottom=b;grid.width=r-l;grid.height=t-b;
        InterfaceParametricProjection::Options options;options.materialSegments=64;options.relativeH2=.001*std::pow(16./64,2);
        const auto solve=[&](const InterfaceParametricProjection::Source& source,const char* layer,bool keep){
            auto selected=options;if(keep)selected.minimumMaterialKnots=retained;
            const auto started=std::chrono::steady_clock::now();
            std::cout<<std::setprecision(17)<<"MOMENT_FREE_RESOLUTION_BEGIN step="<<step<<" layer="<<layer<<" retained="<<keep<<std::endl;
            try{const auto resolved=InterfaceParametricProjection::resolveMaterialSpace(source,grid,selected);
                std::cout<<std::setprecision(17)<<"MOMENT_FREE_RESOLUTION step="<<step<<" layer="<<layer<<" retained="<<keep
                    <<" rejected=0 source_arcs="<<source.arcs.size()<<" segments="<<resolved.size()-1;}
            catch(const InterfaceParametricProjection::MaterialResolutionLimit& error){std::cout<<"MOMENT_FREE_RESOLUTION step="<<step
                <<" layer="<<layer<<" retained="<<keep<<" rejected=1 reason="<<error.what();}
            std::cout<<" relative_h2_budget="<<selected.relativeH2<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()
                <<" accepted=0 phase_flux_evaluated=0 path="<<path<<std::endl;
        };
        InterfaceParametricProjection::Source stationary;stationary.boundary=curve.boundary;
        for(const auto& arc:curve.arcs){P::PolynomialArc<9> lifted;lifted.a=arc.a;lifted.b=arc.b;
            std::copy(arc.p.begin(),arc.p.end(),lifted.p.begin());stationary.arcs.push_back(lifted);}
        solve(stationary,"accepted",false);
        const auto map=HydroInterfaceTransport::frozenParametricFlow(field);
        for(bool whole:{false,true}){
            InterfaceParametricTransport::SourceReport sourceReport;
            const auto source=InterfaceParametricTransport::transportSource(curve,grid,whole?dt:.5*dt,map,sourceReport);
            solve(source,whole?"full":"half",false);solve(source,whole?"full":"half",true);
        }
        if(const char* output=std::getenv("QUANTUM_MOMENT_PARTITION_OUTPUT")){
            std::ofstream out(output);require(bool(out),"cannot open accepted geometry diagnostic");out<<std::setprecision(17)<<"s,x,y,curvature,arc_width,time\n";
            const auto& rule=P::gaussEight();
            for(const auto& arc:curve.arcs)for(double node:rule.x){const double q=.5*(1+node);const auto sample=P::sampleSurface(arc,q);
                out<<arc.a+(arc.b-arc.a)*q<<','<<sample.point.x<<','<<sample.point.y<<','<<sample.curvature<<','<<arc.b-arc.a<<','<<begin<<'\n';}
            require(bool(out),"cannot write accepted geometry diagnostic");
        }
        return;
    }
    double speed=0,interfaceSpeed=0,normalSpeed=0,interfaceGradient=0,wallNormal=0;
    for(const auto& cell:field.cellSpeedSamples())speed=std::max(speed,cell.fineMaximum);
    for(const auto& arc:curve.arcs)for(int i=0;i<17;i++){
        const auto at=P::sampleSurface(arc,double(i)/16);const auto jet=field.firstJet({at.point.x,at.point.y});
        interfaceSpeed=std::max(interfaceSpeed,std::hypot(jet.value.x,jet.value.y));
        normalSpeed=std::max(normalSpeed,std::abs(jet.value.x*at.normal.x+jet.value.y*at.normal.y));
        interfaceGradient=std::max(interfaceGradient,std::sqrt(jet.dx.x*jet.dx.x+jet.dx.y*jet.dx.y+jet.dy.x*jet.dy.x+jet.dy.y*jet.dy.y));
    }
    for(int i=0;i<=2*std::max(nx,ny);i++){
        const double q=double(i)/(2*std::max(nx,ny)),x=l+(r-l)*q,y=b+(t-b)*q;
        wallNormal=std::max({wallNormal,std::abs(field.firstJet({l,y}).value.x),std::abs(field.firstJet({r,y}).value.x),
            std::abs(field.firstJet({x,b}).value.y),std::abs(field.firstJet({x,t}).value.y)});
    }
    require(wallNormal<1e-10*std::max(1.,speed),"saved candidate does not satisfy the declared free-slip normal data");
    std::cout<<std::setprecision(17)<<"MOMENT_COMPLETE_PHASE_PROBE step="<<step<<" time="<<begin<<" dt="<<dt<<" cells="<<cells<<" coordinates="<<count
        <<" coordinate_norm="<<wind.norm()<<" speed_sample="<<speed<<" interface_speed_sample="<<interfaceSpeed<<" interface_normal_sample="<<normalSpeed
        <<" interface_gradient_sample="<<interfaceGradient<<" dt_interface_gradient="<<dt*interfaceGradient<<" wall_normal="<<wallNormal
        <<" minimum_knots="<<knots<<" minimum_material_width="<<minimum<<" accepted=0 path="<<path<<std::endl;
}
void hydroMomentRZFlowReplay(const char* path){
    namespace P=Quantum::Geometry::Parametric;
    std::ifstream in(path);std::string magic;in>>magic;require(magic=="MOMENT_PHASE_CANDIDATE_V1","invalid moment phase checkpoint");
    int rz,uniform,cells,count,step;double begin,dt,l,b,r,t,hx,hy;
    in>>rz>>uniform>>cells>>count>>step>>begin>>dt>>l>>b>>r>>t>>hx>>hy;
    const int n=int(std::lround(1/hx));require(rz&&uniform&&l==0&&b==0&&r==1&&t==1&&hx==hy&&n*n==cells,"RZ replay expects the saved uniform unit domain");
    Eigen::VectorXd wind(count);for(int i=0;i<count;i++)in>>wind[i];P::Curve curve;std::size_t arcs;bool axisEnds;in>>axisEnds>>arcs;if(axisEnds)curve.boundary=P::CurveBoundary::axisMeridian();
    for(std::size_t i=0;i<arcs;i++){P::Arc arc;in>>arc.a>>arc.b;for(auto& p:arc.p)in>>p.x>>p.y;curve.arcs.push_back(arc);}require(bool(in)&&wind.allFinite(),"incomplete moment phase checkpoint");
    QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
    const auto spline=Spline::build(mesh,domain,space,wind,begin+.5*dt);double maximumSpeed=0,maximumGradient=0,join=0;
    for(int j=0;j<64;j++)for(int i=0;i<64;i++){const auto q=spline.firstJet({(i+.5)/64,(j+.5)/64});
        maximumSpeed=std::max(maximumSpeed,std::hypot(q.value.x,q.value.y));maximumGradient=std::max(maximumGradient,std::max(std::abs(q.dx.x)+std::abs(q.dy.x),std::abs(q.dx.y)+std::abs(q.dy.y)));}
    for(int line=1;line<n;line++)for(int k=0;k<37;k++)for(int axis=0;axis<2;axis++){
        Vector2D a{std::nextafter(double(line)/n,0.),(k+.37)/37},b0{std::nextafter(double(line)/n,1.),a.y};if(axis){std::swap(a.x,a.y);std::swap(b0.x,b0.y);}
        const auto p=spline.firstJet(a),q=spline.firstJet(b0);join=std::max({join,std::hypot(p.value.x-q.value.x,p.value.y-q.value.y),std::hypot(p.dx.x-q.dx.x,p.dx.y-q.dx.y),std::hypot(p.dy.x-q.dy.x,p.dy.y-q.dy.y)});
    }
    // Degree 21 radial / 14 axial physical kinetic integrands: twelve by
    // eight Gauss points integrate this represented polynomial exactly.
    double rx[12],rw[12];const double pi=std::acos(-1.);const auto& gy=P::gaussEight();
    kineticQuadrature(rx,rw);
    double represented=0;
    std::vector<double> cellKinetic(mesh.cellIndexNum(),0);
    for(int j=0;j<n;j++)for(int i=0;i<n;i++)for(int a=0;a<12;a++)for(int q=0;q<8;q++){
        const Vector2D p{(i+.5+.5*rx[a])/n,(j+.5+.5*gy.x[q])/n};const auto jet=spline.firstJet(p);const auto u=jet.value;
        const double contribution=.25*hx*hy*rw[a]*gy.w[q]*2*pi*p.x*(u.x*u.x+u.y*u.y);
        represented+=contribution;cellKinetic[jet.nativeCell]+=contribution;}
    const double moment=wind.dot(space.mass()*wind);
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_REPLAY_FIELD step="<<step<<" dt="<<dt<<" wind_norm="<<wind.norm()<<" sampled_speed="<<maximumSpeed<<" sampled_gradient="<<maximumGradient
        <<" sampled_dt_gradient="<<dt*maximumGradient<<" C1_join="<<join<<" reconstructed_mass_norm="<<represented<<" owned_mass_norm="<<moment<<" kinetic_ratio="<<represented/moment<<std::endl;
    // Separate the C1 extension from the pre-existing local minimum-L2
    // solenoidal reconstruction. They share owned moments but have different
    // trace spaces, so this comparison is diagnostic, not a common norm proof.
    using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;
    const auto transport=Transport::build(mesh,domain,space);
    double localKinetic=0,affineKinetic=0,largestRatio=0;int largestCell=-1;
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& block=space.cells()[c];Eigen::VectorXd local(block.rows.size());
        for(int i=0;i<local.size();i++)local[i]=wind[block.rows[i]];
        const double owned=local.dot(block.mass*local);double localEnergy=0,affineEnergy=0;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto center=cell.centerCoord(),size=cell.size();
        for(const auto& piece:transport.cells()[c].pieces){const Eigen::VectorXd coefficients=piece.coefficients*local;
            const int panels=piece.axisRegular?1:4;
            for(int panel=0;panel<panels;panel++)for(int j=0;j<8;j++)for(int i=0;i<8;i++){
                const Vector2D q{-.5+(panel+.5+.5*gy.x[i])/panels,.5*gy.x[j]},p=piece.center+Vector2D{piece.size.x*q.x,piece.size.y*q.y};
                const double weight=.25*piece.size.x*piece.size.y*gy.w[i]*gy.w[j]*2*pi*p.x/panels;
                const Eigen::Vector2d v=Transport::physicalBasis(piece,q)*coefficients;
                const Eigen::Vector2d a=Space::affineVelocityBasis({(p.x-center.x)/size.x,(p.y-center.y)/size.y},size)*block.affineProjection*local;
                localEnergy+=weight*v.squaredNorm();affineEnergy+=weight*a.squaredNorm();
            }}
        localKinetic+=localEnergy;affineKinetic+=affineEnergy;
        const double ratio=cellKinetic[c]/owned;if(ratio>largestRatio){largestRatio=ratio;largestCell=c;}
        if(cellKinetic[c]>.05*represented)std::cout<<std::setprecision(17)<<"MOMENT_RZ_REPLAY_CELL cell="<<c<<" radius="<<center.x<<" z="<<center.y
            <<" C1_energy="<<cellKinetic[c]<<" local_minimum_energy="<<localEnergy<<" affine_energy="<<affineEnergy<<" native_energy="<<owned<<std::endl;
    }
    std::cout<<std::setprecision(17)<<"MOMENT_RZ_REPLAY_NORMS local_minimum="<<localKinetic<<" affine="<<affineKinetic<<" local_to_native="<<localKinetic/moment
        <<" C1_to_local="<<represented/localKinetic<<" largest_cell_ratio="<<largestRatio<<" largest_cell="<<largestCell<<std::endl;
    const double parameters[]={.37334809826284554,.37442453928687264,.37467823127896349,.37540694182951778,.62366427949646241,.62668477494472008};
    for(int blocks:{8,16,32})for(double duration:{.5*dt,dt}){
        const auto map=HydroInterfaceTransport::frozenParametricFlow(spline,blocks);double jump=0,finiteDifference=0,metric=0,mapDifference=0;int rejected=0;
        const auto finer=HydroInterfaceTransport::frozenParametricFlow(spline,2*blocks);
        for(double s:parameters)try{
            const auto p=curve.sample(s),pnext=curve.sample(std::nextafter(s,1.)),d=curve.sample(s,1);const auto a=map({p.x,p.y},duration),b0=map({pnext.x,pnext.y},duration),reference=finer({p.x,p.y},duration);
            jump=std::max({jump,std::hypot(a.dx.x-b0.dx.x,a.dx.y-b0.dx.y),std::hypot(a.dy.x-b0.dy.x,a.dy.y-b0.dy.y)});
            const auto left=curve.sample(s-1e-7),right=curve.sample(s+1e-7);const auto lm=map({left.x,left.y},duration),rm=map({right.x,right.y},duration);
            const Vector2D derivative{a.dx.x*d.x+a.dy.x*d.y,a.dx.y*d.x+a.dy.y*d.y};const auto fd=(rm.point-lm.point)*(1/(2e-7));
            finiteDifference=std::max(finiteDifference,std::hypot(fd.x-derivative.x,fd.y-derivative.y)/std::max(1.,derivative.mod()));
            metric=std::max(metric,std::abs(a.point.x/p.x*(a.dx.x*a.dy.y-a.dx.y*a.dy.x)-1));
            mapDifference=std::max({mapDifference,(a.point-reference.point).mod(),(a.dx-reference.dx).mod(),(a.dy-reference.dy).mod()});
        }catch(const std::exception& error){rejected++;std::cout<<"MOMENT_RZ_REPLAY_REJECTION blocks="<<blocks<<" duration="<<duration<<" reason="<<error.what()<<std::endl;}
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_REPLAY_MAP blocks="<<blocks<<" duration="<<duration<<" neighboring_jet_jump="<<jump<<" finite_difference_error="<<finiteDifference
            <<" weighted_jacobian_error="<<metric<<" doubled_blocks_difference="<<mapDifference<<" rejected="<<rejected<<std::endl;
    }
}
void hydroMomentRZKineticContract(){
    namespace P=Quantum::Geometry::Parametric;
    for(int n:{4,8}){
        QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,true);const auto space=Space::build(mesh,domain);
        const int nv=int(space.mass().rows()),nd=int(space.divergence().rows());std::vector<int> walls;
        for(const auto face:domain.faces(mesh))if(face.outwardSign&&space.faces()[face.nativeRow].first>=0){walls.push_back(space.faces()[face.nativeRow].first);walls.push_back(walls.back()+1);}
        Eigen::MatrixXd constraints=Eigen::MatrixXd::Zero(nd+walls.size(),nv);constraints.topRows(nd)=space.divergence();
        for(int i=0;i<nd;i++)constraints.row(i)/=constraints.row(i).norm();
        for(int i=0;i<int(walls.size());i++)constraints(nd+i,walls[i])=1;
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(constraints,Eigen::ComputeFullV);int rank=0;
        for(int i=0;i<svd.singularValues().size();i++)if(svd.singularValues()[i]>1e-12*svd.singularValues()[0])rank++;
        Eigen::MatrixXd basis=svd.matrixV().rightCols(nv-rank);
        // The wall is an exact essential trace, not an SVD-roundoff velocity.
        for(int row:walls)basis.row(row).setZero();
        const double constraintError=(constraints*basis).norm();require(constraintError<2e-11,"RZ kinetic nullspace violates owned constraints");
        const auto base=Spline::build(mesh,domain,space,Eigen::VectorXd::Zero(nv),0);
        double rx[12],rw[12];kineticQuadrature(rx,rw);const auto& gy=P::gaussEight();const double pi=std::acos(-1.);
        Eigen::MatrixXd values(2*n*n*12*8,basis.cols());
        for(int column=0;column<basis.cols();column++){
            const auto field=base.withVelocity(mesh,domain,space,basis.col(column),0);int row=0;
            for(int j=0;j<n;j++)for(int i=0;i<n;i++)for(int a=0;a<12;a++)for(int b=0;b<8;b++){
                const Vector2D p{(i+.5+.5*rx[a])/n,(j+.5+.5*gy.x[b])/n};const auto u=field.firstJet(p).value;
                const double root=std::sqrt(.25/(n*n)*rw[a]*gy.w[b]*2*pi*p.x);
                values(row++,column)=root*u.x;values(row++,column)=root*u.y;
            }
        }
        const Eigen::MatrixXd actual=values.transpose()*values,owned=basis.transpose()*space.mass()*basis;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigen(actual,owned,Eigen::EigenvaluesOnly);
        require(eigen.info()==Eigen::Success&&eigen.eigenvalues().allFinite()&&eigen.eigenvalues()[0]>0,"RZ common velocity kinetic norm lost injectivity");
        std::cout<<std::setprecision(17)<<"MOMENT_RZ_KINETIC_SPECTRUM n="<<n<<" solenoidal_dofs="<<basis.cols()<<" minimum="<<eigen.eigenvalues()[0]
            <<" maximum="<<eigen.eigenvalues()[basis.cols()-1]<<" constraints="<<constraintError<<" affine_consistency="<<space.report().maximumConsistencyDefect
            <<" hodge_diagonal_scaled_min="<<space.report().minimumScaledEigenvalue<<" hodge_diagonal_scaled_max="<<space.report().maximumScaledEigenvalue<<std::endl;
    }
}
void hydroMomentMaterialSpaceContract(){
    constexpr int n=8;constexpr double dt=.02;
    QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const auto velocity=space.interpolateVelocity(mesh,domain,[](Vector2D p){const auto v=exact(p).value;return Eigen::Vector2d(v.x,v.y);});
    const auto spline=Spline::build(mesh,domain,space,velocity,0);const auto map=HydroInterfaceTransport::frozenParametricFlow(spline);const auto flux=spline.normalFlux();
    auto old=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.43,.57},{.13,.13},24));
    const InterfaceStepTime time{1,0,dt,dt};InterfaceParametricProjection::Options options;
    std::atomic<std::size_t> queries{0};const InterfaceParametricTransport::Flow counted=[&](Vector2D p,double t){queries.fetch_add(1,std::memory_order_relaxed);return map(p,t);};
    auto prepared=InterfaceParametricMotion::prepare(mesh,domain,old,time,dt,counted,flux,options);
    const auto before=prepared.state.parametricGeometry().geometryHandle();const auto inventory=prepared.state.parametricGeometry().inventory();const auto volumes=prepared.state.parametricFluxes();
    std::vector<double> knots;for(const auto& arc:before->arcs){knots.push_back(arc.a);knots.push_back(.5*(arc.a+arc.b));}knots.push_back(1);
    const auto count=queries.load();const auto start=std::chrono::steady_clock::now();
    InterfaceParametricMotion::refineMaterialSpace(prepared,knots);
    const double reusedSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    require(queries.load()==count,"material-space refinement retraced a frozen trajectory");
    require(prepared.state.parametricGeometry().inventory()==inventory&&prepared.state.parametricFollows(old,time),"material refinement changed inventory or lineage");
    options.minimumMaterialKnots=knots;const auto freshStart=std::chrono::steady_clock::now();
    auto fresh=InterfaceParametricMotion::prepare(mesh,domain,old,time,dt,counted,flux,options);
    const double freshSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-freshStart).count();
    const auto& a=prepared.state.parametricGeometry().geometry();const auto& b=fresh.state.parametricGeometry().geometry();
    require(a.arcs.size()==b.arcs.size()&&a.arcs.size()>=2*before->arcs.size(),"material refinement lost its nested partition");
    for(std::size_t i=0;i<a.arcs.size();i++){
        require(a.arcs[i].a==b.arcs[i].a&&a.arcs[i].b==b.arcs[i].b,"reused material partition differs from fresh preparation");
        for(int k=0;k<5;k++)require(a.arcs[i].p[k].x==b.arcs[i].p[k].x&&a.arcs[i].p[k].y==b.arcs[i].p[k].y,"reused reconstruction differs from fresh geometry coefficients");
    }
    require(fresh.state.parametricGeometry().inventory()==inventory,"fresh motion has a different independent inventory");
    const auto& reusedFlux=prepared.state.parametricFluxes();const auto& freshFlux=fresh.state.parametricFluxes();
    for(int direction=0;direction<2;direction++)for(std::size_t face=0;face<volumes.faces[direction].size();face++){
        const auto& oldFace=volumes.faces[direction][face];
        for(const auto* actual:{&reusedFlux.faces[direction][face],&freshFlux.faces[direction][face]})
            require(actual->total==oldFace.total&&actual->liquid==oldFace.liquid&&actual->totalPositive==oldFace.totalPositive&&actual->totalNegative==oldFace.totalNegative
                &&actual->liquidPositive==oldFace.liquidPositive&&actual->liquidNegative==oldFace.liquidNegative&&actual->totalHigher==oldFace.totalHigher&&actual->liquidHigher==oldFace.liquidHigher,"material refinement changed shared phase flux moments");
    }
    const auto acceptedGeometry=prepared.state.parametricGeometry().geometryHandle();bool invalid=false,moved=false;
    try{InterfaceParametricMotion::refineMaterialSpace(prepared,{std::numeric_limits<double>::quiet_NaN()});}catch(const std::invalid_argument&){invalid=true;}
    require(invalid&&prepared.state.parametricGeometry().geometryHandle()==acceptedGeometry&&prepared.state.parametricGeometry().inventory()==inventory,"failed refinement changed its prepared state");
    auto accepted=std::move(prepared.state);
    try{InterfaceParametricMotion::refineMaterialSpace(prepared,knots);}catch(const std::invalid_argument&){moved=true;}
    require(moved&&accepted.parametricGeometry().geometryHandle()==acceptedGeometry,"refinement accepted a consumed motion result");
    std::cout<<std::setprecision(17)<<"MOMENT_MATERIAL_SPACE source_queries="<<count<<" extra_refine_queries=0 refined_arcs="<<a.arcs.size()
        <<" coefficients_bitwise=1 inventory_bitwise=1 phase_moments_bitwise=1 refined_seconds="<<reusedSeconds<<" fresh_seconds="<<freshSeconds<<" invalid_rejected="<<invalid<<" consumed_rejected="<<moved<<std::endl;
}
void hydroMomentCapillaryContract(){
    using Traction=Quantum::Hydrodynamics::HydroCompatibleMomentTraction;using Material=Quantum::Hydrodynamics::HydroCompatibleMomentMaterial;
    namespace P=Quantum::Geometry::Parametric;
    QuadtreeMesh mesh(input(8));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const auto wind=space.interpolateVelocity(mesh,domain,[](Vector2D p){const auto v=exact(p).value;return Eigen::Vector2d(v.x,v.y);});
    const auto spline=Spline::build(mesh,domain,space,wind,0);
    auto state=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.43,.57},{.23,.17},24));
    const auto material=Material::buildFrozen(mesh,domain,space,state.parametricGeometry().geometryHandle(),0,{1,.1,.01,.001});
    const auto constant=Traction::build(mesh,domain,space,spline,material,[](const P::SurfaceSample&){return Traction::Components{3,0};});
    const auto load=Traction::build(mesh,domain,space,spline,material,[](const P::SurfaceSample& q){return Traction::Components{q.curvature,0};},1);
    long double work=0;const auto& rule=P::gaussEight();
    for(const auto& arc:state.parametricGeometry().geometry().arcs)for(int part=0;part<128;part++)for(int k=0;k<8;k++){
        const double t=(part+.5+.5*rule.x[k])/128;const auto q=P::sampleSurface(arc,t);const auto u=spline.firstJet({q.point.x,q.point.y}).value;
        work+=(.5L*rule.w[k]/128)*(q.curvature-load.report.meanNormal)*(-q.derivative.y*u.x+q.derivative.x*u.y);
    }
    const double virtualWork=std::abs(wind.dot(load.variation)-double(work)),constantVariation=constant.variation.norm();
    require(virtualWork<2e-10&&constantVariation<2e-11,"moment sharp traction differs from independent surface work or constant pressure range");
    std::cout<<std::setprecision(17)<<"MOMENT_CAPILLARY_WORK defect="<<virtualWork<<" constant_variation="<<constantVariation<<" pressure_work="<<wind.dot(space.pressureLoad()*constant.pressurePotential)<<std::endl;
    double previousSpeed=0;
    for(int segments:{16,32}){
        auto circle=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.43,.57},{.2,.2},segments));
        auto flow=Flow::createTwoPhase(mesh,domain,space,std::move(circle),{1,.1,.01,.001},Eigen::VectorXd::Zero(wind.size()),1);
        const HydroStepTime time{1,0,.001};const auto zero=flow.makeLoad(mesh,domain,Eigen::VectorXd::Zero(wind.size()),.0005);Flow::Options options;
        options.interfaceProjection.materialSegments=segments;options.interfaceProjection.relativeH2=.001*std::pow(16./segments,2);
        auto candidate=flow.prepareNavierStokesStep(mesh,domain,time,zero,options);const auto report=candidate.report();flow.accept(mesh,domain,std::move(candidate));
        const auto velocity=*flow.state().physicalVelocity;double speed=0;
        for(int j=0;j<40;j++)for(int i=0;i<40;i++){const auto v=velocity.firstJet({(i+.5)/40,(j+.5)/40}).value;speed=std::max(speed,std::hypot(v.x,v.y));}
        const auto surface=Traction::measure(flow.interfaceState().parametricGeometry().geometry());
        const double curvature=std::sqrt(std::max(0.,surface.curvatureSquaredIntegral/surface.length-10*surface.meanCurvature+25))/5;
        std::cout<<std::setprecision(17)<<"MOMENT_STATIC_DROP segments="<<segments<<" time="<<flow.state().time<<" speed="<<speed<<" curvature_relative_rms="<<curvature
            <<" pressure_mean_jump="<<surface.meanCurvature<<" iterations="<<report.nonlinearIterations<<" equation="<<report.equationResidual<<" kinetic_budget="<<report.energyDefect<<" surface_work_defect="<<report.surfaceWorkDefect<<std::endl;
        require(speed<.001&&std::abs(surface.meanCurvature-5)<.002&&report.equationResidual<1e-10,"numerical-curvature static drop failed");
        if(previousSpeed>0)require(speed<.3*previousSpeed,"static-drop spurious velocity did not decrease with interface refinement");previousSpeed=speed;
    }
}
void hydroMomentTwoPhaseContract(){
    constexpr int n=8,steps=6;constexpr double dt=.02;
    QuadtreeMesh mesh(input(n));FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);const auto space=Space::build(mesh,domain);
    const auto initial=space.interpolateVelocity(mesh,domain,[](Vector2D p){const auto v=exact(p).value;return Eigen::Vector2d(v.x,v.y);});
    auto interface=InterfaceParametricMotion::initialize(mesh,domain,InterfacePrescribedData::parametricEllipse({.43,.57},{.13,.13},24));
    const auto first=interface.parametricGeometry().geometry();const double initialVolume=first.area();
    auto hydro=Flow::createTwoPhase(mesh,domain,space,std::move(interface),{5,1,.02,.01},initial);
    const Space::Matrix firstMass=hydro.kineticMass();const double firstEnergy=hydro.kineticEnergy();
    std::ofstream output,history;if(const char* prefix=std::getenv("QUANTUM_MOMENT_TWO_PHASE_OUTPUT")){
        output.open(std::string(prefix)+"_curves.csv");output<<std::setprecision(17)<<"step,time,parameter,x,y\n";
        history.open(std::string(prefix)+"_history.csv");history<<std::setprecision(17)<<"step,time,volume_drift,fraction_error,kinetic_energy,equation_residual,temporal_mass_work,material_metric_work,energy_defect\n";
    }
    const auto save=[&](int step){if(output)for(int k=0;k<=512;k++){const double s=double(k)/512;const auto p=hydro.interfaceState().parametricGeometry().geometry().sample(s);output<<step<<','<<hydro.state().time<<','<<s<<','<<p.x<<','<<p.y<<'\n';}};save(0);
    double maxVolume=0,maxFraction=0,maxEquation=0,maxEnergy=0,maxMass=0,maxMetric=0,maxLocal=0,temporalWork=0,metricWork=0,viscousWork=0;
    int iterations=0;bool rollback=false,duplicate=false;
    for(int step=1;step<=steps;step++){
        HydroStepTime time;time.index=step;time.begin=hydro.state().time;time.dt=dt;
        const auto oldGeometry=hydro.interfaceState().parametricGeometry().geometryHandle();const auto oldVelocity=hydro.state().velocity;
        const auto oldMass=hydro.kineticMass();const auto oldTrace=hydro.activeTangentialRows();
        const auto load=hydro.makeLoad(mesh,domain,Eigen::VectorXd::Zero(initial.size()),time.begin+.5*dt);
        if(step==1){Flow::Options fail;fail.maximumIterations=1;
            try{auto rejected=hydro.prepareNavierStokesStep(mesh,domain,time,load,fail);}catch(const std::invalid_argument&){rollback=true;}
            require(rollback&&hydro.state().acceptedSteps==0&&hydro.interfaceState().parametricGeometry().geometryHandle()==oldGeometry
                &&hydro.state().velocity==oldVelocity&&(hydro.kineticMass()-oldMass).norm()==0&&hydro.activeTangentialRows()==oldTrace,"moving midpoint rejection changed accepted state");}
        auto candidate=hydro.prepareNavierStokesStep(mesh,domain,time,load);const auto report=candidate.report();
        require(candidate.interfaceState().parametricFollows(hydro.interfaceState(),{step,time.begin,time.begin+dt,dt}),"moving midpoint candidate lost interface lineage");
        require(hydro.interfaceState().parametricGeometry().geometryHandle()==oldGeometry&&hydro.state().velocity==oldVelocity,"moving midpoint preparation mutated accepted data");
        hydro.accept(mesh,domain,std::move(candidate));
        if(step==1){try{hydro.accept(mesh,domain,std::move(candidate));}catch(const std::invalid_argument&){duplicate=true;}require(duplicate,"moving midpoint accepted a consumed candidate");}
        const auto& geometry=hydro.interfaceState().parametricGeometry();const auto geometric=geometry.geometricVolumes();double volume=0,fraction=0;
        for(int k=0;k<int(geometry.inventory().size());k++){volume+=geometry.inventory()[k];fraction=std::max(fraction,std::abs(geometric[k]-geometry.inventory()[k])*n*n);}
        maxVolume=std::max(maxVolume,std::abs(volume-initialVolume));maxFraction=std::max(maxFraction,fraction);maxEquation=std::max(maxEquation,report.equationResidual);
        maxEnergy=std::max(maxEnergy,std::abs(report.energyDefect)/firstEnergy);maxMass=std::max(maxMass,report.phaseMassClosure);maxMetric=std::max(maxMetric,report.massMetricDefect);maxLocal=std::max(maxLocal,report.localMomentumDefect);
        temporalWork+=report.temporalMassWork;metricWork+=report.materialMetricWork;viscousWork+=report.viscousDissipation;iterations+=report.nonlinearIterations;
        require(report.transportDissipation==0&&report.viscousDissipation>=0,"moving midpoint added hidden dissipation or negative viscosity work");
        if(history)history<<step<<','<<hydro.state().time<<','<<volume-initialVolume<<','<<fraction<<','<<report.kineticAfter<<','<<report.equationResidual<<','<<report.temporalMassWork<<','<<report.materialMetricWork<<','<<report.energyDefect<<'\n';
        if(step==3||step==steps)save(step);
        std::cout<<std::setprecision(17)<<"MOMENT_TWO_PHASE_STEP step="<<step<<" time="<<hydro.state().time<<" rho_l=5 rho_g=1 iterations="<<report.nonlinearIterations
            <<" energy="<<report.kineticAfter<<" equation="<<report.equationResidual<<" divergence="<<report.divergenceResidual<<" trace="<<report.viscousTraceResidual
            <<" mass_closure="<<report.phaseMassClosure<<" geometric_closure="<<report.geometricMassClosure<<" metric_defect="<<report.massMetricDefect
            <<" temporal_work="<<report.temporalMassWork<<" material_work="<<report.materialMetricWork<<" energy_defect="<<report.energyDefect<<std::endl;
    }
    const double massChange=(hydro.kineticMass()-firstMass).norm()/firstMass.norm();double motion=0;
    for(int k=0;k<512;k++){const auto a=first.sample(double(k)/512),b=hydro.interfaceState().parametricGeometry().geometry().sample(double(k)/512);motion=std::max(motion,std::hypot(a.x-b.x,a.y-b.y));}
    const double finalEnergy=hydro.kineticEnergy();
    const double balance=(finalEnergy-firstEnergy+viscousWork-temporalWork+metricWork)/firstEnergy;
    std::cout<<std::setprecision(17)<<"MOMENT_TWO_PHASE_COMPLETE steps="<<steps<<" time="<<hydro.state().time<<" iterations="<<iterations<<" volume_drift="<<maxVolume<<" fraction_error="<<maxFraction
        <<" equation="<<maxEquation<<" relative_energy_defect="<<maxEnergy<<" mass_closure="<<maxMass<<" metric_defect="<<maxMetric<<" local_momentum="<<maxLocal<<" mass_matrix_change="<<massChange
        <<" material_displacement="<<motion<<" accumulated_energy_balance="<<balance<<" temporal_work="<<temporalWork<<" material_work="<<metricWork
        <<" rollback="<<rollback<<" duplicate_rejected="<<duplicate<<" variable_mass=1 capillary_feedback=0\n";
    require(maxVolume<2e-13&&maxFraction<1e-10&&maxMass<1e-12&&maxMetric<2e-10&&maxEnergy<2e-10&&std::abs(balance)<1e-9&&maxLocal<1e-10&&massChange>.001&&motion>.005,
        "moving two-phase midpoint trajectory or physical budgets failed");
}
