#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentTransport.h"
#include "Modules/Hydrodynamics/Hydrodynamics.h"
#include "Modules/Hydrodynamics/HydroSplineVelocity.h"
#include "Modules/Hydrodynamics/HydroInterfaceTransport.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <numeric>

namespace {
void require(bool valid,const char* why){if(!valid)throw std::runtime_error(why);}
// u_x=a+b*y, u_y=0. Polynomial coefficients are the actual native primitive,
// not an interpolation fitted to queries of an unknown field.
InterfaceNormalFluxView affineFlux(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,double a,double b,bool products=false) {
    using namespace Quantum::Geometry;
    InterfaceNormalFluxView::Queries q;q.retainedMomentDegree=2;q.retainedVelocityProducts=products;
    q.chart=[&mesh,&domain](int row){const auto face=domain.face(mesh,row);const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto lo=mesh.verts(edge.edgeToVert[0]).coord(),hi=mesh.verts(edge.edgeToVert[1]).coord();return InterfaceNormalFluxView::FaceChart{lo,(lo+hi)*.5,hi};};
    q.polynomial=[chart=q.chart,&mesh,&domain,a,b](int row,int half){UnitPolynomialTerm term;const auto c=chart(row);
        const double lower=half?c.midpoint.y:c.begin.y,width=domain.face(mesh,row).horizontal?0:(c.end.y-c.begin.y)*.5;
        term.factor=width;term.coefficients[0]=a+b*lower;term.coefficients[1]=b*width;return std::vector<UnitPolynomialTerm>{term};};
    q.integral=[polynomial=q.polynomial](int row,int half,double a0,double a1,double s0,double s1){return integrateUnitPolynomial(polynomial(row,half),a0,a1,s0,s1);};
    q.density=[polynomial=q.polynomial](int row,int half,bool end,long double offset){const auto p=polynomial(row,half)[0];return p.factor*(p.coefficients[0]+p.coefficients[1]*(end?1-offset:offset));};
    q.phaseIntegral=[](int,int,const LevelSetMappedSegment&)->long double{throw std::logic_error("unused mapped scalar integral in parameter control");};
    return InterfaceNormalFluxView(std::move(q));
}
void phaseMass(bool reversal,bool products=false) {
    using namespace Quantum::Geometry::Parametric;
    using namespace InterfaceParametricTransport;
    using Quantum::Hydrodynamics::HydroConservativeTransport;
    constexpr int n=8;const double dt=reversal?.2:.1;
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);
    const auto initial=InterfacePrescribedData::parametricEllipse({.45,.5},{.15,.15},16);
    auto old=InterfaceParametricMotion::initialize(mesh,domain,initial);const auto inventory=old.parametricGeometry().inventory();
    const double a=reversal?.7:-.7*.53,b=reversal?0:.7;
    const auto first=affineFlux(mesh,domain,a,b,products),last=affineFlux(mesh,domain,reversal?-a:a,b,products);
    const Flow flow=[=](Vector2D p,double time){const double factor=reversal?time*(1-time/dt):time;
        const Vector2D move{(a+b*p.y)*factor,0};return MapJet{p+move,{1,0},{b*factor,1},{(a+b*p.y)*(reversal?1-2*time/dt:1),0},move};};
    const InterfaceStepTime time{1,0,dt,dt};InterfaceParametricProjection::Options options;
    auto next=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,old,time,flow,first,last,options);
    require(old.progress().completedSteps==0 && old.parametricGeometry().inventory()==inventory && next.state.parametricFollows(old,time),
        "parameter mass publication changed the source or lost lineage");
    const auto& flux=next.state.parametricFluxes();require(flux.hasDirections,"parameter motion did not publish directional volumes");
    auto final=initial;if(!reversal)for(auto& arc:final.arcs){arc.p[0].x+=dt*(a+b*arc.p[0].y);for(int k=1;k<5;++k)arc.p[k].x+=dt*b*arc.p[k].y;}
    PhaseMeasure before(initial),after(final);double error=0,netError=0,gross=0,balance=0,momentError=0;
    require(flux.momentDegree==2,"actual phase update did not retain its momentum flux moments");
    const auto area=[&](PhaseMeasure& phase,double x,double lo,double hi){return x>0 && hi>lo?phase.simpleRectangle(0,x,lo,hi):0.;};
    for(int line=0;line<=n;++line)for(int strip=0;strip<n;++strip){const auto& f=flux.faces[0][line*n+strip];
        const double x=double(line)/n,lo=double(strip)/n,hi=double(strip+1)/n;double positive=0,negative=0;
        if(reversal){const double move=a*dt/4;if(x>=move)positive=before.simpleRectangle(x-move,x,lo,hi);negative=positive;}
        else {const double split=std::clamp(.53,lo,hi);
            positive=area(before,x,split,hi)-area(after,x,split,hi);
            negative=area(after,x,lo,split)-area(before,x,lo,split);}
        error=std::max({error,std::abs(f.liquidPositive-positive),std::abs(f.liquidNegative-negative)});
        netError=std::max(netError,std::abs(f.liquid-(positive-negative)));gross+=f.liquidPositive+f.liquidNegative;
        // y is unchanged by this affine shear. Its first/second material
        // moments in x<face give an independent space-time flux reference.
        if(x>0){const auto firstMoments=before.rectangleMoments(0,x,lo,hi),lastMoments=after.rectangleMoments(0,x,lo,hi);
            const double m0=firstMoments[0]-lastMoments[0],m1=firstMoments[2]-lastMoments[2],m2=firstMoments[5]-lastMoments[5];
            momentError=std::max({momentError,std::abs(f.liquidHigher[0]-(m1+.5*m0)),std::abs(f.liquidHigher[1]-(m2+m1+.25*m0))});}
        for(int power=1;power<=2;power++){const double expected=reversal?0:dt*(hi-lo)*((a+b*lo)/(power+1)+b*(hi-lo)/(power+2));
            momentError=std::max(momentError,std::abs(f.totalHigher[power-1]-expected));}
        require(f.liquidPositive>=0 && f.liquidPositive<=f.totalPositive && f.liquidNegative>=0 && f.liquidNegative<=f.totalNegative,
            "published direction is not a phase subvolume");
    }
    for(std::size_t k=0;k<inventory.size();++k)balance+=next.state.parametricGeometry().inventory()[k]-inventory[k];
    require(error<3e-13 && netError<3e-13 && momentError<3e-13 && gross>1e-4 && std::abs(balance)<3e-15,
        "gross phase volumes disagree with independent swept/affinely transformed geometry");
    const auto mass=HydroConservativeTransport::build(mesh,domain,old,next.state,time,1000,1);
    const auto space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace::build(mesh,domain);
    const auto momentMass=Quantum::Hydrodynamics::HydroCompatibleMomentTransport::phaseMassFlux(mesh,domain,space,old,next.state,time,1000,1);
    if(products){
        constexpr int degree=FaceVolumes::productDegree;double totalError=0,directionError=0,phaseError=0;
        const auto& g=gaussEight();
        for(int line=0;line<=n;line++)for(int strip=0;strip<n;strip++){
            const auto& f=flux.faces[0][line*n+strip];const double lo=double(strip)/n,hi=double(strip+1)/n;
            for(int k=0;k<=degree;k++){
                const double exact=reversal?0:dt*(hi-lo)*((a+b*lo)/(degree+1)+b*(hi-lo)*(k+1)/((degree+1)*(degree+2)));
                totalError=std::max(totalError,std::abs(f.totalProducts[k]-exact));
                if(reversal&&(line==0||line==n))directionError=std::max(directionError,std::abs(f.totalPositiveProducts[k]-dt*(hi-lo)*a/(4*(degree+1))));
            }
            if(line!=n/2||(strip!=n/2&&strip!=n/2-1))continue;
            for(int k:{0,3,10,17,degree}){
                const auto integrand=[&](double t){const double y=lo+(hi-lo)*t,x=double(line)/n;
                    const auto width=[&](PhaseMeasure& phase){double value=0;for(const auto& ab:phase.intervals(y))value+=std::max(0.,std::min(x,ab[1])-std::max(0.,ab[0]));return value;};
                    return (hi-lo)*FaceVolumes::productWeight(k,t,false,false,lo,hi)*(width(before)-width(after));};
                const auto rule=[&](double left,double right){double result=0;for(int q=0;q<8;q++)result+=.5*(right-left)*g.w[q]*integrand(left+.5*(right-left)*(1+g.x[q]));return result;};
                std::function<double(double,double,double,int)> integral;
                integral=[&](double left,double right,double budget,int depth){const double m=.5*(left+right),coarse=rule(left,right),fine=rule(left,m)+rule(m,right);
                    if(std::abs(fine-coarse)<=budget)return fine;require(depth<24,"independent weighted swept-area integral unresolved");return integral(left,m,.5*budget,depth+1)+integral(m,right,.5*budget,depth+1);};
                phaseError=std::max(phaseError,std::abs(f.liquidProducts[k]-integral(0,1,2e-15,0)));
            }
        }
        std::cout<<std::setprecision(17)<<"MOMENT_WEIGHTED_PHASE mode="<<(reversal?"temporal_reversal":"spatial_counterflow")<<" total="<<totalError<<" exterior_direction="<<directionError
            <<" swept_phase="<<phaseError<<" shared_constant="<<momentMass.productConstantDefect<<" mass_closure="<<momentMass.maximumRelativeClosure<<std::endl;
        require(flux.hasVelocityProducts&&totalError<3e-13&&directionError<3e-13&&phaseError<3e-13&&momentMass.productConstantDefect<2e-12,
            "weighted flux lost exact total, directional flow, swept phase geometry or its original inventory functional");
    }
    auto unrelatedMoment=InterfaceParametricMotion::initialize(mesh,domain,initial);bool momentRejected=false;
    try{Quantum::Hydrodynamics::HydroCompatibleMomentTransport::phaseMassFlux(mesh,domain,space,unrelatedMoment,next.state,time,1000,1);}catch(const std::invalid_argument&){momentRejected=true;}
    require(momentRejected&&momentMass.maximumRelativeClosure<1e-13,"moment phase flux lost source identity or inventory mass closure");
    std::vector<double> velocity(domain.nativeFaceCount(),.23);
    const auto momentum=mass.upwind(velocity,[](const auto&,int){return .23;});double constantError=0;
    for(std::size_t k=0;k<velocity.size();++k)constantError=std::max(constantError,
        static_cast<double>(std::abs(.23L*(mass.massAfter()[k]-mass.massBefore()[k])+momentum.outward[k])));
    auto unrelated=InterfaceParametricMotion::initialize(mesh,domain,initial);bool rejected=false;
    try{HydroConservativeTransport::build(mesh,domain,unrelated,next.state,time,1000,1);}catch(const std::exception&){rejected=true;}
    require(rejected && constantError<1e-11,"mass bridge accepted an unrelated source or lost constant momentum");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_MASS mode="<<(reversal?"temporal_reversal":"spatial_counterflow")
        <<" directional_volume_error="<<error<<" net_error="<<netError<<" liquid_gross="<<gross<<" inventory_balance="<<balance
        <<" dual_mass_relative_closure="<<mass.maxRelativeClosure()<<" constant_momentum_error="<<constantError
        <<" phase_moment_error="<<momentError<<" moment_mass_closure="<<momentMass.maximumRelativeClosure
        <<" directional_budget_adjustment="<<next.flux.directionalBoundAdjustment<<" wrong_source_rejected="<<rejected<<'\n';
}
void solvedMassSteps() {
    using Quantum::Hydrodynamics::HydroConservativeTransport;
    constexpr int n=8;constexpr double dt=.005;
    DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=n;input.maxCellNum=1024;input.maxInitialIterations=0;
    QuadtreeMesh mesh(input);FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,false);
    VEdgeField<double> ux(mesh,0);HEdgeField<double> uy(mesh,0);const double pi=std::acos(-1.);
    const auto stream=[&](Vector2D p){return .1*std::sin(pi*p.x)*std::sin(pi*p.y)/pi;};
    for(const auto face:domain.faces(mesh)){const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto lo=mesh.verts(edge.edgeToVert[0]).coord(),hi=mesh.verts(edge.edgeToVert[1]).coord();
        (face.horizontal?uy.value:ux.value)[edge.index()]=face.outwardSign?0:(face.horizontal?-1:1)*(stream(hi)-stream(lo))/face.area;}
    const auto initial=InterfacePrescribedData::parametricEllipse({.45,.57},{.15,.15},16);
    auto state=InterfaceParametricMotion::initialize(mesh,domain,initial);
    const auto kind=HydroStepGeometry::complete(false);Hydrodynamics hydro(mesh);hydro.initializeState(mesh,kind,{nullptr,&ux,&uy});
    HydroSolveOptions options;options.compatibleNativePressure=true;InterfaceParametricProjection::Options geometry;
    double maxMass=0,maxFraction=0,maxMomentum=0,maxVolume=0,maxDirection=0;
    for(int step=1;step<=4;++step){const InterfaceStepTime time{step,(step-1)*dt,step*dt,dt};
        hydro.updateMaterials(mesh,{kind,time.begin,HydroMaterialEvolution::TimeDependent,HydroParametricMaterial{&domain,&state,1000,1,.02,.001}});
        const auto velocity=hydro.splineVelocity(mesh,domain,HydroVelocityTimeLevel::Current);
        const auto flux=HydroInterfaceTransport::normalFlux(velocity);const auto flow=HydroInterfaceTransport::frozenParametricFlow(velocity);
        auto next=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,state,time,flow,flux,flux,geometry);
        const auto mass=HydroConservativeTransport::build(mesh,domain,state,next.state,time,1000,1);
        HydroForceDensity force;force.evaluationTime=time.begin;
        const auto solved=hydro.advance(mesh,{kind,{step,time.begin,dt,HydroTimePolicy::ConservativeOneStep},options,force,nullptr,&mass});
        if(!solved.success)throw std::runtime_error("parameter mass coupled step: "+solved.stage+": "+solved.reason);
        require(solved.compatibleNativePressureUsed && solved.conservativeMomentumUsed && solved.coupledUsedSuperLU,
            "parameter mass did not enter the actual compatible Hydro solve");
        maxMass=std::max(maxMass,mass.maxRelativeClosure());maxFraction=std::max(maxFraction,next.state.parametricGeometry().report().maximumInventoryFractionError);
        maxMomentum=std::max(maxMomentum,solved.momentumBalanceRelativeResidual);maxVolume=std::max(maxVolume,std::abs(next.inventoryBalance));
        maxDirection=std::max(maxDirection,next.flux.directionalQuadratureError);state=std::move(next.state);
    }
    double motion=0,velocityChange=0;
    for(int k=0;k<64;++k)motion=std::max(motion,Quantum::Geometry::Parametric::norm(state.parametricGeometry().geometry().sample(k/64.)-initial.sample(k/64.)));
    const auto solved=hydro.state();for(int row=0;row<ux.size();++row)velocityChange=std::max(velocityChange,std::abs(solved.velocityX[row]-ux.value[row]));
    require(state.progress().completedSteps==4 && solved.progress.completedSteps==4 && motion>1e-6 && velocityChange>1e-8 && maxFraction<3e-12 && maxVolume<3e-15,
        "parameter mass dynamics did not advance one shared material/velocity history");
    std::cout<<std::setprecision(17)<<"HYDRO_PARAMETRIC_MASS_DYNAMIC steps=4 end_time="<<4*dt<<" density_ratio=1000 interface_motion="<<motion
        <<" velocity_change="<<velocityChange<<" max_fraction_error="<<maxFraction<<" dual_mass_relative_closure="<<maxMass
        <<" momentum_relative_residual="<<maxMomentum<<" inventory_balance="<<maxVolume<<" directional_quadrature_error="<<maxDirection
        <<" capillary_feedback=0 time_order_claim=0\n";
}
}
void hydroParametricMeasureContract(){
    using namespace Quantum::Geometry::Parametric;
    // At a floating-point pole, the meridian and its virtual axial closure
    // must share exactly one endpoint incidence. The adjacent representable
    // planes exercise the gap produced by reconstructing b as a+(b-a).
    int closureTraces=0;
    for(const auto& curve:{InterfacePrescribedData::parametricSpheroid(.5,{.35,.3},64),
                          InterfacePrescribedData::parametricNeckedSpheroid(.5,{.35,.3},.008/.35,64)}){
        PhaseMeasure phase(curve);const double bottom=curve.sample(0).y,top=curve.sample(1).y;
        for(double pole:{bottom,top})for(double y:{std::nextafter(pole,-std::numeric_limits<double>::infinity()),pole,
                                                 std::nextafter(pole,std::numeric_limits<double>::infinity())})
        for(bool fromLower:{false,true}){
            double length=0;for(const auto& interval:phase.intervals(1,y,fromLower)){
                require(interval[1]>=interval[0]&&interval[0]>=0,"axis closure produced an inverted or negative-radius phase interval");
                length+=interval[1]-interval[0];}
            require(std::isfinite(length)&&(y>=bottom&&y<=top||length==0),"pole closure invented material outside the owned meridian");++closureTraces;
        }
        const auto& interior=phase.intervals(0,0,false);const auto& exterior=phase.intervals(0,0,true);
        require(interior.size()==1&&interior[0][0]==bottom&&interior[0][1]==top&&exterior.empty(),
            "one-sided axial phase traces lost their owned poles");closureTraces+=2;
    }
    std::cout<<"PARAMETRIC_CLOSURE_TRACES samples="<<closureTraces<<" pole_neighbors=1 one_sided_axis=1\n";
    double maximumFractionDefect=0,maximumTotalError=0,maximumMomentDifference=0,maximumStationaryDisplacement=0;
    for(int n:{16,32})for(int segments:{64,128})for(double shift:{0.,std::numeric_limits<double>::epsilon(),-std::numeric_limits<double>::epsilon(),.013}){
        const auto curve=InterfacePrescribedData::parametricEllipse({.5+shift,.5+shift},{.25,.25},segments);
        // The Hysing circle has grid-aligned tangencies on a rectangular
        // domain. Also exercise shifted intersections and the transported
        // polynomial representation, without changing the physical curve.
        PolynomialCurve<9> transported;
        for(const auto& arc:curve.arcs){PolynomialArc<9> next;next.a=arc.a;next.b=arc.b;
            std::copy(arc.p.begin(),arc.p.end(),next.p.begin());transported.arcs.push_back(next);}
        PhaseMeasure phase(curve);PolynomialPhaseMeasure<9> moved(transported);
        const double h=1./n,capacity=h*h;double total=0;
        for(int i=0;i<n;i++)for(int j=0;j<2*n;j++){
            const double x=i*h,y=j*h,v=phase.simpleRectangle(x,x+h,y,y+h);
            const double material=phase.rectanglePhaseMoments(x,x+h,y,y+h).values[0];
            const double advected=moved.simpleRectangle(x,x+h,y,y+h);
            maximumFractionDefect=std::max({maximumFractionDefect,-v/capacity,v/capacity-1});
            maximumMomentDifference=std::max({maximumMomentDifference,std::abs(v-material)/capacity,std::abs(v-advected)/capacity});
            total+=v;
        }
        maximumTotalError=std::max(maximumTotalError,std::abs(total-std::acos(-1.)/16));
        require(maximumFractionDefect<2e-12&&maximumMomentDifference<4e-12&&maximumTotalError<2e-14,
            "grid-aligned tangent circle lost physical phase measure or global area");
        const double bottom=.25+shift,cap=phase.simpleRectangle(.4+shift,.6+shift,bottom-.0625,bottom+1e-4);
        require(cap>1e-7&&cap<1e-5,"pure-cell classification removed a resolved shallow cap");
        // A representable stationary curve is the zero-error constrained
        // minimizer. Its roundoff-scale changes of chart coefficients must
        // not create tangent-cell area defects or a finite geometric motion.
        InterfaceParametricProjection::Grid grid{n,2*n};grid.height=2;
        const auto initial=InterfaceParametricProjection::initialize(curve,grid);
        InterfaceParametricProjection::Options options;options.materialSegments=segments;
        const auto stationary=InterfaceParametricProjection::prepare(transported,grid,initial.inventory(),options);
        for(int k=0;k<256;k++)maximumStationaryDisplacement=std::max(maximumStationaryDisplacement,
            norm(stationary.geometry().sample((k+.37)/256)-curve.sample((k+.37)/256)));
        if(!(maximumStationaryDisplacement<1e-12&&stationary.report().maximumInventoryFractionError<=3e-12)){
            std::cerr<<std::setprecision(17)<<"PARAMETRIC_STATIONARY_FAILURE n="<<n<<" segments="<<segments<<" shift="<<shift
                <<" displacement="<<maximumStationaryDisplacement<<" fraction="<<stationary.report().maximumInventoryFractionError<<'\n';
            require(false,"stationary tangent geometry changed under its own conservative projection");}
    }
    std::cout<<std::setprecision(17)<<"PARAMETRIC_TANGENT_MEASURE maximum_fraction_defect="<<maximumFractionDefect
        <<" maximum_total_error="<<maximumTotalError<<" maximum_moment_difference="<<maximumMomentDifference
        <<" stationary_displacement="<<maximumStationaryDisplacement<<'\n';
}
void hydroParametricMassContract(){
    hydroParametricMeasureContract();
    Quantum::Geometry::UnitPolynomialTerm term;term.factor=1;term.coefficients[2]=1;
    const double endpoint=1.215197711193733e-32;const auto exact=static_cast<long double>(endpoint)*endpoint*endpoint/3;
    const auto integral=Quantum::Geometry::integrateUnitPolynomial({term},0,1,0,endpoint);
    require(std::abs(integral-exact)<=4*std::numeric_limits<double>::epsilon()*exact,
        "quadratic primitive formed unsupported unused high powers at a tiny endpoint");
    Quantum::Geometry::UnitPolynomialTerm large;large.factor=std::ldexp(1.,500);
    for(int i=0;i<7;++i)large.coefficients[i]=i+1;
    auto opposite=large;opposite.weight=-1;
    Quantum::Geometry::UnitPolynomialTerm small;small.factor=1e-8;
    small.coefficients[0]=.25;small.coefficients[1]=-1;small.coefficients[2]=1;
    const std::vector<Quantum::Geometry::UnitPolynomialTerm> cancelling{large,small,opposite};
    const auto bound=Quantum::Geometry::encloseUnitPolynomial(cancelling,.49,.51);
    require(bound[0]<=0 && bound[0]>-1.001e-12 && bound[1]>=1e-12 && bound[1]<1.001e-12,
        "retained polynomial range lost exact cancellation or interval restriction");
    for(int i=0;i<=64;++i){const double s=.49+.02*i/64;
        const long double value=1e-8L*(static_cast<long double>(s)-.5L)*(static_cast<long double>(s)-.5L);
        require(value>=bound[0] && value<=bound[1],"Bernstein range excludes an independent polynomial value");}
    const auto zero=Quantum::Geometry::encloseUnitPolynomial({large,opposite});
    require(zero[0]==0 && zero[1]==0,"exactly cancelled polynomial has a nonzero range");
    Quantum::Geometry::UnitPolynomialTerm high;high.factor=1;high.coefficients[11]=1;high.coefficients[0]=-1./2048;
    const auto highRoots=Quantum::Geometry::isolateUnitPolynomial({high});
    require(highRoots.roots.size()==1&&std::abs(highRoots.roots[0].point-.5)<2e-13,"degree-eleven retained trace lost its actual root");
    const auto highIntegral=Quantum::Geometry::integrateUnitPolynomial({high},0,1,0,1);
    require(std::abs(highIntegral-(1.L/12-1.L/2048))<2e-17L,"degree-eleven retained trace integral is inconsistent");
    const auto highRange=Quantum::Geometry::encloseUnitPolynomial({high});
    require(highRange[0]<=-1./2048&&highRange[1]>=1-1./2048,"degree-eleven Bernstein range excludes its endpoints");
    phaseMass(true);phaseMass(false);solvedMassSteps();
}
void hydroParametricWeightedPhaseContract(){phaseMass(true,true);phaseMass(false,true);}

void hydroParametricMaterialSupportContract(){
    namespace IP=InterfaceParametricProjection;
    using namespace Quantum::Geometry::Parametric;
    // Retain a fine material space on an unchanged physical boundary. This
    // exercises more than the former 2048-arc ceiling without changing its
    // independent volume, endpoint conditions or one-sided curvature.
    for(int kind=0;kind<3;kind++){
        IP::Grid grid;grid.nx=grid.ny=8;grid.axisymmetric=kind!=0;
        const auto initial=kind==0?InterfacePrescribedData::parametricEllipse({.47,.51},{.22,.17},16):
            kind==1?InterfacePrescribedData::parametricSpheroid(.5,{.22,.3},16):InterfacePrescribedData::parametricSupplyCap(.3,.4,16);
        const auto reference=IP::initialize(initial,grid);IP::Source source;source.boundary=initial.boundary;
        for(const auto& arc:initial.arcs){PolynomialArc<9> a;a.a=arc.a;a.b=arc.b;std::copy(arc.p.begin(),arc.p.end(),a.p.begin());source.arcs.push_back(a);}
        IP::Options options;for(int i=0;i<=4096;i++)options.minimumMaterialKnots.push_back(double(i)/4096);
        auto bounded=options;bounded.maximumMaterialSegments=2048;bool refused=false;
        try{IP::prepare(source,grid,reference.inventory(),bounded);}catch(const std::runtime_error& e){refused=std::string(e.what()).find("material-space budget")!=std::string::npos;}
        require(refused,"declared material-space resource bound was ignored");
        const auto fitted=IP::prepare(source,grid,reference.inventory(),options);const auto& curve=fitted.geometry();
        double position=0,curvature=0;
        for(int i=0;i<1024;i++){const double s=(i+.37)/1024;
            position=std::max(position,norm(curve.sample(s)-initial.sample(s)));
            const auto k=[&](const Curve& c){const auto p=c.sample(s),d=c.sample(s,1),dd=c.sample(s,2);const double speed=norm(d);
                return cross(d,dd)/(speed*speed*speed)+(kind?d.y/(speed*p.x):0);};
            curvature=std::max(curvature,std::abs(k(curve)-k(initial))/std::max(1.,std::abs(k(initial))));
        }
        const auto& r=fitted.report();
        const auto assembled=IP::combine({fitted},{1});require(assembled.report().basisCoefficients==r.basisCoefficients&&assembled.report().denseBasisCoefficients==r.denseBasisCoefficients,"material assembly lost basis-storage diagnostics");
        std::cout<<std::setprecision(17)<<"PARAMETRIC_MATERIAL_SUPPORT kind="<<kind<<" arcs="<<curve.arcs.size()<<" position="<<position<<" curvature="<<curvature
            <<" fraction="<<r.maximumInventoryFractionError<<" h2="<<r.relativeH2<<" basis_coefficients="<<r.basisCoefficients<<" dense_equivalent="<<r.denseBasisCoefficients
            <<" original_inventory_unchanged=1 budget_strict=1"<<std::endl;
        // Use the existing translated-geometry position criterion (1e-9),
        // not an additional exact-interpolation claim for a constrained H2
        // projection. Its original H2 and inventory budgets remain strict.
        require(curve.arcs.size()==4096&&position<1e-9&&curvature<2e-7&&r.relativeH2<=options.relativeH2&&r.maximumInventoryFractionError<3e-12
            &&fitted.inventory()==reference.inventory()&&r.basisCoefficients*20<r.denseBasisCoefficients,"compact material support changed geometry, curvature, volume or local storage");
    }
}

void hydroParametricAmrContract(){
    using namespace Quantum::Geometry;
    using namespace Quantum::Geometry::Parametric;
    using namespace InterfaceParametricTransport;
    using Transport=Quantum::Hydrodynamics::HydroCompatibleMomentTransport;
    using Space=Quantum::Hydrodynamics::HydroCompatibleMomentSpace;
    std::ofstream curves,cells,history;
    if(const char* prefix=std::getenv("QUANTUM_PARAMETRIC_AMR_OUTPUT")){
        curves.open(std::string(prefix)+"_curves.csv");cells.open(std::string(prefix)+"_cells.csv");history.open(std::string(prefix)+"_history.csv");
        require(curves&&cells&&history,"cannot open AMR control output");curves<<std::setprecision(17)<<"rz,step,time,parameter,x,y\n";
        cells<<std::setprecision(17)<<"rz,cell,left,right,bottom,top\n";history<<std::setprecision(17)<<"rz,step,time,volume_drift,inventory_error,face_moment_error,mass_closure,position_error,curvature_error\n";
    }
    for(bool rz:{false,true}){
        DictionaryMesh input;input.minCoord={0,0};input.maxCoord={1,1};input.Nx=input.Ny=4;input.maxCellNum=4096;input.maxInitialIterations=0;
        QuadtreeMesh mesh(input);DictionaryMesh::RefineBox outer,inner;
        outer.min=rz?Vector2D{0,.25}:Vector2D{.25,0};outer.max=rz?Vector2D{.5,1}:Vector2D{1,1};outer.level=1;
        inner.min=rz?Vector2D{0,.5}:Vector2D{.5,.375};inner.max=rz?Vector2D{.25,.75}:Vector2D{.75,.625};inner.level=2;
        mesh.applyRefineBoxes({outer,inner});FaceSegmentGraph graph;graph.build(mesh);const auto domain=NativeLeafDomain::build(mesh,graph,rz);
        require(!domain.uniformLeaves(),"AMR control must use mixed native leaves");const auto space=Space::build(mesh,domain);
        const auto initial=rz?InterfacePrescribedData::parametricSpheroid(.3,{.18,.08},16):InterfacePrescribedData::parametricEllipse({.3,.5},{.08,.11},16);
        auto state=InterfaceParametricMotion::initialize(mesh,domain,initial);const auto grid=state.parametricGeometry().grid();
        require(grid.leaves&&grid.size()==domain.cellIndexByRow().size(),"AMR interface does not own actual leaf slots");
        if(cells)for(std::size_t k=0;k<grid.size();k++){const auto c=grid.rectangle(k);cells<<rz<<','<<k<<','<<c.left<<','<<c.right<<','<<c.bottom<<','<<c.top<<'\n';}
        const auto writeCurve=[&](int step){if(curves)for(int k=0;k<=256;k++){const double s=k/256.;const auto p=state.parametricGeometry().geometry().sample(s);
            curves<<rz<<','<<step<<','<<state.progress().endTime<<','<<s<<','<<p.x<<','<<p.y<<'\n';}};writeCurve(0);
        const double volume=std::accumulate(state.parametricGeometry().inventory().begin(),state.parametricGeometry().inventory().end(),0.);
        double maxInventory=0,maxFace=0,maxMass=0,maxPosition=0,maxCurvature=0,maxVolume=0;int coarseFineWet=0;
        for(int step=1;step<=8;step++){
            constexpr double dt=.2;const double speed=step<=4?.5:-.5;const Vector2D velocity=rz?Vector2D{0,speed}:Vector2D{speed,0};
            InterfaceNormalFluxView::Queries q;q.retainedMomentDegree=2;
            q.chart=[&](int row){const auto f=domain.face(mesh,row);const auto& e=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
                const auto a=mesh.verts(e.edgeToVert[0]).coord(),b=mesh.verts(e.edgeToVert[1]).coord();return InterfaceNormalFluxView::FaceChart{a,(a+b)*.5,b};};
            q.polynomial=[chart=q.chart,&mesh,&domain,rz,velocity](int row,int half){const auto f=domain.face(mesh,row);const auto c=chart(row);
                const auto a=half?c.midpoint:c.begin,b=half?c.end:c.midpoint;const double width=f.horizontal?b.x-a.x:b.y-a.y;
                UnitPolynomialTerm t;t.factor=width*(f.horizontal?velocity.y:velocity.x);
                t.coefficients[0]=rz?2*std::acos(-1.)*a.x:1;t.coefficients[1]=rz?2*std::acos(-1.)*(b.x-a.x):0;return std::vector<UnitPolynomialTerm>{t};};
            q.integral=[polynomial=q.polynomial](int row,int half,double a0,double a1,double s0,double s1){return integrateUnitPolynomial(polynomial(row,half),a0,a1,s0,s1);};
            q.density=[polynomial=q.polynomial](int row,int half,bool end,long double s){return evaluateUnitPolynomial(polynomial(row,half),double(end?1-s:s));};
            q.phaseIntegral=[](int,int,const LevelSetMappedSegment&)->long double{throw std::logic_error("unused scalar query in AMR parameter control");};
            const InterfaceNormalFluxView flux(std::move(q));const Flow flow=[=](Vector2D p,double t){const Vector2D move=velocity*t;return MapJet{p+move,{1,0},{0,1},velocity,move};};
            const double begin=state.progress().endTime;const InterfaceStepTime time{step,begin,begin+dt,dt};InterfaceParametricProjection::Options options;
            auto next=InterfaceParametricMotion::prepareTimeLinear(mesh,domain,state,time,flow,flux,flux,options);
            const auto mass=Transport::phaseMassFlux(mesh,domain,space,state,next.state,time,5,.1);maxMass=std::max(maxMass,mass.maximumRelativeClosure);
            // Independent swept prism: exact translation turns every face
            // moment into a volume moment without the face-time algorithm.
            PhaseMeasure oldPhase(state.parametricGeometry().geometry());
            for(const auto f:domain.faces(mesh)){
                if(f.horizontal!=rz)continue;const double plane=rz?f.center.y:f.center.x,lo=(rz?f.center.x:f.center.y)-.5*f.tangentSpacing,hi=lo+f.tangentSpacing;
                const double a=std::min(plane,plane-speed*dt),b=std::max(plane,plane-speed*dt);
                const auto m=rz?oldPhase.rectanglePhaseMoments(lo,hi,a,b,true).values:oldPhase.rectanglePhaseMoments(a,b,lo,hi,false).values;
                const double exact[3]={m[0],m[rz?1:2]+.5*m[0],m[rz?3:5]+m[rz?1:2]+.25*m[0]};
                const auto& actual=InterfaceParametricMotion::faceVolumes(next.state,f.nativeRow);double capacity=std::numeric_limits<double>::infinity();
                for(int c:{f.negativeCell,f.positiveCell})if(c>=0)capacity=std::min(capacity,domain.cellVolumes()[domain.cellRowByIndex()[c]]);
                for(int k=0;k<3;k++)maxFace=std::max(maxFace,std::abs((k?actual.liquidHigher[k-1]:actual.liquid)-std::copysign(1.,speed)*exact[k])/capacity);
                if(f.negativeCell>=0&&f.positiveCell>=0&&std::abs(actual.liquid)>1e-10){
                    const auto width=[&](int c){const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);return mesh.verts(cell.cellToVert[2]).coord().x-mesh.verts(cell.cellToVert[0]).coord().x;};
                    if(width(f.negativeCell)!=width(f.positiveCell))coarseFineWet++;}
            }
            const double offset=.1*(step<=4?step:8-step);auto expected=initial;for(auto& arc:expected.arcs)arc.p[0]+=rz?Point{0,offset}:Point{offset,0};PhaseMeasure phase(expected);
            const auto& accepted=next.state.parametricGeometry();
            for(std::size_t k=0;k<grid.size();k++){const auto c=grid.rectangle(k);const double exact=rz?phase.rectanglePhaseMoments(c.left,c.right,c.bottom,c.top,true).values[0]:phase.simpleRectangle(c.left,c.right,c.bottom,c.top);
                maxInventory=std::max(maxInventory,std::abs(accepted.inventory()[k]-exact)/grid.capacity(k));}
            for(int k=0;k<=256;k++){const double s=k/256.;const auto p=accepted.geometry().sample(s),e=expected.sample(s);maxPosition=std::max(maxPosition,norm(p-e));
                const auto curvature=[&](const Curve& c){const auto v=c.sample(s,1),a=c.sample(s,2),x=c.sample(s);const double meridian=(v.x*a.y-v.y*a.x)/std::pow(norm(v),3);return rz?((s==0||s==1)?2*meridian:meridian+v.y/(norm(v)*x.x)):meridian;};
                const double reference=curvature(expected),error=std::abs(curvature(accepted.geometry())-reference)/std::max(1.,std::abs(reference));require(std::isfinite(error),"nonfinite translated curvature");maxCurvature=std::max(maxCurvature,error);}
            const double drift=std::accumulate(accepted.inventory().begin(),accepted.inventory().end(),0.)-volume;maxVolume=std::max(maxVolume,std::abs(drift));
            if(history)history<<rz<<','<<step<<','<<time.end<<','<<drift<<','<<maxInventory<<','<<maxFace<<','<<maxMass<<','<<maxPosition<<','<<maxCurvature<<'\n';
            state=std::move(next.state);writeCurve(step);
        }
        std::cout<<std::setprecision(17)<<"PARAMETRIC_AMR metric="<<(rz?"RZ":"XY")<<" cells="<<grid.size()<<" steps=8 time="<<state.progress().endTime
            <<" coarse_fine_wet="<<coarseFineWet<<" volume_drift="<<maxVolume<<" inventory_error="<<maxInventory<<" face_moment_error="<<maxFace
            <<" mass_closure="<<maxMass<<" position_error="<<maxPosition<<" curvature_error="<<maxCurvature<<std::endl;
        require(coarseFineWet>0&&maxInventory<3e-12&&maxFace<3e-12&&maxMass<1e-12&&maxPosition<1e-9&&maxCurvature<1e-6,"native AMR transport exceeds its unchanged physical budgets");
    }
}
