#include "Modules/Electricity/Detail/DielectricAxisTransmission.h"
#include "Modules/Electricity/Detail/DielectricCapTransmission.h"
#include "Modules/Electricity/FittedChargeEvolution.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <iostream>
#include <iomanip>

void axisTransmissionGlobalContract();

namespace {
void capTransmissionAlgebraContract()
{
    namespace Cap=Quantum::Electricity::Detail::DielectricCapTransmission;
    const auto check=[](bool ok,const char* why){if(!ok)throw std::runtime_error(why);};
    double maxPotential=0,maxField=0,maxConstraint=0,maxCondition=0,maxFlux=0,maxTinyFlux=0,maxTraction=0;
    for(double h:{1.,.0625})for(double radius:{.4,.25,.1})for(double em:{1.,2.,10.}) {
        const double ep=1,origin=3*h,E=.025,Q=.07*radius*radius/ep;
        const double A=E*radius*radius*radius*(em-ep)/(em+2*ep),B=3*ep/(em+2*ep);
        const double A2=.011*std::pow(radius,4)/(3*ep+2*em),A3=-.004*std::pow(radius,5)/(4*ep+3*em);
        const auto potential=[&](Vector2D p,int liquid) {
            const double r=p.x/h,z=(p.y-origin)/h,rho=std::hypot(r,z);
            const double scale=liquid<0?rho:radius;
            return 2.7+(liquid<0?-E*z+A*z/(rho*rho*rho)+Q/rho:-B*E*z+Q/radius)
                +A2*(z*z-.5*r*r)/std::pow(scale,5)+A3*(z*z*z-1.5*z*r*r)/std::pow(scale,7);
        };
        const auto electric=[&](Vector2D p,int liquid) {
            const double r=p.x/h,z=(p.y-origin)/h,rho=std::hypot(r,z);
            Vector2D result=liquid>=0?Vector2D{0,B*E/h}:
                Vector2D{(3*A*z*r/std::pow(rho,5)+Q*r/std::pow(rho,3))/h,
                    (E-A*(1/std::pow(rho,3)-3*z*z/std::pow(rho,5))+Q*z/std::pow(rho,3))/h};
            for(int l:{2,3}){
                const double polynomial=l==2?z*z-.5*r*r:z*z*z-1.5*z*r*r;
                const double pr=l==2?-r:-3*z*r,pz=l==2?2*z:3*z*z-1.5*r*r;
                const double scale=liquid<0?rho:radius,amplitude=l==2?A2:A3;
                result.x-=amplitude/std::pow(scale,2*l+1)*(pr-(liquid<0?(2*l+1)*polynomial*r/(rho*rho):0))/h;
                result.y-=amplitude/std::pow(scale,2*l+1)*(pz-(liquid<0?(2*l+1)*polynomial*z/(rho*rho):0))/h;
            }
            return result;
        };
        const Vector2D n{.8,.6},point{.8*radius*h,origin+.6*radius*h};
        const std::vector<Cap::Anchor> anchors{{point,n,1/(radius*h),0,origin}};
        std::vector<Cap::Sample> samples;std::vector<double> values;
        for(int j=-3;j<3;++j)for(int i=0;i<5;++i){const Vector2D p{(i+.5)*h,origin+(j+.5)*h};
            samples.push_back({p,-1});values.push_back(potential(p,-1));}
        const Cap::Reconstruction fit(anchors,samples,h,em,ep);
        // P2(u)=(3u^2-1)/2, P3(u)=(5u^3-3u)/2 on the sphere.
        const double u=n.y,amplitude2=.011,amplitude3=-.004;
        const double polynomial2=(3*u*u-1)/2,polynomial3=(5*u*u*u-3*u)/2;
        const double first2=3*u,first3=(15*u*u-3)/2,second2=3,second3=15*u;
        const std::vector<Cap::Jump> q{{(.07+amplitude2*polynomial2+amplitude3*polynomial3)/h,
            (amplitude2*first2+amplitude3*first3)*n.x/(radius*h*h),
            (amplitude2*(second2*n.x*n.x-first2*u)+amplitude3*(second3*n.x*n.x-first3*u))/(radius*radius*h*h*h)}};
        const auto trace=fit.evaluate(0,values,q);const auto inner=electric(point,0),outer=electric(point,-1);
        maxPotential=std::max(maxPotential,std::abs(trace.interfacePotential-potential(point,0)));
        maxField=std::max({maxField,h*std::abs(trace.normalDisplacement/em-(inner.x*n.x+inner.y*n.y)),
            h*std::abs((trace.normalDisplacement+q[0].value)/ep-(outer.x*n.x+outer.y*n.y))});
        const Vector2D t{-n.y,n.x};const double enp=outer.x*n.x+outer.y*n.y,en=inner.x*n.x+inner.y*n.y;
        const double et=inner.x*t.x+inner.y*t.y,etp=outer.x*t.x+outer.y*t.y;
        const double exactTraction=.5*(ep*(enp*enp-etp*etp)-em*(en*en-et*et));
        const double actualTraction=.5*((trace.normalDisplacement+q[0].value)*(trace.normalDisplacement+q[0].value)/ep
            -trace.normalDisplacement*trace.normalDisplacement/em-(ep-em)*trace.tangentialElectricField*trace.tangentialElectricField);
        maxTraction=std::max(maxTraction,std::abs(actualTraction-exactTraction)/std::max(1e-12/(h*h),std::abs(exactTraction)));
        for(int liquid:{-1,0})for(const auto p:liquid<0?std::vector<Vector2D>{{0,origin+1.3*h},{.6*h,origin-.2*h},{1.1*h,origin+1.7*h}}:
            std::vector<Vector2D>{{0,origin},{.3*radius*h,origin+.2*radius*h}}) {
            const auto exact=electric(p,liquid);const double epsilon=liquid<0?ep:em;
            maxPotential=std::max(maxPotential,std::abs(fit.potential(p,liquid,values,q)-potential(p,liquid)));
            maxField=std::max({maxField,h*std::abs(fit.displacementNormal(p,liquid,{1,0},values,q)/epsilon-exact.x),
                h*std::abs(fit.displacementNormal(p,liquid,{0,1},values,q)/epsilon-exact.y)});
        }
        maxConstraint=std::max(maxConstraint,fit.constraintConditionEstimate());maxCondition=std::max(maxCondition,fit.conditionEstimate());
        using Phase=Quantum::Electricity::Detail::DielectricCubicTrace::Phase;
        const auto stream=[&](double r,double z){const double rho=std::hypot(r,z);
            return -E*r*r/2-A*r*r/std::pow(rho,3)+Q*z/rho
                -1.5*A2*z*r*r/std::pow(rho,5)+.5*A3*r*r*(r*r-4*z*z)/std::pow(rho,7);};
        const Cap::Piece radial{{{.6*h,origin-.8*h},{.6*h,origin+.9*h},Phase::Positive},-1};
        const double radialFlux=CartesianMeasure::twoPi*ep*h*(stream(.6,.9)-stream(.6,-.8));
        maxFlux=std::max(maxFlux,std::abs(fit.integratedDisplacement(radial,{1,0},values,q)-radialFlux)/h);
        const double zf=.31*radius,rc=std::sqrt(radius*radius-zf*zf),faceZ=origin+zf*h;
        const std::vector<Cap::Piece> pieces{{{{0,faceZ},{rc*h,faceZ},Phase::Negative},0},
            {{{rc*h,faceZ},{h,faceZ},Phase::Positive},-1}};
        const double innerStream=-B*E*rc*rc/2+A2/std::pow(radius,5)*zf*rc*rc
            +A3/std::pow(radius,7)*(1.5*zf*zf*rc*rc-.375*std::pow(rc,4));
        const double exact=-CartesianMeasure::twoPi*h*(em*innerStream+ep*(stream(1,zf)-stream(rc,zf)));
        const Vector2D owner{.5*h,faceZ-.5*h},neighbor{.5*h,faceZ+.5*h};
        const auto affine=fit.affineCorrection(owner,-1,neighbor,-1,{0,1},h,h,1.3,pieces,q);
        long double assembled=affine.constant;
        for(std::size_t i=1;i<values.size();++i)assembled+=affine.weights[i]*(values[i]-values[0]);
        assembled+=1.3*(CartesianMeasure::twoPi/2)*h*(potential(owner,-1)-potential(neighbor,-1));
        maxFlux=std::max(maxFlux,static_cast<double>(std::abs(assembled-exact)/h));
        double integrated=0;for(const auto& piece:pieces)integrated+=fit.integratedDisplacement(piece,{0,1},values,q);
        maxFlux=std::max(maxFlux,std::abs(integrated-exact)/h);
        auto tiny=radial;tiny.geometry.retainedRange=true;tiny.geometry.offset=.5L;tiny.geometry.fraction=1e-25L;
        const Vector2D midpoint{.6*h,origin+.05*h};
        const double tinyExact=ep*electric(midpoint,-1).x*CartesianMeasure::twoPi*.6*h*1.7*h*1e-25;
        const double tinyActual=fit.integratedDisplacement(tiny,{1,0},values,q);
        maxTinyFlux=std::max(maxTinyFlux,std::abs(tinyActual/tinyExact-1));
        auto shifted=values;for(auto& value:shifted)value+=1024;
        check(h*std::abs(fit.evaluate(0,shifted,q).normalDisplacement-trace.normalDisplacement)<1e-8,
            "cap transmission lost potential gauge invariance");
    }
    check(maxPotential<1e-9 && maxField<1e-8 && maxFlux<1e-9 && maxTinyFlux<1e-7 && maxTraction<1e-7,
        "cap transmission failed exact sphere potential, field or retained physical-face flux");
    std::cout<<"CAP_TRANSMISSION sphere_controls=18 max_potential="<<maxPotential<<" max_hE="<<maxField
        <<" max_flux_over_h="<<maxFlux<<" max_tiny_flux_relative="<<maxTinyFlux<<" max_Maxwell_relative="<<maxTraction
        <<" constraint_condition="<<maxConstraint<<" sample_condition="<<maxCondition<<" global_enabled=0\n";

    // Two independently charged planar traces bounding ONE common gas slab.
    // This isolates joint charge response; it is not a curved two-drop case.
    for(double h:{1.,.0625}) {
        constexpr double em=2,ep=1,a=.025,b=.017,d=-.003;
        const double origin=3*h;const std::array<double,2> plane{{-.2,.2}},slope{{.013,-.021}};
        const auto polynomial=[&](double r,double z){return a*z+b*(z*z-.5*r*r)+d*(z*z*z-1.5*z*r*r);};
        const auto potential=[&](Vector2D p,int liquid){const double r=p.x/h,z=(p.y-origin)/h;
            return 2.7+polynomial(r,z)+(liquid<0?0:slope[liquid]*(z-plane[liquid]));};
        std::vector<Cap::Anchor> anchors;
        for(int i=0;i<2;++i)anchors.push_back({{.3*h,origin+plane[i]*h},{0,i? -1.:1.},0,0,origin+(i?1.:-1.)*h});
        std::vector<Cap::Sample> samples;std::vector<double> values;
        for(int j=0;j<6;++j)for(int i=0;i<5;++i){const Vector2D p{(.2+.35*i)*h,origin+(-.15+.06*j)*h};
            samples.push_back({p,-1});values.push_back(potential(p,-1));}
        const Cap::Reconstruction fit(anchors,samples,h,em,ep);std::vector<Cap::Jump> jumps;
        for(int i=0;i<2;++i){const double r=.3,z=plane[i],nz=anchors[i].normal.y,tx=-nz;
            jumps.push_back({nz*((em-ep)*(a+2*b*z+d*(3*z*z-1.5*r*r))+em*slope[i])/h,
                nz*(em-ep)*(-3*d*r)*tx/(h*h),nz*(em-ep)*(-3*d)/(h*h*h)});}
        double error=0;
        for(int i=0;i<2;++i){const auto state=fit.evaluate(i,values,jumps);const double z=plane[i],r=.3;
            const double en=-(a+2*b*z+d*(3*z*z-1.5*r*r)+slope[i])*anchors[i].normal.y/h;
            error=std::max({error,std::abs(state.interfacePotential-potential(anchors[i].point,i)),
                h*std::abs(state.normalDisplacement/em-en)});
        }
        check(error<1e-8,"coupled cap operator failed independent two-trace charge jets");
        using Phase=Quantum::Electricity::Detail::DielectricCubicTrace::Phase;
        const std::vector<Cap::Piece> pieces{{{{.3*h,origin-.6*h},{.3*h,origin-.2*h},Phase::Negative},0},
            {{{.3*h,origin-.2*h},{.3*h,origin+.2*h},Phase::Positive},-1},
            {{{.3*h,origin+.2*h},{.3*h,origin+.6*h},Phase::Negative},1}};
        double flux=0;for(const auto& piece:pieces)flux+=fit.integratedDisplacement(piece,{1,0},values,jumps);
        const double expected=CartesianMeasure::twoPi*h*.09*b*(em*.8+ep*.4);
        check(std::abs(flux-expected)/h<1e-9,"common gas physical face lost one liquid or reversed flux");
        std::cout<<"CAP_TRANSMISSION common_gas_controls=2 h="<<h<<" max_scaled_error="<<error
            <<" flux_over_h_error="<<std::abs(flux-expected)/h<<" sample_condition="<<fit.conditionEstimate()<<" global_enabled=0\n";
    }
}
}

void axisTransmissionContract()
{
    capTransmissionAlgebraContract();
    namespace Axis=Quantum::Electricity::Detail::DielectricAxisTransmission;
    using Phase=Axis::Phase;
    const auto check=[](bool ok,const char* why){if(!ok)throw std::runtime_error(why);};
    double maxPotential=0,maxField=0,maxTraction=0,maxFlux=0,maxCondition=0,maxConstraintCondition=0;
    for(double h:{1.,.0625})for(double radiusCells:{.4,.25,.1})for(double em:{1.,2.,10.}) {
        const double ep=1,R=radiusCells*h,origin=3*h,Q0=.07,Q1=-.03,A2=.001,A3=-.0004;
        const auto potential=[&](Vector2D p,Phase phase) {
            const double r=p.x/h,z=(p.y-origin)/h;
            double v=2.7+.025*z+.017*(z*z-.5*r*r)-.003*(z*z*z-1.5*z*r*r);
            if(phase==Phase::Positive) {
                const double l=std::log(p.x/R);
                v-=(Q0+Q1*z)*radiusCells/ep*l;
                v-=A2*((z*z-.5*r*r)*l+.5*(r*r-radiusCells*radiusCells));
                v-=A3*((z*z*z-1.5*z*r*r)*l+1.5*z*(r*r-radiusCells*radiusCells));
            }
            return v;
        };
        const auto electric=[&](Vector2D p,Phase phase) {
            const double r=p.x/h,z=(p.y-origin)/h;
            Vector2D e{(.017*r-.009*z*r)/h,(-.025-.034*z+.009*z*z-.0045*r*r)/h};
            if(phase==Phase::Positive) {
                e.x+=(Q0+Q1*z)*radiusCells/(ep*r*h);
                e.y+=Q1*radiusCells/(ep*h)*std::log(p.x/R);
                const double l=std::log(p.x/R);
                e.x+=A2/h*(-r*l+z*z/r+.5*r)+A3/h*(-3*z*r*l+z*z*z/r+1.5*z*r);
                e.y+=A2/h*2*z*l+A3/h*((3*z*z-1.5*r*r)*l+1.5*(r*r-radiusCells*radiusCells));
            }
            return e;
        };
        std::vector<Axis::SamplePoint> samples;std::vector<double> values;
        for(int j=-3;j<=2;++j)for(int i=0;i<5;++i) {
            const Vector2D p{(i+.5)*h,origin+(j+.5)*h};
            check(p.x>R,"analytic thin column unexpectedly has a liquid center");
            samples.push_back({p,Phase::Positive});values.push_back(potential(p,Phase::Positive));
        }
        const Vector2D anchor{R,origin};
        const Axis::Jump jump{(Q0-(em-ep)*.017*radiusCells+.5*ep*A2*radiusCells)/h,
            (Q1+(em-ep)*.009*radiusCells+1.5*ep*A3*radiusCells)/(h*h),2*ep*A2/(radiusCells*h*h*h)};
        const Axis::Reconstruction fit(anchor,{1,0},0,0,em,ep,samples,h);
        const auto state=fit.evaluate(values,jump);
        const auto inner=electric(anchor,Phase::Negative),outer=electric(anchor,Phase::Positive);
        maxPotential=std::max(maxPotential,std::abs(state.interfacePotential-potential(anchor,Phase::Negative)));
        maxField=std::max({maxField,h*std::abs(state.normalDisplacement/em-inner.x),
            h*std::abs((state.normalDisplacement+jump.value)/ep-outer.x),h*std::abs(state.tangentialElectricField-inner.y)});
        const double actual=.5*(std::pow(state.normalDisplacement+jump.value,2)/ep
            -state.normalDisplacement*state.normalDisplacement/em-(ep-em)*state.tangentialElectricField*state.tangentialElectricField);
        const double expected=.5*(ep*(outer.x*outer.x-outer.y*outer.y)-em*(inner.x*inner.x-inner.y*inner.y));
        maxTraction=std::max(maxTraction,std::abs(actual-expected)/std::max(1e-12,std::abs(expected)));
        const double zf=.13,constant=-.025-.034*zf+.009*zf*zf;
        const double baseInner=constant*radiusCells*radiusCells/2-.0045*std::pow(radiusCells,4)/4;
        const double logIntegral=.5*std::log(1/radiusCells)-.25+.25*radiusCells*radiusCells;
        const double cubicLogIntegral=.25*std::log(1/radiusCells)-.0625+.0625*std::pow(radiusCells,4);
        const double exactFlux=CartesianMeasure::twoPi*h*(em*baseInner
            +ep*(constant/2-.0045/4-baseInner)+Q1*radiusCells*logIntegral+ep*2*zf*A2*logIntegral
            +ep*A3*(3*zf*zf*logIntegral-1.5*cubicLogIntegral+1.5*(1-std::pow(radiusCells,4))/4
                -1.5*radiusCells*radiusCells*(1-radiusCells*radiusCells)/2));
        const double faceZ=origin+zf*h;
        const std::vector<Axis::SegmentPiece> pieces{{{0,faceZ},{R,faceZ},Phase::Negative},
            {{R,faceZ},{h,faceZ},Phase::Positive}};
        const Vector2D owner{.5*h,faceZ-.5*h},neighbor{.5*h,faceZ+.5*h};
        const auto correction=fit.affineCorrection(owner,Phase::Positive,neighbor,Phase::Positive,{0,1},
            h,h,1.3,pieces,jump);
        long double reconstructed=correction.constant;
        for(std::size_t i=1;i<values.size();++i)reconstructed+=correction.weights[i]*(values[i]-values[0]);
        reconstructed+=1.3*(CartesianMeasure::twoPi/2)*h*(potential(owner,Phase::Positive)-potential(neighbor,Phase::Positive));
        maxFlux=std::max(maxFlux,static_cast<double>(std::abs(reconstructed-exactFlux)/h));
        for(Phase phase:{Phase::Negative,Phase::Positive})for(double r:phase==Phase::Negative?
            std::vector<double>{0.,.3*R,.8*R}:std::vector<double>{R,1.2*R,.5*h,2.1*h}) {
            const Vector2D p{r,origin+.31*h};const auto e=electric(p,phase);const double eps=phase==Phase::Negative?em:ep;
            maxPotential=std::max(maxPotential,std::abs(fit.potential(p,phase,values,jump)-potential(p,phase)));
            maxField=std::max({maxField,h*std::abs(fit.displacementNormal(p,phase,{1,0},values,jump)/eps-e.x),
                h*std::abs(fit.displacementNormal(p,phase,{0,1},values,jump)/eps-e.y)});
        }
        // A gauge shift preserves fields up to rounding in the supplied data.
        auto shifted=values;for(auto& v:shifted)v+=1024;
        const auto gauge=fit.evaluate(shifted,jump);
        check(h*std::abs(gauge.normalDisplacement-state.normalDisplacement)<1e-8
            && std::abs(gauge.interfacePotential-state.interfacePotential-1024)<1e-8,"axis transmission lost gauge cancellation");
        maxCondition=std::max(maxCondition,fit.conditionEstimate());
        maxConstraintCondition=std::max(maxConstraintCondition,fit.constraintConditionEstimate());
    }
    check(maxPotential<1e-10 && maxField<1e-10 && maxTraction<1e-8 && maxFlux<1e-10,
        "axis transmission failed independent analytic thin-cylinder potential/field/Maxwell stress");
    {
        // A common exact RZ quadratic through a CURVED interface isolates
        // normal rotation, meridional curvature and its arclength derivative.
        const double h=.0625,a=.017,em=2,ep=1,k=.4/h,ks=.1/(h*h);
        const Vector2D anchor{.3*h,3*h},n{.8,.6},t{-.6,.8};
        const auto E=[&](Vector2D p){return Vector2D{a*p.x/(h*h),-2*a*(p.y-anchor.y)/(h*h)};};
        const auto J=[&](Vector2D v){return Vector2D{a*v.x/(h*h),-2*a*v.y/(h*h)};};
        const auto dot=[](Vector2D u,Vector2D v){return u.x*v.x+u.y*v.y;};
        const auto electric=E(anchor);const double en=dot(electric,n),et=dot(electric,t);
        const Axis::Jump jump{(ep-em)*en,(ep-em)*(dot(J(t),n)+k*et),
            (ep-em)*(-k*dot(J(n),n)+2*k*dot(J(t),t)+ks*et-k*k*en)};
        std::vector<Axis::SamplePoint> samples;std::vector<double> values;
        for(int j=-2;j<=2;++j)for(int i=0;i<5;++i) {
            Vector2D p{(i+.5)*h,anchor.y+j*h};samples.push_back({p,Phase::Positive});
            values.push_back(3+a*(std::pow((p.y-anchor.y)/h,2)-.5*std::pow(p.x/h,2)));
        }
        const Axis::Reconstruction fit(anchor,n,k,ks,em,ep,samples,h);
        const auto state=fit.evaluate(values,jump);
        check(h*std::abs(state.normalDisplacement-em*en)<1e-10 && h*std::abs(state.tangentialElectricField-et)<1e-10
            && h*h*std::abs(state.normalDisplacementDerivative-em*(dot(J(t),n)+k*et))<1e-10,
            "axis transmission lost the curved-interface harmonic jump jets");
    }
    std::cout<<std::setprecision(17)<<"AXIS_TRANSMISSION thin_cylinder_cases=18 no_negative_samples=1"
        <<" max_potential_error="<<maxPotential<<" max_h_scaled_field_error="<<maxField
        <<" max_Maxwell_relative_error="<<maxTraction<<" max_sample_condition="<<maxCondition
        <<" max_constraint_condition="<<maxConstraintCondition<<" max_scaled_mixed_flux_error="<<maxFlux
        <<" curved_jump_jets=checked local_analytic_control=1\n";
    axisTransmissionGlobalContract();
}

void axisTransmissionGlobalContract()
{
    using namespace Quantum::Electricity;
    using Solver=PolarizationPoissonFVSolverCellCentered;
    const auto check=[](bool ok,const char* why){if(!ok)throw std::runtime_error(why);};
    DictionaryMesh input;input.Nx=16;input.Ny=24;input.maxCoord={1,1.5};
    input.maxCellNum=1000;input.maxInitialIterations=0;QuadtreeMesh mesh(input);
    constexpr double radius=.025,density=.001;
    Quantum::Geometry::ActiveHermitePatch patch;patch.activeCells.assign(mesh.cellIndexNum(),1);patch.nodes.resize(mesh.vertIndexNum());
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();patch.nodes[v]={p.x*p.x-radius*radius,2*p.x,0,0};
    }
    auto geometry=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,patch,true));
    auto discovered=LevelSetInterfacePanels::discoverNativeComponents(mesh,geometry,true,{},{{1,0},{1,1.5}});
    check(discovered.curves.size()==1,"thin-cylinder global fixture lost its one reflected surface");
    const LevelSetInterfacePanels panels(std::move(discovered.curves.front()));
    FaceSegmentGraph graph;graph.build(mesh);std::vector<Solver::FittedBoundary> boundary(graph.segments().size());
    for(std::size_t i=0;i<boundary.size();++i)if(graph.segments()[i].L<0) {
        boundary[i].type=graph.segments()[i].side==FaceSide::East?Solver::FittedBoundaryType::Dirichlet:
            Solver::FittedBoundaryType::HomogeneousNeumann;boundary[i].value=0;
    }
    DielectricConfig config;config.model="leakyDielectric";config.interfaceDiscretization="sharpHermiteP3";
    config.axisymmetric=true;config.permittivityNegative=2;config.permittivityPositive=1;config.linearResidualTolerance=1e-10;
    const DielectricConductivity conductivity{0,.01};CellField<double> initial(mesh,0.);
    SurfaceChargeField charge(mesh,panels,density);FittedSharpDielectricSolver solver;
    bool strictRejected=false;
    try {solver.solveNative(mesh,graph,boundary,geometry,panels,config,initial,charge.source(),&conductivity);}
    catch(const std::exception& e){strictRejected=std::string(e.what()).find("both phases")!=std::string::npos;}
    check(strictRejected,"default global P3 unexpectedly fabricated a liquid-center sample");
    config.nativeAxisNeckHarmonics=true;
    const auto field=solver.solveNative(mesh,graph,boundary,geometry,panels,config,initial,charge.source(),&conductivity);
    check(field.report.success && field.report.axisNeckHarmonicFaces>0 && field.report.conductionCurrent,
        "explicit thin-neck model did not enter the actual global FV/current solve");
    double potentialError=0,potentialScale=0;
    for(int c=0;c<mesh.cellIndexNum();++c) {
        const auto p=mesh.cells(mesh.cellIndexToLabel[c]).centerCoord();
        check(p.x>radius,"thin-cylinder global fixture has a resolved liquid potential");
        const double exact=density*radius*std::log(1/p.x);
        potentialError=std::max(potentialError,std::abs(field.potential[c]-exact));potentialScale=std::max(potentialScale,std::abs(exact));
    }
    std::vector<InterfaceTraceGeometry::Point> points;
    for(const auto& panel:panels.panels())for(const auto& q:panel.quadratureNodes)points.push_back({q.frame,0,0});
    const auto queries=std::make_shared<const InterfaceTraceGeometry>(geometry,std::move(points));
    const auto traction=reconstructNativeFittedSharpDielectricTractions(mesh,queries,panels,config,field.potential,charge.source());
    double stressError=0,modelResidual=0;
    for(const auto& trace:traction.interfaceTraces) {
        check(trace.axisNeckHarmonicModel && !trace.axisCapExteriorModel,"thin-neck stress used a different electric model");
        stressError=std::max(stressError,std::abs(trace.normalTraction/(.5*density*density)-1));
        modelResidual=std::max(modelResidual,trace.modelRelativeResidual);
    }
    const auto repeated=solver.solveNative(mesh,graph,boundary,geometry,panels,config,initial,charge.source(),&conductivity);
    check(repeated.geometryAssemblyReused,"thin-neck model did not retain immutable prepared geometry");
    double cacheError=0;for(int c=0;c<mesh.cellIndexNum();++c)cacheError=std::max(cacheError,std::abs(field.potential[c]-repeated.potential[c]));
    auto strict=config;strict.nativeAxisNeckHarmonics=false;bool changedPolicyRejected=false;
    try {solver.solveNative(mesh,graph,boundary,geometry,panels,strict,initial,charge.source(),&conductivity);}
    catch(const std::exception& e){changedPolicyRejected=std::string(e.what()).find("both phases")!=std::string::npos;}
    check(changedPolicyRejected && cacheError==0,"thin-neck cache ignored model policy or changed the global potential");
    // A subcell sphere has gas on the axis within the SAME fixed sampling
    // window. It must not acquire the column's logarithmic/regularity model.
    auto dropPatch=patch;
    for(int v=0;v<mesh.vertIndexNum();++v) {
        const auto p=mesh.verts(mesh.vertIndexToLabel[v]).coord();const double z=p.y-.743;
        dropPatch.nodes[v]={p.x*p.x+z*z-radius*radius,2*p.x,2*z,0};
    }
    auto dropGeometry=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,dropPatch,true));
    auto dropComponents=LevelSetInterfacePanels::discoverNativeComponents(mesh,dropGeometry,true);
    check(dropComponents.curves.size()==1,"axis-window rejection fixture lost its subcell sphere");
    const LevelSetInterfacePanels dropPanels(std::move(dropComponents.curves.front()));
    const SurfaceChargeField dropCharge(mesh,dropPanels,density);bool gasAxisRejected=false;
    try {solver.solveNative(mesh,graph,boundary,dropGeometry,dropPanels,config,initial,dropCharge.source(),&conductivity);}
    catch(const std::exception& e){gasAxisRejected=std::string(e.what()).find("both phases")!=std::string::npos;}
    check(gasAxisRejected,"axis-neck model leaked through a gas-axis/detached-cap window");
    FittedChargeProblem problem{graph,boundary,{config,conductivity,1e-8}};
    auto state=initializeFittedChargeState(mesh,problem,geometry,panels,{},initial,density);
    double closure=0;
    for(int step=1;step<=2;++step) {
        FittedChargeMotion motion{geometry,&panels,{},[](Vector2D p){return p;},.001,(step-1)*.001,step*.001};
        auto next=prepareFittedChargeStep(mesh,problem,state,motion);closure=std::max(closure,next.report.relativeChargeClosure);
        check(next.state.completedStep()==step && state.completedStep()==step-1,"thin-neck charge advance lost its clocks");
        state=std::move(next.state);
    }
    std::cout<<"AXIS_TRANSMISSION global_cells="<<mesh.cellIndexNum()<<" harmonic_faces="<<field.report.axisNeckHarmonicFaces
        <<" relative_potential_error="<<potentialError/potentialScale<<" relative_Maxwell_error="<<stressError
        <<" maximum_fit_residual="<<modelResidual<<" charge_steps=2 max_charge_closure="<<closure
        <<" cache_error="<<cacheError<<" strict_policy_rejection_preserved=1 gas_axis_rejected=1\n";
    check(potentialError/potentialScale<.05 && stressError<.1 && closure<1e-8,
        "global thin-cylinder field/stress or current closure exceeds the development control budget");
}
