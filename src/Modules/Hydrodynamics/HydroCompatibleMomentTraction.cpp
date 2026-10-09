#include "Modules/Hydrodynamics/HydroCompatibleMomentTraction.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineBasis.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace Quantum::Hydrodynamics {
namespace {
namespace P=Geometry::Parametric;
using Report=HydroCompatibleMomentTraction::Report;
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
double tangentBound(const P::Arc& arc){double b=0;for(int k=1;k<5;k++)b+=k*P::norm(arc.p[k]);return b;}
P::SurfaceSample surfaceSample(const P::Curve& curve,const P::Arc& arc,double t,bool rz){
    return rz?P::sampleSurfaceRZ(arc,t,curve.boundary):P::sampleSurface(arc,t);
}
template<std::size_t N,class Function,class Accept>
std::array<long double,N> integrate(const Function& fn,const Accept& accept,double lo,double hi,
    const std::array<double,N>& scale,Report& report,int depth=0){
    const auto rule=[&](const auto& x,const auto& w){std::array<long double,N> sum{};
        for(std::size_t k=0;k<x.size();k++){const auto v=fn(.5*(lo+hi+(hi-lo)*x[k]));report.samples++;
            for(std::size_t j=0;j<N;j++){require(std::isfinite(v[j]),"nonfinite sharp traction quadrature");sum[j]+=.5L*(hi-lo)*w[k]*v[j];}}
        return sum;
    };
    const auto coarse=rule(P::gx,P::gw),fine=rule(P::gaussEight().x,P::gaussEight().w);double error=0;
    for(std::size_t k=0;k<N;k++)error=std::max(error,double(std::abs(coarse[k]-fine[k])/std::max(std::abs(fine[k]),static_cast<long double>((hi-lo)*scale[k]))));
    if(error<=512*std::numeric_limits<double>::epsilon()){report.quadratureError=std::max(report.quadratureError,error);accept(lo,hi);return fine;}
    const double middle=.5*(lo+hi);require(depth<24&&middle>lo&&middle<hi,"sharp moment traction quadrature unresolved");
    auto result=integrate<N>(fn,accept,lo,middle,scale,report,depth+1);const auto right=integrate<N>(fn,accept,middle,hi,scale,report,depth+1);
    for(std::size_t k=0;k<N;k++)result[k]+=right[k];return result;
}
using SurfacePartitions=std::vector<std::vector<std::array<double,2>>>;
std::array<std::vector<double>,2> supportPlanes(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const HydroCompatibleMomentComplement* complement=nullptr){
    std::array<std::vector<double>,2> planes;
    for(int index:domain.cellIndexByRow()){
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);const auto p=cell.centerCoord(),h=cell.size();
        for(double x:{p.x-.5*h.x,p.x+.5*h.x})if(x>mesh.minCoord.x&&x<mesh.maxCoord.x)planes[0].push_back(x);
        for(double y:{p.y-.5*h.y,p.y+.5*h.y})if(y>mesh.minCoord.y&&y<mesh.maxCoord.y)planes[1].push_back(y);
    }
    if(complement)for(const auto& c:complement->cells())for(int axis=0;axis<2;axis++)
        for(double line:c.partition[axis])if(line>(axis?mesh.minCoord.y:mesh.minCoord.x)&&line<(axis?mesh.maxCoord.y:mesh.maxCoord.x))planes[axis].push_back(line);
    for(auto& lines:planes){std::sort(lines.begin(),lines.end());lines.erase(std::unique(lines.begin(),lines.end()),lines.end());}
    return planes;
}
std::vector<double> arcCuts(const P::Arc& arc,const std::vector<std::array<double,2>>& intervals,const std::array<std::vector<double>,2>& planes){
    std::vector<double> cuts;for(const auto& ab:intervals){cuts.push_back(ab[0]);cuts.push_back(ab[1]);}if(cuts.empty())return cuts;
    for(int axis=0;axis<2;axis++)for(double line:planes[axis]){auto equation=arc.coordinate(axis);equation[0]-=line;for(double t:P::roots(equation,0,1))cuts.push_back(t);}
    std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());return cuts;
}
SurfacePartitions partitions(const P::Curve& curve,bool rz,const P::SurfaceRestriction* restriction=nullptr){
    SurfacePartitions result;result.reserve(curve.arcs.size());
    for(std::size_t i=0;i<curve.arcs.size();i++){
        std::vector<std::array<double,2>> pieces;const auto critical=P::surfacePartition(curve.arcs[i],curve.boundary,rz);
        const auto windows=restriction?restriction->intervals[i]:std::vector<std::array<double,2>>{{0,1}};
        for(const auto& ab:windows){std::vector<double> cuts{ab[0]};for(double t:critical)if(t>ab[0]&&t<ab[1])cuts.push_back(t);cuts.push_back(ab[1]);
            for(std::size_t k=1;k<cuts.size();k++)pieces.push_back({cuts[k-1],cuts[k]});}
        result.push_back(std::move(pieces));}
    return result;
}
Report measureCurve(const P::Curve& curve,bool rz,const SurfacePartitions& cuts){
    Report report;std::array<long double,3> sum{};
    const double twoPi=2*std::acos(-1.);
    for(std::size_t i=0;i<curve.arcs.size();i++){const auto& arc=curve.arcs[i];
      const double b=tangentBound(arc)*(rz?twoPi*P::range(arc.coordinate(0))[1]:1.),k=std::max(1.,std::abs(surfaceSample(curve,arc,.5,rz).curvature));
      for(const auto& ab:cuts[i]){
        const auto part=integrate<3>([&](double t){const auto q=surfaceSample(curve,arc,t,rz);const double w=q.speed*(rz?twoPi*q.point.x:1.);return std::array<double,3>{w,q.curvature*w,q.curvature*q.curvature*w};},
            [](double,double){},ab[0],ab[1],{b,b*k,b*k*k},report);for(int j=0;j<3;j++)sum[j]+=part[j];}}
    require(sum[0]>=0,"negative moment traction surface");report.length=double(sum[0]);report.meanCurvature=sum[0]>0?double(sum[1]/sum[0]):0;report.curvatureSquaredIntegral=double(sum[2]);return report;
}
}
HydroCompatibleMomentTraction::Report HydroCompatibleMomentTraction::measure(const P::Curve& curve,bool rz){
    return measureCurve(curve,rz,partitions(curve,rz));
}
HydroCompatibleMomentTraction::Result HydroCompatibleMomentTraction::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentSpline& spline,const HydroCompatibleMomentMaterial& material,const Law& law,double sigma){
    if(material.boundaryHandle()->components().empty()){
        material.validate(mesh,domain,space);require(bool(law)&&std::isfinite(sigma)&&sigma>=0,"invalid empty-surface traction law");
        require(std::abs(spline.time()-material.geometryTime())<=16*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(spline.time())),"traction material/velocity time mismatch");
        return {spline.transposePointForces({}),Eigen::VectorXd::Zero(space.divergence().rows()),{}};
    }
    Result result;bool first=true;
    for(std::size_t i=0;i<material.boundaryHandle()->components().size();i++){
        auto part=buildComponent(mesh,domain,space,spline,material,i,law,sigma,false,nullptr);
        if(first){result=std::move(part);first=false;continue;}
        const double length=result.report.length+part.report.length;
        result.report.meanCurvature=length>0?(result.report.length*result.report.meanCurvature+part.report.length*part.report.meanCurvature)/length:0;
        result.report.meanNormal=length>0?(result.report.length*result.report.meanNormal+part.report.length*part.report.meanNormal)/length:0;
        result.report.clipped=result.report.clipped||part.report.clipped;
        result.report.length=length;result.report.curvatureSquaredIntegral+=part.report.curvatureSquaredIntegral;
        result.report.quadratureError=std::max(result.report.quadratureError,part.report.quadratureError);result.report.pieces+=part.report.pieces;result.report.samples+=part.report.samples;
        result.variation+=part.variation;result.pressurePotential+=part.pressurePotential;
    }
    return result;
}
HydroCompatibleMomentTraction::Report HydroCompatibleMomentTraction::measure(const P::Boundary& boundary,bool rz){
    Report result;for(const auto& c:boundary.components()){const auto r=measure(*c.curve,rz);
        result.meanCurvature+=r.length*r.meanCurvature;result.length+=r.length;result.curvatureSquaredIntegral+=r.curvatureSquaredIntegral;result.samples+=r.samples;result.quadratureError=std::max(result.quadratureError,r.quadratureError);}
    if(result.length>0)result.meanCurvature/=result.length;return result;
}
HydroCompatibleMomentTraction::Report HydroCompatibleMomentTraction::measure(const P::Boundary& boundary,const QuadtreeMesh& mesh,const NativeLeafDomain& domain){
    const P::DomainRectangle box{mesh.minCoord.x,mesh.maxCoord.x,mesh.minCoord.y,mesh.maxCoord.y};Report result;
    for(const auto& c:boundary.components()){
        const auto restriction=P::restrictSurface(*c.curve,box);const auto r=measureCurve(*c.curve,domain.axisymmetric(),partitions(*c.curve,domain.axisymmetric(),&restriction));
        result.clipped=result.clipped||restriction.clipped;result.length+=r.length;result.meanCurvature+=r.length*r.meanCurvature;
        result.curvatureSquaredIntegral+=r.curvatureSquaredIntegral;result.samples+=r.samples;result.quadratureError=std::max(result.quadratureError,r.quadratureError);}
    if(result.length>0)result.meanCurvature/=result.length;return result;
}
template<class Velocity> HydroCompatibleMomentTraction::Result HydroCompatibleMomentTraction::buildComponent(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const Velocity& spline,const HydroCompatibleMomentMaterial& material,std::size_t component,const Law& law,double sigma,bool complete,const HydroCompatibleMomentComplement* complement){
    material.validate(mesh,domain,space);require(bool(law),"sharp moment traction requires a law");const bool rz=domain.axisymmetric();const double twoPi=2*std::acos(-1.);
    require(std::isfinite(sigma)&&sigma>=0,"invalid weak surface tension");
    require(std::abs(spline.time()-material.geometryTime())<=16*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(spline.time())),"traction material/velocity time mismatch");
    const auto& curve=*material.boundaryHandle()->components().at(component).curve;const bool open=spline.openBoundary();
    const P::DomainRectangle box{mesh.minCoord.x,mesh.maxCoord.x,mesh.minCoord.y,mesh.maxCoord.y};P::validateContainedSurface(curve,box,rz);const auto restriction=P::restrictSurface(curve,box);
    const auto geometryCuts=partitions(curve,rz,&restriction);Result result;result.report=measureCurve(curve,rz,geometryCuts);result.report.clipped=restriction.clipped;long double normal=0,squared=0;
    const bool generalTrace=spline.usesGeneralTraceSpace();
    for(int end=0;end<2;end++)if(curve.boundary.fixed(end)){
        const auto p=curve.sample(end);const auto v=spline.firstJet({p.x,p.y}).value;
        require(v.x==0&&v.y==0,"fixed free-surface endpoint has nonzero common velocity");
    }
    for(std::size_t index=0;index<curve.arcs.size();index++){const auto& arc=curve.arcs[index];
        if(geometryCuts[index].empty())continue;
        const auto interval=geometryCuts[index].front();
        const double b=tangentBound(arc);const auto probe=law(surfaceSample(curve,arc,.5*(interval[0]+interval[1]),rz));const double scale=std::max(1.,std::hypot(probe.normal,probe.tangential));
        const double metric=rz?twoPi*P::range(arc.coordinate(0))[1]:1.;
        for(const auto& ab:geometryCuts[index]){
          const auto values=integrate<2>([&](double t){const auto q=surfaceSample(curve,arc,t,rz);const auto j=law(q);const double w=q.speed*(rz?twoPi*q.point.x:1.);
                return std::array<double,2>{j.normal*w,(j.normal*j.normal+j.tangential*j.tangential)*w};},
            [](double,double){},ab[0],ab[1],{metric*b*scale,metric*b*scale*scale},result.report);normal+=values[0];squared+=values[1];
        }
    }
    result.report.meanNormal=result.report.length>0?double(normal/result.report.length):0;
    result.pressurePotential=Eigen::VectorXd::Zero(space.divergence().rows());double totalLiquid=0,totalVolume=0;
    for(int c=0;c<mesh.cellIndexNum();c++){
        const auto& m=material.cells()[c];const auto& liquid=m.componentLiquid.at(component);const Eigen::Vector3d alpha=m.volume.ldlt().solve(liquid.col(0));
        result.pressurePotential.segment<3>(3*c)=result.report.meanNormal*alpha;totalLiquid+=liquid(0,0);totalVolume+=m.volume(0,0);
    }
    if(!open)for(int c=0;c<mesh.cellIndexNum();c++)result.pressurePotential[3*c]-=result.report.meanNormal*totalLiquid/totalVolume;
    // Only INTERNAL support boundaries split an interface arc. A physical
    // axis/wall is already an endpoint condition; root-finding its plane
    // again can manufacture a roundoff-sized interval beyond that endpoint.
    const auto planes=supportPlanes(mesh,domain,complement);
    const double stressScale=std::max(1.,result.report.length>0?std::sqrt(double(squared/result.report.length)):0.);
    std::vector<HydroCompatibleMomentSpline::PointForce> forces;
    for(std::size_t index=0;index<curve.arcs.size();index++){const auto& arc=curve.arcs[index];const auto cuts=arcCuts(arc,geometryCuts[index],planes);if(cuts.empty())continue;
        const auto force=[&](double t){const auto q=surfaceSample(curve,arc,t,rz);const auto j=law(q);const double variation=j.normal-sigma*q.curvature-result.report.meanNormal;
            const auto d=q.derivative;const double metric=rz?twoPi*q.point.x:1.,w=-sigma*metric/q.speed;
            // -sigma * integral 2*pi*(r*t.grad(v).t + v_r) ds.
            // The hoop contribution is finite and has no division by r.
            return HydroCompatibleMomentSpline::PointForce{{q.point.x,q.point.y},{metric*(-variation*d.y+j.tangential*d.x)-(rz?twoPi*sigma*q.speed:0.),metric*(variation*d.x+j.tangential*d.y)},
                {w*d.x*d.x,w*d.y*d.x},{w*d.x*d.y,w*d.y*d.y}};};
        for(std::size_t part=1;part<cuts.size();part++){
            const double lo=cuts[part-1],hi=cuts[part],fraction=.5*(lo+hi);
            if(std::none_of(restriction.intervals[index].begin(),restriction.intervals[index].end(),[&](const auto& ab){return fraction>=ab[0]&&fraction<=ab[1];}))continue;
            const auto mid=arc.sampleFraction(fraction);
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[spline.firstJet({mid.x,mid.y}).nativeCell]);const auto h=cell.size(),center=cell.centerCoord();
            const double x0=center.x-.5*h.x,y0=center.y-.5*h.y;
            // Open off-axis RZ uses physical Psi. Value/gradient bases have
            // a common r^2 denominator; monitor that actual rational space.
            // Axis cells divide the exact r^2 factor analytically instead.
            const bool rational=generalTrace&&rz&&x0>0;
            const auto quadrature=[&](auto degreeX,auto degreeY){constexpr int nx=decltype(degreeX)::value,ny=decltype(degreeY)::value,count=6*nx*ny;
            std::array<double,count> scale;scale.fill(tangentBound(arc)*std::max(stressScale,sigma)*(rz?twoPi:1.)/(rational?x0*x0:1.));
            integrate<count>([&](double t){const auto f=force(t);const double x=(f.point.x-x0)/h.x,y=(f.point.y-y0)/h.y;
                const double components[6]={f.measure.x,f.measure.y,f.dx.x,f.dx.y,f.dy.x,f.dy.y};std::array<double,6*decltype(degreeX)::value*decltype(degreeY)::value> v{};
                double px=1/(rational?f.point.x*f.point.x:1.);for(int a=0;a<nx;a++){double py=1;for(int b=0;b<ny;b++){for(int j=0;j<6;j++)v[nx*ny*j+ny*a+b]=components[j]*px*py;py*=y;}px*=x;}return v;},
                [&](double a,double b){const auto& rule=P::gaussEight();for(int k=0;k<8;k++){auto f=force(.5*(a+b+(b-a)*rule.x[k]));const double w=.5*(b-a)*rule.w[k];f.measure=f.measure*w;f.dx=f.dx*w;f.dy=f.dy*w;forces.push_back(f);}},lo,hi,scale,result.report);
            };
            if(complete||generalTrace||!domain.uniformLeaves()){
                if(rz)quadrature(std::integral_constant<int,10>{},std::integral_constant<int,9>{});
                else quadrature(std::integral_constant<int,9>{},std::integral_constant<int,9>{});
            }else if(rz)quadrature(std::integral_constant<int,11>{},std::integral_constant<int,8>{});
            else quadrature(std::integral_constant<int,5>{},std::integral_constant<int,5>{});
            result.report.pieces++;
        }
    }
    // Integrating sigma*kappa*n by parts on a TRUNCATED surface leaves
    // sigma*m*t.v at its exposed ends (m=1 or 2*pi*r). These are control-
    // volume intersections, not a prescribed contact angle or a wall force.
    for(const auto& end:restriction.ends)if(sigma>0){
        const auto q=surfaceSample(curve,curve.arcs[end.arc],end.fraction,rz);auto point=q.point;
        const double tolerance=128*std::numeric_limits<double>::epsilon()*std::max({box.right-box.left,box.top-box.bottom,std::abs(point.x),std::abs(point.y)});
        require(point.x>=box.left-tolerance&&point.x<=box.right+tolerance&&point.y>=box.bottom-tolerance&&point.y<=box.top+tolerance,"surface intersection exceeds coordinate roundoff");
        // Only the evaluation coordinate of a computed intersection is
        // enclosed on the physical face; polynomial coefficients stay fixed.
        point.x=std::clamp(point.x,box.left,box.right);point.y=std::clamp(point.y,box.bottom,box.top);
        const double w=end.sign*sigma*(rz?twoPi*point.x:1.)/q.speed;
        forces.push_back({{point.x,point.y},{w*q.derivative.x,w*q.derivative.y}});
    }
    // C alpha includes the wet physical outer boundary. Its complementary
    // work uses the INTERIOR trace of the same region, including a loop that
    // crosses an outlet. A virtual phase closure contributes no surface area.
    if(open){P::PhaseMeasure phase(curve);
      for(const auto face:domain.faces(mesh))if(face.outwardSign&&face.area>0){
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        const int axis=face.horizontal?1:0;const double plane=axis?a.y:a.x,first=axis?a.x:a.y,last=axis?b.x:b.y;
        for(const auto& wet:phase.intervals(axis,plane,face.outwardSign>0)){
            const double lo=std::max(first,wet[0]),hi=std::min(last,wet[1]);if(!(hi>lo))continue;
            std::vector<double> cuts{lo,hi};for(double line:planes[1-axis])if(line>lo&&line<hi)cuts.push_back(line);
            std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
            for(std::size_t i=1;i<cuts.size();i++){
                const double lower=cuts[i-1],upper=cuts[i],width=upper-lower;
                const bool rational=generalTrace&&rz&&(axis?lower:plane)>0;
                const double radius=axis?upper:plane,minRadius=axis?lower:plane;
                std::array<double,12> scale;scale.fill(std::max(1.,std::abs(result.report.meanNormal))*(rz?twoPi*radius:1.)/(rational?minRadius*minRadius:1.));
                const auto measure=[&](double s){return result.report.meanNormal*face.outwardSign*(rz?twoPi*(axis?s:plane):1.);};
                integrate<12>([&](double s){std::array<double,12> values{};const double r=axis?s:plane;
                    double power=measure(s)/(rational?r*r:1.);for(auto& value:values){value=power;power*=(s-lower)/width;}return values;},
                    [&](double u,double v){const auto& rule=P::gaussEight();for(int k=0;k<8;k++){
                        const double s=.5*(u+v+(v-u)*rule.x[k]),w=.5*(v-u)*rule.w[k]*measure(s);
                        forces.push_back({axis?Vector2D{s,plane}:Vector2D{plane,s},axis?Vector2D{0,w}:Vector2D{w,0}});}},
                    lower,upper,scale,result.report);
            }
        }
      }
    }
    result.variation=spline.transposePointForces(forces);
    if(complement){const Eigen::VectorXd eta=complement->transposePointForces(forces);const int nq=int(result.variation.size());result.variation.conservativeResize(nq+eta.size());result.variation.tail(eta.size())=eta;}
    require(result.variation.allFinite()&&result.pressurePotential.allFinite(),"nonfinite integrated moment traction");return result;
}
HydroCompatibleMomentTraction::Result HydroCompatibleMomentTraction::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentSpline::NativeVelocity& velocity,const HydroCompatibleMomentComplement& complement,
    const HydroCompatibleMomentMaterial& material,const Law& law,double sigma){
    velocity.validate(mesh,domain,space);complement.validate(mesh,domain,space);Result result;bool first=true;
    if(material.boundaryHandle()->components().empty()){
        material.validate(mesh,domain,space);require(bool(law)&&std::isfinite(sigma)&&sigma>=0,"invalid empty-surface traction law");
        require(std::abs(velocity.time()-material.geometryTime())<=16*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(velocity.time())),"traction material/velocity time mismatch");
        return {Eigen::VectorXd::Zero(velocity.coordinates().size()+complement.basisMap().cols()),Eigen::VectorXd::Zero(space.divergence().rows()),{}};
    }
    for(std::size_t i=0;i<material.boundaryHandle()->components().size();i++){
        auto part=buildComponent(mesh,domain,space,velocity,material,i,law,sigma,true,&complement);
        if(first){result=std::move(part);first=false;continue;}const double length=result.report.length+part.report.length;
        result.report.meanCurvature=length>0?(result.report.length*result.report.meanCurvature+part.report.length*part.report.meanCurvature)/length:0;
        result.report.meanNormal=length>0?(result.report.length*result.report.meanNormal+part.report.length*part.report.meanNormal)/length:0;
        result.report.clipped=result.report.clipped||part.report.clipped;
        result.report.length=length;result.report.curvatureSquaredIntegral+=part.report.curvatureSquaredIntegral;
        result.report.quadratureError=std::max(result.report.quadratureError,part.report.quadratureError);result.report.pieces+=part.report.pieces;result.report.samples+=part.report.samples;
        result.variation+=part.variation;result.pressurePotential+=part.pressurePotential;
    }return result;
}
HydroCompatibleMomentSpace::Matrix HydroCompatibleMomentTraction::iterationMetric(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroCompatibleMomentSpace& space,const HydroCompatibleMomentSpline::NativeVelocity& velocity,const HydroCompatibleMomentMaterial& material,double sigma){
    velocity.validate(mesh,domain,space);material.validate(mesh,domain,space);
    require(std::isfinite(sigma)&&sigma>=0,"invalid capillary iteration metric");
    require(velocity.time()==material.geometryTime(),"capillary metric changed the accepted geometry time");
    std::vector<HydroCompatibleMomentSpline::PointForce> rows;const bool rz=domain.axisymmetric();
    const auto planes=supportPlanes(mesh,domain);const MomentSplineBasis::GaussTwelve rule;
    const P::DomainRectangle box{mesh.minCoord.x,mesh.maxCoord.x,mesh.minCoord.y,mesh.maxCoord.y};
    if(sigma>0)for(const auto& component:material.boundaryHandle()->components()){
        const auto& curve=*component.curve;const auto restriction=P::restrictSurface(curve,box);const auto geometryCuts=partitions(curve,rz,&restriction);
        for(std::size_t index=0;index<curve.arcs.size();index++){
            const auto& arc=curve.arcs[index];const auto cuts=arcCuts(arc,geometryCuts[index],planes);
            for(std::size_t part=1;part<cuts.size();part++){
                const double lo=cuts[part-1],hi=cuts[part],mid=.5*(lo+hi);
                if(std::none_of(restriction.intervals[index].begin(),restriction.intervals[index].end(),[&](const auto& ab){return mid>=ab[0]&&mid<=ab[1];}))continue;
                // Fixed positive quadrature of a solver approximation. The
                // original adaptive sharp force and its budgets are unchanged.
                for(int k=0;k<12;k++){
                    const auto at=surfaceSample(curve,arc,mid+.5*(hi-lo)*rule.x[k],rz);
                    const double weight=.5*(hi-lo)*rule.w[k]*sigma*at.speed*(rz?2*std::acos(-1.)*at.point.x:1.);
                    require(std::isfinite(weight)&&weight>=0,"invalid capillary metric measure");
                    const double w=std::sqrt(weight),tx=at.derivative.x/at.speed,ty=at.derivative.y/at.speed;
                    // d_s(u.n) = n . (grad u)t + k_meridian*(u.t).
                    // The principal curvature response is -Delta_Gamma(u.n),
                    // not a Laplacian of every Cartesian velocity component.
                    // Axisymmetric normal speed has no azimuthal derivative;
                    // the hoop curvature is a lower-order shape response.
                    const double curvature=rz?(at.point.x>0?at.curvature-at.normal.x/at.point.x:.5*at.curvature):at.curvature;
                    const double nx=at.normal.x,ny=at.normal.y;
                    rows.push_back({{at.point.x,at.point.y},{w*curvature*tx,w*curvature*ty},
                        {w*tx*nx,w*tx*ny},{w*ty*nx,w*ty*ny}});
                }
            }
        }
    }
    return velocity.pointWorkGram(rows);
}
}
