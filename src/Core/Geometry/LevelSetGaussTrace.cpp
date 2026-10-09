#include "Core/Geometry/LevelSetGaussTrace.h"
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace LevelSetGaussTrace {
using Jet=VelocityJet;
using Matrix=Eigen::Matrix2d;
namespace {
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
Matrix matrix(const Jet& j){Matrix out;out<<j.dx.x,j.dy.x,j.dx.y,j.dy.y;return out;}
Vector2D vector(const Eigen::Vector2d& v){return {v.x(),v.y()};}
Eigen::Vector2d vector(Vector2D v){return {v.x,v.y};}
}
struct SmoothStep {Vector2D point,timeDerivative,displacement;Matrix derivative;std::array<Vector2D,2> stages;};
SmoothStep gaussStep(const Query& query,Vector2D start,double h,Vector2D initial) {
    constexpr double root=.28867513459481288225;
    const double a[2][2]={{.25,.25-root},{.25+root,.25}},c[2]={.5-root,.5+root};
    // The autonomous start velocity already selected the current patch.
    // Reuse it for every bracket trial from this same point; the stages,
    // Newton solve and differentiated map still query their actual points.
    std::array<Vector2D,2> z;
    for(int i=0;i<2;++i)z[i]=start+initial*(c[i]*h);
    Eigen::Matrix4d lhs;std::array<Jet,2> jet;
    for(int iteration=0;iteration<25;++iteration) {
        for(int i=0;i<2;++i)jet[i]=query(z[i]);
        Eigen::Vector4d residual;lhs=Eigen::Matrix4d::Identity();
        for(int i=0;i<2;++i) {
            Vector2D r=z[i]-start;
            for(int j=0;j<2;++j){r-=jet[j].value*(h*a[i][j]);lhs.block<2,2>(2*i,2*j)-=h*a[i][j]*matrix(jet[j]);}
            residual.segment<2>(2*i)=vector(r);
        }
        if(residual.norm()<12*std::numeric_limits<double>::epsilon()*std::max(1.,start.mod()))break;
        const Eigen::Vector4d update=lhs.fullPivLu().solve(residual);
        for(int i=0;i<2;++i)z[i]-=vector(update.segment<2>(2*i));
        require(iteration<24,"crossing Gauss stage Newton did not converge");
    }
    Eigen::Matrix<double,4,3> rhs;
    for(int i=0;i<2;++i){rhs.block<2,2>(2*i,0)=Matrix::Identity();rhs.block<2,1>(2*i,2)=vector(jet[0].value*a[i][0]+jet[1].value*a[i][1]);}
    const Eigen::Matrix<double,4,3> stageJet=lhs.fullPivLu().solve(rhs);
    SmoothStep out;out.stages=z;out.displacement=(jet[0].value+jet[1].value)*(.5*h);out.point=start+out.displacement;
    out.derivative=Matrix::Identity();out.timeDerivative=(jet[0].value+jet[1].value)*.5;
    for(int i=0;i<2;++i){out.derivative+=.5*h*matrix(jet[i])*stageJet.block<2,2>(2*i,0);out.timeDerivative+=vector(matrix(jet[i])*stageJet.block<2,1>(2*i,2))*(.5*h);}
    return out;
}
Map advance(const Query& source,const PatchQuery& patchQuery,Vector2D begin,double dt,int blocks) {
    require(bool(source) && bool(patchQuery) && blocks>0 && blocks<=4096 && std::isfinite(dt) && std::isfinite(begin.x) && std::isfinite(begin.y),"invalid Gauss characteristic input");
    // Each fixed macro interval has a fixed end time. Internal event durations
    // depend on the initial point and are differentiated by the event equation.
    // The third column differentiates the ENTIRE duration, including all blocks.
    const double sign=dt<0?-1:1,total=std::abs(dt),block=total/blocks;
    const Query query=[&](Vector2D p){auto j=source(p);j.value*=sign;j.dx*=sign;j.dy*=sign;return j;};
    Eigen::Matrix<double,2,3> jet=Eigen::Matrix<double,2,3>::Zero();jet.block<2,2>(0,0)=Matrix::Identity();
    Matrix naive=Matrix::Identity();Vector2D point=begin,displacement{},compensation{};Map out;
    for(int k=0;k<blocks;++k) {
        double elapsed=0;Eigen::RowVector3d elapsedJet=Eigen::RowVector3d::Zero();
        for(int count=0;;++count) {
            require(count<64,"too many native seam crossings in one macro interval");
            const auto velocity=query(point).value;const auto patch=patchQuery(point,velocity);
            double h=block-elapsed;auto step=gaussStep(query,point,h,velocity);
            // The collocation stages must belong to one native branch. A solved
            // macro endpoint can straddle a seam by its stage-solve roundoff.
            // Retain that endpoint unchanged; record the numerical incidence,
            // instead of solving a fictitious positive segment beyond its time.
            const double enclosure=32*std::numeric_limits<double>::epsilon()*std::max(1.,point.mod());
            out.coordinateResolution=std::max(out.coordinateResolution,enclosure);
            const bool endpointInside=contains(patch,step.point);
            const bool endpointEnclosed=step.point.x>=patch.lower.x-enclosure && step.point.x<=patch.upper.x+enclosure && step.point.y>=patch.lower.y-enclosure && step.point.y<=patch.upper.y+enclosure;
            const bool full=endpointEnclosed && contains(patch,step.stages[0]) && contains(patch,step.stages[1]);
            if(full && !endpointInside)++out.endpointIncidences;
            Eigen::RowVector3d hJet;
            if(full) {hJet=-elapsedJet;hJet(2)+=1./blocks;}
            else {
                // The Euler ray does not determine event ordering near a
                // corner. Locate each bracketed exit of the SAME Gauss map,
                // then select the earliest, before differentiating that event.
                int axis=-1;double boundary=0,best=std::numeric_limits<double>::infinity();SmoothStep firstExit;
                const double remaining=block-elapsed;
                for(int d=0;d<2;++d)for(int side:{-1,1}) {
                    const double value=d?point.y:point.x,bound=side<0?(d?patch.lower.y:patch.lower.x):(d?patch.upper.y:patch.upper.x);
                    double hi=remaining;SmoothStep bracket=step;
                    if(side*((d?bracket.point.y:bracket.point.x)-bound)<=0)continue;
                    if(value==bound && side*(d?velocity.y:velocity.x)<=0)continue;
                    const double initialSpeed=d?velocity.y:velocity.x;
                    double lo=0,time=initialSpeed==0?.5*remaining:std::clamp((bound-value)/initialSpeed,0.,remaining);
                    if(!(time>lo && time<hi))time=.5*(lo+hi);SmoothStep candidate;
                    for(int iteration=0;iteration<48;++iteration) {
                        candidate=gaussStep(query,point,time,velocity);const double residual=(d?candidate.point.y:candidate.point.x)-bound;
                        if(std::abs(residual)<16*std::numeric_limits<double>::epsilon()*std::max(1.,std::abs(bound)))break;
                        if(side*residual>0)hi=time;else lo=time;
                        const double speed=d?candidate.timeDerivative.y:candidate.timeDerivative.x;
                        double next=time-residual/speed;if(!(next>lo && next<hi))next=.5*(lo+hi);
                        require(next>0 && iteration<47,"native seam event root did not converge");time=next;
                    }
                    if(time<best){best=time;axis=d;boundary=bound;firstExit=candidate;}
                }
                require(axis>=0,"native seam crossing has no bracketed transverse exit");h=best;step=firstExit;
                if(!(h>0 && h<=block-elapsed+32*std::numeric_limits<double>::epsilon())) {
                    std::cerr<<std::setprecision(17)<<"CROSSING_EVENT_DIAGNOSTIC point="<<point.x<<","<<point.y<<" h="<<h<<" remaining="<<block-elapsed<<" axis="<<axis<<" boundary="<<boundary
                        <<" patch="<<patch.lower.x<<","<<patch.lower.y<<":"<<patch.upper.x<<","<<patch.upper.y<<" endpoint="<<step.point.x<<","<<step.point.y<<std::endl;
                    throw std::runtime_error("native seam root is outside the current time interval");
                }
                require(contains(patch,step.stages[0]) && contains(patch,step.stages[1]),"accepted event stages crossed another native branch");
                const double speed=axis?step.timeDerivative.y:step.timeDerivative.x;
                require(std::abs(speed)>1e-12,"native seam exit is grazing or unresolved");
                hJet=-(step.derivative.row(axis)*jet)/speed;
                const double residual=(axis?step.point.y:step.point.x)-boundary;out.eventResidual=std::max(out.eventResidual,std::abs(residual));
                // The coordinate is the solved event constraint, not a tolerance
                // based change of trajectory or a shifted point for side selection.
                if(axis){step.point.y=boundary;step.displacement.y=boundary-point.y;}
                else{step.point.x=boundary;step.displacement.x=boundary-point.x;}
                if(!contains(patch,step.point)) {
                    // The same stage-solve coordinate enclosure already used
                    // for macro endpoints also identifies a corner incidence.
                    // Retain the other computed coordinate unchanged. This
                    // does not move a point or permit a resolved earlier exit.
                    const bool cornerEnclosed=step.point.x>=patch.lower.x-enclosure && step.point.x<=patch.upper.x+enclosure && step.point.y>=patch.lower.y-enclosure && step.point.y<=patch.upper.y+enclosure;
                    if(!cornerEnclosed) {
                        std::cerr<<std::setprecision(17)<<"CROSSING_CORNER_DIAGNOSTIC start="<<point.x<<","<<point.y<<" end="<<step.point.x<<","<<step.point.y
                            <<" axis="<<axis<<" duration="<<h<<" patch="<<patch.lower.x<<","<<patch.lower.y<<":"<<patch.upper.x<<","<<patch.upper.y<<std::endl;
                        throw std::runtime_error("event hit another patch side first");
                    }
                    ++out.cornerIncidences;
                }
                ++out.events;
            }
            jet=(step.derivative*jet+vector(step.timeDerivative)*hJet).eval();
            naive=step.derivative*naive;point=step.point;
            const Vector2D increment=step.displacement-compensation,next=displacement+increment;
            compensation=(next-displacement)-increment;displacement=next;
            if(full)break;
            elapsed+=h;elapsedJet+=hJet;
            require(elapsed<block,"event at macro endpoint needs explicit two-sided treatment");
        }
    }
    out.point=begin+displacement;out.displacement=displacement;out.derivative=jet.block<2,2>(0,0);out.timeDerivative=vector(jet.col(2))*sign;out.naiveDerivative=naive;return out;
}

Map advanceTimeDependent(const TimeQuery& query,const SmoothPatch& domain,Vector2D begin,double duration,int blocks) {
    require(bool(query) && blocks>0 && blocks<=4096 && std::isfinite(duration) && duration>=0 && contains(domain,begin),
        "invalid time-dependent Gauss characteristic input");
    constexpr double root=.28867513459481288225;
    const double a[2][2]={{.25,.25-root},{.25+root,.25}},c[2]={.5-root,.5+root};
    Map out;out.point=begin;out.derivative=Matrix::Identity();out.timeDerivative={};out.displacement={};
    const double h=duration/blocks;Vector2D compensation{};
    for(int k=0;k<blocks;++k) {
        const Vector2D start=out.point;const double time=duration*(double(k)/blocks);
        const auto initial=query(start,time).value;std::array<Vector2D,2> z;
        for(int i=0;i<2;++i)z[i]=start+initial*(c[i]*h);
        Eigen::Matrix4d lhs;std::array<TimeVelocityJet,2> jets;
        const auto gradient=[](const TimeVelocityJet& j){Matrix m;m<<j.dx.x,j.dy.x,j.dx.y,j.dy.y;return m;};
        for(int iteration=0;iteration<25;++iteration) {
            for(int i=0;i<2;++i){require(contains(domain,z[i]),"time-dependent Gauss stage left the physical domain");
                jets[i]=query(z[i],duration*((k+c[i])/blocks));}
            Eigen::Vector4d residual;lhs=Eigen::Matrix4d::Identity();
            for(int i=0;i<2;++i){Vector2D r=z[i]-start;
                for(int j=0;j<2;++j){r-=jets[j].value*(h*a[i][j]);lhs.block<2,2>(2*i,2*j)-=h*a[i][j]*gradient(jets[j]);}
                residual.segment<2>(2*i)=vector(r);}
            if(residual.norm()<12*std::numeric_limits<double>::epsilon()*std::max(1.,start.mod()))break;
            const Eigen::Vector4d update=lhs.fullPivLu().solve(residual);
            for(int i=0;i<2;++i)z[i]-=vector(update.segment<2>(2*i));
            require(iteration<24,"time-dependent Gauss stage Newton did not converge");
        }
        // Differentiate the same stage equations with respect to the original
        // position and ENTIRE elapsed duration. h'=1/blocks and each stage's
        // physical time derivative is (k+c_i)/blocks, not merely c_i/blocks.
        Eigen::Matrix<double,4,3> rhs;
        for(int i=0;i<2;++i){rhs.block<2,2>(2*i,0)=out.derivative;Vector2D r=out.timeDerivative;
            for(int j=0;j<2;++j)r+=(jets[j].value*(1./blocks)+jets[j].dt*(h*(k+c[j])/blocks))*a[i][j];
            rhs.block<2,1>(2*i,2)=vector(r);}
        const Eigen::Matrix<double,4,3> stage=lhs.fullPivLu().solve(rhs);
        for(int i=0;i<2;++i){out.derivative+=.5*h*gradient(jets[i])*stage.block<2,2>(2*i,0);
            out.timeDerivative+=jets[i].value*(.5/blocks)+jets[i].dt*(.5*h*(k+c[i])/blocks)
                +vector(gradient(jets[i])*stage.block<2,1>(2*i,2))*(.5*h);}
        const Vector2D increment=(jets[0].value+jets[1].value)*(.5*h)-compensation,next=out.displacement+increment;
        compensation=(next-out.displacement)-increment;out.displacement=next;out.point=begin+out.displacement;
        require(contains(domain,out.point),"time-dependent Gauss endpoint left the physical domain");
        out.coordinateResolution=std::max(out.coordinateResolution,32*std::numeric_limits<double>::epsilon()*std::max(1.,start.mod()));
    }
    out.naiveDerivative=out.derivative;return out;
}
} // namespace LevelSetGaussTrace
