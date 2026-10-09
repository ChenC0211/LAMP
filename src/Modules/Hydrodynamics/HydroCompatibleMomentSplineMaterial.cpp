#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineNative.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineAssembly.h"
#include "Modules/Hydrodynamics/detail/HydroCompatibleMomentSplineBasis.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentMaterial.h"
#include "Modules/Hydrodynamics/HydroCompatibleMomentComplement.h"
#include "Core/Mesh/CartesianMeasure.h"
#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <exception>
#include <chrono>
#include <cstdlib>
#include <iostream>

namespace Quantum::Hydrodynamics {
namespace {
namespace P=Geometry::Parametric;
constexpr int terms=81;
using Jets=Eigen::Matrix<double,6,Eigen::Dynamic>;
using Strain=Eigen::Matrix<double,4,Eigen::Dynamic>;
using Form=HydroCompatibleMomentSpline::NativeMaterialForms::Cell;
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
Strain strain(const Jets& v,Vector2D p,bool rz){
    Strain e(4,v.cols());e.row(0)=v.row(2);e.row(1)=v.row(5);e.row(2).setZero();
    if(rz){if(p.x)e.row(2)=v.row(0)/p.x;else e.row(2)=v.row(2);}
    e.row(3)=(v.row(4)+v.row(3))/std::sqrt(2.);return e;
}
struct Integral {
    Form form;
    Eigen::VectorXd load;
    double sourceNormSquared=0;
    double estimate=0;
    std::size_t points=0;
    explicit Integral(int n){form.mass=Eigen::MatrixXd::Zero(n,n);form.viscosity=form.mass;form.momentum=Eigen::MatrixXd::Zero(2,n);load=Eigen::VectorXd::Zero(n);
        for(auto& phase:form.phaseVelocityIntegral)phase=Eigen::MatrixXd::Zero(2,n);}
    void clear(){
        form.mass.setZero();form.viscosity.setZero();form.momentum.setZero();load.setZero();
        for(auto& phase:form.phaseVelocityIntegral)phase.setZero();
        form.liquid.fill(0);form.gas.fill(0);sourceNormSquared=estimate=0;points=0;
    }
    void add(const Integral& b,double scale=1){
        form.mass+=scale*b.form.mass;form.viscosity+=scale*b.form.viscosity;form.momentum+=scale*b.form.momentum;
        for(int phase=0;phase<2;phase++)form.phaseVelocityIntegral[phase]+=scale*b.form.phaseVelocityIntegral[phase];
        load+=scale*b.load;sourceNormSquared+=scale*b.sourceNormSquared;
        for(int k=0;k<6;k++){form.liquid[k]+=scale*b.form.liquid[k];form.gas[k]+=scale*b.form.gas[k];}
        estimate+=std::abs(scale)*b.estimate;points+=b.points;
    }
};
}

struct HydroCompatibleMomentSpline::NativeLayout::MaterialCache {
    struct Cell {Eigen::VectorXd massScale,strainScale;std::array<std::optional<Integral>,2> pure;};
    std::optional<HydroCompatibleMomentComplement> complement;
    std::array<double,4> properties{};
    std::shared_ptr<const HydroCompatibleMomentMaterial::ViscosityField> viscosity;
    double tolerance=0;
    std::vector<Cell> cells;
};

HydroCompatibleMomentSpline::NativeMaterialForms HydroCompatibleMomentSpline::NativeLayout::materialForms(
    const HydroCompatibleMomentMaterial& material,const Eigen::MatrixXd& coefficients,double tolerance,const HydroCompatibleMomentComplement* complement,ProjectionRequest* projection) const{
    const auto started=std::chrono::steady_clock::now();const bool profile=std::getenv("QUANTUM_MOMENT_SETUP_PROFILE")!=nullptr;
    const auto timing=[&](const char* stage,Eigen::Index entries=0){if(profile)std::cout<<"MOMENT_MATERIAL_SETUP stage="<<stage<<" cells="<<cells.size()
        <<" entries="<<entries<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<std::endl;};
    require(std::isfinite(tolerance)&&tolerance>=64*std::numeric_limits<double>::epsilon()&&tolerance<1e-3,"invalid native material quadrature tolerance");
    require(coefficients.rows()==terms&&coefficients.cols()==int(cells.size())&&coefficients.allFinite(),"invalid native material field");
    const auto properties=material.properties();const double rhoMax=std::max(properties.liquidDensity,properties.gasDensity),muMax=material.viscosityUpperBound();
    // The cache owns its complement, so the basis identity cannot dangle.
    // Geometry/phase classification is still evaluated on EVERY invocation.
    // Only support, physical properties and already resolved pure-cell forms
    // are reused; a cut cell never receives a previous material snapshot.
    const auto assembly=coordinateAssembly(complement);
    std::lock_guard<std::mutex> lock(mutex);
    const std::array<double,4> key{properties.liquidDensity,properties.gasDensity,properties.liquidViscosity,properties.gasViscosity};
    if(!materialCache||materialCache->properties!=key||materialCache->viscosity!=material.prescribedViscosity()||materialCache->tolerance!=tolerance
        ||bool(materialCache->complement)!=bool(complement)
        ||(complement && &materialCache->complement->basisMap()!=&complement->basisMap())){
        auto fresh=std::make_shared<MaterialCache>();fresh->properties=key;fresh->viscosity=material.prescribedViscosity();fresh->tolerance=tolerance;
        if(complement)fresh->complement=*complement;
        fresh->cells.resize(cells.size());materialCache=std::move(fresh);
    }
    auto& cache=*materialCache;
    std::vector<P::PhaseMeasure> sourceMeasures;for(const auto& c:material.boundaryHandle()->components())sourceMeasures.emplace_back(*c.curve);
    const MomentSplineBasis::GaussTwelve rule;NativeMaterialForms result;
    result.solenoidalCoordinates=int(polynomial.cols());result.complementCoordinates=complement?int(complement->basisMap().cols()):0;
    const auto& offsets=assembly->offsets;const auto& transform=assembly->transform;
    timing("local_coordinate_map",transform.nonZeros());
    Eigen::MatrixXd flatMomentum=Eigen::MatrixXd::Zero(2,offsets.back());
    std::array<Eigen::MatrixXd,2> flatPhaseVelocity;
    for(auto& phase:flatPhaseVelocity)phase=Eigen::MatrixXd::Zero(2,offsets.back());
    Eigen::VectorXd flatLoad=Eigen::VectorXd::Zero(offsets.back());
    std::vector<Eigen::VectorXd> projectionLoads(cells.size());std::vector<double> sourceNorms(cells.size(),0);
    std::vector<NativeMaterialForms> cellResults(cells.size());
    struct Observation {
        std::size_t evaluatedPoints=0,sliceQueries=0,panels=0,axialRefinements=0,radialRefinements=0;
        int maximumAxialDepth=0;
        double sliceWorkerSeconds=0,classificationSeconds=0,partitionSeconds=0;
    };
    std::vector<Observation> observations(profile?cells.size():0);
    const auto integrate=[&](int slot){
        auto& result=cellResults[slot];
        Observation observation;
        // PhaseMeasure owns a mutable horizontal-slice cache. Each cell gets
        // its own copy; only the immutable curves and fixed basis are shared.
        auto measures=sourceMeasures;
        const auto& cell=cells[slot];const auto lo=cell.lower,h=cell.size;const Vector2D hi=lo+h,center=lo+h*.5;
        const double volume=h.x*h.y*(rz?CartesianMeasure::twoPi*center.x:1.);
        const int localTerms=offsets[slot+1]-offsets[slot];
        const auto jets=[&](Vector2D p){Jets out(6,localTerms);out.leftCols(terms)=cellJets(slot,p);
            if(complement)out.rightCols(localTerms-terms)=complement->cellJets(cell.index,p);return out;};
        std::vector<double> xs{lo.x,hi.x},ys{lo.y,hi.y};
        if(complement){xs=complement->cells().at(cell.index).partition[0];ys=complement->cells().at(cell.index).partition[1];}
        // Every donor polynomial boundary is an integration boundary. A
        // target coarse leaf never samples across an unresolved old patch.
        if(projection){for(int index:projection->source.locator.intersectingBox(lo,hi)){
            const auto& old=projection->source.cells[projection->source.slotByCell.at(index)];const auto end=old.lower+old.size;
            if(end.x<=lo.x||old.lower.x>=hi.x||end.y<=lo.y||old.lower.y>=hi.y)continue;
            for(double x:{old.lower.x,end.x})if(x>lo.x&&x<hi.x)xs.push_back(x);
            for(double y:{old.lower.y,end.y})if(y>lo.y&&y<hi.y)ys.push_back(y);}
            for(auto* values:{&xs,&ys}){std::sort(values->begin(),values->end());values->erase(std::unique(values->begin(),values->end()),values->end());}}
        const auto donor=[&](Vector2D p){const auto u=projection->source.firstJet(projection->coefficients,p).value;return Eigen::Vector2d(u.x,u.y);};
        // Per-coordinate physical norms scale the error estimator, including
        // tiny interior basis functions. No clipping of a matrix eigenvalue.
        auto& cached=cache.cells[slot];auto& massScale=cached.massScale;auto& strainScale=cached.strainScale;
        const bool prepareNorm=massScale.size()==0;double sourceScaleSquared=0;
        if(prepareNorm){massScale=Eigen::VectorXd::Zero(localTerms);strainScale=massScale;++result.preparedNormCells;}
        if(prepareNorm||projection){
        for(std::size_t x=1;x<xs.size();x++)for(std::size_t y=1;y<ys.size();y++)
        for(int panel=0;panel<2;panel++)for(int i=0;i<12;i++)for(int j=0;j<12;j++){
            const double hx=xs[x]-xs[x-1],hy=ys[y]-ys[y-1];
            const Vector2D p{xs[x-1]+.5*hx*(panel+.5+.5*rule.x[i]),ys[y-1]+.5*hy*(1+rule.x[j])};
            const double w=.125*hx*hy*rule.w[i]*rule.w[j]*(rz?CartesianMeasure::twoPi*p.x:1.);
            if(prepareNorm){const auto v=jets(p);const auto e=strain(v,p,rz);
                massScale+=w*rhoMax*v.topRows<2>().colwise().squaredNorm().transpose();
                if(muMax>0)strainScale+=2*w*muMax*e.colwise().squaredNorm().transpose();}
            if(projection)sourceScaleSquared+=w*rhoMax*donor(p).squaredNorm();
        }
        if(prepareNorm){massScale=massScale.cwiseSqrt();strainScale=strainScale.cwiseSqrt();}
        }
        const auto error=[&](const Integral& a,const Integral& b){double maximum=0;
            for(int i=0;i<localTerms;i++){
                for(int j=0;j<localTerms;j++){
                    const double m=massScale[i]*massScale[j],k=strainScale[i]*strainScale[j];
                    if(m>0)maximum=std::max(maximum,std::abs(a.form.mass(i,j)-b.form.mass(i,j))/m);
                    if(k>0)maximum=std::max(maximum,std::abs(a.form.viscosity(i,j)-b.form.viscosity(i,j))/k);
                }
                if(massScale[i]>0)for(int d=0;d<2;d++)maximum=std::max(maximum,std::abs(a.form.momentum(d,i)-b.form.momentum(d,i))/(massScale[i]*std::sqrt(rhoMax*volume)));
                if(massScale[i]>0)for(int phase=0;phase<2;phase++)for(int d=0;d<2;d++)maximum=std::max(maximum,
                    std::abs(a.form.phaseVelocityIntegral[phase](d,i)-b.form.phaseVelocityIntegral[phase](d,i))/(massScale[i]*std::sqrt(volume/rhoMax)));
                if(projection&&massScale[i]>0&&sourceScaleSquared>0)maximum=std::max(maximum,std::abs(a.load[i]-b.load[i])/(massScale[i]*std::sqrt(sourceScaleSquared)));
            }
            for(int k=0;k<6;k++)maximum=std::max({maximum,std::abs(a.form.liquid[k]-b.form.liquid[k])/volume,std::abs(a.form.gas[k]-b.form.gas[k])/volume});
            if(projection&&sourceScaleSquared>0)maximum=std::max(maximum,std::abs(a.sourceNormSquared-b.sourceNormSquared)/sourceScaleSquared);
            return maximum;
        };
        // Per-cell scratch belongs to one worker. A rule finishes before
        // recursion begins, so its samples/jet are never retained by a child.
        // Cartesian lines also reuse the local result instead of allocating
        // two dense Gram matrices for every twelve-point interval.
        const int extra=localTerms-terms;const bool dense=rz||extra>0;
        Eigen::MatrixXd velocitySamples,strainSamples;
        if(dense){velocitySamples.resize(24,localTerms);if(muMax>0)strainSamples.resize(48,localTerms);}
        Jets v(6,localTerms);Integral radialScratch(localTerms);
        const auto fillRadialRule=[&](double a,double b,double y,bool liquid,Integral& out){
            out.clear();if(dense&&muMax>0)strainSamples.setZero();
            if(profile)observation.evaluatedPoints+=12;
            const double rho=liquid?properties.liquidDensity:properties.gasDensity;
            // Contract the Cartesian native tensor block along this exact
            // phase/support interval. Reuse its fixed y basis at all twelve
            // positive nodes; no cross-slice sample buffers are needed.
            // RZ and the non-tensor pressure complement retain direct Grams.
            std::optional<MomentSplineBasis::CartesianGram> native;
            if(!rz)native.emplace(cartesianAxis(slot,1,y));
            bool viscous=false;
            for(int i=0;i<12;i++){
                const Vector2D p{a+.5*(b-a)*(1+rule.x[i]),y};const double w=.5*(b-a)*rule.w[i]*(rz?CartesianMeasure::twoPi*p.x:1.);
                const double mu=material.viscosityAt(cell.index,{(p.x-lo.x)/h.x,(p.y-lo.y)/h.y},liquid);
                if(rz)v=jets(p);
                else{
                    const auto x=cartesianAxis(slot,0,p.x);const auto& z=native->axial;native->add(x,w*rho,w*mu);
                    for(int a=0;a<9;a++)for(int b=0;b<9;b++){
                        const int k=9*a+b;v(0,k)=x(0,a)*z(1,b);v(1,k)=-x(1,a)*z(0,b);
                        v(2,k)=x(1,a)*z(1,b);v(3,k)=-x(2,a)*z(0,b);
                        v(4,k)=x(0,a)*z(2,b);v(5,k)=-x(1,a)*z(1,b);
                    }
                    if(extra)v.rightCols(extra)=complement->cellJets(cell.index,p);
                }
                if(dense){velocitySamples.middleRows(2*i,2)=std::sqrt(w*rho)*v.topRows<2>();
                    if(mu>0){strainSamples.middleRows(4*i,4)=std::sqrt(2*w*mu)*strain(v,p,rz);viscous=true;}}
                out.form.momentum.noalias()+=w*rho*v.topRows<2>();
                out.form.phaseVelocityIntegral[liquid?0:1].noalias()+=w*v.topRows<2>();
                if(projection){const Eigen::Vector2d source=donor(p);out.load.noalias()+=w*rho*v.topRows<2>().transpose()*source;out.sourceNormSquared+=w*rho*source.squaredNorm();}
                const double x=(p.x-center.x)/h.x,z=(y-center.y)/h.y;const double powers[6]={1,x,z,x*x,x*z,z*z};
                auto& phase=liquid?out.form.liquid:out.form.gas;for(int k=0;k<6;k++)phase[k]+=w*powers[k];++out.points;
            }
            if(rz){out.form.mass.noalias()=velocitySamples.transpose()*velocitySamples;
                if(viscous)out.form.viscosity.noalias()=strainSamples.transpose()*strainSamples;}
            else{
                native->addTo(out.form.mass,out.form.viscosity);
                const auto complete=[&](Eigen::MatrixXd& form,const Eigen::MatrixXd& values){
                    form.topRightCorner(terms,extra).noalias()=values.leftCols(terms).transpose()*values.rightCols(extra);
                    form.bottomLeftCorner(extra,terms)=form.topRightCorner(terms,extra).transpose();
                    form.bottomRightCorner(extra,extra).noalias()=values.rightCols(extra).transpose()*values.rightCols(extra);
                };
                if(extra){complete(out.form.mass,velocitySamples);if(viscous)complete(out.form.viscosity,strainSamples);}
            }
        };
        const auto radialRule=[&](double a,double b,double y,bool liquid){Integral out(localTerms);fillRadialRule(a,b,y,liquid,out);return out;};
        std::function<Integral(double,double,double,bool,double,int,const Integral&)> radial;
        radial=[&](double a,double b,double y,bool liquid,double budget,int depth,const Integral& coarse){
            const double m=.5*(a+b);const auto left=radialRule(a,m,y,liquid),right=radialRule(m,b,y,liquid);
            auto fine=left;fine.add(right);
            const double difference=error(coarse,fine);
            if(difference<=budget){fine.estimate+=difference;return fine;}
            if(profile)++observation.radialRefinements;
            if(!(depth<12&&m>a&&m<b)){std::ostringstream message;message.precision(17);message<<"native material radial quadrature unresolved: cell="<<cell.index<<" x="<<a<<':'<<b<<" y="<<y<<" difference="<<difference<<" budget="<<budget;throw std::invalid_argument(message.str());}
            // Each child starts with the exact rule already evaluated on
            // this parent's half interval. Preserve nodes, budgets and the
            // left-then-right sum; only redundant evaluations disappear.
            auto out=radial(a,m,y,liquid,.5*budget,depth+1,left);out.add(radial(m,b,y,liquid,.5*budget,depth+1,right));return out;
        };
        const auto addRadial=[&](Integral& sum,double a,double b,double y,bool liquid,double budget){
            // In XY, each partition contains one degree-eight tensor basis
            // and degree-four complement patches. Quadratic prescribed mu
            // raises the largest Gram-integrand degree to eighteen. The
            // twelve-point rule is exact through degree twenty-three, also
            // for a partitioned degree-eight projection donor. RZ jets may
            // contain 1/r factors and retain monitored integration below.
            if(!rz){fillRadialRule(a,b,y,liquid,radialScratch);sum.add(radialScratch);}
            else sum.add(radial(a,b,y,liquid,budget,0,radialRule(a,b,y,liquid)));
        };
        bool cut=false;int full=0;
        std::vector<double> cuts=ys;
        for(auto& measure:measures){
            std::chrono::steady_clock::time_point classificationStarted;
            if(profile)classificationStarted=std::chrono::steady_clock::now();
            const auto classified=measure.rectanglePhaseMoments(lo.x,hi.x,lo.y,hi.y,rz);cut=cut||classified.phase==P::RectanglePhase::Cut;full+=classified.phase==P::RectanglePhase::Full;
            if(profile)observation.classificationSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-classificationStarted).count();
            if(classified.phase!=P::RectanglePhase::Cut)continue;
            std::chrono::steady_clock::time_point partitionStarted;
            if(profile)partitionStarted=std::chrono::steady_clock::now();
            const auto add=[&](const auto& arc){
                auto ys=arc.coordinate(1);std::array<double,5> derivative{};for(int k=1;k<5;k++)derivative[k-1]=k*ys[k];
                cuts.push_back(arc.sampleFraction(0).y);cuts.push_back(arc.sampleFraction(1).y);
                for(double t:P::roots(derivative,0,1))cuts.push_back(arc.sampleFraction(t).y);
                for(double x:xs){auto equation=arc.coordinate(0);equation[0]-=x;for(double t:P::roots(equation,0,1))cuts.push_back(arc.sampleFraction(t).y);}
            };
            for(std::size_t k=0;k<measure.curve.arcs.size();k++){
                const auto& b=measure.bounds[k];if(b.x[1]<lo.x||b.x[0]>hi.x||b.y[1]<lo.y||b.y[0]>hi.y)continue;
                add(measure.curve.arcs[k]);cuts.push_back(measure.sharedEnds[k].y);
            }
            for(const auto& arc:measure.closureEdges)add(arc);
            if(profile)observation.partitionSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-partitionStarted).count();
        }
        require(full<=1&&!(full&&cut),"native material components overlap a cell interior");
        const auto line=[&](double y){Integral out(localTerms);
            if(!cut){for(std::size_t x=1;x<xs.size();x++)addRadial(out,xs[x-1],xs[x],y,full==1,.25*tolerance/(h.y*(xs.size()-1)));return out;}
            std::chrono::steady_clock::time_point sliceStarted;
            if(profile)sliceStarted=std::chrono::steady_clock::now();
            std::vector<std::array<double,2>> intervals;
            for(auto& measure:measures)for(const auto& ab:measure.intervals(y)){
                const double a=std::max(lo.x,ab[0]),b=std::min(hi.x,ab[1]);if(b>a)intervals.push_back({a,b});
            }
            std::sort(intervals.begin(),intervals.end());double begin=lo.x;
            if(profile){++observation.sliceQueries;observation.sliceWorkerSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-sliceStarted).count();}
            // Divide one line budget among its actual phase pieces. A tiny
            // physical piece is retained, but does not demand a vanishing
            // relative tolerance merely because two roots are very close.
            const double pieceBudget=.25*tolerance/(h.y*(2*intervals.size()+1)*(xs.size()-1));
            const auto part=[&](double a,double b,bool liquid){for(std::size_t x=1;x<xs.size();x++){
                const double first=std::max(a,xs[x-1]),last=std::min(b,xs[x]);if(last>first)addRadial(out,first,last,y,liquid,pieceBudget);}};
            for(const auto& ab:intervals){require(ab[0]>=begin,"native material components overlap on a phase slice");part(begin,ab[0],false);part(ab[0],ab[1],true);begin=ab[1];}
            part(begin,hi.x,false);return out;
        };
        cuts.erase(std::remove_if(cuts.begin(),cuts.end(),[&](double y){return y<lo.y||y>hi.y;}),cuts.end());
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        if(profile)observation.panels=cuts.size()-1;
        Integral total(localTerms);const double pi=std::acos(-1.);
        const bool reused=!projection&&!cut&&bool(cached.pure[full]);
        if(reused){total=*cached.pure[full];total.points=0;++result.reusedPureCells;}
        for(std::size_t k=1;!reused&&k<cuts.size();k++){
            const double a=cuts[k-1],b=cuts[k];if(!(b>a))continue;
            const auto axialRule=[&](double first,double last){Integral out(localTerms);
                for(int q=0;q<12;q++){
                    const double t=first+.5*(last-first)*(1+rule.x[q]);
                    // Same monotone-endpoint map as the existing phase-area
                    // integrator. All retained volume weights stay positive.
                    const double y=cut?a+(b-a)*std::pow(std::sin(.5*pi*t),2):a+(b-a)*t;
                    const double jacobian=cut?.5*pi*(b-a)*std::sin(pi*t):b-a;
                    out.add(line(y),.5*(last-first)*rule.w[q]*jacobian);
                }return out;
            };
            std::function<Integral(double,double,double,int,const Integral&)> axial;
            axial=[&](double first,double last,double budget,int depth,const Integral& coarse){
                if(profile)observation.maximumAxialDepth=std::max(observation.maximumAxialDepth,depth);
                const double mid=.5*(first+last);const auto left=axialRule(first,mid),right=axialRule(mid,last);auto fine=left;fine.add(right);
                const double difference=error(coarse,fine);
                if(difference<=budget){fine.estimate+=difference;return fine;}
                if(profile)++observation.axialRefinements;
                if(!(depth<14&&mid>first&&mid<last)){std::ostringstream message;message.precision(17);message<<"native material axial quadrature unresolved: cell="<<cell.index<<" z="<<a<<':'<<b<<" local="<<first<<':'<<last<<" difference="<<difference<<" budget="<<budget;throw std::invalid_argument(message.str());}
                auto out=axial(first,mid,.5*budget,depth+1,left);out.add(axial(mid,last,.5*budget,depth+1,right));return out;
            };
            // The sum of all initial panel budgets is the same cell budget.
            // Allocating it by panel count avoids sub-roundoff relative
            // requirements on small but valid geometric intervals.
            // Only a cut slice has geometry-dependent integration limits.
            // An uncut XY panel is polynomial in y as well as x; all basis
            // and donor discontinuities already belong to xs/ys.
            total.add(!rz&&!cut?axialRule(0,1):axial(0,1,.75*tolerance/(cuts.size()-1),0,axialRule(0,1)));
        }
        if(!projection&&!cut&&!reused)cached.pure[full]=total;
        if(projection){projectionLoads[slot]=total.load;sourceNorms[slot]=total.sourceNormSquared;}
        auto& f=total.form;f.index=cell.index;
        const auto& exact=material.cells().at(cell.index);constexpr int i0[6]={0,0,0,1,1,2},i1[6]={0,1,2,1,2,2};
        for(int k=0;k<6;k++)result.maximumPhaseMomentError=std::max({result.maximumPhaseMomentError,
            std::abs(f.liquid[k]-exact.liquid(i0[k],i1[k]))/volume,
            std::abs(f.liquid[k]+f.gas[k]-exact.volume(i0[k],i1[k]))/volume});
        require(result.maximumPhaseMomentError<=8*tolerance,"native material quadrature differs from the owned Hermite phase moments");
        Eigen::VectorXd u=Eigen::VectorXd::Zero(localTerms);u.head(terms)=coefficients.col(slot);result.kinetic+=.5*u.dot(f.mass*u);result.dissipation+=u.dot(f.viscosity*u);
        const Eigen::Vector2d momentum=f.momentum*u;result.physicalMomentum.x+=momentum[0];result.physicalMomentum.y+=momentum[1];
        result.cells.push_back(std::move(f));
        result.quadratureEstimate=std::max(result.quadratureEstimate,total.estimate);result.quadraturePoints+=total.points;
        if(profile)observations[slot]=observation;
    };
    std::vector<std::exception_ptr> failures(cells.size());
#pragma omp parallel for schedule(dynamic, 1) if(cells.size() > 1)
    for(int slot=0;slot<int(cells.size());slot++){
        try{integrate(slot);}catch(...){failures[slot]=std::current_exception();}
    }
    timing("quadrature");
    // Workers never sum into a shared physical quantity. Merge in the same
    // native order as serial assembly, including the first reported failure.
    for(int slot=0;slot<int(cells.size());slot++){
        if(failures[slot])std::rethrow_exception(failures[slot]);
        auto& local=cellResults[slot];auto& f=local.cells.front();const int localTerms=offsets[slot+1]-offsets[slot];
        result.kinetic+=local.kinetic;result.dissipation+=local.dissipation;
        result.physicalMomentum.x+=local.physicalMomentum.x;result.physicalMomentum.y+=local.physicalMomentum.y;
        result.maximumPhaseMomentError=std::max(result.maximumPhaseMomentError,local.maximumPhaseMomentError);
        result.quadratureEstimate=std::max(result.quadratureEstimate,local.quadratureEstimate);result.quadraturePoints+=local.quadraturePoints;
        result.preparedNormCells+=local.preparedNormCells;result.reusedPureCells+=local.reusedPureCells;
        flatMomentum.middleCols(offsets[slot],localTerms)=f.momentum;
        for(int phase=0;phase<2;phase++)flatPhaseVelocity[phase].middleCols(offsets[slot],localTerms)=f.phaseVelocityIntegral[phase];
        result.cells.push_back(std::move(f));
        if(projection)flatLoad.segment(offsets[slot],localTerms)=projectionLoads[slot];
    }
    if(profile){
        Observation total;std::size_t maximumPoints=0,arcs=0;int maximumCell=-1;
        for(const auto& c:material.boundaryHandle()->components())arcs+=c.curve->arcs.size();
        for(std::size_t slot=0;slot<observations.size();slot++){
            const auto& o=observations[slot];total.evaluatedPoints+=o.evaluatedPoints;total.sliceQueries+=o.sliceQueries;
            total.panels+=o.panels;total.axialRefinements+=o.axialRefinements;total.radialRefinements+=o.radialRefinements;
            total.maximumAxialDepth=std::max(total.maximumAxialDepth,o.maximumAxialDepth);total.sliceWorkerSeconds+=o.sliceWorkerSeconds;
            total.classificationSeconds+=o.classificationSeconds;total.partitionSeconds+=o.partitionSeconds;
            if(o.evaluatedPoints>maximumPoints){maximumPoints=o.evaluatedPoints;maximumCell=cells[slot].index;}
        }
        std::cout<<"MOMENT_MATERIAL_WORK geometry_time="<<material.geometryTime()<<" arcs="<<arcs
            <<" evaluated_points="<<total.evaluatedPoints<<" retained_points="<<result.quadraturePoints
            <<" panels="<<total.panels<<" slice_queries="<<total.sliceQueries<<" slice_worker_seconds="<<total.sliceWorkerSeconds
            <<" classification_worker_seconds="<<total.classificationSeconds<<" partition_worker_seconds="<<total.partitionSeconds
            <<" axial_refinements="<<total.axialRefinements<<" radial_refinements="<<total.radialRefinements
            <<" maximum_axial_depth="<<total.maximumAxialDepth<<" maximum_cell_points="<<maximumPoints<<" maximum_cell="<<maximumCell<<std::endl;
    }
    timing("local_forms");
    result.mass=assembly->assemble([&](int slot)->const Eigen::MatrixXd&{return result.cells[slot].mass;});timing("global_mass",result.mass.nonZeros());
    result.viscosity=assembly->assemble([&](int slot)->const Eigen::MatrixXd&{return result.cells[slot].viscosity;});timing("global_viscosity",result.viscosity.nonZeros());
    result.momentum=flatMomentum*transform;
    for(int phase=0;phase<2;phase++)result.phaseVelocityIntegral[phase]=flatPhaseVelocity[phase]*transform;
    if(projection){projection->load=transform.transpose()*flatLoad;projection->sourceNormSquared=0;for(double value:sourceNorms)projection->sourceNormSquared+=value;}
    return result;
}
}
