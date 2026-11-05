// Independent native-flux, polynomial and point-divergence checks of the
// production regular source. This test supplies no reference to its owner.
#include "Modules/Hydrodynamics/HydroSpatialAdvector.h"
#include "Modules/Hydrodynamics/Detail/HydroRegularAdvectorGeometry.h"
#include "interface_first_jet_fixture.h"
#include <map>
#include <tuple>

namespace {
using namespace LAMP::Hydrodynamics;
using namespace LAMP::Hydrodynamics::Detail;
constexpr double nodes[5]={.046910077030668,.230765344947158,.5,.769234655052841,.953089922969332};
constexpr double weights[5]={.118463442528095,.239314335249683,.284444444444444,.239314335249683,.118463442528095};
double axis(Vector2D p,int n){return n==0?p.x:p.y;}
HydroPhysicalDivergence metrics(const QuadtreeMesh& mesh){
    FaceSegmentGraph graph;graph.build(mesh);
    return HydroPhysicalDivergence::build(mesh,HydroIntegralDivergence::build(mesh,
        HydroDofLayout::build(mesh,NativeLeafDomain::build(mesh,graph,false))));
}
double primitive(Vector2D p){
    if(p.x==0||p.x==1||p.y==0||p.y==1)return 0.;
    const double x=std::sin(double(pi)*p.x),y=std::sin(double(pi)*p.y);
    return x*x*y*y;
}
LevelSetBoundaryTrace::VelocityJet polynomial(Vector2D p,int degree){
    const double x=p.x,y=p.y;
    if(degree==2)return {{x*x+y,x-y*y},{2*x,1,1,-2*y}};
    return {{x*x*x+x*y*y-2*y*y*y,-2*x*x*x+x*x*y+y*y*y},
        {3*x*x+y*y,2*x*y-6*y*y,-6*x*x+2*x*y,x*x+3*y*y}};
}
double power(const std::array<long double,25>& c,Vector2D unit,int dx,int dy,Vector2D h){
    long double sum=0;
    for(int y=dy;y<5;++y)for(int x=dx;x<5;++x){
        long double f=1;for(int k=0;k<dx;++k)f*=x-k;
        for(int k=0;k<dy;++k)f*=y-k;
        sum+=c[5*y+x]*f*std::pow(unit.x,x-dx)*std::pow(unit.y,y-dy)/
            std::pow(h.x,dx)/std::pow(h.y,dy);
    }
    return double(sum);
}
// Differentiate the owner's retained piece polynomial analytically. At a
// seam, use each piece's own closed chart; no finite-difference limit,
// independently fitted field, or common interpolated Hessian is substituted.
struct C2Continuity {
    double native=0,quarter=0,lowerOrder=0;
    std::size_t nativeSeams=0,quarterSeams=0;
};
void checkC2(const QuadtreeMesh& mesh,const HydroRegularAdvectorGeometry& regular,
    const HydroRegularCoefficients& coefficients,const std::vector<HydroRegularPiecePolynomial>& pieces,
    C2Continuity& peak,double& sourcePeak){
    using Edge=std::tuple<int,double,double,double>; // normal axis, line, tangent bounds
    std::map<Edge,std::vector<std::size_t>> edges;
    for(std::size_t k=0;k<pieces.size();++k){
        const auto& p=pieces[k];const auto h=p.upper-p.lower;
        const Vector2D unit{.37,.61},point=p.lower+Vector2D{h.x*unit.x,h.y*unit.y};
        const auto actual=regular.cellJet(coefficients,p.pressure,point).jet;
        const std::array<std::array<double,3>,2> actualFirstJet{{
            {{actual.value.x,actual.derivative.xx,actual.derivative.xy}},
            {{actual.value.y,actual.derivative.yx,actual.derivative.yy}}}};
        for(int component=0;component<2;++component)for(int d=0;d<3;++d){
            const double gap=std::abs(power(p.powers[component],unit,d==1?1:0,d==2?1:0,h)-
                                      actualFirstJet[component][d]);
            check(std::isfinite(gap),"nonfinite actual piece/source first-jet difference");
            sourcePeak=std::max(sourcePeak,gap);
        }
        for(double x:{p.lower.x,p.upper.x})edges[{0,x,p.lower.y,p.upper.y}].push_back(k);
        for(double y:{p.lower.y,p.upper.y})edges[{1,y,p.lower.x,p.upper.x}].push_back(k);
    }
    std::vector<std::array<unsigned,2>> seen(regular.cells().size());
    constexpr std::array<std::array<int,2>,6> derivative{{{{0,0}},{{1,0}},{{0,1}},{{2,0}},{{1,1}},{{0,2}}}};
    for(const auto& entry:edges){
        const int normal=std::get<0>(entry.first);const double line=std::get<1>(entry.first);
        const double lo=std::get<2>(entry.first),hi=std::get<3>(entry.first);
        const auto& owners=entry.second;
        if(owners.size()==1){
            check(line==axis(mesh.minCoord,normal)||line==axis(mesh.maxCoord,normal),
                "actual source piece has an unmatched interior seam");
            continue;
        }
        check(owners.size()==2,"actual source piece seam has more than two owners");
        const auto& a=pieces[owners[0]];const auto& b=pieces[owners[1]];
        const bool native=a.pressure!=b.pressure;
        if(native)++peak.nativeSeams;
        else{
            ++peak.quarterSeams;
            const auto& box=regular.cells()[a.pressure];
            const double fraction=(line-axis(box.lower,normal))/axis(box.upper-box.lower,normal);
            const int knot=static_cast<int>(std::lround(4*fraction));
            check(knot>=1&&knot<=3&&std::abs(fraction-.25*knot)<=1e-12,
                "actual compatible source uses an unexpected internal support knot");
            seen[a.pressure][normal]|=1u<<(knot-1);
        }
        const auto ha=a.upper-a.lower,hb=b.upper-b.lower;
        for(double t:nodes){
            const Vector2D x=normal==0?Vector2D{line,lo+(hi-lo)*t}:Vector2D{lo+(hi-lo)*t,line};
            const Vector2D ua{(x.x-a.lower.x)/ha.x,(x.y-a.lower.y)/ha.y},
                           ub{(x.x-b.lower.x)/hb.x,(x.y-b.lower.y)/hb.y};
            for(int component=0;component<2;++component)for(std::size_t d=0;d<derivative.size();++d){
                const int dx=derivative[d][0],dy=derivative[d][1];
                const double gap=std::abs(power(a.powers[component],ua,dx,dy,ha)-
                                          power(b.powers[component],ub,dx,dy,hb));
                check(std::isfinite(gap),"nonfinite actual source continuity difference");
                if(d<3)peak.lowerOrder=std::max(peak.lowerOrder,gap);
                else if(native)peak.native=std::max(peak.native,gap);
                else peak.quarter=std::max(peak.quarter,gap);
            }
        }
    }
    check(peak.nativeSeams>0&&peak.quarterSeams>0,"C2 check did not exercise both seam kinds");
    for(const auto& knots:seen)check(knots[0]==7u&&knots[1]==7u,
        "actual source did not expose every quarter knot in both directions");
}
void run(int nx,int ny){
    DictionaryMesh dictionary;dictionary.Nx=nx;dictionary.Ny=ny;
    dictionary.maxCellNum=16384;dictionary.maxInitialIterations=0;
    QuadtreeMesh mesh(dictionary);
    const auto geometry=HydroMomentumGeometry::build(mesh,metrics(mesh));
    const auto reconstruction=HydroAdvectorReconstruction::prepare(mesh,geometry);
    const auto& regular=reconstruction.regularGeometry();const auto& layout=geometry.layout();
    double peakFlux=0,peakDivergence=0,peakPolynomial=0,peakSymbolic=0,peakContinuity=0,peakHessian=0;
    C2Continuity c2;
    check(regular.usesCompatibleTensor(),"complete uniform grid lost the compatible tensor source");
    // A two-cell direction has only three normal sample nodes. Do not claim
    // recovery of an unobservable cubic mode; retain both arbitrary native
    // coefficients and the divergence-free closed-wall field on this grid.
    const std::vector<int> scenarios=nx==2&&ny==2?std::vector<int>{0,1}:std::vector<int>{0,1,2,3};
    for(int scenario:scenarios){
        std::vector<double> u(layout.velocityCount()),p(layout.pressureCount());std::vector<int> wall;
        for(int f=0;f<layout.velocityCount();++f){
            const auto& face=regular.faces()[f];
            if(scenario==0){
                u[f]=(face.normalAxis==0?1:-1)*(primitive(face.upper)-primitive(face.lower))/double(face.area);
                if(face.external)wall.push_back(f);
            }else if(scenario==1){
                // Deliberately rough, non-polynomial native FACE MEANS.
                // These coefficients are not treated as point velocities.
                u[f]=std::sin(11.37*(f+1))+.3*std::cos(7.19*(f+1)*(face.normalAxis+1));
            }else{
                for(int q=0;q<5;++q)u[f]+=weights[q]*axis(polynomial(face.lower+(face.upper-face.lower)*nodes[q],scenario).value,face.normalAxis);
            }
        }
        std::optional<HydroAdvectorBoundary> law;
        if(scenario==0)law=HydroAdvectorBoundary::impermeable(mesh,geometry,wall);
        const HydroNativeValues values(layout,u,p);const auto coefficients=regular.apply(mesh,values,law?&*law:nullptr);
        const auto polynomials=regular.piecePolynomials(coefficients,reconstruction.pieces());
        if(scenario<=1)checkC2(mesh,regular,coefficients,polynomials,c2,peakSymbolic);
        for(int f=0;f<layout.velocityCount();++f){
            const auto& face=regular.faces()[f];const int owner=face.negativePressure>=0?face.negativePressure:face.positivePressure;
            double mean=0;
            for(int quarter=0;quarter<4;++quarter){
                const auto begin=face.lower+(face.upper-face.lower)*(quarter*.25),
                           end=face.lower+(face.upper-face.lower)*((quarter+1)*.25);
                const auto normal=regular.normalPolynomial(coefficients,owner,face.normalAxis,begin,end);
                for(int q=0;q<5;++q){
                    const auto x=begin+(end-begin)*nodes[q];
                    const double sample=axis(regular.cellJet(coefficients,owner,x).jet.value,face.normalAxis);
                    mean+=.25*weights[q]*sample;
                    peakSymbolic=std::max(peakSymbolic,std::abs(sample-double(normal.value(false,nodes[q]*normal.length))));
                    if(scenario==0&&face.external)check(std::abs(sample)<=2e-12,"nonzero reconstructed impermeable normal trace");
                }
            }
            peakFlux=std::max(peakFlux,std::abs(mean-u[f])/(1+std::abs(u[f])));
            if(face.negativePressure>=0&&face.positivePressure>=0){
                for(double t:{.17,.5,.83}){
                    const auto x=face.lower+(face.upper-face.lower)*t;
                    const auto a=regular.cellJet(coefficients,face.negativePressure,x).jet,
                               b=regular.cellJet(coefficients,face.positivePressure,x).jet;
                    for(double gap:{a.value.x-b.value.x,a.value.y-b.value.y,
                        a.derivative.xx-b.derivative.xx,a.derivative.xy-b.derivative.xy,
                        a.derivative.yx-b.derivative.yx,a.derivative.yy-b.derivative.yy})
                        peakContinuity=std::max(peakContinuity,std::abs(gap));
                }
            }
        }
        for(int cell=0;cell<layout.pressureCount();++cell){
            const auto& box=regular.cells()[cell];const auto h=box.upper-box.lower;
            const auto stencil=regular.cellStencil(mesh,cell,law?&*law:nullptr);
            for(auto fraction:{Vector2D{.23,.31},Vector2D{.67,.72},Vector2D{.5,.5}}){
                const auto x=box.lower+Vector2D{h.x*fraction.x,h.y*fraction.y};
                const auto actual=regular.cellJet(coefficients,cell,x).jet;
                const double div=actual.derivative.xx+actual.derivative.yy;
                const auto exact=scenario<=1?LevelSetBoundaryTrace::VelocityJet{}:polynomial(x,scenario);
                const double expectedDiv=scenario==0?0:exact.derivative.xx+exact.derivative.yy;
                if(scenario!=1)peakDivergence=std::max(peakDivergence,std::abs(div-expectedDiv));
                if(scenario>=2){
                    for(double gap:{actual.value.x-exact.value.x,actual.value.y-exact.value.y,
                        actual.derivative.xx-exact.derivative.xx,actual.derivative.xy-exact.derivative.xy,
                        actual.derivative.yx-exact.derivative.yx,actual.derivative.yy-exact.derivative.yy})
                        peakPolynomial=std::max(peakPolynomial,std::abs(gap));
                }
                std::array<long double,6> symbolic{};
                for(const auto& row:regular.cellJetWeights(stencil,x))
                    for(int k=0;k<6;++k)symbolic[k]+=row.jet[k]*u[row.velocityDof];
                const std::array<double,6> numerical{{actual.value.x,actual.value.y,
                    actual.derivative.xx,actual.derivative.xy,actual.derivative.yx,actual.derivative.yy}};
                for(int k=0;k<6;++k)peakSymbolic=std::max(peakSymbolic,std::abs(double(symbolic[k])-numerical[k]));
            }
        }
        if(scenario==3){
            for(const auto& piece:polynomials){
                const auto h=piece.upper-piece.lower,x=piece.lower+Vector2D{h.x*.37,h.y*.61};
                const std::array<std::array<double,3>,2> exact{{{{6*x.x,2*x.y,2*x.x-12*x.y}},{{-12*x.x+2*x.y,2*x.x,6*x.y}}}};
                for(int c=0;c<2;++c)for(int k=0;k<3;++k){
                    const int dx=k==0?2:k==1?1:0,dy=k==0?0:k==1?1:2;
                    peakHessian=std::max(peakHessian,std::abs(power(piece.powers[c],{.37,.61},dx,dy,h)-exact[c][k]));
                }
            }
        }
    }
    check(c2.lowerOrder<=2e-10,"actual piece velocity or first derivative is discontinuous");
    check(c2.native<=2e-9,"compatible velocity Hessian is not continuous across native edges");
    check(c2.quarter<=2e-9,"compatible velocity Hessian is not continuous across quarter knots");
    check(peakFlux<=2e-12,"compatible source lost a native face mean");
    check(peakDivergence<=2e-10,"point divergence does not commute with native flux divergence");
    check(peakPolynomial<=2e-10&&peakHessian<=2e-9,"compatible source lost cubic polynomial or Hessian reproduction");
    check(peakSymbolic<=2e-10,"symbolic/polynomial query describes another source");
    check(peakContinuity<=2e-10,"compatible velocity is not C1 across native edges");
    std::cout<<std::setprecision(17)<<"COMPATIBLE_SOURCE,"<<nx<<','<<ny<<','
        <<peakFlux<<','<<peakDivergence<<','<<peakPolynomial<<','<<peakHessian<<','
        <<peakSymbolic<<','<<peakContinuity<<",PASS\n";
    std::cout<<"COMPATIBLE_C2,"<<nx<<','<<ny<<','<<c2.native<<','<<c2.quarter<<','
        <<c2.lowerOrder<<','<<c2.nativeSeams<<','<<c2.quarterSeams<<",PASS\n";
}
}
int main(){
    try{run(8,8);run(8,6);run(2,2);std::cout<<"COMPATIBLE_SOURCE_SUMMARY,PASS\n";return 0;}
    catch(const std::exception& error){std::cerr<<"COMPATIBLE_SOURCE_FAIL: "<<error.what()<<'\n';return 1;}
}
