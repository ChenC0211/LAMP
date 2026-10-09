#include "Modules/Hydrodynamics/HydroConservativeTransport.h"
#include "Modules/Interface/Compatibility/InterfaceConservativePhase.h"
#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/CartesianLeafLocator.h"
#include "Modules/Hydrodynamics/HydroStep.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Hydrodynamics
{
namespace
{
void require(bool valid,const char* reason)
{if(!valid) throw std::invalid_argument(std::string("[HydroConservativeTransport] ")+reason);}
long double lowerFraction(const QuadtreeMesh& mesh,int index,int axis,bool rz)
{
    if(!rz || axis==1) return .5L;
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);
    const long double lo=mesh.verts(cell.cellToVert[0]).coord().x,hi=mesh.verts(cell.cellToVert[2]).coord().x;
    const long double mid=cell.centerCoord().x;
    return (mid-lo)*(mid+lo)/((hi-lo)*(hi+lo));
}
int cellFace(const QuadtreeMesh& mesh,int index,int axis,int side)
{
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[index]);
    return axis?mesh.vEdgeIndexNum()+mesh.hEdges(cell.cellToHEdge[side][0]).index():
        mesh.vEdges(cell.cellToVEdge[side][0]).index();
}

using Duals=HydroConservativeTransport::Duals;
using MassFlux=HydroConservativeTransport::MassFlux;
long double coordinate(Vector2D p,int axis) {return axis?p.y:p.x;}
long double measure(Vector2D lo,Vector2D hi,bool rz)
{
    return (static_cast<long double>(hi.x)-lo.x)*(static_cast<long double>(hi.y)-lo.y)
        *(rz?CartesianMeasure::twoPi/2*(static_cast<long double>(hi.x)+lo.x):1.L);
}
long double area(int normal,long double plane,long double lo,long double hi,bool rz)
{return (hi-lo)*(rz?CartesianMeasure::twoPi*(normal?(hi+lo)/2:plane):1.L);}
struct Piece {int native;Vector2D lo,hi;long double fraction;};
using Pieces=std::array<std::vector<Piece>,2>;
Pieces cellPieces(const QuadtreeMesh& mesh,int c,const Duals& duals,bool rz)
{
    const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
    const auto lo=mesh.verts(cell.cellToVert[0]).coord(),hi=mesh.verts(cell.cellToVert[2]).coord();
    const auto capacity=measure(lo,hi,rz);Pieces pieces;
    for(int component=0;component<2;++component) {
        long double sum=0;
        for(int side=0;side<2;++side)for(int part=0;part<2;++part) {
            const int label=component?cell.cellToHEdge[side][part]:cell.cellToVEdge[side][part];
            if(label<0)continue;
            const int row=component?mesh.vEdgeIndexNum()+mesh.hEdges(label).index():mesh.vEdges(label).index();
            const auto& v=duals.volumes()[row];require(v.active,"inactive momentum piece");
            Piece p{row,{std::max(lo.x,v.lower.x),std::max(lo.y,v.lower.y)},
                {std::min(hi.x,v.upper.x),std::min(hi.y,v.upper.y)},0};
            require(p.hi.x>p.lo.x && p.hi.y>p.lo.y,"empty native half-cell intersection");
            p.fraction=measure(p.lo,p.hi,rz)/capacity;sum+=p.fraction;
            pieces[component].push_back(p);
        }
        require(pieces[component].size()>=2 && pieces[component].size()<=4
            && std::abs(sum-1)<64*std::numeric_limits<double>::epsilon(),"native pieces do not partition cell mass");
        if(pieces[component].size()==2) {
            // Match the existing uniform NS half-mass arithmetic exactly.
            pieces[component][0].fraction=lowerFraction(mesh,c,component,rz);
            pieces[component][1].fraction=1-pieces[component][0].fraction;
        }
    }
    return pieces;
}

// Lift accepted primal fluxes into each physical half-cell partition. For
// every boundary-flux basis vector, require
//   B J_internal = fraction * div_primal - J_external.
// The connected 2--4 piece graph has one redundant row. A local weighted
// graph Laplacian selects its minimum-energy flux, with one multiplier pinned.
// This is a geometric linear map, not a new global flow solve. Applying every
// signed coefficient separately to gross +/- flux retains temporal reversal;
// a negative coefficient swaps directions. Added counterflow is numerical
// diffusion, not a separately certified subcell physical velocity.
std::size_t mixedFluxes(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const Duals& duals,
    const std::vector<Pieces>& cells,const std::vector<MassFlux>& primal,std::vector<MassFlux>& fluxes)
{
    const bool rz=domain.axisymmetric();std::size_t split=0;
    std::vector<std::vector<int>> internal(mesh.cellIndexNum());
    for(int i=0;i<static_cast<int>(duals.segments().size());++i) {
        const auto& s=duals.segments()[i];if(s.ownerCell>=0)internal[s.ownerCell].push_back(i);
    }
    for(int c:domain.cellIndexByRow()) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        struct Boundary {int native,normal,sign;long double plane,lo,hi,area;};
        std::vector<Boundary> boundary;
        for(int normal=0;normal<2;++normal)for(int side=0;side<2;++side)for(int part=0;part<2;++part) {
            const int label=normal?cell.cellToHEdge[side][part]:cell.cellToVEdge[side][part];if(label<0)continue;
            const auto& edge=normal?mesh.hEdges(label):mesh.vEdges(label);
            const int native=edge.index()+(normal?mesh.vEdgeIndexNum():0);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const auto plane=coordinate(a,normal),lo=coordinate(a,1-normal),hi=coordinate(b,1-normal);
            const auto weight=area(normal,plane,lo,hi,rz);
            if(weight==0){require(primal[native].positive==0 && primal[native].negative==0,"nonzero RZ axis mass flux");continue;}
            require(weight>0,"invalid native boundary measure");
            boundary.push_back({native,normal,side?1:-1,plane,lo,hi,weight});
        }
        require(boundary.size()<=8,"native cell has more than eight boundary segments");
        for(int component=0;component<2;++component) {
            const auto& pieces=cells[c][component];const int n=static_cast<int>(pieces.size()),m=static_cast<int>(boundary.size());
            if(n>2)++split;
            struct Edge {int segment,negative,positive;long double weight;};
            std::vector<Edge> edges;long double maxWeight=0;
            const auto index=[&](int native) {
                for(int p=0;p<n;++p)if(pieces[p].native==native)return p;
                throw std::invalid_argument("internal segment leaves its native half-cell partition");
            };
            for(int id:internal[c]) {
                const auto& s=duals.segments()[id];
                if((s.negativeNative>=mesh.vEdgeIndexNum())!=static_cast<bool>(component))continue;
                const int a=index(s.negativeNative),b=index(s.positiveNative);
                const auto x=(pieces[a].lo+pieces[a].hi)*.5,y=(pieces[b].lo+pieces[b].hi)*.5;
                const auto distance=coordinate(y,s.normalAxis)-coordinate(x,s.normalAxis);
                require(distance>0 && s.measure>0,"invalid internal momentum metric");
                const long double weight=s.measure/distance;maxWeight=std::max(maxWeight,weight);
                edges.push_back({id,a,b,weight});
            }
            require(edges.size()>=static_cast<std::size_t>(n-1) && maxWeight>0,"disconnected native momentum partition");
            std::array<std::array<long double,8>,4> rhs{};
            for(int j=0;j<m;++j) {
                const auto& f=boundary[j];long double sum=0;
                for(int p=0;p<n;++p) {
                    const auto& piece=pieces[p];long double share=0;
                    if(coordinate(f.sign>0?piece.hi:piece.lo,f.normal)==f.plane) {
                        const auto lo=std::max(f.lo,coordinate(piece.lo,1-f.normal));
                        const auto hi=std::min(f.hi,coordinate(piece.hi,1-f.normal));
                        if(hi>lo)share=area(f.normal,f.plane,lo,hi,rz)/f.area;
                    }
                    rhs[p][j]=f.sign*(piece.fraction-share);sum+=rhs[p][j];
                }
                require(std::abs(sum)<128*std::numeric_limits<double>::epsilon(),"incompatible local mass-flux partition");
            }
            // Remove the last redundant equation, keeping its multiplier zero.
            std::array<std::array<long double,11>,3> system{};
            for(auto& e:edges) {
                e.weight/=maxWeight;
                if(e.negative<n-1)system[e.negative][e.negative]+=e.weight;
                if(e.positive<n-1)system[e.positive][e.positive]+=e.weight;
                if(e.negative<n-1 && e.positive<n-1){system[e.negative][e.positive]-=e.weight;system[e.positive][e.negative]-=e.weight;}
            }
            for(int p=0;p<n-1;++p)for(int j=0;j<m;++j)system[p][n-1+j]=rhs[p][j];
            for(int k=0;k<n-1;++k) {
                int pivot=k;for(int p=k+1;p<n-1;++p)if(std::abs(system[p][k])>std::abs(system[pivot][k]))pivot=p;
                require(std::isfinite(system[pivot][k]) && system[pivot][k]!=0,"singular local mass-flux partition");
                if(pivot!=k)std::swap(system[pivot],system[k]);
                const auto diagonal=system[k][k];for(int j=k;j<n-1+m;++j)system[k][j]/=diagonal;
                for(int p=0;p<n-1;++p)if(p!=k){const auto value=system[p][k];for(int j=k;j<n-1+m;++j)system[p][j]-=value*system[k][j];}
            }
            std::array<std::array<long double,8>,4> divergence{};
            for(const auto& e:edges) {
                auto& flux=fluxes[e.segment];
                for(int j=0;j<m;++j) {
                    const auto negative=e.negative<n-1?system[e.negative][n-1+j]:0;
                    const auto positive=e.positive<n-1?system[e.positive][n-1+j]:0;
                    const auto coefficient=e.weight*(negative-positive);
                    require(std::isfinite(coefficient),"nonfinite local mass-flux coefficient");
                    divergence[e.negative][j]+=coefficient;divergence[e.positive][j]-=coefficient;
                    const auto& source=primal[boundary[j].native];
                    flux.positive+=coefficient>=0?coefficient*source.positive:-coefficient*source.negative;
                    flux.negative+=coefficient>=0?coefficient*source.negative:-coefficient*source.positive;
                }
            }
            for(int p=0;p<n;++p)for(int j=0;j<m;++j)
                require(std::abs(divergence[p][j]-rhs[p][j])<4096*std::numeric_limits<double>::epsilon(),
                    "local mass-flux map does not commute with the native partition");
        }
    }
    return split;
}
}

HydroConservativeTransport HydroConservativeTransport::singlePhase(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const HydroStepTime& time,double density,const std::vector<double>& velocity)
{
    domain.validateMeshIdentity(mesh);
    require(domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"single phase requires a complete domain");
    require(std::isfinite(density) && density>0 && time.index>0 && std::isfinite(time.begin) && time.begin>=0
        && std::isfinite(time.dt) && time.dt>0 && std::isfinite(time.begin+time.dt) && time.begin+time.dt>time.begin,
        "invalid single-phase density or interval");
    require(velocity.size()==static_cast<std::size_t>(domain.nativeFaceCount()),"single-phase velocity layout mismatch");
    std::vector<long double> mass(mesh.cellIndexNum()),net(mesh.cellIndexNum(),0),gross(net.size(),0);
    std::vector<MassFlux> primal(domain.nativeFaceCount());
    for(int c:domain.cellIndexByRow())mass[c]=density*static_cast<long double>(domain.cellVolumes()[domain.cellRowByIndex()[c]]);
    for(const auto face:domain.faces(mesh)){
        const double u=velocity[face.nativeRow];require(std::isfinite(u),"nonfinite single-phase face velocity");
        require(!face.outwardSign || u==0,"single-phase closed boundary has nonzero normal velocity");
        const long double flux=static_cast<long double>(density)*time.dt*face.area*u;
        primal[face.nativeRow]={std::max(flux,0.L),std::max(-flux,0.L)};
        if(face.negativeCell>=0){net[face.negativeCell]+=flux;gross[face.negativeCell]+=std::abs(flux);}
        if(face.positiveCell>=0){net[face.positiveCell]-=flux;gross[face.positiveCell]+=std::abs(flux);}
    }
    // Use the same extensive inventory/throughput scale as the dual balance
    // below. A pure relative-to-throughput test becomes singular near a
    // stagnation point and can reject an accepted incompressible solve even
    // when its full-step mass defect is below inventory roundoff.
    double primalClosure=0;
    for(int c:domain.cellIndexByRow()){
        const auto residual=std::abs(net[c])/(2*mass[c]+gross[c]);
        require(std::isfinite(residual) && residual<=4096*std::numeric_limits<double>::epsilon(),
            "single-phase velocity does not preserve local cell mass");
        primalClosure=std::max(primalClosure,static_cast<double>(residual));}
    auto result=fromCellMasses(mesh,domain,{time.index,time.begin,time.begin+time.dt,time.dt},mass,mass,primal);
    result.maxRelativeClosure_=std::max(result.maxRelativeClosure_,primalClosure);
    result.singlePhaseVelocity_=velocity;
    return result;
}
void HydroConservativeTransport::validateSinglePhaseMidpoint(const std::vector<double>& velocity) const
{
    require(isSinglePhase() && velocity==singlePhaseVelocity_,"single-phase momentum and mass use different midpoint velocities");
}
HydroConservativeTransport HydroConservativeTransport::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& previous,const InterfaceConservativePhase& next,
    const InterfaceStepTime& time,double liquidDensity,double gasDensity)
{
    previous.validate(mesh,domain);next.validate(mesh,domain);
    require(domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "requires a complete native domain");
    require(next.follows(previous,time),"inventory step lineage or time mismatch");
    require(std::isfinite(liquidDensity) && liquidDensity>0 && std::isfinite(gasDensity) && gasDensity>0,
        "phase densities must be finite and positive");
    std::vector<long double> before(mesh.cellIndexNum()),after(mesh.cellIndexNum());
    for(int c:domain.cellIndexByRow()){
        const auto capacity=previous.capacities()[c];
        before[c]=liquidDensity*previous.liquidVolumes()[c]+gasDensity*(capacity-previous.liquidVolumes()[c]);
        after[c]=liquidDensity*next.liquidVolumes()[c]+gasDensity*(capacity-next.liquidVolumes()[c]);
    }
    const int count=domain.nativeFaceCount();std::vector<MassFlux> primal(count);
    require(next.directionalFluxes().size()==static_cast<std::size_t>(count),"missing accepted directional fluxes");
    for(int row=0;row<count;++row){const auto& f=next.directionalFluxes()[row];
        primal[row]={liquidDensity*f.liquidPositive+gasDensity*(f.totalPositive-f.liquidPositive),
            liquidDensity*f.liquidNegative+gasDensity*(f.totalNegative-f.liquidNegative)};}
    return fromCellMasses(mesh,domain,time,before,after,primal);
}
HydroConservativeTransport HydroConservativeTransport::build(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& previous,const InterfaceState& next,const InterfaceStepTime& time,double liquidDensity,double gasDensity,
    const InterfaceState* midpoint)
{
    InterfaceParametricMotion::validateDomain(mesh,domain,previous);InterfaceParametricMotion::validateDomain(mesh,domain,next);
    require(next.parametricFollows(previous,time),"parametric inventory step lineage or time mismatch");
    require(std::isfinite(liquidDensity) && liquidDensity>0 && std::isfinite(gasDensity) && gasDensity>0,
        "phase densities must be finite and positive");
    const auto& flux=next.parametricFluxes();require(flux.hasDirections,"parametric step has no retained gross directional fluxes");
    std::vector<long double> before(mesh.cellIndexNum()),after(mesh.cellIndexNum());
    const auto& cells=InterfaceParametricMotion::cellIndices(previous);
    for(std::size_t slot=0;slot<cells.size();++slot){const int c=cells[slot];const auto capacity=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        const auto first=static_cast<long double>(previous.parametricGeometry().inventory()[slot]),last=static_cast<long double>(next.parametricGeometry().inventory()[slot]);
        before[c]=liquidDensity*first+gasDensity*(capacity-first);after[c]=liquidDensity*last+gasDensity*(capacity-last);
    }
    std::vector<MassFlux> primal(domain.nativeFaceCount());
    for(const auto face:domain.faces(mesh)){
        const auto& f=InterfaceParametricMotion::faceVolumes(next,face.nativeRow);
        primal[face.nativeRow]={liquidDensity*f.liquidPositive+gasDensity*(f.totalPositive-f.liquidPositive),
            liquidDensity*f.liquidNegative+gasDensity*(f.totalNegative-f.liquidNegative)};
    }
    auto result=fromCellMasses(mesh,domain,time,before,after,primal);
    result.parameterBefore_=previous.parametricGeometry().geometryHandle();
    result.parameterAfter_=next.parametricGeometry().geometryHandle();
    if(midpoint){InterfaceParametricMotion::validateDomain(mesh,domain,*midpoint);
        const InterfaceStepTime half{time.index,time.begin,time.begin+.5*time.dt,.5*time.dt};
        require(midpoint->parametricFollows(previous,half),"midpoint geometry is not from the actual old interface and half interval");
        result.parameterMiddle_=midpoint->parametricGeometry().geometryHandle();}
    return result;
}
HydroConservativeTransport HydroConservativeTransport::fromCellMasses(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceStepTime& time,const std::vector<long double>& cellBefore,const std::vector<long double>& cellAfter,
    const std::vector<MassFlux>& primal)
{
    HydroConservativeTransport result;result.duals_=Duals::build(mesh,domain);result.radialRows_=mesh.vEdgeIndexNum();
    result.identity_=domain.meshIdentity();result.axisymmetric_=domain.axisymmetric();
    result.uniformLeaves_=domain.uniformLeaves();
    result.step_=time.index;result.begin_=time.begin;result.dt_=time.dt;result.sourceDensity_.resize(mesh.cellIndexNum());
    const auto count=domain.nativeFaceCount();result.before_.assign(count,0);result.after_.assign(count,0);
    const bool rz=domain.axisymmetric();
    std::vector<Pieces> pieces;
    if(!domain.uniformLeaves()) {
        pieces.resize(mesh.cellIndexNum());
        for(int c:domain.cellIndexByRow())pieces[c]=cellPieces(mesh,c,result.duals_,rz);
    }
    for(int c:domain.cellIndexByRow()) {
        const auto capacity=domain.cellVolumes()[domain.cellRowByIndex()[c]];
        const auto before=cellBefore[c],after=cellAfter[c];
        result.sourceDensity_[c]=static_cast<double>(before/capacity);
        require(std::isfinite(before) && std::isfinite(after) && before>0 && after>0,"invalid cell mixture mass");
        if(!domain.uniformLeaves()) {
            for(const auto& component:pieces[c])for(const auto& piece:component) {
                result.before_[piece.native]+=piece.fraction*before;result.after_[piece.native]+=piece.fraction*after;
            }
        } else for(int axis=0;axis<2;++axis) {
            const auto lower=lowerFraction(mesh,c,axis,rz);
            for(int side=0;side<2;++side) {
                const auto fraction=side?1-lower:lower;const int row=cellFace(mesh,c,axis,side);
                result.before_[row]+=fraction*before;result.after_[row]+=fraction*after;
            }
        }
    }
    for(int row=0;row<count;++row) {
        require(std::isfinite(primal[row].positive) && std::isfinite(primal[row].negative)
            && primal[row].positive>=0 && primal[row].negative>=0,"invalid accepted phase mass flux");
    }
    std::vector<long double> net(count,0),gross(count,0);
    result.fluxes_.resize(result.duals_.segments().size());
    if(!domain.uniformLeaves())result.splitCellMappings_=mixedFluxes(mesh,domain,result.duals_,pieces,primal,result.fluxes_);
    for(std::size_t i=0;i<result.duals_.segments().size();++i) {
        const auto& segment=result.duals_.segments()[i];auto& flux=result.fluxes_[i];
        if(segment.ownerCell>=0) {
            if(domain.uniformLeaves()) {
            const int cell=segment.ownerCell,axis=segment.normalAxis;
            const auto lower=lowerFraction(mesh,cell,axis,rz);
            const auto& low=primal[cellFace(mesh,cell,axis,0)];const auto& high=primal[cellFace(mesh,cell,axis,1)];
            // Opposite half weights are required at the cell-centre plane:
            // J=(1-a)*F_low+a*F_high. Its divergence equals the mapped
            // primal mass divergence, also on nonuniform radial measures.
            flux={(1-lower)*low.positive+lower*high.positive,(1-lower)*low.negative+lower*high.negative};
            }
        } else {
            const auto face=domain.face(mesh,segment.primalNative);const auto& source=primal[segment.primalNative];
            if(face.area==0) require(source.positive==0 && source.negative==0,"nonzero mass flux through the RZ axis");
            const long double share=face.area>0?static_cast<long double>(segment.measure)/face.area:0;
            require(share>=0 && share<=1+64*std::numeric_limits<double>::epsilon(),"invalid primal-face partition");
            flux={share*source.positive,share*source.negative};
        }
        require(std::isfinite(flux.positive) && std::isfinite(flux.negative) && flux.positive>=0 && flux.negative>=0,
            "invalid mapped directional mass flux");
        for(int side=0;side<2;++side) {
            const int row=side?segment.positiveNative:segment.negativeNative;if(row<0) continue;
            net[row]+=(side?-1:1)*(flux.positive-flux.negative);gross[row]+=flux.positive+flux.negative;
        }
    }
    for(int row=0;row<count;++row) {
        const auto scale=result.before_[row]+result.after_[row]+gross[row];
        const auto residual=std::abs(result.after_[row]-result.before_[row]+net[row])/scale;
        require(std::isfinite(residual) && residual<=4096*std::numeric_limits<double>::epsilon(),
            "dual mass update is inconsistent with accepted phase fluxes");
        result.maxRelativeClosure_=std::max(result.maxRelativeClosure_,static_cast<double>(residual));
    }
    return result;
}
void HydroConservativeTransport::validateParametricForce(const InterfaceState& source,HydroForceTimeLevel level) const
{
    require(source.hasParametricGeometry(),"parameter force has another geometry representation");
    const bool stepEnd=level==HydroForceTimeLevel::StepEnd,midpoint=level==HydroForceTimeLevel::StepMidpoint;
    require(stepEnd || midpoint || level==HydroForceTimeLevel::StepBegin,"invalid parameter force layer");
    const auto& expected=midpoint?parameterMiddle_:(stepEnd?parameterAfter_:parameterBefore_);
    require(expected && source.parametricGeometry().geometryHandle()==expected,"parameter force is not from the mass step's actual geometry");
    const auto progress=source.progress();
    const double expectedTime=begin_+(midpoint?.5:(stepEnd?1:0))*dt_;
    require(progress.completedSteps==(stepEnd || midpoint?step_:step_-1) &&
        std::abs(progress.endTime-expectedTime)<=64*std::numeric_limits<double>::epsilon()*std::max({dt_,std::abs(progress.endTime),std::abs(expectedTime)}),
        "parameter force belongs to another physical time layer");
}
void HydroConservativeTransport::validateFor(const QuadtreeMesh& mesh,bool axisymmetric,const HydroStepTime& time,
    const CellField<double>& density) const
{
    require(bool(identity_),"transport has no retained mesh identity");
    identity_->validate(mesh);
    require(mesh.Nx>=2 && mesh.Ny>=2,"conservative NS requires a two-dimensional native background grid");
    require(axisymmetric==axisymmetric_ && time.index==step_ && time.begin==begin_ && time.dt==dt_,
        "Hydro and phase transport have different coordinates or physical steps");
    require(density.value && density.size()==static_cast<int>(sourceDensity_.size()),"source density layout changed");
    for(int i=0;i<density.size();++i) require(density.label(i)==mesh.cellIndexToLabel[i] && std::isfinite(density.value[i])
        && std::abs(density.value[i]-sourceDensity_[i])<=128*std::numeric_limits<double>::epsilon()*sourceDensity_[i],
        "Hydro source density differs from the conserved phase mixture");
}
HydroConservativeTransport::MomentumFlux HydroConservativeTransport::upwind(const std::vector<double>& velocity,
    const std::function<double(const Duals::Segment&,int)>& incoming) const
{
    require(velocity.size()==before_.size(),"wrong native velocity layout");
    for(double v:velocity) require(std::isfinite(v),"nonfinite native velocity");
    MomentumFlux result;result.outward.assign(before_.size(),0);
    for(std::size_t i=0;i<fluxes_.size();++i) {
        const auto& segment=duals_.segments()[i];const auto& flux=fluxes_[i];
        const int owner=segment.negativeNative>=0?segment.negativeNative:segment.positiveNative;
        const int component=owner<radialRows_?0:1;
        const auto donor=[&](int row,long double amount) {
            if(amount==0) return 0.;if(row>=0) return velocity[row];
            require(bool(incoming),"physical incoming momentum velocity is missing");
            const auto value=incoming(segment,component);require(std::isfinite(value),"nonfinite incoming momentum velocity");
            return value;
        };
        const auto value=flux.positive*donor(segment.negativeNative,flux.positive)
            -flux.negative*donor(segment.positiveNative,flux.negative);
        require(std::isfinite(value),"overflowed momentum flux");
        if(segment.negativeNative>=0) result.outward[segment.negativeNative]+=value;
        else result.boundaryOutward[component]-=value;
        if(segment.positiveNative>=0) result.outward[segment.positiveNative]-=value;
        else result.boundaryOutward[component]+=value;
    }
    for(auto value:result.outward) require(std::isfinite(value),"overflowed momentum divergence");
    for(auto value:result.boundaryOutward) require(std::isfinite(value),"overflowed boundary momentum balance");
    return result;
}

HydroConservativeTransport::MomentumFlux HydroConservativeTransport::centered(const std::vector<double>& velocity) const
{
    require(uniformLeaves_ && velocity.size()==before_.size(),"centered momentum requires uniform XY/RZ native rows");
    for(double u:velocity)require(std::isfinite(u),"nonfinite centered momentum velocity");
    MomentumFlux result;result.outward.assign(velocity.size(),0);
    for(std::size_t k=0;k<fluxes_.size();++k){const auto& face=duals_.segments()[k];const auto& f=fluxes_[k];
        if(face.negativeNative<0 || face.positiveNative<0){require(f.positive==0 && f.negative==0,"centered momentum lacks closed boundary mass flux");continue;}
        const long double momentum=(f.positive-f.negative)*(.5L*(velocity[face.negativeNative]+velocity[face.positiveNative]));
        result.outward[face.negativeNative]+=momentum;result.outward[face.positiveNative]-=momentum;}
    return result;
}

HydroConservativeTransport::MomentumRemap HydroConservativeTransport::remapMomentum(
    const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,const InterfaceConservativePhase& old,
    const std::vector<double>& oldVelocity,const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceConservativePhase& next,double liquidDensity,double gasDensity)
{
    old.validate(oldMesh,oldDomain);next.validate(mesh,domain);
    require(next.transferredFrom(old),"momentum remap requires the actual target inventory lineage");
    require(oldDomain.axisymmetric()==domain.axisymmetric() && oldMesh.Nx==mesh.Nx && oldMesh.Ny==mesh.Ny
        && oldMesh.minCoord.x==mesh.minCoord.x && oldMesh.minCoord.y==mesh.minCoord.y
        && oldMesh.maxCoord.x==mesh.maxCoord.x && oldMesh.maxCoord.y==mesh.maxCoord.y,
        "momentum remap requires unchanged background lattice and measure");
    require(oldDomain.cellIndexByRow().size()==static_cast<std::size_t>(oldMesh.cellIndexNum())
        && domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "momentum remap requires complete native domains");
    require(std::isfinite(liquidDensity) && liquidDensity>0 && std::isfinite(gasDensity) && gasDensity>0,
        "momentum remap phase densities must be finite and positive");
    require(oldVelocity.size()==static_cast<std::size_t>(oldDomain.nativeFaceCount()),"source velocity layout mismatch");
    for(double u:oldVelocity)require(std::isfinite(u),"nonfinite source momentum velocity");
    const bool rz=domain.axisymmetric();const auto oldDuals=Duals::build(oldMesh,oldDomain),duals=Duals::build(mesh,domain);
    const auto mass=[&](const InterfaceConservativePhase& phase,int c) {
        const auto liquid=phase.liquidVolumes()[c];
        const auto value=liquidDensity*liquid+gasDensity*(phase.capacities()[c]-liquid);
        require(std::isfinite(value) && value>0,"invalid remapped mixture mass");return value;
    };
    const auto box=[](const QuadtreeMesh& m,int c) {
        const auto& cell=m.cells(m.cellIndexToLabel[c]);
        return std::pair<Vector2D,Vector2D>{m.verts(cell.cellToVert[0]).coord(),m.verts(cell.cellToVert[2]).coord()};
    };
    const auto overlap=[](Vector2D a,Vector2D b,Vector2D c,Vector2D d) {
        return std::pair<Vector2D,Vector2D>{{std::max(a.x,c.x),std::max(a.y,c.y)},
            {std::min(b.x,d.x),std::min(b.y,d.y)}};
    };
    std::vector<Pieces> oldPieces(oldMesh.cellIndexNum()),pieces(mesh.cellIndexNum());
    std::vector<long double> oldMass(oldDomain.nativeFaceCount(),0),consumed(oldDomain.nativeFaceCount(),0);
    MomentumRemap result;result.mass.assign(domain.nativeFaceCount(),0);result.momentum.assign(domain.nativeFaceCount(),0);
    result.velocity.resize(domain.nativeFaceCount());
    std::vector<long double> supplied(domain.nativeFaceCount(),0);std::vector<int> exclusive(domain.nativeFaceCount(),-2);
    for(int c:oldDomain.cellIndexByRow()) {
        oldPieces[c]=cellPieces(oldMesh,c,oldDuals,rz);
        for(const auto& component:oldPieces[c])for(const auto& p:component)oldMass[p.native]+=p.fraction*mass(old,c);
    }
    for(int c:domain.cellIndexByRow()) {
        pieces[c]=cellPieces(mesh,c,duals,rz);
        for(const auto& component:pieces[c])for(const auto& p:component)result.mass[p.native]+=p.fraction*mass(next,c);
    }
    const CartesianLeafLocator locator(mesh);
    std::vector<int> owners(mesh.cellIndexNum(),0);
    const auto transfer=[&](int source,int target,long double amount) {
        require(std::isfinite(amount) && amount>=0,"invalid local momentum transfer weight");
        if(amount==0)return;
        consumed[source]+=amount;supplied[target]+=amount;result.momentum[target]+=amount*oldVelocity[source];
        if(exclusive[target]==-2)exclusive[target]=source;else if(exclusive[target]!=source)exclusive[target]=-1;
    };
    for(int c:oldDomain.cellIndexByRow()) {
        const auto parent=box(oldMesh,c);std::vector<int> children;long double sumMass=0;
        for(int target:locator.intersectingBox(parent.first,parent.second)) {
            const auto child=box(mesh,target),cut=overlap(parent.first,parent.second,child.first,child.second);
            if(cut.second.x<=cut.first.x || cut.second.y<=cut.first.y)continue;
            require(child.first.x>=parent.first.x && child.first.y>=parent.first.y
                && child.second.x<=parent.second.x && child.second.y<=parent.second.y,
                "conservative momentum remap currently requires refinement, not coarsening");
            ++owners[target];children.push_back(target);sumMass+=mass(next,target);
        }
        require(!children.empty(),"momentum remap lost a source parent");
        const auto parentMass=mass(old,c);
        require(std::abs(sumMass-parentMass)<=1024*std::numeric_limits<double>::epsilon()*parentMass,
            "target phase masses do not conserve the actual source parent");
        if(children.size()>1) {
            ++result.refinedParents;
            for(int target:children)for(int component=0;component<2;++component)
                for(const auto& t:pieces[target][component])for(const auto& s:oldPieces[c][component])
                    transfer(s.native,t.native,s.fraction*t.fraction*mass(next,target));
        } else {
            ++result.retainedParents;const int target=children.front();
            const auto capacity=measure(parent.first,parent.second,rz);
            for(int component=0;component<2;++component)for(const auto& t:pieces[target][component])
                for(const auto& s:oldPieces[c][component]) {
                    const auto cut=overlap(s.lo,s.hi,t.lo,t.hi);
                    if(cut.second.x<=cut.first.x || cut.second.y<=cut.first.y)continue;
                    // Retain the original arithmetic when the two pieces are
                    // identical; otherwise split only at the new native face.
                    const auto fraction=cut.first.x==s.lo.x && cut.first.y==s.lo.y
                        && cut.second.x==s.hi.x && cut.second.y==s.hi.y?s.fraction:measure(cut.first,cut.second,rz)/capacity;
                    transfer(s.native,t.native,parentMass*fraction);
                }
        }
    }
    for(int n:owners)require(n==1,"target momentum cell has missing or multiple source parents");
    const auto closure=[&](long double actual,long double expected) {
        require(expected>0 && std::isfinite(actual),"invalid momentum transfer mass closure");
        const auto residual=std::abs(actual-expected)/expected;
        require(residual<=2048*std::numeric_limits<double>::epsilon(),"local momentum transfer mass marginal failed");
        result.maxRelativeMassClosure=std::max(result.maxRelativeMassClosure,static_cast<double>(residual));
    };
    long double absoluteMomentum[2]{};
    for(int row=0;row<oldDomain.nativeFaceCount();++row) {
        closure(consumed[row],oldMass[row]);const int component=row<oldMesh.vEdgeIndexNum()?0:1;
        const auto momentum=oldMass[row]*oldVelocity[row];result.sourceMomentum[component]+=momentum;
        absoluteMomentum[component]+=std::abs(momentum);
        result.sourceKineticEnergy[component]+=.5L*momentum*oldVelocity[row];
    }
    for(int row=0;row<domain.nativeFaceCount();++row) {
        closure(supplied[row],result.mass[row]);const int component=row<mesh.vEdgeIndexNum()?0:1;
        result.velocity[row]=exclusive[row]>=0?oldVelocity[exclusive[row]]:static_cast<double>(result.momentum[row]/result.mass[row]);
        require(std::isfinite(result.velocity[row]),"nonfinite remapped momentum velocity");
        // Report the actually returned binary64 velocity, including its final
        // rounding; extensive momentum is never a separate acceptance field.
        result.momentum[row]=result.mass[row]*result.velocity[row];
        result.targetMomentum[component]+=result.momentum[row];
        result.targetKineticEnergy[component]+=.5L*result.momentum[row]*result.velocity[row];
    }
    for(int component=0;component<2;++component) {
        result.sourceAbsoluteMomentum[component]=absoluteMomentum[component];
        const auto difference=std::abs(result.targetMomentum[component]-result.sourceMomentum[component]);
        const auto scale=absoluteMomentum[component];
        require(difference<=4096*std::numeric_limits<double>::epsilon()*scale,"remapped component momentum did not close");
        if(scale>0)result.maxRelativeMomentumClosure=std::max(result.maxRelativeMomentumClosure,static_cast<double>(difference/scale));
        require(result.targetKineticEnergy[component]<=result.sourceKineticEnergy[component]*(1+4096*std::numeric_limits<double>::epsilon()),
            "positive momentum remap increased kinetic energy beyond roundoff");
    }
    return result;
}
}
