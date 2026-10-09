#include "Modules/Interface/Compatibility/InterfacePlicTransport.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Core/Geometry/ExactPointKey.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <map>
#include <sstream>
#include <iomanip>

namespace InterfacePlicTransport
{
namespace
{
namespace Cut=Quantum::Geometry::PlanarCellCut;
void require(bool value,const char* message)
{ if (!value) throw std::invalid_argument(std::string("[InterfacePlicTransport] ")+message); }
namespace Curve=Quantum::Geometry::QuadraticCellCut;
struct Donor { Cut::Rectangle box; Cut::Reconstruction plane; long double fraction=0; bool uniform=false;
    std::optional<Curve::Reconstruction> curved; };
Donor reconstruct(Cut::Rectangle box,CellData data,bool rz)
{
    const bool uniform=std::min(data.fraction,1-data.fraction)<=Cut::relativeVolumeResolution;
    // Below geometric inverse precision there is no resolved plane direction.
    // Transport the retained cell average over the same donor overlap instead.
    // reconstruct still validates fraction/box; its dummy plane is NOT used
    // for a uniform donor's liquid subvolume.
    const auto normal=uniform?Vector2D{1,0}:data.normal;
    Donor result{box,Cut::reconstruct(box,normal,data.fraction,rz),data.fraction,uniform};
    if (!uniform && data.curved) result.curved=Curve::reconstruct(box,*data.curved,data.fraction,rz);
    return result;
}
Vector2D facePoint(Vector2D begin,Vector2D midpoint,Vector2D end,double fraction)
{
    if (fraction==0) return begin;
    if (fraction==.5) return midpoint;
    if (fraction==1) return end;
    return fraction<.5?begin+(midpoint-begin)*(2*fraction):midpoint+(end-midpoint)*(2*fraction-1);
}
}
std::vector<std::vector<FaceSweep>> traceSweeps(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceStepTime& time,const InterfaceMotionVelocityLayers& velocities,
    const InterfaceNormalFluxView& previousFlux,const InterfaceNormalFluxView& currentFlux,
    const std::vector<std::vector<double>>& partitions,TraceReport* report)
{
    domain.validateMeshIdentity(mesh);
    require(domain.uniformLeaves() && domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "velocity sweeps require the complete uniform rectangle");
    require(time.index>0 && std::isfinite(time.begin) && std::isfinite(time.end) && time.end>time.begin
        && std::isfinite(time.dt) && time.dt>0 && std::abs(time.end-time.begin-time.dt)
        <=16*std::numeric_limits<double>::epsilon()*std::max({time.dt,std::abs(time.begin),std::abs(time.end)}),
        "invalid external sweep time interval");
    require(velocities.current && velocities.previous,"both immutable velocity layers are required");
    require(partitions.empty() || partitions.size()==static_cast<std::size_t>(domain.nativeFaceCount()),
        "sweep partitions have wrong native layout");
    const LevelSetBoundaryTrace::Domain bounds{mesh.minCoord,mesh.maxCoord,domain.spacing(),
        domain.axisymmetric() && mesh.minCoord.x==0};
    std::map<Quantum::Geometry::ExactPointKey,Vector2D> displacements;
    const auto displacement=[&](Vector2D point) {
        const auto key=Quantum::Geometry::exactPointKey(point);const auto found=displacements.find(key);
        if (found!=displacements.end()) return found->second;
        const auto traced=LevelSetBoundaryTrace::trace(bounds,point,time.dt,velocities.current,velocities.previous);
        require(!traced.boundaryHit,"PLIC backtrace needs an explicit physical inflow/exterior trajectory closure");
        const auto value=point-traced.point;displacements.emplace(key,value);return value;
    };
    const auto same=[](Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;};
    TraceReport quality;
    std::vector<std::vector<FaceSweep>> result(domain.nativeFaceCount());
    for (const auto face:domain.faces(mesh)) {
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if (face.horizontal?a.x>b.x:a.y>b.y) std::swap(a,b);
        for (const auto* view:{&previousFlux,&currentFlux}) {
            const auto chart=view->faceChart(face.nativeRow);
            require(same(chart.begin,a) && same(chart.midpoint,face.center) && same(chart.end,b),
                "velocity flux chart differs from its actual native face");
        }
        auto cuts=partitions.empty()?std::vector<double>{0,.5,1}:partitions[face.nativeRow];
        if (partitions.empty()) for (bool current:{false,true}) {
            const auto& view=current?currentFlux:previousFlux;
            if (!view.hasNormalFluxRoots()) continue;
            for (int half:{0,1}) {
                const auto roots=view.normalFluxRoots(face.nativeRow,half);
                if (roots.identicallyZero) ++quality.zeroTraceHalves;
                double previousPosition=-1,previousUpper=-1;
                require(!roots.identicallyZero || roots.roots.empty(),"zero velocity trace has isolated roots");
                for (const auto& root:roots.roots) {
                    require(std::isfinite(root.point) && root.lower>=0 && root.upper<=1
                        && root.lower<=root.point && root.point<=root.upper,"invalid velocity root enclosure");
                    require(root.lower>previousUpper,"velocity roots have overlapping or unordered enclosures");
                    previousUpper=root.upper;
                    const double position=.5*(half+root.point);
                    require(position>previousPosition,
                        "distinct velocity roots collapse in the shared face chart; use an endpoint-relative subdivision");
                    previousPosition=position;
                    quality.roots.push_back({face.nativeRow,half,current,root});
                    if (root.point==0 || root.point==1) continue;
                    require(position>.5*half && position<.5*(half+1),
                        "velocity root is unresolved in the shared face chart; use an endpoint-relative subdivision");
                    cuts.push_back(position);
                }
            }
        }
        if (partitions.empty()) {
            std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
            quality.automaticInteriorCuts+=cuts.size()-3;
        }
        require(cuts.size()>=3 && cuts.front()==0 && cuts.back()==1
            && std::find(cuts.begin(),cuts.end(),.5)!=cuts.end(),"sweep partition must contain 0, 1/2 and 1");
        result[face.nativeRow].reserve(cuts.size()-1);
        for (std::size_t i=1;i<cuts.size();++i) {
            const double lo=cuts[i-1],hi=cuts[i];
            require(std::isfinite(lo) && std::isfinite(hi) && lo>=0 && hi<=1 && hi>lo,"invalid sweep partition");
            const int half=lo>=.5?1:0;const double first=2*lo-half,last=2*hi-half;
            const auto old=previousFlux.integratedNormalFlux(face.nativeRow,half,first,last);
            const auto current=currentFlux.integratedNormalFlux(face.nativeRow,half,first,last);
            require(std::isfinite(old) && std::isfinite(current) && (old==0 || current==0 || (old>0)==(current>0)),
                "temporal face-flux reversal requires physical substeps");
            for (double fraction:{0.,.25,.5,.75,1.}) {
                const long double t=first+static_cast<long double>(last-first)*fraction;
                const auto u=previousFlux.normalFluxDensity(face.nativeRow,half,false,t);
                const auto v=currentFlux.normalFluxDensity(face.nativeRow,half,false,t);
                require(std::isfinite(u) && std::isfinite(v) && (u==0 || v==0 || (u>0)==(v>0)),
                    "sampled temporal trace reversal requires physical substeps");
            }
            const long double volume=(old+current)*(.5L*time.dt);
            require(std::isfinite(volume),"time-integrated face volume overflow");
            result[face.nativeRow].push_back({lo,hi,volume,displacement(facePoint(a,face.center,b,lo)),
                displacement(facePoint(a,face.center,b,hi))});
        }
    }
    if (report) *report=std::move(quality);
    return result;
}
static std::vector<InterfaceConservativePhase::FaceFlux> sweepImpl(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& inventory,
    const std::vector<Vector2D>& normals,const std::vector<std::vector<FaceSweep>>& faceSweeps,
    const ExteriorCell& exterior,bool correct,SweepReport* report,const CurvedDonors& curved)
{
    inventory.validate(mesh,domain);
    require(domain.uniformLeaves() && domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),
        "sweeping requires the complete uniform Cartesian rectangle");
    require(normals.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"normal array has wrong native layout");
    require(curved.empty() || curved.size()==normals.size(),"curved donor array has wrong native layout");
    const auto h=domain.spacing();
    require(faceSweeps.size()==static_cast<std::size_t>(domain.nativeFaceCount()),"sweep array has wrong native layout");
    const bool rz=domain.axisymmetric();
    const UniformLeafGrid grid(mesh,CartesianGridScope::UniformLeaves);
    std::vector<Donor> donors(mesh.cellIndexNum());
    std::vector<std::array<int,2>> logical(mesh.cellIndexNum());
    for (int j=0;j<grid.ny();++j) for (int i=0;i<grid.nx();++i) {
        const auto& cell=mesh.cells(grid.cells()[j*grid.nx()+i]);
        const int c=cell.index();logical[c]={{i,j}};
        const Cut::Rectangle box{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
        donors[c]=reconstruct(box,{inventory.liquidVolumes()[c]/inventory.capacities()[c],normals[c],
            curved.empty()?std::nullopt:curved[c]},rz);
    }
    std::vector<InterfaceConservativePhase::FaceFlux> result(domain.nativeFaceCount());
    SweepReport quality;
    for (const auto face:domain.faces(mesh)) {
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto first=mesh.verts(edge.edgeToVert[0]).coord(),last=mesh.verts(edge.edgeToVert[1]).coord();
        if (face.horizontal?first.x>last.x:first.y>last.y) std::swap(first,last);
        const auto coordinate=[&](double t) {return facePoint(first,face.center,last,t);};
        InterfaceConservativePhase::DirectionalFlux directions;bool correctedFace=false;
        double lastEnd=0;Vector2D lastDisplacement;
        const auto& pieces=faceSweeps[face.nativeRow];
        for (const auto& piece:pieces) {
            require(std::isfinite(piece.begin) && std::isfinite(piece.end) && piece.begin==lastEnd
                && piece.end>piece.begin && piece.end<=1,"subface sweeps must partition the face without gaps/overlap");
            require(piece.begin==0 || (piece.beginDisplacement.x==lastDisplacement.x
                && piece.beginDisplacement.y==lastDisplacement.y),"shared subface endpoint has two different trajectories");
            lastDisplacement=piece.endDisplacement;
            lastEnd=piece.end;
            const long double total=piece.totalVolume,magnitude=std::abs(total);
            require(std::isfinite(total),"nonfinite prescribed face volume");
            const auto a=coordinate(piece.begin),b=coordinate(piece.end);
            require(face.horizontal?a.x<b.x:a.y<b.y,"unresolved physical subface interval");
            const auto da=piece.beginDisplacement,db=piece.endDisplacement;
            for (const auto d:{da,db}) require(std::isfinite(d.x) && std::isfinite(d.y)
                && std::abs(d.x)/h.x+std::abs(d.y)/h.y<=1,"sweep endpoint CFL exceeds one or is nonfinite");
            const double na=face.horizontal?da.y:da.x,nb=face.horizontal?db.y:db.x;
            if (total==0) {
                require(na==0 && nb==0,"zero net face volume has unresolved normal endpoint motion");
                continue;
            }
            if (!(face.area>0 && total*na>=0 && total*nb>=0 && (na!=0 || nb!=0))) {
                std::ostringstream message;
                message<<std::setprecision(17)<<"[InterfacePlicTransport] sweep normal directions disagree with prescribed flux or contain spatial counterflow"
                    <<"; face="<<face.nativeRow<<", interval=["<<piece.begin<<','<<piece.end<<"]"
                    <<", flux="<<total<<", normal displacements=("<<na<<','<<nb<<')';
                throw std::invalid_argument(message.str());
            }
            // Work in local coordinates throughout; scaling a displacement does
            // not round a new global departure point back to the mesh lattice.
            const Cut::Point origin{a.x,a.y};
            const auto polygonAt=[&](long double scale) {
                const Cut::Point p{0,0};
                const Cut::Point q{static_cast<long double>(b.x)-a.x,static_cast<long double>(b.y)-a.y};
                return std::vector<Cut::Point>{p,q,{q.x-db.x*(face.horizontal?1:scale),q.y-db.y*(face.horizontal?scale:1)},
                    {p.x-da.x*(face.horizontal?1:scale),p.y-da.y*(face.horizontal?scale:1)}};
            };
            long double scale=1;
            if (correct) {
                const Cut::Rectangle support{{rz?std::max(0.,std::min(a.x,b.x)-h.x):std::min(a.x,b.x)-h.x,
                    std::min(a.y,b.y)-h.y},{std::max(a.x,b.x)+h.x,std::max(a.y,b.y)+h.y}};
                const auto measure=[&](long double factor) {
                    return Cut::overlapInFrame(support,polygonAt(factor),origin,{1,0},0,rz).total;
                };
                const long double raw=measure(1);
                require(raw>0,"zero swept measure for nonzero face flux");
                const auto discrepancy=std::abs(raw-magnitude)/magnitude;
                require(discrepancy<=.25L,"swept geometry and prescribed volume differ by more than 25 percent");
                quality.maxRelativeVolumeCorrection=std::max(quality.maxRelativeVolumeCorrection,discrepancy);
                long double upper=std::numeric_limits<long double>::infinity();
                for (const auto d:{da,db}) {
                    const double normal=face.horizontal?d.y:d.x,tangent=face.horizontal?d.x:d.y;
                    if (normal!=0) upper=std::min(upper,(1-std::abs(static_cast<long double>(tangent))/(face.horizontal?h.x:h.y))
                        *(face.horizontal?h.y:h.x)/std::abs(normal));
                    if (rz && !face.horizontal && d.x>0) upper=std::min(upper,static_cast<long double>(a.x)/d.x);
                }
                require(upper>0 && std::isfinite(upper) && measure(upper)>=magnitude,
                    "volume correction needs a larger stencil/CFL or crosses the RZ axis");
                if (!rz) scale=magnitude/raw;
                else {
                    long double lo=0,hi=upper;
                    for (int iteration=0;iteration<96;++iteration) {
                        scale=lo+(hi-lo)/2;const auto value=measure(scale);
                        if (value==magnitude || scale==lo || scale==hi) break;
                        if (value<magnitude) lo=scale;else hi=scale;
                    }
                }
                require(scale>0 && scale<=upper && std::abs(measure(scale)-magnitude)
                    <=256*std::numeric_limits<double>::epsilon()*magnitude,"unresolved swept-volume correction");
                quality.maxNormalScaleChange=std::max(quality.maxNormalScaleChange,std::abs(scale-1));
                if (discrepancy>256*std::numeric_limits<double>::epsilon()) correctedFace=true;
            }
            ++quality.sweptPieces;
            const int owner=face.negativeCell>=0?face.negativeCell:face.positiveCell;
            const auto ij=logical[owner];
            long double swept=0,liquid=0;
            const auto polygon=polygonAt(scale);
            // With outflow CFL<=1 the true face preimage meets only this 3x3
            // donor neighbourhood (including one actual exterior ghost layer).
            for (int j=std::max(-1,ij[1]-1);j<=std::min(grid.ny(),ij[1]+1);++j)
                for (int i=std::max(-1,ij[0]-1);i<=std::min(grid.nx(),ij[0]+1);++i) {
                    const bool inside=i>=0 && j>=0 && i<grid.nx() && j<grid.ny();
                    Donor donor;
                    if (inside) donor=donors[mesh.cells(grid.cells()[j*grid.nx()+i]).index()];
                    else {
                        // Continue the nearest actual boundary rectangle by its
                        // own width. No independent regenerated nominal grid.
                        const int ni=std::clamp(i,0,grid.nx()-1),nj=std::clamp(j,0,grid.ny()-1);
                        auto box=donors[mesh.cells(grid.cells()[nj*grid.nx()+ni]).index()].box;
                        const auto size=box.upper-box.lower;
                        box.lower=box.lower+Vector2D{(i-ni)*size.x,(j-nj)*size.y};
                        box.upper=box.upper+Vector2D{(i-ni)*size.x,(j-nj)*size.y};
                        // Physical RZ sweeps cannot use negative-radius ghosts.
                        if (rz && box.lower.x<0) continue;
                        donor.box=box;
                    }
                    const auto& box=donor.box;
                    if (!inside) {
                        const auto cap=Cut::overlapInFrame(box,polygon,origin,{1,0},0,rz).total;
                        if (cap==0) continue;
                        require(static_cast<bool>(exterior),"swept inflow needs explicit exterior phase geometry");
                        donor=reconstruct(box,exterior(box),rz);
                    }
                    const auto part=donor.curved?Curve::overlapInFrame(box,polygon,origin,*donor.curved,rz):
                        Cut::overlapInFrame(box,polygon,origin,donor.plane.normal,donor.plane.offset,rz);
                    swept+=part.total;
                    liquid+=donor.uniform?donor.fraction*part.total:part.liquid;
                }
            if (!(std::isfinite(swept) && swept>0 && std::abs(swept-magnitude)<=
                256*std::numeric_limits<double>::epsilon()*magnitude)) {
                std::ostringstream message;
                message<<std::setprecision(17)<<"[InterfacePlicTransport] swept polygon does not cover the face volume"
                    <<"; face="<<face.nativeRow<<", interval=["<<piece.begin<<','<<piece.end<<"]"
                    <<", target="<<magnitude<<", donor sum="<<swept<<", normal scale="<<scale;
                throw std::invalid_argument(message.str());
            }
            require(std::isfinite(liquid) && liquid>=0 && liquid<=swept,"invalid swept liquid subvolume");
            // Match the retained face metric (double) to the same geometric ratio.
            // Only the roundoff-size metric discrepancy admitted above is allowed.
            const auto phase=magnitude*(liquid/swept);
            if (total>0) {directions.totalPositive+=magnitude;directions.liquidPositive+=phase;}
            else {directions.totalNegative+=magnitude;directions.liquidNegative+=phase;}
        }
        require(pieces.empty() || lastEnd==1,"subface sweeps leave an uncovered face interval");
        if (directions.totalPositive>0 || directions.totalNegative>0) ++quality.sweptFaces;
        if (correctedFace) ++quality.correctedFaces;
        auto& flux=result[face.nativeRow];
        if (directions.totalPositive>0 && directions.totalNegative>0)
            flux=InterfaceConservativePhase::FaceFlux::fromDirectional(directions);
        else {
            flux.totalVolume=directions.totalPositive-directions.totalNegative;
            flux.highOrderLiquidVolume=directions.liquidPositive-directions.liquidNegative;
        }
        if (face.negativeCell<0 && directions.totalPositive>0)
            flux.boundaryLiquidFraction=directions.liquidPositive/directions.totalPositive;
        if (face.positiveCell<0 && directions.totalNegative>0)
            flux.boundaryLiquidFraction=directions.liquidNegative/directions.totalNegative;
    }
    if (report) *report=quality;
    return result;
}
std::vector<InterfaceConservativePhase::FaceFlux> swept(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& inventory,const std::vector<Vector2D>& normals,
    const std::vector<long double>& faceVolumes,const std::vector<Vector2D>& vertexDisplacements,
    const ExteriorCell& exterior,SweepReport* report,const CurvedDonors& curved)
{
    domain.validateMeshIdentity(mesh);
    require(faceVolumes.size()==static_cast<std::size_t>(domain.nativeFaceCount())
        && vertexDisplacements.size()==static_cast<std::size_t>(mesh.vertIndexNum()),"sweep arrays have wrong native layout");
    std::vector<std::vector<FaceSweep>> pieces(domain.nativeFaceCount());
    for (const auto face:domain.faces(mesh)) {
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        int a=mesh.verts(edge.edgeToVert[0]).index(),b=mesh.verts(edge.edgeToVert[1]).index();
        const auto p=mesh.verts(edge.edgeToVert[0]).coord(),q=mesh.verts(edge.edgeToVert[1]).coord();
        if (face.horizontal?p.x>q.x:p.y>q.y) std::swap(a,b);
        pieces[face.nativeRow].push_back({0,1,faceVolumes[face.nativeRow],vertexDisplacements[a],vertexDisplacements[b]});
    }
    return sweepImpl(mesh,domain,inventory,normals,pieces,exterior,true,report,curved);
}
std::vector<InterfaceConservativePhase::FaceFlux> sweptSubfaces(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& inventory,const std::vector<Vector2D>& normals,
    const std::vector<std::vector<FaceSweep>>& faceSweeps,const ExteriorCell& exterior,SweepReport* report,const CurvedDonors& curved)
{ return sweepImpl(mesh,domain,inventory,normals,faceSweeps,exterior,true,report,curved); }
std::vector<InterfaceConservativePhase::FaceFlux> translation(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& inventory,
    const std::vector<Vector2D>& normals,Vector2D displacement,const ExteriorCell& exterior,const CurvedDonors& curved)
{
    require(!domain.axisymmetric() || displacement.x==0,"constant radial translation is not RZ incompressible");
    domain.validateMeshIdentity(mesh);
    std::vector<std::vector<FaceSweep>> pieces(domain.nativeFaceCount());
    for (const auto face:domain.faces(mesh)) pieces[face.nativeRow].push_back({0,1,
        static_cast<long double>(face.horizontal?displacement.y:displacement.x)*face.area,displacement,displacement});
    return sweepImpl(mesh,domain,inventory,normals,pieces,exterior,false,nullptr,curved);
}
}
