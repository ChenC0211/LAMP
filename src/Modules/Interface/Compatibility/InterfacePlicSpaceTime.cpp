#include "Modules/Interface/Compatibility/InterfacePlicTransport.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Core/Geometry/ExactPointKey.h"
#include "Core/Geometry/QuadraticChordEvents.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/NativeLeafCoverage.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include <exception>
#include <mutex>

namespace InterfacePlicTransport
{
namespace
{
namespace Cut=Quantum::Geometry::PlanarCellCut;
namespace Curve=Quantum::Geometry::QuadraticCellCut;
using Rates=std::array<long double,5>; // +/- total, +/- liquid, spatial error estimate
void require(bool ok,const char* reason)
{ if (!ok) throw std::invalid_argument(std::string("[InterfacePlicSpaceTime] ")+reason); }
Rates add(Rates a,const Rates& b) {for (std::size_t i=0;i<a.size();++i) a[i]+=b[i];return a;}
Rates scale(Rates a,long double w) {for (auto& x:a) x*=w;return a;}
struct Donor {Cut::Rectangle box;Cut::Reconstruction plane;long double fraction;bool uniform;
    std::optional<Curve::Reconstruction> curved;};
struct Range {double begin,end;long double weight=1;};
// Intersect a LINE parameter with a half-space. No point-phase sampling or
// polygon area subtraction: shared rectangle boundaries cut the same chord.
bool clip(long double value,long double slope,long double bound,long double& lo,long double& hi)
{
    if (slope==0) return value<=bound;
    const auto crossing=(bound-value)/slope;
    if (slope>0) hi=std::min(hi,crossing);else lo=std::max(lo,crossing);
    return lo<hi;
}
}
std::vector<InterfaceConservativePhase::FaceFlux> spaceTime(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const InterfaceConservativePhase& inventory,
    const std::vector<Vector2D>& normals,const InterfaceStepTime& time,
    const InterfaceMotionVelocityLayers& velocities,const InterfaceNormalFluxView& previousFlux,
    const InterfaceNormalFluxView& currentFlux,const SpaceTimeOptions& options,SpaceTimeReport* report,const CurvedDonors& curved,
    const SinglePhasePorts& ports)
{
    inventory.validate(mesh,domain);
    require(normals.size()==static_cast<std::size_t>(mesh.cellIndexNum()),"wrong normal layout");
    require(curved.empty() || curved.size()==normals.size(),"wrong curved donor layout");
    require(time.index>0 && std::isfinite(time.dt) && time.dt>0 && std::isfinite(time.begin) && std::isfinite(time.end)
        && time.end>time.begin && std::abs(time.end-time.begin-time.dt)<=16*std::numeric_limits<double>::epsilon()
        *std::max({time.dt,std::abs(time.begin),std::abs(time.end)}),"invalid physical step time");
    require(velocities.current && velocities.previous && previousFlux.canInterpolateRoots(currentFlux),
        "physical velocity layers and a shared flux interpolation/root authority required");
    require(!options.parallelFaces || velocities.independentQueries,
        "parallel face transport requires independent velocity queries");
    require(std::isfinite(options.relativeFluxTolerance) && options.relativeFluxTolerance>0
        && options.relativeFluxTolerance<1 && std::isfinite(options.mapToleranceCells) && options.mapToleranceCells>0
        && options.mapToleranceCells<1 && options.maxSpaceDepth>=1 && options.maxSpaceDepth<=24
        && options.maxSpaceLeaves>=2 && options.maxSpaceLeaves<=65536
        && options.maxTimeDepth>=1 && options.maxTimeDepth<=20,"invalid bounded quadrature options");
    const NativeLeafCoverage coverage(mesh,domain);
    Vector2D spacing{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
    for(int c:domain.cellIndexByRow()) {
        const auto size=mesh.cells(mesh.cellIndexToLabel[c]).size();
        spacing.x=std::min(spacing.x,size.x);spacing.y=std::min(spacing.y,size.y);
    }
    const double h=std::min(spacing.x,spacing.y);
    require(std::isfinite(h) && h>0,"empty or invalid selected leaf spacing");
    LevelSetBoundaryTrace::Domain bounds{mesh.minCoord,mesh.maxCoord,spacing,domain.axisymmetric() && mesh.minCoord.x==0};
    if (domain.cellIndexByRow().size()!=static_cast<std::size_t>(mesh.cellIndexNum()))
        bounds.activeSegment=[&](Vector2D a,Vector2D b){return coverage.covers(a,b);};
    std::vector<Donor> donors(mesh.cellIndexNum());
    for (int c:domain.cellIndexByRow()) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
        const Cut::Rectangle box{mesh.verts(cell.cellToVert[0]).coord(),mesh.verts(cell.cellToVert[2]).coord()};
        const auto f=inventory.liquidVolumes()[c]/inventory.capacities()[c];
        const bool uniform=std::min(f,1-f)<=Cut::relativeVolumeResolution;
        donors[c]={box,Cut::reconstruct(box,uniform?Vector2D{1,0}:normals[c],f,domain.axisymmetric()),f,uniform};
        if (!uniform && !curved.empty() && curved[c])
            donors[c].curved=Curve::reconstruct(box,*curved[c],f,domain.axisymmetric());
    }
    struct Port {int phase;Boundary side;bool horizontal;Vector2D lower,upper;};
    std::vector<Port> portFaces;std::vector<int> portPhase(domain.nativeFaceCount(),0);
    if (ports) for (const auto face:domain.faces(mesh)) {
        const int phase=ports(face.nativeRow);
        require(phase>=-1 && phase<=1,"invalid single-phase port declaration");
        if (phase==0) continue;
        require((face.negativeCell<0)!=(face.positiveCell<0),"single-phase port is not an exterior face");
        require(face.area>0,"zero-area axis is not a physical port");
        const int c=face.negativeCell<0?face.positiveCell:face.negativeCell;
        require(inventory.liquidVolumes()[c]==(phase<0?inventory.capacities()[c]:0),
            "single-phase port has no matching pure donor layer; refine or use a phase-history closure");
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if (face.horizontal?a.x>b.x:a.y>b.y) std::swap(a,b);
        require(face.horizontal?(a.y==mesh.minCoord.y || a.y==mesh.maxCoord.y):
            (a.x==mesh.minCoord.x || a.x==mesh.maxCoord.x),"internal solid boundary cannot be declared an outer port");
        const auto side=face.horizontal?(face.negativeCell<0?Boundary::south:Boundary::north):
            (face.negativeCell<0?Boundary::west:Boundary::east);
        portPhase[face.nativeRow]=phase;portFaces.push_back({phase,side,face.horizontal,a,b});
    }
    SpaceTimeReport quality;std::vector<InterfaceConservativePhase::FaceFlux> result(domain.nativeFaceCount());
    // A donor's corners and interface/edge intersections depend only on this
    // call's already reconstructed OLD inventory. Neighboring half faces and
    // their trace arrangements may reuse those exact vertices. Keep the
    // original ordering/arithmetic and build only when an arrangement asks;
    // time-dependent crossings/tangencies still run on each actual trace.
    std::vector<std::vector<Cut::Point>> donorVertices(donors.size());
    std::vector<std::once_flag> donorVerticesReady(donors.size());
    const auto arrangementVertices=[&](int cell,SpaceTimeReport& quality) -> const std::vector<Cut::Point>& {
        ++quality.arrangementVertexRequests;
        std::call_once(donorVerticesReady[cell],[&] {
        // Publish a complete immutable donor only. A failed construction is
        // retriable by call_once and must leave no partially populated cache.
        std::vector<Cut::Point> vertices;
        ++quality.arrangementVertexBuilds;
        const auto& d=donors[cell];const auto& box=d.box;
        const std::array<Cut::Point,4> corners{{{box.lower.x,box.lower.y},{box.upper.x,box.lower.y},
            {box.upper.x,box.upper.y},{box.lower.x,box.upper.y}}};
        vertices.assign(corners.begin(),corners.end());
        if(d.curved) for(int edge=0;edge<4;++edge) {
            const auto x=corners[edge],y=corners[(edge+1)%4];
            // Retain isolated gas-edge contacts as arrangement vertices too.
            ++quality.arrangementEdgeRootQueries;
            for(const auto root:Quantum::Geometry::QuadraticChordEvents::intersections(box,*d.curved,x,y).roots) {
                const long double t=root.point;
                if(t>0 && t<1) vertices.push_back({std::fma(t,y.x-x.x,x.x),std::fma(t,y.y-x.y,x.y)});
            }
        } else if(!d.uniform) {
            const auto ends=Cut::interfaceSegment(box,d.plane.normal,d.plane.offset);
            const long double cx=(static_cast<long double>(box.lower.x)+box.upper.x)/2;
            const long double cy=(static_cast<long double>(box.lower.y)+box.upper.y)/2;
            if(ends) for(auto p:*ends) vertices.push_back({p.x+cx,p.y+cy});
        }
        donorVertices[cell]=std::move(vertices);
        });
        return donorVertices[cell];
    };
    const auto timeRoots=previousFlux.interpolationTimeRoots(currentFlux);
    quality.temporalRootEnclosures=timeRoots.roots.size();
    for (const auto& root:timeRoots.roots)
        quality.maxTemporalRootEnclosureWidth=std::max(quality.maxTemporalRootEnclosureWidth,root.upper-root.lower);
    std::vector<double> timeCuts{0,1};
    for (const auto& root:timeRoots.roots) if (root.point>0 && root.point<1) timeCuts.push_back(root.point);
    std::sort(timeCuts.begin(),timeCuts.end());
    const auto faces=domain.faces(mesh);
    std::vector<SpaceTimeReport> faceReports(faces.size());
    const auto integrateFace=[&](std::size_t index) {
        const auto face=faces[static_cast<int>(index)];
        auto& quality=faceReports[index];
        const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
        auto first=mesh.verts(edge.edgeToVert[0]).coord(),last=mesh.verts(edge.edgeToVert[1]).coord();
        if (face.horizontal?first.x>last.x:first.y>last.y) std::swap(first,last);
        const auto same=[](Vector2D a,Vector2D b){return a.x==b.x && a.y==b.y;};
        for (const auto* view:{&previousFlux,&currentFlux}) {
            const auto chart=view->faceChart(face.nativeRow);
            require(same(chart.begin,first) && same(chart.midpoint,face.center) && same(chart.end,last),"native face chart mismatch");
        }
        const bool inactive=face.negativeCell<0 && face.positiveCell<0;
        const bool internalWall=((face.negativeCell<0)!=(face.positiveCell<0)) &&
            !(face.horizontal?(first.y==mesh.minCoord.y || first.y==mesh.maxCoord.y):
                (first.x==mesh.minCoord.x || first.x==mesh.maxCoord.x));
        if (inactive || internalWall) {
            // Zero net flux is insufficient: exclude hidden spatial or
            // temporal counterflow through solid faces in either time layer.
            for (const auto* view:{&previousFlux,&currentFlux}) for (int half:{0,1})
                require(view->interpolatedRoots(*view,face.nativeRow,half,0).identicallyZero,
                    "inactive or internal solid face has nonzero normal-flux density");
            result[face.nativeRow]=InterfaceConservativePhase::FaceFlux::fromDirectional({});
            return;
        }
        const int owner=face.negativeCell>=0?face.negativeCell:face.positiveCell;
        // Every checked RK displacement is <= h/2. A full h halo, rounded
        // outward, includes all actual donor leaves of this native subface,
        // including a coarse owner and its several fine neighbors.
        const double infinity=std::numeric_limits<double>::infinity();
        auto nearby=coverage.locator().intersectingBox(
            {std::nextafter(first.x-h,-infinity),std::nextafter(first.y-h,-infinity)},
            {std::nextafter(last.x+h,infinity),std::nextafter(last.y+h,infinity)});
        nearby.erase(std::remove_if(nearby.begin(),nearby.end(),[&](int c){return domain.cellRowByIndex()[c]<0;}),nearby.end());
        Rates faceVolume{};
        for (int half:{0,1}) {
            const auto a=half==0?first:face.center,b=half==0?face.center:last;
            std::map<std::pair<double,double>,Vector2D> departures;
            const LevelSetBoundaryTrace::Velocity oldVelocity=[&](const Vector2D& p) {
                const auto u=velocities.previous(p);
                require(std::isfinite(u.x) && std::isfinite(u.y) && std::hypot(u.x,u.y)*time.dt<=.5*h,"full-step trace CFL exceeded");return u;
            };
            const auto tracePoint=[&](double fraction,double t) {
                const LevelSetBoundaryTrace::Velocity atTime=[&](const Vector2D& p) {
                    const auto u=velocities.current(p),v=oldVelocity(p);
                    require(std::isfinite(u.x) && std::isfinite(u.y) && std::hypot(u.x,u.y)*time.dt<=.5*h,"full-step trace CFL exceeded");
                    return v*(1-fraction)+u*fraction;
                };

                const auto found=departures.find({fraction,t});if (found!=departures.end()) return found->second;
                const auto p=t==0?a:t==1?b:a+(b-a)*t;
                const auto trace=LevelSetBoundaryTrace::trace(bounds,p,time.dt*fraction,atTime,oldVelocity);
                if (trace.boundaryHit) {
                    int phase=0;
                    for (const auto& port:portFaces) {
                        if (!trace.side || *trace.side!=port.side) continue;
                        const auto p=trace.point;
                        const bool contains=port.horizontal?p.y==port.lower.y && p.x>=port.lower.x && p.x<=port.upper.x:
                            p.x==port.lower.x && p.y>=port.lower.y && p.y<=port.upper.y;
                        if (!contains) continue;
                        require(phase==0 || phase==port.phase,"inflow event has conflicting incident port histories");
                        phase=port.phase;
                    }
                    require(phase!=0,"physical inflow needs an explicit single-phase port history");
                    ++quality.inflowTraceEvents;
                    // The event remains at the physical boundary. Its
                    // declared constant phase equals that entire old donor
                    // cell's phase, so the existing chord partition may
                    // use this boundary value. No continuation point or
                    // velocity outside the domain is constructed.
                }
                ++quality.tracedPoints;departures.emplace(std::make_pair(fraction,t),trace.point);return trace.point;
            };
            // Seed time quadrature with crossings of the OLD donor geometry.
            // Interior Gauss nodes can otherwise miss a short phase interval
            // at either end while coarse/fine estimates both report zero.
            // Use the same characteristic callback as the actual phase slices;
            // no exact benchmark geometry or inventory alteration is involved.
            std::vector<double> phaseTimeCuts=timeCuts;
            bool mixed=false;
            for(int c:nearby) mixed=mixed || !donors[c].uniform || donors[c].fraction!=donors[owner].fraction;
            const bool movingHalf=!previousFlux.interpolatedRoots(currentFlux,face.nativeRow,half,0).identicallyZero
                || !previousFlux.interpolatedRoots(currentFlux,face.nativeRow,half,1).identicallyZero;
            if(mixed && movingHalf) {
                namespace Events=Quantum::Geometry::QuadraticChordEvents;
                const auto point=[](Vector2D p){return Cut::Point{p.x,p.y};};
                std::vector<const std::vector<Cut::Point>*> vertices(nearby.size());
                for(std::size_t k=0;k<nearby.size();++k) {
                    vertices[k]=&arrangementVertices(nearby[k],quality);
                }
                // A short liquid pocket can be born/die INSIDE the half face
                // without crossing any of its three sampled trajectories.
                // The old donor arrangement changes when a mapped chord meets
                // a vertex (cell corner/interface edge root), or is tangent to
                // a quadratic donor. Seed both types on a bilinear trace map.
                const auto arrangement=[&](const auto& self,double lo,double hi,double left,double right,int td,int sd) -> void {
                    const double middle=lo+(hi-lo)*.5,centre=left+(right-left)*.5;
                    const auto a0=tracePoint(lo,left),b0=tracePoint(lo,right),a1=tracePoint(hi,left),b1=tracePoint(hi,right);
                    const auto am=tracePoint(middle,left),bm=tracePoint(middle,right);
                    const auto m0=tracePoint(lo,centre),m1=tracePoint(hi,centre),mm=tracePoint(middle,centre);
                    const auto distance=[](Vector2D p,Vector2D q){return std::hypot(p.x-q.x,p.y-q.y);};
                    const double spaceError=std::max(distance(m0,(a0+b0)*.5),distance(m1,(a1+b1)*.5))/h;
                    const double timeError=std::max(distance(am,(a0+a1)*.5),distance(bm,(b0+b1)*.5))/h;
                    const double centreError=distance(mm,(a0+b0+a1+b1)*.25)/h;
                    if(std::max({spaceError,timeError,centreError})>options.mapToleranceCells) {
                        if(timeError>=spaceError) {
                            require(td<options.maxTimeDepth,"temporal arrangement map exhausted its bound");
                            self(self,lo,middle,left,right,td+1,sd);self(self,middle,hi,left,right,td+1,sd);
                        } else {
                            require(sd<options.maxSpaceDepth,"spatial arrangement map exhausted its bound");
                            self(self,lo,hi,left,centre,td,sd+1);self(self,lo,hi,centre,right,td,sd+1);
                        }
                        return;
                    }
                    phaseTimeCuts.insert(phaseTimeCuts.end(),{lo,middle,hi});
                    const Events::Sweep sweep{point(a0),point(b0),point(a1),point(b1)};
                    const double pad=options.mapToleranceCells*h;
                    const Vector2D lower{std::min({a0.x,b0.x,a1.x,b1.x})-pad,std::min({a0.y,b0.y,a1.y,b1.y})-pad};
                    const Vector2D upper{std::max({a0.x,b0.x,a1.x,b1.x})+pad,std::max({a0.y,b0.y,a1.y,b1.y})+pad};
                    const auto addEvents=[&](const Quantum::Geometry::UnitPolynomialRoots& events,const auto& value) {
                        std::vector<double> seeds{0,1};for(const auto r:events.roots) if(r.point>0 && r.point<1) seeds.push_back(r.point);
                        std::sort(seeds.begin(),seeds.end());seeds.erase(std::unique(seeds.begin(),seeds.end()),seeds.end());
                        for(std::size_t k=1;k+1<seeds.size();++k) {
                            double l=lo+(hi-lo)*(seeds[k-1]+seeds[k])*.5,r=lo+(hi-lo)*(seeds[k]+seeds[k+1])*.5;
                            auto f=value(l),g=value(r);double root=lo+(hi-lo)*seeds[k];
                            if(f==0) root=l;else if(g==0) root=r;
                            else if((f<0)!=(g<0)) {
                                for(int j=0;j<52;++j) {const double m=l+(r-l)*.5;if(m==l || m==r) break;
                                    const auto v=value(m);if(v==0){l=r=m;break;}
                                    if((v<0)==(f<0)){l=m;f=v;}else r=m;}
                                root=l+(r-l)*.5;
                            }
                            if(root>lo && root<hi) phaseTimeCuts.push_back(root);
                        }
                    };
                    for(std::size_t k=0;k<nearby.size();++k) {
                        const auto& d=donors[nearby[k]];const auto& box=d.box;
                        if(box.upper.x<lower.x || box.lower.x>upper.x || box.upper.y<lower.y || box.lower.y>upper.y) continue;
                        for(auto vertex:*vertices[k]) {
                            const auto roots=Events::vertex(sweep,vertex);quality.temporalVertexSeeds+=roots.roots.size();
                            addEvents(roots,[&](double t){return Events::vertexValue(point(tracePoint(t,left)),point(tracePoint(t,right)),vertex);});
                        }
                        if(d.curved) {
                            const auto roots=Events::tangencies(box,*d.curved,sweep);quality.temporalTangencySeeds+=roots.roots.size();
                            addEvents(roots,[&](double t){return Events::tangencyValue(box,*d.curved,point(tracePoint(t,left)),point(tracePoint(t,right)));});
                        }
                    }
                };
                for(std::size_t t=1;t<timeCuts.size();++t) arrangement(arrangement,timeCuts[t-1],timeCuts[t],0,1,0,0);
            }
            if(mixed && movingHalf) for(double tangent:{0.,.5,1.}) {
                const auto events=[&](double lo,double hi,Vector2D x,Vector2D y) {
                    for(int c:nearby) {
                        const auto& d=donors[c];const auto& box=d.box;
                        if(std::max(x.x,y.x)<box.lower.x || std::min(x.x,y.x)>box.upper.x
                            || std::max(x.y,y.y)<box.lower.y || std::min(x.y,y.y)>box.upper.y) continue;
                        const auto addRoots=[&](std::vector<double> roots,const auto& value) {
                            roots.push_back(0);roots.push_back(1);std::sort(roots.begin(),roots.end());
                            roots.erase(std::unique(roots.begin(),roots.end()),roots.end());
                            for(std::size_t k=1;k+1<roots.size();++k) {
                                double left=lo+(hi-lo)*(roots[k-1]+roots[k])*.5;
                                double right=lo+(hi-lo)*(roots[k]+roots[k+1])*.5;
                                auto f=value(tracePoint(left,tangent)),g=value(tracePoint(right,tangent));
                                double root=lo+(hi-lo)*roots[k];
                                // Polish simple crossings on the actual trace.
                                // Unbracketed chord tangencies remain seeds,
                                // not a claim of certified temporal isolation.
                                if(f==0) root=left;
                                else if(g==0) root=right;
                                else if((f<0)!=(g<0)) {
                                    for(int j=0;j<52;++j) {
                                        const double middle=left+(right-left)*.5;
                                        if(middle==left || middle==right) break;
                                        const auto m=value(tracePoint(middle,tangent));
                                        if(m==0) {left=right=middle;break;}
                                        if((m<0)==(f<0)) {left=middle;f=m;}else right=middle;
                                    }
                                    root=left+(right-left)*.5;
                                }
                                if(root>lo && root<hi) phaseTimeCuts.push_back(root);
                            }
                        };
                        const auto linear=[&](const auto& value) {
                            const auto u=value(x),v=value(y);std::vector<double> roots;
                            if(u!=v) {const auto r=u/(u-v);if(r>0 && r<1) roots.push_back(static_cast<double>(r));}
                            addRoots(std::move(roots),value);
                        };
                        // Cell boundaries also matter: neighboring PLICs can
                        // have different one-sided phase traces there.
                        for(double bound:{box.lower.x,box.upper.x}) linear([=](Vector2D p){return static_cast<long double>(p.x)-bound;});
                        for(double bound:{box.lower.y,box.upper.y}) linear([=](Vector2D p){return static_cast<long double>(p.y)-bound;});
                        if(d.uniform) continue;
                        const long double cx=(static_cast<long double>(box.lower.x)+box.upper.x)/2;
                        const long double cy=(static_cast<long double>(box.lower.y)+box.upper.y)/2;
                        if(d.curved) {
                            std::vector<double> roots;
                            for(const auto r:Curve::intervals(box,*d.curved,{x.x,x.y},{y.x,y.y}))
                                for(auto t:r) if(t>0 && t<1) roots.push_back(static_cast<double>(t));
                            const auto curve=*d.curved;
                            addRoots(std::move(roots),[=](Vector2D p){return Curve::implicitValue(box,curve,{p.x,p.y});});
                        } else linear([&](Vector2D p){return d.plane.normal.x*(p.x-cx)+d.plane.normal.y*(p.y-cy)-d.plane.offset;});
                    }
                };
                const auto partition=[&](const auto& self,double lo,double hi,int depth) -> void {
                    const double middle=lo+(hi-lo)*.5;
                    const auto x=tracePoint(lo,tangent),y=tracePoint(hi,tangent),m=tracePoint(middle,tangent);
                    const double deviation=std::hypot(m.x-(x.x+y.x)*.5,m.y-(x.y+y.y)*.5)/h;
                    if(deviation>options.mapToleranceCells) {
                        require(depth<options.maxTimeDepth,"temporal event trajectory exhausted its map bound");
                        self(self,lo,middle,depth+1);self(self,middle,hi,depth+1);
                    } else {
                        // Retain trace-leaf boundaries too: an event can lie
                        // exactly on a shared endpoint of the chord pieces.
                        phaseTimeCuts.insert(phaseTimeCuts.end(),{lo,middle,hi});
                        events(lo,middle,x,m);events(middle,hi,m,y);
                    }
                };
                for(std::size_t t=1;t<timeCuts.size();++t) partition(partition,timeCuts[t-1],timeCuts[t],0);
            }
            std::sort(phaseTimeCuts.begin(),phaseTimeCuts.end());
            phaseTimeCuts.erase(std::unique(phaseTimeCuts.begin(),phaseTimeCuts.end()),phaseTimeCuts.end());
            quality.temporalGeometryCuts+=phaseTimeCuts.size()-timeCuts.size();
            std::map<double,Rates> slices;
            const auto slice=[&](double fraction) -> Rates {
                const auto found=slices.find(fraction);if (found!=slices.end()) return found->second;
                ++quality.timeSlices;
                const auto roots=previousFlux.interpolatedRoots(currentFlux,face.nativeRow,half,fraction);quality.rootEnclosures+=roots.roots.size();
                if (roots.identicallyZero) {slices.emplace(fraction,Rates{});return Rates{};}
                std::vector<double> cuts{0,1};
                for (const auto& root:roots.roots) {
                    quality.maxRootEnclosureWidth=std::max(quality.maxRootEnclosureWidth,root.upper-root.lower);
                    if (root.point>0 && root.point<1) cuts.push_back(root.point);
                }
                std::sort(cuts.begin(),cuts.end());
                // The spatial coarse/fine and phase-partition checks query
                // the same anchored primitive repeatedly. Reuse only exact
                // parameter tuples within this frozen half-face/time slice;
                // keep every partition, quadrature estimate and comparison.
                // Release the cache before moving to another time or face.
                std::map<std::array<double,4>,long double> primitives;
                const auto integral=[&](double lo,double hi,double partLo=0,double partHi=1) {
                    ++quality.primitiveRequests;
                    const std::array<double,4> key{lo,hi,partLo,partHi};
                    const auto found=primitives.find(key);
                    if(found!=primitives.end()) return found->second;
                    const auto value=previousFlux.interpolatedIntegral(currentFlux,face.nativeRow,half,fraction,lo,hi,partLo,partHi);
                    ++quality.primitiveEvaluations;primitives.emplace(key,value);return value;
                };
                const auto departure=[&](double t){return tracePoint(fraction,t);};
                const auto chord=[&](double lo,double hi,Vector2D x,Vector2D y) {
                    require(coverage.covers(x,y),"departure chord crosses outside selected fluid cells");
                    // The CFL-enclosing donor halo has one exactly constant
                    // pure phase. Its characteristic geometry cannot change
                    // the phase integral: integrate the retained native flux
                    // primitive directly. Keep all trace, domain coverage,
                    // sign-root and map-refinement checks. Tiny nonzero tails
                    // or differing donor fractions still use the full split.
                    if(!mixed && (donors[owner].fraction==0 || donors[owner].fraction==1)) {
                        const auto total=integral(lo,hi);
                        require(std::isfinite(total) && total!=0,"constant phase chord has no finite native volume primitive");
                        ++quality.constantPhaseChords;
                        return donors[owner].fraction==0?0.L:total;
                    }
                    const long double dx=static_cast<long double>(y.x)-x.x,dy=static_cast<long double>(y.y)-x.y;
                    struct Span {int cell;long double lo,hi;};
                    std::vector<Span> spans;std::vector<long double> partition{0,1};
                    for (int c:nearby) {
                        const auto& box=donors[c].box;
                        long double left=0,right=1;
                        if (!clip(0,dx,static_cast<long double>(box.upper.x)-x.x,left,right)
                            || !clip(0,-dx,static_cast<long double>(x.x)-box.lower.x,left,right)
                            || !clip(0,dy,static_cast<long double>(box.upper.y)-x.y,left,right)
                            || !clip(0,-dy,static_cast<long double>(x.y)-box.lower.y,left,right)) continue;
                        spans.push_back({c,left,right});partition.push_back(left);partition.push_back(right);
                    }
                    std::sort(partition.begin(),partition.end());
                    partition.erase(std::unique(partition.begin(),partition.end()),partition.end());
                    std::vector<Range> liquid;
                    for(std::size_t k=1;k<partition.size();++k) {
                        long double left=partition[k-1],right=partition[k];int selected=-1;
                        // Split at every real leaf boundary before assigning
                        // ownership. A coarse upper edge may face two fine
                        // lower edges (or only a partly selected neighbor).
                        // Positive-side ownership is per interval, never a
                        // whole-coarse-cell skip or a rounded midpoint test.
                        for(const auto& span:spans) if(span.lo<=left && span.hi>=right) {
                            if(selected<0) {selected=span.cell;continue;}
                            const auto& old=donors[selected].box;const auto& next=donors[span.cell].box;
                            const bool vertical=dx==0 && (x.x==old.lower.x || x.x==old.upper.x)
                                && (x.x==next.lower.x || x.x==next.upper.x);
                            const bool horizontal=dy==0 && (x.y==old.lower.y || x.y==old.upper.y)
                                && (x.y==next.lower.y || x.y==next.upper.y);
                            require(vertical || horizontal,"overlapping interior donor intervals");
                            if((vertical && next.lower.x>old.lower.x)
                                || ((!vertical || next.lower.x==old.lower.x) && horizontal && next.lower.y>old.lower.y))
                                selected=span.cell;
                        }
                        require(selected>=0,"departure chord not covered once by its donor cells");
                        const auto& d=donors[selected];const auto& box=d.box;
                        if (d.fraction==0) continue;
                        if (d.curved) {
                            for (auto range:Curve::intervals(box,*d.curved,{x.x,x.y},{y.x,y.y},left,right))
                                liquid.push_back({static_cast<double>(range[0]),static_cast<double>(range[1]),1});
                            continue;
                        }
                        const long double cx=static_cast<long double>(box.lower.x)+(static_cast<long double>(box.upper.x)-box.lower.x)/2;
                        const long double cy=static_cast<long double>(box.lower.y)+(static_cast<long double>(box.upper.y)-box.lower.y)/2;
                        if (!d.uniform && !clip(d.plane.normal.x*(static_cast<long double>(x.x)-cx)
                            +d.plane.normal.y*(static_cast<long double>(x.y)-cy),d.plane.normal.x*dx+d.plane.normal.y*dy,
                            d.plane.offset,left,right)) continue;
                        liquid.push_back({static_cast<double>(left),static_cast<double>(right),d.uniform?d.fraction:1});
                    }
                    std::sort(liquid.begin(),liquid.end(),[](auto x,auto y){return x.begin<y.begin;});
                    std::vector<Range> merged;
                    for (const auto range:liquid) {
                        if (!merged.empty() && range.begin==merged.back().end && range.weight==merged.back().weight)
                            merged.back().end=range.end;
                        else {require(merged.empty() || range.begin>=merged.back().end,"overlapping liquid donor intervals");merged.push_back(range);}
                    }
                    // Use one common partition for both phases. Independently
                    // summing several liquid intervals can exceed a separately
                    // evaluated whole primitive by roundoff near full liquid.
                    // Reconcile that metric difference, never clip a phase.
                    const auto total=integral(lo,hi);long double liquidRate=0,gasRate=0;double cursor=0;
                    const auto deposit=[&](double begin,double end,long double liquidFraction) {
                        if (begin==end) return;
                        const auto piece=integral(lo,hi,begin,end);
                        require(std::isfinite(piece) && (total>0?piece>=0:piece<=0),"phase partition crossed an unisolated flux zero");
                        liquidRate+=liquidFraction*piece;gasRate+=(1-liquidFraction)*piece;
                    };
                    for (const auto range:merged) {
                        deposit(cursor,range.begin,0);deposit(range.begin,range.end,range.weight);cursor=range.end;
                    }
                    deposit(cursor,1,0);const auto sum=liquidRate+gasRate;
                    if(!(sum!=0 && std::abs(sum-total)<=256*std::numeric_limits<double>::epsilon()*std::abs(total))) {
                        std::ostringstream message;message<<std::setprecision(17)
                            <<"[InterfacePlicSpaceTime] liquid/gas partition changed its native volume primitive"
                            <<"; face="<<face.nativeRow<<", half="<<half<<", time fraction="<<fraction
                            <<", interval=["<<lo<<','<<hi<<"], total="<<total<<", split="<<sum
                            <<", liquid="<<liquidRate<<", gas="<<gasRate<<", ranges="<<merged.size();
                        for(const auto range:merged) message<<" ["<<range.begin<<','<<range.end<<']';
                        throw std::invalid_argument(message.str());
                    }
                    const auto fraction=liquidRate/sum;
                    require(std::isfinite(fraction) && fraction>=0 && fraction<=1,"invalid liquid/gas partition ratio");
                    return total*fraction;
                };
                using Spatial=std::pair<long double,long double>; // phase, observed error
                struct SpaceInterval {double lo,hi,deviation;int depth;long double total,phase,error;};
                const auto consistent=[&](long double total,long double splitTotal,long double phase) {
                    if (!(splitTotal!=0 && std::abs(splitTotal-total)<=256*std::numeric_limits<double>::epsilon()*std::abs(total))) {
                        std::ostringstream message;message<<std::setprecision(17)
                            <<"[InterfacePlicSpaceTime] spatial partition changed its native volume primitive"
                            <<"; face="<<face.nativeRow<<", half="<<half<<", time fraction="<<fraction
                            <<", total="<<total<<", split="<<splitTotal<<", phase="<<phase;
                        throw std::invalid_argument(message.str());
                    }
                    const auto fraction=phase/splitTotal;
                    require(std::isfinite(fraction) && fraction>=0 && fraction<=1,"invalid geometric phase fraction in spatial partition");
                    return total*fraction; // Metric reconciliation only; no phase clipping.
                };
                const auto evaluateSpace=[&](double lo,double hi,int depth) {
                    const double middle=lo+(hi-lo)*.5;require(middle>lo && middle<hi,"unresolved spatial quadrature interval");
                    const auto x=departure(lo),y=departure(hi),m=departure(middle);
                    const auto total=integral(lo,hi),splitTotal=integral(lo,middle)+integral(middle,hi);
                    const auto coarse=chord(lo,hi,x,y),fine=consistent(total,splitTotal,chord(lo,middle,x,m)+chord(middle,hi,m,y));
                    const double deviation=std::hypot(m.x-(x.x+y.x)*.5,m.y-(x.y+y.y)*.5)/h;
                    return SpaceInterval{lo,hi,deviation,depth,total,fine,std::abs(fine-coarse)};
                };
                const auto spatial=[&](double lo,double hi) -> Spatial {
                    const auto total=integral(lo,hi),gross=std::abs(total);
                    const auto budget=(options.relativeFluxTolerance/8)*gross;
                    std::vector<SpaceInterval> intervals{evaluateSpace(lo,hi,0)};
                    for (;;) {
                        long double sum=0,phase=0,error=0,depthError=0;double deviation=0;int depth=0;
                        std::size_t refine=intervals.size();bool needsMap=false,depthMap=false;
                        for (const auto& interval:intervals) {
                            sum+=interval.total;phase+=interval.phase;error+=interval.error;
                            deviation=std::max(deviation,interval.deviation);
                            depth=std::max(depth,interval.depth);
                            needsMap=needsMap || interval.deviation>options.mapToleranceCells;
                            if(interval.depth>=options.maxSpaceDepth) {
                                depthError+=interval.error;
                                depthMap=depthMap || interval.deviation>options.mapToleranceCells;
                            }
                        }
                        const auto normalized=consistent(total,sum,phase);error*=std::abs(total/sum);
                        if (!needsMap && error<=budget) {
                            quality.spaceLeaves+=intervals.size();
                            quality.maxAcceptedSpaceLeaves=std::max(quality.maxAcceptedSpaceLeaves,intervals.size());
                            quality.maxAcceptedSpaceDepth=std::max(quality.maxAcceptedSpaceDepth,depth);
                            quality.maxAcceptedMapDeviationCells=std::max(quality.maxAcceptedMapDeviationCells,deviation);
                            quality.maxAcceptedSpatialRelativeError=std::max(quality.maxAcceptedSpatialRelativeError,static_cast<double>(error/gross));
                            return {normalized,error};
                        }
                        for (std::size_t i=0;i<intervals.size();++i) {
                            const auto& interval=intervals[i];
                            if (interval.depth>=options.maxSpaceDepth || (needsMap && interval.deviation<=options.mapToleranceCells)) continue;
                            if (refine==intervals.size() || (needsMap?interval.deviation>intervals[refine].deviation:interval.error>intervals[refine].error)) refine=i;
                        }
                        // Frozen leaves cannot lose error by subdividing a
                        // different leaf. The metric reconciliation above
                        // permits only 256 eps of relative change in sum;
                        // twice that allowance makes this rejection robust
                        // to rounding. Avoid filling the tree with irrelevant
                        // subdivisions after its error budget is unattainable.
                        const bool depthBudget=depthError>budget*(1+512*std::numeric_limits<double>::epsilon());
                        if (depthMap || depthBudget || refine==intervals.size() || intervals.size()>=options.maxSpaceLeaves) {
                            std::ostringstream message;message<<std::setprecision(17)<<"[InterfacePlicSpaceTime] spatial phase/map quadrature exhausted its bound"
                                <<"; face="<<face.nativeRow<<", half="<<half<<", time fraction="<<fraction
                                <<", error="<<error<<", budget="<<budget<<", deviation/h="<<deviation
                                <<", depth="<<depth<<'/'<<options.maxSpaceDepth
                                <<", leaves="<<intervals.size()<<'/'<<options.maxSpaceLeaves
                                <<", frozen_error="<<depthError<<", frozen_map="<<depthMap;
                            throw std::invalid_argument(message.str());
                        }
                        const auto old=intervals[refine];const double middle=old.lo+(old.hi-old.lo)*.5;
                        intervals[refine]=evaluateSpace(old.lo,middle,old.depth+1);
                        intervals.insert(intervals.begin()+refine+1,evaluateSpace(middle,old.hi,old.depth+1));
                    }
                };
                Rates value{};
                for (std::size_t i=1;i<cuts.size();++i) {
                    const double lo=cuts[i-1],hi=cuts[i];require(hi>lo,"collapsed density roots");
                    const auto total=integral(lo,hi);if (total==0) continue;
                    const auto phase=spatial(lo,hi);const bool positive=total>0;
                    const auto magnitude=std::abs(total),liquid=positive?phase.first:-phase.first;
                    if (!(std::isfinite(liquid) && liquid>=0 && liquid<=magnitude)) {
                        std::ostringstream message;message<<std::setprecision(17)<<"[InterfacePlicSpaceTime] unresolved signed liquid phase interval"
                            <<"; face="<<face.nativeRow<<", half="<<half<<", time fraction="<<fraction
                            <<", interval=["<<lo<<','<<hi<<"], total="<<magnitude<<", liquid="<<liquid;
                        throw std::invalid_argument(message.str());
                    }
                    value[positive?0:1]+=magnitude;value[positive?2:3]+=liquid;
                    value[4]+=phase.second;
                }
                slices.emplace(fraction,value);return value;
            };
            const auto quadrature=[&](double lo,double hi) {
                const double middle=(lo+hi)*.5,offset=(hi-lo)*.5/std::sqrt(3.);
                return scale(add(slice(middle-offset),slice(middle+offset)),(hi-lo)*.5L);
            };
            struct TimeInterval {double lo,hi;int depth;Rates left,right,fine;long double error;};
            const auto evaluateTime=[&](double lo,double hi,const Rates& coarse,int depth) {
                const double middle=lo+(hi-lo)*.5;
                if(!(middle>lo && middle<hi)) {
                    // Keep adjacent representable event times distinct. With
                    // no interior floating point, refinement cannot observe
                    // a phase transition. Charge its ENTIRE possible phase
                    // volume to the original half-face error budget instead.
                    // The retained flux interpolates linearly in time, so its
                    // spatial gross norm is convex: endpoint trapezoids bound
                    // every directional total and liquid integral, independent
                    // of the interface motion or any missed phase pocket.
                    const auto a=slice(lo),b=slice(hi);
                    const long double width=static_cast<long double>(hi)-lo;
                    const auto grossBound=width*(a[0]+a[1]+b[0]+b[1])/2;
                    const auto bound=std::max(grossBound,coarse[0]+coarse[1])
                        *(1+32*std::numeric_limits<double>::epsilon());
                    require(std::isfinite(bound) && bound>=0,"invalid indivisible-time volume bound");
                    ++quality.indivisibleTimeIntervals;quality.indivisibleTimeVolumeBound+=time.dt*bound;
                    return TimeInterval{lo,hi,options.maxTimeDepth,{},{},coarse,bound};
                }
                const auto left=quadrature(lo,middle),right=quadrature(middle,hi),fine=add(left,right);
                long double error=0;for (int i=0;i<4;++i) error=std::max(error,std::abs(fine[i]-coarse[i]));
                return TimeInterval{lo,hi,depth,left,right,fine,error};
            };
            std::vector<TimeInterval> intervals;
            for (std::size_t t=1;t<phaseTimeCuts.size();++t) {
                const double lo=phaseTimeCuts[t-1],hi=phaseTimeCuts[t];
                require(hi>lo,"collapsed temporal density roots");
                intervals.push_back(evaluateTime(lo,hi,quadrature(lo,hi),0));
            }
            // One integrated half-face budget. A geometric crossing may put
            // most error into a short time interval; do not force each leaf
            // to spend exactly its length-proportional share. Refine the
            // largest observed contribution, retaining all original bounds.
            for (;;) {
                Rates sum{};long double error=0;std::size_t refine=intervals.size();
                for (std::size_t i=0;i<intervals.size();++i) {
                    const auto& interval=intervals[i];sum=add(sum,interval.fine);error+=interval.error;
                    if (interval.depth<options.maxTimeDepth && (refine==intervals.size() || interval.error>intervals[refine].error)) refine=i;
                }
                const auto budget=options.relativeFluxTolerance*(sum[0]+sum[1]);
                if (error+sum[4]<=budget) {
                    quality.timeLeaves+=intervals.size();quality.estimatedTemporalError+=time.dt*error;
                    faceVolume=add(faceVolume,scale(sum,time.dt));break;
                }
                if (refine==intervals.size()) {
                    std::ostringstream message;message<<std::setprecision(17)<<"[InterfacePlicSpaceTime] temporal phase quadrature exhausted its bound"
                        <<"; face="<<face.nativeRow<<", half="<<half<<", leaves="<<intervals.size()
                        <<", temporal error="<<error<<", spatial error="<<sum[4]<<", budget="<<budget;
                    throw std::invalid_argument(message.str());
                }
                const auto old=intervals[refine];const double middle=old.lo+(old.hi-old.lo)*.5;
                intervals[refine]=evaluateTime(old.lo,middle,old.left,old.depth+1);
                intervals.insert(intervals.begin()+refine+1,evaluateTime(middle,old.hi,old.right,old.depth+1));
            }
        }
        const auto net=faceVolume[0]-faceVolume[1];long double exact=0;
        quality.estimatedSpatialError+=faceVolume[4];
        for (int half:{0,1}) exact+=.5L*time.dt*(previousFlux.integratedNormalFlux(face.nativeRow,half,0,1)
            +currentFlux.integratedNormalFlux(face.nativeRow,half,0,1));
        const auto difference=std::abs(net-exact);quality.maxNetFluxDifference=std::max(quality.maxNetFluxDifference,difference);
        require(difference<=256*std::numeric_limits<double>::epsilon()*(faceVolume[0]+faceVolume[1]),"space-time quadrature changed native net volume");
        // Geometry quadrature determines the directional phase fractions,
        // not a replacement total-volume flux. Restore the native two-layer
        // integral after its discrepancy has passed the roundoff check. Keep
        // counterflow and rescale only the larger direction, preserving its
        // phase fraction (including exact pure phases). Otherwise independent
        // temporal subdivisions can destroy GCL in a full liquid cell.
        const int direction=exact>=0?0:1;
        const auto magnitude=faceVolume[1-direction]+std::abs(exact),old=faceVolume[direction];
        require(old>0 || magnitude==0,"native net flux has no quadrature direction");
        if (old>0) faceVolume[direction+2]=magnitude*(faceVolume[direction+2]/old);
        faceVolume[direction]=magnitude;
        const int incoming=face.negativeCell<0?0:face.positiveCell<0?1:-1;
        long double boundaryFraction=std::numeric_limits<long double>::quiet_NaN();
        if (incoming>=0 && faceVolume[incoming]>0) {
            require(portPhase[face.nativeRow]!=0,"physical inflow needs phase history");
            boundaryFraction=portPhase[face.nativeRow]<0?1:0;
            require(faceVolume[incoming+2]==boundaryFraction*faceVolume[incoming],
                "incoming phase integral disagrees with the declared pure port");
            quality.boundaryInflowVolume+=faceVolume[incoming];quality.boundaryInflowLiquid+=faceVolume[incoming+2];
        }
        result[face.nativeRow]=InterfaceConservativePhase::FaceFlux::fromDirectional(
            {faceVolume[0],faceVolume[1],faceVolume[2],faceVolume[3]},boundaryFraction);
    };
    // Keep the first active query on the caller for the strict geometry/FP
    // contract. Each worker owns one native face, its quadrature and report.
    std::size_t first=0;
    while(first<faces.size() && faces[static_cast<int>(first)].negativeCell<0
        && faces[static_cast<int>(first)].positiveCell<0) integrateFace(first++);
    if(first<faces.size()) integrateFace(first++);
    std::vector<std::exception_ptr> failures(faces.size());
#pragma omp parallel for schedule(dynamic, 8) if(options.parallelFaces && faces.size()>128)
    for(std::ptrdiff_t index=static_cast<std::ptrdiff_t>(first);index<static_cast<std::ptrdiff_t>(faces.size());++index) {
        try {integrateFace(static_cast<std::size_t>(index));}
        catch(...) {failures[index]=std::current_exception();}
    }
    // Replay the original native-face failure order; publish nothing on a
    // failed call. Only diagnostics are summed here, in a fixed serial order.
    for(const auto& failure:failures) if(failure) std::rethrow_exception(failure);
    for(const auto& q:faceReports) {
        quality.tracedPoints+=q.tracedPoints;quality.timeSlices+=q.timeSlices;
        quality.spaceLeaves+=q.spaceLeaves;quality.timeLeaves+=q.timeLeaves;quality.rootEnclosures+=q.rootEnclosures;
        quality.maxAcceptedSpaceLeaves=std::max(quality.maxAcceptedSpaceLeaves,q.maxAcceptedSpaceLeaves);
        quality.maxAcceptedSpaceDepth=std::max(quality.maxAcceptedSpaceDepth,q.maxAcceptedSpaceDepth);
        quality.constantPhaseChords+=q.constantPhaseChords;
        quality.primitiveRequests+=q.primitiveRequests;quality.primitiveEvaluations+=q.primitiveEvaluations;
        quality.arrangementVertexRequests+=q.arrangementVertexRequests;
        quality.arrangementVertexBuilds+=q.arrangementVertexBuilds;
        quality.arrangementEdgeRootQueries+=q.arrangementEdgeRootQueries;
        quality.temporalGeometryCuts+=q.temporalGeometryCuts;
        quality.indivisibleTimeIntervals+=q.indivisibleTimeIntervals;
        quality.indivisibleTimeVolumeBound+=q.indivisibleTimeVolumeBound;
        quality.temporalVertexSeeds+=q.temporalVertexSeeds;quality.temporalTangencySeeds+=q.temporalTangencySeeds;
        quality.estimatedTemporalError+=q.estimatedTemporalError;quality.estimatedSpatialError+=q.estimatedSpatialError;
        quality.maxNetFluxDifference=std::max(quality.maxNetFluxDifference,q.maxNetFluxDifference);
        quality.maxAcceptedMapDeviationCells=std::max(quality.maxAcceptedMapDeviationCells,q.maxAcceptedMapDeviationCells);
        quality.maxAcceptedSpatialRelativeError=std::max(quality.maxAcceptedSpatialRelativeError,q.maxAcceptedSpatialRelativeError);
        quality.maxRootEnclosureWidth=std::max(quality.maxRootEnclosureWidth,q.maxRootEnclosureWidth);
        quality.inflowTraceEvents+=q.inflowTraceEvents;quality.boundaryInflowVolume+=q.boundaryInflowVolume;
        quality.boundaryInflowLiquid+=q.boundaryInflowLiquid;
    }
    if (report) *report=quality;
    return result;
}
}
