#include "Modules/Interface/Compatibility/InterfaceMeshTransfer.h"
#include "Modules/Interface/Compatibility/InterfaceMotion.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Geometry/ActiveHermiteTransfer.h"
#include "Core/Geometry/InterfaceGeometryResolution.h"
#include "Core/Geometry/ClosedCurveVolumeFit.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <array>
#include <stdexcept>

namespace {
// Follow the actual native fluid boundary with fluid on the left. Intermediate
// vertices only: the open interface's endpoints already belong to its trace.
std::vector<Quantum::Geometry::ClosedCurveVolumeFit::Point> boundaryClosure(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,Vector2D from,Vector2D to)
{
    using Point=Quantum::Geometry::ClosedCurveVolumeFit::Point;
    using Key=std::array<double,2>;
    struct Segment {Vector2D a,b;};std::vector<Segment> edges;std::map<Key,int> outgoing;
    const auto key=[](Vector2D p){return Key{p.x,p.y};};
    const auto on=[](const Segment& e,Vector2D p){return (e.a.x==e.b.x?p.x==e.a.x:p.y==e.a.y)
        && p.x>=std::min(e.a.x,e.b.x) && p.x<=std::max(e.a.x,e.b.x)
        && p.y>=std::min(e.a.y,e.b.y) && p.y<=std::max(e.a.y,e.b.y);};
    for(const auto f:domain.faces(mesh))if((f.negativeCell>=0)!=(f.positiveCell>=0)){
        const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
        auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
        if(f.horizontal?a.x>b.x:a.y>b.y)std::swap(a,b);
        if(f.horizontal?f.negativeCell>=0:f.positiveCell>=0)std::swap(a,b);
        if(!outgoing.emplace(key(a),static_cast<int>(edges.size())).second)
            throw std::invalid_argument("contact closure has a nonmanifold fluid boundary");
        edges.push_back({a,b});
    }
    int start=-1;
    for(std::size_t i=0;i<edges.size();++i)if(on(edges[i],from) && key(from)!=key(edges[i].b)){
        if(start>=0)throw std::invalid_argument("contact closure has an ambiguous starting face");start=static_cast<int>(i);}
    if(start<0 || key(from)==key(to))throw std::invalid_argument("contact closure endpoint is not on the native boundary");
    std::vector<Point> result;int current=start;Vector2D point=from;
    for(std::size_t i=0;i<=edges.size();++i){const auto& edge=edges[current];
        if(on({point,edge.b},to))return result;
        result.push_back({edge.b.x,edge.b.y});const auto next=outgoing.find(key(edge.b));
        if(next==outgoing.end())throw std::invalid_argument("contact closure encountered an open domain boundary");
        current=next->second;point=edge.b;
        if(current==start){if(on({edge.b,from},to))return result;break;}
    }
    throw std::invalid_argument("contact endpoints belong to different native boundary cycles");
}
std::vector<InterfaceMeshTransfer::ConservativeReport::SharedCurve> sharedReferences(
    const QuadtreeMesh& oldMesh,const InterfaceState& source,const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,const std::vector<std::vector<int>>& sourceCells,
    std::vector<long double>& reference,bool synchronizeGeometry)
{
    namespace F=Quantum::Geometry::ClosedCurveVolumeFit;
    namespace P=Quantum::Geometry::PlanarCellCut;
    using Report=InterfaceMeshTransfer::ConservativeReport::SharedCurve;
    const auto box=[](const QuadtreeMesh& m,int c){const auto& cell=m.cells(m.cellIndexToLabel[c]);
        return P::Rectangle{m.verts(cell.cellToVert[0]).coord(),m.verts(cell.cellToVert[2]).coord()};};
    std::vector<Report> reports;
    const auto& inventory=source.conservativePhase();const bool rz=domain.axisymmetric();
    std::vector<int> incidence(oldMesh.cellIndexNum(),0);
    for(const auto& surface:source.topology().surfaces())for(const auto& graph:surface.geometry.nativeChain().graphs())
        ++incidence[graph.nativeCellIndex];
    for(const auto& surface:source.topology().surfaces()) {
        Report report;report.surfaceId=surface.id;const auto& chain=surface.geometry.nativeChain();
        const auto attempt=[&]() -> std::string {
            using Kind=LevelSetInterfacePanels::EndpointKind;
            const bool axisToAxis=rz && !chain.closed() && surface.geometry.beginEndpointKind()==Kind::Axis
                && surface.geometry.endEndpointKind()==Kind::Axis;
            const bool contact=!chain.closed() && (surface.geometry.beginEndpointKind()==Kind::DeclaredBoundary
                || surface.geometry.endEndpointKind()==Kind::DeclaredBoundary);
            const bool axisToContact=contact && rz && (surface.geometry.beginEndpointKind()==Kind::Axis || surface.geometry.endEndpointKind()==Kind::Axis);
            if(!chain.closed() && !axisToAxis && !contact)return "open endpoint has no constrained representation";
            if(chain.graphs().size()<8)return "insufficient closed branch support";
            std::vector<F::Point> trace;std::set<int> parents;double h=std::numeric_limits<double>::infinity();
            for(std::size_t g=0;g<chain.graphs().size();++g) {
                const int c=chain.graphs()[g].nativeCellIndex;
                if(incidence[c]!=1)return "multiple branch incidences in one parent";
                if(!parents.insert(c).second)return "multiple arcs in one parent";
                const auto p=chain.frameAt(static_cast<int>(g),.5).point;trace.push_back({p.x,p.y});
                const auto size=oldMesh.cells(oldMesh.cellIndexToLabel[c]).size();h=std::min({h,size.x,size.y});
                const auto amount=inventory.liquidVolumes()[c],cap=inventory.capacities()[c];
                if(std::min(amount,cap-amount)<=P::relativeVolumeResolution*cap)return "unresolved parent inventory";
            }
            std::vector<F::Point> closure;
            if(axisToAxis || contact) {
                const auto first=chain.frameAt(0,0).point,last=chain.frameAt(static_cast<int>(chain.graphs().size())-1,1).point;
                trace.insert(trace.begin(),{first.x,first.y});trace.push_back({last.x,last.y});
                if(contact)closure=boundaryClosure(mesh,domain,last,first);
            }
            std::vector<int> targets;bool refined=false;
            for(int c:domain.cellIndexByRow()) if(sourceCells[c].size()==1 && parents.count(sourceCells[c][0])) {
                targets.push_back(c);refined|=mesh.cells(mesh.cellIndexToLabel[c]).level()>
                    oldMesh.cells(oldMesh.cellIndexToLabel[sourceCells[c][0]]).level();
            }
            if(!refined)return "no mixed parent refinement";
            auto polygon=trace;polygon.insert(polygon.end(),closure.begin(),closure.end());
            long double orientation=0;const auto origin=polygon.front();
            for(std::size_t i=0;i<polygon.size();++i){const auto a=polygon[i],b=polygon[(i+1)%polygon.size()];
                orientation+=(a.x-origin.x)*(b.y-origin.y)-(a.y-origin.y)*(b.x-origin.x);}
            auto tangent=F::Point{trace[1].x-trace.back().x,trace[1].y-trace.back().y};
            if(axisToAxis || contact){const auto a=chain.frameAt(0,.25).point,b=chain.frameAt(0,.75).point;tangent={b.x-a.x,b.y-a.y};}
            const auto normal=chain.frameAt(0,.5).normal;
            const bool inside=orientation*(normal.x*tangent.y-normal.y*tangent.x)>0;
            std::vector<F::Constraint> constraints;
            for(int c:parents)constraints.push_back({box(oldMesh,c),inside?inventory.liquidVolumes()[c]:
                inventory.capacities()[c]-inventory.liquidVolumes()[c]});
            const auto fitted=axisToAxis?F::fitAxisToAxis(trace,constraints,h):axisToContact?F::fitAxisToContact(trace,closure,constraints,h):
                contact?F::fitOpen(trace,closure,constraints,h,rz):F::fit(trace,constraints,h,rz);
            report.controls=fitted.controls;report.rank=fitted.rank;report.iterations=fitted.iterations;
            report.initialSupportRepairs=fitted.initialSupportRepairs;report.maxInitialDisplacement=fitted.maxInitialDisplacement;
            report.parentFractionResidual=fitted.maxFractionResidual;report.maxNormalDisplacement=fitted.maxNormalDisplacement;
            if(!fitted.converged)return "curve fit or polygon admissibility failed";
            std::vector<F::Point> coarse;for(std::size_t i=0;i<fitted.polygon.size();++i)if(i%2==0 ||
                std::find(fitted.fixedPolygonVertices.begin(),fitted.fixedPolygonVertices.end(),i)!=fitted.fixedPolygonVertices.end())coarse.push_back(fitted.polygon[i]);
            std::vector<F::Rectangle> rectangles;for(int c:targets)rectangles.push_back(box(mesh,c));
            const auto fineMeasures=F::measureMany(fitted.polygon,rectangles,rz),coarseMeasures=F::measureMany(coarse,rectangles,rz);
            std::vector<long double> values;values.reserve(targets.size());
            for(std::size_t i=0;i<targets.size();++i) {
                const int c=targets[i];const long double cap=domain.cellVolumes()[domain.cellRowByIndex()[c]];
                auto value=fineMeasures[i];const auto half=coarseMeasures[i];
                report.maxIntegrationDifference=std::max(report.maxIntegrationDifference,std::abs(value-half)/cap);
                if(!inside)value=cap-value;
                const auto tolerance=P::relativeVolumeResolution*cap;
                if(value< -tolerance || value>cap+tolerance)return "curve measure outside capacity";
                value=std::clamp(value,0.L,cap); // Only the checked arithmetic margin.
                if(reference[c]==0 || reference[c]==cap) {
                    if(std::abs(value-reference[c])>tolerance) {
                        if(!synchronizeGeometry)return "curve conflicts with geometric pure child";
                        ++report.changedPureChildren;
                    } else value=reference[c];
                }
                values.push_back(value);
            }
            // A refinement observation, not a certified quadrature bound.
            // Keep it two orders below the existing 1e-3 geometry/C budget.
            if(report.maxIntegrationDifference>1e-5)return "polygon integration needs refinement";
            for(std::size_t i=0;i<targets.size();++i)reference[targets[i]]=values[i];
            report.applied=true;return "applied";
        };
        report.reason=attempt();reports.push_back(std::move(report));
    }
    return reports;
}
}

void InterfaceMeshTransfer::protectScalarSupport(const QuadtreeMesh& mesh,
    const InterfaceState& source, const NativeLeafDomain& domain,
    std::vector<std::uint8_t>& protectedCells)
{
    source.validateMesh(mesh);
    domain.validateMeshIdentity(mesh);
    if (source.stage() != InterfaceState::Stage::Ready || !source.geometry()
        || protectedCells.size() != static_cast<std::size_t>(mesh.cellIndexNum()))
        throw std::invalid_argument("Interface mesh protection requires ready geometry and a native cell mask");
    for (int i : domain.cellIndexByRow())
        if (!source.geometry()->closedCellSign(mesh.cells(mesh.cellIndexToLabel[i]).centerCoord()))
            protectedCells[i] = 1;
}

InterfaceMeshTransfer::Result InterfaceMeshTransfer::preparePinnedToAxis(const QuadtreeMesh& sourceMesh,
    const InterfaceState& source, const QuadtreeMesh& target, const NativeLeafDomain& targetDomain,
    const std::vector<std::uint8_t>& protectedSourceCells, InterfaceState& destination)
{
    using Patch = Quantum::Geometry::ActiveQuarticHermitePatch;
    if (source.hasConservativePhase())
        throw std::invalid_argument("conservative inventory requires a conservative mesh transfer");
    source.validateMesh(sourceMesh);
    targetDomain.validateMeshIdentity(target);
    if (source.stage() != InterfaceState::Stage::Ready || &source == &destination
        || destination.stage() != InterfaceState::Stage::Empty || destination.vertexCount() != target.vertIndexNum()
        || !source.geometry() || !source.panels() || !source.geometry()->usesNativeLeaves()
        || !std::holds_alternative<Patch>(source.scalarPatch()) || !targetDomain.axisymmetric()
        || targetDomain.uniformLeaves() || source.pinnedCoordinates().size() != 1)
        throw std::invalid_argument("Interface mesh transfer requires native Q4, pinned RZ and a fresh destination");
    auto& next = InterfaceStateAccess::get(destination);
    InterfacePreparation preparation(next);
    std::vector<std::uint8_t> active(target.cellIndexNum(), 0);
    for (int i : targetDomain.cellIndexByRow()) active[i] = 1;
    auto transferred = ActiveQuarticHermiteTransfer::prepare(sourceMesh, std::get<Patch>(source.scalarPatch()),
        target, active, true, protectedSourceCells);
    next.pinnedCoordinates = source.pinnedCoordinates();
    next.initialContactJets = source.initialContactJets();
    for (int old : source.pinnedVertices()) {
        if (old < 0 || old >= static_cast<int>(transferred.sourceToTargetVertex.size())
            || transferred.sourceToTargetVertex[old] < 0)
            throw std::invalid_argument("Interface mesh transfer lost a pinned contact vertex");
        next.pinnedVertexIndices.push_back(transferred.sourceToTargetVertex[old]);
    }
    auto scalar = std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(target, transferred.patch, true));
    auto chain = LevelSetInterfacePanels::discoverNativeGraphs(target, scalar, true, next.pinnedCoordinates);
    auto panels = std::make_unique<LevelSetInterfacePanels>(
        source.panels()->rebindNative(std::move(chain), transferred.sourceToTargetCell));
    InterfaceMotion::validatePinnedToAxisChain(*panels, next.pinnedCoordinates.front());
    next.activeVertex.assign(target.vertIndexNum(), 0);
    for (int i : targetDomain.cellIndexByRow())
        for (int label : target.cells(target.cellIndexToLabel[i]).cellToVert)
            next.activeVertex[target.verts(label).index()] = 1;
    for (int i = 0; i < target.vertIndexNum(); ++i)
        if (next.activeVertex[i])
            next.freeLevelSet.value[i] = target.verts(target.vertIndexToLabel[i]).supCell >= 0
                ? scalar->value(target.verts(target.vertIndexToLabel[i]).coord()) : transferred.patch.nodes[i].value;
    next.scalarPatch = std::move(transferred.patch);
    next.meshIdentity = targetDomain.meshIdentity();
    next.retainGeometry(std::move(scalar), std::move(panels));
    next.progress = source.progress();
    preparation.complete(InterfaceState::Stage::Ready);
    return {std::move(transferred.sourceToTargetCell), std::move(transferred.sourceToTargetVertex), transferred.report};
}

// Implementation-only access to topology identity. Geometry extraction still
// has a single owner; this only binds its remeshed components to their history.
void InterfaceMeshTransferTopology::retain(const InterfaceTopology& old,InterfaceTopology& next,
    const std::vector<std::vector<int>>& oldCellsByTarget)
{
    const auto require=[](bool ok,const char* why) {
        if(!ok) throw Quantum::Geometry::InterfaceGeometryResolution(
            Quantum::Geometry::InterfaceGeometryResolution::Reason::GraphPairing,-1,why);
    };
        require(old.symmetryPlanes_==next.symmetryPlanes_,"zero-time transfer changed symmetry sections");
        require(old.liquids_.size()==next.liquids_.size() && old.surfaces_.size()==next.surfaces_.size(),
            "zero-time mesh transfer changed component counts");
        std::vector<int> oldLiquidByCell(old.cellPhase_.size(),-1);
        for(std::size_t i=0;i<old.liquids_.size();++i) for(int c:old.liquids_[i].pureCells) {
            require(oldLiquidByCell.at(c)<0,"old pure cell belongs to multiple liquids");
            oldLiquidByCell[c]=static_cast<int>(i);
        }
        std::vector<int> usedLiquids(old.liquids_.size(),0);
        for(auto& liquid:next.liquids_) {
            std::set<int> parents;
            for(int c:liquid.pureCells) for(int previous:oldCellsByTarget.at(c))
                if(oldLiquidByCell.at(previous)>=0) parents.insert(oldLiquidByCell[previous]);
            require(parents.size()==1,"zero-time liquid lineage lacks a unique pure-overlap witness");
            const int parent=*parents.begin();
            require(++usedLiquids[parent]==1,"zero-time mesh transfer would split a liquid identity");
            require(liquid.touchesBoundary==old.liquids_[parent].touchesBoundary,
                "zero-time mesh transfer changed liquid boundary attachment");
            liquid.id=old.liquids_[parent].id;
        }
        require(std::all_of(usedLiquids.begin(),usedLiquids.end(),[](int n){return n==1;}),
            "zero-time mesh transfer lost a liquid identity");
        std::vector<std::set<int>> oldSurfaceCells(old.surfaces_.size());
        for(std::size_t i=0;i<old.surfaces_.size();++i)
            for(const auto& panel:old.surfaces_[i].geometry.panels()) oldSurfaceCells[i].insert(panel.cellIndex);
        // Surface liquidId initially names the fresh extraction's region.
        // Its ordered region vector has not changed; rebind via region index.
        std::vector<int> usedSurfaces(old.surfaces_.size(),0);
        for(auto& surface:next.surfaces_) {
            const int region=surface.geometry.nativeChain().negativeRegion();
            require(region>=0 && static_cast<std::size_t>(region)<next.liquids_.size(),"remeshed surface lost its region");
            surface.liquidId=next.liquids_[region].id;
            std::set<int> sourceCells;
            for(const auto& panel:surface.geometry.panels())
                for(int c:oldCellsByTarget.at(panel.cellIndex)) sourceCells.insert(c);
            std::vector<int> parents;
            for(std::size_t i=0;i<old.surfaces_.size();++i) {
                const auto& previous=old.surfaces_[i];
                if(previous.liquidId!=surface.liquidId || previous.geometry.closed()!=surface.geometry.closed()
                    || previous.geometry.beginEndpointKind()!=surface.geometry.beginEndpointKind()
                    || previous.geometry.endEndpointKind()!=surface.geometry.endEndpointKind()) continue;
                if(std::any_of(sourceCells.begin(),sourceCells.end(),[&](int c){return oldSurfaceCells[i].count(c)!=0;}))
                    parents.push_back(static_cast<int>(i));
            }
            require(parents.size()==1,"zero-time surface lineage is ambiguous or absent");
            const int parent=parents.front();
            require(++usedSurfaces[parent]==1,"zero-time mesh transfer would split a surface identity");
            surface.id=old.surfaces_[parent].id;
            surface.parents=old.surfaces_[parent].parents; // Same physical frame, not a new ancestry edge.
        }
        require(std::all_of(usedSurfaces.begin(),usedSurfaces.end(),[](int n){return n==1;}),
            "zero-time mesh transfer lost a surface identity");
        next.nextLiquidId_=old.nextLiquidId_;next.nextSurfaceId_=old.nextSurfaceId_;next.events_=old.events_;
}

InterfaceMeshTransfer::ConservativeReport InterfaceMeshTransfer::prepareConservative(
    const QuadtreeMesh& oldMesh,const NativeLeafDomain& oldDomain,const InterfaceState& source,
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const std::vector<Vector2D>& normals,
    const std::vector<std::optional<Quantum::Geometry::QuadraticCellCut::Shape>>& curves,
    const ConservativeLimits& limits,InterfaceState& destination,const InterfaceConservativeGeometry::BoundaryInput* targetBoundaries)
{
    const auto require=[](bool ok,const char* why){if(!ok) throw std::invalid_argument(why);};
    require(&source!=&destination && destination.stage()==InterfaceState::Stage::Empty
        && destination.vertexCount()==mesh.vertIndexNum(),"conservative remap requires a fresh target");
    InterfacePreparation preparation(InterfaceStateAccess::get(destination));
    require(source.componentInventories().empty(),
        "split-child material inventory mesh transfer is not yet implemented");
    auto provisional=limits.synchronizeGeometry?std::make_unique<InterfaceState>(mesh):nullptr;
    auto& geometryTarget=provisional?*provisional:destination;
    auto& next=InterfaceStateAccess::get(geometryTarget);
    source.validateMesh(oldMesh);oldDomain.validateMeshIdentity(oldMesh);domain.validateMeshIdentity(mesh);
    require(source.hasTopology() && source.hasConservativePhase() && source.geometry()->usesNativeLeaves()
        && std::holds_alternative<Quantum::Geometry::ActiveHermitePatch>(source.scalarPatch())
        && oldDomain.axisymmetric()==domain.axisymmetric(),"conservative remap requires native cubic topology and inventory");
    require(std::isfinite(limits.scalarDifference) && limits.scalarDifference>=0
        && std::isfinite(limits.geometryInventoryL1) && limits.geometryInventoryL1>=0,"invalid mesh transfer budgets");
    require(!limits.synchronizeGeometry || limits.sharedCurveReference,
        "joint geometry synchronization requires the shared-curve inventory representation");
    require(source.pinnedCoordinates().empty(),"legacy contact metadata requires an explicit remap");
    std::vector<Vector2D> contacts;
    for(const auto& surface:source.topology().surfaces()) {
        using Kind=LevelSetInterfacePanels::EndpointKind;
        const auto& chain=surface.geometry.nativeChain();
        if(surface.geometry.beginEndpointKind()==Kind::DeclaredBoundary)contacts.push_back(chain.frameAt(0,0).point);
        if(surface.geometry.endEndpointKind()==Kind::DeclaredBoundary)contacts.push_back(chain.frameAt(static_cast<int>(chain.graphs().size())-1,1).point);
    }
    require((contacts.empty() && source.topology().symmetryPlanes().empty()) || targetBoundaries,
        "declared contact/section remap requires explicit target boundary laws");
    if(targetBoundaries){
        require(targetBoundaries->domain.geometry().meshIdentity()==domain.meshIdentity()
            && targetBoundaries->domain.geometry().axisymmetric()==domain.axisymmetric()
            && targetBoundaries->domain.geometry().cellIndexByRow()==domain.cellIndexByRow(),"contact remap boundary domain mismatch");
        for(const auto& section:source.topology().symmetryPlanes()) {
            const bool horizontal=section.normalAxis==1;
            require(section.coordinate==(horizontal?mesh.minCoord.y:mesh.minCoord.x)
                || section.coordinate==(horizontal?mesh.maxCoord.y:mesh.maxCoord.x),
                "section remap requires the same physical outer reflection plane");
            bool found=false;
            for(const auto face:domain.faces(mesh))if(face.horizontal==horizontal
                && (horizontal?face.center.y:face.center.x)==section.coordinate) {
                found=true;
                require(face.outwardSign!=0 && targetBoundaries->domain.boundary(face.nativeRow)==InterfaceTransportDomain::Boundary::Symmetry
                    && (!targetBoundaries->portPhase || targetBoundaries->portPhase(face.nativeRow)==0),
                    "section remap requires an impermeable symmetry law on every target face");
            }
            require(found,"section remap lost its physical target boundary");
        }
        for(auto p:contacts){bool fixed=false;
            for(const auto face:domain.faces(mesh))if((face.negativeCell>=0)!=(face.positiveCell>=0)
                && targetBoundaries->domain.boundary(face.nativeRow)==InterfaceTransportDomain::Boundary::StationaryWall){
                const auto& edge=face.horizontal?mesh.hEdges(face.edgeLabel):mesh.vEdges(face.edgeLabel);
                const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
                fixed|=(face.horizontal?p.y==a.y:p.x==a.x) && p.x>=std::min(a.x,b.x) && p.x<=std::max(a.x,b.x)
                    && p.y>=std::min(a.y,b.y) && p.y<=std::max(a.y,b.y);
            }
            require(fixed,"contact remap endpoint has no stationary target wall");
        }
    }
    InterfacePhaseGeometry::validateLayout(oldMesh,oldDomain,*source.geometry());
    std::vector<std::uint8_t> oldAssigned(oldMesh.cellIndexNum(),0);
    for(const auto& liquid:source.topology().liquids()) for(int c:liquid.cells) {
        require(!oldAssigned[c],"source inventory cell contains multiple liquids; refine before single-donor remap");
        oldAssigned[c]=1;
    }
    const auto clock=source.progress();const auto& oldInventory=source.conservativePhase();
    require(oldInventory.progress().completedSteps==clock.completedSteps && oldInventory.progress().endTime==clock.endTime,
        "source geometry/inventory clocks differ");
    std::vector<std::uint8_t> active(mesh.cellIndexNum(),0);for(int c:domain.cellIndexByRow()) active[c]=1;
    auto scalar=ActiveHermiteTransfer::prepare(oldMesh,std::get<Quantum::Geometry::ActiveHermitePatch>(source.scalarPatch()),mesh,active,domain.axisymmetric());
    require((limits.reextendPureSectors?scalar.maxNonExtensionDifferenceBound:scalar.maxScalarDifferenceBound)
        <=limits.scalarDifference,"scalar remap exceeds caller error budget");
    std::vector<int> eventCells,eventSource(mesh.cellIndexNum(),-1);
    std::vector<unsigned char> modeled(oldMesh.cellIndexNum(),0),retained(oldMesh.cellIndexNum(),0);
    for(int c:source.eventCurvatureCells())modeled.at(c)=1;
    for(int c:domain.cellIndexByRow())for(int old:scalar.sourceCellsByTarget[c])if(modeled[old]) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto& oldCell=oldMesh.cells(oldMesh.cellIndexToLabel[old]);
        const auto lower=mesh.verts(cell.cellToVert[0]).coord(),upper=mesh.verts(cell.cellToVert[2]).coord();
        const auto box=source.geometry()->nativeCellBounds(old);
        require(scalar.sourceCellsByTarget[c].size()==1 && cell.level()==oldCell.level()
            && lower.x==box[0].x && lower.y==box[0].y && upper.x==box[1].x && upper.y==box[1].y,
            "terminal model AMR must retain the physical finest leaf; only surrounding support may change");
        require(!retained[old],"terminal model leaf has multiple transfer targets");
        retained[old]=1;eventCells.push_back(c);eventSource[c]=old;
    }
    for(int c:source.eventCurvatureCells())require(retained[c],"terminal model leaf was lost during AMR");
    std::vector<int> extensions;
    for(int c:domain.cellIndexByRow())if(scalar.pureSectorExtension[c]) {
        require(targetBoundaries,"wet-sector remap requires explicit target boundary laws");
        for(int old:scalar.sourceCellsByTarget[c])require(oldInventory.liquidVolumes()[old]==
            (scalar.pureSectorExtension[c]<0?oldInventory.capacities()[old]:0),
            "pure geometric sector has non-pure independent inventory");
        extensions.push_back(c);
    }
    InterfaceGeometryReconstruction::initialize(mesh,std::move(scalar.patch),domain.axisymmetric(),contacts,geometryTarget,
        source.topology().symmetryPlanes());
    InterfaceMeshTransferTopology::retain(source.topology(),*next.topology,scalar.sourceCellsByTarget);
    next.eventCurvatureCells=std::move(eventCells);
    InterfacePhaseGeometry::MappedFacePartitions partitions(mesh,domain,geometryTarget.geometry());
    auto measure=InterfaceCellMeasure::evaluateVolumes(mesh,domain,partitions,geometryTarget);
    // The retained zero set supplies subcell shape, never a replacement for
    // independently transported parent C. Canonicalize only certified pure
    // geometric references; raw integrals remain in the discrepancy report.
    auto reference=measure.negativeVolumeByCell;
    for(int c:domain.cellIndexByRow()) {
        const auto phase=next.topology->cellPhase()[c];
        if(phase==1) reference[c]=0;
        else if(phase==-1) reference[c]=domain.cellVolumes()[domain.cellRowByIndex()[c]];
    }
    // A whole-surface auxiliary fit must not reinterpret the terminal cap.
    // The native scalar and its original error budget remain authoritative.
    auto shared=limits.sharedCurveReference && source.eventCurvatureCells().empty()?sharedReferences(oldMesh,source,mesh,domain,scalar.sourceCellsByTarget,reference,limits.synchronizeGeometry):
        std::vector<ConservativeReport::SharedCurve>{};
    auto inventory=InterfaceConservativePhase::transferMesh(oldMesh,oldDomain,oldInventory,mesh,domain,normals,curves,reference);
    next.progress=clock;next.conservativePhase.emplace(std::move(inventory));
    for(int c:next.eventCurvatureCells)require(next.conservativePhase->liquidVolumes()[c]==oldInventory.liquidVolumes()[eventSource[c]],
        "terminal model leaf inventory changed in an identity transfer");
    const bool synchronize=provisional && std::any_of(shared.begin(),shared.end(),[](const auto& curve){return curve.applied;});
    if(synchronize) InterfaceConservativeGeometry::resynchronize(mesh,domain,geometryTarget,destination,4,targetBoundaries);
    else if(provisional) InterfaceStateAccess::get(destination)=std::move(next);
    const auto& finalInventory=destination.conservativePhase();
    const auto& finalTopology=destination.topology();
    ConservativeReport report;report.sourceGeometry=source.geometry();
    report.sharedCurves=std::move(shared);
    report.scalarDifferenceBound=std::move(scalar.scalarDifferenceBound);report.maxScalarDifferenceBound=scalar.maxScalarDifferenceBound;
    report.pureSectorExtensionCells=std::move(extensions);report.maxNonExtensionDifferenceBound=scalar.maxNonExtensionDifferenceBound;
    report.geometricVolumeByCell=std::move(measure.negativeVolumeByCell);report.geometricVolume=measure.summedNegativeVolume;
    if(synchronize) {
        const auto& synchronization=destination.conservativeGeometry();
        report.geometrySynchronized=true;report.maxSynchronizationScalarChange=synchronization.maxNodalScalarChange;
        report.maxSynchronizationInterfaceShift=synchronization.maxZeroTimeInterfaceShift;
        report.geometricVolumeByCell=synchronization.geometricVolumeByCell;report.geometricVolume=synchronization.geometricVolume;
    }
    report.retainedLiquids=finalTopology.liquids().size();report.retainedSurfaces=finalTopology.surfaces().size();
    report.retainedEventCells=destination.eventCurvatureCells().size();
    std::vector<std::uint8_t> assigned(mesh.cellIndexNum(),0);
    for(const auto& liquid:finalTopology.liquids()) {
        ConservativeReport::Component component;component.liquidId=liquid.id;
        for(int c:liquid.cells) {
            require(!assigned[c],"coarsening put multiple liquid components in one inventory cell");assigned[c]=1;
            component.inventoryVolume+=finalInventory.liquidVolumes()[c];
            component.geometricVolume+=report.geometricVolumeByCell[c];
            component.geometryInventoryL1+=std::abs(report.geometricVolumeByCell[c]-finalInventory.liquidVolumes()[c]);
        }
        report.components.push_back(component);
    }
    for(int c:domain.cellIndexByRow()) {
        const auto volume=finalInventory.liquidVolumes()[c],difference=std::abs(report.geometricVolumeByCell[c]-volume);
        require(std::isfinite(difference),"nonfinite remapped geometry/inventory discrepancy");
        report.inventoryVolume+=volume;report.geometryInventoryL1+=difference;
        report.maxFractionDifference=std::max(report.maxFractionDifference,difference/finalInventory.capacities()[c]);
        if(!assigned[c]) report.unassignedInventoryVolume+=volume;
    }
    require(report.geometryInventoryL1<=limits.geometryInventoryL1,"remapped geometry/inventory mismatch exceeds caller budget");
    preparation.complete(InterfaceState::Stage::Ready);
    return report;
}
