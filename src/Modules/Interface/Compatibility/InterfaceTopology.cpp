#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Core/Geometry/InterfaceGeometryResolution.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <sstream>

namespace
{
bool overlaps(const std::vector<int>& a, const std::vector<int>& b)
{
    auto i=a.begin(), j=b.begin();
    while (i!=a.end() && j!=b.end()) {
        if (*i==*j) return true;
        if (*i<*j) ++i; else ++j;
    }
    return false;
}
void require(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
}

void InterfaceGeometryReconstruction::initialize(const QuadtreeMesh& mesh, InterfaceScalarPatch patch,
    bool axisymmetric, const std::vector<Vector2D>& contacts, InterfaceState& destination,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes)
{ prepare(mesh, std::move(patch), axisymmetric, contacts, nullptr, nullptr, destination,symmetryPlanes); }

void InterfaceGeometryReconstruction::prepareManufactured(const QuadtreeMesh& mesh, InterfaceScalarPatch patch,
    bool axisymmetric, const std::vector<Vector2D>& contacts, const InterfaceState& previous,
    const InterfaceStepTime& time, InterfaceState& destination)
{
    require(!previous.hasConservativePhase(),
        "manufactured geometry cannot discard conservative inventory");
    prepare(mesh, std::move(patch), axisymmetric, contacts, &previous, &time, destination);
}

void InterfaceGeometryReconstruction::prepare(const QuadtreeMesh& mesh, InterfaceScalarPatch patch,
    bool axisymmetric, const std::vector<Vector2D>& contacts, const InterfaceState* previous,
    const InterfaceStepTime* time, InterfaceState& destination,
    const std::vector<LevelSetInterfacePanels::SymmetryPlane>& symmetryPlanes,
    std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> retainedScalar)
{
    require(destination.stage()==InterfaceState::Stage::Empty, "topology reconstruction requires a fresh destination");
    if (previous) {
        previous->validateMesh(mesh);
        require(previous->geometry() && previous->geometry()->radialAxisSymmetry()==axisymmetric,
            "manufactured sequence cannot change its XY/RZ interpretation");
        const auto progress=previous->progress();
        require(previous->hasTopology() && time && time->index==progress.completedSteps+1 &&
            time->begin==progress.endTime && std::isfinite(time->end) && time->end>time->begin &&
            std::isfinite(time->dt) && time->dt>0 &&
            std::abs((time->end-time->begin)-time->dt)<=16*std::numeric_limits<double>::epsilon()*
                std::max({std::abs(time->dt),std::abs(time->begin),std::abs(time->end)}),
            "manufactured sequence requires consistent external clock and prior complete topology");
    }
    auto& state=InterfaceStateAccess::get(destination);
    InterfacePreparation preparation(state);
    auto source=retainedScalar?std::move(retainedScalar):std::visit([&](const auto& data) {
        return std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
            LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,data,axisymmetric));
    },patch);
    require(source->usesNativeLeaves() && source->radialAxisSymmetry()==axisymmetric &&
        source->nativeCellCount()==static_cast<std::size_t>(mesh.cellIndexNum()),"prepared scalar layout does not match topology");
    const auto& retainedPlanes=previous?previous->topology().symmetryPlanes():symmetryPlanes;
    auto discovered=LevelSetInterfacePanels::discoverNativeComponents(mesh,source,axisymmetric,contacts,retainedPlanes);
    auto topology=std::make_unique<InterfaceTopology>();
    topology->symmetryPlanes_=retainedPlanes;
    topology->cellPhase_=std::move(discovered.cellPhase);
    const auto* old=previous ? &previous->topology() : nullptr;
    if (old) {
        topology->nextLiquidId_=old->nextLiquidId_; topology->nextSurfaceId_=old->nextSurfaceId_;
        topology->events_=old->events_; // Ancestry remains module-owned after the old handle is replaced.
    }
    for (auto& region:discovered.negativeRegions)
        topology->liquids_.push_back({0,std::move(region.cells),std::move(region.pureCells),region.touchesNonAxisBoundary,{}});

    if (!old) for (auto& liquid:topology->liquids_) liquid.id=topology->nextLiquidId_++;
    else {
        // Bipartite overlap groups permit 1:N, N:1 and N:M. No array position
        // is used as persistent identity, and no global volume match is used.
        const int n=static_cast<int>(old->liquids_.size()), m=static_cast<int>(topology->liquids_.size());
        std::vector<int> parent(n+m), degree(n+m,0); std::iota(parent.begin(),parent.end(),0);
        const auto root=[&](int i) { while (parent[i]!=i) i=parent[i]; return i; };
        for (int i=0;i<n;++i) for (int j=0;j<m;++j)
            if (overlaps(old->liquids_[i].pureCells,topology->liquids_[j].pureCells)) {
                parent[root(n+j)]=root(i); ++degree[i]; ++degree[n+j];
            }
        if(n==m && previous->hasConservativePhase()
            && std::any_of(degree.begin(),degree.end(),[](int d){return d==0;})
            && (!previous->eventCurvatureCells().empty()
                || std::any_of(old->events_.begin(),old->events_.end(),[](const auto& event){return !event.interiorNecks.empty();}))) {
            // Explicit modeled children may have no pure liquid cell. They
            // still have distinct, connected intervals on the RZ axis.
            // Unique positive overlap retains identity only for an unchanged
            // set of components; it cannot create a split, merge or new ID.
            using Interval=std::pair<double,double>;
            const auto interval=[](const auto& chain)->std::optional<Interval> {
                if(!chain.axisymmetric() || chain.closed())return {};
                std::array<double,2> z{};
                for(bool end:{false,true}) {
                    const auto* plane=chain.endpointSymmetryPlane(end);
                    if((end?chain.endEndpointKind():chain.beginEndpointKind())!=LevelSetInterfacePanels::EndpointKind::Axis
                        && !(plane && plane->normalAxis==1))return {};
                    z[end]=chain.frameAt(end?static_cast<int>(chain.graphs().size())-1:0,end?1:0).point.y;
                }
                if(z[0]==z[1])return {};return std::minmax(z[0],z[1]);
            };
            std::vector<std::optional<Interval>> from(n),to(m);bool witnessed=true;
            for(const auto& surface:old->surfaces_) {
                const auto liquid=std::find_if(old->liquids_.begin(),old->liquids_.end(),
                    [&](const auto& value){return value.id==surface.liquidId;});
                const int i=static_cast<int>(liquid-old->liquids_.begin());
                if(i==n || from[i]){witnessed=false;break;}from[i]=interval(surface.geometry.nativeChain());
            }
            for(const auto& chain:discovered.curves) {
                const int j=chain.negativeRegion();
                if(j<0 || j>=m || to[j]){witnessed=false;break;}to[j]=interval(chain);
            }
            witnessed=witnessed && std::all_of(from.begin(),from.end(),[](const auto& x){return x.has_value();})
                && std::all_of(to.begin(),to.end(),[](const auto& x){return x.has_value();});
            std::vector<int> match(n,-1),targets(m,0);
            if(witnessed) {
                for(int i=0;i<n;++i)for(int j=0;j<m;++j)
                    if(std::min(from[i]->second,to[j]->second)>std::max(from[i]->first,to[j]->first)) {
                        if(match[i]>=0)witnessed=false;match[i]=j;++targets[j];
                    }
                for(int i=0;i<n;++i) {
                    const int j=match[i];
                    if(j<0 || targets[j]!=1){witnessed=false;continue;}
                    if((degree[i] || degree[n+j]) && (degree[i]!=1 || degree[n+j]!=1 || root(i)!=root(n+j)))witnessed=false;
                    for(int k=i+1;k<n;++k)
                        if(std::min(from[i]->second,from[k]->second)>std::max(from[i]->first,from[k]->first))witnessed=false;
                }
                for(int j=0;j<m;++j)for(int k=j+1;k<m;++k)
                    if(std::min(to[j]->second,to[k]->second)>std::max(to[j]->first,to[k]->first))witnessed=false;
            }
            if(witnessed)for(int i=0;i<n;++i)if(degree[i]==0) {
                parent[root(n+match[i])]=root(i);degree[i]=degree[n+match[i]]=1;
            }
        }
        if (std::any_of(degree.begin(),degree.end(),[](int d){return d==0;})) {
            std::ostringstream detail;detail<<"lineage lacks shared pure-cell witnesses; previous="<<n<<", candidate="<<m;
            for (int i=0;i<n;++i) if (degree[i]==0) detail<<"; unmatched previous="<<i
                <<", cells="<<old->liquids_[i].cells.size()<<", pure="<<old->liquids_[i].pureCells.size();
            for (int j=0;j<m;++j) if (degree[n+j]==0) {
                const auto& liquid=topology->liquids_[j];
                detail<<"; unmatched candidate="<<j<<", cells="<<liquid.cells.size()<<", pure="<<liquid.pureCells.size()<<", first cells=[";
                for (std::size_t k=0;k<std::min<std::size_t>(liquid.cells.size(),12);++k) detail<<liquid.cells[k]<<',';
                detail<<']';
            }
            throw Quantum::Geometry::InterfaceGeometryResolution(
                Quantum::Geometry::InterfaceGeometryResolution::Reason::GraphPairing,-1,
                detail.str());
        }
        std::map<int,std::pair<std::vector<int>,std::vector<int>>> groups;
        for (int i=0;i<n;++i) groups[root(i)].first.push_back(i);
        for (int j=0;j<m;++j) groups[root(n+j)].second.push_back(j);
        for (const auto& entry:groups) {
            const auto& parents=entry.second.first; const auto& children=entry.second.second;
            if (parents.size()==1 && children.size()==1) {
                topology->liquids_[children.front()].id=old->liquids_[parents.front()].id; continue;
            }
            InterfaceTopologyEvent event{parents.size()==1 ? InterfaceTopologyEvent::Kind::Split :
                children.size()==1 ? InterfaceTopologyEvent::Kind::Merge : InterfaceTopologyEvent::Kind::Reconnection,
                time->begin,time->end,{},{}};
            for (int i:parents) event.parents.push_back(old->liquids_[i].id);
            for (int j:children) {
                auto& liquid=topology->liquids_[j]; liquid.id=topology->nextLiquidId_++;
                event.children.push_back(liquid.id);
            }
            topology->events_.push_back(std::move(event));
        }
    }
    std::vector<long double> volumes(topology->liquids_.size(),0);
    for (std::size_t i=0;i<volumes.size();++i) volumes[i]=discovered.negativeRegions[i].fixedBoundaryVolume;
    for (auto& chain:discovered.curves) {
        const int region=chain.negativeRegion();
        require(region>=0 && static_cast<std::size_t>(region)<topology->liquids_.size(),"curve has no negative region");
        const auto liquidId=topology->liquids_[region].id;
        LevelSetInterfacePanels panels(std::move(chain));
        long double area=0, volume=0;
        const bool controlled=std::all_of(panels.nativeChain().graphs().begin(),panels.nativeChain().graphs().end(),
            [](const auto& graph){return graph.controlledIntegration;});
        for (const auto& panel:panels.panels()) for (const auto& q:panel.quadratureNodes) {
            area+=q.physicalWeight;
            const long double r=q.frame.point.x;
            const long double contrib=(axisymmetric ? std::acos(-1.0L)*r*r : r)*q.frame.normal.x*q.weight;
            if (!controlled) {
                volume+=contrib;
            }
        }
        if (controlled) for (const auto& graph:panels.nativeChain().graphs()) {
            const auto moments=LevelSetMappedInterface::nativeMoments(*panels.nativeChain().source(),graph,axisymmetric);
            volume+=moments.radial;
        }
        volumes[region]+=volume;
        topology->surfaces_.push_back({0,liquidId,std::move(panels),static_cast<double>(area),{}});
    }
    // Establish both directions BEFORE assigning surface IDs. An old curve
    // splitting into two must not donate its ID to whichever is visited first.
    std::vector<std::vector<int>> predecessors(topology->surfaces_.size());
    std::vector<int> successorCount(old ? old->surfaces_.size() : 0, 0);
    if (old) for (std::size_t j=0;j<topology->surfaces_.size();++j) {
        const auto& target=topology->surfaces_[j];
        std::vector<std::uint64_t> parentLiquids{target.liquidId};
        for (const auto& event:topology->events_)
            if (event.end==time->end &&
                std::find(event.children.begin(),event.children.end(),target.liquidId)!=event.children.end())
                parentLiquids=event.parents;
        std::vector<int> eligible;
        for (std::size_t i=0;i<old->surfaces_.size();++i) {
            const auto& surface=old->surfaces_[i];
            if (std::find(parentLiquids.begin(),parentLiquids.end(),surface.liquidId)==parentLiquids.end()) continue;
            if (surface.geometry.closed()!=target.geometry.closed() ||
                surface.geometry.beginEndpointKind()!=target.geometry.beginEndpointKind() ||
                surface.geometry.endEndpointKind()!=target.geometry.endEndpointKind()) continue;
            const auto& from=surface.geometry.nativeChain();const auto& to=target.geometry.nativeChain();
            bool sameSections=true;
            for (bool end:{false,true}) {
                const auto* a=from.endpointSymmetryPlane(end);const auto* b=to.endpointSymmetryPlane(end);
                sameSections &= (!a && !b) || (a && b && *a==*b);
            }
            if (!sameSections) continue;
            eligible.push_back(static_cast<int>(i));
        }
        const bool soleTarget=std::count_if(topology->surfaces_.begin(),topology->surfaces_.end(),
            [&](const auto& surface){return surface.liquidId==target.liquidId;})==1;
        const bool onePerParent=std::all_of(parentLiquids.begin(),parentLiquids.end(),[&](auto id) {
            return std::count_if(old->surfaces_.begin(),old->surfaces_.end(),
                [&](const auto& surface){return surface.liquidId==id;})==1;
        });
        if (eligible.size()==1 || (soleTarget && onePerParent)) predecessors[j]=eligible;
        else {
            std::vector<int> cells;
            for (const auto& panel:target.geometry.panels()) cells.push_back(panel.cellIndex);
            std::sort(cells.begin(),cells.end());
            for (int i:eligible)
                for (const auto& panel:old->surfaces_[i].geometry.panels())
                    if (std::binary_search(cells.begin(),cells.end(),panel.cellIndex)) {
                        predecessors[j].push_back(i); break;
                    }
        }
        if (predecessors[j].empty())
            throw Quantum::Geometry::InterfaceGeometryResolution(
                Quantum::Geometry::InterfaceGeometryResolution::Reason::GraphPairing,-1,
                "manufactured surface lineage is unresolved; refine or sample more closely");
        for (int i:predecessors[j]) ++successorCount[i];
    }
    if (old && std::any_of(successorCount.begin(),successorCount.end(),[](int count){return count==0;}))
        throw Quantum::Geometry::InterfaceGeometryResolution(
            Quantum::Geometry::InterfaceGeometryResolution::Reason::GraphPairing,-1,
            "manufactured surface disappeared without resolved lineage");
    for (std::size_t j=0;j<topology->surfaces_.size();++j) {
        auto& surface=topology->surfaces_[j];
        for (int i:predecessors[j]) surface.parents.push_back(old->surfaces_[i].id);
        if (predecessors[j].size()==1 && successorCount[predecessors[j].front()]==1 &&
            old->surfaces_[predecessors[j].front()].liquidId==surface.liquidId)
            surface.id=old->surfaces_[predecessors[j].front()].id;
        else surface.id=topology->nextSurfaceId_++;
    }
    for (std::size_t i=0;i<volumes.size();++i) {
        require(std::isfinite(static_cast<double>(volumes[i])) && volumes[i]>0,
            "negative region has nonpositive or nonfinite volume");
        topology->liquids_[i].volume=static_cast<double>(volumes[i]);
    }
    state.meshIdentity=std::make_shared<const NativeMeshIdentity>(mesh);
    state.activeVertex.assign(mesh.vertIndexNum(),0);
    std::visit([&](const auto& data) {
        const std::array<Vector2D,4> corners{{{0,0},{1,0},{1,1},{0,1}}};
        for (int c=0;c<mesh.cellIndexNum();++c) if (data.activeCells[c]) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);
            for (int k=0;k<4;++k) {
                const int v=mesh.verts(cell.cellToVert[k]).index();
                if (state.activeVertex[v]) continue;
                state.activeVertex[v]=1;
                // Hanging-node payloads may be unused/NaN. Read the actual
                // constrained scalar, not the caller's ignored input value.
                state.freeLevelSet.value[v]=source->scalarJetAtNativeCellChart(c,
                    Quantum::Geometry::cellChartFromLocal(corners[k])).value;
            }
        }
    },patch);
    state.scalarPatch=std::move(patch); state.source=std::move(source); state.topology=std::move(topology);
    if (time) state.progress={time->index,time->end};
    preparation.complete(InterfaceState::Stage::Ready);
}
