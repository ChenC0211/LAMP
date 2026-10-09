#include "Modules/Interface/Compatibility/InterfaceNeckCap.h"
#include "Modules/Interface/Compatibility/InterfaceTopology.h"
#include "Modules/Interface/InterfaceCellMeasure.h"
#include "Modules/Interface/InterfacePhaseGeometry.h"
#include "Modules/Interface/Compatibility/InterfaceCurvature.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Core/Geometry/LevelSetMappedSegment.h"
#include "Core/Geometry/CellChartPoint.h"
#include "Core/Geometry/Detail/ExactExpansion.h"
#include "Core/Mesh/CartesianMeasure.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace
{
void require(bool ok,const char* why)
{if(!ok)throw std::runtime_error(std::string("[InterfaceNeckCap] ")+why);}
struct TerminalEnd { bool end,detached; };
TerminalEnd terminalEnd(const InterfaceState& state,double upper)
{
    require(state.hasConservativePhase() && state.hasTopology() && state.topology().surfaces().size()==1,
        "terminal model requires one conservative column surface");
    const auto& chain=state.topology().surfaces()[0].geometry.nativeChain();
    require(!chain.closed() && chain.axisymmetric(),"terminal model requires an open RZ column or cap");
    const int count=static_cast<int>(chain.graphs().size());
    const bool end=chain.frameAt(count-1,1).point.y>chain.frameAt(0,0).point.y;
    const auto* section=chain.endpointSymmetryPlane(end);bool detached=false;
    for(const auto& event:state.topology().events())if(event.kind==InterfaceTopologyEvent::Kind::SectionDetachment)
        for(const auto& plane:event.detachedSections)detached|=plane.normalAxis==1 && plane.coordinate==upper;
    require(section?(section->normalAxis==1 && section->coordinate==upper && !detached):
        (detached && (end?chain.endEndpointKind():chain.beginEndpointKind())==LevelSetInterfacePanels::EndpointKind::Axis),
        "terminal model lost its declared section or reattached a detached cap");
    return {end,detached};
}
double sectionRadius(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,int cell,double z)
{
    const auto box=source->nativeCellBounds(cell);
    // The native mapped port retains the shared-edge restriction. The legacy
    // logical-grid edge API deliberately does not accept native snapshots.
    LevelSetMappedSegment section(source,{box[0].x,z},{box[1].x,z});
    require(!section.ranges().empty()&&section.ranges().front().negative,"source section does not start in axial liquid");
    long double radius=0;int negativeIntervals=0;
    for(const auto& phase:section.ranges())if(phase.negative) {
        ++negativeIntervals;
        radius=box[1].x*phase.length;
    }
    require(negativeIntervals==1,"section must have one radial liquid interval");
    const Vector2D p{static_cast<double>(radius),z};
    require(p.x>0&&p.x<box[1].x,"source neck is not strictly inside the first radial cell");
    const auto jet=source->scalarJetInNativeCell(cell,p);
    require(jet.dx>0,"source neck has the wrong liquid orientation");
    return p.x;
}
// Independent observation of the sampled native scalar, rather than merely
// reporting the model's polynomial integrals. Split at its ACTUAL axis root,
// then integrate actual radial phase ranges with GL8 in z. This is a numerical
// consistency check, not an interval bound on the complete 2D zero set.
long double nativeVolume(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,int cell,int subdivisions,bool reverse=false)
{
    constexpr double nodes[]{-.96028985649753623168,-.79666647741362673959,-.52553240991632898582,-.18343464249564980494,
        .18343464249564980494,.52553240991632898582,.79666647741362673959,.96028985649753623168};
    constexpr double weights[]{.10122853629037625915,.22238103445337447054,.31370664587788728734,.36268378337836198297,
        .36268378337836198297,.31370664587788728734,.22238103445337447054,.10122853629037625915};
    const auto box=source->nativeCellBounds(cell);const double hr=box[1].x,h=box[1].y-box[0].y;
    LevelSetMappedSegment axis(source,{0,box[0].y},{0,box[1].y});
    require(!axis.ranges().empty()&&(reverse?axis.ranges().back().negative:axis.ranges().front().negative),
        "local cap does not start in axial liquid");
    long double volume=0;int negativeIntervals=0;
    for(const auto& interval:axis.ranges())if(interval.negative) {
        ++negativeIntervals;
        // Ordered complete phase coverage plus a single initial negative
        // range establishes this endpoint identity. Do not recover exact zero
        // by subtracting separately rounded from-end offsets and lengths.
        const long double lo=0;
        for(int sub=0;sub<subdivisions;++sub)for(int k=0;k<8;++k) {
            const long double t=lo+interval.length*(sub+(nodes[k]+1)/2)/subdivisions;
            const double z=static_cast<double>(reverse?box[1].y-h*t:box[0].y+h*t);
            LevelSetMappedSegment radial(source,{0,z},{hr,z});
            require(!radial.ranges().empty()&&radial.ranges().front().negative,"local slice does not start in axial liquid");
            long double area=0;int radialIntervals=0;
            for(const auto& phase:radial.ranges())if(phase.negative) {
                ++radialIntervals;
                const long double a=0;
                area+=(CartesianMeasure::twoPi/2)*hr*hr*phase.length*(2*a+phase.length);
            }
            require(radialIntervals==1,"local cap lacks one liquid radial slice");
            volume+=weights[k]*area*h*interval.length/(2*subdivisions);
        }
    }
    require(negativeIntervals==1,"local cap lacks one liquid axial interval");
    return volume;
}
}
InterfaceNeckCap::InterfaceNeckCap(const InterfaceState& source,Quantum::Geometry::AxisymmetricNeckCap model)
    :source_(source.geometry()),
     sourcePanels_(std::make_shared<LevelSetInterfacePanels>(source.topology().surfaces().front().geometry)),
     sourceSurfaceId_(source.topology().surfaces().front().id),liquidId_(source.topology().surfaces().front().liquidId),
     inventory_(source.conservativePhase()),model_(std::move(model)){}

InterfaceNeckCap::MaterialTransfer InterfaceNeckCap::materialTransfer() const
{
    require(complete_ && sourcePanels_,"material closure requires an unmoved cap proposal");
    namespace Exact=Quantum::Geometry::Detail::ExactExpansion;
    MaterialTransfer result;result.source=sourcePanels_;
    const auto& surface=complete_->topology().surfaces().front();
    result.target=std::make_shared<LevelSetInterfacePanels>(surface.geometry);
    result.cells=cells_;result.sourceSurfaceId=sourceSurfaceId_;result.targetSurfaceId=surface.id;
    result.liquidId=liquidId_;result.clock=complete_->progress();result.beginZ=beginZ_;
    result.upperSection=source_->nativeCellBounds(cells_[1])[1].y;
    require(surface.liquidId==liquidId_ && surface.parents==std::vector<std::uint64_t>{sourceSurfaceId_},
        "material cap lost its unique parent liquid/surface");
    const auto& old=result.source->panels();const auto& next=result.target->panels();
    std::vector<std::size_t> local[2],remote[2];
    const auto partition=[&](const auto& panels,int side) {
        for(std::size_t i=0;i<panels.size();++i) {
            const auto& panel=panels[i];
            require(std::isfinite(panel.physicalMeasure) && panel.physicalMeasure>0,"invalid event material area");
            const bool modeled=panel.cellIndex==cells_[0] || panel.cellIndex==cells_[1];
            (modeled?local[side]:remote[side]).push_back(i);
        }
        require(!local[side].empty() && !remote[side].empty(),"material event must be a proper local surface interval");
        require(local[side].back()-local[side].front()+1==local[side].size()
            && (local[side].front()==0 || local[side].back()+1==panels.size()),
            "material cap is not one endpoint interval");
    };
    partition(old,0);partition(next,1);
    require((local[0].front()==0)==(local[1].front()==0) && remote[0].size()==remote[1].size(),
        "material cap changed the remote partition or orientation");
    result.weights.resize(old.size());
    for(std::size_t k=0;k<remote[0].size();++k) {
        const auto i=remote[0][k],j=remote[1][k];const auto& a=old[i];const auto& b=next[j];
        require(a.cellIndex==b.cellIndex && a.coordinateIsX==b.coordinateIsX
            && a.measure==b.measure && a.physicalMeasure==b.physicalMeasure
            && a.beginFrame.hasCellChartPoint && b.beginFrame.hasCellChartPoint
            && a.endFrame.hasCellChartPoint && b.endFrame.hasCellChartPoint
            && Quantum::Geometry::cellChartEqual(a.beginFrame.cellChartPoint,b.beginFrame.cellChartPoint)
            && Quantum::Geometry::cellChartEqual(a.endFrame.cellChartPoint,b.endFrame.cellChartPoint),
            "material cap changed a remote native panel");
        result.weights[i].push_back({j,1});
    }
    const auto cumulative=[](const auto& panels,const auto& indices) {
        std::vector<Exact::Expansion> area(1);
        for(auto i:indices)area.push_back(Exact::add(area.back(),Exact::fromDouble(panels[i].physicalMeasure)));
        return area;
    };
    const auto a=cumulative(old,local[0]),b=cumulative(next,local[1]);
    result.sourceLocalArea=Exact::estimate(a.back());result.targetLocalArea=Exact::estimate(b.back());
    // Compare A_i/A_total with B_j/B_total by exact cross-products. A tiny
    // retained pole panel must not disappear when its cumulative area rounds
    // to the total; only the final positive overlap ratio is rounded.
    const auto less=[](const auto& x,const auto& y){return Exact::sign(Exact::subtract(x,y))<0;};
    std::vector<bool> covered(next.size(),false);
    for(auto j:remote[1])covered[j]=true;
    for(std::size_t i=0;i<local[0].size();++i) {
        const auto lo=Exact::multiply(a[i],b.back()),hi=Exact::multiply(a[i+1],b.back());
        const auto width=Exact::subtract(hi,lo);long double sum=0;
        for(std::size_t j=0;j<local[1].size();++j) {
            const auto start=Exact::multiply(b[j],a.back()),end=Exact::multiply(b[j+1],a.back());
            const auto left=less(lo,start)?start:lo,right=less(hi,end)?hi:end;
            const auto overlap=Exact::subtract(right,left);
            if(Exact::sign(overlap)<=0)continue;
            const long double share=static_cast<long double>(Exact::estimate(overlap))/Exact::estimate(width);
            require(std::isfinite(share) && share>0 && share<=1,"unrepresentable local material overlap");
            result.weights[local[0][i]].push_back({local[1][j],share});sum+=share;covered[local[1][j]]=true;
        }
        require(std::abs(sum-1)<=128*std::numeric_limits<long double>::epsilon()*result.weights[local[0][i]].size(),
            "local event material row does not conserve its integral");
    }
    require(std::all_of(covered.begin(),covered.end(),[](bool x){return x;}),"event material closure lost a target panel");
    return result;
}

std::vector<int> InterfaceNeckCap::terminalResolutionCells(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& state,double upper,int maximumLevel,double maxCurvatureCells)
{
    state.validateMesh(mesh);domain.validateMeshIdentity(mesh);
    require(domain.axisymmetric() && upper==mesh.maxCoord.y && maximumLevel>=0 && maximumLevel<=24
        && std::isfinite(maxCurvatureCells) && maxCurvatureCells>0,"invalid terminal resolution scope");
    const auto terminal=terminalEnd(state,upper);const bool end=terminal.end;
    const auto& chain=state.topology().surfaces()[0].geometry.nativeChain();
    const int count=static_cast<int>(chain.graphs().size());std::vector<int> cells;
    for(int g=end?count-1:0;g>=0 && g<count;g+=end?-1:1) {
        const int c=chain.graphs()[g].nativeCellIndex;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto size=cell.size();
        require(cell.level()<=maximumLevel,"terminal source exceeds the declared finest level");
        const double finestWidth=std::ldexp(std::max(size.x,size.y),static_cast<int>(cell.level())-maximumLevel);
        bool outside=true;
        for(double f:{0.,.5,1.}) {
            const auto frame=chain.frameAt(g,f);
            require(terminal.detached || frame.normal.x>=0,"terminal resolution corridor is not a radial liquid graph");
            outside&=frame.point.x>=finestWidth/maxCurvatureCells;
        }
        if(outside)break;
        cells.push_back(c);
    }
    std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    return cells;
}

std::unique_ptr<InterfaceState> InterfaceNeckCap::prepareTerminalCurvature(const QuadtreeMesh& mesh,
    const NativeLeafDomain& domain,std::unique_ptr<InterfaceState> state,double upper,int maximumLevel,bool includeInteriorNecks)
{
    require(state && state->stage()==InterfaceState::Stage::Ready,"terminal curvature requires a complete candidate");
    state->validateMesh(mesh);domain.validateMeshIdentity(mesh);
    require(domain.axisymmetric() && upper==mesh.maxCoord.y && state->hasConservativePhase(),
        "terminal curvature requires conservative RZ column surfaces");
    if(includeInteriorNecks) {
        auto& data=InterfaceStateAccess::get(*state);data.eventCurvatureCells.clear();
        const bool afterEvent=!state->topology().events().empty();
        InterfaceCurvature resolved(mesh,domain,*state,InterfaceCurvature::Regularization::None);
        std::vector<int> cells;
        for(std::size_t s=0;s<state->topology().surfaces().size();++s) {
            const auto& chain=state->topology().surfaces()[s].geometry.nativeChain();
            require(chain.axisymmetric() && !chain.closed(),"neck models require open RZ axis/section components");
            for(bool end:{false,true}) {
                const auto* section=chain.endpointSymmetryPlane(end);
                require((end?chain.endEndpointKind():chain.beginEndpointKind())==LevelSetInterfacePanels::EndpointKind::Axis
                    || (section && section->normalAxis==1),"neck model has an unsupported endpoint");
            }
            const int count=static_cast<int>(chain.graphs().size());
            std::vector<bool> available(count,true),selected(count,false);
            for(int g=0;g<count;++g) {
                const int c=chain.graphs()[g].nativeCellIndex;
                const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const double h=cell.size().x;
                bool seed=false,haveNormal=false;Vector2D normal;
                for(double f:{0.,.5,1.}) {
                    const auto frame=chain.frameAt(g,f);const auto knot=resolved.evaluatePointwise(s,g,f);
                    seed=seed || frame.point.x<2*h;
                    available[g]=available[g] && frame.point.x>=2*h && knot.valid();
                    if(haveNormal)available[g]=available[g] && normal.x*knot.normal.x+normal.y*knot.normal.y>0;
                    normal=knot.normal;haveNormal=true;
                }
                if(seed) {
                    require(cell.level()<=maximumLevel && (!afterEvent || cell.level()==maximumLevel),
                        "interior neck/cap requires its declared finest level after an event");
                    selected[g]=cell.level()==maximumLevel;
                }
            }
            for(int g=0;g<count;++g)if(selected[g])for(int direction:{-1,1})
                for(int j=g+direction;j>=0 && j<count && !available[j];j+=direction)selected[j]=true;
            // A model tag owns the whole native cell, including subcharts.
            for(int g=0;g<count;++g)if(selected[g]) {
                const int c=chain.graphs()[g].nativeCellIndex;
                require(mesh.cells(mesh.cellIndexToLabel[c]).level()==maximumLevel,
                    "interior model scope exceeds its declared finest region");cells.push_back(c);
            }
        }
        std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
        data.eventCurvatureCells=std::move(cells);
        if(data.conservativeGeometry){data.conservativeGeometry->curvature.clear();
            data.conservativeGeometry->validHeightSamples=data.conservativeGeometry->invalidHeightSamples=0;
            data.conservativeGeometry->validLocalFitSamples=data.conservativeGeometry->unavailableCurvatureSamples=0;}
        return state;
    }
    require(state->topology().surfaces().size()==1,"terminal curvature requires one conservative RZ column surface");
    const auto& chain=state->topology().surfaces()[0].geometry.nativeChain();
    const auto terminal=terminalEnd(*state,upper);
    const int count=static_cast<int>(chain.graphs().size());
    const bool end=terminal.end,detached=terminal.detached;
    auto& data=InterfaceStateAccess::get(*state);
    // This owned candidate may carry the previous model's transported scope.
    // Model availability is not evidence that an ordinary curvature stencil
    // has become resolved. Probe the current geometry/inventory with those
    // tags removed; otherwise the first old model cell above r=2h can remove
    // its own tag and publish an unavailable ordinary force provider.
    // Scalar geometry, C, physical time and event lineage remain unchanged.
    data.eventCurvatureCells.clear();
    InterfaceCurvature resolved(mesh,domain,*state,InterfaceCurvature::Regularization::None);
    std::vector<int> cells;
    for(int g=end?count-1:0;g>=0 && g<count;g+=end?-1:1) {
        const int c=chain.graphs()[g].nativeCellIndex;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const double h=cell.size().x;
        bool resolvedSection=true,haveNormal=false;Vector2D previousNormal;
        for(double f:{0.,.5,1.}) {
            const auto frame=chain.frameAt(g,f);
            // A radial neck is required to construct the section-detachment
            // model. Its transported implicit cap may subsequently turn;
            // signed native curvature remains defined on that regular curve.
            // Topology, finest-level, curvature and force-projection gates
            // still apply. Do not clamp the normal or alter C/geometry here.
            require(detached || frame.normal.x>=0,"terminal neck is not a single radial liquid graph");
            // These are the actual ordinary-side endpoint/midpoint knots.
            // Interpolated availability before tagging would not certify
            // the new boundary knot introduced by this model partition.
            const auto knot=resolved.evaluatePointwise(0,g,f);
            resolvedSection&=frame.point.x>=2*h && knot.valid();
            if(haveNormal)resolvedSection&=previousNormal.x*knot.normal.x+previousNormal.y*knot.normal.y>0;
            previousNormal=knot.normal;haveNormal=true;
        }
        // Tags own whole native cells, including every retained subchart.
        if(resolvedSection && std::find(cells.begin(),cells.end(),c)==cells.end())break;
        // Below the declared finest level the ordinary AMR path must resolve
        // the neck. This policy does not silently waive a refinement request.
        if(cell.level()!=maximumLevel) {
            require(!detached,"detached event region left its declared finest level");
            break;
        }
        cells.push_back(c);
    }
    std::sort(cells.begin(),cells.end());cells.erase(std::unique(cells.begin(),cells.end()),cells.end());
    data.eventCurvatureCells=std::move(cells);
    // Existing pointwise diagnostics were made before this model selection.
    // Do not present them as the actual force provider after tagging.
    if(data.conservativeGeometry){data.conservativeGeometry->curvature.clear();
        data.conservativeGeometry->validHeightSamples=data.conservativeGeometry->invalidHeightSamples=0;
        data.conservativeGeometry->validLocalFitSamples=data.conservativeGeometry->unavailableCurvatureSamples=0;}
    return state;
}

std::unique_ptr<InterfaceState> InterfaceNeckCap::takeUnchargedCandidate(const InterfaceStepTime& time) &&
{return takeCandidate(time);}
std::unique_ptr<InterfaceState> InterfaceNeckCap::takeMaterialCandidate(const InterfaceStepTime& time,const MaterialTransfer& transfer) &&
{
    require(complete_ && transfer.source==sourcePanels_ && transfer.target && transfer.target->usesSource(*complete_->geometry())
        && transfer.cells==cells_ && transfer.sourceSurfaceId==sourceSurfaceId_ && transfer.liquidId==liquidId_
        && transfer.clock.completedSteps==time.index && transfer.clock.endTime==time.end,
        "material candidate requires the matching local event transfer");
    return takeCandidate(time);
}
std::unique_ptr<InterfaceState> InterfaceNeckCap::takeCandidate(const InterfaceStepTime& time)
{
    require(complete_!=nullptr,"candidate was already moved");
    const auto progress=complete_->progress();
    require(progress.completedSteps==time.index && progress.endTime==time.end && time.end>time.begin
        && time.dt>0 && std::isfinite(time.dt) && std::abs(time.end-time.begin-time.dt)
            <=16*std::numeric_limits<double>::epsilon()*std::max({std::abs(time.begin),std::abs(time.end),time.dt}),
        "physical candidate requires its actual completed step interval");
    auto& event=InterfaceStateAccess::get(*complete_).topology->events_.back();
    event.begin=time.begin;event.end=time.end;
    return std::move(complete_);
}

bool InterfaceNeckCap::prepareInterior(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& previous,std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> predictor,
    const InterfaceConservativePhase& inventory,const InterfaceStepTime& time,
    Quantum::Geometry::ActiveHermitePatch& patch,double cutoff,int maximumLevel,InterfaceState& destination,
    std::shared_ptr<const InterfaceInteriorMaterialTransfer>* material)
{
    previous.validateMesh(mesh);domain.validateMeshIdentity(mesh);inventory.validate(mesh,domain);
    require(destination.stage()==InterfaceState::Stage::Empty && predictor && predictor->usesNativeLeaves()
        && inventory.follows(previous.conservativePhase(),time),"interior event requires its actual source/predictor/inventory step");
    require(std::isfinite(cutoff) && cutoff>0 && cutoff<1 && maximumLevel>=0 && maximumLevel<=24,
        "invalid interior-neck cutoff or finest level");
    if(!domain.axisymmetric() || mesh.minCoord.x!=0 || previous.topology().surfaces().size()!=1
        || previous.topology().liquids().size()!=1)return false;
    const auto& sourceChain=previous.topology().surfaces().front().geometry.nativeChain();
    require(!sourceChain.closed() && patch.cornerSectors.empty(),"interior neck requires one open RZ scalar without corner sectors");
    std::vector<int> axial;
    for(int c:domain.cellIndexByRow())if(previous.geometry()->nativeCellBounds(c)[0].x==0)axial.push_back(c);
    std::sort(axial.begin(),axial.end(),[&](int a,int b){
        return previous.geometry()->nativeCellBounds(a)[0].y<previous.geometry()->nativeCellBounds(b)[0].y;});
    const auto radius=[&](const auto& scalar,int c,double z)->std::optional<double> {
        const auto box=scalar->nativeCellBounds(c);const double hr=box[1].x;
        if(!(scalar->scalarJetInNativeCell(c,{0,z}).value<0)
            || !(scalar->scalarJetInNativeCell(c,{hr,z}).value>0))return {};
        LevelSetMappedSegment section(scalar,{0,z},{hr,z});
        if(section.ranges().size()!=2 || !section.ranges().front().negative || section.ranges().back().negative)return {};
        const double r=static_cast<double>(hr*section.ranges().front().length);
        if(!(r>0 && r<hr) || !(scalar->scalarJetInNativeCell(c,{r,z}).dx>0))return {};
        return r;
    };
    std::size_t selected=axial.size();double sourceRadius=0,predictedRadius=0,best=1;
    for(std::size_t k=2;k+2<=axial.size();++k) {
        const auto box=previous.geometry()->nativeCellBounds(axial[k]);const double cut=box[0].y,hr=box[1].x;
        const auto old=radius(previous.geometry(),axial[k],cut),next=radius(predictor,axial[k],cut);
        if(!old || !next || !(*next<*old) || *next>cutoff*hr)continue;
        const auto lo=previous.geometry()->nativeCellBounds(axial[k-1]);
        const auto hi=previous.geometry()->nativeCellBounds(axial[k]);
        if(!(predictor->scalarJetInNativeCell(axial[k-1],{0,lo[0].y}).value<0)
            || !(predictor->scalarJetInNativeCell(axial[k],{0,hi[1].y}).value<0))continue;
        const auto left=radius(predictor,axial[k-1],lo[0].y),right=radius(predictor,axial[k],hi[1].y);
        // This is a declared native SECTION cutoff. Its neighbours must be
        // wider; neither a pole nor a monotonically tapering tip is a neck.
        if((left && *next>=*left) || (right && *next>=*right))continue;
        if(*next/hr<best){selected=k;sourceRadius=*old;predictedRadius=*next;best=*next/hr;}
    }
    if(selected==axial.size())return false;
    const std::array<int,4> cells{{axial[selected-2],axial[selected-1],axial[selected],axial[selected+1]}};
    const auto lower=previous.geometry()->nativeCellBounds(cells[0]);
    const auto upper=previous.geometry()->nativeCellBounds(cells[3]);
    const double cut=previous.geometry()->nativeCellBounds(cells[2])[0].y,hr=lower[1].x;
    const auto size=mesh.cells(mesh.cellIndexToLabel[cells[0]]).size();const double h=size.y;
    const auto eventCell=[&](int c){return std::find(cells.begin(),cells.end(),c)!=cells.end();};
    for(int c:cells) {
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[c]);const auto box=previous.geometry()->nativeCellBounds(c);
        require(cell.level()==maximumLevel && cell.size().x==size.x && cell.size().y==size.y && box[1].x==hr,
            "interior event requires four equal leaves at the declared finest level");
        for(const auto& scalar:{previous.geometry(),predictor}) {
            LevelSetMappedSegment axis(scalar,{0,box[0].y},{0,box[1].y});
            require(!axis.ranges().empty() && std::all_of(axis.ranges().begin(),axis.ranges().end(),
                [](const auto& interval){return interval.negative;}),"interior source/predictor is already disconnected in its model band");
        }
    }
    for(const auto& graph:sourceChain.graphs())if(eventCell(graph.nativeCellIndex))
        require(!graph.coordinateIsX,"interior source neck is not a radial graph throughout the event band");
    auto synchronized=std::make_shared<InterfaceState>(mesh);
    InterfaceGeometryReconstruction::initialize(mesh,patch,true,{},*synchronized,previous.topology().symmetryPlanes());
    const auto base=synchronized->geometry();
    for(int c:domain.cellIndexByRow()) {
        const auto box=base->nativeCellBounds(c);
        if(eventCell(c) || box[1].y<=lower[0].y || box[0].y>=upper[1].y)continue;
        require(synchronized->topology().cellPhase()[c]==1 && inventory.liquidVolumes()[c]
            <=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution*inventory.capacities()[c],
            "interior model band contains resolved liquid beyond its four axis cells");
    }
    using Model=Quantum::Geometry::AxisymmetricNeckCap;
    std::array<Model::Input,2> inputs;std::array<Model::LeftTrace,2> traces;
    std::array<double,2> radii{},scales{};
    const std::array<double,2> seam{{lower[0].y,upper[1].y}};
    for(int side=0;side<2;++side) {
        const int c=side?cells[3]:cells[0];const auto r=radius(base,c,seam[side]);
        require(r.has_value(),"interior cap has no resolved remote radial seam");radii[side]=*r;
        const auto jet=base->scalarJetInNativeCell(c,{*r,seam[side]});
        scales[side]=jet.dx*hr*hr/(2* *r);
        require(std::isfinite(scales[side]) && scales[side]>0,"invalid interior seam scalar scale");
    }
    const double scale=(scales[0]+scales[1])/2,gas=cutoff*cutoff;
    require(std::isfinite(scale) && scale>0,"unrepresentable common interior gas trace scale");
    for(int side=0;side<2;++side) {
        const int a=side?cells[3]:cells[0],b=side?cells[2]:cells[1];const double orientation=side?-1:1;
        const double z=seam[side],r=radii[side];const auto q=base->scalarJetInNativeCell(a,{r,z});
        const auto axis=base->scalarJetInNativeCell(a,{0,z}),outer=base->scalarJetInNativeCell(a,{hr,z});
        require(axis.dx==0 && axis.dxy==0,"interior source seam lacks radial parity");
        inputs[side]={hr,h,r*r,-orientation*2*r*q.dy/q.dx,
            {{inventory.liquidVolumes()[a],inventory.liquidVolumes()[b]}}};
        traces[side]={{axis.value/scale,orientation*axis.dy*h/scale},
            {outer.value/scale,orientation*outer.dy*h/scale},outer.dx*hr/scale,orientation*outer.dxy*hr*h/scale};
    }
    const std::array<Model,2> models{{Model::reconstructMatchedGap(inputs[0],traces[0],gas),
        Model::reconstructMatchedGap(inputs[1],traces[1],gas)}};
    auto completePatch=patch;std::vector<bool> changed(mesh.vertIndexNum(),false);
    for(int i=0;i<4;++i) {
        const int side=i<2?0:1;const double orientation=side?-1:1;
        const auto middle=base->nativeCellBounds(side?cells[2]:cells[1]);
        const double middleZ=side?middle[1].y:middle[0].y;
        for(int label:mesh.cells(mesh.cellIndexToLabel[cells[i]]).cellToVert) {
            const auto& vertex=mesh.verts(label);const auto p=vertex.coord();const int v=vertex.index();
            require(vertex.supCell<0,"interior cap requires independent shared nodes");
            if(p.y==seam[side])continue; // Keep actual remote C1 DOFs exactly.
            require(p.y==middleZ || p.y==cut,"unexpected interior model node");
            const auto q=models[side].scalar(p.x,p.y==cut?2*h:h);
            const Quantum::Geometry::ActiveHermitePatch::NodalDof next{scale*q.value,scale*q.dr,
                orientation*scale*q.dz,orientation*scale*q.drz};
            auto& node=completePatch.nodes[v];
            if(changed[v])require(node.value==next.value && node.dx==next.dx && node.dy==next.dy && node.dxy==next.dxy,
                "opposed caps lost their exact common C1 gas trace");
            node=next;changed[v]=true;
        }
    }
    InterfaceGeometryReconstruction::initialize(mesh,std::move(completePatch),true,{},destination,previous.topology().symmetryPlanes());
    auto& data=InterfaceStateAccess::get(destination);auto& target=*data.topology;
    require(target.liquids_.size()==2 && target.surfaces_.size()==2,"interior cap did not produce exactly two children");
    for(int c:domain.cellIndexByRow())if(!eventCell(c)) {
        require(synchronized->topology().cellPhase()[c]==target.cellPhase_[c],"interior cap changed a remote cell phase");
        bool touched=false;for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)touched=touched || changed[mesh.verts(label).index()];
        require(!touched || (target.cellPhase_[c]==1 && inventory.liquidVolumes()[c]
            <=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution*inventory.capacities()[c]),
            "interior cap changed resolved inventory outside its local cells");
    }
    for(int i=0;i<4;++i) {
        const auto volume=inventory.liquidVolumes()[cells[i]];auto before=nativeVolume(destination.geometry(),cells[i],1,i>=2);
        bool converged=false;
        for(int n=2;n<=64;n*=2) {
            const auto measured=nativeVolume(destination.geometry(),cells[i],n,i>=2);
            if(std::abs(measured-before)<=2048*std::numeric_limits<double>::epsilon()*volume) {
                require(std::abs(measured-volume)<=16384*std::numeric_limits<double>::epsilon()*volume,
                    "interior native cap failed its separate cell inventory");converged=true;break;
            }
            before=measured;
        }
        require(converged,"interior native cap volume quadrature did not converge");
    }
    const auto& old=previous.topology();target.nextLiquidId_=old.nextLiquidId_;target.nextSurfaceId_=old.nextSurfaceId_;
    std::array<int,2> regions{{-1,-1}};
    for(std::size_t s=0;s<target.surfaces_.size();++s) {
        const auto& chain=target.surfaces_[s].geometry.nativeChain();
        const double a=chain.frameAt(0,0).point.y,b=chain.frameAt(static_cast<int>(chain.graphs().size())-1,1).point.y;
        const bool below=std::max(a,b)<cut,above=std::min(a,b)>cut;
        require(below || above,"interior child crosses the declared gas section");const int side=above?1:0;
        require(regions[side]<0,"interior event has two children on one side");
        const auto region=std::find_if(target.liquids_.begin(),target.liquids_.end(),
            [&](const auto& liquid){return liquid.id==target.surfaces_[s].liquidId;});
        require(region!=target.liquids_.end(),"interior surface has no liquid region");
        regions[side]=static_cast<int>(region-target.liquids_.begin());
    }
    require(regions[0]>=0 && regions[1]>=0 && regions[0]!=regions[1],"interior event lost one child");
    std::vector<std::uint64_t> children;std::array<std::uint64_t,2> temporaryIds{};
    for(int side=0;side<2;++side) {
        auto& region=target.liquids_[regions[side]];temporaryIds[side]=region.id;
        region.id=target.nextLiquidId_++;children.push_back(region.id);
    }
    for(auto& surface:target.surfaces_) {
        const auto where=std::find(temporaryIds.begin(),temporaryIds.end(),surface.liquidId);
        require(where!=temporaryIds.end(),"interior child lost its temporary region binding");
        surface.liquidId=children[where-temporaryIds.begin()];surface.id=target.nextSurfaceId_++;
        surface.parents={old.surfaces_[0].id};
    }
    target.events_=old.events_;
    InterfaceTopologyEvent event{InterfaceTopologyEvent::Kind::Split,time.begin,time.end,{old.liquids_[0].id},children,{}};
    event.interiorNecks.push_back({cut,hr,h,cutoff,sourceRadius,predictedRadius,gas,{cells.begin(),cells.end()}});
    target.events_.push_back(std::move(event));data.progress={time.index,time.end};
    // The declared cut partitions the sole parent's material at event birth.
    // It labels every retained C volume, including unresolved remote tails.
    // RZ capacities are linear in axial extent; no nearest-surface assignment.
    data.componentInventories={{children[0],std::vector<long double>(mesh.cellIndexNum())},
        {children[1],std::vector<long double>(mesh.cellIndexNum())}};
    for(int c:domain.cellIndexByRow()) {
        const auto box=destination.geometry()->nativeCellBounds(c);
        const long double below=box[1].y<=cut?1:box[0].y>=cut?0:
            (static_cast<long double>(cut)-box[0].y)/(static_cast<long double>(box[1].y)-box[0].y);
        const auto volume=inventory.liquidVolumes()[c];
        data.componentInventories[0].liquidVolumes[c]=volume*below;
        data.componentInventories[1].liquidVolumes[c]=volume-data.componentInventories[0].liquidVolumes[c];
    }
    data.eventCurvatureCells.assign(cells.begin(),cells.end());
    if(material) {
        // Publish no coupled state here. Retain the synchronized, still
        // connected geometry as the endpoint of the ordinary material leg.
        auto& physical=InterfaceStateAccess::get(*synchronized);
        require(physical.topology->surfaces_.size()==1 && physical.topology->liquids_.size()==1,
            "interior material source must remain connected before closure");
        physical.topology->surfaces_[0].id=old.surfaces_[0].id;
        physical.topology->surfaces_[0].liquidId=old.liquids_[0].id;
        physical.topology->surfaces_[0].parents={old.surfaces_[0].id};
        physical.topology->liquids_[0].id=old.liquids_[0].id;
        physical.topology->nextSurfaceId_=old.nextSurfaceId_;physical.topology->nextLiquidId_=old.nextLiquidId_;
        physical.topology->events_=old.events_;physical.progress={time.index,time.end};
        physical.conservativePhase=inventory;physical.conservativeGeometry.emplace();
        physical.conservativeGeometry->predictor=predictor;
        auto closure=std::make_shared<InterfaceInteriorMaterialTransfer>();
        closure->sourceState=synchronized;closure->time=time;closure->cells=cells;closure->section=cut;
        for(const auto& surface:target.surfaces_)
            closure->targets.push_back({surface.id,surface.liquidId,std::make_shared<LevelSetInterfacePanels>(surface.geometry)});
        const auto& before=physical.topology->surfaces_[0].geometry.panels();closure->weights.resize(before.size());
        std::array<std::vector<std::size_t>,2> sourceLocal,targetLocal;
        std::array<std::size_t,2> targetSide{{0,0}};
        std::vector<std::vector<bool>> covered;
        for(std::size_t s=0;s<closure->targets.size();++s) {
            const auto& panels=closure->targets[s].panels->panels();covered.emplace_back(panels.size(),false);
            const int side=panels.front().frame.point.y>cut?1:0;targetSide[side]=s;
            for(std::size_t j=0;j<panels.size();++j)if(eventCell(panels[j].cellIndex))targetLocal[side].push_back(j);
        }
        for(std::size_t i=0;i<before.size();++i) {
            const auto& a=before[i];const auto box=base->nativeCellBounds(a.cellIndex);
            require(box[1].y<=cut || box[0].y>=cut,"event material panel straddles its native section");
            const int side=box[0].y>=cut?1:0;const auto s=targetSide[side];
            if(eventCell(a.cellIndex)){sourceLocal[side].push_back(i);continue;}
            const auto& panels=closure->targets[s].panels->panels();bool found=false;
            for(std::size_t j=0;j<panels.size();++j) {
                const auto& b=panels[j];
                if(a.cellIndex!=b.cellIndex || a.coordinateIsX!=b.coordinateIsX || a.measure!=b.measure
                    || a.physicalMeasure!=b.physicalMeasure
                    || !Quantum::Geometry::cellChartEqual(a.beginFrame.cellChartPoint,b.beginFrame.cellChartPoint)
                    || !Quantum::Geometry::cellChartEqual(a.endFrame.cellChartPoint,b.endFrame.cellChartPoint))continue;
                require(!found && !covered[s][j],"interior event duplicated a remote material panel");
                found=true;covered[s][j]=true;closure->weights[i].push_back({s,j,1});
            }
            require(found,"interior event changed a remote material panel");
        }
        namespace Exact=Quantum::Geometry::Detail::ExactExpansion;
        const auto cumulative=[](const auto& panels,const auto& indices) {
            std::vector<Exact::Expansion> result(1);
            for(auto i:indices)result.push_back(Exact::add(result.back(),Exact::fromDouble(panels[i].physicalMeasure)));
            return result;
        };
        const auto less=[](const auto& a,const auto& b){return Exact::sign(Exact::subtract(a,b))<0;};
        for(int side=0;side<2;++side) {
            const auto s=targetSide[side];const auto& panels=closure->targets[s].panels->panels();
            require(!sourceLocal[side].empty() && !targetLocal[side].empty(),"interior material half-neck has no area");
            const auto a=cumulative(before,sourceLocal[side]),b=cumulative(panels,targetLocal[side]);
            for(std::size_t i=0;i<sourceLocal[side].size();++i) {
                const auto lo=Exact::multiply(a[i],b.back()),hi=Exact::multiply(a[i+1],b.back());
                const auto width=Exact::subtract(hi,lo);long double sum=0;
                auto& row=closure->weights[sourceLocal[side][i]];
                for(std::size_t j=0;j<targetLocal[side].size();++j) {
                    const auto start=Exact::multiply(b[j],a.back()),end=Exact::multiply(b[j+1],a.back());
                    const auto overlap=Exact::subtract(less(hi,end)?hi:end,less(lo,start)?start:lo);
                    if(Exact::sign(overlap)<=0)continue;
                    const long double share=static_cast<long double>(Exact::estimate(overlap))/Exact::estimate(width);
                    require(std::isfinite(share) && share>0 && share<=1,"unrepresentable interior material overlap");
                    row.push_back({s,targetLocal[side][j],share});sum+=share;covered[s][targetLocal[side][j]]=true;
                }
                require(std::abs(sum-1)<=128*std::numeric_limits<long double>::epsilon()*row.size(),
                    "interior material row does not conserve its integral");
            }
        }
        for(const auto& child:covered)require(std::all_of(child.begin(),child.end(),[](bool v){return v;}),
            "interior material closure lost a child panel");
        *material=std::move(closure);
    }
    return true;
}

InterfaceNeckCap InterfaceNeckCap::prepare(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& source,double upper,double cutoff)
{
    source.validateMesh(mesh);domain.validateMeshIdentity(mesh);
    require(domain.axisymmetric()&&source.hasTopology()&&source.hasConservativePhase(),"requires RZ topology and independent conservative inventory");
    source.conservativePhase().validate(mesh,domain);
    const auto* sourcePatch=std::get_if<Quantum::Geometry::ActiveHermitePatch>(&source.scalarPatch());
    require(sourcePatch&&sourcePatch->cornerSectors.empty(),"matched cap requires the source cubic scalar without corner sectors");
    const auto time=source.progress(),inventoryTime=source.conservativePhase().progress();
    require(time.completedSteps==inventoryTime.completedSteps&&time.endTime==inventoryTime.endTime,"geometry and inventory clocks differ");
    require(std::isfinite(upper)&&upper==mesh.maxCoord.y&&mesh.minCoord.x==0,"only an upper axial reflection section is supported");
    require(std::isfinite(cutoff)&&cutoff>0&&cutoff<1,"cutoff must be an explicit fraction in (0,1)");
    const auto& topology=source.topology();
    require(topology.surfaces().size()==1&&topology.liquids().size()==1,"requires one connected source column");
    const auto& chain=topology.surfaces()[0].geometry.nativeChain();
    require(!chain.closed()&&chain.axisymmetric(),"requires an open axisymmetric source column");
    bool upperEndpoint=false;
    for(bool end:{false,true}) {
        const auto* section=chain.endpointSymmetryPlane(end);
        require(section&&section->normalAxis==1,"source ends must be axial reflection sections");
        upperEndpoint|=section->coordinate==upper;
    }
    require(upperEndpoint,"source does not end at the requested reflection section");
    require(source.geometry()->usesNativeLeaves(),"requires native source geometry");
    std::array<int,2> cells{{-1,-1}};
    for(int c:domain.cellIndexByRow()) {
        const auto box=source.geometry()->nativeCellBounds(c);
        if(box[0].x==0&&box[1].y==upper) {require(cells[1]<0,"ambiguous axis cell");cells[1]=c;}
    }
    require(cells[1]>=0,"no native axis cell at the requested section");
    const auto endBox=source.geometry()->nativeCellBounds(cells[1]);
    const double hr=endBox[1].x,h=endBox[1].y-endBox[0].y;
    for(int c:domain.cellIndexByRow()) {
        const auto box=source.geometry()->nativeCellBounds(c);
        if(box[0].x==0&&box[1].y==endBox[0].y) {require(cells[0]<0,"ambiguous penultimate axis cell");cells[0]=c;}
    }
    require(cells[0]>=0,"no adjacent penultimate axis cell");
    const auto firstBox=source.geometry()->nativeCellBounds(cells[0]);
    const auto& firstCell=mesh.cells(mesh.cellIndexToLabel[cells[0]]);
    const auto& endCell=mesh.cells(mesh.cellIndexToLabel[cells[1]]);
    // Equal dyadic leaves share their nominal widths. Subtracting rounded
    // physical endpoints can differ by an ulp on a non-binary domain (pi/k),
    // even on the uniform root mesh. Keep the exact mesh-level/size condition;
    // the native-volume and C1 seam checks below still use actual coordinates.
    require(firstCell.level()==endCell.level()&&firstCell.size().x==endCell.size().x
        &&firstCell.size().y==endCell.size().y&&firstBox[1].x==hr,
        "event cells must have equal native widths");
    const double z0=firstBox[0].y;
    for(int c:domain.cellIndexByRow()) {
        const auto box=source.geometry()->nativeCellBounds(c);
        if(box[1].y<=z0||box[0].y>=upper||c==cells[0]||c==cells[1])continue;
        require(topology.cellPhase()[c]==1&&source.conservativePhase().liquidVolumes()[c]
            <=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution*source.conservativePhase().capacities()[c],
            "event band contains resolved liquid outside the two axis cells");
    }
    for(const auto& graph:chain.graphs())if(graph.nativeCellIndex==cells[0]||graph.nativeCellIndex==cells[1])
        require(!graph.coordinateIsX,"source event band is not a single radius graph");
    const double neck=sectionRadius(source.geometry(),cells[1],upper);
    require(neck<=cutoff*hr,"source neck has not reached the declared cutoff");
    const double r0=sectionRadius(source.geometry(),cells[0],z0);
    const auto jet=source.geometry()->scalarJetInNativeCell(cells[0],{r0,z0});
    Quantum::Geometry::AxisymmetricNeckCap::Input input;
    input.radialWidth=hr;input.axialWidth=h;input.leftRadiusSquared=r0*r0;
    input.leftRadiusSquaredDerivative=-2*r0*jet.dy/jet.dx;
    for(int i=0;i<2;++i)input.liquidVolumes[i]=source.conservativePhase().liquidVolumes()[cells[i]];
    const double scale=jet.dx*hr*hr/(2*r0);
    require(std::isfinite(scale)&&scale>0,"unrepresentable source scalar normalization");
    const auto axisJet=source.geometry()->scalarJetInNativeCell(cells[0],{0,z0});
    const auto outerJet=source.geometry()->scalarJetInNativeCell(cells[0],{hr,z0});
    require(axisJet.dx==0&&axisJet.dxy==0,"source scalar lacks an even radial extension at the seam");
    const Quantum::Geometry::AxisymmetricNeckCap::LeftTrace left{
        {axisJet.value/scale,axisJet.dy*h/scale},{outerJet.value/scale,outerJet.dy*h/scale},
        outerJet.dx*hr/scale,outerJet.dxy*hr*h/scale};
    InterfaceNeckCap result(source,Quantum::Geometry::AxisymmetricNeckCap::reconstructMatched(input,left));
    result.cells_=cells;result.beginZ_=z0;result.neckRadius_=neck;result.cutoffRadiusCells_=cutoff;
    result.scalarScale_=scale;
    const double nan=std::numeric_limits<double>::quiet_NaN();
    result.patch_.activeCells.assign(mesh.cellIndexNum(),0);
    result.patch_.nodes.assign(mesh.vertIndexNum(),{nan,nan,nan,nan});
    for(int i=0;i<2;++i) {
        result.patch_.activeCells[cells[i]]=1;
        const auto& cell=mesh.cells(mesh.cellIndexToLabel[cells[i]]);
        for(int label:cell.cellToVert) {
            const auto& vertex=mesh.verts(label);const auto p=vertex.coord();
            require(vertex.supCell<0,"event band requires independent seam/cap nodes");
            const double localZ=p.y==z0?0:p.y==endBox[0].y?h:2*h;
            const auto q=result.model_.scalar(p.x,localZ);
            auto& node=result.patch_.nodes[vertex.index()];
            node={scale*q.value,scale*q.dr,scale*q.dz,scale*q.drz};
            if(localZ==0) {
                // Retain actual stored shared DOFs, not a divide/multiply or
                // cache-evaluation round trip. Cells below the seam are exact.
                node=sourcePatch->nodes[vertex.index()];
            }
        }
    }
    result.local_=std::make_shared<const LevelSetHermiteGeometry::Snapshot>(
        LevelSetHermiteGeometry::Snapshot::fromNativeLeaves(mesh,result.patch_,true));
    for(double u:{0.,.25,.5,.75,1.}) {
        const Vector2D p{hr*u,z0};
        const auto old=source.geometry()->scalarJetInNativeCell(cells[0],p);
        const auto next=result.local_->scalarJetInNativeCell(cells[0],p);
        result.seamScalarDefect_=std::max({result.seamScalarDefect_,std::abs(old.value-next.value)/scale,
            hr*std::abs(old.dx-next.dx)/scale,h*std::abs(old.dy-next.dy)/scale,hr*h*std::abs(old.dxy-next.dxy)/scale});
    }
    require(result.seamScalarDefect_<=4096*std::numeric_limits<double>::epsilon()
        *std::max({1.,std::abs(left.axis.value),std::abs(left.outer.value),std::abs(left.axis.derivative),std::abs(left.outer.derivative)}),
        "native scalar does not retain the source C1 seam");
    for(int i=0;i<2;++i) {
        auto previous=nativeVolume(result.local_,cells[i],1);bool converged=false;
        for(int subdivisions=2;subdivisions<=64;subdivisions*=2) {
            const auto current=nativeVolume(result.local_,cells[i],subdivisions);
            if(std::abs(current-previous)<=2048*std::numeric_limits<double>::epsilon()*input.liquidVolumes[i]) {
                result.nativeVolumes_[i]=current;converged=true;break;
            }
            previous=current;
        }
        require(converged,"native bridge volume quadrature did not converge");
        require(std::abs(result.nativeVolumes_[i]-input.liquidVolumes[i])
            <=16384*std::numeric_limits<double>::epsilon()*input.liquidVolumes[i],
            "native scalar sampling fails the per-cell inventory check");
    }
    auto completePatch=*sourcePatch;
    std::vector<std::uint8_t> changed(mesh.vertIndexNum(),0);
    for(int c:cells)for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert) {
        const int v=mesh.verts(label).index();const auto& old=completePatch.nodes[v];const auto& next=result.patch_.nodes[v];
        changed[v]=old.value!=next.value||old.dx!=next.dx||old.dy!=next.dy||old.dxy!=next.dxy;
        completePatch.nodes[v]=next;
    }
    result.complete_=std::make_unique<InterfaceState>(mesh);
    InterfaceGeometryReconstruction::initialize(mesh,std::move(completePatch),true,{},*result.complete_,topology.symmetryPlanes());
    auto& completed=InterfaceStateAccess::get(*result.complete_);auto& target=*completed.topology;
    require(target.liquids_.size()==1&&target.surfaces_.size()==1,"section cap did not retain exactly one domain liquid and surface");
    for(int c:domain.cellIndexByRow())if(c!=cells[0]&&c!=cells[1]) {
        require(topology.cellPhase()[c]==target.cellPhase_[c],"cap changed a cell phase outside its declared event region");
        bool touched=false;for(int label:mesh.cells(mesh.cellIndexToLabel[c]).cellToVert)touched|=changed[mesh.verts(label).index()]!=0;
        require(!touched||(topology.cellPhase()[c]==1&&target.cellPhase_[c]==1
            &&source.conservativePhase().liquidVolumes()[c]
                <=Quantum::Geometry::PlanarCellCut::relativeVolumeResolution*source.conservativePhase().capacities()[c]),
            "cap splice modifies resolved liquid outside its declared event cells");
    }
    const auto& targetChain=target.surfaces_[0].geometry.nativeChain();
    int axisEnds=0,sectionEnds=0;
    for(bool end:{false,true}) {
        const auto kind=end?targetChain.endEndpointKind():targetChain.beginEndpointKind();
        axisEnds+=kind==LevelSetInterfacePanels::EndpointKind::Axis;
        if(const auto* section=targetChain.endpointSymmetryPlane(end)) {
            ++sectionEnds;require(section->coordinate!=upper,"cap retained the detached upper section");
            const auto* a=chain.endpointSymmetryPlane(false);const auto* b=chain.endpointSymmetryPlane(true);
            require((a&&*section==*a)||(b&&*section==*b),"cap changed the surviving source reflection section");
        }
    }
    require(axisEnds==1&&sectionEnds==1,"cap did not replace exactly one reflection endpoint by an axis endpoint");
    target.nextLiquidId_=topology.nextLiquidId_;target.nextSurfaceId_=topology.nextSurfaceId_;
    target.liquids_[0].id=topology.liquids_[0].id;
    target.surfaces_[0].liquidId=target.liquids_[0].id;
    target.surfaces_[0].id=target.nextSurfaceId_++;
    target.surfaces_[0].parents={topology.surfaces_[0].id};
    target.events_=topology.events_;
    target.events_.push_back({InterfaceTopologyEvent::Kind::SectionDetachment,time.endTime,time.endTime,
        {topology.liquids_[0].id},{target.liquids_[0].id},{{1,upper}}});
    completed.progress=time;completed.conservativePhase.emplace(result.inventory_);
    completed.eventCurvatureCells.assign(cells.begin(),cells.end());
    // Rebuild current geometric diagnostics. The pre-event surface transfer
    // and curvature samples cannot be copied onto a different material surface.
    InterfaceConservativeGeometryReport report;report.predictor=source.geometry();
    if(source.hasConservativeGeometry())report.distanceBandWidth=source.conservativeGeometry().distanceBandWidth;
    InterfacePhaseGeometry::MappedFacePartitions partitions(mesh,domain,completed.source);
    auto measured=InterfaceCellMeasure::evaluateVolumes(mesh,domain,partitions,*result.complete_);
    report.geometricVolume=measured.summedNegativeVolume;
    report.geometricVolumeByCell=std::move(measured.negativeVolumeByCell);
    for(int c:domain.cellIndexByRow()) {
        const auto volume=result.inventory_.liquidVolumes()[c],delta=report.geometricVolumeByCell[c]-volume;
        report.inventoryVolume+=volume;report.signedDifference+=delta;report.l1Difference+=std::abs(delta);
        report.maxFractionDifference=std::max(report.maxFractionDifference,std::abs(delta)/result.inventory_.capacities()[c]);
    }
    InterfaceConservativeGeometryReport::ComponentDifference component;component.liquidId=target.liquids_[0].id;
    std::vector<bool> assigned(mesh.cellIndexNum(),false);
    for(int c:target.liquids_[0].cells) {
        assigned[c]=true;const auto volume=result.inventory_.liquidVolumes()[c],geometric=report.geometricVolumeByCell[c];
        component.inventoryVolume+=volume;component.geometricVolume+=geometric;
        component.signedDifference+=geometric-volume;component.l1Difference+=std::abs(geometric-volume);
    }
    for(int c:domain.cellIndexByRow())if(!assigned[c])report.unassignedInventoryVolume+=result.inventory_.liquidVolumes()[c];
    report.components.push_back(component);completed.conservativeGeometry.emplace(std::move(report));
    result.complete_=prepareTerminalCurvature(mesh,domain,std::move(result.complete_),upper,
        mesh.cells(mesh.cellIndexToLabel[cells[1]]).level());
    return result;
}
