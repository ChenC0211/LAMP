#include "Modules/Interface/Compatibility/InterfaceParametricMotion.h"
#include "Modules/Interface/Compatibility/Detail/InterfaceStateInternal.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <map>

namespace {
void require(bool ok,const char* message){if(!ok)throw std::invalid_argument(message);}
InterfaceParametricProjection::Grid makeNativeGrid(const QuadtreeMesh& mesh,const NativeLeafDomain& domain) {
    domain.validateMeshIdentity(mesh);
    require(domain.cellIndexByRow().size()==static_cast<std::size_t>(mesh.cellIndexNum()),"parametric interface requires the complete native domain");

    InterfaceParametricProjection::Grid phaseGrid;phaseGrid.axisymmetric=domain.axisymmetric();
    phaseGrid.left=mesh.minCoord.x;phaseGrid.bottom=mesh.minCoord.y;
    phaseGrid.width=mesh.maxCoord.x-mesh.minCoord.x;phaseGrid.height=mesh.maxCoord.y-mesh.minCoord.y;
    require(std::isfinite(phaseGrid.left)&&std::isfinite(phaseGrid.bottom)&&std::isfinite(mesh.maxCoord.x)&&std::isfinite(mesh.maxCoord.y)
        &&phaseGrid.width>0&&phaseGrid.height>0&&(!domain.axisymmetric()||phaseGrid.left>=0),"invalid physical parametric domain");
    // Uniform and mixed native meshes share this one physical leaf/face
    // binding. No unit-domain numbering or globally finest shadow mesh.
    {
        auto leaves=std::make_shared<InterfaceParametricProjection::Grid::Leaves>(mesh);
        leaves->slotByCell.assign(mesh.cellIndexNum(),-1);long double area=0;
        auto ordered=domain.cellIndexByRow();
        std::sort(ordered.begin(),ordered.end(),[&](int a,int b){
            const auto& ca=mesh.cells(mesh.cellIndexToLabel[a]);const auto& cb=mesh.cells(mesh.cellIndexToLabel[b]);
            const auto pa=mesh.verts(ca.cellToVert[0]).coord(),pb=mesh.verts(cb.cellToVert[0]).coord();
            return pa.x<pb.x||(pa.x==pb.x&&pa.y<pb.y);
        });
        const auto resolution=[](double length,double h){const double ratio=length/h;
            require(std::isfinite(ratio)&&ratio>=1&&ratio<=65536,"invalid parametric leaf scale");
            const int n=int(std::lround(ratio));require(std::abs(ratio-n)<=128*std::numeric_limits<double>::epsilon()*n,"parametric leaf is not a dyadic rectangle");return n;};
        for(int cell:ordered){
            const auto& native=mesh.cells(mesh.cellIndexToLabel[cell]);const auto a=mesh.verts(native.cellToVert[0]).coord(),b=mesh.verts(native.cellToVert[2]).coord();
            const double hx=b.x-a.x,hy=b.y-a.y;
            require(a.x>=mesh.minCoord.x&&a.y>=mesh.minCoord.y&&b.x<=mesh.maxCoord.x&&b.y<=mesh.maxCoord.y&&hx>0&&hy>0,"parametric leaf is outside the physical rectangle");
            phaseGrid.nx=std::max(phaseGrid.nx,resolution(phaseGrid.width,hx));phaseGrid.ny=std::max(phaseGrid.ny,resolution(phaseGrid.height,hy));
            area+=static_cast<long double>(hx)*hy;leaves->slotByCell[cell]=int(leaves->cells.size());
            leaves->cells.push_back({a.x,b.x,a.y,b.y});
        }
        const long double expected=static_cast<long double>(phaseGrid.width)*phaseGrid.height;
        require(std::abs(area-expected)<=128*std::numeric_limits<double>::epsilon()*expected,"parametric leaves do not cover the physical rectangle");
        phaseGrid.leaves=leaves;
        std::array<std::map<double,std::vector<InterfaceParametricProjection::Grid::Face>>,2> lines;
        for(const auto f:domain.faces(mesh)){
            const auto& edge=f.horizontal?mesh.hEdges(f.edgeLabel):mesh.vEdges(f.edgeLabel);
            const auto a=mesh.verts(edge.edgeToVert[0]).coord(),b=mesh.verts(edge.edgeToVert[1]).coord();
            const double lower=f.horizontal?a.x:a.y,upper=f.horizontal?b.x:b.y;
            require(lower<upper,"native parametric face orientation is reversed");
            lines[f.horizontal][f.horizontal?f.center.y:f.center.x].push_back({lower,upper,
                f.negativeCell<0?-1:leaves->slotByCell.at(f.negativeCell),f.positiveCell<0?-1:leaves->slotByCell.at(f.positiveCell),f.nativeRow});
        }
        for(int axis=0;axis<2;axis++){
            int offset=0;
            for(auto& item:lines[axis]){
                auto& faces=item.second;std::sort(faces.begin(),faces.end(),[](const auto& a,const auto& b){return a.lower<b.lower;});
                double capacity=std::numeric_limits<double>::infinity(),previous=-std::numeric_limits<double>::infinity();
                for(const auto& f:faces){require(f.lower>=previous,"native parametric face segments overlap");previous=f.upper;
                    for(int c:{f.negative,f.positive})if(c>=0)capacity=std::min(capacity,phaseGrid.capacity(c));
                    ++offset;}
                require(capacity>0&&std::isfinite(capacity),"native face has no inventory neighbor");
                leaves->planes[axis].push_back({item.first,capacity,offset-int(faces.size()),std::move(faces)});
            }
            const double lo=axis?mesh.minCoord.y:mesh.minCoord.x,hi=axis?mesh.maxCoord.y:mesh.maxCoord.x;
            require(!leaves->planes[axis].empty()&&leaves->planes[axis].front().coordinate==lo&&leaves->planes[axis].back().coordinate==hi,"parametric native planes do not close the physical domain");
        }
    }
    return phaseGrid;
}
InterfaceParametricTransport::NativeFlux bindNativeFlux(const InterfaceParametricProjection::Grid& grid,const InterfaceNormalFluxView& spatial) {
    InterfaceParametricTransport::NativeFlux flux;flux.grid=grid;flux.retainedMomentDegree=spatial.retainedMomentDegree();
    flux.retainedVelocityProducts=spatial.retainedVelocityProducts();
    const auto layout=flux;
    // Check every original half chart before the synchronous calculation.
    for(bool horizontal:{false,true})for(int line=0;line<layout.planeCount(horizontal);++line)for(int strip=0;strip<layout.stripCount(horizontal,line);++strip) {
        const int row=layout.face(horizontal,line,strip).nativeRow;const auto chart=spatial.faceChart(row);
        const auto face=layout.face(horizontal,line,strip);const double coordinate=layout.planeCoordinate(horizontal,line);
        const auto normal=[&](Vector2D p){return horizontal?p.y:p.x;};const auto tangent=[&](Vector2D p){return horizontal?p.x:p.y;};
        require(normal(chart.begin)==coordinate && normal(chart.midpoint)==coordinate && normal(chart.end)==coordinate &&
            tangent(chart.begin)==face.lower && tangent(chart.end)==face.upper && tangent(chart.midpoint)==.5*(face.lower+face.upper),
            "parametric flux primitive is not on the bound native half charts");
    }
    flux.integral=[&spatial,layout](bool horizontal,int line,int strip,double a,double b) {
        const int row=layout.face(horizontal,line,strip).nativeRow;const auto chart=spatial.faceChart(row);
        const Vector2D points[3]={chart.begin,chart.midpoint,chart.end};long double out=0;
        for(int half=0;half<2;++half){const double lo=horizontal?points[half].x:points[half].y,hi=horizontal?points[half+1].x:points[half+1].y;
            out+=spatial.integratedNormalFlux(row,half,std::clamp((a-lo)/(hi-lo),0.,1.),std::clamp((b-lo)/(hi-lo),0.,1.));}
        return static_cast<double>(out);
    };
    flux.density=[&spatial,layout](bool horizontal,int line,int strip,double coordinate) {
        const int row=layout.face(horizontal,line,strip).nativeRow;const auto chart=spatial.faceChart(row);
        const double middle=horizontal?chart.midpoint.x:chart.midpoint.y;const int half=coordinate<=middle?0:1;
        const double a=half?middle:(horizontal?chart.begin.x:chart.begin.y),b=half?(horizontal?chart.end.x:chart.end.y):middle;
        return static_cast<double>(spatial.normalFluxDensity(row,half,false,std::clamp((coordinate-a)/(b-a),0.,1.))/(b-a));
    };
    auto bounds=std::make_shared<std::array<std::vector<double>,2>>();
    for(bool horizontal:{false,true}){auto& values=(*bounds)[horizontal];values.resize(layout.faceCount(horizontal));
        for(int line=0;line<layout.planeCount(horizontal);++line)for(int strip=0;strip<layout.stripCount(horizontal,line);++strip){
            const int slot=layout.faceOffset(horizontal,line)+strip,row=layout.face(horizontal,line,strip).nativeRow;const auto face=layout.face(horizontal,line,strip);
            values[slot]=static_cast<double>(2/(face.upper-face.lower)*std::max(spatial.densityMagnitudeBound(row,0),spatial.densityMagnitudeBound(row,1)));
            if(values[slot]>0)values[slot]=std::nextafter(values[slot],std::numeric_limits<double>::infinity());}}
    flux.densityMagnitudeBound=[bounds,layout](bool horizontal,int line,int strip){return (*bounds)[horizontal].at(layout.faceOffset(horizontal,line)+strip);};
    if(spatial.hasNormalFluxPolynomial())flux.polynomial=[&spatial,layout](bool horizontal,int line,int strip,int half) {
        return spatial.normalFluxPolynomial(layout.face(horizontal,line,strip).nativeRow,half);
    };
    return flux;
}
}

std::vector<InterfaceParametricTransport::FaceVolumes::Face> InterfaceParametricMotion::backgroundFaceVolumes(
    const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceNormalFluxView& spatial,double duration){
    require(std::isfinite(duration)&&duration>0,"background face flux needs a positive interval");
    const auto grid=makeNativeGrid(mesh,domain);const auto flux=bindNativeFlux(grid,spatial);
    InterfaceParametricTransport::FluxReport report;InterfaceParametricTransport::FaceVolumes volumes;
    InterfaceParametricTransport::backgroundFlux({flux,flux},duration,report,volumes);
    std::vector<InterfaceParametricTransport::FaceVolumes::Face> result(domain.nativeFaceCount());
    for(int axis=0;axis<2;axis++)for(const auto& plane:grid.leaves->planes[axis])
        for(int i=0;i<int(plane.faces.size());i++)result.at(plane.faces[i].nativeRow)=volumes.faces[axis].at(plane.offset+i);
    return result;
}

InterfaceState InterfaceParametricMotion::initialize(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    InterfaceParametricProjection::Curve curve) {
    return bind(mesh,domain,[&](const auto& grid){return InterfaceParametricProjection::initialize(std::move(curve),grid);},{});
}
InterfaceState InterfaceParametricMotion::initialize(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    std::vector<InterfaceParametricProjection::Curve> curves) {
    return bind(mesh,domain,[&](const auto& grid){return InterfaceParametricProjection::initialize(std::move(curves),grid);},{});
}
InterfaceParametricMotion::MeshTransfer InterfaceParametricMotion::prepareMeshTransfer(const QuadtreeMesh& sourceMesh,
    const NativeLeafDomain& sourceDomain,const InterfaceState& source,const QuadtreeMesh& targetMesh,
    const NativeLeafDomain& targetDomain,const InterfaceParametricProjection::Options& options){
    validateDomain(sourceMesh,sourceDomain,source);
    require(sourceDomain.axisymmetric()==targetDomain.axisymmetric(),"parametric mesh transfer changed physical measure");
    MeshTransfer result;
    result.state=bind(targetMesh,targetDomain,[&](const auto& grid){
        return InterfaceParametricProjection::transferMesh(source.parametricGeometry(),grid,options,result.report);
    },source.progress());result.state.parametric_->neckEvents=source.parametric_->neckEvents;return result;
}
InterfaceState InterfaceParametricMotion::prepareNeckSplit(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& source,const InterfaceParametricProjection::NeckSplitOptions& options){
    auto assessment=assessNeckSplit(mesh,domain,source,options);require(bool(assessment.state),"neck event is not geometrically ready");return std::move(*assessment.state);
}
InterfaceParametricMotion::NeckAssessment InterfaceParametricMotion::assessNeckSplit(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& source,const InterfaceParametricProjection::NeckSplitOptions& options){
    validateDomain(mesh,domain,source);
    auto assessment=InterfaceParametricProjection::assessNeckSplit(source.parametricGeometry(),options);
    NeckAssessment result{assessment.readiness,assessment.report,{}};if(!assessment.split)return result;
    auto split=std::move(*assessment.split);
    InterfaceState next;const auto clock=source.progress();
    next.parametric_=std::make_unique<InterfaceState::ParametricData>(source.parametric_->binding,std::move(split.candidate),clock);
    next.parametric_->neckEvents=source.parametric_->neckEvents;
    next.parametric_->neckEvents.push_back({std::move(split.report),clock.completedSteps,clock.endTime});
    result.state=std::move(next);return result;
}
InterfaceState InterfaceParametricMotion::bind(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const std::function<InterfaceParametricProjection::Candidate(const InterfaceParametricProjection::Grid&)>& make,
    InterfaceStepProgress progress) {
    const auto phaseGrid=makeNativeGrid(mesh,domain);
    auto binding=std::make_shared<InterfaceState::ParametricData::Binding>();binding->mesh=&mesh;binding->domain=&domain;binding->identity=domain.meshIdentity();
    binding->faceSlots.resize(domain.nativeFaceCount());binding->cellIndices.resize(phaseGrid.size());
    for(int cell=0;cell<mesh.cellIndexNum();cell++)binding->cellIndices.at(phaseGrid.leaves->slotByCell.at(cell))=cell;
    for(int axis=0;axis<2;axis++){
        auto& rows=axis?binding->horizontalRows:binding->verticalRows;
        for(const auto& plane:phaseGrid.leaves->planes[axis])for(const auto& face:plane.faces){
            binding->faceSlots.at(face.nativeRow)={axis,int(rows.size())};rows.push_back(face.nativeRow);
        }
    }
    auto candidate=make(phaseGrid);InterfaceState out;
    out.parametric_=std::make_unique<InterfaceState::ParametricData>(std::move(binding),std::move(candidate),progress);return out;
}

const std::vector<int>& InterfaceParametricMotion::cellIndices(const InterfaceState& state) {
    require(state.hasParametricGeometry(),"parametric cell indices require the parametric state representation");return state.parametric_->binding->cellIndices;
}
const InterfaceParametricTransport::FaceVolumes::Face& InterfaceParametricMotion::faceVolumes(const InterfaceState& state,int row){
    const auto& flux=state.parametricFluxes();const auto slot=state.parametric_->binding->faceSlots.at(row);
    return flux.faces[slot.first].at(slot.second);
}
void InterfaceParametricMotion::validateDomain(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,const InterfaceState& state) {
    require(state.hasParametricGeometry(),"parametric consumer requires a parametric InterfaceState");
    state.validateMesh(mesh);domain.validateMeshIdentity(mesh);
    require(state.parametric_->binding->domain==&domain && state.parametric_->binding->identity==domain.meshIdentity(),
        "parametric consumer domain ownership changed");
    require(state.parametric_->candidate.grid().axisymmetric==domain.axisymmetric(),"parametric inventory/domain measure mismatch");
}

void InterfaceParametricMotion::refineMaterialSpace(Result& prepared,const std::vector<double>& minimumKnots) {
    require(prepared.state.hasParametricGeometry()
        &&prepared.state.parametricGeometry().boundaryHandle()==prepared.projected_
        &&prepared.sources_.size()==prepared.state.parametricGeometry().componentCount(),"material refinement lost its prepared motion identity");
    auto options=prepared.projection_;auto& knots=options.minimumMaterialKnots;
    for(double s:minimumKnots){require(std::isfinite(s)&&s>=0&&s<=1,"invalid retained material knot");knots.push_back(s);}
    for(const auto& c:prepared.projected_->components())for(const auto& arc:c.curve->arcs)knots.push_back(arc.a);
    knots.push_back(1);std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
    auto& state=*prepared.state.parametric_;
    const auto local=InterfaceParametricProjection::activeOptions(state.candidate,options);
    if(prepared.sources_.empty()){prepared.projection_=std::move(options);return;}
    std::vector<InterfaceParametricProjection::Candidate> components;std::vector<std::uint64_t> ids;
    for(std::size_t i=0;i<prepared.sources_.size();i++){
        const auto part=state.candidate.component(i);
        components.push_back(InterfaceParametricProjection::prepare(prepared.sources_[i],state.candidate.grid(),part.inventory(),local));ids.push_back(prepared.projected_->components()[i].id);
    }
    auto candidate=InterfaceParametricProjection::combine(state.candidate.grid(),std::move(components),ids,state.candidate.retiredComponents());
    require(candidate.report().maximumInventoryFractionError<=options.inventoryFractionUncertainty+options.measureFractionReserve,
        "material refinement and retired inventory exceed the common fraction budget");
    require(std::abs(candidate.report().totalGeometryInventoryDifference)<=options.measureFractionReserve*state.candidate.grid().totalCapacity(),
        "material refinement and retired inventory exceed the common total-volume reserve");
    // Only the successfully reconstructed candidate is replaced. The native
    // phase flux, clock, binding and old-geometry lineage remain identical.
    state.candidate=std::move(candidate);prepared.projected_=state.candidate.boundaryHandle();prepared.projection_=std::move(options);
}

void InterfaceParametricMotion::retireExteriorComponents(Result& prepared){
    require(prepared.state.hasParametricGeometry()&&prepared.state.parametricGeometry().boundaryHandle()==prepared.projected_
        &&prepared.sources_.size()==prepared.projected_->components().size(),"component retirement lost its prepared motion identity");
    auto& state=*prepared.state.parametric_;
    auto candidate=InterfaceParametricProjection::retireExterior(state.candidate,state.progress.endTime,prepared.projection_);
    if(candidate.boundaryHandle()==prepared.projected_)return;
    std::vector<InterfaceParametricProjection::Source> sources;std::size_t old=0;
    for(const auto& c:candidate.boundaryHandle()->components()){
        while(old<prepared.projected_->components().size()&&prepared.projected_->components()[old].id!=c.id)++old;
        require(old<prepared.sources_.size(),"component retirement changed active identity order");sources.push_back(prepared.sources_[old]);
    }
    state.candidate=std::move(candidate);prepared.sources_=std::move(sources);prepared.projected_=state.candidate.boundaryHandle();
}

InterfaceParametricMotion::Result InterfaceParametricMotion::prepare(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& old,const InterfaceStepTime& time,double duration,const InterfaceParametricTransport::Flow& flow,
    const InterfaceNormalFluxView& spatial,const InterfaceParametricProjection::Options& options) {
    return prepareInterval(mesh,domain,old,time,duration,flow,spatial,spatial,options);
}
InterfaceParametricMotion::Result InterfaceParametricMotion::prepareTimeLinear(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& old,const InterfaceStepTime& time,const InterfaceParametricTransport::Flow& flow,
    const InterfaceNormalFluxView& beginFlux,const InterfaceNormalFluxView& endFlux,const InterfaceParametricProjection::Options& options) {
    return prepareInterval(mesh,domain,old,time,time.dt,flow,beginFlux,endFlux,options);
}
InterfaceParametricMotion::Result InterfaceParametricMotion::prepareInterval(const QuadtreeMesh& mesh,const NativeLeafDomain& domain,
    const InterfaceState& old,const InterfaceStepTime& time,double duration,const InterfaceParametricTransport::Flow& flow,
    const InterfaceNormalFluxView& beginFlux,const InterfaceNormalFluxView& endFlux,const InterfaceParametricProjection::Options& options) {
    validateDomain(mesh,domain,old);
    const auto binding=old.parametric_->binding;
    const auto progress=old.progress();
    require(time.index==progress.completedSteps+1 && time.begin==progress.endTime && time.end>time.begin && time.dt>0 &&
        std::isfinite(time.begin) && std::isfinite(time.end) && std::isfinite(time.dt) &&
        std::abs((time.end-time.begin)-time.dt)<=32*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(time.begin),std::abs(time.end)}),
        "parametric motion time does not follow the accepted InterfaceState");
    const auto& grid=old.parametricGeometry().grid();
    const InterfaceParametricTransport::TimeLinearFlux flux{bindNativeFlux(grid,beginFlux),bindNativeFlux(grid,endFlux)};
    Result out;
    auto step=InterfaceParametricTransport::prepare(old.parametric_->candidate,duration,flow,flux,options,
        [&](const InterfaceParametricProjection::Source& source,const InterfaceParametricProjection::Candidate&){
            out.sources_.push_back(source);
        });
    InterfaceState next;next.parametric_=std::make_unique<InterfaceState::ParametricData>(binding,std::move(step.candidate),InterfaceStepProgress{time.index,time.end});
    next.parametric_->stepTime=time;next.parametric_->previousGeometry=old.parametric_->candidate.boundaryHandle();
    next.parametric_->fluxes=std::move(step.faceVolumes);
    next.parametric_->neckEvents=old.parametric_->neckEvents;
    out.projected_=next.parametric_->candidate.boundaryHandle();out.projection_=options;
    out.state=std::move(next);out.transport=step.transport;out.flux=step.flux;out.inventoryBalance=step.inventoryBalance;return out;
}
