#include "IO/HydroOutputComponents.h"
#include "Modules/Hydrodynamics/HydroVofStep.h"
#include "Modules/Interface/GeometricVofMeshTransfer.h"
#include "Core/Geometry/PlanarCellCut.h"
#include "Core/IO/OutputProvenance.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace AppIO::HydroOutput {
using namespace LAMP::Hydrodynamics;
using LAMP::IO::stream;
namespace {
const LAMP::Interface::GeometricVofPublication& publication(const HydroTwoPhaseState& state){
    if(state.model()!=HydroPhaseGeometry::Model::GeometricVof)throw std::invalid_argument("VOF output requires a genuine VOF state");
    const auto& p=state.phaseGeometry().vofPublication();
    if(p.geometryKind()!=LAMP::Interface::GeometricVofGeometryKind::LocalPlic)throw std::invalid_argument("VOF output requires LocalPlic");
    return p;
}
bool sameClock(LAMP::TimeLayer a,LAMP::TimeLayer b){
    return a.step==b.step&&a.begin==b.begin&&a.end==b.end&&a.position==b.position;
}
void validateStep(const HydroTwoPhaseState& state,const HydroVofStep* step,
    const LAMP::Interface::GeometricVofMeshTransfer* transfer){
    if(!transfer){
        if(step&&!step->candidateState().sameIdentity(state))throw std::invalid_argument("VOF output step differs from accepted state");
        return;
    }
    if(!step)throw std::invalid_argument("VOF remapped output requires the original physical step");
    const auto& source=publication(step->candidateState());const auto& target=publication(state);
    const auto& candidate=transfer->candidate();const auto clock=state.timeLayer();
    if(!transfer->report().localPlic||candidate.geometryKind()!=LAMP::Interface::GeometricVofGeometryKind::LocalPlic||
        transfer->sourceIdentity()!=source.identity()||candidate.identity()!=target.identity()||
        clock.position!=LAMP::TimeLayer::Position::End||!sameClock(clock,source.timeLayer())||
        !sameClock(clock,candidate.timeLayer())||!sameClock(clock,step->candidateState().timeLayer()))
        throw std::invalid_argument("VOF remapped output requires authenticated source/target identities and same complete End clock");
}
template<class T> std::string scalarText(T value){auto out=stream();out<<value;return out.str();}
struct Sum {long double value=0,tail=0;void add(long double x){const auto a=x-tail,t=value+a;tail=(t-value)-a;value=t;}};
struct ForceRange {bool defined=false;double lower=0,upper=0;};
ForceRange forceRange(const HydroVofStep* step){
    ForceRange out;if(!step||!step->hasCapillaryLoad())return out;
    const auto& force=step->capillaryLoad();const auto& k=force.heightFunctionFaceCurvatures();const auto& support=force.supportReports();
    if(k.size()!=support.size())throw std::logic_error("VOF HF output support dimensions differ");
    for(std::size_t u=0;u<k.size();++u)if(support[u].curvatureDefined&&std::isfinite(k[u])){
        if(!out.defined){out.lower=out.upper=k[u];out.defined=true;}else{out.lower=std::min(out.lower,k[u]);out.upper=std::max(out.upper,k[u]);}}
    return out;
}
void appendVofMomentumBudget(std::ostream& out,const HydroVofStep* step){
    constexpr std::size_t count=29;
    const auto* report=step?&step->report():nullptr;
    if(!report||report->stages.empty()){for(std::size_t i=0;i<count;++i)out<<",NA";return;}
    // Predictor and corrector are alternative states, not two time advances.
    const auto& stage=report->stages.back();const auto& energy=stage.energy;
    const auto scalar=[&](long double value){out<<',';if(std::isfinite(value))out<<value;else out<<"NA";};
    const std::array<long double,18> work{{
        stage.kineticChange,
        stage.capillaryWork,
        energy.bodyWork,
        energy.viscousWork,
        energy.pressureWork,
        energy.halfRWork,
        energy.mappedPrimalHalfRWork,
        energy.geometryTransferHalfRWork,
        energy.spatialExchangeHalfRWork,
        energy.timeMassWork,
        energy.roundedPhaseMassWork,
        energy.operatorRateRoundoffWork,
        energy.rawMinusCanonicalOutflowHalfWork,
        energy.halfRDecompositionResidual,
        energy.momentumEquationWork,
        energy.kineticBudgetResidual,
        energy.closureResidual,
        energy.closureRelativeResidual
    }};
    for(const auto value:work)scalar(value);
    out<<','<<energy.stableGasDofs;
    if(energy.stableGasDofs>0){scalar(energy.stableGasHalfRWork);scalar(energy.stableGasSpatialExchangeHalfRWork);scalar(energy.stableGasMaxAbsIntegratedDefectOverMidMass);}
    else out<<",NA,NA,NA";
    const bool peak=energy.stableGasDofs>0&&energy.peakGasNativeFace>=0;
    out<<',';if(peak)out<<energy.peakGasNativeFace;else out<<"NA";
    const std::array<long double,6> peakValues{{
        energy.peakGasAbsEndNormal,
        energy.peakGasIntegratedDefectOverMidMass,
        energy.peakGasMassChangeOverMidMass,
        energy.peakGasOutwardIntegralOverMidMass,
        energy.peakGasGeometryTransferOverMidMass,
        energy.peakGasSpatialExchangeOverMidMass
    }};
    for(const auto value:peakValues){if(peak)scalar(value);else out<<",NA";}
}
}
Observation observeVof(const HydroTwoPhaseState& state,const HydroVofStep* step,const LAMP::Interface::GeometricVofMeshTransfer* transfer){
    publication(state);validateStep(state,step,transfer);Observation out;out.scalarRepresentation="localPlic";out.gradientUpdate="conservativeAlpha";
    out.capillaryDiscretization="localPlicHeightFunction";out.heightFunctionForce=step&&step->hasCapillaryLoad();const auto r=forceRange(step);
    out.faceCurvatureDefined=r.defined;out.faceCurvatureMin=r.lower;out.faceCurvatureMax=r.upper;return out;
}
std::string vofBuildInfo(const HydrodynamicsInput& input,const Observation&){
    const bool projected=input.vof.momentumDensityRepresentation==HydroMomentumMass::DensityRepresentation::CellAverageProjection;
    return "compiler="+LAMP::IO::compiler()+"\nbuild_type="+LAMP::IO::buildType()+"\ncompiled="+__DATE__+" "+__TIME__+
        "\ndensity_representation="+std::string(projected?"cellAverageProjection":"sharpGeometry")+
        "\nrepresented_dual_phase_volumes="+std::string(projected?"cellwise uniform P0 from independently retained liquid/gas inventories":"sharp actual LocalPlic subdual phase geometry")+
        "\nunits=code\ninterface_method=vof\nprimary=alpha (independent liquid/gas stocks)\nscalar_representation=localPlic\n"
        "volume_output=tecplotDat\nframe_index=output/meta/frames.csv\n"
        "transport=conservative local split PLIC; retained accepted P0 face flux and boundary-preserving volume-coordinate Coons material extension\n"
        "capillary_model=pressure-compatible alpha height function\ncurvature=actual midpoint HF force; unavailable initially; not PLIC line curvature\n"
        "height_support=coarse fallback uses actual PLIC segment midpoint and minimum native leaf spacing\n"
        "height_fallback=native9/endpoint/refined HF first; fixed 5h same-owner segment-midpoint native-face cluster, segment-length-weighted quadratic single QR only after all HF fail; no curvature clip or physical-parameter change\n"
        "height_function_support=original preferred-axis 3x9 slices then alternate-axis 3x9; only if both fail, same-owner batched missing-pure-end extension shared by all three columns, at most four extra rings per end; original face curvature averaging retained\n"
        "height_function_extensions=height_extended_cells counts cells whose final curvature uses extended support; height_extended_slices counts all additionally measured virtual slices; zero-sigma behavior unchanged\n"
        "pressure=actual accepted midpoint native saddle pressure\nvelocity=accepted native P0-side RT0/Piola cell-center diagnostic; not a new primary velocity\nvof_native_velocity=first four groups retain exact adjacent inventory classes; appended geometry_pure_gas uses accepted End actual measured PLIC phase; all normal DOFs are actual, force columns retain midpoint clock\n"
        "interface=independent actual two-endpoint local PLIC segments; no global connection or smoothing\n"
        "liquid_p0_moments=raw extensive first/second native-cell inventory moments at accepted clock; uniform alpha in each native rectangle, not a PLIC subcell moment\n"
        "liquid_p0_moment_coordinates=XY x/y or RZ r/z; RZ integrates actual 2*pi*r measure; nonfinite derived diagnostic is NA\n"
        "vof_momentum_work=step-integrated actual final AC stage only; initial/nonfinite diagnostics NA; predictor and corrector are not summed; original energy-series meaning unchanged\n"
        "vof_kinetic_budget=deltaK=capillary+body-viscous-pressure+halfR+timeMass+roundedPhaseMass+operatorRateRoundoff+kineticBudgetResidual; no surface or total-energy identity\n"
        "vof_half_R_budget=halfR=mappedPrimalHalfR+geometryTransferHalfR+spatialExchangeHalfR+halfRDecompositionResidual\n"
        "vof_work_closure=closureResidual=kineticBudgetResidual-momentumEquationWork; momentumEquationWork=dt*actual momentum power residual; closureRelativeResidual uses gross absolute work; storage/midpoint/matrix evaluation roundoff may remain, with no new gate\n"
        "vof_operator_rate_roundoff=sum((Q/2-r-dt*hii)*ubar^2); Q=canonical integrated dual outward flux, r=integrated defect, hii=actual frozen H diagonal; not the full velocity matrix\n"
        "vof_raw_outflow_half_work=rawMinusCanonicalOutflowHalfWork=sum((dt*rawOutwardRate-Q)*ubar^2/2), a side diagnostic not added to the kinetic budget\n"
        "vof_stable_gas=begin/mid/end represented dual liquid volumes all exactly zero for the selected density representation; Sharp uses actual PLIC phase geometry, P0 uses independent inventories projected uniformly; positive micro tails may remove this exact-zero group; count0 leaves gas work/ratios/peak unavailable\n"
        "vof_geometry_end_gas=adjacent accepted End PLIC measured liquid volumes all exactly zero; independent inventory tails do not alter this End-only classification; it is not a three-clock stable-gas or no-interface-throughout-step claim\n"
        "vof_peak_gas=maximum abs accepted-end normal velocity among stable-gas duals; actual native face; all R/mass-change/outflow/geometry/dt*spatial-exchange ratios refer to the same final-stage midpoint mass; R is integrated mass defect\n"
        "dynamic_amr_output=new accepted stock/moments/divergence/velocity/retained-M kinetic use target mesh; physical work/dual-R/HF scalar rows retain original physical-step owner\n"
        "dynamic_amr_force=after actual remap, new native force/kappa columns are unavailable; field placeholders are zero+Available0, original physical force clock/identity retained; no End HF recomputation\n"
        "dynamic_amr_energy=accepted kinetic includes remap change; physical work budget uses original step and excludes remap; regrid.dat native deltaP/deltaKE are separate remap diagnostics, never external work\n"
        "dynamic_amr_pressure=physical midpoint sample volume-averaged and regauged onto target; remap multiplier is not physical pressure\n"
        "dynamic_amr_native_face_ids=vof_work peak face/R diagnostics retain physical-step source layout; vof_native_velocity peaks use current target layout\n"
        "capillary_curvature=HF fields retain physical HF/fit estimates; effective fields retain actual force closure; method3 is unresolved isolated zero closure, not physical curvature or a resolved droplet\n"
        "surface_area=not_computed\nenergy_area=not_available\nenergy=actual kinetic and native work diagnostics only; surface/total/closure unavailable\n";
}
std::string vofNativeVelocityCsvHeader(){
    auto out=stream();out<<"step,accepted_time,geometry_identity,force_time,force_geometry_identity";
    for(const char* group:{"global","pure_liquid","pure_gas","interface","geometry_pure_gas"})for(const char* field:{"abs_normal_velocity","normal_velocity","velocity_dof","native_face","normal_axis","face_x","face_y","negative_cell","positive_cell","negative_alpha","positive_alpha","midpoint_capillary_load","midpoint_force_kappa"})out<<','<<group<<'_'<<field;
    out<<'\n';return out.str();
}
std::string vofNativeVelocityRow(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroVofStep* step,const LAMP::Interface::GeometricVofMeshTransfer* transfer){
    state.validateFor(mesh);const auto& pub=publication(state);validateStep(state,step,transfer);const auto& layout=state.values().layout();const auto& velocity=state.values().normalVelocity();
    std::array<int,5> peak{{-1,-1,-1,-1,-1}};
    const auto keep=[&](int group,int u){if(peak[group]<0||std::abs(velocity[u])>std::abs(velocity[peak[group]]))peak[group]=u;};
    for(int u=0;u<layout.velocityCount();++u){const auto f=layout.domain().face(mesh,layout.nativeFace(u));bool liquid=true,gas=true,geometryGas=true;
        for(int native:{f.negativeCell,f.positiveCell})if(native>=0){const auto& cell=pub.plicCells().at(native);const auto& inventory=cell.inventory;
            liquid=liquid&&inventory[1]==0;gas=gas&&inventory[0]==0;geometryGas=geometryGas&&cell.measured[0]==0;}
        keep(0,u);keep(liquid?1:gas?2:3,u);if(geometryGas)keep(4,u);
    }
    const auto* force=step&&step->hasCapillaryLoad()?&step->capillaryLoad():nullptr;
    auto out=stream();out<<state.timeLayer().step<<','<<state.timeLayer().time()<<','<<pub.identity().value()<<',';
    if(force)out<<force->phaseGeometry().timeLayer().time()<<','<<force->phaseGeometry().vofIdentity().value();else out<<"NA,NA";
    for(int u:peak){if(u<0){for(int field=0;field<13;++field)out<<",NA";continue;}
        const auto f=layout.domain().face(mesh,layout.nativeFace(u));out<<','<<std::abs(velocity[u])<<','<<velocity[u]<<','<<u<<','<<f.nativeRow<<','<<(layout.normalAxis(u)==HydroDofLayout::NormalAxis::First?0:1)<<','<<f.center.x<<','<<f.center.y<<','<<f.negativeCell<<','<<f.positiveCell;
        for(int native:{f.negativeCell,f.positiveCell}){out<<',';if(native>=0)out<<pub.alpha().at(native);else out<<"NA";}
        out<<',';if(force&&!transfer)out<<force->momentumRhsLoad().at(u);else out<<"NA";
        out<<',';if(force&&!transfer&&force->supportReports().at(u).curvatureDefined&&std::isfinite(force->heightFunctionFaceCurvatures().at(u)))out<<force->heightFunctionFaceCurvatures()[u];else out<<"NA";
    }
    out<<'\n';return out.str();
}
std::string vofManifestFields(){
    return ",\"primary\":\"alpha\",\"scalar_representation\":\"localPlic\",\"capillary_method\":\"alphaHeightFunction\",\"surface_area_available\":false,\"energy_area_available\":false,\"volume_format\":\"tecplotDat\",\"frame_index\":\"meta/frames.csv\"";
}
TecplotWriter::Zone vofVolumeZone(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroVofStep* step,const LAMP::Interface::GeometricVofMeshTransfer* transfer){
    state.validateFor(mesh);validateStep(state,step,transfer);const auto& pub=publication(state);const auto& layout=state.values().layout();
    const auto velocity=state.sampleCellCenterVelocities(mesh);const auto t=state.timeLayer();
    const std::vector<double>* forceK=!transfer&&step&&step->hasCapillaryLoad()?&step->capillaryLoad().heightFunctionCellCurvatures():nullptr;
    const auto* forceEffective=forceK?&step->capillaryLoad().heightFunctionCellEffectiveCurvatures():nullptr;
    const auto* forceMethods=forceK?&step->capillaryLoad().heightFunctionCellCurvatureMethods():nullptr;
    TecplotWriter::Zone zone;zone.title="accepted native LocalPlic volume";zone.element=TecplotWriter::Element::Quadrilateral;
    zone.coordinateNames=layout.domain().axisymmetric()?std::array<std::string,2>{"R","Z"}:std::array<std::string,2>{"X","Y"};
    zone.solutionTime=t.time();zone.strandId=1;
    zone.auxiliary={{"AcceptedStep",std::to_string(t.step)},{"GeometryIdentity",std::to_string(pub.identity().value())},
        {"Primary","alpha"},{"Representation","localPlic"},{"FieldLocation","actual native cell averages/pressure; no nodal interpolation"},
        {"UnavailableCurvatureValue","0 is only a serialization placeholder when kappa_HF_forceAvailable is 0"}};
    // Retain every former VTU FieldData value as text, including exact identities
    // and the pressure/HF clocks. Unavailable force time remains explicitly nan.
    zone.auxiliary["RemappedForceUnavailable"]=transfer?"old physical force arrays belong to source mesh; new native values unavailable":"not_remapped";
    zone.auxiliary["CapillaryCurvatureMethods"]="0=unavailable;1=heightFunction;2=segmentFit;3=underresolvedIsolatedFragmentZero;4=mixedHeightFit;5=twoSegmentCircleFit";
    zone.auxiliary["CapillaryCurvatureMeaning"]="HF field is physical HF/fit estimate; effective field is actual force closure; method3 has physical unavailable, effective0 available";
    zone.auxiliary["VTK_TimeValue"]=scalarText(t.time());
    zone.auxiliary["VTK_geometry_identity"]=std::to_string(pub.identity().value());
    zone.auxiliary["VTK_accepted_step"]=std::to_string(t.step);
    zone.auxiliary["VTK_pressure_sample_time"]=scalarText(t.step?t.begin+.5*(t.end-t.begin):0.);
    zone.auxiliary["VTK_hf_force_time"]=step&&step->hasCapillaryLoad()?scalarText(step->capillaryLoad().phaseGeometry().timeLayer().time()):"nan";
    zone.auxiliary["VTK_hf_force_geometry_identity"]=std::to_string(step&&step->hasCapillaryLoad()?step->capillaryLoad().phaseGeometry().vofIdentity().value():0);
    for(const char* name:{"cell_id","cell_capacity","phase_volume_fraction","liquid_volume","gas_volume","plic_normal_X","plic_normal_Y","plic_normal_Z","plic_offset","plic_resolved","velocity_X","velocity_Y","velocity_Z","pressure","kappa_HF_force","kappa_HF_forceAvailable","kappa_effective_force","kappa_effective_forceAvailable","capillary_curvature_method"})
        zone.fields.push_back({name,TecplotWriter::Location::CellCentered,{}});
    for(int v=0;v<mesh.vertIndexNum();++v)zone.points.push_back(mesh.verts(mesh.vertIndexToLabel[v]).coord());
    for(int row=0;row<layout.pressureCount();++row){
        const auto native=layout.nativeCell(row);const auto& cell=mesh.cells(mesh.cellIndexToLabel[native]);const auto& p=pub.plicCells().at(native);
        std::vector<int> corners;for(auto v:cell.cellToVert)corners.push_back(mesh.verts(v).index());zone.connectivity.push_back(std::move(corners));
        const double k=forceK&&static_cast<std::size_t>(row)<forceK->size()?(*forceK)[row]:std::numeric_limits<double>::quiet_NaN();
        const double effective=forceEffective?forceEffective->at(row):std::numeric_limits<double>::quiet_NaN();
        const double method=forceMethods?static_cast<double>(forceMethods->at(row)):0.;
        const std::array<double,19> values{{static_cast<double>(native),static_cast<double>(pub.cells().at(native).capacity),pub.alpha().at(native),
            static_cast<double>(p.inventory[0]),static_cast<double>(p.inventory[1]),p.normal.x,p.normal.y,0.,static_cast<double>(p.offset),p.resolved?1.:0.,
            velocity.at(row).x,velocity.at(row).y,0.,state.values().pressure().at(row),std::isfinite(k)?k:0.,std::isfinite(k)?1.:0.,std::isfinite(effective)?effective:0.,std::isfinite(effective)?1.:0.,method}};
        for(std::size_t c=0;c<values.size();++c)zone.fields[c].values.push_back(values[c]);
    }
    zone.auxiliary["ForceGeometryIdentity"]=step&&step->hasCapillaryLoad()?std::to_string(step->capillaryLoad().phaseGeometry().vofIdentity().value()):"unavailable";
    return zone;
}
TecplotWriter::Zone actualForceGeometryZone(const QuadtreeMesh& sourceMesh,const HydroBalancedCsfLoad& force){
    const auto& phase=force.phaseGeometry();const auto& pub=phase.vofPublication();const auto t=pub.timeLayer();
    if(pub.geometryKind()!=LAMP::Interface::GeometricVofGeometryKind::LocalPlic)
        throw std::invalid_argument("actual force geometry diagnostic requires LocalPlic");
    force.validateVofFor(sourceMesh,force.pressureLoad(),pub.identity(),t);
    const auto& layout=force.layout();const auto& physical=force.heightFunctionCellCurvatures();
    const auto& effective=force.heightFunctionCellEffectiveCurvatures();const auto& methods=force.heightFunctionCellCurvatureMethods();
    if(!sameClock(t,phase.timeLayer())||!sameClock(t,force.timeLayer())||pub.identity()!=force.vofIdentity()
        ||pub.plicCells().size()!=static_cast<std::size_t>(sourceMesh.cellIndexNum())
        ||pub.cells().size()!=pub.plicCells().size()||physical.size()!=static_cast<std::size_t>(layout.pressureCount())
        ||effective.size()!=physical.size()||methods.size()!=physical.size())
        throw std::invalid_argument("actual force geometry diagnostic source/clock/native arrays differ");
    TecplotWriter::Zone zone;zone.title="diagnostic actual capillary force owner geometry";zone.element=TecplotWriter::Element::Quadrilateral;
    zone.coordinateNames=pub.domain().axisymmetric()?std::array<std::string,2>{"R","Z"}:std::array<std::string,2>{"X","Y"};
    zone.solutionTime=t.time();zone.strandId=1;
    const char* position=t.position==LAMP::TimeLayer::Position::Initial?"Initial":t.position==LAMP::TimeLayer::Position::Midpoint?"Midpoint":"End";
    zone.auxiliary={{"DiagnosticOnly","true; excluded from accepted frame index"},
        {"GeometryIdentity",std::to_string(pub.identity().value())},{"ForceGeometryIdentity",std::to_string(pub.identity().value())},
        {"ForceStep",std::to_string(t.step)},{"ForceBegin",scalarText(t.begin)},{"ForceEnd",scalarText(t.end)},
        {"ForcePosition",position},{"ForcePositionCode",std::to_string(static_cast<int>(t.position))},{"ForceTime",scalarText(t.time())},
        {"SourceMeshNativeCells",std::to_string(sourceMesh.cellIndexNum())},{"SourceMeshNativeVertices",std::to_string(sourceMesh.vertIndexNum())},
        {"GeometryMeaning","actual force publication and retained source mesh; never accepted End or remapped target geometry"},
        {"PlaneCoordinates","cell-centred local offset; rectangle bounds and nodes are actual source-mesh coordinates"},
        {"CurvatureMeaning","physical HF/fit estimate and separate effective force coefficient; availability flags distinguish placeholders"},
        {"CapillaryCurvatureMethods","0=unavailable;1=heightFunction;2=segmentFit;3=underresolvedIsolatedFragmentZero;4=mixedHeightFit;5=twoSegmentCircleFit"},
        {"UnavailableCurvatureValue","0 is a serialization placeholder only when the corresponding Available value is 0"},
        {"Units","code; raw independent stocks and actual measured PLIC volumes are unchanged"}};
    for(const char* name:{"cell_id","publication_active","force_pressure_dof","cell_level","rectangle_lower_X","rectangle_lower_Y","rectangle_upper_X","rectangle_upper_Y",
        "cell_capacity","liquid_volume","gas_volume","measured_liquid_volume","measured_gas_volume","plic_normal_X","plic_normal_Y","plic_offset","plic_resolved",
        "kappa_HF_force","kappa_HF_forceAvailable","kappa_effective_force","kappa_effective_forceAvailable","capillary_curvature_method"})
        zone.fields.push_back({name,TecplotWriter::Location::CellCentered,{}});
    for(int v=0;v<sourceMesh.vertIndexNum();++v)zone.points.push_back(sourceMesh.verts(sourceMesh.vertIndexToLabel[v]).coord());
    for(int native=0;native<sourceMesh.cellIndexNum();++native){
        const auto& cell=sourceMesh.cells(sourceMesh.cellIndexToLabel[native]);const auto& p=pub.plicCells().at(native);
        const auto lower=sourceMesh.verts(cell.cellToVert[0]).coord(),upper=sourceMesh.verts(cell.cellToVert[2]).coord();
        std::vector<int> corners;for(auto v:cell.cellToVert)corners.push_back(sourceMesh.verts(v).index());zone.connectivity.push_back(std::move(corners));
        const int row=layout.pressureDof(native);const bool active=pub.domain().cellRowByIndex().at(native)>=0;
        const double k=row>=0?physical.at(row):std::numeric_limits<double>::quiet_NaN();
        const double e=row>=0?effective.at(row):std::numeric_limits<double>::quiet_NaN();
        const std::array<double,22> values{{static_cast<double>(native),active?1.:0.,static_cast<double>(row),static_cast<double>(cell.level()),
            lower.x,lower.y,upper.x,upper.y,static_cast<double>(pub.cells().at(native).capacity),
            static_cast<double>(p.inventory[0]),static_cast<double>(p.inventory[1]),static_cast<double>(p.measured[0]),static_cast<double>(p.measured[1]),
            p.normal.x,p.normal.y,static_cast<double>(p.offset),p.resolved?1.:0.,std::isfinite(k)?k:0.,std::isfinite(k)?1.:0.,
            std::isfinite(e)?e:0.,std::isfinite(e)?1.:0.,row>=0?static_cast<double>(methods.at(row)):0.}};
        for(std::size_t c=0;c<values.size();++c)zone.fields[c].values.push_back(values[c]);
    }
    return zone;
}
TecplotWriter::Zone vofInterfaceZone(const HydroTwoPhaseState& state,const HydroVofStep* step,bool radial,const LAMP::Interface::GeometricVofMeshTransfer* transfer){
    const auto& pub=publication(state);validateStep(state,step,transfer);TecplotWriter::Zone zone;zone.title="accepted independent LocalPlic segments";zone.element=TecplotWriter::Element::Line;
    zone.coordinateNames=radial?std::array<std::string,2>{"R","Z"}:std::array<std::string,2>{"X","Y"};
    zone.fields={{radial?"NR":"NX",TecplotWriter::Location::Nodal,{}},{radial?"NZ":"NY",TecplotWriter::Location::Nodal,{}},{"KappaHFForce",TecplotWriter::Location::Nodal,{}},{"KappaHFForceAvailable",TecplotWriter::Location::Nodal,{}},{"KappaEffectiveForce",TecplotWriter::Location::Nodal,{}},{"KappaEffectiveForceAvailable",TecplotWriter::Location::Nodal,{}},{"CapillaryCurvatureMethod",TecplotWriter::Location::Nodal,{}}};
    zone.solutionTime=state.timeLayer().time();zone.strandId=1;zone.auxiliary={{"AcceptedStep",std::to_string(state.timeLayer().step)},{"GeometryIdentity",std::to_string(pub.identity().value())},{"Primary","alpha"},{"Representation","localPlic"},{"CurvatureMeaning","midpoint physical HF/fit estimates and separate effective force closure; not PLIC geometric curvature"}};
    const auto& layout=state.values().layout();const std::vector<double>* forceK=!transfer&&step&&step->hasCapillaryLoad()?&step->capillaryLoad().heightFunctionCellCurvatures():nullptr;
    const auto* forceEffective=forceK?&step->capillaryLoad().heightFunctionCellEffectiveCurvatures():nullptr;
    const auto* forceMethods=forceK?&step->capillaryLoad().heightFunctionCellCurvatureMethods():nullptr;
    zone.auxiliary["CapillaryCurvatureMethods"]="0=unavailable;1=heightFunction;2=segmentFit;3=underresolvedIsolatedFragmentZero;4=mixedHeightFit;5=twoSegmentCircleFit";
    zone.auxiliary["CapillaryCurvatureMeaning"]="HF field is physical HF/fit estimate; effective field is actual force closure; method3 has physical unavailable, effective0 available";
    std::size_t patches=0,unresolved=0;for(int row=0;row<layout.pressureCount();++row){const auto& p=pub.plicCells().at(layout.nativeCell(row));if(p.inventory[0]==0||p.inventory[1]==0)continue;
        if(!p.resolved){++unresolved;continue;}const auto segment=LAMP::Geometry::PlanarCellCut::interfaceSegment(p.rectangle,p.normal,p.offset);if(!segment){++unresolved;continue;}
        const long double cx=(static_cast<long double>(p.rectangle.lower.x)+p.rectangle.upper.x)/2,cy=(static_cast<long double>(p.rectangle.lower.y)+p.rectangle.upper.y)/2;
        const auto norm=std::hypot(p.normal.x,p.normal.y);const double k=forceK&&static_cast<std::size_t>(row)<forceK->size()?(*forceK)[row]:std::numeric_limits<double>::quiet_NaN();const auto base=static_cast<int>(zone.points.size());
        const double effective=forceEffective?forceEffective->at(row):std::numeric_limits<double>::quiet_NaN();
        const double method=forceMethods?static_cast<double>(forceMethods->at(row)):0.;
        for(const auto point:*segment){zone.points.push_back({static_cast<double>(cx+point.x),static_cast<double>(cy+point.y)});zone.fields[0].values.push_back(p.normal.x/norm);zone.fields[1].values.push_back(p.normal.y/norm);zone.fields[2].values.push_back(std::isfinite(k)?k:0);zone.fields[3].values.push_back(std::isfinite(k)?1:0);zone.fields[4].values.push_back(std::isfinite(effective)?effective:0);zone.fields[5].values.push_back(std::isfinite(effective)?1:0);zone.fields[6].values.push_back(method);}
        zone.connectivity.push_back({base,base+1});++patches;
    }
    zone.auxiliary["UnavailableCurvatureValue"]="0 is a serialization placeholder only when KappaHFForceAvailable is 0";
    zone.auxiliary["RemappedForceUnavailable"]=transfer?"old physical force arrays belong to source mesh; new native values unavailable":"not_remapped";
    zone.auxiliary["SegmentCount"]=std::to_string(patches);zone.auxiliary["UnresolvedMixedCells"]=std::to_string(unresolved);
    zone.auxiliary["ForceGeometryIdentity"]=step&&step->hasCapillaryLoad()?std::to_string(step->capillaryLoad().phaseGeometry().vofIdentity().value()):"unavailable";
    zone.auxiliary["ForceTime"]=step&&step->hasCapillaryLoad()?scalarText(step->capillaryLoad().phaseGeometry().timeLayer().time()):"unavailable";return zone;
}
Rows vofRows(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroVofStep* step,const HydrodynamicsInput& input,const Totals& totals,double compute,std::ios_base::fmtflags flags,std::streamsize precision,bool final,const LAMP::Interface::GeometricVofMeshTransfer* transfer){
    state.validateFor(mesh);const auto& pub=publication(state);validateStep(state,step,transfer);Rows out;out.next=totals;const auto& t=state.timeLayer();const auto* r=step?&step->report():nullptr;
    Sum liquid,mass;std::array<Sum,5> moments;long double closure=0;
    for(std::size_t i=0;i<pub.plicCells().size();++i){const auto& p=pub.plicCells()[i];const auto c=pub.cells()[i].capacity;liquid.add(p.inventory[0]);mass.add(input.vof.liquidDensity*p.inventory[0]+input.vof.gasDensity*p.inventory[1]);closure=std::max(closure,std::abs(p.inventory[0]+p.inventory[1]-c)/c);
        // Constant inventory density on the actual native rectangle. Preserve
        // extensive moments so centroid/shape choices stay in postprocessing.
        const long double dx=static_cast<long double>(p.rectangle.upper.x)-p.rectangle.lower.x,dy=static_cast<long double>(p.rectangle.upper.y)-p.rectangle.lower.y;
        const long double x=static_cast<long double>(p.rectangle.lower.x)+dx/2,y=static_cast<long double>(p.rectangle.lower.y)+dy/2;
        const long double meanX=input.radial?x+dx*dx/(12*x):x;
        const long double meanXX=x*x+dx*dx/(input.radial?4:12),meanYY=y*y+dy*dy/12;
        const std::array<long double,5> average{{meanX,y,meanXX,meanYY,meanX*y}};
        for(std::size_t k=0;k<moments.size();++k)moments[k].add(p.inventory[0]*average[k]);
    }
    if(!step)out.next.initialVolume=liquid.value;const auto drift=out.next.initialVolume>0?(liquid.value-out.next.initialVolume)/out.next.initialVolume:0;
    double umax=0,div=0;for(double u:state.values().normalVelocity())umax=std::max(umax,std::abs(u));for(double d:state.phaseGeometry().metrics().apply(mesh,state.values()))div=std::max(div,std::abs(d));
    long double defect=0;if(step)for(const auto& c:step->frozenMidpoint().material().report().components)defect=std::max(defect,c.maxAbsDefectRate);
    const auto range=forceRange(step);const auto kmax=std::max(std::abs(range.lower),std::abs(range.upper));
    if(t.step==0||t.step==1||final||t.step%input.run.logEverySteps==0){
        std::ostringstream p;p.flags(flags);p.precision(precision);
        p<<"[Step "<<t.step<<"] t="<<t.time()<<" dt="<<(r?r->deltaTime:0)
            <<" cells="<<mesh.cellIndexNum()<<" AC="<<(r?r->linearSolves:0)
            <<" div="<<div<<" dV/V="<<drift<<" umax="<<umax;
        if(input.run.detailedLog){
            p<<" R="<<defect;
            if(range.defined)p<<" kappa_force="<<kmax;
            if(step){const auto& q=step->transport().report();
                p<<" plic_sweeps="<<q.localDirectionalSweeps<<" plic_clips="<<q.localStripClips
                    <<" unresolved_donors="<<q.localUnresolvedDonors<<" phase_inventory_residual="<<q.maximumInventoryResidual;}
            if(step&&step->hasCapillaryLoad())p<<" force_time="<<step->capillaryLoad().phaseGeometry().timeLayer().time()
                <<" force_geometry="<<step->capillaryLoad().phaseGeometry().vofIdentity().value()
                <<" unresolved_capillary_fragments="<<step->capillaryLoad().curvatureReport().unresolvedIsolatedFragmentCells;
        }
        p<<" wall="<<compute<<"s";out.progress=p.str();
    }
    auto row=stream();row<<t.step<<','<<t.time()<<','<<(r?r->deltaTime:0)<<','<<pub.identity().value()<<','<<liquid.value<<','<<out.next.initialVolume<<','<<drift<<','<<mass.value<<','<<div<<',';if(step)row<<defect;else row<<"NA";row<<','<<umax<<',';if(range.defined)row<<kmax;else row<<"NA";row<<','<<closure<<',';if(state.retainedMass())row<<state.retainedMass()->report().maxRelativeLedgerDefect;else row<<"NA";for(const auto& moment:moments){row<<',';if(std::isfinite(moment.value))row<<moment.value;else row<<"NA";}appendVofMomentumBudget(row,step);row<<'\n';out.conservation=row.str();
    row=stream();row<<t.step<<','<<t.time()<<','<<(step?"AC":"initialization")<<",accepted,"<<(r?r->linearSolves:0)<<",NA,";if(step){const auto& s=step->frozenMidpoint().solveResult().report();row<<s.linearRefinements<<','<<s.algebraicRelativeResidual<<','<<s.freeMomentumDefect;}else row<<"NA,NA,NA";row<<",\n";out.solver=row.str();
    long double kinetic=0,body=0,viscous=0,capillary=0;bool kineticAvailable=state.retainedMass()!=nullptr;if(kineticAvailable){const auto& m=state.retainedMass()->masses();const auto& u=state.values().normalVelocity();for(std::size_t i=0;i<u.size();++i)kinetic+=.5L*m.at(i)*u[i]*u[i];}else if(umax==0)kineticAvailable=true;
    if(step){const auto& frozen=step->frozenMidpoint();const auto& u=step->solvedMidpointValues().normalVelocity();const auto& f=frozen.bodyLoads();const auto& k=frozen.viscosityMatrix();for(std::size_t i=0;i<u.size();++i)body+=u[i]*f.at(i);for(Eigen::Index j=0;j<k.outerSize();++j)for(HydroFrozenMidpoint::Matrix::InnerIterator it(k,j);it;++it)viscous+=static_cast<long double>(u[it.row()])*it.value()*u[it.col()];
        if(step->hasCapillaryLoad()){const auto& force=step->capillaryLoad().momentumRhsLoad();for(std::size_t i=0;i<u.size();++i)capillary+=u[i]*force.at(i);}out.next.bodyWork+=r->deltaTime*body;out.next.viscousLoss+=r->deltaTime*viscous;}
    row=stream();row<<t.step<<','<<t.time()<<',';if(kineticAvailable)row<<kinetic;else row<<"NA";row<<",NA,NA,";if(step)row<<body;else row<<0;row<<','<<out.next.bodyWork<<','<<out.next.viscousLoss<<',';out.energyPrefix=row.str();
    row=stream();row<<",NA,";if(step&&step->hasCapillaryLoad())row<<capillary<<",true,"<<capillary<<",true\n";else row<<"NA,false,NA,false\n";out.energySuffix=row.str();
    if(step&&step->hasCapillaryLoad()){const auto& force=step->capillaryLoad();const auto& h=force.curvatureReport();row=stream();row<<t.step<<','<<t.time()<<','<<force.phaseGeometry().timeLayer().time()<<','<<force.phaseGeometry().vofIdentity().value()<<',';if(range.defined)row<<range.lower<<','<<range.upper;else row<<"NA,NA";std::size_t loaded=0;for(double x:force.momentumRhsLoad())loaded+=x!=0;row<<','<<loaded<<','<<h.heightFunctionCells<<','<<h.heightFunctionColumns<<','<<h.heightFunctionCellQueries<<','<<h.heightFunctionAlternateDirections<<','<<h.heightFunctionExtendedCells<<','<<h.heightFunctionExtendedSlices<<','<<h.heightFunctionRefinedCells<<','<<h.heightFunctionRefinedSlices<<','<<h.heightFunctionFitCells<<','<<h.heightFunctionFitSegments<<','<<h.unresolvedIsolatedFragmentCells<<'\n';out.capillaryFile="vof_csf.csv";out.capillary=row.str();}
    return out;
}
std::string vofPerformance(const QuadtreeMesh& mesh,const HydroTwoPhaseState& state,const HydroVofStep* step,double compute,double io,double elapsed,const LAMP::Interface::GeometricVofMeshTransfer* transfer){
    state.validateFor(mesh);validateStep(state,step,transfer);
    const auto& physical=step?step->candidateState():state;
    auto out=stream();const auto append=[&](const char* name,double seconds,bool original){
        const auto& layout=original?physical.values().layout():state.values().layout();
        const auto cells=original?physical.phaseGeometry().vofPublication().plicCells().size():static_cast<std::size_t>(mesh.cellIndexNum());
        out<<state.timeLayer().step<<','<<state.timeLayer().time()<<','<<name<<','<<seconds<<','<<cells<<','<<layout.velocityCount()<<','<<layout.pressureCount()<<'\n';};
    append(step?"VOF_AC":"initialization",compute,true);append("io",io,false);append("step_total",elapsed,false);
    if(step){std::array<double,12> seconds{};for(const auto& s:step->report().stages){const auto& t=s.timing;const std::array<double,12> v{{t.transport,t.material,t.capillary,t.maxwell,t.frozenNs,t.diagnostics,t.phaseGeometry,t.materialTransport,t.dualMass,t.convection,t.viscosity,t.total}};for(std::size_t i=0;i<v.size();++i)seconds[i]+=v[i];}
        constexpr std::array<const char*,12> names{{"vof_transport","vof_material","vof_hf","vof_maxwell","vof_frozen_ns","vof_post_diagnostics","vof_phase_geometry","vof_material_transport","vof_dual_mass","vof_convection","vof_viscosity","vof_ac_stages"}};
        for(std::size_t i=0;i<seconds.size();++i)append(names[i],seconds[i],true);append("vof_other_and_retries",compute-seconds.back(),true);}
    return out.str();
}
}
