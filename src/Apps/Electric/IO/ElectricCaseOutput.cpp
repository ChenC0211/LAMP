#include "Electric/IO/ElectricCaseOutput.h"
#include "Core/IO/OutputProvenance.h"
#include "Modules/Electric/IO/ElectricTecplotFieldView.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
using namespace LAMP::Electric;
using Layer=LAMP::Interface::HermiteGalsTimeLayer;
namespace fs=std::filesystem;

void require(bool condition,const std::string& reason)
{
    if(!condition)throw std::runtime_error("[Electric/Output] "+reason);
}
bool sameClock(const Layer& a,const Layer& b)
{
    return a.step==b.step && a.begin==b.begin && a.end==b.end && a.position==b.position;
}
void text(const fs::path& path,const std::string& contents,bool append=false)
{
    std::ofstream out;out.exceptions(std::ios::failbit|std::ios::badbit);
    try {
        out.open(path,std::ios::binary|(append?std::ios::app:std::ios::trunc));
        out<<contents;out.close();
    } catch(const std::exception& error) {
        throw std::runtime_error("[Electric/Output] "+path.generic_string()+": "+error.what());
    }
}
std::string number(double value)
{
    auto out=LAMP::IO::stream();out<<value;return out.str();
}
std::string fileName(const std::string& prefix,const char* kind,int step)
{
    std::ostringstream out;out<<prefix<<'_'<<kind<<"_step"
        <<std::setw(6)<<std::setfill('0')<<step<<".dat";return out.str();
}
std::array<double,3> basis(double coordinate)
{
    return {{1,coordinate,.5*(3*coordinate*coordinate-1)}};
}
double polynomial(const std::array<double,3>& coefficients,double coordinate)
{
    const auto b=basis(coordinate);
    return coefficients[0]*b[0]+coefficients[1]*b[1]+coefficients[2]*b[2];
}
double polynomial(const std::array<double,4>& coefficients,double coordinate)
{
    return polynomial(std::array<double,3>{{coefficients[0],coefficients[1],coefficients[2]}},coordinate)
        +coefficients[3]*.5*(5*coordinate*coordinate*coordinate-3*coordinate);
}
void sameArithmetic(long double defect,long double scale,double tolerance,const char* reason)
{
    require(std::isfinite(defect) && std::isfinite(scale)
        && (scale>0 ? defect<=tolerance*scale : defect==0),reason);
}
TecplotWriter::Zone zone(const ElectricPhaseResult& field,const char* title,int strand)
{
    TecplotWriter::Zone result;result.title=title;
    result.element=TecplotWriter::Element::OrderedPoints;
    result.coordinateNames=field.axisymmetric()?std::array<std::string,2>{{"r [code]","z [code]"}}
        :std::array<std::string,2>{{"x [code]","y [code]"}};
    const auto physical=field.timeLayer(),source=field.geometry().timeLayer();
    result.solutionTime=physical.time();result.strandId=strand;
    result.auxiliary={{"Units","code"},{"Method","HHO"},
        {"Backend","Eigen::SimplicialLDLT (SPD)"},{"InterfaceMotion","fixed"},
        {"AcceptedStep",std::to_string(physical.step)},
        {"ElectricClockPosition",physical.position==Layer::Position::Initial?"Initial":"End"},
        {"ElectricClockBegin",number(physical.begin)},{"ElectricClockEnd",number(physical.end)},
        {"GeometryIdentity",std::to_string(field.geometry().geometryIdentity().value())},
        {"GeometrySourceStep",std::to_string(source.step)},
        {"GeometrySourceClockPosition","Initial"},
        {"GeometrySourceClockBegin",number(source.begin)},{"GeometrySourceClockEnd",number(source.end)},
        {"CoordinateSystem",field.axisymmetric()?"axisymmetricRZ":"cartesianXY"}};
    return result;
}
TecplotWriter::Zone phaseZone(const ElectricPhaseResult& field,bool liquid,int phaseId,double sigma)
{
    auto output=zone(field,liquid?"electric liquid phase samples":"electric gas phase samples",liquid?1:2);
    output.auxiliary["PhaseId"]=std::to_string(phaseId);
    output.auxiliary["PhaseSign"]=liquid?"negative phi":"positive phi";
    output.auxiliary["Sampling"]="one retained positive-weight volume quadrature point per phase element; nearest element center; one-sided P3 sample, not an average";
    const bool radial=field.axisymmetric();
    output.fields={{"V [code]",TecplotWriter::Location::Nodal,{}},
        {"epsilon [code]",TecplotWriter::Location::Nodal,{}},
        {"sigma [code]",TecplotWriter::Location::Nodal,{}},
        {radial?"Er [code]":"Ex [code]",TecplotWriter::Location::Nodal,{}},
        {radial?"Ez [code]":"Ey [code]",TecplotWriter::Location::Nodal,{}}};
    const auto& elements=field.geometry().elements();
    for(std::size_t element=0;element<elements.size();++element) {
        if(field.liquid(element)!=liquid)continue;
        const auto& volume=elements[element].geometry.volume;
        require(!volume.empty(),"phase element has no retained physical volume quadrature");
        const auto center=field.elementCenter(element);
        const auto distance=[&](const LAMP::ElectricDiscretization::CellQuadraturePoint& q) {
            const auto delta=q.point-center;return delta.x*delta.x+delta.y*delta.y;
        };
        const auto selected=std::min_element(volume.begin(),volume.end(),
            [&](const auto& a,const auto& b){return distance(a)<distance(b);});
        const auto point=selected->point;
        const auto gradient=field.gradientInElement(element,point);
        output.points.push_back(point);
        output.fields[0].values.push_back(field.potentialInElement(element,point));
        output.fields[1].values.push_back(field.epsilonByPhase()[liquid?0:1]);
        output.fields[2].values.push_back(sigma);
        output.fields[3].values.push_back(-gradient.x);
        output.fields[4].values.push_back(-gradient.y);
    }
    require(!output.points.empty(),"phase has no retained physical volume samples");
    return output;
}
TecplotWriter::Zone interfaceZone(const ElectricPhaseResult& field,const ElectricMaxwellTrace& maxwell)
{
    auto output=zone(field,"electric interface electric and Maxwell samples",3);
    output.auxiliary["Sampling"]="original Maxwell quadrature on the shared interface skeleton; no refit";
    output.auxiliary["NormalOrientation"]="negative phi to positive phi";
    output.auxiliary["TractionMeaning"]="(T_plus - T_minus) n; tangent=(-ny,nx)";
    output.auxiliary["ChargeMeaning"]="q_h from the original integrated Q via the shared weighted P2 load";
    output.auxiliary["PotentialMeaning"]="same-source shared P2 trace with retained weighted cubic completion";
    const bool radial=field.axisymmetric();
    output.fields={{radial?"nr":"nx",TecplotWriter::Location::Nodal,{}},
        {radial?"nz":"ny",TecplotWriter::Location::Nodal,{}},
        {"V_Gamma [code]",TecplotWriter::Location::Nodal,{}},
        {"q_h [code]",TecplotWriter::Location::Nodal,{}},
        {"Et [code]",TecplotWriter::Location::Nodal,{}},
        {"Dn_liquid [code]",TecplotWriter::Location::Nodal,{}},
        {"Dn_gas [code]",TecplotWriter::Location::Nodal,{}},
        {"Maxwell_normal [code]",TecplotWriter::Location::Nodal,{}},
        {"Maxwell_tangent [code]",TecplotWriter::Location::Nodal,{}},
        {radial?"Maxwell_r [code]":"Maxwell_x [code]",TecplotWriter::Location::Nodal,{}},
        {radial?"Maxwell_z [code]":"Maxwell_y [code]",TecplotWriter::Location::Nodal,{}}};
    for(const auto face:maxwell.interfaceFaces()) {
        const auto potential=maxwell.potentialPolynomialP3(face);
        for(const auto& q:maxwell.quadrature(face)) {
            output.points.push_back(q.point);
            const std::array<double,11> values{{q.normal.x,q.normal.y,
                polynomial(potential,q.coordinate),q.chargeDensity,q.tangentialField,
                q.normalDisplacement[0],q.normalDisplacement[1],q.normalTraction,
                q.tangentialTraction,q.tractionJump.x,q.tractionJump.y}};
            for(std::size_t variable=0;variable<values.size();++variable)
                output.fields[variable].values.push_back(values[variable]);
        }
    }
    require(!output.points.empty(),"interface has no retained Maxwell samples");
    return output;
}
} // namespace

ElectricCaseOutput::ElectricCaseOutput(const ElectricCaseInput& input,int finalStep)
    :run_(input.run),finalStep_(finalStep),radial_(input.radial),
    phaseIds_{{input.negativePhase,input.positivePhase}},epsilon_(input.electric.epsilonByPhase),
    conductivity_(input.conductivity),sourceTolerance_(input.linear.relativeConservationTolerance),
    projectionTolerance_(input.maxwell.relativeProjectionTolerance)
{
    require(finalStep>0,"requires a positive final step from the App");
    require(!fs::exists("output") || (fs::is_directory("output") && fs::is_empty("output")),
        "output already exists; use a fresh run directory");
    for(const auto* part:{"fields/volume/samples","fields/interface","series","meta/input","meta/effective"})
        fs::create_directories(fs::path("output")/part);
    text("output/fields/README.txt",
        "MAIN VIEW: volume/"+run_.outputPrefix+"_fields_stepXXXXXX.dat\n"
        "Load one accepted-step file in Tecplot. Both phase zones are in this file.\n"
        "Use 2D Cartesian, Contour Flood, V [code] or E_mag [code]; turn Scatter off.\n"
        "Electric field components are -grad of the actual one-sided electric P3 potential.\n"
        "The visualization triangles follow the same retained interface with linear chords;\n"
        "nodes remain local to each phase element, preserving physical/numerical jumps.\n"
        "The display mesh is not the electric integration geometry or a conservation diagnostic.\n"
        "XY axes: x/y. Axisymmetric axes: r/z (meridional view, no mirror added).\n\n"
        "VERIFICATION SAMPLES: volume/samples/*_phase_*_samples_stepXXXXXX.dat\n"
        "One original positive volume quadrature node per phase element; scatter samples, not contours.\n"
        "INTERFACE DIAGNOSTICS: interface/*_interface_samples_stepXXXXXX.dat\n"
        "Original q_h, common Et, two-sided equilibrated Dn and Maxwell quadrature; scatter samples.\n"
        "Accepted-step/time/source links: ../meta/frames.csv. Time histories: ../series/.\n");
    if(run_.writeFields)
        std::cout<<"[OUTPUT] Tecplot contours: output/fields/volume/"<<run_.outputPrefix
            <<"_fields_stepXXXXXX.dat (V, E_mag; both phases)"<<std::endl;
    auto effective=LAMP::IO::stream(),hashes=LAMP::IO::stream(),trace=LAMP::IO::stream();
    hashes<<"{\n";int count=0;
    for(const auto& document:input.documents) {
        const auto name=fs::path(document.file).filename();
        text(fs::path("output/meta/input")/name,document.sourceText);
        text(fs::path("output/meta/effective")/name,document.effectiveText);
        effective<<"// "<<name.generic_string()<<'\n'<<document.effectiveText<<'\n';
        if(count++)hashes<<",\n";
        hashes<<LAMP::IO::quoted(document.file)<<": {\"source_sha256\": "
            <<LAMP::IO::quoted(LAMP::IO::sha256(document.sourceText))
            <<", \"effective_sha256\": "<<LAMP::IO::quoted(LAMP::IO::sha256(document.effectiveText))<<'}';
    }
    hashes<<"\n}\n";text("output/meta/input_hashes.json",hashes.str());
    text("output/meta/effective_config.dict",effective.str());
    trace<<"field,file,line,raw_value,effective_value,units,consumer,defaulted,used,unused_reason\n";
    for(const auto& entry:input.trace.entries)
        trace<<LAMP::IO::csv(entry.field)<<','<<LAMP::IO::csv(entry.file)<<','<<entry.line<<','
            <<LAMP::IO::csv(entry.rawValue)<<','<<LAMP::IO::csv(LAMP::IO::value(entry.value))<<','
            <<LAMP::IO::csv(entry.units)<<','<<LAMP::IO::csv(entry.consumer)<<','
            <<(entry.defaulted?"true":"false")<<','<<(entry.used?"true":"false")<<','
            <<LAMP::IO::csv(entry.unusedReason)<<'\n';
    text("output/meta/configuration_trace.csv",trace.str());
    text("output/meta/build_info.txt","compiler="+LAMP::IO::compiler()
        +"\nbuild_type="+LAMP::IO::buildType()+"\ncompiled="+__DATE__+" "+__TIME__
        +"\nunits=code\nmethod=HHO\nbackend=Eigen::SimplicialLDLT (SPD)"
        "\ncharge_integration=BE startup; variable BDF2 thereafter\ninterface_motion=fixed"
        "\ngeometry_source_clock=original Initial FirstJet publication"
        "\nphysical_output_clock=matching Q/V/I/Maxwell Initial or accepted End"
        "\nphase_contour=two same-source phase triangle zones in one file; element-local P3 nodal samples; no phase or cell averaging"
        "\ncontour_geometry=linear chord visualization of retained boundaries, not integration/acceptance geometry"
        "\nphase_field=one actual positive-weight volume node nearest element center per phase element; one-sided P3 point sample, not an average"
        "\ninterface_field=retained weighted cubic trace and original Maxwell quadrature"
        "\nnormal_displacement=equilibrated electric flux; not a raw gradient replacement\n"
        +"series_profile="+(run_.diagnosticSeries?"diagnostic":"physical")+"\n");
    text("output/meta/frames.csv","step,time,geometry_identity,liquid,gas,interface,cloud\n");
    history_=std::make_unique<TecplotWriter::History>(
        fs::path("output/series")/(run_.outputPrefix+"_electric.dat"),
        "Accepted sharp-interface electrostatics physical states",
        std::vector<std::string>{"step","time","geometry_source_step","geometry_source_time",
            "total_Q","panel_current","exterior_current","bulk_current","current_balance",
            "current_balance_relative","linear_residual","local_Gauss_relative","global_Gauss",
            "global_Gauss_relative","time_update_present","temporal_order","dt","previous_dt",
            "a0","a1","a2","Schur_iterations","Schur_relative","panel_equation_relative",
            "charge_equation","charge_equation_relative","Maxwell_charge_jump_relative"},
        std::map<std::string,std::string>{{"Units","code"},{"Method","HHO"},
            {"PhysicalClock","same accepted Q/V/I/Maxwell layer"},
            {"GeometryClock","original Initial publication; fixed geometry"},
            {"CurrentMeaning","Ohmic panel influx and exterior outflow; no emission law"},
            {"InitialTimeEquation","time_update_present=0; time-equation columns are undefined placeholders"}});
}

void ElectricCaseOutput::write(const QuadtreeMesh& mesh,const ElectricPhaseResult& field,
    const ElectricPanelCurrentOperator& current,const ElectricMaxwellTrace& maxwell,int step,
    const ElectricChargeTimeReport* report)
{
    const auto physical=field.timeLayer(),source=field.geometry().timeLayer();
    const auto& publication=field.geometry().publication();
    using Clock=std::chrono::steady_clock;
    const auto elapsed=[](Clock::time_point started) {
        return std::chrono::duration<double>(Clock::now()-started).count();
    };
    auto started=Clock::now();outputTiming_={};
    const auto identity=field.geometry().geometryIdentity();
    const std::string where="step="+std::to_string(step)+" geometry="+std::to_string(identity.value())+": ";
    (void)physical.time();(void)source.time();
    require(step>=0 && step<=finalStep_ && step==lastWrittenStep_+1,
        where+"output requires the next actual accepted step, without duplicate files");
    require(physical.step==step && (step==0
        ?physical.position==Layer::Position::Initial
        :physical.position==Layer::Position::End),where+"actual electric physical clock is not Initial/End");
    require(source.position==Layer::Position::Initial,where+"fixed geometry source clock must remain Initial");
    require((step==0)==(report==nullptr),where+"time report must belong to each advanced physical layer only");
    if(report)require(report->dt==physical.end-physical.begin,
        where+"time report interval differs from the accepted electric clock");
    require(field.axisymmetric()==radial_ && field.epsilonByPhase()==epsilon_,
        where+"physical coordinate measure or material coefficients differ from the configured output");
    field.validateSource(mesh,publication);
    current.geometry().validateSource(mesh,publication);
    current.momentOperator().geometry().validateSource(mesh,publication);
    maxwell.validateSource(mesh,publication,physical);
    require(current.geometry().geometryIdentity()==identity
        && current.momentOperator().geometry().geometryIdentity()==identity
        && maxwell.geometry().geometryIdentity()==identity
        && sameClock(current.geometry().timeLayer(),source)
        && sameClock(current.momentOperator().timeLayer(),source)
        && sameClock(maxwell.geometrySourceTimeLayer(),source)
        && sameClock(current.timeLayer(),physical) && sameClock(maxwell.timeLayer(),physical),
        where+"Q/V/I/Maxwell geometry identity or full clock differs");
    require(current.epsilonByPhase()==epsilon_ && current.conductivityByPhase()==conductivity_,
        where+"panel current uses different physical material coefficients");
    const auto* load=field.surfaceChargeLoad();
    require(load!=nullptr,where+"formal sharp-interface electrostatics output requires the retained integrated-charge load");
    load->validateGeometry(field.geometry());
    const auto& charge=load->panelCharge();
    require(charge.size()==current.panelCurrent().size()
        && charge.size()==current.momentOperator().panelCount(),
        where+"integrated Q and current panel ordering differ");
    require(load->momentOperator().geometry().geometryIdentity()==identity
        && sameClock(load->momentOperator().timeLayer(),source),
        where+"integrated Q uses a different P_Q source");
    // A current has no public source-result identity. This read-only affine
    // query checks that the caller supplied the current for this actual V.
    const auto applied=current.apply(field.skeletonSolution());
    require(applied.size()==current.panelCurrent().size(),where+"current output size mismatch");
    long double currentDefect=0,currentScale=0,totalCharge=0;
    for(std::size_t panel=0;panel<charge.size();++panel) {
        totalCharge+=charge[panel];
        currentDefect+=std::abs(static_cast<long double>(applied[panel])-current.panelCurrent()[panel]);
        currentScale+=std::abs(static_cast<long double>(current.offset()[panel]))
            +std::abs(current.panelCurrent()[panel]);
    }
    const auto& matrix=current.matrix();const auto& skeleton=field.skeletonSolution();
    for(int column=0;column<matrix.outerSize();++column)
        for(ElectricPanelCurrentOperator::Matrix::InnerIterator entry(matrix,column);entry;++entry)
            currentScale+=std::abs(static_cast<long double>(entry.value())
                *skeleton[static_cast<std::size_t>(entry.col())]);
    sameArithmetic(currentDefect,currentScale,sourceTolerance_,
        "actual potential and panel-current affine receipt differ");
    // Maxwell's clock/geometry checks do not expose its originating result.
    // Compare its retained original-node flux/charge/shared trace receipts.
    // This neither prepares a new trace nor alters physical acceptance.
    long double displacementDefect=0,displacementScale=0;
    for(const auto face:maxwell.interfaceFaces()) {
        const auto minus=field.normalDisplacementPolynomialP2(face,true);
        const auto plus=field.normalDisplacementPolynomialP2(face,false);
        const auto shared=field.skeletalFaceTrace(face);
        const auto completed=maxwell.potentialPolynomialP3(face);
        const auto chargeMoments=field.interfaceChargeMoments(face);
        std::array<long double,3> projectedCharge{{0,0,0}},projectionDefect{{0,0,0}};
        long double projectionScale=0,chargeScale=0;
        for(const auto& node:maxwell.quadrature(face)) {
            const auto b=basis(node.coordinate);
            const std::array<double,2> expected{{polynomial(minus,node.coordinate),
                polynomial(plus,node.coordinate)}};
            for(std::size_t phase=0;phase<2;++phase) {
                displacementDefect+=node.physicalWeight*std::abs(
                    static_cast<long double>(expected[phase])-node.normalDisplacement[phase]);
                displacementScale+=node.physicalWeight*(std::abs(static_cast<long double>(expected[phase]))
                    +std::abs(node.normalDisplacement[phase]));
            }
            const double v=polynomial(completed,node.coordinate),vShared=polynomial(shared,node.coordinate);
            for(std::size_t mode=0;mode<3;++mode) {
                projectedCharge[mode]+=static_cast<long double>(node.physicalWeight)*b[mode]*node.chargeDensity;
                chargeScale+=std::abs(static_cast<long double>(node.physicalWeight)*b[mode]*node.chargeDensity);
                projectionDefect[mode]+=static_cast<long double>(node.physicalWeight)*b[mode]*(v-vShared);
                projectionScale+=node.physicalWeight*std::abs(b[mode])
                    *(std::abs(static_cast<long double>(v))+std::abs(vShared));
            }
        }
        long double chargeDefect=0,traceDefect=0;
        for(std::size_t mode=0;mode<3;++mode) {
            chargeDefect=std::max(chargeDefect,std::abs(projectedCharge[mode]-chargeMoments[mode]));
            chargeScale+=std::abs(chargeMoments[mode]);
            traceDefect=std::max(traceDefect,std::abs(projectionDefect[mode]));
        }
        sameArithmetic(chargeDefect,chargeScale,sourceTolerance_,
            "integrated Q and retained Maxwell charge moments differ");
        sameArithmetic(traceDefect,projectionScale,projectionTolerance_,
            "actual potential and retained Maxwell shared-trace projection differ");
    }
    sameArithmetic(displacementDefect,displacementScale,sourceTolerance_,
        "actual potential and retained Maxwell normal displacement differ");
    outputTiming_.validation=elapsed(started);
    if(run_.writesFieldsAtStep(step,finalStep_)) {
        started=Clock::now();
        const auto liquid="volume/samples/"+fileName(run_.outputPrefix,"phase_liquid_samples",step);
        const auto gas="volume/samples/"+fileName(run_.outputPrefix,"phase_gas_samples",step);
        const auto surface="interface/"+fileName(run_.outputPrefix,"interface_samples",step);
        const auto cloud="volume/"+fileName(run_.outputPrefix,"fields",step);
        auto contours=LAMP::Electric::buildElectricTecplotFieldView(mesh,field,phaseIds_);
        for(auto& phase:contours)phase.auxiliary["InterfaceMotion"]="fixed";
        outputTiming_.cloudBuild=elapsed(started);started=Clock::now();
        TecplotWriter::write(fs::path("output/fields")/cloud,contours);
        outputTiming_.cloudWrite=elapsed(started);started=Clock::now();
        TecplotWriter::write(fs::path("output/fields")/liquid,
            phaseZone(field,true,phaseIds_[0],conductivity_[0]));
        TecplotWriter::write(fs::path("output/fields")/gas,
            phaseZone(field,false,phaseIds_[1],conductivity_[1]));
        TecplotWriter::write(fs::path("output/fields")/surface,interfaceZone(field,maxwell));
        outputTiming_.samples=elapsed(started);
        auto frame=LAMP::IO::stream();frame<<step<<','<<physical.time()<<','<<identity.value()<<','
            <<LAMP::IO::csv(liquid)<<','<<LAMP::IO::csv(gas)<<','<<LAMP::IO::csv(surface)<<','<<LAMP::IO::csv(cloud)<<'\n';
        text("output/meta/frames.csv",frame.str(),true);
    }
    started=Clock::now();
    const auto& electric=field.report();const auto& currentReport=current.report();
    history_->append(TecplotWriter::row(step,physical.time(),source.step,source.time(),totalCharge,
        currentReport.totalPanelCurrent,currentReport.exteriorCurrent,currentReport.bulkCurrent,
        currentReport.globalCurrentBalance,currentReport.relativeGlobalCurrentBalance,
        electric.skeleton.relativeLinearResidual,electric.maxRelativeLocalGaussResidual,
        electric.globalGaussResidual,electric.relativeGlobalGaussResidual,report!=nullptr,
        report?report->temporalOrder:0,report?report->dt:0,report?report->previousDt:0,
        report?report->a0:0,report?report->a1:0,report?report->a2:0,
        report?report->iterations:0,report?report->relativeSchurResidual:0,
        report?report->maxRelativePanelEquationResidual:0,report?report->integratedEquationResidual:0,
        report?report->relativeTotalEquationResidual:0,maxwell.report().maxRelativeChargeJumpResidual));
    lastWrittenStep_=step;
    outputTiming_.series=elapsed(started);
}

void ElectricCaseOutput::writePerformance(int step,double time,
    const ElectricPerformanceReport& app,const LAMP::Electric::ElectricPhaseResult& field,
    const LAMP::Electric::ElectricChargeTimeReport* charge)
{
    if(!run_.diagnosticSeries)return;
    // Parents and children are both recorded; inclusive child clocks must not
    // be added to their parents. These receipts never participate in acceptance.
    const auto& geometry=field.geometry().report();
    const auto& g=geometry.timings;
    const auto& e=field.report().skeleton.timing;
    const LAMP::Electric::ElectricChargeTimeReport empty;
    const auto& q=charge?charge->timing:empty.timing;
    const auto initial=[&](double value){return step==0?value:0.;};
    const double classified=app.inputSeconds+app.outputStartupSeconds+app.meshSeconds
        +app.publicationSeconds+app.geometrySeconds+app.panelChargeSeconds
        +app.momentPrepareSeconds+app.initialSolveSeconds+app.chargePrepareSeconds
        +app.initialCurrentSeconds+app.chargeAdvanceSeconds+app.maxwellSeconds
        +app.commitSeconds+app.outputSeconds;
    const std::vector<std::pair<const char*,double>> values{
        {"step",static_cast<double>(step)},{"time",time},
        {"phase_elements",static_cast<double>(field.phaseElementCount())},
        {"free_skeleton_dofs",static_cast<double>(field.report().skeleton.freeSkeletonDofs)},
        {"panels",static_cast<double>(field.surfaceChargeLoad()->panelCharge().size())},
        {"volume_points",static_cast<double>(geometry.volumePoints)},
        {"surface_points",static_cast<double>(geometry.surfacePoints)},
        {"quadrature_refinement",static_cast<double>(geometry.acceptedRefinement)},
        {"schur_iterations",charge?static_cast<double>(charge->iterations):0},
        {"electric_solves",charge?static_cast<double>(charge->electricSolves):0},
        {"total_s",app.totalSeconds},{"input_s",app.inputSeconds},
        {"output_startup_s",app.outputStartupSeconds},{"mesh_s",app.meshSeconds},
        {"publication_s",app.publicationSeconds},{"phase_geometry_s",app.geometrySeconds},
        {"geometry_phase_quadrature_s",initial(g.phaseQuadrature_s)},
        {"geometry_cartesian_faces_s",initial(g.cartesianFaces_s)},
        {"geometry_interface_faces_s",initial(g.interfaceFaces_s)},
        {"geometry_local_operators_s",initial(g.localOperators_s)},
        {"geometry_matrix_comparison_s",initial(g.matrixComparison_s)},
        {"geometry_remaining_s",initial(g.remaining_s)},
        {"panel_charge_s",app.panelChargeSeconds},{"moment_prepare_s",app.momentPrepareSeconds},
        {"initial_electric_s",app.initialSolveSeconds},
        {"electric_boundary_s",initial(e.boundarySeconds)},
        {"electric_local_operators_s",initial(e.localOperatorsSeconds)},
        {"electric_local_loads_s",initial(e.localLoadsSeconds)},
        {"electric_assembly_s",initial(e.assemblySeconds)},
        {"electric_factor_s",initial(e.factorSeconds)},{"electric_solve_s",initial(e.solveSeconds)},
        {"electric_recovery_s",initial(e.recoverySeconds)},
        {"electric_gauss_affine_s",initial(e.gaussAffineSeconds)},
        {"charge_prepare_s",app.chargePrepareSeconds},{"initial_current_s",app.initialCurrentSeconds},
        {"charge_advance_s",app.chargeAdvanceSeconds},{"charge_schur_s",q.schurSeconds},
        {"charge_schur_inner_electric_s",q.schurElectricSolveSeconds},
        {"charge_moment_load_s",q.momentLoadSeconds},{"charge_final_electric_s",q.finalElectricSolveSeconds},
        {"charge_recovery_s",q.recoverySeconds},{"charge_current_s",q.currentSeconds},
        {"charge_balance_s",q.balanceSeconds},{"maxwell_s",app.maxwellSeconds},
        {"commit_s",app.commitSeconds},{"output_s",app.outputSeconds},
        {"output_validation_s",outputTiming_.validation},{"output_cloud_build_s",outputTiming_.cloudBuild},
        {"output_cloud_write_s",outputTiming_.cloudWrite},{"output_samples_s",outputTiming_.samples},
        {"output_series_s",outputTiming_.series},{"app_unclassified_s",app.totalSeconds-classified}
    };
    auto out=LAMP::IO::stream();
    if(step==0) {
        for(std::size_t j=0;j<values.size();++j)out<<(j?",":"")<<values[j].first;
        out<<'\n';
    }
    for(std::size_t j=0;j<values.size();++j) {
        require(std::isfinite(values[j].second),"nonfinite passive timing receipt");
        out<<(j?",":"")<<values[j].second;
    }
    out<<'\n';text("output/series/performance.csv",out.str(),step!=0);
}
