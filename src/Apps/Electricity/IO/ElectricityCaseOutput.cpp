#include "Electricity/IO/ElectricityCaseOutput.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

ElectricityCaseOutput::ElectricityCaseOutput(const ElectricityCaseInput& input,int finalStep)
    : run_(input.run), finalStep_(finalStep), radial_(input.domain.electric.axisymmetric)
{
    if(finalStep<=0)throw std::invalid_argument("Electricity output requires a positive final step from the App");
    namespace fs=std::filesystem;
    if (fs::exists("output") && (!fs::is_directory("output") || !fs::is_empty("output")))
        throw std::runtime_error("Electricity output already exists; use a fresh run directory");
    for (const auto* part : {"fields","series","meta"}) fs::create_directories(fs::path("output")/part);
    fs::copy("system","output/meta/input",fs::copy_options::recursive);
    std::ofstream effective("output/meta/electricity_effective.txt");
    effective<<std::setprecision(17)<<"model="<<input.domain.electric.model
        <<"\ninterfaceDiscretization="<<input.domain.electric.interfaceDiscretization
        <<"\ncoordinateSystem="<<(radial_?"axisymmetricRZ":"cartesianXY")
        <<"\nunits=code\nbackend=SuperLU\ninterfaceMotion=fixed"
        <<"\nlinearResidualTolerance="<<input.domain.electric.linearResidualTolerance
        <<"\ninterfaceHalfWidthConsumed="<<(input.domain.electric.interfaceDiscretization=="diffuseTanh")
        <<"\nchargeClosureTolerance="<<input.domain.electric.chargeClosureTolerance.value_or(16e-9)
        <<"\nchargeClosureConsumed="<<(input.domain.electric.model=="leakyDielectric")
        <<"\noutputTime=coeval electric field and charge\n";
    if(input.domain.electric.interfaceDiscretization=="diffuseTanh")
        effective<<"interfaceHalfWidth="<<input.domain.diffuseThickness
            <<"\nboundaryFlux="<<input.domain.electric.boundaryFluxTolerance.value_or(1e-9)<<'\n';
    else {
        const auto checks=input.domain.electric.sharpFluxChecks.value_or(DielectricFluxChecks());
        effective<<"reconstructedPoisson="<<checks.reconstructedPoisson<<"\ngaussClosure="<<checks.gaussClosure<<'\n';
    }
    if(input.domain.ignoredSharpHalfWidth)
        effective<<"deprecatedIgnoredSharpHalfWidth="<<*input.domain.ignoredSharpHalfWidth<<'\n';
    effective.close();
    if (!effective) throw std::runtime_error("cannot write effective electric configuration");
    history_=std::make_unique<TecplotWriter::History>(
        fs::path("output/series")/(run_.outputPrefix+"_electric.dat"),"Accepted standalone electric states",
        std::vector<std::string>{"step","time","superlu","linear_residual","charge_present","total_charge",
            "charge_update_present","current_source_time","old_boundary_current","new_boundary_current",
            "electrode_charge_residual","max_electric_field","electric_energy"},
        std::map<std::string,std::string>{{"Units","code"},{"FieldTime","coeval potential and accepted charge"},{"CurrentMeaning","boundary conduction, not emission"}});
}
void ElectricityCaseOutput::write(const QuadtreeMesh& mesh,
    const Quantum::Electricity::DielectricEvaluation& evaluation,
    int step,const Quantum::Electricity::DielectricChargeUpdate* update)
{
    const double time=step*run_.dt;
    const double previous=(step-1)*run_.dt; // Source time in the charge-update history.
    const auto& report=evaluation.electric.report;
    if(!report.success || !evaluation.electric.fields || evaluation.time.step!=step+1 || evaluation.time.forceTime!=time)
        throw std::runtime_error("electric output requires the matching successful field time layer");
    const auto& charge=evaluation.charge;
    const auto& e=*evaluation.electric.fields;
    e.meshIdentity.validate(mesh);
    if(e.axisymmetric!=radial_) throw std::runtime_error("electric result/output coordinate measures disagree");
    if(charge) {
        charge->meshIdentity.validate(mesh);
        if(charge->completedStep!=step || charge->time!=time)
            throw std::runtime_error("electric output requires coeval field and charge");
    }
    if(update && (update->charge!=charge || update->sourceTime.step!=step || update->sourceTime.dt!=evaluation.time.dt))
        throw std::runtime_error("charge update does not belong to the published electric state");
    long double total=0.; if(charge) for(double q:charge->integratedCharge) total+=q;
    if (run_.writesFieldsAtStep(step,finalStep_))
    {
        TecplotWriter::Zone zone;
        zone.title="Standalone electric field"; zone.solutionTime=time; zone.strandId=1;
        zone.auxiliary={{"Units","code"},{"CoordinateSystem",radial_?"axisymmetricRZ":"cartesianXY"},
            {"FieldTime","electric field and charge at this time"},{"InterfaceMotion","fixed"}};
        for(int i=0;i<mesh.vertIndexNum();++i) zone.points.push_back(mesh.verts(mesh.vertIndexToLabel[i]).coord());
        for(int i=0;i<mesh.cellIndexNum();++i) {
            const auto& cell=mesh.cells(mesh.cellIndexToLabel[i]);
            std::vector<int> indices; for(int label:cell.cellToVert) indices.push_back(mesh.verts(label).index());
            zone.connectivity.push_back(std::move(indices));
        }
        const auto field=[&](const char* name,const std::vector<double>& values) {
            TecplotWriter::Field f; f.name=name; f.location=TecplotWriter::Location::CellCentered;
            f.values=values; zone.fields.push_back(std::move(f));
        };
        field("potential",e.potential); field("epsilon",e.permittivity);
        field("Ex",e.ex); field("Ey",e.ey);
        TecplotWriter::Field phi; phi.name="levelSet";
        phi.values=e.levelSet; zone.fields.push_back(std::move(phi));
        if(charge) {
            TecplotWriter::Field q; q.name="interface_panel_mean_charge"; q.location=TecplotWriter::Location::CellCentered;
            q.values.assign(mesh.cellIndexNum(),0.);
            const auto density=charge->density();
            const auto& panels=charge->geometry->panels();
            for(std::size_t i=0;i<panels.size();++i) q.values[panels[i].cellIndex]=density[i];
            zone.fields.push_back(std::move(q));
            zone.auxiliary["ChargeMeaning"]="panel surface density on interface cells, zero elsewhere; not volume charge";
        }
        std::ostringstream file; file<<run_.outputPrefix<<"_step"<<std::setw(6)<<std::setfill('0')<<step<<".dat";
        TecplotWriter::write(std::filesystem::path("output/fields")/file.str(),zone);
    }
    history_->append(TecplotWriter::row(step,time,report.usedSuperLU,report.relativeResidual,
        bool(charge),total,bool(update),update?previous:time,
        update?update->current->netBoundaryCurrent:0.,
        report.conductionCurrent?report.conductionCurrent->netBoundaryCurrent:0.,
        update?update->electrodeResidual:0.,report.maxField,report.fieldEnergy));
}
