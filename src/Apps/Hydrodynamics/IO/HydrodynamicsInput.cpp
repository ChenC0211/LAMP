#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include "IO/TwoPhaseInputComponents.h"
#include <stdexcept>

LAMP::ScalarExpression::Jet HydrodynamicsInput::initialJet(Vector2D p) const
{
    if(!expression_)throw std::logic_error("Hydro initial data have not been validated");
    return expression_->sample(p.x,p.y);
}
HydrodynamicsInput readHydrodynamicsInput(const std::string& systemDir,std::optional<int> threadsOverride)
{
    const AppIO::CaseDocuments documents(systemDir,
        {"run.dict","geometry.dict","mesh.dict","interface.dict","materials.dict","hydro.dict"},
        "Hydrodynamics");
    const auto& method=documents.root("interface.dict").block("Interface").get("method");
    if(method.word()!="vof")method.fail("vof for standalone Hydrodynamics production");
    for(const auto& phase:documents.root("materials.dict").block("Materials").block("Phases").entries)
        phase.allow({"phaseId","density","dynamicViscosity"});
    HydrodynamicsInput input;AppIO::ConstantsUsage usage;
    const auto components=AppIO::readTwoPhaseInputComponents(input,documents,usage,threadsOverride);
    usage.validate(documents.root("geometry.dict").block("Geometry"),input.constants,"a constant used by Geometry.Initial.phi");
    AppIO::recordTwoPhaseInputDefaults(documents,input,input.trace);
    const std::array<const char*,6> files{"run.dict","geometry.dict","mesh.dict","interface.dict","materials.dict","hydro.dict"};
    for(std::size_t i=0;i<files.size();++i)
        input.documents.push_back({documents.path(files[i]),documents.sourceText(files[i]),components.effective[i]});
    return input;
}
