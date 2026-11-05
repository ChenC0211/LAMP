#include "IO/CaseDocuments.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace AppIO
{
namespace
{
std::string readText(const std::string& path)
{
    std::ifstream in(path,std::ios::binary);
    if(!in)throw std::runtime_error("cannot open Hydro input: "+path);
    std::ostringstream out;out<<in.rdbuf();
    if(in.bad())throw std::runtime_error("cannot read Hydro input: "+path);
    return out.str();
}
}
CaseDocuments::CaseDocuments(const std::string& systemDir,std::initializer_list<const char*> files,
    const std::string& expectedApp,const std::string& contract):app_(expectedApp)
{
    namespace fs=std::filesystem;const fs::path directory(systemDir);
    if(!fs::is_directory(directory))throw std::invalid_argument("missing "+contract+" system directory: "+directory.string());
    std::set<std::string> names(files.begin(),files.end());
    if(names.size()!=files.size())throw std::logic_error("duplicate declared case document");
    for(const auto* name:files)if(!fs::is_regular_file(directory/name)||fs::is_symlink(directory/name))
        throw std::invalid_argument("missing or invalid "+contract+" input: "+(directory/name).string()+":1");
    std::string expected;for(const auto* name:files){if(!expected.empty())expected+='/';expected+=name;}
    for(const auto* name:files){
        const auto path=(directory/name).generic_string();const auto before=readText(path);
        auto document=std::make_unique<DictionaryInput::Document>(path);
        if(readText(path)!=before)throw std::runtime_error(path+":1: case source changed during Core dictionary parse");
        sources_.push_back({name,path,before,std::move(document)});
    }
    // Only an explicit local declaration extends the strict required set.
    // This asset is legacy Core grammar, not a StrictDictionary document.
    if((expectedApp=="Hydrodynamics"||expectedApp=="VofElectrohydrodynamics"||expectedApp=="Electrospray")&&names.count("geometry.dict")){
        const auto& initial=root("geometry.dict").block("Geometry").block("Initial");
        if(const auto* asset=initial.find("csgFile")){
            const auto name=asset->word();
            if(!std::regex_match(name,std::regex("[A-Za-z0-9][A-Za-z0-9_.-]*"))||name=="."||name==".."||names.count(name))
                asset->fail("a distinct simple local CSG basename without paths or traversal");
            const auto path=directory/name;
            if(!fs::is_regular_file(path)||fs::is_symlink(path))asset->fail("a regular nonsymlink declared CSG source file");
            const auto before=readText(path.generic_string());
            if(readText(path.generic_string())!=before)throw std::runtime_error(path.string()+":1: CSG source changed while reading");
            names.insert(name);expected+='/'+name;sources_.push_back({name,path.generic_string(),before,{}});
        }
    }
    for(const auto& entry:fs::directory_iterator(directory))
        if(!entry.is_regular_file()||entry.is_symlink()||!names.count(entry.path().filename().string()))
            throw std::invalid_argument("unused "+contract+" input: "+entry.path().string()+":1; expected only "+expected);
    const auto& dispatch=root("run.dict").block("Run");dispatch.allow({"app"});
    dispatch.get("app").choice({app_.c_str()});
}
const CaseDocuments::Source& CaseDocuments::source(const std::string& name) const
{
    const auto found=std::find_if(sources_.begin(),sources_.end(),[&](const Source& s){return s.name==name;});
    if(found==sources_.end())throw std::logic_error("case document was not declared: "+name);return *found;
}
const DictionaryInput::Node& CaseDocuments::root(const std::string& name) const
{const auto& entry=source(name);if(!entry.document)throw std::logic_error("raw CSG asset has no StrictDictionary root: "+name);return entry.document->root();}
const std::string& CaseDocuments::path(const std::string& name) const{return source(name).path;}
const std::string& CaseDocuments::sourceText(const std::string& name) const{return source(name).text;}
void ConstantsUsage::include(const std::string& formula)
{
    const std::regex token("(?:^|[^A-Za-z_0-9.])([A-Za-z_][A-Za-z_0-9]*)");
    for(auto it=std::sregex_iterator(formula.begin(),formula.end(),token);it!=std::sregex_iterator();++it)
        identifiers_.insert((*it)[1].str());
}
void ConstantsUsage::validate(const DictionaryInput::Node& fields,const LAMP::ScalarExpression::Constants& constants,
    const char* expected) const
{
    for(const auto& item:constants)if(!identifiers_.count(item.first))fields.block("Constants").get(item.first).fail(expected);
}
}
