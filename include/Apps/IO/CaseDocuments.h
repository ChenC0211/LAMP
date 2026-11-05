#pragma once
#include "Core/IO/StrictDictionary.h"
#include "Core/Base/ScalarExpression.h"
#include <initializer_list>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace AppIO
{
// Provenance and effective text only; no physical state or configuration owner.
struct CaseDocument {std::string file,sourceText,effectiveText;};
// One caller-defined strict document set. Grammar remains owned by Core.
class CaseDocuments
{
public:
    CaseDocuments(const std::string& systemDir,std::initializer_list<const char*> files,
        const std::string& expectedApp,const std::string& contract="formal Hydro");
    const DictionaryInput::Node& root(const std::string& name) const;
    const std::string& path(const std::string& name) const;
    const std::string& sourceText(const std::string& name) const;
    const std::string& app() const noexcept{return app_;}
private:
    struct Source {std::string name,path,text;std::unique_ptr<DictionaryInput::Document> document;};
    const Source& source(const std::string& name) const;
    std::string app_;
    std::vector<Source> sources_;
};
// The caller validates usage only after all of its expressions are parsed.
struct ConstantsUsage
{
    void include(const std::string& formula);
    void validate(const DictionaryInput::Node& fields,const LAMP::ScalarExpression::Constants&,
        const char* expected) const;
private:
    std::set<std::string> identifiers_;
};
}
