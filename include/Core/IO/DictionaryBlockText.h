#pragma once

#include "Core/IO/Dictionary.h"

// Shared lexical block discovery. Callers own block names, source precedence,
// parsing of values and any materialized files. No solver or physics policy.
class DictionaryBlockText : public Dictionary
{
public:
    enum class Scope { AnyDepth, TopLevel };

    // Preserve the existing line-comment rules of Dictionary.
    bool loadFile(const std::string& path, std::string& outText);

    // Match an unquoted block name at a record boundary followed by its brace.
    // Quoted braces and comments do not affect nesting. Return body text only.
    bool extractBlockText(const std::string& text, const std::string& name,
        std::string& body, Scope scope = Scope::AnyDepth) const;
    // Keep alias priority in caller order, then the first matching occurrence.
    bool extractFirstBlockText(const std::string& text, const std::vector<std::string>& names,
        std::string& body, Scope scope = Scope::AnyDepth) const;
};
