#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "Core/IO/DictionaryBlockText.h"

namespace Solvers
{
namespace IO
{

// ----------------------------------------------------------------------------
// KVBlock
//
// A minimal key -> tokens map for a dictionary block.
// Keys are stored in lower-case.
// ----------------------------------------------------------------------------
struct KVBlock
{
    std::unordered_map<std::string, std::vector<std::string>> kv;

    bool has(const std::string& keyLower) const;
    const std::vector<std::string>* get(const std::string& keyLower) const;
};

// ----------------------------------------------------------------------------
// KVBlockParser
//
// A small, reusable helper for parsing QUANTUM-style dictionary files:
//   BlockName { key value ... }
//
// Responsibilities:
//   - Reuse Core's file loading and scoped block discovery.
//   - Parse key/value lines into a KVBlock.
//   - Provide typed getters with defaults.
//
// Notes:
//   - Solver key/value interpretation remains in Solvers/IO.
//   - It is intentionally conservative: it only understands one-line
//     key/value records and ignores unknown formats.
// ----------------------------------------------------------------------------
class KVBlockParser : public ::DictionaryBlockText
{
public:
    KVBlockParser() = default;
    ~KVBlockParser() override = default;

    // Parse key/value lines inside a block.
    KVBlock parseKeyValues(const std::string& blockText) const;

    // Typed getters.
    bool getBool(const KVBlock& block, const std::string& key, bool defaultVal) const;
    int getInt(const KVBlock& block, const std::string& key, int defaultVal) const;
    double getDouble(const KVBlock& block, const std::string& key, double defaultVal) const;
    std::string getString(const KVBlock& block, const std::string& key, const std::string& defaultVal) const;
    std::vector<std::string> getList(const KVBlock& block, const std::string& key) const;

    // Token helpers.
    static std::string toLowerCopy(std::string s);
    static bool parseBoolToken(const std::string& s, bool& out);

private:
    static std::string stripTrailingPunct(std::string s);
    static std::vector<std::string> stripParenTokens(const std::vector<std::string>& toks);
};

} // namespace IO
} // namespace Solvers
