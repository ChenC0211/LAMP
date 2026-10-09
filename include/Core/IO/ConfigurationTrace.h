#pragma once

#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "Core/IO/StrictDictionary.h"

namespace DictionaryInput
{
using ConfigurationValue = std::variant<double, int, bool, std::string, std::vector<int>>;

struct ConfigurationEntry
{
    std::string field, file, rawValue, units, consumer;
    int line = 0;
    ConfigurationValue value;
    bool defaulted = false;
    bool used = true;
    std::string unusedReason;
};

struct ConfigurationTrace
{
    std::vector<ConfigurationEntry> entries;

    template<class Value> Value record(const Node& node, Value value,
        const std::string& units, const std::string& consumer)
    {
        entries.push_back({node.path, node.file, node.raw(), units, consumer,
            node.line, value, false, true, {}});
        return value;
    }
    template<class Value> Value defaultValue(const Node& block, const std::string& name,
        Value value, const std::string& units, const std::string& consumer)
    {
        entries.push_back({block.path + '.' + name, block.file, "<omitted>", units, consumer,
            block.line, value, true, true, {}});
        return value;
    }
};
}
