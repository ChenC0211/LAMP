#pragma once

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "Core/IO/DictionaryTokenReader.h"
#include "Core/Units/UnitParser.h"

// Versioned physics dictionaries use one scalar entry per line or ';'. Blocks
// may open on the following line. Collections use named child blocks or '( )'
// lists. Geometry/mesh/PIC keep their existing readers and grammars.
namespace DictionaryInput
{
struct Node
{
    std::string file, path, name;
    int line = 1;
    bool isBlock = false;
    std::vector<std::string> values;
    std::vector<Node> entries;

    std::string raw() const
    {
        if (isBlock) return "{...}";
        std::string result;
        for (const auto& value : values) { if (!result.empty()) result += ' '; result += value; }
        return result;
    }
    [[noreturn]] void fail(const std::string& expected) const
    {
        throw std::runtime_error(file + ":" + std::to_string(line) + ": " + path
            + " = '" + raw() + "'; expected " + expected);
    }
    const Node* find(const std::string& key) const
    {
        for (const auto& entry : entries) if (entry.name == key) return &entry;
        return nullptr;
    }
    const Node& get(const std::string& key) const
    {
        if (!isBlock) fail("a block");
        if (const auto* entry = find(key)) return *entry;
        fail("required field/block " + (path.empty() ? key : path + '.' + key));
    }
    const Node& block(const std::string& key) const
    {
        const auto& entry = get(key);
        if (!entry.isBlock) entry.fail("a block");
        return entry;
    }
    void allow(std::initializer_list<const char*> keys) const
    {
        if (!isBlock) fail("a block");
        for (const auto& entry : entries)
            if (std::none_of(keys.begin(), keys.end(), [&](const char* key) { return entry.name == key; }))
                entry.fail("a supported field/block (unknown input is not ignored)");
    }
    std::string word() const
    {
        if (isBlock || values.size() != 1) fail("exactly one value");
        return values.front();
    }
    std::string choice(std::initializer_list<const char*> choices) const
    {
        const auto value = word();
        for (const char* choice : choices) if (value == choice) return value;
        std::string expected;
        for (const char* choice : choices) { if (!expected.empty()) expected += " | "; expected += choice; }
        fail(expected);
    }
    double number() const
    {
        const auto value = word();
        std::size_t used = 0;
        double result = 0.0;
        try { result = std::stod(value, &used); }
        catch (...) { fail("a finite number with no trailing text"); }
        if (used != value.size() || !std::isfinite(result)) fail("a finite number with no trailing text");
        return result;
    }
    int integer() const
    {
        const auto value = word();
        std::size_t used = 0;
        long long result = 0;
        try { result = std::stoll(value, &used); }
        catch (...) { fail("an integer with no trailing text"); }
        if (used != value.size() || result < std::numeric_limits<int>::min()
            || result > std::numeric_limits<int>::max()) fail("an integer with no trailing text");
        return static_cast<int>(result);
    }
    bool boolean() const { return choice({"true", "false"}) == "true"; }
    double quantity(UnitParser::Dimension dimension, const char* expectedUnit) const
    {
        if (isBlock || values.size() != 2) fail(std::string("number + unit, e.g. 1 ") + expectedUnit);
        UnitParser::Result result;
        if (!UnitParser::TryParseQuantityToSI(raw(), dimension, result))
            fail(std::string("a finite physical quantity with dimensions of ") + expectedUnit);
        return result.value_si;
    }
    std::vector<int> ids() const
    {
        if (isBlock || values.size() < 3 || values.front() != "(" || values.back() != ")")
            fail("a nonempty parenthesized list of positive 1-based IDs");
        std::vector<int> result;
        for (std::size_t i = 1; i + 1 < values.size(); ++i)
        {
            Node item = *this;
            item.values = {values[i]}; item.entries.clear();
            const int id = item.integer();
            if (id <= 0 || std::find(result.begin(), result.end(), id) != result.end())
                fail("distinct positive 1-based IDs");
            result.push_back(id);
        }
        return result;
    }
};

class Document
{
public:
    explicit Document(const std::string& file) : reader_(file, true)
    {
        root_.file = file; root_.isBlock = true;
        readEntries(root_, 0);
    }
    const Node& root() const noexcept { return root_; }

private:
    Reader reader_;
    Node root_;

    void readEntries(Node& parent, unsigned depth)
    {
        if (depth > 64) parent.fail("block nesting depth <= 64");
        while (!reader_.done())
        {
            if (reader_.peek() == ";") { reader_.next(); continue; }
            if (reader_.peek() == "}")
            {
                if (!depth) parent.fail("no unmatched closing brace");
                reader_.next(); return;
            }
            Node node;
            node.file = parent.file; node.line = reader_.line(); node.name = reader_.next();
            node.path = parent.path.empty() ? node.name : parent.path + '.' + node.name;
            if (node.name.empty() || node.name.find_first_of("{}();,") != std::string::npos)
                node.fail("a field or block name");
            if (parent.find(node.name)) node.fail("a unique field/block name (duplicate input)");
            if (reader_.peek() == "{")
            {
                reader_.next(); node.isBlock = true; readEntries(node, depth + 1);
            }
            else
            {
                int parentheses = 0;
                while (!reader_.done())
                {
                    const auto token = reader_.peek();
                    if (!parentheses && (token == ";" || token == "}" || reader_.line() != reader_.previousLine())) break;
                    if (token == "{" || token == "}" || token == ";") node.fail("a scalar/list followed by a line end or ';'");
                    node.values.push_back(reader_.next());
                    if (token == "(") ++parentheses;
                    if (token == ")" && --parentheses < 0) node.fail("balanced list parentheses");
                }
                if (parentheses || node.values.empty()) node.fail("a complete scalar value or list");
            }
            parent.entries.push_back(std::move(node));
        }
        if (depth) parent.fail("a closing brace");
    }
};
}
