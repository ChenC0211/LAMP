#pragma once

#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Shared dictionary lexer extracted from the strict Hydro reader. Existing
// Hydro callers retain token, numeric and error behavior. Versioned documents
// retain semicolons and use token lines to delimit scalar entries.
namespace DictionaryInput
{
class Reader
{
public:
    explicit Reader(const std::string& path, bool retainSemicolons = false) : path_(path)
    {
        std::ifstream input(path);
        if (!input) fail("cannot open dictionary");
        std::string line;
        int lineNumber = 0;
        while (std::getline(input, line))
        {
            ++lineNumber;
            for (std::size_t i = 0; i < line.size();)
            {
                const char c = line[i];
                if (std::isspace(static_cast<unsigned char>(c))) { ++i; continue; }
                if (c == '#' || (c == '/' && i + 1 < line.size() && line[i + 1] == '/')) break;
                if (c == ';')
                {
                    if (retainSemicolons) tokens_.push_back({";", lineNumber});
                    ++i; continue;
                }
                if (c == '"' || c == '\'')
                {
                    const auto end = line.find(c, i + 1);
                    if (end == std::string::npos) failAt(lineNumber, "unterminated quoted string");
                    tokens_.push_back({line.substr(i + 1, end - i - 1), lineNumber});
                    i = end + 1;
                    continue;
                }
                if (isPunctuation(c))
                {
                    tokens_.push_back({std::string(1, c), lineNumber});
                    ++i;
                    continue;
                }
                const auto start = i;
                while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i])) &&
                       !isPunctuation(line[i]) && line[i] != ';' && line[i] != '#' &&
                       !(line[i] == '/' && i + 1 < line.size() && line[i + 1] == '/')) ++i;
                tokens_.push_back({line.substr(start, i - start), lineNumber});
            }
        }
    }

    const std::string& path() const noexcept { return path_; }
    int line() const { return done() ? (tokens_.empty() ? 1 : tokens_.back().line) : tokens_[position_].line; }
    int previousLine() const { return position_ ? tokens_[position_ - 1].line : 1; }

    bool done() const { return position_ == tokens_.size(); }
    std::string peek() const { return done() ? std::string() : tokens_[position_].text; }
    std::string next()
    {
        if (done()) fail("unexpected end of dictionary");
        return tokens_[position_++].text;
    }
    void expect(const std::string& expected)
    {
        if (done() || next() != expected) fail("expected '" + expected + "'");
    }
    double number(const std::string& key)
    {
        const std::string token = next();
        std::size_t used = 0;
        double value = 0.0;
        try { value = std::stod(token, &used); }
        catch (...) { fail(key + " must be a finite numeric value in code units"); }
        if (used != token.size() || !std::isfinite(value))
            fail(key + " must be a finite numeric value in code units");
        return value;
    }
    int integer(const std::string& key)
    {
        const std::string token = next();
        std::size_t used = 0;
        long long value = 0;
        try { value = std::stoll(token, &used); }
        catch (...) { fail(key + " must be an integer"); }
        if (used != token.size() || value < std::numeric_limits<int>::min() ||
            value > std::numeric_limits<int>::max()) fail(key + " must be an integer");
        return static_cast<int>(value);
    }
    bool boolean(const std::string& key)
    {
        const std::string token = next();
        if (token == "true") return true;
        if (token == "false") return false;
        fail(key + " must be true or false");
    }
    [[noreturn]] void fail(const std::string& reason) const
    {
        const int line = done() ? (tokens_.empty() ? 1 : tokens_.back().line) : tokens_[position_].line;
        failAt(line, reason);
    }

private:
    struct Token { std::string text; int line; };
    std::string path_;
    std::vector<Token> tokens_;
    std::size_t position_ = 0;
    static bool isPunctuation(char c)
    {
        return c == '{' || c == '}' || c == '(' || c == ')' || c == ',';
    }
    [[noreturn]] void failAt(int line, const std::string& reason) const
    {
        throw std::runtime_error(path_ + ":" + std::to_string(line) + ": " + reason);
    }
};
}
