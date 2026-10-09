#pragma once

#include <initializer_list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

enum class PicFieldMode { Explicit, Implicit };

// One sequence per species-field group. Retains the legacy per-base suffix
// counter; it does not resolve collisions with an independently suffixed base.
class PicFieldKeys
{
public:
    PicFieldKeys(PicFieldMode mode, std::size_t expected);
    std::string next(const std::string& rawName);
private:
    PicFieldMode mode_;
    std::unordered_map<std::string, int> counts_;
};

// Output metadata only: no dictionary reader, fields, physical state or IO.
// Explicit defaults/mandatory densities and fixed order differ from implicit
// explicit selection and input-token ordering; retain those policies here.
class PicFieldSelection
{
public:
    struct Fields
    {
        bool phi=false, rho=false, psi=false, psiComposite=false, psiAll=false;
        bool ex=false, ey=false, n=false, te=false, meanEnergy=false;
        bool muCell=false, muFace=false;
        bool weightMean=false, weightSum=false, weightMax=false, macroCount=false;
    };

    PicFieldSelection(PicFieldMode mode, const std::vector<std::string>& tokens);
    PicFieldMode mode() const { return mode_; }
    const Fields& fields() const { return fields_; }
    bool configured() const { return !orderedTokens_.empty(); }
    std::string key(const std::string& rawName) const;
    PicFieldKeys keys(std::size_t expected) const { return {mode_, expected}; }

    // Only implicit output is reordered. Mismatched name/value counts retain
    // the old no-op behavior; the serializer still owns its own validation.
    template<class Field>
    void order(std::vector<std::string>& names, std::vector<Field>& fields) const
    {
        if (mode_ != PicFieldMode::Implicit || orderedTokens_.empty() || names.size() != fields.size())
            return;
        const auto perm = permutation(names);
        std::vector<std::string> namesOrdered;
        namesOrdered.reserve(names.size());
        std::vector<Field> fieldsOrdered;
        fieldsOrdered.reserve(fields.size());
        for (const std::size_t index : perm)
        {
            namesOrdered.push_back(std::move(names[index]));
            fieldsOrdered.push_back(std::move(fields[index]));
        }
        names = std::move(namesOrdered);
        fields = std::move(fieldsOrdered);
    }

private:
    using Aliases = std::initializer_list<const char*>;
    bool has(Aliases aliases) const;
    int first(Aliases aliases) const;
    int priority(const std::string& outputName) const;
    std::vector<std::size_t> permutation(const std::vector<std::string>& names) const;
    PicFieldMode mode_;
    Fields fields_;
    std::unordered_set<std::string> tokens_;
    std::vector<std::string> orderedTokens_;
};
