#pragma once

#include <string>
#include <functional>

#include "Core/IO/Dictionary.h"
#include "Modules/Plasma/ExplicitPIC/PicConfig.h"

class DictionaryPIC : public Dictionary
{
public:
    DictionaryPIC() = default;
    ~DictionaryPIC() = default;

    // Optional legacy input records owned by the caller. Invoked synchronously
    // with the existing brace/token interpretation; never retained. Return true
    // only for a consumed record. Ordinary PIC inputs need no additional reader.
    using AdditionalRecord = std::function<bool(const std::string& block, const std::string& parent,
        const std::string& line, const std::vector<std::string>& tokens)>;

    // Return false if the file cannot be opened.
    // If units!=nullptr, any value+unit scalars (e.g. "1 ns", "500 V") will be
    // converted to code values using the UnitsContext reference scales.
    bool readDict(const std::string& picDictPath = "./system/pic.dict", const UnitsContext* units = nullptr,
                  const AdditionalRecord& additionalRecord = {});

    const PicConfig& config() const { return config_; }

private:
    PicConfig config_{};

    static std::string toLowerCopy(std::string s);
    static bool parseBoolToken(const std::string& s, bool& out);
    static PicBCType parseBCTypeToken(const std::string& s, bool& ok);
};
