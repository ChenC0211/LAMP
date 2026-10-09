#include "Core/IO/TecplotWriter.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <stdexcept>

namespace TecplotWriter
{
namespace
{
void require(bool condition, const char* reason)
{
    if (!condition) throw std::invalid_argument(std::string("TecplotWriter: ") + reason);
}

std::string quoted(const std::string& input)
{
    std::string result = "\"";
    for (unsigned char c : input)
    {
        require(c >= 32 && c != 127, "strings cannot contain control characters");
        if (c == '"' || c == '\\') result += '\\';
        result += static_cast<char>(c);
    }
    result += '"';
    // Each string is emitted on a single line; quoted strings cannot continue.
    require(result.size() <= 31900, "quoted string exceeds the ASCII line limit");
    return result;
}

bool identifier(const std::string& name)
{
    const auto alpha = [](unsigned char c) {
        return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '_';
    };
    if (name.empty() || !alpha(static_cast<unsigned char>(name.front()))) return false;
    for (unsigned char c : name)
        if (!alpha(c) && !(c >= '0' && c <= '9')) return false;
    return true;
}

struct Header
{
    std::string title;
    std::vector<std::string> variables;
    std::vector<std::size_t> cellVariables;
    std::vector<std::string> auxiliary;
    const char* zoneType = nullptr;
};

Header validate(const Zone& zone)
{
    Header header;
    header.title = quoted(zone.title);
    require(zone.title.size() <= 128, "zone title exceeds 128 bytes");
    require(std::isfinite(zone.solutionTime), "solution time must be finite");
    require(zone.strandId > 0, "strand ID must be positive");
    const auto maxInt = static_cast<std::size_t>(std::numeric_limits<int>::max());
    require(!zone.points.empty() && zone.points.size() <= maxInt, "invalid node count");
    require(!zone.connectivity.empty() && zone.connectivity.size() <= maxInt, "invalid element count");
    require(zone.fields.size() <= maxInt - 2, "too many variables");
    std::size_t nodesPerElement = 0;
    switch (zone.element)
    {
    case Element::Line: nodesPerElement = 2; header.zoneType = "FELINESEG"; break;
    case Element::Quadrilateral: nodesPerElement = 4; header.zoneType = "FEQUADRILATERAL"; break;
    default: require(false, "unsupported element type");
    }
    for (const auto& point : zone.points)
        require(std::isfinite(point.x) && std::isfinite(point.y), "nonfinite coordinate");
    for (const auto& element : zone.connectivity)
    {
        require(element.size() == nodesPerElement, "incorrect connectivity size");
        for (std::size_t k = 0; k < element.size(); ++k)
        {
            require(element[k] >= 0 && static_cast<std::size_t>(element[k]) < zone.points.size(),
                "connectivity node is out of range");
            for (std::size_t j = 0; j < k; ++j)
                require(element[j] != element[k], "element repeats a node");
        }
    }
    std::set<std::string> names;
    const auto addName = [&](const std::string& name) {
        require(!name.empty(), "variable name cannot be empty");
        require(names.insert(name).second, "variable names must be unique, including coordinates");
        header.variables.push_back(quoted(name));
    };
    for (const auto& name : zone.coordinateNames) addName(name);
    for (const auto& field : zone.fields)
    {
        addName(field.name);
        std::size_t count = 0;
        switch (field.location)
        {
        case Location::Nodal: count = zone.points.size(); break;
        case Location::CellCentered:
            count = zone.connectivity.size();
            header.cellVariables.push_back(header.variables.size());
            break;
        default: require(false, "unsupported field location");
        }
        require(field.values.size() == count, "field value count does not match its location");
        for (double value : field.values) require(std::isfinite(value), "nonfinite field value");
    }
    for (const auto& entry : zone.auxiliary)
    {
        require(identifier(entry.first), "invalid AUXDATA identifier");
        const std::string line = "AUXDATA " + entry.first + "=" + quoted(entry.second);
        require(line.size() <= 32000, "AUXDATA exceeds the ASCII line limit");
        header.auxiliary.push_back(line);
    }
    return header;
}

template<class Value>
void block(std::ostream& out, std::size_t count, Value value)
{
    for (std::size_t k = 0; k < count; ++k)
        out << value(k) << ((k + 1) % 6 == 0 || k + 1 == count ? '\n' : ' ');
}
}

History::History(const std::filesystem::path& path, const std::string& title,
    const std::vector<std::string>& columns,
    const std::map<std::string, std::string>& auxiliary)
    : path_(path), columns_(columns.size())
{
    require(!path.empty() && !std::filesystem::exists(path), "history requires a new output path");
    require(!columns.empty() && title.size() <= 128, "invalid history columns or title");
    std::ostringstream header;
    header.imbue(std::locale::classic());
    header << "VARIABLES=\n";
    std::set<std::string> unique;
    for (const auto& column : columns)
    {
        require(!column.empty() && unique.insert(column).second, "empty or duplicate history column");
        header << quoted(column) << '\n';
    }
    header << "ZONE T=" << quoted(title) << ", DATAPACKING=POINT\nDT=(\n";
    for (std::size_t k = 0; k < columns_; ++k) header << "DOUBLE\n";
    header << ")\n";
    for (const auto& entry : auxiliary)
    {
        require(identifier(entry.first), "invalid history AUXDATA identifier");
        const std::string line = "AUXDATA " + entry.first + "=" + quoted(entry.second);
        require(line.size() <= 32000, "history AUXDATA exceeds the ASCII line limit");
        header << line << '\n';
    }
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(path_, std::ios::binary);
    out << header.str();
    out.close();
    bytes_ = std::filesystem::file_size(path_);
}

void History::append(const std::vector<long double>& values)
{
    require(!failed_, "cannot append after history I/O failure");
    require(values.size() == columns_, "history row has incorrect column count");
    std::ostringstream record;
    record.imbue(std::locale::classic());
    record << std::setprecision(17) << std::defaultfloat;
    for (std::size_t k = 0; k < values.size(); ++k)
    {
        require(std::isfinite(values[k]), "history values must be finite; use an explicit validity flag");
        record << (k ? " " : "") << values[k];
    }
    record << '\n';
    try
    {
        require(std::filesystem::file_size(path_) == bytes_, "history file was modified outside its owner");
        std::ofstream out;
        out.exceptions(std::ios::failbit | std::ios::badbit);
        out.open(path_, std::ios::binary | std::ios::app);
        out << record.str();
        out.close();
        bytes_ = std::filesystem::file_size(path_);
    }
    catch (const std::exception& error)
    {
        failed_ = true;
        throw std::runtime_error("TecplotWriter: history append failed for " + path_.string() + ": " + error.what());
    }
}

void write(const std::filesystem::path& path, const Zone& zone)
{
    const Header header = validate(zone);
    require(!path.empty(), "output path cannot be empty");
    std::ofstream out;
    out.imbue(std::locale::classic());
    out.exceptions(std::ios::failbit | std::ios::badbit);
    try
    {
        out.open(path, std::ios::out | std::ios::trunc);
        out << std::setprecision(17) << std::defaultfloat;
        out << "VARIABLES=\n";
        for (const auto& name : header.variables) out << name << '\n';
        out << "ZONE T=" << header.title << '\n'
            << "NODES=" << zone.points.size() << ", ELEMENTS=" << zone.connectivity.size()
            << ", ZONETYPE=" << header.zoneType << ", DATAPACKING=BLOCK\n"
            << "SOLUTIONTIME=" << zone.solutionTime << ", STRANDID=" << zone.strandId << '\n';
        out << "DT=(\n";
        for (std::size_t k = 0; k < header.variables.size(); ++k) out << "DOUBLE\n";
        out << ")\n";
        if (!header.cellVariables.empty())
        {
            out << "VARLOCATION=([\n";
            for (std::size_t k = 0; k < header.cellVariables.size(); ++k)
                out << (k ? "," : "") << header.cellVariables[k] << '\n';
            out << "]=CELLCENTERED)\n";
        }
        for (const auto& entry : header.auxiliary) out << entry << '\n';
        block(out, zone.points.size(), [&](std::size_t k) { return zone.points[k].x; });
        block(out, zone.points.size(), [&](std::size_t k) { return zone.points[k].y; });
        for (const auto& field : zone.fields)
            block(out, field.values.size(), [&](std::size_t k) { return field.values[k]; });
        for (const auto& element : zone.connectivity)
        {
            for (std::size_t k = 0; k < element.size(); ++k)
                out << element[k] + 1 << (k + 1 == element.size() ? '\n' : ' ');
        }
        out.close();
    }
    catch (const std::ios_base::failure& error)
    {
        throw std::runtime_error("TecplotWriter: I/O failure for " + path.string() + ": " + error.what());
    }
}
}
