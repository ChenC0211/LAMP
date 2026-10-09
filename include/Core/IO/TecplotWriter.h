#pragma once

#include "Core/Base/Vector2D.h"
#include <array>
#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

// Serialization only: callers supply geometry, field locations, units and time.
// A valid file does not certify mesh orientation or physical correctness.
namespace TecplotWriter
{
enum class Location { Nodal, CellCentered };
enum class Element { Line, Quadrilateral };

struct Field
{
    std::string name;
    Location location = Location::Nodal;
    std::vector<double> values;
};

struct Zone
{
    std::string title;
    std::array<std::string, 2> coordinateNames{"x", "y"};
    std::vector<Vector2D> points;
    Element element = Element::Quadrilateral;
    // Zero-based input node indices; 2 for Line, 4 for Quadrilateral.
    std::vector<std::vector<int>> connectivity;
    std::vector<Field> fields;
    double solutionTime = 0.0;
    int strandId = 1;
    // Keys use [A-Za-z_][A-Za-z0-9_]*; values are quoted and escaped.
    std::map<std::string, std::string> auxiliary;
};

// Single-zone ASCII FE/BLOCK file; all floating variables explicitly DOUBLE.
// Validates the entire zone before opening/truncating the requested path.
// Does not create directories. An I/O failure throws and may leave a partial
// file; callers decide publication/retention policy and choose unique paths.
void write(const std::filesystem::path& path, const Zone& zone);

// Ordered POINT history for XY Line. Each append is a complete numeric row,
// flushed and closed before returning. I is deliberately omitted: Tecplot
// infers the current number of rows, including while the simulation is running.
// No time integration or resampling; column units and constant metadata belong
// to the caller. Undefined quantities require an explicit caller-owned flag.
class History
{
public:
    History(const std::filesystem::path& path, const std::string& title,
        const std::vector<std::string>& columns,
        const std::map<std::string, std::string>& auxiliary = {});
    History(const History&) = delete;
    History& operator=(const History&) = delete;
    void append(const std::vector<long double>& values);
private:
    std::filesystem::path path_;
    std::size_t columns_;
    std::uintmax_t bytes_ = 0;
    bool failed_ = false;
};

// Preserve evaluation in each argument's original type, then serialize without
// narrowing long-double conservation diagnostics or changing their arithmetic.
template<class... Values>
std::vector<long double> row(Values... values)
{
    return {static_cast<long double>(values)...};
}
}
