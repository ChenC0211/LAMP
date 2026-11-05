#pragma once

#include "Core/Base/Vector2D.h"
#include <array>
#include <cstdint>
#include <cstddef>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

// Serialization only: callers supply geometry, field locations, units and time.
// A valid file does not certify mesh orientation or physical correctness.
namespace TecplotWriter
{
enum class Location { Nodal, CellCentered };
enum class Element { Line, Quadrilateral, OrderedPoints, Triangle };

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
    // Zero-based input node indices; 2 for Line, 3 for Triangle, 4 for Quadrilateral.
    // OrderedPoints requires no FE connectivity and only nodal fields.
    std::vector<std::vector<int>> connectivity;
    std::vector<Field> fields;
    double solutionTime = 0.0;
    int strandId = 1;
    // Keys use [A-Za-z_][A-Za-z0-9_]*; values are quoted and escaped.
    std::map<std::string, std::string> auxiliary;
};

// Single-zone ASCII FE or ordered-points BLOCK file; all floating variables DOUBLE.
// Validates the entire zone before opening/truncating the requested path.
// Does not create directories. An I/O failure throws and may leave a partial
// file; callers decide publication/retention policy and choose unique paths.
void write(const std::filesystem::path& path, const Zone& zone);

// Shared VARIABLES followed by multiple zones. Coordinate names and field
// names/order must match; each zone may have its own field locations and topology.
// Validates every zone before opening/truncating the requested path.
void write(const std::filesystem::path& path, const std::vector<Zone>& zones);

// Read-only CSV companion export, intended after the owner flushes/closes its
// accepted/failure history. Numeric tokens retain their original precision.
// Each metric/category/finite continuous segment is one ordered POINT zone;
// missing values are omitted, never represented by zero or a nonfinite token.
// Returns false without touching datPath when no curve has a finite sample.
// Parses and validates before publication. An existing DAT may be replaced only
// when its Generator marker and canonical SourceCsv identify this exporter.
// Does not create directories; malformed input and I/O failures throw.
bool exportCsvHistory(const std::filesystem::path& csvPath,
    const std::filesystem::path& datPath);
// Recursively export *.csv below an existing series directory to sibling DATs.
// Symlinks are skipped; existing CSVs and History/field DAT semantics are intact.
// Returns the number of nonempty companion files published.
std::size_t exportCsvHistories(const std::filesystem::path& seriesDirectory);

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

// Explicit wide-table schema: no inference from an initial NA or category value.
enum class TableColumnKind { Numeric, Category };
struct TableColumn
{
    std::string name;
    TableColumnKind kind = TableColumnKind::Numeric;
};

// Buffered, append-only POINT table, created lazily at the first nonempty batch.
// append takes complete CSV records in memory (no header or CSV file). Numeric
// names/order and finite tokens are retained; each numeric column also has an
// appended <name>Available variable. Missing/nonfinite values are 0 placeholders
// with Available=0: callers must use those flags for XY value blanking.
// Categories are zone AUXDATA (readable display plus exact original byte hex).
// Category changes or a nonincreasing finite clock start a new zone; I is omitted.
// Each batch is fully validated before writing; a validation failure preserves
// the accepted prefix and permits another batch. An I/O failure poisons the
// writer and may leave a partial physical write. No existing file is replaced.
// flush/close expose I/O errors; destruction only releases the buffered stream.
// One owner only: do not modify the file while this writer is alive.
class TableHistory
{
public:
    TableHistory(const std::filesystem::path& path, const std::string& title,
        const std::vector<TableColumn>& columns,
        const std::map<std::string, std::string>& auxiliary = {});
    ~TableHistory();
    TableHistory(const TableHistory&) = delete;
    TableHistory& operator=(const TableHistory&) = delete;
    void append(const std::string& completeCsvRecords);
    void flush();
    void close();
private:
    struct Data;
    std::unique_ptr<Data> data_;
};

// Preserve evaluation in each argument's original type, then serialize without
// narrowing long-double conservation diagnostics or changing their arithmetic.
template<class... Values>
std::vector<long double> row(Values... values)
{
    return {static_cast<long double>(values)...};
}
}
