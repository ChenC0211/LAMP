#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Core/Base/Enumeration.h"
#include "Core/Base/Structures.h"
#include "Core/Base/Vector2D.h"
#include "Core/Mesh/Cell.h"
#include "Core/Mesh/Edge.h"
#include "Core/Mesh/Vert.h"

namespace MeshCacheIO
{

// Opt-in byte budget for a seekable cache whose file size is already known.
// Bounds are checked before allocations, not only after a truncated read.
class ReadBudget
{
public:
    ReadBudget(std::istream& input, std::uint64_t bytes) : input_(input), remaining_(bytes)
    {
        previous_ = input_.pword(slot());
        input_.pword(slot()) = this;
    }
    ~ReadBudget() { input_.pword(slot()) = previous_; }
    ReadBudget(const ReadBudget&) = delete;
    ReadBudget& operator=(const ReadBudget&) = delete;
    std::uint64_t remaining() const { return remaining_; }
    static void requireBytes(std::istream& input, std::uint64_t bytes)
    {
        auto* budget = static_cast<ReadBudget*>(input.pword(slot()));
        if (bytes > static_cast<std::uint64_t>((std::numeric_limits<std::streamsize>::max)()) ||
            (budget && bytes > budget->remaining_))
            throw std::runtime_error("mesh cache is truncated or declares a size exceeding remaining file bytes");
    }
    static void consumed(std::istream& input, std::uint64_t bytes)
    {
        auto* budget = static_cast<ReadBudget*>(input.pword(slot()));
        if (budget) budget->remaining_ -= bytes;
    }
private:
    static int slot() { static const int index = std::ios_base::xalloc(); return index; }
    std::istream& input_;
    std::uint64_t remaining_;
    void* previous_ = nullptr;
};

inline void readBytes(std::istream& input, char* destination, std::uint64_t bytes)
{
    ReadBudget::requireBytes(input, bytes);
    if (bytes != 0 && !input.read(destination, static_cast<std::streamsize>(bytes)))
        throw std::runtime_error("mesh cache read failed or file is truncated");
    ReadBudget::consumed(input, bytes);
}

template<typename T>
inline void writeRaw(std::ostream& os, const T& value)
{
    os.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

inline void writeDiscCoeff(std::ostream& os, const DiscCoeff& value)
{
    static_assert(std::is_standard_layout<DiscCoeff>::value &&
                  std::is_trivially_copyable<DiscCoeff>::value, "native mesh coefficient layout required");
    // QMC2/BMC2 retain their native-layout record size and member offsets.
    // Serialize members into zeroed storage rather than leaking C++ padding.
    // Even a zero-count stencil may carry a meaningful RHS: preserve all bits.
    std::array<char, sizeof(DiscCoeff)> encoded{};
    std::memcpy(encoded.data() + offsetof(DiscCoeff, discNum), &value.discNum, sizeof(value.discNum));
    std::memcpy(encoded.data() + offsetof(DiscCoeff, discOpr), value.discOpr, sizeof(value.discOpr));
    std::memcpy(encoded.data() + offsetof(DiscCoeff, rhs), &value.rhs, sizeof(value.rhs));
    os.write(encoded.data(), static_cast<std::streamsize>(encoded.size()));
}

template<typename T>
inline void readRaw(std::istream& is, T& value)
{
    if constexpr (std::is_same<T, bool>::value)
    {
        static_assert(sizeof(bool) == 1, "QMC2 stores native one-byte booleans");
        unsigned char encoded = 0;
        readBytes(is, reinterpret_cast<char*>(&encoded), 1);
        if (encoded > 1) throw std::runtime_error("mesh cache contains an invalid boolean");
        value = encoded != 0;
    }
    else readBytes(is, reinterpret_cast<char*>(&value), sizeof(T));
}

template<typename T>
inline void writeRawArray(std::ostream& os, const T* data, std::size_t count)
{
    if (count == 0)
        return;
    os.write(reinterpret_cast<const char*>(data), sizeof(T) * count);
}

template<typename T>
inline void readRawArray(std::istream& is, T* data, std::size_t count)
{
    if (count == 0)
        return;
    if (count > (std::numeric_limits<std::uint64_t>::max)() / sizeof(T))
        throw std::runtime_error("mesh cache array byte size overflows");
    readBytes(is, reinterpret_cast<char*>(data), static_cast<std::uint64_t>(sizeof(T)) * count);
}

template<typename T>
inline void writeVector(std::ostream& os, const std::vector<T>& values)
{
    const int n = static_cast<int>(values.size());
    writeRaw(os, n);
    if (n > 0)
        writeRawArray(os, values.data(), static_cast<std::size_t>(n));
}

template<typename T>
inline void readVector(std::istream& is, std::vector<T>& values)
{
    int n = 0;
    readRaw(is, n);
    if (n < 0) throw std::runtime_error("mesh cache contains a negative vector size");
    ReadBudget::requireBytes(is, static_cast<std::uint64_t>(n) * sizeof(T));
    values.resize(static_cast<std::size_t>(n));
    if (n > 0)
        readRawArray(is, values.data(), static_cast<std::size_t>(n));
}

inline void writeVert(std::ostream& os, const Vert& v, int levelSetCount)
{
    const std::size_t n = static_cast<std::size_t>(std::max(1, levelSetCount));

    const int label = v.label();
    const int index = v.index();
    const Vector2D coord = v.coord();
    const bool needUpdate = v.needUpdate();
    const Boundary boundary = v.boundary();
    const TVertType type = v.type();
    const int phase = v.phase();

    writeRaw(os, label);
    writeRaw(os, index);
    writeRaw(os, coord);
    writeRaw(os, needUpdate);
    writeRaw(os, boundary);
    writeRaw(os, type);
    writeRaw(os, phase);

    writeRawArray(os, v.vertToHEdge, 2);
    writeRawArray(os, v.vertToVEdge, 2);
    writeRawArray(os, v.vertToCell, 4);
    writeRawArray(os, v.neiVertices, 4);
    writeRawArray(os, v.supVert, 2);
    writeRaw(os, v.supCell);
    writeRawArray(os, v.son, 1);
    writeRaw(os, v.parent);

    writeRaw(os, v.inPhaseNum);
    writeRaw(os, v.relatedVertsNum);
    writeRawArray(os, v.relatedVerts, 6);

    writeDiscCoeff(os, v.lapCoeff);
    writeDiscCoeff(os, v.selfCoeff);
    writeDiscCoeff(os, v.nCoeff);
    writeDiscCoeff(os, v.tCoeff);
    writeDiscCoeff(os, v.xCoeff);
    writeDiscCoeff(os, v.yCoeff);
    writeDiscCoeff(os, v.nnCoeff);
    writeDiscCoeff(os, v.ttCoeff);

    writeRawArray(os, v.onCutEdge.data(), n);
    writeRawArray(os, v.onCutCell.data(), n);
    writeRawArray(os, v.nearestNodeCoord.data(), n);
    writeRawArray(os, v.nearestNodeLabel.data(), n);
}

inline void readVert(std::istream& is, Vert& v, int levelSetCount)
{
    const std::size_t n = static_cast<std::size_t>(std::max(1, levelSetCount));

    int label = -1;
    int index = -1;
    Vector2D coord;
    bool needUpdate = false;
    Boundary boundary = Boundary::non;
    TVertType type = TVertType::non;
    int phase = -1;

    readRaw(is, label);
    readRaw(is, index);
    readRaw(is, coord);
    readRaw(is, needUpdate);
    readRaw(is, boundary);
    readRaw(is, type);
    readRaw(is, phase);

    v.setLabel(label);
    v.setIndex(index);
    v.setCoord(coord);
    v.setNeedUpdate(needUpdate);
    v.setBoundary(boundary);
    v.setType(type);
    v.setPhase(phase);

    readRawArray(is, v.vertToHEdge, 2);
    readRawArray(is, v.vertToVEdge, 2);
    readRawArray(is, v.vertToCell, 4);
    readRawArray(is, v.neiVertices, 4);
    readRawArray(is, v.supVert, 2);
    readRaw(is, v.supCell);
    readRawArray(is, v.son, 1);
    readRaw(is, v.parent);

    readRaw(is, v.inPhaseNum);
    readRaw(is, v.relatedVertsNum);
    readRawArray(is, v.relatedVerts, 6);

    readRaw(is, v.lapCoeff);
    readRaw(is, v.selfCoeff);
    readRaw(is, v.nCoeff);
    readRaw(is, v.tCoeff);
    readRaw(is, v.xCoeff);
    readRaw(is, v.yCoeff);
    readRaw(is, v.nnCoeff);
    readRaw(is, v.ttCoeff);

    v.ensureLevelSetSize(n);
    readRawArray(is, v.onCutEdge.data(), n);
    readRawArray(is, v.onCutCell.data(), n);
    readRawArray(is, v.nearestNodeCoord.data(), n);
    readRawArray(is, v.nearestNodeLabel.data(), n);
}

inline void writeEdge(std::ostream& os, const Edge& e, int levelSetCount)
{
    const std::size_t n = static_cast<std::size_t>(std::max(1, levelSetCount));

    const int label = e.label();
    const int index = e.index();
    const Vector2D midpoint = e.midpointCoord();
    const double length = e.length();
    const bool isHorizontal = e.isHorizontal();
    const int preIndex = e.preIndex();
    const Boundary boundary = e.boundary();
    const bool needUpdate = e.needUpdate();
    const bool cut = e.cut();
    const double cutLength = e.cutLength();

    writeRawArray(os, e.isCut.data(), n);
    writeRawArray(os, e.edgeToVert, 2);
    writeRaw(os, e.masterCell);
    writeRaw(os, e.slaveCell);
    writeRaw(os, e.parent);
    writeRawArray(os, e.son, 2);
    writeVector(os, e.relatedPreIndex);
    writeVector(os, e.relatedPreCoord);

    writeRaw(os, label);
    writeRaw(os, index);
    writeRaw(os, midpoint);
    writeRaw(os, length);
    writeRaw(os, isHorizontal);
    writeRaw(os, preIndex);
    writeRaw(os, boundary);
    writeRaw(os, needUpdate);
    writeRaw(os, cut);
    writeRaw(os, cutLength);
}

inline void readEdge(std::istream& is, Edge& e, int levelSetCount)
{
    const std::size_t n = static_cast<std::size_t>(std::max(1, levelSetCount));

    e.ensureLevelSetSize(n);
    readRawArray(is, e.isCut.data(), n);
    readRawArray(is, e.edgeToVert, 2);
    readRaw(is, e.masterCell);
    readRaw(is, e.slaveCell);
    readRaw(is, e.parent);
    readRawArray(is, e.son, 2);
    readVector(is, e.relatedPreIndex);
    readVector(is, e.relatedPreCoord);

    int label = -1;
    int index = -1;
    Vector2D midpoint;
    double length = 0.0;
    bool isHorizontal = true;
    int preIndex = -1;
    Boundary boundary = Boundary::non;
    bool needUpdate = false;
    bool cut = false;
    double cutLength = 0.0;

    readRaw(is, label);
    readRaw(is, index);
    readRaw(is, midpoint);
    readRaw(is, length);
    readRaw(is, isHorizontal);
    readRaw(is, preIndex);
    readRaw(is, boundary);
    readRaw(is, needUpdate);
    readRaw(is, cut);
    readRaw(is, cutLength);

    e.setLabel(label);
    e.setIndex(index);
    e.setMidpointCoord(midpoint);
    e.setLength(length);
    e.setIsHorizontal(isHorizontal);
    e.setPreIndex(preIndex);
    e.setBoundary(boundary);
    e.setNeedUpdate(needUpdate);
    e.setCut(cut);
    e.setCutLength(cutLength);
}

inline void writeCell(std::ostream& os, const Cell& c, int levelSetCount)
{
    const std::size_t n = static_cast<std::size_t>(std::max(1, levelSetCount));

    const int label = c.label();
    const int index = c.index();
    const int preIndex = c.preIndex();
    const Vector2D center = c.centerCoord();
    const int level = c.level();
    const Vector2D size = c.size();
    const bool needUpdate = c.needUpdate();
    const Boundary boundary = c.boundary();

    writeRawArray(os, c.isCut.data(), n);
    writeRawArray(os, c.cellToVert, 4);
    writeRawArray(os, &c.cellToHEdge[0][0], 4);
    writeRawArray(os, &c.cellToVEdge[0][0], 4);
    writeRawArray(os, &c.neiCell[0][0], 8);
    writeVector(os, c.relatedPreIndex);
    writeVector(os, c.relatedPreCoord);
    writeRawArray(os, c.son, 4);
    writeRaw(os, c.parent);
    writeRaw(os, c.supCellIndex);
    writeRaw(os, c.coarsenFlag);
    writeRaw(os, c.refineFlag);
    writeRaw(os, c.needCoarsen);
    writeRaw(os, c.needRefine);

    writeRaw(os, label);
    writeRaw(os, index);
    writeRaw(os, preIndex);
    writeRaw(os, center);
    writeRaw(os, level);
    writeRaw(os, size);
    writeRaw(os, needUpdate);
    writeRaw(os, boundary);
}

inline void readCell(std::istream& is, Cell& c, int levelSetCount)
{
    const std::size_t n = static_cast<std::size_t>(std::max(1, levelSetCount));

    c.ensureLevelSetSize(n);
    readRawArray(is, c.isCut.data(), n);
    readRawArray(is, c.cellToVert, 4);
    readRawArray(is, &c.cellToHEdge[0][0], 4);
    readRawArray(is, &c.cellToVEdge[0][0], 4);
    readRawArray(is, &c.neiCell[0][0], 8);
    readVector(is, c.relatedPreIndex);
    readVector(is, c.relatedPreCoord);
    readRawArray(is, c.son, 4);
    readRaw(is, c.parent);
    readRaw(is, c.supCellIndex);
    readRaw(is, c.coarsenFlag);
    readRaw(is, c.refineFlag);
    readRaw(is, c.needCoarsen);
    readRaw(is, c.needRefine);

    int label = -1;
    int index = -1;
    int preIndex = -1;
    Vector2D center;
    int level = 0;
    Vector2D size;
    bool needUpdate = false;
    Boundary boundary = Boundary::non;

    readRaw(is, label);
    readRaw(is, index);
    readRaw(is, preIndex);
    readRaw(is, center);
    readRaw(is, level);
    readRaw(is, size);
    readRaw(is, needUpdate);
    readRaw(is, boundary);

    c.setLabel(label);
    c.setIndex(index);
    c.setPreIndex(preIndex);
    c.setCenterCoord(center);
    c.setLevel(level);
    c.setSize(size);
    c.setNeedUpdate(needUpdate);
    c.setBoundary(boundary);
}

inline bool magicEquals(const char a[4], const char b[4])
{
    return std::equal(a, a + 4, b);
}

} // namespace MeshCacheIO
