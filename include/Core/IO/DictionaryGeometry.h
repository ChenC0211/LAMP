#pragma once

#include <iostream>
#include <unordered_map>
#include "Core/IO/Dictionary.h"
#include "Core/Geometry/Line.h"
#include "Core/Base/Vector2D.h"
#include "Core/Base/Enumeration.h"
#include "Core/Base/Structures.h"

class DictionaryGeometry : public Dictionary
{
public:
    Vector2D domainMin = Vector2D(0.0, 0.0);
    Vector2D domainMax = Vector2D(1.0, 1.0);
    bool hasDomain = false;

    // Optional domain dimension/closure metadata.
    // dimension: 2 (2D slab) or 3 (fully 3D).
    int dimension = 2;
    bool hasDimension = false;
    // Out-of-plane thickness for 2D slab closure (SI meters).
    double slabThickness_m = 1.0;
    bool hasSlabThickness = false;

    std::vector<Vector2D> points;
    std::vector<Line> lines;

    // Level set info
    std::vector<std::string> levelSetName;
    std::vector<Vector2D> innerPoints;
    std::vector<bool> moving;
    std::vector<std::vector<SignInt>> linesOnInterface;
    // Optional directed 1-based Line IDs in the dictionary, stored 0-based.
    // Appended only to the initialization ray-parity path, never the actual
    // distance/cut/force interface. Empty entries retain the legacy closure.
    std::vector<std::vector<SignInt>> signClosureLines;

    // Phase info
    std::vector<std::string> phaseName;
    std::vector<std::vector<LevelSetOperation>> levelSetOperation;
    std::vector<std::vector<SignInt>> levelSetIndexA;
    std::vector<std::vector<SignInt>> levelSetIndexB;
    std::vector<std::vector<SignInt>> boundaryLines;

    struct BoundaryEntity
    {
        int id = -1; // 1-based ID from input
        std::string name;
        GeometricBoundary kind = GeometricBoundary::domain;
        std::vector<int> lineIds; // 0-based line indices
        int levelSetIndex = -1;   // 0-based
        int phase = -1;           // 0-based
        int opposedPhase = -1;    // 0-based
        int adjacentPhase = -1;   // 0-based (domainBoundary)
        bool immersed = true;
        bool ignored = false;
    };
    std::vector<BoundaryEntity> boundaries;
    std::unordered_map<int, BoundaryEntity> boundaryById;

    std::vector<int> immersedBoundaryIds() const;
    bool inferLevelSetInsideIsNegative(int levelSetIndex, bool& insideIsNegative) const;
    bool inferBoundaryInsideIsNegative(int boundaryId, bool& insideIsNegative) const;
    bool inferEmbeddedInsideIsNegative(bool& insideIsNegative) const;

    // Periodic boundary pairing (optional).
    // Each pair connects two domain boundaries and specifies a translation shift that maps
    // points on boundaryA to their periodic counterparts on boundaryB.
    //
    // NOTE: In PB0 we only parse and store these; subsequent PB steps consume them in
    // Poisson / E-field / PIC deposition and particle wrap.
    struct PeriodicPair
    {
        int boundaryA = -1;
        int boundaryB = -1;
        Vector2D shiftAtoB = Vector2D(0.0, 0.0);
        std::string mode = "strict"; // strict|loose (future)
    };

    std::vector<PeriodicPair> periodicPairs;
    // Lookup tables: boundaryId -> mate boundaryId / shift to mate (A->B uses +shift, B->A uses -shift).
    std::unordered_map<int, int> periodicMateByBoundary;
    std::unordered_map<int, Vector2D> periodicShiftToMate;

    inline bool isPeriodicBoundary(const int boundaryId) const
    {
        return periodicMateByBoundary.find(boundaryId) != periodicMateByBoundary.end();
    }

    inline int periodicMate(const int boundaryId) const
    {
        auto it = periodicMateByBoundary.find(boundaryId);
        return (it == periodicMateByBoundary.end()) ? -1 : it->second;
    }

    inline Vector2D periodicShift(const int boundaryId) const
    {
        auto it = periodicShiftToMate.find(boundaryId);
        return (it == periodicShiftToMate.end()) ? Vector2D(0.0, 0.0) : it->second;
    }

    DictionaryGeometry();
    ~DictionaryGeometry();

    // File/parse failures throw; the caller reports them. Parsing is in place:
    // discard the record after failure; no rollback of partially read data.
    void readDict(const UnitsContext* units = nullptr);
    // Explicit source path; does not change the process working directory.
    void readDict(const std::string& path, const UnitsContext* units = nullptr);

    // Uninterpreted extension sections belong to their consuming modules.
    const std::string* extensionBlock(const std::string& name) const
    {
        const auto found = extensionBlocks_.find(name);
        return found == extensionBlocks_.end() ? nullptr : &found->second;
    }
    bool extensionEnabled(const std::string& name, bool defaultWhenPresent = true) const;

private:
    std::unordered_map<std::string, std::string> extensionBlocks_;
    void captureExtension(std::istream& file, const std::string& name);
    const std::string geometryDictPath = "./system/geometry.dict";
};
