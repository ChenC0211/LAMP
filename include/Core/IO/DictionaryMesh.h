#pragma once

#include <iostream>
#include <unordered_map>
#include <vector>

#include "Core/IO/Dictionary.h"

class UnitsContext;

class DictionaryMesh : public Dictionary
{
public:
    Vector2D minCoord = Vector2D(0.0, 0.0);
    Vector2D maxCoord = Vector2D(1.0, 1.0);
    int Nx = 1;
    int Ny = 1;
    int levelSetCount = 1;

    // 自适应网格信息
    int maxCellNum = 0;

    int maxInitialIterations = 3;
    std::vector<int> minLevel_neg;
    std::vector<int> minLevel_pos;
    std::vector<int> maxLevel_local;

    double C_Lipschitz = 2.6;

    // ---------------------------------------------------------------------
    // Static refinement directives (demo support)
    // ---------------------------------------------------------------------
    struct RefineBox
    {
        Vector2D min;
        Vector2D max;
        int level = 0;
    };
    std::vector<RefineBox> refineBoxes;

    // ---------------------------------------------------------------------
    // V3 geometry/physics refinement settings
    // ---------------------------------------------------------------------
    struct GeometryRefinementConfig
    {
        enum class RegionShape
        {
            kUnknown = 0,
            kBox,
            kCircle,
            kPolygon
        };

        bool present = false;
        bool enabled = false;

        struct BoundaryRule
        {
            int boundaryId = -1;   // 1-based in input
            double cGamma = 1.0;   // C_Gamma coefficient
            int maxLevel = -1;     // max refinement level
        };

        struct PhaseRule
        {
            int phase = -1;        // 1-based in input
            int minLevel = -1;
            int maxLevel = -1;
        };

        struct RegionRule
        {
            std::string name;
            RegionShape shape = RegionShape::kUnknown;
            Vector2D min;
            Vector2D max;
            int maxLevel = -1;
        };

        std::vector<BoundaryRule> boundaryRules;
        std::vector<PhaseRule> phaseRules;
        std::vector<RegionRule> regionRules;
    };
    GeometryRefinementConfig geometryRefinement;

    struct PhysicsRefinementConfig
    {
        bool present = false;
        bool enabled = false;
        std::vector<std::string> indicators; // e.g. DebyeLength, ElectricField
        int adaptEverySteps = 0;
        int indicatorEverySteps = 0;
        int reportEvery = 0;
        double hysteresis = -1.0;
        int persistenceSteps = 0;
        std::unordered_map<std::string, int> maxLevelCaps; // indicator -> maxLevel
        // Explicit SI thresholds for current Hydro fields (not PIC references).
        double velocityThresholdSI = -1.0;
        double velocityGradientThresholdSI = -1.0;
    };
    PhysicsRefinementConfig physicsRefinement;

    struct GeometryRefreshConfig
    {
        bool present = false;
        // full | adaptive | skip
        std::string mode = "full";
        double changeRel = 0.2;
        bool verbose = false;
        bool allowWithMovingEb = false;
        // full | exact_refresh (experimental, env-gated in ImplicitPIC)
        std::string movingEbRefreshMode = "full";
    };
    GeometryRefreshConfig geometryRefresh;

    // Explicit time-loop policy, separate from initial mesh preparation and
    // physics indicators. The calling workflow must opt in as a consumer.
    struct RuntimeAdaptationConfig
    {
        bool present = false, enabled = false;
        int everySteps = 1;
        int maxSweeps = 2;
    };
    RuntimeAdaptationConfig runtimeAdaptation;

    DictionaryMesh();
    ~DictionaryMesh();

    // File/parse failures throw; the caller reports them. Parsing is in place:
    // discard the record after failure; no rollback of partially read data.
    void readDict(const UnitsContext* units = nullptr, bool acceptsRuntimeAdaptation = false);
    // Both input paths are explicit: the mesh reader also consumes Domain
    // from geometry. This does not change the process working directory.
    void readDict(const std::string& meshPath, const std::string& geometryPath,
        const UnitsContext* units = nullptr, bool acceptsRuntimeAdaptation = false);
    void ensureLevelSetCount(int n);

private:
    const std::string meshDictPath = "./system/mesh.dict";
    const std::string geometryDictPath = "./system/geometry.dict";
};
