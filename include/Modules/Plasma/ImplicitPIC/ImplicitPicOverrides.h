#pragma once
#include <optional>
#include <string>

// Typed legacy overrides supplied by the caller, independent of process state.
// An absent value preserves the dictionary/default policy at its consumer.
struct ImplicitPicOverrides
{
    struct Polarization
    {
        std::optional<std::string> mode;
        std::optional<bool> temporalGuess, reuseIterativePattern;
        std::optional<double> temporalOmega;
    } polarization;
    bool segmenterTrace = false;
    bool segmenterIgnoreHostHints = false;
    std::optional<std::string> initGeomSnapshotPath;
    std::optional<bool> seeDebug;
    std::optional<int> packedFlushEvery;
    std::optional<double> amrGateEFieldRel;
    std::optional<double> amrGateDebyeRel;
    std::optional<int> amrMaxChanges;
    std::optional<bool> initGeomSnapshot;
    std::optional<bool> periodicDiagDebug;
    std::optional<double> consoleFallbackFracWarn;
    std::optional<int> consoleFallbackDebounceSteps;
    std::optional<bool> consoleFallbackPrintRecover;
    std::optional<double> consoleDqrelWarn;
    std::optional<double> consoleDqrelRecover;
    std::optional<int> consoleDqrelDebounceSteps;
    std::optional<int> consoleDqrelEventCooldownSteps;
    std::optional<bool> consoleDqrelPrintRecover;
    std::optional<double> consoleRgL2Warn;
    std::optional<int> consoleAmrSkipStreakWarn;
    std::optional<bool> consolePhysicsEnabled;
    std::optional<int> consolePhysicsEvery;
    std::optional<bool> consolePhysicsShowCharge;
    std::optional<bool> consolePhysicsShowCurrent;
    std::optional<bool> consolePhysicsShowWall;
    std::optional<bool> consolePhysicsShowMcc;
    std::optional<bool> consolePhysicsShowPerf;
    std::optional<int> depositFillOmpMinParticles;
    std::optional<bool> orbitAvgFused;
    std::optional<bool> resetTransientsOnAmr;
    std::optional<int> cutcellStateDumpStep;
    std::optional<std::string> cutcellStateDumpPath;
    std::optional<bool> cutcellStateDumpOnlyActive;
    std::optional<int> movingTraceStep;
    std::optional<bool> movingTraceAll;
    std::optional<std::string> movingTracePath;
    std::optional<bool> outputFastDensity;
    std::optional<bool> outputCellConstant;
    std::optional<bool> outputDisableCellConstant;
    std::optional<bool> outputMaskDiag;
    std::optional<bool> packedBatchWithAmr;
    std::optional<bool> amrFastRemap;
    std::optional<std::string> amrGeomRefreshMode;
    std::optional<double> amrGeomRefreshChangeRel;
    std::optional<bool> amrGeomRefreshVerbose;
    std::optional<bool> amrGeomRefreshAllowWithMovingEb;
    std::optional<bool> amrAllowExactMovingEbRefresh;
    std::optional<std::string> amrMovingEbRefreshMode;
    std::optional<bool> amrReinitActiveMovingOnly;
    std::optional<bool> amrReinitHybridSensitiveInactive;
    std::optional<bool> amrReinitHybridUseDelta;
    std::optional<int> amrReinitHybridNearZeroMin;
    std::optional<int> amrReinitHybridNarrowBandMin;
    std::optional<int> amrReinitHybridNearZeroDeltaMin;
    std::optional<int> amrReinitHybridNarrowBandDeltaMin;
    std::optional<bool> amrReinitHybridIncludeMaxIter;
    std::optional<int> amrReinitHybridMaxIterThreshold;
    std::optional<int> amrReinitActiveMovingMaxStreak;
    std::optional<bool> orbitAvgSpeciesPar;
    std::optional<int> segmenterMaxSegments;
    std::optional<int> segmenterMaxLevel;
};
