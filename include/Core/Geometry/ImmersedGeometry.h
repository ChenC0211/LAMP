#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <queue>
#include <unordered_set>
#include <optional>
#include <array>
#include <functional>
#include "Core/Base/Base.h"
#include "Core/Base/Vector2D.h"
#include "Core/Geometry/Line.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/Geometry/BSpline.h"
#include "Core/Geometry/LevelSet.h"
#include "Core/Geometry/Phase.h"
#include "Core/Geometry/LevelSetBoundaryConditions.h"
#include "Core/Base/Utility.h"
#include "Core/Discrete/FiniteDifference.h"

namespace ImmersedGeometryDetail { struct InitialCurve; }

template<typename Mesh>
class ImmersedGeometry
{
public:
    struct TraceDiagnostics
    {
        int clampCount = 0;
        int fallbackCellCount = 0;
        int longTraceCount = 0;
        int interpolationFallbackCount = 0;
        int nonFiniteInterpolationCount = 0;
        double maxTraceRatio = 0.0;
        int inflowBoundaryCount = 0;
        int prescribedBoundaryVertexCount = 0;
    };

    struct LevelSetUpdateDiagnostics
    {
        double minValue = 0.0;
        double maxValue = 0.0;
        int nearZeroCount = 0;
        int nonFiniteCount = 0;
        int narrowBandCount = 0;
        int reinitIterations = 0;
        double reinitFineResidual = 0.0;
        double reinitCoarseResidual = 0.0;
        bool reinitialized = false;
    };

    struct CutInfoDiagnostics
    {
        int cutCellCount = 0;
        int cutHEdgeCount = 0;
        int cutVEdgeCount = 0;
    };

    std::vector<Vector2D> points;
    std::vector<Line> lines;

    // 水平集信息(当前)
    std::vector<LevelSet> levelSets;
    // 水平集信息(上一步)，随网格更新插值更新
    std::vector<VertField<double>> levelSetsPre;
    std::vector<VertField<double>> levelSetsPre_dx;
    std::vector<VertField<double>> levelSetsPre_dy;
    // 水平集信息(上上一步)，不需随网格插值更新，利用上一步赋值即可
    std::vector<VertField<double>> levelSetsPrePre;
    std::vector<VertField<double>> levelSetsPrePre_dx;
    std::vector<VertField<double>> levelSetsPrePre_dy;

    // 相信息
    std::vector<Phase> phases;

    ImmersedGeometry(DictionaryGeometry& dictGeometry, Mesh& mesh);
    ~ImmersedGeometry();

    std::vector<Vector2D> discreteLines(const int levelSetI, const int localLineI, int numSamples = 50);
    // 初始化水平集
    // Diagnostic scatter is an explicit caller choice; it never controls
    // geometry calculations. Default preserves existing non-spray callers.
    void initializeLevelSets(Mesh& mesh, double time = 0.0, bool writeInitialScatter = true);
    // Distance to ONLY the selected physical lines, using the same exact
    // straight/arc and sampled-spline primitives as level-set initialization.
    // No sign closure, zero displacement, phase marking or history mutation.
    // An open interface needs separate active-domain/sign/boundary data before
    // this unsigned field can be used as a transported scalar.
    // Optional endpoint tangents point INTO each oriented initial primitive;
    // spline tangents belong to the same sampled polyline, not another curve.
    void initializeUnsignedDistanceToLines(const Mesh& mesh,
        const std::vector<SignInt>& selectedLines, VertField<double>& distance,
        std::vector<std::array<Vector2D, 2>>* inwardEndpointTangents = nullptr);
    struct InitialCurveSample { double unsignedDistance; bool negative; };
    // Owned t=0 primitive sampler for local reconstruction queries, not an
    // interpolated evolving field. Distance uses only distanceLines; sign uses
    // the existing complete level-set path/closure and innerPoint convention.
    // No legacy zero displacement. A physical solid wall still needs the
    // caller's fluid-side phase trace; parity is not a wall boundary condition.
    std::function<InitialCurveSample(Vector2D)> initialCurveSampler(const Mesh& mesh,
        const std::vector<SignInt>& distanceLines, int signLevelSet);
    // Opt-in, fixed-uniform nonperiodic scalar boundaries. Configure once per
    // level set before initialization/advancement. No field is changed here.
    // Mesh changes, legacy cut geometry/cache and reinitialization are rejected
    // while this mode is active; the shared Hermite geometry consumes phi.
    void setLevelSetBoundaryConditions(const Mesh& mesh, int levelSetI,
        LevelSetBoundaryConditions conditions);
    // Scalar derivatives/history storage only; does not construct legacy cut
    // lengths or change phi. Safe with declared exact-zero boundary vertices.
    void refreshLevelSetDerivatives(const Mesh& mesh);

    // 相标记
    void markPhases(Mesh& mesh);

    // 重初始化水平集
    void reinitializeLevelSets(const Mesh& mesh, const double& t);
    void pdeReinitializeMovingLevelSets(const Mesh& mesh, const double& t);
    void pdeReinitializeMovingLevelSets(const Mesh& mesh, const double& t, const std::vector<int>& activeLevelSetIndices);
    void exactRedistanceMovingLevelSets(const Mesh& mesh);
    void exactRedistanceMovingLevelSets(const Mesh& mesh, const std::vector<int>& activeLevelSetIndices);
    void refreshLevelSetScalarDiagnostics(const Mesh& mesh);
    void refreshLevelSetDerivedData(const Mesh& mesh, const double& t, const bool reinitializeMovingLevelSets);
    void dumpAllLevelSetsScatter(const Mesh& mesh, const std::string& tag);

    // 更新水平集函数

    // previousTime is required only when explicit scalar boundaries are active.
    // A boundary event samples phi at previousTime+dt-elapsed, not at t_old.
    // Optional velocity callbacks require explicit scalar boundaries and
    // replace only velocity interpolation, not RK2/time/history handling.
    // A retained old scalar interpolant can likewise replace the default
    // point interpolation, without modifying any nodal source or boundary data.
    void updateLevelSets(const Mesh& mesh, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, bool radialAxisSymmetry = false,
        double previousTime = std::numeric_limits<double>::quiet_NaN(),
        const LevelSetAdvectionVelocity* advectionVelocity = nullptr,
        const LevelSetAdvectionScalar* advectionScalar = nullptr);
    Vector2D getPreCoord(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, bool radialAxisSymmetry = false);
    // Forward material point obtained by inverting the SAME numerical
    // backtrace as updateLevelSets. Vx/Vy are new, preVx/preVy are old.
    // No geometry/diagnostic state is written. Requires finite dt >= 0,
    // interior unclipped trajectories, sampled speed*dt and displacement
    // <= half the local minimum cell size. Throws on invalid input or if
    // y <- y + coord - backtrace(y) does not converge within 32 evaluations.
    // Residual tolerance: max(128*epsilon*coordinateScale, 1e-12*cellSize).
    // Explicit RZ mode requires min.r=0 and permits the exact axis in strict
    // traces. Axis radial velocity must be exactly zero; negative-r stages
    // are rejected, never clipped/reflected. Other strict walls stay interior.
    Vector2D getNextCoord(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, bool radialAxisSymmetry = false) const;
    Vector2D getPrePreCoord(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, const double& preDt);

    void checkLevelSet(const Mesh& mesh, const int levelSetI);
    void checkCutInfo(Mesh& mesh, const bool isInit);

    inline VertField<double>& LS_dx(const int i) { return LSs_dx_[i]; }
    inline VertField<double>& LS_dy(const int i) { return LSs_dy_[i]; }
    inline const VertField<double>& LS_dx(const int i) const { return LSs_dx_[i]; }
    inline const VertField<double>& LS_dy(const int i) const { return LSs_dy_[i]; }
    inline Field<double>& LS_dxx(const int i) { return LSs_dxx_[i]; }
    inline Field<double>& LS_dyy(const int i) { return LSs_dyy_[i]; }
    inline Field<double>& LS_dxy(const int i) { return LSs_dxy_[i]; }
    inline const Field<double>& LS_dxx(const int i) const { return LSs_dxx_[i]; }
    inline const Field<double>& LS_dyy(const int i) const { return LSs_dyy_[i]; }
    inline const Field<double>& LS_dxy(const int i) const { return LSs_dxy_[i]; }
    inline Field<double>& LS_dx_backward(const int i) { return LSs_dx_backward_[i]; }
    inline Field<double>& LS_dy_backward(const int i) { return LSs_dy_backward_[i]; }
    inline Field<double>& LS_dx_forward(const int i) { return LSs_dx_forward_[i]; }
    inline Field<double>& LS_dy_forward(const int i) { return LSs_dy_forward_[i]; }
    inline const TraceDiagnostics& traceDiagnostics() const { return traceDiagnostics_; }
    inline const std::vector<LevelSetUpdateDiagnostics>& levelSetUpdateDiagnostics() const { return levelSetUpdateDiagnostics_; }
    inline const CutInfoDiagnostics& cutInfoDiagnostics() const { return cutInfoDiagnostics_; }

    void wirteCache(const std::string& filePath) const;
    void readCache(Mesh& mesh, const std::string& filePath);

private:
    std::vector<VertField<double>> LSs_dx_;
    std::vector<VertField<double>> LSs_dy_;
    std::vector<Field<double>> LSs_dxx_;
    std::vector<Field<double>> LSs_dyy_;
    std::vector<Field<double>> LSs_dxy_;

    std::vector<Field<double>> LSs_dx_backward_;
    std::vector<Field<double>> LSs_dy_backward_;
    std::vector<Field<double>> LSs_dx_forward_;
    std::vector<Field<double>> LSs_dy_forward_;
    std::vector<std::vector<Vector2D>> periodicDistanceShiftsByLevelSet_;
    // Owned, sign-only initialization closures. These references are never
    // inserted into LevelSet::lines(), so distance/cut/forces remain physical.
    std::vector<std::vector<SignInt>> signClosureLines_;
    ImmersedGeometryDetail::InitialCurve makeInitialCurve(SignInt reference, double hMin);
    std::vector<Vector2D> discreteLine(SignInt reference, int numSamples);
    TraceDiagnostics traceDiagnostics_;
    std::vector<LevelSetUpdateDiagnostics> levelSetUpdateDiagnostics_;
    CutInfoDiagnostics cutInfoDiagnostics_;
    std::vector<LevelSetBoundaryConditions> levelSetBoundaryConditions_;
    struct BoundaryMeshLayout
    {
        Vector2D lower, upper, spacing;
        int nx = 0, ny = 0;
        std::vector<Vector2D> vertices; // Exact coordinates in solver-index order.
        struct CellLayout
        {
            int label;
            Vector2D center;
            std::array<int, 4> vertices;
        };
        std::vector<CellLayout> cells;
    };
    std::optional<BoundaryMeshLayout> boundaryMeshLayout_;
    bool periodicGeometry_ = false;
    bool boundaryHistoryInitialized_ = false;

    struct TraceQueryDiagnostics
    {
        int clampCount = 0;
        int fallbackCellCount = 0;
        int longTraceCount = 0;
        double traceRatio = 0.0;
        int inflowBoundaryCount = 0;
        int prescribedBoundaryVertexCount = 0;
    };

    bool hasLevelSetBoundaryConditions() const;
    void validateBoundaryMesh(const Mesh& mesh) const;
    bool prescribedBoundaryVertex(const Mesh& mesh, int levelSetI, Vector2D point) const;
    std::optional<double> prescribedBoundaryValue(const Mesh& mesh, int levelSetI,
        Vector2D point, double time, std::optional<Boundary> side = std::nullopt) const;
    void stageBoundaryLevelSet(const Mesh& mesh, int levelSetI,
        const VertField<double>& Vx, const VertField<double>& Vy,
        const VertField<double>& preVx, const VertField<double>& preVy,
        double dt, bool radialAxisSymmetry, double previousTime,
        std::vector<double>& values, std::vector<TraceQueryDiagnostics>& diagnostics,
        const LevelSetAdvectionVelocity* advectionVelocity,
        const LevelSetAdvectionScalar* advectionScalar) const;

    void resetUpdateDiagnostics();
    void syncHistoricalLevelSetStorageToCurrent(const Mesh& mesh, const int levelSetI);
    void summarizeLevelSetDiagnostics(const Mesh& mesh, const int levelSetI);
    void prepareLevelSetDerivedStorage(const Mesh& mesh, const int levelSetI);
    void computeSecondDerivatives(const Mesh& mesh, const int levelSetI);
    void computeFirstDerivativesAndMixed(const Mesh& mesh, const int levelSetI);
    Vector2D getPreCoordDetailed(const Mesh& mesh, const Vector2D& coord, const VertField<double>& Vx, const VertField<double>& Vy, const VertField<double>& preVx, const VertField<double>& preVy, const double& dt, TraceQueryDiagnostics* diag, bool strict = false, bool radialAxisSymmetry = false) const;

    void diffForReinitLS(const Mesh& mesh, const int levelSetI);

    void searchCutEdges(const Mesh& mesh, const int levelSetI);
    void computeCutLength(const Mesh& mesh, const int levelSetI);

    void reinitLevelSet(const Mesh& mesh, const int levelSetI, const double& t);
    void exactRedistanceLevelSet(const Mesh& mesh, const int levelSetI);
    template<typename Field1, typename Field2>
    void updateLevelSet(const Mesh& mesh, const int levelSetI, const Field1& Vx, const Field2& Vy, const Field1& preVx, const Field2& preVy, const double& dt, bool radialAxisSymmetry = false,
        const std::vector<double>* stagedValues = nullptr,
        const std::vector<TraceQueryDiagnostics>* stagedDiagnostics = nullptr);

};

# include "Core/Geometry/ImmersedGeometry_impl.h"
