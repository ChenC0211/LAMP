#pragma once

#include "Solvers/Poisson/CutCellMetrics.h"
#include <memory>

class DictionaryGeometry;
class UnitsContext;
class NormalizationContext;
class MovingEBController;

// Owns PIC's selected solid boundaries, composite level-set fields and cut-cell
// caches. Geometry/mesh/motion are borrowed synchronously for each operation.
// The single-level-set view borrows the geometry field: rebuild after geometry
// or topology changes before querying views. Composite views belong to this
// object. All views expire on rebuild, reselection or destruction.
class PicEmbeddedGeometry
{
public:
    enum class Mode { Explicit, Implicit };
    using Geometry = ImmersedGeometry<QuadtreeMesh>;
    using Metrics = Solvers::Poisson::CutCellMetrics;
    explicit PicEmbeddedGeometry(Mode mode);
    ~PicEmbeddedGeometry();
    PicEmbeddedGeometry(const PicEmbeddedGeometry&) = delete;
    PicEmbeddedGeometry& operator=(const PicEmbeddedGeometry&) = delete;
    static constexpr int alphaSampleCount = 6;

    void resolveSelection(const DictionaryGeometry& geometry, int levelSetCount);
    void rebuildComposite(const QuadtreeMesh& mesh, const Geometry* geometry);
    void rebuildMetrics(const QuadtreeMesh& mesh, Geometry* geometry, const FaceSegmentGraph& graph,
        const std::vector<int>& lineToBoundaryId, const MovingEBController* moving = nullptr);
    bool activeCell(const QuadtreeMesh& mesh, const Geometry* geometry, const Cell& cell) const;
    bool containsPoint(const QuadtreeMesh& mesh, const Geometry* geometry, const Vector2D& point, int hint) const;
    int boundaryIdAt(const QuadtreeMesh& mesh, Geometry* geometry,
        const std::vector<int>& lineToBoundaryId, const Vector2D& point, int hint) const;
    double surfaceAreaM2(const Geometry* geometry, const DictionaryGeometry& input,
        const std::vector<int>& lineToBoundaryId, const UnitsContext* units, const NormalizationContext* normalization) const;

    bool ready() const { return embeddedGeomReady; }
    bool absorbing() const { return embeddedAbsorb; }
    bool insideIsNegative() const { return embeddedInsideIsNegative; }
    int levelSetIndex() const { return embeddedLevelSetIndex; }
    const std::vector<int>& levelSetIndices() const { return embeddedLevelSetIndices; }
    const std::vector<int>& boundaryIds() const { return embeddedBoundaryIds; }
    bool usesComposite() const { return embeddedUseComposite; }
    const VertField<double>* levelSet() const { return embeddedLevelSetVert; }
    const std::vector<SignInt>& compositeLines() const { return embeddedCompositeLines; }
    bool metricsReady() const { return cutMetricsReady; }
    const Metrics& metrics() const { return cutMetrics; }
    const std::vector<double>& alphaByCell() const { return cutAlphaByCellIndex; }
    const std::vector<int>& boundaryIdsByCell() const { return embeddedBoundaryIdByCellIndex; }
private:
    void clearComposite();
    const Mode mode_;
    bool embeddedGeomReady = false, embeddedAbsorb = false, embeddedInsideIsNegative = true;
    int embeddedLevelSetIndex = -1;
    std::vector<int> embeddedLevelSetIndices, embeddedBoundaryIds;
    bool embeddedUseComposite = false;
    const VertField<double>* embeddedLevelSetVert = nullptr;
    std::unique_ptr<VertField<double>> embeddedPsiComposite, embeddedDxComposite, embeddedDyComposite;
    std::vector<SignInt> embeddedCompositeLines;
    Metrics cutMetrics;
    bool cutMetricsReady = false;
    std::vector<double> cutAlphaByCellIndex;
    std::vector<int> embeddedBoundaryIdByCellIndex;
};
