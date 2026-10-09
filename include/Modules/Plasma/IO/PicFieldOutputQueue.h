#pragma once

#include "Core/Field/CellField.h"
#include "Core/Field/VertField.h"
#include "Modules/Plasma/IO/PicFieldSelection.h"
#include "Modules/Plasma/Physics/PicCurrentWeights.h"
#include <cstddef>
#include <string>
#include <vector>

class QuadtreeMesh;
class UnitsContext;
class NormalizationContext;
struct ImplicitPicOverrides;
class FaceSegmentGraph;
namespace Solvers { namespace Poisson { class CutCellMetrics; } }
template<class Mesh> class ImmersedGeometry;

// Owns implicit-PIC output policy, coordinates and prepared frames. Physical
// inputs are borrowed synchronously while preparing a frame; only the mesh and
// immutable unit context must outlive this object.
// The caller must flush before changing that geometry and at normal completion.
// Destruction drops pending frames; it never writes or retries failed IO.
class PicFieldOutputQueue
{
public:
    struct Options {
        int flushEvery = 1;
        bool fastDensity = true, forceCells = false, disableCells = true;
        bool maskDiagnostics = false, allowBatchingWithAmr = false;
        static Options fromImplicit(const ImplicitPicOverrides& overrides);
    };
    struct Input {
        const FaceSegmentGraph& graph;
        const PicEmbeddedGeometry& embedded;
        const ImmersedGeometry<QuadtreeMesh>* geometry;
        const std::vector<std::string>& levelSetNames;
        const Solvers::Poisson::CutCellMetrics* metrics;
        bool useCutMetrics;
        const CellField<double>& potential;
        const CellField<double>& density;
        const CellField<Vector2D>& electric;
        const std::vector<CellField<double>>& numberAmounts;
        const std::vector<ParticleType>& species;
        const std::vector<double>& weights;
        const NormalizationContext* normalization;
        const PicCurrentWeights::Policy& coupling;
        const std::vector<double>& muCell;
        const std::vector<double>& kappaCell;
        const std::vector<double>& kappaSegment;
        bool periodicX, periodicY;
    };
    struct Timing { double compute_s = 0.0, io_s = 0.0; };

    PicFieldOutputQueue(const QuadtreeMesh& mesh, const UnitsContext* units,
        const std::vector<std::string>& fields, Options options);
    // Already-prepared frames may supply their own coordinate snapshot.
    PicFieldOutputQueue(const QuadtreeMesh& mesh, const std::vector<double>& x,
        const std::vector<double>& y, int flushEveryOutputs);
    PicFieldOutputQueue(const PicFieldOutputQueue&) = delete;
    PicFieldOutputQueue& operator=(const PicFieldOutputQueue&) = delete;

    // Called at the original stage where AMR availability is resolved.
    void configureBatching(bool amrEnabled);
    bool batchingEnabled() const { return batchingEnabled_; }
    // Prepare the selected positive completed step. Run cadence is supplied by
    // PicOutputCadence before this call; empty selection/null units still skip IO.
    Timing writeImplicit(int completedStep, const Input& input);
    // Caller flushes before AMR, then refreshes at the original post-rebuild
    // callback. Wrong mesh or pending frames reject before changing coordinates.
    void refreshCoordinates(const QuadtreeMesh& mesh);

    // Returns only actual write time. Queued frames take ownership of names/
    // fields and copy the mask; the immediate path writes the supplied frame.
    double submitVertices(const std::string& name, std::vector<std::string>&& names,
        std::vector<VertField<double>>&& fields, const std::vector<unsigned char>& mask, bool batch);
    double submitCells(const std::string& name, std::vector<std::string>&& names,
        std::vector<CellField<double>>&& fields, const std::vector<unsigned char>& mask, bool batch);
    double flush();
    std::size_t pending() const noexcept { return jobs_.size(); }

private:
    struct Job
    {
        std::string nameNoExt;
        std::vector<std::string> names;
        std::vector<VertField<double>> fields;
        std::vector<CellField<double>> cellFields;
        std::vector<unsigned char> cellMask;
        bool cellConstant = false;
    };
    template<class Field>
    double submit(const std::string& name, std::vector<std::string>&& names,
        std::vector<Field>&& fields, const std::vector<unsigned char>& mask, bool batch);

    const QuadtreeMesh& mesh_;
    std::vector<double> x_, y_;
    const UnitsContext* units_;
    const PicFieldSelection selection_;
    const Options options_;
    const int flushEveryOutputs_;
    bool batchingEnabled_ = false;
    std::vector<Job> jobs_;
};
