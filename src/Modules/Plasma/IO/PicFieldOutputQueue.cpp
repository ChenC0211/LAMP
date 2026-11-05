#include "Modules/Plasma/IO/PicFieldOutputQueue.h"
#include "Modules/Plasma/IO/PicFieldOutput.h"
#include "Modules/Plasma/IO/PicFieldFrame.h"
#include "Modules/Plasma/ImplicitPIC/ImplicitPicOverrides.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Units/UnitsContext.h"
#include <algorithm>
#include <chrono>
#include <type_traits>
#include <utility>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace {
using Clock = std::chrono::steady_clock;
double secondsSince(const Clock::time_point& start)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - start).count();
}
}

PicFieldOutputQueue::Options PicFieldOutputQueue::Options::fromImplicit(const ImplicitPicOverrides& overrides)
{
    return {std::max(1, overrides.packedFlushEvery.value_or(1)),
        overrides.outputFastDensity.value_or(true), overrides.outputCellConstant.value_or(false),
        overrides.outputDisableCellConstant.value_or(true), overrides.outputMaskDiag.value_or(false),
        overrides.packedBatchWithAmr.value_or(false)};
}

PicFieldOutputQueue::PicFieldOutputQueue(const QuadtreeMesh& mesh, const UnitsContext* units,
    const std::vector<std::string>& fields, Options options)
    : mesh_(mesh), units_(units), selection_(PicFieldMode::Implicit, fields), options_(options),
      flushEveryOutputs_(std::max(1, options.flushEvery))
{
    refreshCoordinates(mesh);
    if (flushEveryOutputs_ > 1)
        jobs_.reserve(static_cast<std::size_t>(flushEveryOutputs_));
}

PicFieldOutputQueue::PicFieldOutputQueue(const QuadtreeMesh& mesh, const std::vector<double>& x,
    const std::vector<double>& y, int flushEveryOutputs)
    : PicFieldOutputQueue(mesh, nullptr, {}, Options{flushEveryOutputs})
{
    x_ = x;
    y_ = y;
}

void PicFieldOutputQueue::refreshCoordinates(const QuadtreeMesh& mesh)
{
    if (&mesh != &mesh_) throw std::invalid_argument("PIC output coordinates belong to another mesh");
    if (!jobs_.empty()) throw std::logic_error("flush PIC output before changing coordinates");
    if (!units_) return;
    x_.assign(static_cast<std::size_t>(mesh.vertIndexNum()), 0.0);
    y_.assign(static_cast<std::size_t>(mesh.vertIndexNum()), 0.0);
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const auto& v = mesh.verts(mesh.vertIndexToLabel[i]);
        x_[static_cast<std::size_t>(i)] = units_->LengthCodeToSI(v.coord().x);
        y_[static_cast<std::size_t>(i)] = units_->LengthCodeToSI(v.coord().y);
    }
}

void PicFieldOutputQueue::configureBatching(bool amrEnabled)
{
    batchingEnabled_ = (flushEveryOutputs_ > 1) && (!amrEnabled || options_.allowBatchingWithAmr);
    if (options_.fastDensity)
        std::cout << "[ImplicitPIC][output] fast density mode enabled (LAMP_IMPLICIT_OUTPUT_FAST_DENSITY=1)" << std::endl;
    if (flushEveryOutputs_ > 1)
        std::cout << "[ImplicitPIC][output] packed-field flush every outputs=" << flushEveryOutputs_
                  << " (batching=" << (batchingEnabled_ ? "enabled" : "disabled")
                  << ", amr=" << (amrEnabled ? "on" : "off") << ")" << std::endl;
    if ((flushEveryOutputs_ > 1) && !batchingEnabled_)
        std::cout << "[ImplicitPIC][output] packed-field batching disabled (AMR active)" << std::endl;
    else if ((flushEveryOutputs_ > 1) && amrEnabled && options_.allowBatchingWithAmr)
        std::cout << "[ImplicitPIC][output] packed-field batching forced with AMR (LAMP_IMPLICIT_PACKED_BATCH_WITH_AMR=1)" << std::endl;
}

PicFieldOutputQueue::Timing PicFieldOutputQueue::writeImplicit(int stepForOutput, const Input& input)
{
    if (stepForOutput <= 0) throw std::invalid_argument("PIC output requires a positive completed step");
    Timing outTiming;
    const auto tCompute0 = Clock::now();
    if (!units_) return outTiming;
    if (!selection_.configured()) return outTiming;

    const auto& embedded = input.embedded;
    const auto& alpha = embedded.alphaByCell();
    const bool useCells = input.metrics && (options_.forceCells || !options_.disableCells);
    PicFieldFrame frame(mesh_, selection_,
        PicOutputReconstruction::buildActiveCellMask(mesh_, input.graph,
            input.metrics, input.useCutMetrics ? &alpha : nullptr), input.potential, useCells);
    if (options_.maskDiagnostics)
        PicFieldOutput::writeMaskDiagnostics(mesh_, stepForOutput, frame.validCells, input.potential, input.graph,
            input.metrics, input.useCutMetrics ? &alpha : nullptr,
            embedded.ready() && input.geometry ? &input.geometry->levelSets : nullptr,
            embedded.levelSetIndices(), embedded.insideIsNegative());

    frame.addSusceptibility(input.graph, input.muCell, input.kappaCell, input.kappaSegment);
    frame.addElectric(input.electric, *units_);
    frame.addElectrostatic(input.potential, input.density, *units_, *input.normalization);
    frame.addLevelSets(*units_, embedded.ready() ? embedded.levelSet() : nullptr,
        embedded.ready() && input.geometry ? &input.geometry->levelSets : nullptr,
        embedded.levelSetIndices(), input.levelSetNames);
    const bool hasAlpha = input.useCutMetrics && alpha.size() == static_cast<std::size_t>(mesh_.cellIndexNum());
    auto buildCache = [&](std::size_t si, CouplingWeightsCache& cache) {
        PicCurrentWeights::buildDiagnostic(input.coupling, embedded, input.species[si], cache, input.useCutMetrics);
    };
    frame.addParticles(input.species, input.weights, *input.normalization,
        hasAlpha ? &alpha : nullptr, buildCache, options_.fastDensity, &input.numberAmounts);
    if (frame.names.empty()) return outTiming;
    frame.finish(input.periodicX, input.periodicY);

    std::ostringstream tag;
    tag << std::setw(6) << std::setfill('0') << stepForOutput;
    const std::string nameNoExt = std::string("fields/fields_") + tag.str();
    outTiming.compute_s += secondsSince(tCompute0);
    if (useCells)
        outTiming.io_s += submitCells(nameNoExt, std::move(frame.names), std::move(frame.cells), frame.validCells, batchingEnabled_);
    else
        outTiming.io_s += submitVertices(nameNoExt, std::move(frame.names), std::move(frame.vertices), frame.validCells, batchingEnabled_);
    return outTiming;
}

template<class Field>
double PicFieldOutputQueue::submit(const std::string& name, std::vector<std::string>&& names,
    std::vector<Field>&& fields, const std::vector<unsigned char>& mask, bool batch)
{
    constexpr bool cellConstant = std::is_same_v<Field, CellField<double>>;
    if (batch)
    {
        Job job;
        job.nameNoExt = name;
        job.names = std::move(names);
        job.cellConstant = cellConstant;
        if constexpr (cellConstant)
            job.cellFields = std::move(fields);
        else
            job.fields = std::move(fields);
        job.cellMask = mask;
        jobs_.push_back(std::move(job));
        if (static_cast<int>(jobs_.size()) >= flushEveryOutputs_)
            return flush();
        return 0.0;
    }
    const auto tIo0 = Clock::now();
    if constexpr (cellConstant)
        PicFieldOutput::writeImplicitCells(mesh_, name, x_, y_, names, fields, &mask);
    else
        PicFieldOutput::writeImplicitVertices(mesh_, name, x_, y_, names, fields, &mask);
    return secondsSince(tIo0);
}

double PicFieldOutputQueue::submitVertices(const std::string& name, std::vector<std::string>&& names,
    std::vector<VertField<double>>&& fields, const std::vector<unsigned char>& mask, bool batch)
{
    return submit(name, std::move(names), std::move(fields), mask, batch);
}

double PicFieldOutputQueue::submitCells(const std::string& name, std::vector<std::string>&& names,
    std::vector<CellField<double>>&& fields, const std::vector<unsigned char>& mask, bool batch)
{
    return submit(name, std::move(names), std::move(fields), mask, batch);
}

double PicFieldOutputQueue::flush()
{
    if (jobs_.empty()) return 0.0;
    const auto tIo0 = Clock::now();
    for (auto& job : jobs_)
    {
        if (job.cellConstant)
            PicFieldOutput::writeImplicitCells(mesh_, job.nameNoExt, x_, y_, job.names, job.cellFields, &job.cellMask);
        else
            PicFieldOutput::writeImplicitVertices(mesh_, job.nameNoExt, x_, y_, job.names, job.fields, &job.cellMask);
    }
    // Preserve the old failure behavior: an exception leaves all jobs pending,
    // including any files already written. The caller stops the run.
    jobs_.clear();
    return secondsSince(tIo0);
}
