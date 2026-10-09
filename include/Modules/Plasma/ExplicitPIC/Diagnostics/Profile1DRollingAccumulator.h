#pragma once

#include <cstddef>
#include <vector>

namespace ExplicitPIC
{
namespace Diagnostics
{

// Rolling accumulator for 1D profile quantities.
//
// Layout for a single sample is field-major:
//   sample[field * nBins + bin]
//
// It maintains per-entry sum/sum2/count and supports a fixed-size sliding
// window through an internal ring buffer.
class Profile1DRollingAccumulator
{
public:
    Profile1DRollingAccumulator() = default;

    // Configure sizes and reset state.
    // windowSteps==0 disables the accumulator.
    void reset(std::size_t nBins, std::size_t nFields, std::size_t windowSteps);

    // Clear current state but keep sizes/window.
    void clear();

    bool enabled() const;

    std::size_t nBins() const { return nBins_; }
    std::size_t nFields() const { return nFields_; }
    std::size_t windowSteps() const { return windowSteps_; }
    std::size_t count() const { return count_; }
    std::size_t sampleStride() const { return nBins_ * nFields_; }

    // Push a new sample of length nBins*nFields in field-major layout.
    void pushSample(const std::vector<double>& sample);
    void pushSample(const double* sample, std::size_t len);

    // Compute mean/std for each (field,bin) entry; outputs are size nBins*nFields.
    void mean(std::vector<double>& outMean) const;
    void stddev(std::vector<double>& outStd) const;
    void meanStd(std::vector<double>& outMean, std::vector<double>& outStd) const;

    // On-demand accessors.
    double meanAt(std::size_t field, std::size_t bin) const;
    double stddevAt(std::size_t field, std::size_t bin) const;

private:
    std::size_t nBins_ = 0;
    std::size_t nFields_ = 0;
    std::size_t windowSteps_ = 0;

    std::size_t count_ = 0; // number of samples currently kept (<=windowSteps_)
    std::size_t head_ = 0;  // next write slot in ring buffer

    std::vector<double> sum_;
    std::vector<double> sum2_;
    std::vector<double> ring_; // windowSteps_ * sampleStride()

    std::size_t off(std::size_t field, std::size_t bin) const { return field * nBins_ + bin; }

    const double* slotPtr(std::size_t slot) const;
    double* slotPtr(std::size_t slot);

    void addSample(const double* s);
    void subSample(const double* s);
};

} // namespace Diagnostics
} // namespace ExplicitPIC
