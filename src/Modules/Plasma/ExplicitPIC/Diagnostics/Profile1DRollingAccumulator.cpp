#include "Modules/Plasma/ExplicitPIC/Diagnostics/Profile1DRollingAccumulator.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace ExplicitPIC
{
namespace Diagnostics
{

bool Profile1DRollingAccumulator::enabled() const
{
    return windowSteps_ > 0 && nBins_ > 0 && nFields_ > 0;
}

void Profile1DRollingAccumulator::reset(std::size_t nBins, std::size_t nFields, std::size_t windowSteps)
{
    nBins_ = nBins;
    nFields_ = nFields;
    windowSteps_ = windowSteps;

    count_ = 0;
    head_ = 0;

    const std::size_t stride = sampleStride();
    sum_.assign(stride, 0.0);
    sum2_.assign(stride, 0.0);

    if (windowSteps_ > 0)
    {
        ring_.assign(windowSteps_ * stride, 0.0);
    }
    else
    {
        ring_.clear();
    }
}

void Profile1DRollingAccumulator::clear()
{
    count_ = 0;
    head_ = 0;
    std::fill(sum_.begin(), sum_.end(), 0.0);
    std::fill(sum2_.begin(), sum2_.end(), 0.0);
    std::fill(ring_.begin(), ring_.end(), 0.0);
}

const double* Profile1DRollingAccumulator::slotPtr(std::size_t slot) const
{
    const std::size_t stride = sampleStride();
    return ring_.data() + slot * stride;
}

double* Profile1DRollingAccumulator::slotPtr(std::size_t slot)
{
    const std::size_t stride = sampleStride();
    return ring_.data() + slot * stride;
}

void Profile1DRollingAccumulator::addSample(const double* s)
{
    const std::size_t stride = sampleStride();
    for (std::size_t i = 0; i < stride; ++i)
    {
        const double x = s[i];
        sum_[i] += x;
        sum2_[i] += x * x;
    }
}

void Profile1DRollingAccumulator::subSample(const double* s)
{
    const std::size_t stride = sampleStride();
    for (std::size_t i = 0; i < stride; ++i)
    {
        const double x = s[i];
        sum_[i] -= x;
        sum2_[i] -= x * x;
    }
}

void Profile1DRollingAccumulator::pushSample(const std::vector<double>& sample)
{
    pushSample(sample.data(), sample.size());
}

void Profile1DRollingAccumulator::pushSample(const double* sample, std::size_t len)
{
    if (!enabled())
    {
        return;
    }

    const std::size_t stride = sampleStride();
    if (len != stride)
    {
        throw std::runtime_error("Profile1DRollingAccumulator::pushSample: sample length mismatch");
    }

    // head_ always points to the next insertion slot.
    // When the buffer is full, that slot contains the oldest sample.
    if (count_ == windowSteps_)
    {
        const double* old = slotPtr(head_);
        subSample(old);
    }
    else
    {
        ++count_;
    }

    double* dst = slotPtr(head_);
    std::copy(sample, sample + stride, dst);
    addSample(dst);

    head_ = (head_ + 1) % windowSteps_;
}

void Profile1DRollingAccumulator::mean(std::vector<double>& outMean) const
{
    const std::size_t stride = sampleStride();
    outMean.assign(stride, 0.0);

    if (count_ == 0)
    {
        return;
    }

    const double inv = 1.0 / static_cast<double>(count_);
    for (std::size_t i = 0; i < stride; ++i)
    {
        outMean[i] = sum_[i] * inv;
    }
}

void Profile1DRollingAccumulator::stddev(std::vector<double>& outStd) const
{
    std::vector<double> m;
    meanStd(m, outStd);
}

void Profile1DRollingAccumulator::meanStd(std::vector<double>& outMean, std::vector<double>& outStd) const
{
    const std::size_t stride = sampleStride();
    outMean.assign(stride, 0.0);
    outStd.assign(stride, 0.0);

    if (count_ == 0)
    {
        return;
    }

    const double inv = 1.0 / static_cast<double>(count_);
    for (std::size_t i = 0; i < stride; ++i)
    {
        const double mu = sum_[i] * inv;
        // population variance: E[x^2] - (E[x])^2
        double var = sum2_[i] * inv - mu * mu;
        if (var < 0.0)
        {
            // numeric roundoff
            var = 0.0;
        }
        outMean[i] = mu;
        outStd[i] = std::sqrt(var);
    }
}

double Profile1DRollingAccumulator::meanAt(std::size_t field, std::size_t bin) const
{
    if (count_ == 0)
    {
        return 0.0;
    }
    const std::size_t i = off(field, bin);
    assert(i < sum_.size());
    return sum_[i] / static_cast<double>(count_);
}

double Profile1DRollingAccumulator::stddevAt(std::size_t field, std::size_t bin) const
{
    if (count_ == 0)
    {
        return 0.0;
    }
    const std::size_t i = off(field, bin);
    assert(i < sum_.size());
    const double inv = 1.0 / static_cast<double>(count_);
    const double mu = sum_[i] * inv;
    double var = sum2_[i] * inv - mu * mu;
    if (var < 0.0)
    {
        var = 0.0;
    }
    return std::sqrt(var);
}

} // namespace Diagnostics
} // namespace ExplicitPIC
