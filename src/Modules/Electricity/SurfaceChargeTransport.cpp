#include "Modules/Electricity/SurfaceChargeTransport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Quantum::Electricity
{
namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(message);
}

double finite(double value, const char* message)
{
    require(std::isfinite(value), message);
    return value;
}

double product(double left, double right, const char* message)
{
    const double value = finite(left * right, message);
    require(left == 0.0 || right == 0.0 || value != 0.0,
            "nonzero charge/flux product underflow");
    return value;
}

// Compensated totals are important for interfaces carrying both charge signs.
struct Sum
{
    double value = 0.0, correction = 0.0;
    void add(double term)
    {
        const double adjusted = term - correction;
        const double next = finite(value + adjusted, "charge balance sum overflow");
        correction = (next - value) - adjusted;
        value = next;
    }
};

void validateMeasure(const std::vector<double>& measure)
{
    require(measure.size() >= 3, "interface requires at least three panels");
    for (double length : measure)
        require(std::isfinite(length) && length > 0.0,
                "panel measure must be finite and positive");
}

void validateIntegratedCharge(const std::vector<double>& measure, const std::vector<double>& charge)
{
    validateMeasure(measure);
    require(charge.size() == measure.size(), "integrated charge/panel count mismatch");
    for (std::size_t i = 0; i < charge.size(); ++i)
    {
        finite(charge[i], "integrated charge is not finite");
        const double density = finite(charge[i] / measure[i], "integrated charge density overflow");
        require(charge[i] == 0.0 || density != 0.0, "integrated charge density underflow");
    }
}
}

SurfaceChargeTransport::SurfaceChargeTransport(
    const std::vector<double>& panelMeasure, const std::vector<double>& initialDensity, bool closed)
    : measure_(panelMeasure), charge_(initialDensity.size()), closed_(closed)
{
    validateMeasure(measure_);
    require(initialDensity.size() == measure_.size(), "density/panel count mismatch");
    for (std::size_t i = 0; i < charge_.size(); ++i)
    {
        finite(initialDensity[i], "initial charge density is not finite");
        charge_[i] = product(initialDensity[i], measure_[i], "initial charge overflow");
        finite(charge_[i] / measure_[i], "initial stored density overflow");
    }
}

SurfaceChargeTransport::SurfaceChargeTransport(const std::vector<double>& panelMeasure,
    const std::vector<double>& integratedCharge, IntegratedChargeTag, bool closed)
    : measure_(panelMeasure), charge_(integratedCharge), closed_(closed)
{
    validateIntegratedCharge(measure_, charge_);
}

SurfaceChargeTransport SurfaceChargeTransport::fromIntegratedCharge(
    const std::vector<double>& panelMeasure, const std::vector<double>& integratedCharge, bool closed)
{
    return SurfaceChargeTransport(panelMeasure, integratedCharge, IntegratedChargeTag{}, closed);
}

SurfaceChargeStepReport SurfaceChargeTransport::remapWeights(
    const std::vector<std::vector<PanelWeight>>& rows,const std::vector<double>& newMeasure)
{
    SurfaceChargeStepReport report;report.stage="prepared remap";
    try {
        validateMeasure(newMeasure);
        require(rows.size()==charge_.size(),"prepared remap source count mismatch");
        std::vector<Sum> deposited(newMeasure.size());std::vector<bool> covered(newMeasure.size(),false);
        Sum before,scale,after;
        constexpr double roundoff=512*std::numeric_limits<double>::epsilon();
        for (std::size_t i=0;i<rows.size();++i) {
            require(!rows[i].empty(),"prepared remap lost a donor");
            long double total=0;
            for (const auto& weight:rows[i]) {
                require(weight.targetPanel<newMeasure.size() && std::isfinite(weight.share)
                    && weight.share>0 && weight.share<=1,"invalid prepared remap weight");
                total+=weight.share;
                const auto value=static_cast<double>(static_cast<long double>(charge_[i])*weight.share);
                finite(value,"prepared charge deposit overflow");
                require(charge_[i]==0 || value!=0,"prepared charge deposit underflow");
                deposited[weight.targetPanel].add(value);covered[weight.targetPanel]=true;
            }
            require(std::abs(total-1)<=128*std::numeric_limits<long double>::epsilon()*rows[i].size(),
                "prepared remap weights do not sum to one");
            before.add(charge_[i]);scale.add(std::abs(charge_[i]));
        }
        std::vector<double> nextCharge(newMeasure.size()),nextMeasure(newMeasure);
        for (std::size_t j=0;j<nextCharge.size();++j) {
            require(covered[j],"prepared remap omitted a target panel");
            nextCharge[j]=deposited[j].value;after.add(nextCharge[j]);
        }
        validateIntegratedCharge(nextMeasure,nextCharge);
        report.totalChargeBefore=before.value;report.totalChargeAfter=after.value;
        report.conservationScale=scale.value;
        report.conservationResidual=finite(after.value-before.value,"prepared remap balance overflow");
        report.relativeConservationResidual=scale.value==0?0:std::abs(report.conservationResidual)/scale.value;
        require(report.relativeConservationResidual<=roundoff,"prepared remap charge conservation exceeds roundoff");
        for (std::size_t j=0;j<nextCharge.size();++j)
            report.maxAbsDensity=std::max(report.maxAbsDensity,std::abs(nextCharge[j]/nextMeasure[j]));
        report.stage="complete";report.success=true;
        measure_.swap(nextMeasure);charge_.swap(nextCharge);
    } catch (const std::exception& error) {report.reason=error.what();}
    return report;
}

SurfaceChargeStepReport SurfaceChargeTransport::remap(
    const std::vector<double>& mappedOldBoundaries, const std::vector<double>& newBoundaries)
{
    SurfaceChargeStepReport report;
    try
    {
        report.stage = "remap input";
        require(mappedOldBoundaries.size() == charge_.size() + 1 && newBoundaries.size() >= 4,
            "remap boundary/panel count mismatch");
        const auto checkBoundaries = [](const std::vector<double>& boundary)
        {
            for (std::size_t i = 0; i < boundary.size(); ++i)
            {
                finite(boundary[i], "nonfinite remap boundary");
                if (i) require(boundary[i] > boundary[i - 1], "remap boundaries fold or repeat");
            }
        };
        checkBoundaries(mappedOldBoundaries);
        checkBoundaries(newBoundaries);
        require(newBoundaries.front() == 0.0, "new remap boundaries must start at zero");
        const double period = newBoundaries.back();
        require(period > 0.0, "remap period must be positive");
        const double twoPeriods = closed_ ? finite(period + period, "remap period range overflow") : period;
        constexpr double roundoff = 512.0 * std::numeric_limits<double>::epsilon();
        const double span = finite(mappedOldBoundaries.back() - mappedOldBoundaries.front(),
            "mapped old period span overflow");
        require(std::abs(span - period) <= roundoff * period,
            "mapped old boundaries do not span exactly one period within roundoff");

        require(closed_ || (mappedOldBoundaries.front() == 0.0 && mappedOldBoundaries.back() == period),
            "open remap boundaries must share exact zero and final endpoints");
        double origin = closed_ ? finite(std::fmod(mappedOldBoundaries.front(), period), "remap origin reduction failed") : 0.0;
        if (closed_ && origin < 0.0) origin = finite(origin + period, "remap origin reduction overflow");
        require(origin >= 0.0 && origin < period, "unresolved periodic remap origin");
        std::vector<double> oldBoundary(mappedOldBoundaries.size());
        for (std::size_t i = 0; i < oldBoundary.size(); ++i)
        {
            oldBoundary[i] = closed_ ? finite(origin + (mappedOldBoundaries[i] - mappedOldBoundaries.front()),
                "remap boundary reduction overflow") : mappedOldBoundaries[i];
            require(oldBoundary[i] >= 0.0 && oldBoundary[i] <= twoPeriods,
                "reduced remap boundaries leave two-period range");
            if (i) require(oldBoundary[i] > oldBoundary[i - 1], "remap origin reduction collapsed a panel");
        }
        const std::size_t newCount = newBoundaries.size() - 1;
        std::vector<double> nextMeasure(newCount), nextCharge(newCount);
        for (std::size_t i = 0; i < newCount; ++i)
            nextMeasure[i] = finite(newBoundaries[i + 1] - newBoundaries[i], "new remap measure overflow");
        validateMeasure(nextMeasure);
        std::vector<Sum> deposited(newCount), covered(newCount);

        report.stage = "remap overlap";
        std::size_t target = static_cast<std::size_t>(
            std::upper_bound(newBoundaries.begin(), newBoundaries.end(), origin) - newBoundaries.begin() - 1);
        Sum before, scale;
        for (std::size_t old = 0; old < charge_.size(); ++old)
        {
            const double length = finite(mappedOldBoundaries[old + 1] - mappedOldBoundaries[old],
                "mapped old panel measure overflow");
            require(length > 0.0, "mapped old panel measure is not positive");
            Sum donorWeight;
            while (target / newCount < (closed_ ? 2u : 1u))
            {
                const std::size_t cell = closed_ ? target % newCount : target;
                const double shift = !closed_ || target < newCount ? 0.0 : period;
                const double left = newBoundaries[cell] + shift;
                const double right = newBoundaries[cell + 1] + shift;
                require(right > left, "period translation collapsed a new remap panel");
                if (left >= oldBoundary[old + 1]) break;
                const double overlap = std::min(right, oldBoundary[old + 1]) - std::max(left, oldBoundary[old]);
                if (overlap > 0.0)
                {
                    const double fraction = finite(overlap / length, "remap overlap fraction overflow");
                    require(fraction > 0.0, "positive remap overlap fraction underflow");
                    deposited[cell].add(product(charge_[old], fraction, "remapped charge contribution overflow"));
                    covered[cell].add(overlap);
                    donorWeight.add(fraction);
                }
                if (right > oldBoundary[old + 1]) break;
                ++target;
            }
            require(std::abs(donorWeight.value - 1.0) <= roundoff,
                "mapped old panel is not completely covered within roundoff");
            before.add(charge_[old]);
            scale.add(std::abs(charge_[old]));
        }
        Sum after;
        for (std::size_t i = 0; i < newCount; ++i)
        {
            require(std::abs(covered[i].value / nextMeasure[i] - 1.0) <= roundoff,
                "new remap panel has a gap or duplicate coverage");
            nextCharge[i] = deposited[i].value;
            after.add(nextCharge[i]);
        }
        validateIntegratedCharge(nextMeasure, nextCharge);
        report.stage = "remap balance";
        report.totalChargeBefore = before.value;
        report.totalChargeAfter = after.value;
        report.conservationScale = scale.value;
        report.conservationResidual = finite(after.value - before.value, "remap charge balance overflow");
        report.relativeConservationResidual = scale.value == 0.0 ? 0.0 :
            std::abs(report.conservationResidual) / scale.value;
        require(report.relativeConservationResidual <= roundoff, "remap charge conservation exceeds roundoff");
        for (std::size_t i = 0; i < newCount; ++i)
            report.maxAbsDensity = std::max(report.maxAbsDensity, std::abs(nextCharge[i] / nextMeasure[i]));
        report.stage = "complete";
        report.success = true;
        measure_.swap(nextMeasure);
        charge_.swap(nextCharge);
    }
    catch (const std::exception& error)
    {
        report.success = false;
        report.reason = error.what();
    }
    return report;
}

SurfaceChargeStepReport SurfaceChargeTransport::remap(
    const std::vector<PanelCoordinate>& mappedOldBoundaries, const std::vector<double>& newPanelMeasure)
{
    SurfaceChargeStepReport report;
    try
    {
        report.stage = "local remap input";
        require(!closed_, "panel-local remap requires an open chain");
        require(mappedOldBoundaries.size() == charge_.size() + 1,
            "local remap boundary/old panel count mismatch");
        validateMeasure(newPanelMeasure);
        const std::size_t newCount = newPanelMeasure.size();
        std::vector<double> nextMeasure(newPanelMeasure), nextCharge(newCount);
        std::vector<PanelCoordinate> boundary(mappedOldBoundaries);
        for (std::size_t i = 0; i < boundary.size(); ++i)
        {
            auto& point = boundary[i];
            require(point.panelIndex <= newCount && std::isfinite(point.localMeasure)
                && point.localMeasure >= 0.0, "invalid panel-local remap boundary");
            if (point.panelIndex == newCount)
                require(point.localMeasure == 0.0, "final panel-local sentinel must have zero local measure");
            else
            {
                const double measure = nextMeasure[point.panelIndex];
                require(point.localMeasure <= measure, "remap boundary exceeds its actual new panel measure");
                if (point.localMeasure == measure)
                {
                    // A change of shared-boundary identity, not a projection,
                    // tolerance-based merge or numeric adjustment of a root.
                    ++point.panelIndex;
                    point.localMeasure = 0.0;
                }
            }
            if (i != 0)
            {
                const auto& previous = boundary[i-1];
                require(point.panelIndex > previous.panelIndex
                    || (point.panelIndex == previous.panelIndex && point.localMeasure > previous.localMeasure),
                    "panel-local remap boundaries fold or repeat a shared identity");
            }
        }
        require(boundary.front().panelIndex == 0 && boundary.front().localMeasure == 0.0
            && boundary.back().panelIndex == newCount && boundary.back().localMeasure == 0.0,
            "open panel-local remap must cover the exact first and final endpoints");

        constexpr double roundoff = 512.0 * std::numeric_limits<double>::epsilon();
        std::vector<Sum> deposited(newCount), covered(newCount);
        struct Overlap { std::size_t panel; double measure; };
        std::vector<Overlap> overlaps;
        Sum before, scale;
        report.stage = "local remap overlap";
        for (std::size_t old = 0; old < charge_.size(); ++old)
        {
            const auto& begin = boundary[old];
            const auto& end = boundary[old+1];
            overlaps.clear();
            Sum donorMeasure;
            for (std::size_t target = begin.panelIndex;
                target < end.panelIndex || (target == end.panelIndex && end.localMeasure > 0.0); ++target)
            {
                require(target < newCount, "mapped donor leaves the open new chain");
                const double left = target == begin.panelIndex ? begin.localMeasure : 0.0;
                const double right = target == end.panelIndex ? end.localMeasure : nextMeasure[target];
                const double overlap = finite(right-left, "panel-local overlap measure overflow");
                require(overlap > 0.0, "positive panel-local overlap is not representable");
                overlaps.push_back({target,overlap});
                donorMeasure.add(overlap);
            }
            // The finite denominator is an approximation, just as in the old
            // floating overlap rule. A tiny panel need not change its rounded
            // value to be retained: its own positive overlap below still gets
            // an independent fraction and charge contribution. Never recover
            // that overlap by subtracting two ordinary-sized cumulative sums.
            require(std::isfinite(donorMeasure.value) && donorMeasure.value > 0.0,
                "mapped old panel has no representable positive measure");
            Sum donorWeight;
            for (const auto& overlap : overlaps)
            {
                const double fraction = finite(overlap.measure/donorMeasure.value,
                    "panel-local overlap fraction overflow");
                require(fraction > 0.0, "positive panel-local overlap fraction underflow");
                deposited[overlap.panel].add(product(charge_[old],fraction,
                    "panel-local remapped charge contribution overflow"));
                covered[overlap.panel].add(overlap.measure);
                donorWeight.add(fraction);
            }
            require(std::abs(donorWeight.value-1.0) <= roundoff,
                "mapped old panel is not completely covered within roundoff");
            before.add(charge_[old]);
            scale.add(std::abs(charge_[old]));
        }
        Sum after;
        for (std::size_t i = 0; i < newCount; ++i)
        {
            require(std::abs(covered[i].value/nextMeasure[i]-1.0) <= roundoff,
                "new remap panel has a gap or duplicate coverage");
            nextCharge[i] = deposited[i].value;
            after.add(nextCharge[i]);
        }
        validateIntegratedCharge(nextMeasure,nextCharge);
        report.stage = "local remap balance";
        report.totalChargeBefore = before.value;
        report.totalChargeAfter = after.value;
        report.conservationScale = scale.value;
        report.conservationResidual = finite(after.value-before.value,"remap charge balance overflow");
        report.relativeConservationResidual = scale.value == 0.0 ? 0.0
            : std::abs(report.conservationResidual)/scale.value;
        require(report.relativeConservationResidual <= roundoff,"remap charge conservation exceeds roundoff");
        for (std::size_t i = 0; i < newCount; ++i)
            report.maxAbsDensity = std::max(report.maxAbsDensity,std::abs(nextCharge[i]/nextMeasure[i]));
        // All checks and diagnostic allocations precede the nonthrowing commit.
        report.stage = "complete";
        report.success = true;
        measure_.swap(nextMeasure);
        charge_.swap(nextCharge);
    }
    catch (const std::exception& error)
    {
        report.success = false;
        report.reason = error.what();
    }
    return report;
}

std::vector<double> SurfaceChargeTransport::density() const
{
    std::vector<double> result(charge_.size());
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = charge_[i] / measure_[i];
    return result;
}

SurfaceChargeStepReport SurfaceChargeTransport::advance(
    const std::vector<double>& nextPanelMeasure,
    const std::vector<double>& relativeAdvectiveMeasureFlux,
    const std::vector<double>& integratedNormalCurrent, double dt)
{
    return advanceImpl(nextPanelMeasure, relativeAdvectiveMeasureFlux, integratedNormalCurrent, dt, nullptr);
}

SurfaceChargeStepReport SurfaceChargeTransport::advance(
    const std::vector<double>& nextPanelMeasure,
    const std::vector<double>& relativeAdvectiveMeasureFlux,
    const std::vector<double>& integratedNormalCurrent, double dt,
    const SurfaceChargeBoundaryState& boundary)
{
    return advanceImpl(nextPanelMeasure, relativeAdvectiveMeasureFlux, integratedNormalCurrent, dt, &boundary);
}

SurfaceChargeStepReport SurfaceChargeTransport::advanceImpl(
    const std::vector<double>& nextPanelMeasure,
    const std::vector<double>& relativeAdvectiveMeasureFlux,
    const std::vector<double>& integratedNormalCurrent, double dt,
    const SurfaceChargeBoundaryState* boundary)
{
    SurfaceChargeStepReport report;
    try
    {
        report.stage = "input";
        require(std::isfinite(dt) && dt > 0.0, "dt must be finite and positive");
        const std::size_t count = charge_.size();
        const std::size_t fluxCount = count + (closed_ ? 0 : 1);
        require(nextPanelMeasure.size() == count &&
                relativeAdvectiveMeasureFlux.size() == fluxCount &&
                integratedNormalCurrent.size() == count,
                "panel count/order must remain unchanged");
        validateMeasure(nextPanelMeasure);
        for (std::size_t i = 0; i < fluxCount; ++i)
            finite(relativeAdvectiveMeasureFlux[i], "advective measure flux is not finite");
        if (boundary)
        {
            require(!closed_, "external endpoint states require an open interface");
            if (boundary->startDensity) finite(*boundary->startDensity, "external start density is not finite");
            if (boundary->endDensity) finite(*boundary->endDensity, "external end density is not finite");
            require(relativeAdvectiveMeasureFlux.front() <= 0.0 || boundary->startDensity.has_value(),
                "start endpoint inflow requires an explicit external charge density");
            require(relativeAdvectiveMeasureFlux.back() >= 0.0 || boundary->endDensity.has_value(),
                "end endpoint inflow requires an explicit external charge density");
        }
        else
            require(closed_ || (relativeAdvectiveMeasureFlux.front() == 0.0 && relativeAdvectiveMeasureFlux.back() == 0.0),
                "open interface requires exactly zero endpoint measure fluxes");
        for (std::size_t i = 0; i < count; ++i)
        {
            finite(integratedNormalCurrent[i], "normal current is not finite");
        }

        report.stage = "advection";
        report.advectiveChargeFlux.resize(fluxCount);
        constexpr double roundoff = 32.0 * std::numeric_limits<double>::epsilon();
        for (std::size_t i = 0; i < count; ++i)
        {
            const std::size_t previous = closed_ ? (i + count - 1) % count : i;
            const std::size_t outlet = closed_ ? i : i + 1;
            const double outflow = finite(std::max(0.0, relativeAdvectiveMeasureFlux[outlet]) +
                std::max(0.0, -relativeAdvectiveMeasureFlux[previous]), "outflow overflow");
            const double cfl = finite((outflow / measure_[i]) * dt, "advective CFL overflow");
            report.maxAdvectiveCfl = std::max(report.maxAdvectiveCfl, cfl);
            require(cfl <= 1.0 + roundoff, "surface-charge advective CFL exceeds one");
            if (!closed_ && outlet == count) continue; // Endpoint donor is handled separately below.
            const std::size_t donor = relativeAdvectiveMeasureFlux[outlet] >= 0.0 ? i : (closed_ ? (i + 1) % count : i + 1);
            report.advectiveChargeFlux[outlet] = product(relativeAdvectiveMeasureFlux[outlet],
                charge_[donor] / measure_[donor], "advective charge flux overflow");
        }
        if (boundary)
        {
            const double first = relativeAdvectiveMeasureFlux.front();
            const double last = relativeAdvectiveMeasureFlux.back();
            const double firstDensity = first > 0.0 ? *boundary->startDensity : charge_.front() / measure_.front();
            const double lastDensity = last < 0.0 ? *boundary->endDensity : charge_.back() / measure_.back();
            report.advectiveChargeFlux.front() = product(first, firstDensity, "start endpoint charge flux overflow");
            report.advectiveChargeFlux.back() = product(last, lastDensity, "end endpoint charge flux overflow");
            Sum inflow, outflow;
            if (first > 0.0) inflow.add(report.advectiveChargeFlux.front());
            else if (first < 0.0) outflow.add(-report.advectiveChargeFlux.front());
            if (last < 0.0) inflow.add(-report.advectiveChargeFlux.back());
            else if (last > 0.0) outflow.add(report.advectiveChargeFlux.back());
            report.boundaryChargeInflow = inflow.value;
            report.boundaryChargeOutflow = outflow.value;
            report.netBoundaryCurrent = finite(inflow.value - outflow.value, "net endpoint charge current overflow");
        }

        report.stage = "update";
        std::vector<double> nextCharge(count), nextMeasure(nextPanelMeasure);
        Sum before, after, current, balanceScale;
        for (std::size_t i = 0; i < count; ++i)
        {
            const std::size_t previous = closed_ ? (i + count - 1) % count : i;
            const std::size_t outlet = closed_ ? i : i + 1;
            const double source = product(dt, integratedNormalCurrent[i], "current increment overflow");
            const double inflow = product(dt, report.advectiveChargeFlux[previous], "advective increment overflow");
            const double outflow = product(dt, report.advectiveChargeFlux[outlet], "advective increment overflow");
            nextCharge[i] = finite(charge_[i] + source + (inflow - outflow), "updated charge overflow");
            const double nextDensity = finite(nextCharge[i] / nextMeasure[i], "updated density overflow");
            require(nextCharge[i] == 0.0 || nextDensity != 0.0, "updated density underflow");
            report.maxAbsDensity = std::max(report.maxAbsDensity, std::abs(nextDensity));
            before.add(charge_[i]);
            after.add(nextCharge[i]);
            current.add(integratedNormalCurrent[i]);
            balanceScale.add(std::abs(charge_[i]));
            balanceScale.add(std::abs(source));
        }
        report.stage = "balance";
        report.totalChargeBefore = before.value;
        report.totalChargeAfter = after.value;
        report.integratedNormalCurrent = current.value;
        if (boundary)
        {
            balanceScale.add(product(dt, std::abs(report.advectiveChargeFlux.front()), "endpoint balance scale overflow"));
            balanceScale.add(product(dt, std::abs(report.advectiveChargeFlux.back()), "endpoint balance scale overflow"));
        }
        report.conservationScale = balanceScale.value;
        report.conservationResidual = finite(after.value - before.value - dt * current.value,
                                             "charge balance residual overflow");
        if (boundary)
            report.conservationResidual = finite(report.conservationResidual
                - product(dt, report.netBoundaryCurrent, "endpoint current increment overflow"),
                "open-boundary charge balance residual overflow");
        report.relativeConservationResidual = balanceScale.value == 0.0 ? 0.0 :
            std::abs(report.conservationResidual) / balanceScale.value;
        require(report.relativeConservationResidual <= 128.0 * std::numeric_limits<double>::epsilon(),
                "surface-charge conservation residual exceeds roundoff");

        // All validation, diagnostic allocations and strings precede commit.
        report.stage = "complete";
        report.success = true;
        measure_.swap(nextMeasure);
        charge_.swap(nextCharge);
    }
    catch (const std::exception& error)
    {
        report.success = false;
        report.reason = error.what();
    }
    return report;
}
}
