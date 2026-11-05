#pragma once
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace LAMP
{
// Supplied by the App's clock, never advanced by a state container.
struct TimeLayer
{
    enum class Position { Initial, Midpoint, End };
    std::int64_t step = 0;
    double begin = 0.0, end = 0.0;
    Position position = Position::Initial;
    double time() const // validates the interval; Initial permits restart time
    {
        // Preserve the existing validation exception text for all consumers.
        const auto require = [](bool condition, const std::string& message) {
            if (!condition) throw std::invalid_argument(
                std::string("[Interface/HermiteGalsState id=0] ") + message);
        };
        require(step >= 0 && std::isfinite(begin) && std::isfinite(end), "time: invalid step or endpoints");
        if (position == Position::Initial) {
            require(begin == end, "time: initial/restart layer must have zero interval");
            return begin;
        }
        require(position == Position::Midpoint || position == Position::End, "time: unknown position");
        require(step > 0 && end > begin && std::isfinite(end - begin), "time: invalid advancing interval");
        return position == Position::Midpoint ? begin + (end - begin) * 0.5 : end;
    }
};
}
