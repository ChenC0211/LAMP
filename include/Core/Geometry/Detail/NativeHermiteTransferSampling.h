#pragma once
// Implementation-only shared exact-chart sampling. Include after the native
// scalar owner; its TU-local types never cross the public ABI.
namespace LAMP::Geometry::Detail
{
namespace
{
namespace NativeHermiteTransfer
{
namespace Scalar = LevelSetHermiteScalar;
namespace Exact = ExactExpansion;
inline void require(bool valid, const char* reason)
{ if (!valid) throw std::runtime_error(std::string("[NativeHermiteTransfer] ") + reason); }
// Exact physical datum coordinates stay expansions until the old cell and
// its rational local chart have been selected. In particular a target edge
// midpoint is not rounded to an adjacent old vertex before sampling.
template<class Source> class DatumSampler
{
  public:
    explicit DatumSampler(const Source& source) : source_(source) {}
    struct Location
    {
        int cell = -1;
        Scalar::RationalSegmentCoordinate x, y;
    };
    Location locate(Vector2D a, Vector2D b) const
    {
        const auto midpoint = [](double x, double y) {
            return Exact::scalePow2(Exact::add(Exact::fromDouble(x), Exact::fromDouble(y)), -1);
        };
        const auto x = midpoint(a.x, b.x), y = midpoint(a.y, b.y);
        const double px = Exact::estimate(x), py = Exact::estimate(y);
        const double infinity = std::numeric_limits<double>::infinity();
        auto candidates = source_.leafCandidates(
            {std::nextafter(px, -infinity), std::nextafter(py, -infinity)},
            {std::nextafter(px, infinity), std::nextafter(py, infinity)});
        std::sort(candidates.begin(), candidates.end());
        for (int c : candidates)
        {
            if (!source_.cellActive(c)) continue;
            const auto box = source_.leafBounds(c);
            const auto chart = [](const Exact::Expansion& point, double low, double high) {
                Scalar::RationalSegmentCoordinate result;
                result.numeratorBegin = Exact::subtract(point, Exact::fromDouble(low));
                result.numeratorEnd = result.numeratorBegin;
                result.denominator = Exact::subtract(Exact::fromDouble(high), Exact::fromDouble(low));
                return result;
            };
            auto cx = chart(x, box[0].x, box[1].x), cy = chart(y, box[0].y, box[1].y);
            const auto inside = [](const auto& coordinate) {
                return Exact::sign(coordinate.numeratorBegin) >= 0
                    && Exact::sign(Exact::subtract(coordinate.denominator, coordinate.numeratorBegin)) >= 0;
            };
            if (inside(cx) && inside(cy)) return {c, std::move(cx), std::move(cy)};
        }
        throw std::runtime_error("[NativeHermiteTransfer] target datum has no active source chart");
    }
    double value(const Location& point, unsigned dx = 0, unsigned dy = 0) const
    {
        const auto polynomial = exactValue(point, dx, dy);
        const double result = Exact::estimate(polynomial.numerator[0]) / Exact::estimate(polynomial.denominator);
        require(std::isfinite(result), "nonfinite restricted source datum");
        require(result != 0 || Exact::sign(polynomial.numerator[0]) == 0, "restricted source datum underflow");
        return result;
    }
    // Retained homogeneous constant polynomial at the exact physical datum.
    // Consumers that need an exact restriction must keep this ratio intact.
    auto exactValue(const Location& point, unsigned dx = 0, unsigned dy = 0) const
    { return source_.restrictLocal(point.cell, point.x, point.y, dx, dy); }
  private:
    const Source& source_;
};

}
}
}
