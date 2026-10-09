#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "Modules/Electricity/Detail/DielectricInterfaceTrace.h"

// Shared implementation-only ordered face/link geometry. No shape, fit, material law,
// quadrature, or global solver lives here. Actual phase labels are supplied
// by the caller; they are not inferred from a local interface tangent plane.
// Implementation-only: the anonymous ancestor gives every definition and
// type TU-local linkage, so different floating-point modes cannot share these
// engines through COMDAT. All includes stay outside these namespaces; no
// Detail type is a cross-TU public API.
namespace Quantum::Electricity::Detail
{
namespace
{
namespace DielectricTraceFaceGeometry
{
using Phase = DielectricInterfaceTrace::Phase;

enum class LinkGeometry
{
    Orthogonal,
    ProjectedNormalCenters
};

struct SegmentPiece
{
    Vector2D begin;
    Vector2D end;
    Phase phase = Phase::Negative;
    // Native mapped intervals retain width independently of their rounded
    // location. begin/end are then the WHOLE input face, not phase endpoints.
    bool retainedRange = false, fromEnd = false;
    long double offset = 0.0L, fraction = 0.0L;
};

struct LinkCoordinates
{
    std::array<double, 2> owner;
    std::array<double, 2> neighbor;
};

// lengthScale and normalized come from an already validated reconstruction.
// Visit each checked piece before proceeding to the next piece or the final
// whole-face checks. This deliberately preserves the P2 validation/integration
// exception order and its midpoint arithmetic, rather than doing two passes.
// The visitor receives normalized begin/end/midpoint and normalized length.
// Positive pieces must be roundoff-resolved; none are silently discarded.
// errorPrefix preserves each candidate's existing public exception prefix.
template <class Normalize, class VisitPiece>
inline LinkCoordinates visitValidatedFace(const Vector2D &owner, Phase ownerPhase, const Vector2D &neighbor,
                                          Phase neighborPhase, const Vector2D &ownerFaceNormal,
                                          double faceLength, double linkDistance, double twoPointPermittivity,
                                          const std::vector<SegmentPiece> &pieces, double lengthScale,
                                          Normalize normalized, VisitPiece visitPiece,
                                          const char *errorPrefix,
                                          LinkGeometry geometry = LinkGeometry::Orthogonal)
{
    const auto require = [errorPrefix](bool condition, const char *reason) {
        if (!condition)
            throw std::runtime_error(std::string(errorPrefix) + reason);
    };
    const auto checked = [&require](long double value, const char *reason) {
        const double result = static_cast<double>(value);
        require(std::isfinite(value) && std::isfinite(result), reason);
        return result;
    };
    const auto checkPhase = [&require](Phase value) {
        require(value == Phase::Negative || value == Phase::Positive, "invalid phase tag");
    };
    const auto checkUnit = [&require](const Vector2D &value, const char *reason) {
        const double length = std::hypot(value.x, value.y);
        require(std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(length) &&
                    std::abs(length - 1.0) <= 128.0 * std::numeric_limits<double>::epsilon(),
                reason);
    };
    checkPhase(ownerPhase);
    checkPhase(neighborPhase);
    checkUnit(ownerFaceNormal, "owner face normal must have finite unit length");
    require(std::isfinite(faceLength) && faceLength > 0.0 && std::isfinite(linkDistance) &&
                linkDistance > 0.0,
            "face length and link distance must be finite and positive");
    require(std::isfinite(twoPointPermittivity) && twoPointPermittivity > 0.0,
            "two-point permittivity must be finite and positive");
    require(!pieces.empty(), "face partition must contain at least one positive-length piece");
    const auto ownerCoordinate = normalized(owner);
    const auto neighborCoordinate = normalized(neighbor);
    const double scaledLength =
        checked(static_cast<long double>(faceLength) / lengthScale, "unrepresentable normalized face length");
    const double scaledDistance = checked(static_cast<long double>(linkDistance) / lengthScale,
                                          "unrepresentable normalized link distance");
    require(scaledLength > 0.0 && scaledDistance > 0.0, "normalized face/link length underflow");
    const double tolerance = 256.0 * std::numeric_limits<double>::epsilon();
    const double linkScale =
        std::max({1.0, scaledDistance, std::abs(ownerCoordinate[0]), std::abs(ownerCoordinate[1]),
                  std::abs(neighborCoordinate[0]), std::abs(neighborCoordinate[1])});
    const double linkX = checked(static_cast<long double>(neighborCoordinate[0]) - ownerCoordinate[0],
                                 "non-finite normalized link x displacement");
    const double linkY = checked(static_cast<long double>(neighborCoordinate[1]) - ownerCoordinate[1],
                                 "non-finite normalized link y displacement");
    require(geometry == LinkGeometry::Orthogonal || geometry == LinkGeometry::ProjectedNormalCenters,
            "unknown center-link geometry");
    require(geometry == LinkGeometry::Orthogonal
                ? std::hypot(linkX - scaledDistance * ownerFaceNormal.x,
                             linkY - scaledDistance * ownerFaceNormal.y) <= tolerance * linkScale
                : std::abs(static_cast<long double>(linkX) * ownerFaceNormal.x +
                           static_cast<long double>(linkY) * ownerFaceNormal.y - scaledDistance) <=
                      tolerance * linkScale,
            "owner-to-neighbor link must have the supplied distance and outward normal direction");

    const auto visit = [&](const std::array<double, 2> &begin, const std::array<double, 2> &end,
                           const std::array<double, 2> &midpoint, double length, Phase phase,
                           const std::array<long double, 2> &span) {
        if constexpr (std::is_invocable_v<VisitPiece, decltype(begin), decltype(end), decltype(midpoint),
                                          double, Phase, decltype(span)>)
            visitPiece(begin, end, midpoint, length, phase, span);
        else
        {
            require(!pieces.front().retainedRange, "this trace does not consume retained phase widths");
            visitPiece(begin, end, midpoint, length, phase);
        }
    };

    if (pieces.front().retainedRange)
    {
        require(geometry == LinkGeometry::ProjectedNormalCenters,
                "retained native intervals require explicit projected-normal geometry");
        const auto a = normalized(pieces.front().begin), b = normalized(pieces.front().end);
        const std::array<long double, 2> whole{
            {static_cast<long double>(b[0]) - a[0], static_cast<long double>(b[1]) - a[1]}};
        const long double length = std::hypot(whole[0], whole[1]);
        const double scale =
            std::max({1.0, scaledLength, std::abs(a[0]), std::abs(a[1]), std::abs(b[0]), std::abs(b[1])});
        require(length > 0 && std::abs(length - scaledLength) <= tolerance * scale &&
                    std::abs(whole[0] * ownerFaceNormal.x + whole[1] * ownerFaceNormal.y) <=
                        tolerance * scale,
                "retained intervals do not belong to the supplied normal physical face");
        long double total = 0, previous = 0;
        for (const auto &piece : pieces)
        {
            checkPhase(piece.phase);
            require(piece.retainedRange && piece.begin.x == pieces.front().begin.x &&
                        piece.begin.y == pieces.front().begin.y && piece.end.x == pieces.front().end.x &&
                        piece.end.y == pieces.front().end.y && std::isfinite(piece.offset) &&
                        std::isfinite(piece.fraction) && piece.offset >= 0 && piece.fraction > 0 &&
                        piece.offset + piece.fraction <= 1 + tolerance,
                    "invalid retained face interval or source endpoints");
            const long double first = piece.fromEnd ? 1 - piece.offset - piece.fraction : piece.offset;
            require(std::abs(first - previous) <= tolerance,
                    "retained face intervals are not ordered/continuous");
            std::array<double, 2> midpoint{};
            std::array<long double, 2> span{};
            for (int axis = 0; axis < 2; ++axis)
            {
                span[axis] = whole[axis] * piece.fraction;
                midpoint[axis] = checked(
                    piece.fromEnd
                        ? static_cast<long double>(b[axis]) - whole[axis] * piece.offset - .5L * span[axis]
                        : static_cast<long double>(a[axis]) + whole[axis] * piece.offset + .5L * span[axis],
                    "unrepresentable retained face midpoint");
            }
            const double pieceLength =
                checked(length * piece.fraction, "unrepresentable retained phase length");
            require(pieceLength > 0, "positive retained phase length underflow");
            // The visitor uses the independent span for Gauss offsets. Do not
            // subtract two rounded endpoints, even when they coincide.
            visit(midpoint, midpoint, midpoint, pieceLength, piece.phase, span);
            total += piece.fraction;
            previous = piece.fromEnd ? 1 - piece.offset : piece.offset + piece.fraction;
        }
        require(std::abs(total - 1) <= tolerance && std::abs(previous - 1) <= tolerance,
                "retained phase lengths do not cover the complete face");
        const long double normalPosition =
            (.5L * a[0] + .5L * b[0] - ownerCoordinate[0]) * ownerFaceNormal.x +
            (.5L * a[1] + .5L * b[1] - ownerCoordinate[1]) * ownerFaceNormal.y;
        require(normalPosition >= -tolerance * scale && normalPosition <= scaledDistance + tolerance * scale,
                "native face plane does not separate the two cell centers");
        return {ownerCoordinate, neighborCoordinate};
    }

    const auto origin = normalized(pieces.front().begin);
    auto previousEnd = origin;
    std::array<double, 2> direction{{0.0, 0.0}};
    long double totalLength = 0.0L;
    double faceScale = std::max({1.0, scaledLength, std::abs(origin[0]), std::abs(origin[1])});
    for (std::size_t i = 0; i < pieces.size(); ++i)
    {
        checkPhase(pieces[i].phase);
        require(!pieces[i].retainedRange, "mixed endpoint and retained-width face representations");
        const auto begin = normalized(pieces[i].begin);
        const auto end = normalized(pieces[i].end);
        faceScale =
            std::max({faceScale, std::abs(begin[0]), std::abs(begin[1]), std::abs(end[0]), std::abs(end[1])});
        const double dx =
            checked(static_cast<long double>(end[0]) - begin[0], "non-finite piece x displacement");
        const double dy =
            checked(static_cast<long double>(end[1]) - begin[1], "non-finite piece y displacement");
        const double length = std::hypot(dx, dy);
        require(std::isfinite(length) && length > tolerance * faceScale,
                "face piece must have positive, roundoff-resolved length");
        require(std::hypot(begin[0] - previousEnd[0], begin[1] - previousEnd[1]) <= tolerance * faceScale,
                "face pieces must be continuous and ordered without gaps or overlaps");
        if (i == 0)
        {
            direction = {{dx / length, dy / length}};
            require(std::abs(static_cast<long double>(direction[0]) * ownerFaceNormal.x +
                             static_cast<long double>(direction[1]) * ownerFaceNormal.y) <= tolerance,
                    "face-piece line must be perpendicular to the owner face normal");
        }
        const long double along =
            static_cast<long double>(dx) * direction[0] + static_cast<long double>(dy) * direction[1];
        const long double cross =
            static_cast<long double>(dx) * direction[1] - static_cast<long double>(dy) * direction[0];
        const long double offset = (static_cast<long double>(begin[0]) - origin[0]) * direction[1] -
                                   (static_cast<long double>(begin[1]) - origin[1]) * direction[0];
        require(std::isfinite(along) && std::isfinite(cross) && std::isfinite(offset) && along > 0.0L &&
                    std::abs(cross) <= tolerance * faceScale && std::abs(offset) <= tolerance * faceScale,
                "face pieces must share one straight line without backtracking");
        const std::array<double, 2> midpoint{
            {checked(0.5L * begin[0] + 0.5L * end[0], "non-finite normalized piece midpoint x"),
             checked(0.5L * begin[1] + 0.5L * end[1], "non-finite normalized piece midpoint y")}};
        visit(begin, end, midpoint, length, pieces[i].phase,
              {{static_cast<long double>(end[0]) - begin[0], static_cast<long double>(end[1]) - begin[1]}});
        totalLength += length;
        previousEnd = end;
    }
    require(std::isfinite(totalLength) && std::abs(totalLength - scaledLength) <= tolerance * faceScale,
            "face-piece lengths must sum to the supplied face length");
    // The midpoint of the physical face must lie on the orthogonal center
    // link. Permit any position on that link, including either endpoint for
    // boundary values; do not impose the interior-face half-distance rule.
    const long double midpointX = 0.5L * origin[0] + 0.5L * previousEnd[0] - ownerCoordinate[0];
    const long double midpointY = 0.5L * origin[1] + 0.5L * previousEnd[1] - ownerCoordinate[1];
    const long double alongLink = midpointX * ownerFaceNormal.x + midpointY * ownerFaceNormal.y;
    const long double acrossLink = -midpointX * ownerFaceNormal.y + midpointY * ownerFaceNormal.x;
    const double positionTolerance = tolerance * std::max(linkScale, faceScale);
    require(
        std::isfinite(alongLink) && std::isfinite(acrossLink) &&
            (geometry == LinkGeometry::ProjectedNormalCenters || std::abs(acrossLink) <= positionTolerance) &&
            alongLink >= -positionTolerance && alongLink <= scaledDistance + positionTolerance,
        "complete face midpoint must lie on the owner-to-neighbor link segment");
    return {ownerCoordinate, neighborCoordinate};
}
} // namespace DielectricTraceFaceGeometry
} // namespace
} // namespace Quantum::Electricity::Detail
