#pragma once

#include <cstddef>
#include <memory>
#include <vector>
#include "Core/Base/Vector2D.h"
#include "Core/Geometry/ActiveHermiteSensitivity.h"

namespace LevelSetHermiteGeometry { class Snapshot; }
class LevelSetMappedInterface;

// Retained straight-line query of the same cell-local Hermite polynomial,
// explicitly embedded between SHARED canonical cell bounds. Unlike the old
// nominal-h point map, adjacent cells have no exact gap/overlap. Native mesh
// coordinates and nodal DOFs are not changed. This is not a physical-C1 claim.
//
// Grid cuts and composed root positions are retained internally as positive-
// denominator ratios. No global root point is formed, and no short phase is
// snapped away. The original bounded root isolator still rejects unresolved
// roots; arithmetic outside the expansion range also rejects.
// The retained source must remain immutable through every other alias too.
// Cubic and EnrichedQuartic use the same phase/root/refinement algorithm,
// with the actual scalar's seven/nine exact restriction coefficients. Native
// negative-phase sensitivity has separate typed cubic and quartic data ports.
// Native Q4 phase measures also retain algebraically exact isolated multiple
// roots using the square-free polynomial of the SAME exact restriction. Phase
// signs still belong to the original scalar (even roots do not flip phase).
// This does not certify a regular 2D interface chart or a phase derivative;
// negative-phase variations on such charts reject. Positive-length zero
// traces and unresolved near-critical roots still reject.
class LevelSetMappedSegment
{
public:
    class NativeQuarticTrace;
    // Explicit geometry-only port. The normal is the POSITIVE coordinate
    // axis, independent of input traversal. Native Q4 and an axis-aligned
    // segment perpendicular to normalAxis are required; no phase convention
    // is supplied for material transport or legacy variations.
    static NativeQuarticTrace nativeQuarticTrace(
        std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        Vector2D begin,Vector2D end,int normalAxis);
    struct Range
    {
        // The original input segment has parameter t in [0,1]. When fromEnd
        // is false this interval is t=offset+s; otherwise t=1-offset-s,
        // with 0<=s<=length. length is always positive. Integration along the
        // INPUT direction uses |dt|; a fromEnd chart does not reverse the
        // physical co-normal or the input segment's orientation.
        //
        // Terminal floating approximations of a separately retained offset
        // and width, NOT two rounded endpoints. Use translated polynomial
        // integration; do not form offset+length and subtract primitives.
        // These numbers are not interval bounds or exact geometric measures.
        bool fromEnd = false;
        long double offset = 0.0L, length = 0.0L;
        bool negative = false;
    };
    struct Root
    {
        int cellI = -1, cellJ = -1;
        // Opaque retained chart identity, not a global/half-face parameter.
        // The representative and its certified enclosure belong to THAT
        // chart. A tiny distance from a grid-piece endpoint remains small.
        std::size_t chart = 0;
        double parameter = 0.0, lower = 0.0, upper = 0.0;
        bool pointOnly = false;
        // Native leaf provenance is separate from logical Cartesian i/j.
        // Exactly one addressing mode is used by a nonempty root inventory.
        int nativeCellIndex = -1;
    };

    LevelSetMappedSegment(std::shared_ptr<const LevelSetHermiteGeometry::Snapshot> source,
        Vector2D begin, Vector2D end);
    ~LevelSetMappedSegment();
    LevelSetMappedSegment(LevelSetMappedSegment&&) noexcept;
    LevelSetMappedSegment& operator=(LevelSetMappedSegment&&) noexcept;
    LevelSetMappedSegment(const LevelSetMappedSegment&) = delete;
    LevelSetMappedSegment& operator=(const LevelSetMappedSegment&) = delete;

    bool usesSource(const LevelSetHermiteGeometry::Snapshot& source) const;
    Vector2D begin() const;
    Vector2D end() const;
    const std::vector<Range>& ranges() const;
    const std::vector<Root>& roots() const;
    // Exact retained identity in the ORIGINAL input direction. Does not
    // recover t from a rounded physical point or snap a near-endpoint root.
    // rootIndex belongs to roots(); invalid indices and moved-from queries throw.
    bool rootAtInputEndpoint(std::size_t rootIndex, bool atEnd) const;
    std::size_t gridPieceCount() const;

    // Exact proportionality of the two restricted native-Q4 polynomials.
    // Same input segment and one positive-length native chart are required.
    // A true result does not assert equal amplitudes or equal phase signs.
    bool sameNativeQuarticTraceZeros(const LevelSetMappedSegment& other) const;
    bool sameNativeCubicTraceZeros(const LevelSetMappedSegment& other) const;
    struct CubicValueVariation
    {
        double value=0;
        std::vector<Quantum::Geometry::ActiveHermiteNodalSensitivity> independent;
    };
    CubicValueVariation nativeCubicValueVariation(bool fromEnd,long double offset,
        Quantum::Geometry::CornerSectorVariation sectors =
            Quantum::Geometry::CornerSectorVariation::RequireAbsent) const;
    struct QuarticValueVariation
    {
        double value=0;
        std::vector<Quantum::Geometry::ActiveQuarticHermiteSensitivity> independent;
    };
    // Arbitrary input-chart sample, not a root or a reconstructed global point.
    // Used for a space-time crossing coincident with a whole native face.
    QuarticValueVariation nativeQuarticValueVariation(bool fromEnd,long double offset,
        Quantum::Geometry::CornerSectorVariation sectors =
            Quantum::Geometry::CornerSectorVariation::RequireAbsent) const;

    // Refine only roots()[rootIndex] inside its certified enclosure. Returns
    // a new immutable query; the source and input segment do not change.
    // All same-point incidences and adjacent phase ranges use the new root.
    // The original query stays valid, including if refinement fails.
    // rootIndex belongs to this query's roots() inventory; obtain fresh
    // requests after replacement rather than caching stale cut references.
    LevelSetMappedSegment refined(std::size_t rootIndex) const;

    struct NegativePhaseVariation
    {
        // Original INPUT parameter t=offset or 1-offset; never a global point.
        bool fromEnd = false;
        long double offset = 0.0L;
        bool endpoint = false;
        // d integral H(-phi(t)) dt = sum coefficients * physical nodal DOFs.
        // Each distinct root contributes -B/|dphi/dt| exactly once. At an
        // exact input endpoint this returns the symmetric generalized slope
        // (one half). It is not a classical derivative of the clipped range;
        // adjacent halves recover the full interior-root contribution.
        std::vector<Quantum::Geometry::ActiveHermiteNodalSensitivity> independent;
    };
    // Actual native leaves only. Uses the retained root ratios, shared local
    // Hermite basis and hanging constraints. Interior simple transverse roots
    // give the ordinary first variation. Tangency, unresolved derivative and
    // sector-DOF variations reject unless declared jets are explicitly held
    // fixed. No velocity or physical measure is added.
    std::vector<NegativePhaseVariation> nativeNegativePhaseVariation(
        Quantum::Geometry::CornerSectorVariation sectors =
            Quantum::Geometry::CornerSectorVariation::RequireAbsent) const;

    struct QuarticNegativePhaseVariation
    {
        // Same input-root identity, endpoint half weight and derivative meaning
        // as NegativePhaseVariation. Columns also include independent edge
        // midpoint and cell-center data; no reduction to cubic nodal data.
        bool fromEnd = false;
        long double offset = 0.0L;
        bool endpoint = false;
        std::vector<Quantum::Geometry::ActiveQuarticHermiteSensitivity> independent;
    };
    std::vector<QuarticNegativePhaseVariation> nativeQuarticNegativePhaseVariation(
        Quantum::Geometry::CornerSectorVariation sectors =
            Quantum::Geometry::CornerSectorVariation::RequireAbsent) const;

private:
    friend class LevelSetMappedInterface;
    struct Impl;
    explicit LevelSetMappedSegment(std::unique_ptr<const Impl> impl);
    std::unique_ptr<const Impl> impl_;
    const Impl& checked() const;
};

// Same immutable native scalar, with value and positive-length zero traces
// explicitly separated. All terminal ranges use the ORIGINAL input parameter
// and retain offset/width separately, like Range above. They are approximations,
// not certified bounds; do not reconstruct root identity from rounded endpoints.
class LevelSetMappedSegment::NativeQuarticTrace
{
public:
    struct Interval
    {
        bool fromEnd=false;
        long double offset=0.0L,length=0.0L;
    };
    struct Side
    {
        bool available=false;
        int nativeCellIndex=-1,normalOrder=0;
        // Certified sign of the first nonzero normal coefficient, multiplied
        // by side^normalOrder. These ranges cover this coincident chart only.
        std::vector<Range> ranges;
    };
    struct CoincidentRange
    {
        std::size_t chart=0;
        // There is deliberately no negative flag on a zero interval.
        Interval interval;
        Side minus,plus;
    };
    struct Event
    {
        Root root;
        // 0=value root; otherwise a root of this same-source normal derivative.
        int normalOrder=0,normalSide=0;
        std::size_t coincidentChart=0;
    };
    ~NativeQuarticTrace();
    NativeQuarticTrace(NativeQuarticTrace&&) noexcept;
    NativeQuarticTrace& operator=(NativeQuarticTrace&&) noexcept;
    NativeQuarticTrace(const NativeQuarticTrace&)=delete;
    NativeQuarticTrace& operator=(const NativeQuarticTrace&)=delete;
    bool usesSource(const LevelSetHermiteGeometry::Snapshot& source) const;
    Vector2D begin() const;
    Vector2D end() const;
    int normalAxis() const;
    std::size_t gridPieceCount() const;
    const std::vector<Range>& ordinaryRanges() const;
    const std::vector<CoincidentRange>& coincidentRanges() const;
    // Value incidences precede independent side-event inventories. No global
    // root ordering or cross-side equality is implied by these indices.
    const std::vector<Event>& events() const;
    bool rootAtInputEndpoint(std::size_t eventIndex,bool atEnd) const;
    NativeQuarticTrace refined(std::size_t eventIndex) const;
    // Available only when there are NO coincident ranges. Legacy semantics,
    // source identity and variation behavior remain exactly the original.
    const LevelSetMappedSegment& ordinary() const;
private:
    friend class LevelSetMappedSegment;
    struct Impl;
    explicit NativeQuarticTrace(std::unique_ptr<const Impl> impl);
    std::unique_ptr<const Impl> impl_;
    const Impl& checked() const;
};
