#pragma once

#include <array>
#include "Core/Base/Vector2D.h"
#include "Core/Base/TimeLayer.h"
#include <cstdint>
#include <memory>
#include <vector>

class QuadtreeMesh;
namespace LevelSetHermiteGeometry { class Snapshot; }

namespace LAMP::Interface
{
namespace Detail { struct HermiteGalsData; }
class HermiteGalsBubbleResult;
class HermiteGalsPublication;
class HermiteGalsFirstJetLocalPreparation;
struct HermiteGalsPublicationConfig;

// Process-local identity of ONE immutable primary revision. Zero is unbound.
// It is not a restart ID or an accepted-step counter.
class HermiteGalsIdentity
{
public:
    HermiteGalsIdentity() = default;
    std::uint64_t value() const noexcept { return value_; }
    friend bool operator==(HermiteGalsIdentity a, HermiteGalsIdentity b) noexcept
    { return a.value_ == b.value_; }
    friend bool operator!=(HermiteGalsIdentity a, HermiteGalsIdentity b) noexcept
    { return !(a == b); }
private:
    explicit HermiteGalsIdentity(std::uint64_t value) : value_(value) {}
    std::uint64_t value_ = 0;
    friend class HermiteGalsCandidate;
};

// Compatibility name; the clock value has no Hermite state owner.
using HermiteGalsTimeLayer = LAMP::TimeLayer;

enum class HermiteGalsMeasure { CartesianArea, AxisymmetricVolume }; // dx dy / 2*pi*r dr dz
enum class HermiteGalsTransport { ConservativeVolume, FirstJet };
enum class HermiteGalsGradientUpdate { CharacteristicPullback, CompatibleIncrement };

struct HermiteGalsVertexJet
{
    int vertex = -1; // native index, not label; independent master only
    double phi = 0.0, dx = 0.0, dy = 0.0; // physical derivatives
};
struct HermiteGalsEdgeJet
{
    // Full master leaf side, including a coarse side spanning two fine faces.
    // Endpoints ordered in increasing x (horizontal) or y (vertical).
    // dn uses +y on horizontal sides and +x on vertical sides, NOT outward n_K.
    std::array<int, 2> vertices{{-1, -1}};
    double phi = 0.0, dn = 0.0; // at the exact parametric midpoint
};
struct HermiteGalsCellVolume
{
    int cell = -1; // native index, including pure-phase inventory cells
    double liquid = 0.0; // M_K in conservative state; geometric volume in a measurement result
};
struct HermiteGalsPrimaryState
{
    HermiteGalsGradientUpdate gradientUpdate = HermiteGalsGradientUpdate::CharacteristicPullback;
    HermiteGalsMeasure measure = HermiteGalsMeasure::CartesianArea;
    std::vector<HermiteGalsVertexJet> vertices; // ascending native indices
    std::vector<HermiteGalsEdgeJet> edges; // unique ordered endpoint pairs
    std::vector<int> geometryCells; // ascending active polynomial cells
    std::vector<HermiteGalsCellVolume> phaseVolumes; // conservative: covers geometryCells; FirstJet: empty
    HermiteGalsTransport transport = HermiteGalsTransport::ConservativeVolume; // explicit historical default
};

// Derived data, deliberately absent from HermiteGalsPrimaryState.
// I05/I06 produce dxy; c_K enters phi_boundary + c_K*q4(xi)*q4(eta).
// ConservativeVolume uses I10; FirstJet uses the actual traced center value.
// Bubble coefficients depend on the chosen boundary extension. I10/I11 use
// the zero-center cardinal extension, where c equals the center sample.
struct HermiteGalsClosure
{
    std::vector<double> vertexMixedDerivatives; // same order as master vertices
    std::vector<double> cellBubbleCoefficients; // same order as geometryCells
};

// An owning, read-only geometry view. Identity and time cannot be detached
// accidentally by obtaining a naked shared_ptr from this API.
class HermiteGalsGeometryView
{
public:
    HermiteGalsIdentity identity() const;
    HermiteGalsTimeLayer timeLayer() const;
    const LevelSetHermiteGeometry::Snapshot& snapshot() const;
private:
    explicit HermiteGalsGeometryView(std::shared_ptr<const Detail::HermiteGalsData> data);
    std::shared_ptr<const Detail::HermiteGalsData> data_;
    friend class HermiteGalsState;
};

// Immutable completed module candidate. Copying retains the SAME authority.
// App/Coupling owns the accepted handle and replaces it only after all module
// checks pass. finish() below is NOT a coupled commit or an acceptance test.
class HermiteGalsState
{
public:
    HermiteGalsIdentity identity() const;
    HermiteGalsTimeLayer timeLayer() const;
    const HermiteGalsPrimaryState& primary() const;
    const HermiteGalsClosure& closure() const;
    HermiteGalsGeometryView geometry() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
private:
    explicit HermiteGalsState(std::shared_ptr<const Detail::HermiteGalsData> data);
    std::shared_ptr<const Detail::HermiteGalsData> data_;
    friend class HermiteGalsCandidate;
};

class HermiteGalsCandidate
{
public:
    // Copies all input arrays; no writable aliases into retained state.
    HermiteGalsCandidate(const QuadtreeMesh& mesh, HermiteGalsTimeLayer time,
        const HermiteGalsPrimaryState& primary);
    ~HermiteGalsCandidate();
    HermiteGalsCandidate(HermiteGalsCandidate&&) noexcept;
    HermiteGalsCandidate& operator=(HermiteGalsCandidate&&) noexcept;
    HermiteGalsCandidate(const HermiteGalsCandidate&) = delete;
    HermiteGalsCandidate& operator=(const HermiteGalsCandidate&) = delete;

    HermiteGalsIdentity identity() const;
    HermiteGalsTimeLayer timeLayer() const;
    const HermiteGalsPrimaryState& primary() const;
    bool hasGeometry() const;
    void validateMesh(const QuadtreeMesh& mesh) const;
    // Success creates a fresh identity and invalidates ALL derived data.
    // Failure preserves this candidate, including its identity and geometry.
    void replacePrimary(const QuadtreeMesh& mesh, const HermiteGalsPrimaryState& primary);

    // Structural ownership port, not the checked I11 production publication
    // generator. source must identify the unchanged primary used to compute
    // closure and snapshot. Checks structural binding, not polynomial/volume
    // consistency: that numerical prerequisite belongs to I05-I11.
    // Snapshot is moved into a const allocation, never adopted via an externally
    // writable shared_ptr. One binding per primary revision; no fallback.
    void bindGeometry(const QuadtreeMesh& mesh, HermiteGalsIdentity source,
        const HermiteGalsClosure& closure, LevelSetHermiteGeometry::Snapshot&& snapshot);
    HermiteGalsState finish() &&; // requires binding; consumes the candidate
private:
    explicit HermiteGalsCandidate(std::unique_ptr<Detail::HermiteGalsData>) noexcept;
    static HermiteGalsState bindPreparedLocal(const HermiteGalsState&,HermiteGalsTimeLayer,
        const HermiteGalsPrimaryState&,const HermiteGalsClosure&,LevelSetHermiteGeometry::Snapshot&&,
        const std::vector<int>&,const std::vector<Vector2D>&,const std::vector<unsigned char>&);
    friend class HermiteGalsFirstJetLocalPreparation;
    // Private unpublished transaction branch; no second public identity/owner.
    HermiteGalsCandidate publicationCopy() const;
    void consumePublication() noexcept;
    std::unique_ptr<Detail::HermiteGalsData> data_;
    friend class HermiteGalsPublication;
    friend HermiteGalsPublication publishHermiteGalsGeometry(const QuadtreeMesh&,
        HermiteGalsCandidate&&, const HermiteGalsBubbleResult&, const HermiteGalsPublicationConfig&);
};
} // namespace LAMP::Interface
