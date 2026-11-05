#pragma once

#include "Modules/Interface/HermiteGalsPublication.h"
#include <cstddef>
#include <memory>

namespace LAMP::Interface
{
enum class HermiteGalsMeshTransferKind { ExactRefinement, FirstJetProjection };
struct HermiteGalsMeshTransferReport
{
    HermiteGalsMeshTransferKind kind = HermiteGalsMeshTransferKind::ExactRefinement;
    std::size_t refinedCells = 0, coarsenedCells = 0, unchangedCells = 0;
    long double sourceVolume = 0, targetVolume = 0, relativeDomainVolumeChange = 0;
};
class HermiteGalsMeshTransfer
{
public:
    HermiteGalsIdentity sourceIdentity() const;
    const HermiteGalsPublication& publication() const;
    const HermiteGalsMeshTransferReport& report() const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit HermiteGalsMeshTransfer(std::shared_ptr<const Data>);
    const Data& checked() const;
    static HermiteGalsMeshTransfer prepare(const QuadtreeMesh&, const HermiteGalsPublication&,
        const QuadtreeMesh&, HermiteGalsMeshTransferKind, const HermiteGalsPublicationConfig&);
    friend HermiteGalsMeshTransfer transferHermiteGalsFirstJetMesh(const QuadtreeMesh&,
        const HermiteGalsPublication&, const QuadtreeMesh&, HermiteGalsMeshTransferKind,
        const HermiteGalsPublicationConfig&);
};
// Zero physical time, new revision. Source and both external meshes remain
// unchanged. ExactRefinement never projects on failure; FirstJetProjection
// is explicit interpolation, with no claim of general Q4 equality or volume
// conservation. Neither operation supplies an App/Coupling accepted commit.
HermiteGalsMeshTransfer transferHermiteGalsFirstJetMesh(const QuadtreeMesh& sourceMesh,
    const HermiteGalsPublication& source, const QuadtreeMesh& targetMesh,
    HermiteGalsMeshTransferKind kind, const HermiteGalsPublicationConfig& config = {});
}
