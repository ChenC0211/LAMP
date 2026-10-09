#include "Modules/Interface/HermiteGalsVertexTransport.h"

#include "Modules/Interface/detail/HermiteGalsFirstJet.h"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace Quantum::Interface
{
using Detail::requireFirstJet;

HermiteGalsVertexTransportResult::HermiteGalsVertexTransportResult(
    HermiteGalsIdentity source, HermiteGalsTimeLayer sourceTime,
    const HermiteGalsVelocityView& arrival, std::vector<HermiteGalsVertexJet> vertices)
    : source_(source), arrival_(arrival.geometryIdentity()), sourceTime_(sourceTime),
      arrivalTime_(arrival.timeLayer()), vertices_(std::move(vertices)) {}

HermiteGalsVertexTransportResult transportHermiteGalsVertices(
    const QuadtreeMesh& mesh, const HermiteGalsState& source,
    const std::vector<int>& masterVertices, const LevelSetBoundaryTrace::Domain& domain,
    const HermiteGalsVelocityView& previousVelocity, const HermiteGalsVelocityView& arrivalVelocity)
{
    HermiteGalsIdentity sourceId;
    int vertex=-1;
    const char* stage="input";
    try {
        sourceId=source.identity();
        source.validateMesh(mesh);
        requireFirstJet(!masterVertices.empty(), "a nonempty master vertex selection is required");
        const auto sourceTime=source.timeLayer();
        const auto geometry=source.geometry();
        const auto& scalar=geometry.snapshot();
        Detail::validateFirstJetSource(mesh,source,scalar,domain,previousVelocity,arrivalVelocity);
        requireFirstJet(!scalar.hasCornerSectors(), "corner-sector geometry requires explicit one-sided vertex transport");
        const auto& sourceVertices=source.primary().vertices;
        int previous=-1;
        for(int index:masterVertices) {
            vertex=index;
            requireFirstJet(index>previous && index<mesh.vertIndexNum(), "target indices must be ascending, unique and native");
            requireFirstJet(mesh.verts(mesh.vertIndexToLabel[index]).type()==TVertType::non,
                "hanging vertex is not an independent master");
            const auto found=std::lower_bound(sourceVertices.begin(),sourceVertices.end(),index,
                [](const auto& jet,int id){return jet.vertex<id;});
            requireFirstJet(found!=sourceVertices.end() && found->vertex==index, "target is not a source master vertex");
            previous=index;
        }
        std::vector<HermiteGalsVertexJet> values;
        values.reserve(masterVertices.size());
        for(int index:masterVertices) {
            vertex=index;stage="trace";
            const auto point=mesh.verts(mesh.vertIndexToLabel[index]).coord();
            const auto jet=Detail::pullbackFirstJet(scalar,domain,point,previousVelocity,arrivalVelocity,stage);
            values.push_back({index,jet.phi,jet.gradient.x,jet.gradient.y});
        }
        // No partial vector/state is published if any vertex or provider fails.
        stage="final mesh binding";vertex=-1;
        source.validateMesh(mesh);
        return HermiteGalsVertexTransportResult(sourceId,sourceTime,arrivalVelocity,std::move(values));
    } catch(const std::exception& error) {
        throw std::runtime_error("[Interface/HermiteGalsVertexTransport source="+std::to_string(sourceId.value())
            +" arrival="+std::to_string(arrivalVelocity.geometryIdentity().value())+" vertex="+std::to_string(vertex)
            +"] "+stage+": "+error.what());
    }
}
} // namespace Quantum::Interface
