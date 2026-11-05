#include "Core/Mesh/NativeMeshIdentity.h"
#include "Core/Mesh/QuadtreeMesh.h"

#include <cmath>
#include <exception>
#include <cstring>
#include <stdexcept>
#include <string>

namespace
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::invalid_argument(std::string("[NativeMeshIdentity] ") + message);
}

// Keep the unchanged failure construction out of the per-word comparison.
// The existing optimized binary called a separate verifier for every word.
[[noreturn]] void identityChanged()
{
    throw std::invalid_argument(std::string("[NativeMeshIdentity] ") +
        "mesh geometry/topology changed since identity construction");
}
#if defined(_MSC_VER)
#define LAMP_IDENTITY_WORD_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define LAMP_IDENTITY_WORD_INLINE inline __attribute__((always_inline))
#else
#define LAMP_IDENTITY_WORD_INLINE inline
#endif
class IdentityWordVerifier
{
public:
    explicit IdentityWordVerifier(const std::vector<std::uint64_t>& words):words_(words){}
    LAMP_IDENTITY_WORD_INLINE void operator()(std::uint64_t word)
    {
        if(next_>=words_.size()||word!=words_[next_])identityChanged();
        ++next_;
    }
    std::size_t count() const noexcept{return next_;}
private:
    const std::vector<std::uint64_t>& words_;
    std::size_t next_=0;
};
#undef LAMP_IDENTITY_WORD_INLINE

// Exact ordered words, not a fingerprint. Exclude phase/cut flags and the
// mutable supCellIndex solver cache; retain connectivity used by native stencils.
template<class Emit> void visitIdentityWords(const QuadtreeMesh& mesh,Emit&& emit)
{
    require(mesh.cellIndexNum() > 0 && mesh.vertIndexNum() > 0 &&
            mesh.vEdgeIndexNum() > 0 && mesh.hEdgeIndexNum() > 0 &&
            mesh.cellIndexToLabel && mesh.vertIndexToLabel &&
            mesh.vEdgeIndexToLabel && mesh.hEdgeIndexToLabel &&
            mesh.cellLabelCapacity() >= mesh.cellIndexNum() &&
            mesh.vertLabelCapacity() >= mesh.vertIndexNum() &&
            mesh.vEdgeLabelCapacity() >= mesh.vEdgeIndexNum() &&
            mesh.hEdgeLabelCapacity() >= mesh.hEdgeIndexNum(), "missing or undersized indexed mesh storage");
    const auto integer = [&](int value) { emit(static_cast<std::uint64_t>(static_cast<std::int64_t>(value))); };
    const auto scalar = [&](double value) {
        require(std::isfinite(value), "non-finite mesh coordinate/metric");
        std::uint64_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "binary64 mesh identity requires 64-bit double");
        std::memcpy(&bits, &value, sizeof(bits));
        emit(bits);
    };
    const auto point = [&](Vector2D value) { scalar(value.x); scalar(value.y); };
    integer(mesh.cellIndexNum()); integer(mesh.vertIndexNum());
    integer(mesh.vEdgeIndexNum()); integer(mesh.hEdgeIndexNum());
    integer(mesh.Nx); integer(mesh.Ny);
    point(mesh.minCoord); point(mesh.maxCoord);
    bool mixedLevels=false;
    for (int i = 0; i < mesh.cellIndexNum(); ++i)
    {
        const int label = mesh.cellIndexToLabel[i];
        require(label >= 0 && label < mesh.cellLabelCapacity(), "invalid indexed cell label");
        const auto& cell = mesh.cells(label);
        integer(label); integer(cell.label()); integer(cell.index());
        integer(cell.level()); integer(cell.parent); integer(static_cast<int>(cell.boundary()));
        mixedLevels |= cell.level()!=mesh.cells(mesh.cellIndexToLabel[0]).level();
        point(cell.centerCoord()); point(cell.size());
        for (int son : cell.son) integer(son);
        for (int vertex : cell.cellToVert) integer(vertex);
        for (const auto& side : cell.cellToVEdge) for (int edge : side) integer(edge);
        for (const auto& side : cell.cellToHEdge) for (int edge : side) integer(edge);
        for (const auto& side : cell.neiCell) for (int neighbor : side) integer(neighbor);
    }
    // Adaptive RZ rows follow parent.son[] to native sibling faces. Leaf words
    // alone cannot identify that operator after an ancestor-only mutation.
    if (mixedLevels)
    {
        std::vector<unsigned char> seen(mesh.cellLabelCapacity(),0);
        for (int i=0;i<mesh.cellIndexNum();++i)
        {
            int label=mesh.cells(mesh.cellIndexToLabel[i]).parent;
            while (label>=0)
            {
                require(label<mesh.cellLabelCapacity(), "invalid native ancestor label");
                if (seen[label]) break;
                seen[label]=1;
                const auto& cell=mesh.cells(label);
                integer(label); integer(cell.label()); integer(cell.level()); integer(cell.parent);
                point(cell.centerCoord()); point(cell.size());
                for (int son:cell.son) integer(son);
                for (int vertex:cell.cellToVert) integer(vertex);
                for (const auto& side:cell.cellToVEdge) for (int edge:side) integer(edge);
                for (const auto& side:cell.cellToHEdge) for (int edge:side) integer(edge);
                label=cell.parent;
            }
        }
    }
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const int label = mesh.vertIndexToLabel[i];
        require(label >= 0 && label < mesh.vertLabelCapacity(), "invalid indexed vertex label");
        const auto& vertex = mesh.verts(label);
        integer(label); integer(vertex.label()); integer(vertex.index()); point(vertex.coord());
        integer(vertex.parent); integer(vertex.son[0]); integer(static_cast<int>(vertex.boundary()));
        integer(static_cast<int>(vertex.type()));
        integer(vertex.supCell);
        for (int value : vertex.supVert) integer(value);
        for (int value : vertex.vertToCell) integer(value);
        for (int value : vertex.vertToVEdge) integer(value);
        for (int value : vertex.vertToHEdge) integer(value);
        for (int value : vertex.neiVertices) integer(value);
    }
    for (bool horizontal : {false, true})
    {
        const int count = horizontal ? mesh.hEdgeIndexNum() : mesh.vEdgeIndexNum();
        const int capacity = horizontal ? mesh.hEdgeLabelCapacity() : mesh.vEdgeLabelCapacity();
        const auto labels = horizontal ? mesh.hEdgeIndexToLabel : mesh.vEdgeIndexToLabel;
        for (int i = 0; i < count; ++i)
        {
            const int label = labels[i];
            require(label >= 0 && label < capacity, "invalid indexed edge label");
            const auto& edge = horizontal ? mesh.hEdges(label) : mesh.vEdges(label);
            integer(label); integer(edge.label()); integer(edge.index()); integer(edge.parent);
            integer(edge.son[0]); integer(edge.son[1]); integer(edge.masterCell); integer(edge.slaveCell);
            integer(edge.edgeToVert[0]); integer(edge.edgeToVert[1]);
            integer(edge.isHorizontal() ? 1 : 0); integer(static_cast<int>(edge.boundary()));
            point(edge.midpointCoord()); scalar(edge.length());
        }
    }
}

} // namespace

thread_local NativeMeshIdentity::ValidationBatch* NativeMeshIdentity::ValidationBatch::active_=nullptr;

NativeMeshIdentity::ValidationBatch::ValidationBatch(const QuadtreeMesh& mesh) noexcept
    :mesh_(&mesh),previous_(active_),uncaught_(std::uncaught_exceptions())
{active_=this;}
NativeMeshIdentity::ValidationBatch::~ValidationBatch() noexcept
{
    if(std::uncaught_exceptions()>uncaught_)invalidateChain();
    active_=previous_;
}
void NativeMeshIdentity::ValidationBatch::invalidateChain() noexcept
{
    for(auto* batch=this;batch;batch=batch->previous_){
        while(batch->size_)batch->proofs_[--batch->size_].reset();
    }
}
bool NativeMeshIdentity::ValidationBatch::contains(const std::shared_ptr<const unsigned char>& token) const noexcept
{
    if(!token)return false;
    for(std::size_t i=0;i<size_;++i)if(proofs_[i]==token)return true;
    return false;
}
void NativeMeshIdentity::ValidationBatch::remember(const std::shared_ptr<const unsigned char>& token) noexcept
{
    // Capacity exhaustion keeps the original complete scan, never a false hit.
    if(token&&size_<proofs_.size())proofs_[size_++]=token;
}
NativeMeshIdentity::ValidationBatch::Pause::Pause() noexcept:suspended_(active_)
{
    if(suspended_)suspended_->invalidateChain();
    active_=nullptr;
}
NativeMeshIdentity::ValidationBatch::Pause::~Pause() noexcept
{active_=suspended_;}

NativeMeshIdentity::NativeMeshIdentity(const QuadtreeMesh& mesh)
{
    visitIdentityWords(mesh,[&](std::uint64_t word){words_.push_back(word);});
    token_=std::make_shared<const unsigned char>(0);
}

void NativeMeshIdentity::validate(const QuadtreeMesh& mesh) const
{
    ValidationBatch* destination=nullptr;
    for(auto* batch=ValidationBatch::active_;batch;batch=batch->previous_)if(batch->mesh_==&mesh){
        if(batch->contains(token_))return;
        if(!destination)destination=batch;
    }
    // Public calls outside an audited read-only segment always scan all live
    // words; callbacks pause and invalidate every enclosing proof first.
    IdentityWordVerifier verifier(words_);
    visitIdentityWords(mesh,verifier);
    require(verifier.count()==words_.size(),"mesh geometry/topology changed since identity construction");
    if(destination)destination->remember(token_);
}

