#include "Core/IO/FileManagement.h"
#include "Core/IO/OutputLayout.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/DyadicGridCoordinates.h"
#include "Core/Mesh/MeshCacheIO.h"
#include "Core/Mesh/detail/QuadtreeMeshCacheSnapshot.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <cfenv>
#include <cerrno>
#include <exception>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

using namespace std;

namespace
{
// Optional integer-only execution storage must not change the caller's FP/errno.
class NativeMaintenanceEnvironment
{
    std::fenv_t native_{};int error_=errno;bool restore_=true;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr_=_mm_getcsr();
#endif
public:
    bool captured=std::fegetenv(&native_)==0;
    void dismiss() noexcept {restore_=false;}
    ~NativeMaintenanceEnvironment() noexcept
    {
        if(!restore_)return;
        if(captured&&std::fesetenv(&native_)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr_);
#endif
        errno=error_;
    }
};
}
struct QuadtreeMesh::NativeAdaptWorkspace
{
    NativeAdaptationMaintenanceReport* report=nullptr;bool valid=true,transferStarted=false;
    template<class T>struct Pool
    {
        LAMP::Mesh::Detail::EntityPool<T>* records=nullptr;const void* recordOwner=nullptr;int maximum=0,leafCount=0,freeCount=0,dirtyCount=0,priorCount=0,lastLive=-1;
        std::unique_ptr<unsigned char[]> live,leaf,marked;
        std::unique_ptr<int[]> leaves,free,dirty,scratch,prior;
        void allocate(LAMP::Mesh::Detail::EntityPool<T>& objects,int capacity)
        {
            records=&objects;recordOwner=objects.identity();maximum=capacity;const auto n=static_cast<std::size_t>(capacity);
            live=std::make_unique<unsigned char[]>(n);leaf=std::make_unique<unsigned char[]>(n);
            marked=std::make_unique<unsigned char[]>(n);
            leaves=std::unique_ptr<int[]>(new int[n]);free=std::unique_ptr<int[]>(new int[n]);
            dirty=std::unique_ptr<int[]>(new int[n]);scratch=std::unique_ptr<int[]>(new int[n]);prior=std::unique_ptr<int[]>(new int[n]);
        }
        bool inspect(int ptr,int count,const int* indices,const int* labels) noexcept
        {
            if(!records||recordOwner!=records->identity()||records->capacity()!=maximum||!labels||ptr<0||ptr>maximum||count<0||count>maximum||(count&&!indices))return false;
            for(int i=0;i<maximum;++i){const auto& t=records->peek(i);
                if(t.label()!=-1){if(t.label()!=i||i>=ptr)return false;live[i]=1;lastLive=i;
                    if(t.son[0]==-1){if(leafCount>=count||indices[leafCount]!=i||t.index()!=leafCount)return false;
                        leaf[i]=1;leaves[leafCount++]=i;}}
                else{if(labels[freeCount]!=i)return false;free[freeCount++]=i;}}
            return leafCount==count;
        }
        bool touch(int i,NativeAdaptationMaintenanceReport& r) noexcept
        {
            if(i<0||i>=maximum)return false;
            if(!marked[i]){if(dirtyCount>=maximum)return false;marked[i]=1;dirty[dirtyCount++]=i;++r.journaledSlots;}
            return true;
        }
        int labelOf(const void* object) const noexcept
        {return records&&recordOwner==records->identity()?records->labelOf(object):-1;}
        void merge(bool forLeaf)
        {
            auto& list=forLeaf?leaves:free;auto& count=forLeaf?leafCount:freeCount;
            int i=0,j=0,out=0;
            while(i<count||j<dirtyCount){
                const int a=i<count?list[i]:maximum,b=j<dirtyCount?dirty[j]:maximum,key=std::min(a,b);
                const bool changed=b==key;
                const bool keep=changed?(forLeaf?leaf[key]!=0:live[key]==0):true;
                if(keep)scratch[out++]=key;
                if(a==key)++i;if(b==key)++j;
            }
            list.swap(scratch);count=out;
        }
        bool refresh(int ptr,NativeAdaptationMaintenanceReport& r)
        {
            if(ptr<0||ptr>maximum)return false;
            if(dirtyCount==0)return lastLive<ptr;
            std::sort(dirty.get(),dirty.get()+dirtyCount);
            for(int k=0;k<dirtyCount;++k){const int i=dirty[k];const auto& t=records->peek(i);
                if(t.label()!=-1&&t.label()!=i)return false;
                const bool nextLive=t.label()!=-1,nextLeaf=nextLive&&t.son[0]==-1;
                if(nextLive!=(live[i]!=0)||nextLeaf!=(leaf[i]!=0))++r.changedMembershipSlots;
                live[i]=static_cast<unsigned char>(nextLive);leaf[i]=static_cast<unsigned char>(nextLeaf);
                if(nextLive)lastLive=std::max(lastLive,i);
            }
            while(lastLive>=0&&!live[lastLive])--lastLive;
            if(lastLive>=ptr)return false; // Literal statistics handles an unproven binding.
            merge(true);merge(false);return true;
        }
        void exportIndices(int& count,int*& indices,NativeAdaptationMaintenanceReport& r)
        {
            // Same allocation/publication order as the original count().
            count=leafCount;int* next=new int[count];
            for(int i=0;i<count;++i){next[i]=leaves[i];(*records)[leaves[i]].setIndex(i);}
            delete[] indices;indices=next;r.exportedIndexEntries+=static_cast<std::size_t>(count);r.actuallyUsed=true;
        }
        void exportFree(int& ptr,int* labels,NativeAdaptationMaintenanceReport& r) noexcept
        {
            std::copy_n(free.get(),freeCount,labels);if(lastLive>=0)ptr=lastLive+1;
            r.exportedFreeEntries+=static_cast<std::size_t>(freeCount);r.actuallyUsed=true;
        }
        void rememberIndices() noexcept
        {std::copy_n(leaves.get(),leafCount,prior.get());priorCount=leafCount;}
        void resetJournal() noexcept
        {for(int k=0;k<dirtyCount;++k)marked[dirty[k]]=0;dirtyCount=0;}
        void clearTransfer(NativeAdaptationMaintenanceReport& r) noexcept
        {
            if constexpr(!std::is_same_v<T,Vert>){
                const auto clear=[&](int i){auto& t=(*records)[i];t.setPreIndex(-1);t.relatedPreIndex.clear();t.relatedPreCoord.clear();++r.clearedTransferSlots;};
                for(int k=0;k<priorCount;++k){const int i=prior[k];clear(i);marked[i]=0;}
                for(int k=0;k<dirtyCount;++k){const int i=dirty[k];if(marked[i]){clear(i);marked[i]=0;}}
                for(int i=0;i<leafCount;++i)(*records)[leaves[i]].setPreIndex(i);
            }
            resetJournal();rememberIndices();
        }
        void release() noexcept
        {live.reset();leaf.reset();marked.reset();leaves.reset();free.reset();dirty.reset();scratch.reset();prior.reset();}
    };
    Pool<Vert> vertices;Pool<Edge> horizontal,vertical;Pool<Cell> cells;
    ~NativeAdaptWorkspace() noexcept
    {
        NativeMaintenanceEnvironment environment;
        vertices.release();horizontal.release();vertical.release();cells.release();
    }
    void fail() noexcept
    {
        valid=false;if(report)report->fallback=true;
        NativeMaintenanceEnvironment environment;
        vertices.release();horizontal.release();vertical.release();cells.release();
    }
};
std::shared_ptr<QuadtreeMesh::NativeAdaptWorkspace> QuadtreeMesh::beginNativeAdaptationMaintenance(bool requested) noexcept
{
    nativeAdaptationReport_=NativeAdaptationMaintenanceReport{};nativeAdaptationReport_.requested=requested;
    if(!requested)return {};
    NativeMaintenanceEnvironment environment;
    if(!environment.captured){nativeAdaptationReport_.fallback=true;return {};}
    try{
        const int limits[]{maxVertNum_,maxHEdgeNum_,maxVEdgeNum_,maxCellNum_};
        std::size_t bytes=sizeof(NativeAdaptWorkspace)+256u+4u*8u*64u;
        constexpr std::size_t bytesPerSlot=3u*sizeof(unsigned char)+5u*sizeof(int);
        if(bytes>nativeAdaptationReport_.workspaceLimitBytes){nativeAdaptationReport_.fallback=true;return {};}
        // Eligibility includes required storage-owned address tables conservatively;
        // these are shared persistent stock, not another workspace allocation.
        const std::size_t addressMetadata[]{verts_.addressMetadataBytes(),hEdges_.addressMetadataBytes(),
            vEdges_.addressMetadataBytes(),cells_.addressMetadataBytes()};
        for(const auto n:addressMetadata){
            if(n>nativeAdaptationReport_.workspaceLimitBytes-bytes){nativeAdaptationReport_.fallback=true;return {};}
            bytes+=n;
        }
        for(int n:limits){if(n<=0||static_cast<std::size_t>(n)>(nativeAdaptationReport_.workspaceLimitBytes-bytes)/bytesPerSlot){nativeAdaptationReport_.fallback=true;return {};}
            bytes+=static_cast<std::size_t>(n)*bytesPerSlot;}
        auto w=std::make_shared<NativeAdaptWorkspace>();w->report=&nativeAdaptationReport_;
        w->vertices.allocate(verts_,maxVertNum_);w->horizontal.allocate(hEdges_,maxHEdgeNum_);
        w->vertical.allocate(vEdges_,maxVEdgeNum_);w->cells.allocate(cells_,maxCellNum_);
        if(!w->vertices.inspect(vertPtr_,vertIndexNum_,vertIndexToLabel,vertLabelSpace_)
            ||!w->horizontal.inspect(hEdgePtr_,hEdgeIndexNum_,hEdgeIndexToLabel,hEdgeLabelSpace_)
            ||!w->vertical.inspect(vEdgePtr_,vEdgeIndexNum_,vEdgeIndexToLabel,vEdgeLabelSpace_)
            ||!w->cells.inspect(cellPtr_,cellIndexNum_,cellIndexToLabel,cellLabelSpace_)){
            nativeAdaptationReport_.fallback=true;return {};}
        nativeAdaptationReport_.eligible=true;nativeAdaptationReport_.workspaceBytes=bytes;return w;
    }catch(...){nativeAdaptationReport_.fallback=true;return {};}
}
const NativeAdaptationMaintenanceReport& QuadtreeMesh::nativeAdaptationMaintenanceReport() const noexcept
{return nativeAdaptationReport_;}
std::array<NativeEntityPoolStorageReport,4> QuadtreeMesh::entityPoolStorageReport() const noexcept
{
    const auto report=[](const auto& pool){
        using Entity=typename std::decay_t<decltype(pool)>::value_type;
        const auto stock=pool.storageStatistics();
        return NativeEntityPoolStorageReport{stock[0],stock[1],stock[2],stock[1]*sizeof(Entity),pool.denseMode()};
    };
    return {report(verts_),report(hEdges_),report(vEdges_),report(cells_)};
}
void QuadtreeMesh::recordNativeWrite(NativeAdaptWorkspace* w,Cell* t) noexcept
{if(w&&w->valid){const int label=w->cells.labelOf(t);if(label<0||!w->cells.touch(label,*w->report))w->fail();}}
void QuadtreeMesh::recordNativeWrite(NativeAdaptWorkspace* w,Vert* t) noexcept
{if(w&&w->valid){const int label=w->vertices.labelOf(t);if(label<0||!w->vertices.touch(label,*w->report))w->fail();}}
void QuadtreeMesh::recordNativeWrite(NativeAdaptWorkspace* w,Edge* t) noexcept
{
    if(!w||!w->valid)return;
    const int h=w->horizontal.labelOf(t);
    if(h>=0){if(!w->horizontal.touch(h,*w->report))w->fail();}
    else{const int v=w->vertical.labelOf(t);if(v<0||!w->vertical.touch(v,*w->report))w->fail();}
}
void QuadtreeMesh::recordNativeCellEdges(NativeAdaptWorkspace* w,const Cell& c) noexcept
{
    if(!w||!w->valid)return;
    for(int side=0;side<2;++side)for(int half=0;half<2;++half){const int h=c.cellToHEdge[side][half],v=c.cellToVEdge[side][half];
        if(h>=0){if(h>=maxHEdgeNum_)w->fail();else recordNativeWrite(w,&hEdges_[h]);}
        if(v>=0){if(v>=maxVEdgeNum_)w->fail();else recordNativeWrite(w,&vEdges_[v]);}}
}
void QuadtreeMesh::prepareTransferIndices(NativeAdaptWorkspace* w)
{
    if(!w||!w->valid){prepareTransferIndices();return;}
    if(!w->transferStarted){
        prepareTransferIndices();w->vertices.rememberIndices();w->horizontal.rememberIndices();w->vertical.rememberIndices();w->cells.rememberIndices();
        w->transferStarted=true;nativeAdaptationReport_.clearedTransferSlots+=cells_.materializedCount()+hEdges_.materializedCount()+vEdges_.materializedCount();return;
    }
    w->vertices.clearTransfer(nativeAdaptationReport_);w->horizontal.clearTransfer(nativeAdaptationReport_);
    w->vertical.clearTransfer(nativeAdaptationReport_);w->cells.clearTransfer(nativeAdaptationReport_);
    nativeAdaptationReport_.actuallyUsed=true;
}
bool QuadtreeMesh::updateNativeLabelSpaces(NativeAdaptWorkspace* w) noexcept
{
    if(!w||!w->valid)return false;
    w->vertices.exportFree(vertPtr_,vertLabelSpace_,nativeAdaptationReport_);
    w->horizontal.exportFree(hEdgePtr_,hEdgeLabelSpace_,nativeAdaptationReport_);
    w->vertical.exportFree(vEdgePtr_,vEdgeLabelSpace_,nativeAdaptationReport_);
    w->cells.exportFree(cellPtr_,cellLabelSpace_,nativeAdaptationReport_);return true;
}
bool QuadtreeMesh::refreshNativeMembership(NativeAdaptWorkspace* w)
{
    if(!w||!w->valid)return false;
    if(!w->vertices.refresh(vertPtr_,nativeAdaptationReport_)
        ||!w->horizontal.refresh(hEdgePtr_,nativeAdaptationReport_)
        ||!w->vertical.refresh(vEdgePtr_,nativeAdaptationReport_)
        ||!w->cells.refresh(cellPtr_,nativeAdaptationReport_)){w->fail();return false;}
    return true;
}


QuadtreeMesh::QuadtreeMesh(const DictionaryMesh& dictMesh) : 
    BackgroundMesh(dictMesh),
    maxVertNum_(dictMesh.maxCellNum),
    maxHEdgeNum_(dictMesh.maxCellNum),
    maxVEdgeNum_(dictMesh.maxCellNum),
    maxCellNum_(dictMesh.maxCellNum),
    maxInitialIterations_(dictMesh.maxInitialIterations),
    C_Lipschitz(dictMesh.C_Lipschitz)
{
    // BackgroundMesh already preserved the explicit zero-storage constructor opt-in.
    if (dictMesh.levelSetCount != 0) setLevelSetCount(dictMesh.levelSetCount);
    minLevel_neg = dictMesh.minLevel_neg;
    minLevel_pos = dictMesh.minLevel_pos;
    maxLevel_local = dictMesh.maxLevel_local;
    const std::size_t lsCount = static_cast<std::size_t>(levelSetCount_);
    if (minLevel_neg.size() < lsCount) minLevel_neg.resize(lsCount, 0);
    if (minLevel_pos.size() < lsCount) minLevel_pos.resize(lsCount, 0);
    if (maxLevel_local.size() < lsCount) maxLevel_local.resize(lsCount, 0);

    maxLevel_ = 0;
    for (int v : maxLevel_local)
        maxLevel_ = std::max(maxLevel_, v);
    nCellsPerLevel.resize(maxLevel_ + 1, 0);
    nCellsPerLevel[0] = cellIndexNum_;
    
    vertPtr_ = vertIndexNum_;
    hEdgePtr_ = hEdgeIndexNum_;
    vEdgePtr_ = vEdgeIndexNum_;
    cellPtr_ = cellIndexNum_;

    resizeMemory();

    // Optional static refinement boxes (demo support): create coarse--fine
    // interfaces at construction time.
    if (!dictMesh.refineBoxes.empty())
        applyRefineBoxes(dictMesh.refineBoxes);
}

QuadtreeMesh::~QuadtreeMesh() 
{
    delete[] vertLabelSpace_;
    delete[] hEdgeLabelSpace_;
    delete[] vEdgeLabelSpace_;
    delete[] cellLabelSpace_;
}

QuadtreeMesh::QuadtreeMesh(const QuadtreeMesh& source, CopyMetadataTag tag)
    : BackgroundMesh(source, tag), maxLevel_(source.maxLevel_),
      maxVertNum_(std::max({source.maxVertNum_, source.maxHEdgeNum_, source.maxVEdgeNum_, source.maxCellNum_})),
      maxHEdgeNum_(maxVertNum_), maxVEdgeNum_(maxVertNum_), maxCellNum_(maxVertNum_),
      maxInitialIterations_(source.maxInitialIterations_), minLevel_neg(source.minLevel_neg),
      minLevel_pos(source.minLevel_pos), maxLevel_local(source.maxLevel_local),
      C_Lipschitz(source.C_Lipschitz), nCellsPerLevel(source.nCellsPerLevel)
{
    // Preserve the original clone capacity policy, not a new refinement cap.
    vertLabelCapacity_ = maxVertNum_; hEdgeLabelCapacity_ = maxHEdgeNum_;
    vEdgeLabelCapacity_ = maxVEdgeNum_; cellLabelCapacity_ = maxCellNum_;
}

std::unique_ptr<QuadtreeMesh> QuadtreeMesh::clone() const
{
    auto result = std::unique_ptr<QuadtreeMesh>(new QuadtreeMesh(*this, CopyMetadataTag{}));

    // Only live records have a physical/topological identity. Default unused
    // slots in the candidate are deliberately not copies of retired entities.
    // Rebuild the free-label PREFIX from occupancy; the unused suffix of the
    // source free-list may hold legacy/stale entries. Parent records are retained.
    const auto copyEntities = [this](const auto& source, int capacity, auto& destination,
                                   int destinationCapacity, int*& freeLabels)
    {
        using Pool=std::decay_t<decltype(destination)>;
        Pool records;
        if(levelSetCount_==0)records.resetLazy(destinationCapacity,0);
        else records.resetDense(destinationCapacity);
        auto labels = std::make_unique<int[]>(destinationCapacity);
        int freeCount = 0;
        for (int label = 0; label < destinationCapacity; ++label)
        {
            if (label < capacity && source.peek(label).label() != -1)
            {
                if (source.peek(label).label() != label)
                    throw std::runtime_error("cannot clone inconsistent quadtree storage labels");
                records[label] = source.peek(label);
            }
            else
            {
                if(levelSetCount_!=0){
                    records[label].setLabel(-1);
                    records[label].ensureLevelSetSize(static_cast<std::size_t>(levelSetCount_));
                }
                labels[freeCount++] = label;
            }
        }
        // Unique target chunks and complete free prefix publish together.
        destination.swap(records);
        freeLabels = labels.release();
    };
    copyEntities(verts_, maxVertNum_, result->verts_, result->maxVertNum_, result->vertLabelSpace_);
    copyEntities(hEdges_, maxHEdgeNum_, result->hEdges_, result->maxHEdgeNum_, result->hEdgeLabelSpace_);
    copyEntities(vEdges_, maxVEdgeNum_, result->vEdges_, result->maxVEdgeNum_, result->vEdgeLabelSpace_);
    copyEntities(cells_, maxCellNum_, result->cells_, result->maxCellNum_, result->cellLabelSpace_);
    const auto copyIndices = [](const int* source, int count, int*& destination)
    {
        auto indices = std::make_unique<int[]>(count);
        std::copy_n(source, count, indices.get());
        delete[] destination;
        destination = indices.release();
    };
    copyIndices(vertIndexToLabel, vertIndexNum_, result->vertIndexToLabel);
    copyIndices(hEdgeIndexToLabel, hEdgeIndexNum_, result->hEdgeIndexToLabel);
    copyIndices(vEdgeIndexToLabel, vEdgeIndexNum_, result->vEdgeIndexToLabel);
    copyIndices(cellIndexToLabel, cellIndexNum_, result->cellIndexToLabel);
    result->vertIndexNum_ = vertIndexNum_;
    result->hEdgeIndexNum_ = hEdgeIndexNum_;
    result->vEdgeIndexNum_ = vEdgeIndexNum_;
    result->cellIndexNum_ = cellIndexNum_;
    result->vertPtr_ = vertPtr_;
    result->hEdgePtr_ = hEdgePtr_;
    result->vEdgePtr_ = vEdgePtr_;
    result->cellPtr_ = cellPtr_;
    result->maxLevel_ = maxLevel_;
    result->nCellsPerLevel = nCellsPerLevel;
    result->negFaceVerts = negFaceVerts;
    result->posFaceVerts = posFaceVerts;
    return result;
}

int QuadtreeMesh::minimumAdaptationLevel() const
{
    if (minLevel_neg.empty() || minLevel_pos.empty()) return 0;
    int negative = minLevel_neg.front(), positive = minLevel_pos.front();
    for (std::size_t i=1; i<minLevel_neg.size(); ++i) negative = std::min(negative, minLevel_neg[i]);
    for (std::size_t i=1; i<minLevel_pos.size(); ++i) positive = std::min(positive, minLevel_pos[i]);
    return std::min(negative, positive);
}

int QuadtreeMesh::boundedTargetLevel(double value, int minimum) const
{
    int aim = minimum;
    if (std::isfinite(value)) aim = static_cast<int>(std::lround(value));
    return std::max(minimum, std::min(maxLevel_, aim));
}

bool QuadtreeMesh::targetRefinementAllowed(const Cell& cell) const
{
    for (int i=0; i<4; ++i)
        if (cell.neiCell[i][0] != -1 && cells_[cell.neiCell[i][0]].level() < cell.level()) return false;
    return true;
}

bool QuadtreeMesh::targetCoarseningAllowed(const Cell& cell) const
{
    for (int i=0; i<4; ++i)
        if (cell.neiCell[i][0] != -1 && cells_[cell.neiCell[i][0]].level() > cell.level()) return false;
    return true;
}

void QuadtreeMesh::swapStorage(QuadtreeMesh& candidate)
{
    if (minCoord.x != candidate.minCoord.x || minCoord.y != candidate.minCoord.y
        || maxCoord.x != candidate.maxCoord.x || maxCoord.y != candidate.maxCoord.y
        || Nx != candidate.Nx || Ny != candidate.Ny || levelSetCount_ != candidate.levelSetCount_)
        throw std::invalid_argument("quadtree storage exchange requires the same background domain, lattice and level-set count");
    using std::swap;
    verts_.swap(candidate.verts_); hEdges_.swap(candidate.hEdges_);
    vEdges_.swap(candidate.vEdges_); cells_.swap(candidate.cells_);
    swap(vertIndexToLabel, candidate.vertIndexToLabel); swap(hEdgeIndexToLabel, candidate.hEdgeIndexToLabel);
    swap(vEdgeIndexToLabel, candidate.vEdgeIndexToLabel); swap(cellIndexToLabel, candidate.cellIndexToLabel);
    swap(vertIndexNum_, candidate.vertIndexNum_); swap(hEdgeIndexNum_, candidate.hEdgeIndexNum_);
    swap(vEdgeIndexNum_, candidate.vEdgeIndexNum_); swap(cellIndexNum_, candidate.cellIndexNum_);
    swap(vertLabelCapacity_, candidate.vertLabelCapacity_); swap(hEdgeLabelCapacity_, candidate.hEdgeLabelCapacity_);
    swap(vEdgeLabelCapacity_, candidate.vEdgeLabelCapacity_); swap(cellLabelCapacity_, candidate.cellLabelCapacity_);
    swap(maxVertNum_, candidate.maxVertNum_); swap(maxHEdgeNum_, candidate.maxHEdgeNum_);
    swap(maxVEdgeNum_, candidate.maxVEdgeNum_); swap(maxCellNum_, candidate.maxCellNum_);
    swap(vertPtr_, candidate.vertPtr_); swap(hEdgePtr_, candidate.hEdgePtr_);
    swap(vEdgePtr_, candidate.vEdgePtr_); swap(cellPtr_, candidate.cellPtr_);
    swap(vertLabelSpace_, candidate.vertLabelSpace_); swap(hEdgeLabelSpace_, candidate.hEdgeLabelSpace_);
    swap(vEdgeLabelSpace_, candidate.vEdgeLabelSpace_); swap(cellLabelSpace_, candidate.cellLabelSpace_);
    swap(maxLevel_, candidate.maxLevel_); swap(maxInitialIterations_, candidate.maxInitialIterations_);
    swap(C_Lipschitz, candidate.C_Lipschitz);
    minLevel_neg.swap(candidate.minLevel_neg); minLevel_pos.swap(candidate.minLevel_pos);
    maxLevel_local.swap(candidate.maxLevel_local); nCellsPerLevel.swap(candidate.nCellsPerLevel);
    swap(nativeAdaptationReport_,candidate.nativeAdaptationReport_);
    negFaceVerts.swap(candidate.negFaceVerts); posFaceVerts.swap(candidate.posFaceVerts);
}

void QuadtreeMesh::resizeMemory()
{
    LAMP::Mesh::Detail::EntityPool<Vert> tempVerts;
    LAMP::Mesh::Detail::EntityPool<Edge> tempHEdges;
    LAMP::Mesh::Detail::EntityPool<Edge> tempVEdges;
    LAMP::Mesh::Detail::EntityPool<Cell> tempCells;

    if (levelSetCount_ == 0)
    {
        tempVerts.resetLazy(maxVertNum_, 0);
        tempHEdges.resetLazy(maxHEdgeNum_, 0);
        tempVEdges.resetLazy(maxVEdgeNum_, 0);
        tempCells.resetLazy(maxCellNum_, 0);
    }
    else
    {
        tempVerts.resetDense(maxVertNum_);
        tempHEdges.resetDense(maxHEdgeNum_);
        tempVEdges.resetDense(maxVEdgeNum_);
        tempCells.resetDense(maxCellNum_);
    }

    for (int i = 0; i < vertIndexNum_; i++)
    {
        tempVerts[i] = verts_.peek(i);
    }

    for (int i = 0; i < hEdgeIndexNum_; i++)
    {
        tempHEdges[i] = hEdges_.peek(i);
    }

    for (int i = 0; i < vEdgeIndexNum_; i++)
    {
        tempVEdges[i] = vEdges_.peek(i);
    }

    for (int i = 0; i < cellIndexNum_; i++)
    {
        tempCells[i] = cells_.peek(i);
    }

    vertLabelSpace_ = new int[maxVertNum_]{};
    hEdgeLabelSpace_ = new int[maxHEdgeNum_]{};
    vEdgeLabelSpace_ = new int[maxVEdgeNum_]{};
    cellLabelSpace_ = new int[maxCellNum_]{};

    int _remainVertNum = maxVertNum_ - vertPtr_;
    int _remainHEdgeNum = maxHEdgeNum_ - hEdgePtr_;
    int _remainVEdgeNum = maxVEdgeNum_ - vEdgePtr_;
    int _remainCellNum = maxCellNum_ - cellPtr_;
    
    for (int i = 0; i < _remainVertNum; i++)
    {
        vertLabelSpace_[i] = vertPtr_ + i;
    }

    for (int i = 0; i < _remainHEdgeNum; i++)
    {
        hEdgeLabelSpace_[i] = hEdgePtr_ + i;
    }

    for (int i = 0; i < _remainVEdgeNum; i++)
    {
        vEdgeLabelSpace_[i] = vEdgePtr_ + i;
    }

    for (int i = 0; i < _remainCellNum; i++)
    {
        cellLabelSpace_[i] = cellPtr_ + i;
    }


    verts_.swap(tempVerts);
    hEdges_.swap(tempHEdges);
    vEdges_.swap(tempVEdges);
    cells_.swap(tempCells);

    // BackgroundMesh hint-based queries index cells_ by label, so keep an accurate
    // allocated label-space capacity after buffer growth.
    vertLabelCapacity_ = maxVertNum_;
    hEdgeLabelCapacity_ = maxHEdgeNum_;
    vEdgeLabelCapacity_ = maxVEdgeNum_;
    cellLabelCapacity_ = maxCellNum_;

    tempVerts.clear();
    tempHEdges.clear();
    tempVEdges.clear();
    tempCells.clear();

    // Accessing every LS-free label here would materialize the unused chunks.
    if (levelSetCount_ != 0)
    {
        const std::size_t n = static_cast<std::size_t>(std::max(1, levelSetCount_));
        for (int i = 0; i < maxVertNum_; ++i)
            verts_[i].ensureLevelSetSize(n);
        for (int i = 0; i < maxHEdgeNum_; ++i)
            hEdges_[i].ensureLevelSetSize(n);
        for (int i = 0; i < maxVEdgeNum_; ++i)
            vEdges_[i].ensureLevelSetSize(n);
        for (int i = 0; i < maxCellNum_; ++i)
            cells_[i].ensureLevelSetSize(n);
    }
}

// -----------------------------------------------------------------------------
// Static refine boxes (demo support)
// -----------------------------------------------------------------------------
void QuadtreeMesh::prepareTransferIndices()
{
    [[maybe_unused]] const bool useOmp = (cellIndexNum_ > 4096 || maxCellNum_ > 4096);
#ifdef _OPENMP
#pragma omp parallel for if(useOmp)
#endif
    for (int i = 0; i < maxCellNum_; ++i)
    {
        if(auto* record=cells_.materialized(i)){
            record->setPreIndex(-1);
            record->relatedPreIndex.clear();
            record->relatedPreCoord.clear();
        }
    }
#ifdef _OPENMP
#pragma omp parallel for if(useOmp)
#endif
    for (int i = 0; i < maxHEdgeNum_; ++i)
    {
        if(auto* record=hEdges_.materialized(i)){
            record->setPreIndex(-1);
            record->relatedPreIndex.clear();
            record->relatedPreCoord.clear();
        }
    }
#ifdef _OPENMP
#pragma omp parallel for if(useOmp)
#endif
    for (int i = 0; i < maxVEdgeNum_; ++i)
    {
        if(auto* record=vEdges_.materialized(i)){
            record->setPreIndex(-1);
            record->relatedPreIndex.clear();
            record->relatedPreCoord.clear();
        }
    }
#ifdef _OPENMP
#pragma omp parallel for if(useOmp)
#endif
    for (int i = 0; i < cellIndexNum_; ++i) cells_[cellIndexToLabel[i]].setPreIndex(i);
#ifdef _OPENMP
#pragma omp parallel for if(useOmp)
#endif
    for (int i = 0; i < hEdgeIndexNum_; ++i) hEdges_[hEdgeIndexToLabel[i]].setPreIndex(i);
#ifdef _OPENMP
#pragma omp parallel for if(useOmp)
#endif
    for (int i = 0; i < vEdgeIndexNum_; ++i) vEdges_[vEdgeIndexToLabel[i]].setPreIndex(i);
}

void QuadtreeMesh::applyRefineBoxes(const std::vector<DictionaryMesh::RefineBox>& boxes)
{
    if (boxes.empty())
        return;

    int maxTarget = 0;
    for (const auto& b : boxes)
        maxTarget = std::max(maxTarget, b.level);
    if (maxTarget <= 0)
        return;

    // Ensure maxLevel_ / per-level stats can represent the target level.
    if (maxTarget > maxLevel_)
    {
        maxLevel_ = maxTarget;
        nCellsPerLevel.assign(maxLevel_ + 1, 0);
        nCellsPerLevel[0] = cellIndexNum_;
    }

    auto inBox = [](const Vector2D& p, const DictionaryMesh::RefineBox& b) -> bool {
        const double x0 = std::min(b.min.x, b.max.x);
        const double x1 = std::max(b.min.x, b.max.x);
        const double y0 = std::min(b.min.y, b.max.y);
        const double y1 = std::max(b.min.y, b.max.y);
        return (p.x >= x0 && p.x <= x1 && p.y >= y0 && p.y <= y1);
    };

    // Refine progressively: level -> level+1 until reaching maxTarget.
    for (int lvl = 0; lvl < maxTarget; ++lvl)
    {
        bool anyMarked = false;
        for (int i = 0; i < cellIndexNum_; ++i)
        {
            const int cellLabel = cellIndexToLabel[i];
            Cell& c = cells_[cellLabel];
            if (c.level() != lvl)
                continue;
            // Leaf check: only refine leaves.
            const bool isLeaf = (c.son[0] == -1 && c.son[1] == -1 && c.son[2] == -1 && c.son[3] == -1);
            if (!isLeaf)
                continue;

            bool wantsRefine = false;
            const Vector2D cc = c.centerCoord();
            for (const auto& b : boxes)
            {
                if (b.level > lvl && inBox(cc, b))
                {
                    wantsRefine = true;
                    break;
                }
            }
            c.needRefine = wantsRefine;
            if (wantsRefine)
                anyMarked = true;
        }

        if (!anyMarked)
            continue;

        // History gathering distinguishes existing entities from children by
        // preIndex. An unprepared base edge had both preIndex=-1 and parent=-1,
        // so it was incorrectly followed through a nonexistent parent.
        prepareTransferIndices();
        subdivideCell(lvl);
    }
    // Match adapt(): subdivision updates connectivity/indices, but does not
    // classify the new boundary and hanging vertices. Consumers of supCell
    // (notably native Hydro virtual stencils) require these derived identities
    // before the static-refinement mesh is published.
    defineTypes();
}

void QuadtreeMesh::subdivideCell(const int currentLevel)
{subdivideCellWithWorkspace(currentLevel,nullptr);}
void QuadtreeMesh::subdivideCellWithWorkspace(int currentLevel,NativeAdaptWorkspace* workspace)
{
    const LAMP::Mesh::DyadicGridCoordinates coordinates(minCoord, maxCoord, Nx, Ny);
    // Fresh entities have empty LS vectors; enforce the mesh-wide storage
    // contract immediately, including the static refineBoxes path.
    auto assignNewEntity = [this,workspace](auto& destination, const auto& source)
    {
        recordNativeWrite(workspace,&destination);
        destination = source;
        destination.ensureLevelSetSize(static_cast<std::size_t>(levelSetCount_));
    };
    const int newLevel = currentLevel + 1;
    int nNewVert = 0;
    int nNewHEdge = 0;
    int nNewVEdge = 0;
    int nNewCell = 0;
    int nSkipRefine = 0;

    // Conservative per-cell allocation bounds to avoid exceeding max buffers.
    const int kMaxNewCellsPerRefine = 4;
    const int kMaxNewVertsPerRefine = 5;
    const int kMaxNewHEdgesPerRefine = 6;
    const int kMaxNewVEdgesPerRefine = 6;

    nCellsPerLevel[newLevel] = 0;

    for (int c = 0; c < cellIndexNum_; c++)
    {
        const int cellLabel = cellIndexToLabel[c];
        Cell& cell = cells_[cellIndexToLabel[c]];
        if (cell.needRefine)
        {
            const int availCells = maxCellNum_ - (cellPtr_ + nNewCell);
            const int availVerts = maxVertNum_ - (vertPtr_ + nNewVert);
            const int availHEdges = maxHEdgeNum_ - (hEdgePtr_ + nNewHEdge);
            const int availVEdges = maxVEdgeNum_ - (vEdgePtr_ + nNewVEdge);
            if (availCells < kMaxNewCellsPerRefine ||
                availVerts < kMaxNewVertsPerRefine ||
                availHEdges < kMaxNewHEdgesPerRefine ||
                availVEdges < kMaxNewVEdgesPerRefine)
            {
                // Not enough buffer space; skip refining this cell to avoid OOB.
                cell.needRefine = false;
                nSkipRefine++;
                continue;
            }
            cell.needRefine = false;
            recordNativeWrite(workspace,&cell);recordNativeCellEdges(workspace,cell);

            const auto address = LAMP::Mesh::cellAddress(*this, cell);
            if (address.level != static_cast<unsigned>(currentLevel)
                || address.level > std::numeric_limits<unsigned>::max() - 2
                || address.i > (std::numeric_limits<std::uint64_t>::max() - 4) / 4
                || address.j > (std::numeric_limits<std::uint64_t>::max() - 4) / 4)
                throw std::runtime_error("[QuadtreeMesh] invalid or overflowing subdivision address");
            const Vector2D newCellSize = coordinates.cellSize(address.level + 1);
            // Every new location belongs to the parent's integer quarter-grid.
            // Shared points therefore have one coordinate regardless of which
            // neighbor created them or how its parent centers were rounded.
            const std::uint64_t quarterI = 4 * address.i, quarterJ = 4 * address.j;
            const auto quarterPoint = [&](unsigned dx, unsigned dy) {
                return coordinates.point(quarterI + dx, quarterJ + dy, address.level + 2);
            };

            bool newDir[4] = {false, false, false, false};

            // 生成新的单元
            Vector2D cellCenterCoord[4];
            cellCenterCoord[0] = quarterPoint(1, 1);
            cellCenterCoord[1] = quarterPoint(3, 1);
            cellCenterCoord[2] = quarterPoint(3, 3);
            cellCenterCoord[3] = quarterPoint(1, 3);
            Cell newCell[4];
            for (int i = 0; i < 4; i++)
            {
                const int sonCellLabel = cellLabelSpace_[nNewCell];
                newCell[i] = Cell(sonCellLabel, cellCenterCoord[i], newLevel, newCellSize);
                assignNewEntity(cells_[sonCellLabel], newCell[i]);
                cells_[sonCellLabel].setLevel(cell.level() + 1);
                cells_[sonCellLabel].parent = cellLabel;
                cell.son[i] = sonCellLabel;
                cells_[sonCellLabel].cellToVert[i] = cell.cellToVert[i];

                nCellsPerLevel[newLevel]++;
                nNewCell++;
            }
            verts_[cell.cellToVert[0]].vertToCell[2] = cell.son[0];
            verts_[cell.cellToVert[1]].vertToCell[3] = cell.son[1];
            verts_[cell.cellToVert[2]].vertToCell[0] = cell.son[2];
            verts_[cell.cellToVert[3]].vertToCell[1] = cell.son[3];

            findSubdivideRelatedCells(cell,workspace);

            // New neighbor relations
            if (cell.neiCell[0][1] != -1)
            {
                cells_[cell.son[0]].neiCell[0][0] = cell.neiCell[0][0];
                cells_[cell.son[3]].neiCell[0][0] = cell.neiCell[0][1];

                cells_[cell.neiCell[0][0]].neiCell[2][0] = cell.son[0];
                cells_[cell.neiCell[0][1]].neiCell[2][0] = cell.son[3];
            }
            else if (cell.neiCell[0][0] != -1)
            {
                cells_[cell.son[0]].neiCell[0][0] = cell.neiCell[0][0];
                cells_[cell.son[3]].neiCell[0][0] = cell.neiCell[0][0];

                cells_[cell.neiCell[0][0]].neiCell[2][0] = cell.son[0];
                cells_[cell.neiCell[0][0]].neiCell[2][1] = cell.son[3];
            }

            if (cell.neiCell[1][1] != -1)
            {
                cells_[cell.son[0]].neiCell[1][0] = cell.neiCell[1][0];
                cells_[cell.son[1]].neiCell[1][0] = cell.neiCell[1][1];

                cells_[cell.neiCell[1][0]].neiCell[3][0] = cell.son[0];
                cells_[cell.neiCell[1][1]].neiCell[3][0] = cell.son[1];
            }
            else if (cell.neiCell[1][0] != -1)
            {
                cells_[cell.son[0]].neiCell[1][0] = cell.neiCell[1][0];
                cells_[cell.son[1]].neiCell[1][0] = cell.neiCell[1][0];

                cells_[cell.neiCell[1][0]].neiCell[3][0] = cell.son[0];
                cells_[cell.neiCell[1][0]].neiCell[3][1] = cell.son[1];
            }

            if (cell.neiCell[2][1] != -1)
            {
                cells_[cell.son[1]].neiCell[2][0] = cell.neiCell[2][0];
                cells_[cell.son[2]].neiCell[2][0] = cell.neiCell[2][1];

                cells_[cell.neiCell[2][0]].neiCell[0][0] = cell.son[1];
                cells_[cell.neiCell[2][1]].neiCell[0][0] = cell.son[2];
            }
            else if (cell.neiCell[2][0] != -1)
            {
                cells_[cell.son[1]].neiCell[2][0] = cell.neiCell[2][0];
                cells_[cell.son[2]].neiCell[2][0] = cell.neiCell[2][0];

                cells_[cell.neiCell[2][0]].neiCell[0][0] = cell.son[1];
                cells_[cell.neiCell[2][0]].neiCell[0][1] = cell.son[2];
            }

            if (cell.neiCell[3][1] != -1)
            {
                cells_[cell.son[3]].neiCell[3][0] = cell.neiCell[3][0];
                cells_[cell.son[2]].neiCell[3][0] = cell.neiCell[3][1];

                cells_[cell.neiCell[3][0]].neiCell[1][0] = cell.son[3];
                cells_[cell.neiCell[3][1]].neiCell[1][0] = cell.son[2];
            }
            else if (cell.neiCell[3][0] != -1)
            {
                cells_[cell.son[3]].neiCell[3][0] = cell.neiCell[3][0];
                cells_[cell.son[2]].neiCell[3][0] = cell.neiCell[3][0];

                cells_[cell.neiCell[3][0]].neiCell[1][0] = cell.son[3];
                cells_[cell.neiCell[3][0]].neiCell[1][1] = cell.son[2];
            }
            // Internal
            // child cell 0
            {
                cells_[cell.son[0]].neiCell[2][0] = cell.son[1];
                cells_[cell.son[0]].neiCell[3][0] = cell.son[3];
            }
            // child cell 1
            {
                cells_[cell.son[1]].neiCell[0][0] = cell.son[0];
                cells_[cell.son[1]].neiCell[3][0] = cell.son[2];
            }
            // child cell 2
            {
                cells_[cell.son[2]].neiCell[0][0] = cell.son[3];
                cells_[cell.son[2]].neiCell[1][0] = cell.son[1];
            }
            // child cell 3
            {
                cells_[cell.son[3]].neiCell[1][0] = cell.son[0];
                cells_[cell.son[3]].neiCell[2][0] = cell.son[2];
            }

            // 生成新的边
            if (cell.cellToVEdge[0][1] == -1)
            {
                newDir[0] = true;

                Vector2D edgeMidPointCoord[2];               
                edgeMidPointCoord[0] = quarterPoint(0, 1);
                edgeMidPointCoord[1] = quarterPoint(0, 3);

                const int vEdgeLabel0 = vEdgeLabelSpace_[nNewVEdge];
                Edge vEdge0(vEdgeLabel0, edgeMidPointCoord[0], newCellSize.y, false);
                assignNewEntity(vEdges_[vEdgeLabel0], vEdge0);
                vEdges_[vEdgeLabel0].masterCell = cell.son[0];
                vEdges_[vEdgeLabel0].slaveCell = cells_[cell.son[0]].neiCell[0][0];
                vEdges_[vEdgeLabel0].parent = cell.cellToVEdge[0][0];
                vEdges_[vEdges_[vEdgeLabel0].parent].son[0] = vEdgeLabel0;
                verts_[cell.cellToVert[0]].vertToVEdge[1] = vEdgeLabel0;
                cells_[cell.son[0]].cellToVEdge[0][0] = vEdgeLabel0;
                if (cells_[cell.son[0]].neiCell[0][0] != -1)
                {
                    cells_[cells_[cell.son[0]].neiCell[0][0]].cellToVEdge[1][0] = vEdgeLabel0;
                }

                nNewVEdge++;

                const int vEdgeLabel1 = vEdgeLabelSpace_[nNewVEdge];
                Edge vEdge1(vEdgeLabel1, edgeMidPointCoord[1], newCellSize.y, false);
                assignNewEntity(vEdges_[vEdgeLabel1], vEdge1);
                vEdges_[vEdgeLabel1].masterCell = cell.son[3];
                vEdges_[vEdgeLabel1].slaveCell = cells_[cell.son[3]].neiCell[0][0];
                vEdges_[vEdgeLabel1].parent = cell.cellToVEdge[0][0];
                vEdges_[vEdges_[vEdgeLabel1].parent].son[1] = vEdgeLabel1;
                verts_[cell.cellToVert[3]].vertToVEdge[0] = vEdgeLabel1;
                cells_[cell.son[3]].cellToVEdge[0][0] = vEdgeLabel1;
                if (cells_[cell.son[3]].neiCell[0][0] != -1)
                {
                    cells_[cells_[cell.son[3]].neiCell[0][0]].cellToVEdge[1][1] = vEdgeLabel1;
                }

                nNewVEdge++;
            }
            else
            {
                cells_[cell.son[0]].cellToVEdge[0][0] = cell.cellToVEdge[0][0];
                cells_[cell.son[3]].cellToVEdge[0][0] = cell.cellToVEdge[0][1];
                vEdges_[cell.cellToVEdge[0][0]].masterCell = cell.son[0];
                vEdges_[cell.cellToVEdge[0][1]].masterCell = cell.son[3];
            }

            if (cell.cellToHEdge[0][1] == -1)
            {
                newDir[1] = true;

                Vector2D edgeMidPointCoord[2];
                edgeMidPointCoord[0] = quarterPoint(1, 0);
                edgeMidPointCoord[1] = quarterPoint(3, 0);

                const int hEdgeLabel0 = hEdgeLabelSpace_[nNewHEdge];
                Edge hEdge0(hEdgeLabel0, edgeMidPointCoord[0], newCellSize.x, true);
                assignNewEntity(hEdges_[hEdgeLabel0], hEdge0);
                hEdges_[hEdgeLabel0].masterCell = cell.son[0];
                hEdges_[hEdgeLabel0].slaveCell = cells_[cell.son[0]].neiCell[1][0];
                hEdges_[hEdgeLabel0].parent = cell.cellToHEdge[0][0];
                hEdges_[hEdges_[hEdgeLabel0].parent].son[0] = hEdgeLabel0;
                verts_[cell.cellToVert[0]].vertToHEdge[1] = hEdgeLabel0;
                cells_[cell.son[0]].cellToHEdge[0][0] = hEdgeLabel0;
                if (cells_[cell.son[0]].neiCell[1][0] != -1)
                {
                    cells_[cells_[cell.son[0]].neiCell[1][0]].cellToHEdge[1][0] = hEdgeLabel0;
                }

                nNewHEdge++;

                const int hEdgeLabel1 = hEdgeLabelSpace_[nNewHEdge];
                Edge hEdge1(hEdgeLabel1, edgeMidPointCoord[1], newCellSize.x, true);
                assignNewEntity(hEdges_[hEdgeLabel1], hEdge1);
                hEdges_[hEdgeLabel1].masterCell = cell.son[1];
                hEdges_[hEdgeLabel1].slaveCell = cells_[cell.son[1]].neiCell[1][0];
                hEdges_[hEdgeLabel1].parent = cell.cellToHEdge[0][0];
                hEdges_[hEdges_[hEdgeLabel1].parent].son[1] = hEdgeLabel1;
                verts_[cell.cellToVert[1]].vertToHEdge[0] = hEdgeLabel1;
                cells_[cell.son[1]].cellToHEdge[0][0] = hEdgeLabel1;
                if (cells_[cell.son[1]].neiCell[1][0] != -1)
                {
                    cells_[cells_[cell.son[1]].neiCell[1][0]].cellToHEdge[1][1] = hEdgeLabel1;
                }

                nNewHEdge++;
            }
            else
            {
                cells_[cell.son[0]].cellToHEdge[0][0] = cell.cellToHEdge[0][0];
                cells_[cell.son[1]].cellToHEdge[0][0] = cell.cellToHEdge[0][1];
                hEdges_[cell.cellToHEdge[0][0]].masterCell = cell.son[0];
                hEdges_[cell.cellToHEdge[0][1]].masterCell = cell.son[1];
            }

            if (cell.cellToVEdge[1][1] == -1)// 右
            {
                newDir[2] = true;

                Vector2D edgeMidPointCoord[2];
                edgeMidPointCoord[0] = quarterPoint(4, 1);
                edgeMidPointCoord[1] = quarterPoint(4, 3);

                const int vEdgeLabel0 = vEdgeLabelSpace_[nNewVEdge];
                Edge vEdge0(vEdgeLabel0, edgeMidPointCoord[0], newCellSize.y, false);
                assignNewEntity(vEdges_[vEdgeLabel0], vEdge0);
                vEdges_[vEdgeLabel0].slaveCell = cell.son[1];
                vEdges_[vEdgeLabel0].masterCell = cells_[cell.son[1]].neiCell[2][0];
                vEdges_[vEdgeLabel0].parent = cell.cellToVEdge[1][0];
                vEdges_[vEdges_[vEdgeLabel0].parent].son[0] = vEdgeLabel0;
                verts_[cell.cellToVert[1]].vertToVEdge[1] = vEdgeLabel0;
                cells_[cell.son[1]].cellToVEdge[1][0] = vEdgeLabel0;
                if (cells_[cell.son[1]].neiCell[2][0] != -1)
                {
                    cells_[cells_[cell.son[1]].neiCell[2][0]].cellToVEdge[0][0] = vEdgeLabel0;
                }

                nNewVEdge++;

                const int vEdgeLabel1 = vEdgeLabelSpace_[nNewVEdge];
                Edge vEdge1(vEdgeLabel1, edgeMidPointCoord[1], newCellSize.y, false);
                assignNewEntity(vEdges_[vEdgeLabel1], vEdge1);
                vEdges_[vEdgeLabel1].slaveCell = cell.son[2];
                vEdges_[vEdgeLabel1].masterCell = cells_[cell.son[2]].neiCell[2][0];
                vEdges_[vEdgeLabel1].parent = cell.cellToVEdge[1][0];
                vEdges_[vEdges_[vEdgeLabel1].parent].son[1] = vEdgeLabel1;
                verts_[cell.cellToVert[2]].vertToVEdge[0] = vEdgeLabel1;
                cells_[cell.son[2]].cellToVEdge[1][0] = vEdgeLabel1;
                if (cells_[cell.son[2]].neiCell[2][0] != -1)
                {
                    cells_[cells_[cell.son[2]].neiCell[2][0]].cellToVEdge[0][1] = vEdgeLabel1;
                }

                nNewVEdge++;
            }
            else
            {
                cells_[cell.son[1]].cellToVEdge[1][0] = cell.cellToVEdge[1][0];
                cells_[cell.son[2]].cellToVEdge[1][0] = cell.cellToVEdge[1][1];
                vEdges_[cell.cellToVEdge[1][0]].slaveCell = cell.son[1];
                vEdges_[cell.cellToVEdge[1][1]].slaveCell = cell.son[2];
            }

            if (cell.cellToHEdge[1][1] == -1)// 上
            {
                newDir[3] = true;
                
                Vector2D edgeMidPointCoord[2];
                edgeMidPointCoord[0] = quarterPoint(1, 4);
                edgeMidPointCoord[1] = quarterPoint(3, 4);

                const int hEdgeLabel0 = hEdgeLabelSpace_[nNewHEdge];
                Edge hEdge0(hEdgeLabel0, edgeMidPointCoord[0], newCellSize.x, true);
                assignNewEntity(hEdges_[hEdgeLabel0], hEdge0);
                hEdges_[hEdgeLabel0].slaveCell = cell.son[3];
                hEdges_[hEdgeLabel0].masterCell = cells_[cell.son[3]].neiCell[3][0];
                hEdges_[hEdgeLabel0].parent = cell.cellToHEdge[1][0];
                hEdges_[hEdges_[hEdgeLabel0].parent].son[0] = hEdgeLabel0;
                verts_[cell.cellToVert[3]].vertToHEdge[1] = hEdgeLabel0;
                cells_[cell.son[3]].cellToHEdge[1][0] = hEdgeLabel0;
                if (cells_[cell.son[3]].neiCell[3][0] != -1)
                {
                    cells_[cells_[cell.son[3]].neiCell[3][0]].cellToHEdge[0][0] = hEdgeLabel0;
                }

                nNewHEdge++;

                const int hEdgeLabel1 = hEdgeLabelSpace_[nNewHEdge];
                Edge hEdge1(hEdgeLabel1, edgeMidPointCoord[1], newCellSize.x, true);
                assignNewEntity(hEdges_[hEdgeLabel1], hEdge1);
                hEdges_[hEdgeLabel1].slaveCell = cell.son[2];
                hEdges_[hEdgeLabel1].masterCell = cells_[cell.son[2]].neiCell[3][0];
                hEdges_[hEdgeLabel1].parent = cell.cellToHEdge[1][0];
                hEdges_[hEdges_[hEdgeLabel1].parent].son[1] = hEdgeLabel1;
                verts_[cell.cellToVert[2]].vertToHEdge[0] = hEdgeLabel1;
                cells_[cell.son[2]].cellToHEdge[1][0] = hEdgeLabel1;
                if (cells_[cell.son[2]].neiCell[3][0] != -1)
                {
                    cells_[cells_[cell.son[2]].neiCell[3][0]].cellToHEdge[0][1] = hEdgeLabel1;
                }

                nNewHEdge++;
            }
            else
            {
                cells_[cell.son[3]].cellToHEdge[1][0] = cell.cellToHEdge[1][0];
                cells_[cell.son[2]].cellToHEdge[1][0] = cell.cellToHEdge[1][1];
                hEdges_[cell.cellToHEdge[1][0]].slaveCell = cell.son[3];
                hEdges_[cell.cellToHEdge[1][1]].slaveCell = cell.son[2];
            }

            /*内部新生成4条边*/          
            // 中左
            {
                Vector2D edgeMidPointCoord0 = quarterPoint(1, 2);
                
                const int hEdgeLabel0 = hEdgeLabelSpace_[nNewHEdge];
                Edge hEdge0(hEdgeLabel0, edgeMidPointCoord0, newCellSize.x, true);
                assignNewEntity(hEdges_[hEdgeLabel0], hEdge0);
                hEdges_[hEdgeLabel0].masterCell = cell.son[3];
                hEdges_[hEdgeLabel0].slaveCell = cell.son[0];
                cells_[cell.son[0]].cellToHEdge[1][0] = hEdgeLabel0;
                cells_[cell.son[3]].cellToHEdge[0][0] = hEdgeLabel0;

                nNewHEdge++;
            }
            // 中下
            {
                Vector2D edgeMidPointCoord1 = quarterPoint(2, 1);
                
                const int vEdgeLabel0 = vEdgeLabelSpace_[nNewVEdge];
                Edge vEdge0(vEdgeLabel0, edgeMidPointCoord1, newCellSize.y, false);
                assignNewEntity(vEdges_[vEdgeLabel0], vEdge0);
                vEdges_[vEdgeLabel0].masterCell = cell.son[1];
                vEdges_[vEdgeLabel0].slaveCell = cell.son[0];
                cells_[cell.son[0]].cellToVEdge[1][0] = vEdgeLabel0;
                cells_[cell.son[1]].cellToVEdge[0][0] = vEdgeLabel0;

                nNewVEdge++;
            }
            // 中右
            {
                Vector2D edgeMidPointCoord2 = quarterPoint(3, 2);
                
                const int hEdgeLabel1 = hEdgeLabelSpace_[nNewHEdge];
                Edge hEdge1(hEdgeLabel1, edgeMidPointCoord2, newCellSize.x, true);
                assignNewEntity(hEdges_[hEdgeLabel1], hEdge1);
                hEdges_[hEdgeLabel1].masterCell = cell.son[2];
                hEdges_[hEdgeLabel1].slaveCell = cell.son[1];
                cells_[cell.son[2]].cellToHEdge[0][0] = hEdgeLabel1;
                cells_[cell.son[1]].cellToHEdge[1][0] = hEdgeLabel1;

                nNewHEdge++;
            }
            // 中上
            {
                Vector2D edgeMidPointCoord3 = quarterPoint(2, 3);
                
                const int vEdgeLabel1 = vEdgeLabelSpace_[nNewVEdge];
                Edge vEdge1(vEdgeLabel1, edgeMidPointCoord3, newCellSize.y, false);
                assignNewEntity(vEdges_[vEdgeLabel1], vEdge1);
                vEdges_[vEdgeLabel1].masterCell = cell.son[2];
                vEdges_[vEdgeLabel1].slaveCell = cell.son[3];
                cells_[cell.son[2]].cellToVEdge[0][0] = vEdgeLabel1;
                cells_[cell.son[3]].cellToVEdge[1][0] = vEdgeLabel1;

                nNewVEdge++;
            }

            Edge& hEdge00 = hEdges_[cells_[cell.son[0]].cellToHEdge[0][0]];
            Edge& hEdge01 = hEdges_[cells_[cell.son[1]].cellToHEdge[0][0]];
            Edge& hEdge10 = hEdges_[cells_[cell.son[0]].cellToHEdge[1][0]];
            Edge& hEdge11 = hEdges_[cells_[cell.son[1]].cellToHEdge[1][0]];
            Edge& hEdge20 = hEdges_[cells_[cell.son[3]].cellToHEdge[1][0]];
            Edge& hEdge21 = hEdges_[cells_[cell.son[2]].cellToHEdge[1][0]];

            Edge& vEdge00 = vEdges_[cells_[cell.son[0]].cellToVEdge[0][0]];
            Edge& vEdge01 = vEdges_[cells_[cell.son[3]].cellToVEdge[0][0]];
            Edge& vEdge10 = vEdges_[cells_[cell.son[0]].cellToVEdge[1][0]];
            Edge& vEdge11 = vEdges_[cells_[cell.son[3]].cellToVEdge[1][0]];
            Edge& vEdge20 = vEdges_[cells_[cell.son[1]].cellToVEdge[1][0]];
            Edge& vEdge21 = vEdges_[cells_[cell.son[2]].cellToVEdge[1][0]];

            if (cell.neiCell[0][0] != -1)
            {
                findSubdivideRelatedHEdges(cell, hEdge00, hEdge10, hEdge20, 0, false,workspace);
                if (cell.neiCell[2][0] != -1)
                {
                    findSubdivideRelatedHEdges(cell, hEdge01, hEdge11, hEdge21, 2, true,workspace);
                }
                else
                {
                    hEdge01.relatedPreIndex = hEdge00.relatedPreIndex;
                    hEdge01.relatedPreCoord = hEdge00.relatedPreCoord;
                    hEdge11.relatedPreIndex = hEdge10.relatedPreIndex;
                    hEdge11.relatedPreCoord = hEdge10.relatedPreCoord;
                    hEdge21.relatedPreIndex = hEdge20.relatedPreIndex;
                    hEdge21.relatedPreCoord = hEdge20.relatedPreCoord;
                }
            }
            else
            {
                findSubdivideRelatedHEdges(cell, hEdge01, hEdge11, hEdge21, 2, true,workspace);
                hEdge00.relatedPreIndex = hEdge01.relatedPreIndex;
                hEdge00.relatedPreCoord = hEdge01.relatedPreCoord;
                hEdge10.relatedPreIndex = hEdge11.relatedPreIndex;
                hEdge10.relatedPreCoord = hEdge11.relatedPreCoord;
                hEdge20.relatedPreIndex = hEdge21.relatedPreIndex;
                hEdge20.relatedPreCoord = hEdge21.relatedPreCoord;
            }

            if (cell.neiCell[1][0] != -1)
            {
                findSubdivideRelatedVEdges(cell, vEdge00, vEdge10, vEdge20, 1, false,workspace);
                if (cell.neiCell[3][0] != -1)
                {
                    findSubdivideRelatedVEdges(cell, vEdge01, vEdge11, vEdge21, 3, true,workspace);
                }
                else
                {
                    vEdge01.relatedPreIndex = vEdge00.relatedPreIndex;
                    vEdge01.relatedPreCoord = vEdge00.relatedPreCoord;
                    vEdge11.relatedPreIndex = vEdge10.relatedPreIndex;
                    vEdge11.relatedPreCoord = vEdge10.relatedPreCoord;
                    vEdge21.relatedPreIndex = vEdge01.relatedPreIndex;
                    vEdge21.relatedPreCoord = vEdge01.relatedPreCoord;
                }
            }
            else
            {
                findSubdivideRelatedVEdges(cell, vEdge01, vEdge11, vEdge21, 3, true,workspace);
                vEdge00.relatedPreIndex = vEdge01.relatedPreIndex;
                vEdge00.relatedPreCoord = vEdge01.relatedPreCoord;
                vEdge10.relatedPreIndex = vEdge11.relatedPreIndex;
                vEdge10.relatedPreCoord = vEdge11.relatedPreCoord;
                vEdge20.relatedPreIndex = vEdge21.relatedPreIndex;
                vEdge20.relatedPreCoord = vEdge21.relatedPreCoord;
            }

            // 生成新的顶点
            // 中间
            int centerVertLabel = -1;
            {
                const int vertLabel = vertLabelSpace_[nNewVert];
                Vector2D centerVertCoord = quarterPoint(2, 2);
                Vert vert(vertLabel, centerVertCoord);
                assignNewEntity(verts_[vertLabel], vert);
                cells_[cell.son[0]].cellToVert[2] = vertLabel;
                cells_[cell.son[1]].cellToVert[3] = vertLabel;
                cells_[cell.son[2]].cellToVert[0] = vertLabel;
                cells_[cell.son[3]].cellToVert[1] = vertLabel;

                for (int i = 0; i < 4; i++)
                    verts_[vertLabel].vertToCell[i] = cell.son[i];
                
                verts_[vertLabel].vertToHEdge[0] = cells_[cell.son[0]].cellToHEdge[1][0];
                verts_[vertLabel].vertToHEdge[1] = cells_[cell.son[1]].cellToHEdge[1][0];
                verts_[vertLabel].vertToVEdge[0] = cells_[cell.son[0]].cellToVEdge[1][0];
                verts_[vertLabel].vertToVEdge[1] = cells_[cell.son[3]].cellToVEdge[1][0];

                centerVertLabel = vertLabel;

                nNewVert++;
            }

            if (newDir[0])// 检验该边是否是新分裂
            {
                const int vertLabel = vertLabelSpace_[nNewVert];
                Vector2D vertCoord = quarterPoint(0, 2);
                Vert vert(vertLabel, vertCoord);
                assignNewEntity(verts_[vertLabel], vert);
                cells_[cell.son[0]].cellToVert[3] = vertLabel;
                cells_[cell.son[3]].cellToVert[0] = vertLabel;
                verts_[vertLabel].vertToCell[0] = cell.neiCell[0][0];
                verts_[vertLabel].vertToCell[3] = cell.neiCell[0][0];
                verts_[vertLabel].vertToCell[1] = cell.son[0];
                verts_[vertLabel].vertToCell[2] = cell.son[3];

                int vEdgeL0 = cells_[cell.son[0]].cellToVEdge[0][0];
                int hEdgeL1 = cells_[cell.son[0]].cellToHEdge[1][0];
                int vEdgeL1 = cells_[cell.son[3]].cellToVEdge[0][0];
                verts_[vertLabel].neiVertices[1] = cell.cellToVert[0];
                verts_[vertLabel].neiVertices[2] = centerVertLabel;
                verts_[vertLabel].neiVertices[3] = cell.cellToVert[3];
                verts_[cell.cellToVert[0]].neiVertices[3] = vertLabel;
                verts_[centerVertLabel].neiVertices[0] = vertLabel;
                verts_[cell.cellToVert[3]].neiVertices[1] = vertLabel;
                verts_[vertLabel].vertToVEdge[0] = vEdgeL0;
                verts_[vertLabel].vertToHEdge[1] = hEdgeL1;
                verts_[vertLabel].vertToVEdge[1] = vEdgeL1;
                vEdges_[vEdgeL0].edgeToVert[0] = cell.cellToVert[0];
                vEdges_[vEdgeL0].edgeToVert[1] = vertLabel;
                hEdges_[hEdgeL1].edgeToVert[0] = vertLabel;
                hEdges_[hEdgeL1].edgeToVert[1] = centerVertLabel;
                vEdges_[vEdgeL1].edgeToVert[0] = vertLabel;
                vEdges_[vEdgeL1].edgeToVert[1] = cell.cellToVert[3];

                nNewVert++;
            }
            else
            {
                int vertLabel = cells_[cells_[cell.son[0]].neiCell[0][0]].cellToVert[2];
                int hEdgeL1 = cells_[cell.son[0]].cellToHEdge[1][0];
                cells_[cell.son[0]].cellToVert[3] = vertLabel;
                cells_[cell.son[3]].cellToVert[0] = vertLabel;
                verts_[vertLabel].vertToCell[1] = cell.son[0];
                verts_[vertLabel].vertToCell[2] = cell.son[3];
                verts_[vertLabel].neiVertices[2] = centerVertLabel;
                verts_[centerVertLabel].neiVertices[0] = vertLabel;
                verts_[vertLabel].vertToHEdge[1] = hEdgeL1;
                hEdges_[hEdgeL1].edgeToVert[0] = vertLabel;
                hEdges_[hEdgeL1].edgeToVert[1] = centerVertLabel;
            }

            if (newDir[1])// 检验该边是否是新分裂
            {
                const int vertLabel = vertLabelSpace_[nNewVert];
                Vector2D vertCoord = quarterPoint(2, 0);
                Vert vert(vertLabel, vertCoord);
                assignNewEntity(verts_[vertLabel], vert);
                cells_[cell.son[0]].cellToVert[1] = vertLabel;
                cells_[cell.son[1]].cellToVert[0] = vertLabel;
                verts_[vertLabel].vertToCell[0] = cell.neiCell[1][0];
                verts_[vertLabel].vertToCell[1] = cell.neiCell[1][0];
                verts_[vertLabel].vertToCell[2] = cell.son[1];
                verts_[vertLabel].vertToCell[3] = cell.son[0];

                int hEdgeL0 = cells_[cell.son[0]].cellToHEdge[0][0];
                int vEdgeL1 = cells_[cell.son[0]].cellToVEdge[1][0];
                int hEdgeL1 = cells_[cell.son[1]].cellToHEdge[0][0];
                verts_[vertLabel].neiVertices[0] = cell.cellToVert[0];
                verts_[vertLabel].neiVertices[2] = cell.cellToVert[1];
                verts_[vertLabel].neiVertices[3] = centerVertLabel;
                verts_[cell.cellToVert[0]].neiVertices[2] = vertLabel;
                verts_[cell.cellToVert[1]].neiVertices[0] = vertLabel;
                verts_[centerVertLabel].neiVertices[1] = vertLabel;
                verts_[vertLabel].vertToHEdge[0] = hEdgeL0;
                verts_[vertLabel].vertToVEdge[1] = vEdgeL1;
                verts_[vertLabel].vertToHEdge[1] = hEdgeL1;
                hEdges_[hEdgeL0].edgeToVert[0] = cell.cellToVert[0];
                hEdges_[hEdgeL0].edgeToVert[1] = vertLabel;
                vEdges_[vEdgeL1].edgeToVert[0] = vertLabel;
                vEdges_[vEdgeL1].edgeToVert[1] = centerVertLabel;
                hEdges_[hEdgeL1].edgeToVert[0] = vertLabel;
                hEdges_[hEdgeL1].edgeToVert[1] = cell.cellToVert[1];

                nNewVert++;
            }
            else
            {
                int vertLabel = cells_[cells_[cell.son[1]].neiCell[1][0]].cellToVert[3];
                int vEdgeL1 = cells_[cell.son[1]].cellToVEdge[0][0];
                cells_[cell.son[0]].cellToVert[1] = vertLabel;
                cells_[cell.son[1]].cellToVert[0] = vertLabel;
                verts_[vertLabel].vertToCell[2] = cell.son[1];
                verts_[vertLabel].vertToCell[3] = cell.son[0];
                verts_[vertLabel].neiVertices[3] = centerVertLabel;
                verts_[centerVertLabel].neiVertices[1] = vertLabel;
                verts_[vertLabel].vertToVEdge[1] = vEdgeL1;
                vEdges_[vEdgeL1].edgeToVert[0] = vertLabel;
                vEdges_[vEdgeL1].edgeToVert[1] = centerVertLabel;
            }

            if (newDir[2])// 检验该边是否是新分裂
            {
                const int vertLabel = vertLabelSpace_[nNewVert];
                Vector2D vertCoord = quarterPoint(4, 2);
                Vert vert(vertLabel, vertCoord);
                assignNewEntity(verts_[vertLabel], vert);
                cells_[cell.son[1]].cellToVert[2] = vertLabel;
                cells_[cell.son[2]].cellToVert[1] = vertLabel;
                verts_[vertLabel].vertToCell[1] = cell.neiCell[2][0];
                verts_[vertLabel].vertToCell[2] = cell.neiCell[2][0];
                verts_[vertLabel].vertToCell[0] = cell.son[1];
                verts_[vertLabel].vertToCell[3] = cell.son[2];

                int vEdgeL0 = cells_[cell.son[1]].cellToVEdge[1][0];
                int hEdgeL0 = cells_[cell.son[1]].cellToHEdge[1][0];
                int vEdgeL1 = cells_[cell.son[2]].cellToVEdge[1][0];
                verts_[vertLabel].neiVertices[0] = centerVertLabel;
                verts_[vertLabel].neiVertices[1] = cell.cellToVert[1];
                verts_[vertLabel].neiVertices[3] = cell.cellToVert[2];
                verts_[centerVertLabel].neiVertices[2] = vertLabel;
                verts_[cell.cellToVert[1]].neiVertices[3] = vertLabel;
                verts_[cell.cellToVert[2]].neiVertices[1] = vertLabel;
                verts_[vertLabel].vertToVEdge[0] = vEdgeL0;
                verts_[vertLabel].vertToHEdge[0] = hEdgeL0;
                verts_[vertLabel].vertToVEdge[1] = vEdgeL1;
                vEdges_[vEdgeL0].edgeToVert[0] = cell.cellToVert[1];
                vEdges_[vEdgeL0].edgeToVert[1] = vertLabel;
                hEdges_[hEdgeL0].edgeToVert[0] = centerVertLabel;
                hEdges_[hEdgeL0].edgeToVert[1] = vertLabel;
                vEdges_[vEdgeL1].edgeToVert[0] = vertLabel;
                vEdges_[vEdgeL1].edgeToVert[1] = cell.cellToVert[2];

                nNewVert++;
            }
            else
            {
                int vertLabel = cells_[cells_[cell.son[2]].neiCell[2][0]].cellToVert[0];
                int hEdgeL0 = cells_[cell.son[2]].cellToHEdge[0][0];
                cells_[cell.son[1]].cellToVert[2] = vertLabel;
                cells_[cell.son[2]].cellToVert[1] = vertLabel;
                verts_[vertLabel].vertToCell[0] = cell.son[1];
                verts_[vertLabel].vertToCell[3] = cell.son[2];
                verts_[vertLabel].neiVertices[0] = centerVertLabel;
                verts_[centerVertLabel].neiVertices[2] = vertLabel;
                verts_[vertLabel].vertToHEdge[0] = hEdgeL0;
                hEdges_[hEdgeL0].edgeToVert[0] = centerVertLabel;
                hEdges_[hEdgeL0].edgeToVert[1] = vertLabel;
            }

            if (newDir[3])// 检验该边是否是新分裂
            {
                const int vertLabel = vertLabelSpace_[nNewVert];
                Vector2D vertCoord = quarterPoint(2, 4);
                Vert vert(vertLabel, vertCoord);
                assignNewEntity(verts_[vertLabel], vert);
                cells_[cell.son[2]].cellToVert[3] = vertLabel;
                cells_[cell.son[3]].cellToVert[2] = vertLabel;
                verts_[vertLabel].vertToCell[2] = cell.neiCell[3][0];
                verts_[vertLabel].vertToCell[3] = cell.neiCell[3][0];
                verts_[vertLabel].vertToCell[0] = cell.son[3];
                verts_[vertLabel].vertToCell[1] = cell.son[2];

                int hEdgeL0 = cells_[cell.son[3]].cellToHEdge[1][0];
                int vEdgeL0 = cells_[cell.son[3]].cellToVEdge[1][0];
                int hEdgeL1 = cells_[cell.son[2]].cellToHEdge[1][0];
                verts_[vertLabel].neiVertices[0] = cell.cellToVert[3];
                verts_[vertLabel].neiVertices[1] = centerVertLabel;
                verts_[vertLabel].neiVertices[2] = cell.cellToVert[2];
                verts_[cell.cellToVert[3]].neiVertices[2] = vertLabel;
                verts_[centerVertLabel].neiVertices[3] = vertLabel;
                verts_[cell.cellToVert[2]].neiVertices[0] = vertLabel;
                verts_[vertLabel].vertToHEdge[0] = hEdgeL0;
                verts_[vertLabel].vertToVEdge[0] = vEdgeL0;
                verts_[vertLabel].vertToHEdge[1] = hEdgeL1;
                hEdges_[hEdgeL0].edgeToVert[0] = cell.cellToVert[3];
                hEdges_[hEdgeL0].edgeToVert[1] = vertLabel;
                vEdges_[vEdgeL0].edgeToVert[0] = centerVertLabel;
                vEdges_[vEdgeL0].edgeToVert[1] = vertLabel;
                hEdges_[hEdgeL1].edgeToVert[0] = vertLabel;
                hEdges_[hEdgeL1].edgeToVert[1] = cell.cellToVert[2];

                nNewVert++;
            }
            else
            {
                int vertLabel = cells_[cells_[cell.son[3]].neiCell[3][0]].cellToVert[1];
                int vEdgeL0 = cells_[cell.son[3]].cellToVEdge[1][0];
                cells_[cell.son[2]].cellToVert[3] = vertLabel;
                cells_[cell.son[3]].cellToVert[2] = vertLabel;
                verts_[vertLabel].vertToCell[0] = cell.son[3];
                verts_[vertLabel].vertToCell[1] = cell.son[2];
                verts_[vertLabel].neiVertices[1] = centerVertLabel;
                verts_[centerVertLabel].neiVertices[3] = vertLabel;
                verts_[vertLabel].vertToVEdge[0] = vEdgeL0;
                vEdges_[vEdgeL0].edgeToVert[0] = centerVertLabel;
                vEdges_[vEdgeL0].edgeToVert[1] = vertLabel;
            }
        }
    }

    vertPtr_ += nNewVert;
    hEdgePtr_ += nNewHEdge;
    vEdgePtr_ += nNewVEdge;
    cellPtr_ += nNewCell;

    if (nSkipRefine > 0)
    {
        std::cerr << "[Mesh] WARNING: skipped " << nSkipRefine
                  << " refinements at level " << currentLevel
                  << " due to maxCellNum/maxVertNum/maxEdgeNum limits.\n";
    }

    // Refresh label spaces so multi-level refinement doesn't reuse labels.
    const int remainVertNum = std::max(0, maxVertNum_ - vertPtr_);
    const int remainHEdgeNum = std::max(0, maxHEdgeNum_ - hEdgePtr_);
    const int remainVEdgeNum = std::max(0, maxVEdgeNum_ - vEdgePtr_);
    const int remainCellNum = std::max(0, maxCellNum_ - cellPtr_);

    for (int i = 0; i < remainVertNum; ++i)
        vertLabelSpace_[i] = vertPtr_ + i;
    for (int i = 0; i < remainHEdgeNum; ++i)
        hEdgeLabelSpace_[i] = hEdgePtr_ + i;
    for (int i = 0; i < remainVEdgeNum; ++i)
        vEdgeLabelSpace_[i] = vEdgePtr_ + i;
    for (int i = 0; i < remainCellNum; ++i)
        cellLabelSpace_[i] = cellPtr_ + i;

    statistics(workspace);
}

void QuadtreeMesh::submergeCell(const int currentLevel)
{submergeCellWithWorkspace(currentLevel,nullptr);}
void QuadtreeMesh::submergeCellWithWorkspace(int currentLevel,NativeAdaptWorkspace* workspace)
{
    // Recycled slots remain part of fixed-width QMC2 storage. Initialize only
    // the reset entities; never clear surviving/reactivated parent metadata.
    auto assignResetEntity = [this,workspace](auto& destination, const auto& source)
    {
        recordNativeWrite(workspace,&destination);
        destination = source;
        destination.ensureLevelSetSize(static_cast<std::size_t>(levelSetCount_));
    };
    for (int c= 0; c < cellIndexNum_; c++)
    {
        Cell&  cell = cells_[cellIndexToLabel[c]];
        if (cell.needCoarsen)
        {         
            Cell& parentCell = cells_[cell.parent];

            const bool coarsenFlag = coarseningFamilyEligible(parentCell, [&](int child) {
                const bool eligible = cells_[child].needCoarsen;
                cells_[child].needCoarsen = false;
                return eligible;
            });

            if (!coarsenFlag)
            {
                continue;
            }

            recordNativeWrite(workspace,&parentCell);recordNativeCellEdges(workspace,parentCell);
            for(int k=0;k<4;++k)recordNativeCellEdges(workspace,cells_[parentCell.son[k]]);
            findSubmergeRelatedCells(parentCell,workspace);

            Cell& sonCell0 = cells_[parentCell.son[0]];
            Cell& sonCell1 = cells_[parentCell.son[1]];
            Cell& sonCell2 = cells_[parentCell.son[2]];
            Cell& sonCell3 = cells_[parentCell.son[3]];

            if (sonCell0.neiCell[0][0] != -1)// 左
            {
                Cell& neiCell0 = cells_[sonCell0.neiCell[0][0]];
                if (neiCell0.neiCell[2][1] != -1)
                {
                    neiCell0.neiCell[2][1] = -1;
                    neiCell0.neiCell[2][0] = sonCell0.parent;
                    parentCell.neiCell[0][0] = sonCell0.neiCell[0][0];
                    parentCell.neiCell[0][1] = -1;

                    // 合并边
                    int newVEdgeL = vEdges_[sonCell0.cellToVEdge[0][0]].parent;
                    parentCell.cellToVEdge[0][0] = newVEdgeL;
                    parentCell.cellToVEdge[0][1] = -1;
                    neiCell0.cellToVEdge[1][0] = newVEdgeL;
                    neiCell0.cellToVEdge[1][1] = -1;
                    vEdges_[newVEdgeL].masterCell = sonCell0.parent;
                    vEdges_[newVEdgeL].slaveCell = sonCell0.neiCell[0][0];
                    vEdges_[newVEdgeL].son[0] = -1;
                    vEdges_[newVEdgeL].son[1] = -1;

                    vEdges_[newVEdgeL].setLabel(newVEdgeL);

                    findSubmergeRelatedVEdges(vEdges_[newVEdgeL], sonCell0, sonCell3,workspace);

                    Edge edge;
                    assignResetEntity(vEdges_[sonCell0.cellToVEdge[0][0]], edge);
                    assignResetEntity(vEdges_[sonCell3.cellToVEdge[0][0]], edge);

                    // 删除点
                    Vert& deleteVert = verts_[sonCell0.cellToVert[3]];
                    verts_[deleteVert.neiVertices[1]].neiVertices[3] = deleteVert.neiVertices[3];
                    verts_[deleteVert.neiVertices[3]].neiVertices[1] = deleteVert.neiVertices[1];
                    verts_[deleteVert.neiVertices[1]].vertToVEdge[1] = newVEdgeL;
                    verts_[deleteVert.neiVertices[3]].vertToVEdge[0] = newVEdgeL;
                    verts_[deleteVert.neiVertices[1]].vertToCell[2] = sonCell0.parent;
                    verts_[deleteVert.neiVertices[3]].vertToCell[1] = sonCell0.parent;

                    Vert vert;
                    assignResetEntity(deleteVert, vert);
                }
                else if (neiCell0.neiCell[2][0] != -1)
                {
                    neiCell0.neiCell[2][0] = sonCell0.parent;
                    cells_[sonCell3.neiCell[0][0]].neiCell[2][0] = sonCell0.parent;

                    parentCell.neiCell[0][0] = sonCell0.neiCell[0][0];
                    parentCell.neiCell[0][1] = sonCell3.neiCell[0][0];

                    parentCell.cellToVEdge[0][0] = sonCell0.cellToVEdge[0][0];
                    parentCell.cellToVEdge[0][1] = sonCell3.cellToVEdge[0][0];

                    vEdges_[sonCell0.cellToVEdge[0][0]].masterCell = sonCell0.parent;
                    vEdges_[sonCell3.cellToVEdge[0][0]].masterCell = sonCell0.parent;

                    verts_[sonCell0.cellToVert[3]].neiVertices[2] = -1;
                    verts_[sonCell0.cellToVert[3]].vertToHEdge[1] = -1;
                    verts_[sonCell0.cellToVert[3]].vertToCell[1] = sonCell0.parent;
                    verts_[sonCell0.cellToVert[3]].vertToCell[2] = sonCell0.parent;
                    verts_[sonCell0.cellToVert[0]].vertToCell[2] = sonCell0.parent;
                    verts_[sonCell3.cellToVert[3]].vertToCell[1] = sonCell0.parent;
                }
                else
                {
                    cerr << "[QuadtreeMesh] ERROR: neighbor mapping error. Line: " << __LINE__ << endl;
                    exit(-1);
                }
            }
            else
            {
                parentCell.neiCell[0][0] = -1;

                // 合并边
                int newVEdgeL = vEdges_[sonCell0.cellToVEdge[0][0]].parent;
                parentCell.cellToVEdge[0][0] = newVEdgeL;
                parentCell.cellToVEdge[0][1] = -1;
                vEdges_[newVEdgeL].masterCell = sonCell0.parent;
                vEdges_[newVEdgeL].son[0] = -1;
                vEdges_[newVEdgeL].son[1] = -1;

                vEdges_[newVEdgeL].setLabel(newVEdgeL);

                findSubmergeRelatedVEdges(vEdges_[newVEdgeL], sonCell0, sonCell3,workspace);

                Edge edge;
                assignResetEntity(vEdges_[sonCell0.cellToVEdge[0][0]], edge);
                assignResetEntity(vEdges_[sonCell3.cellToVEdge[0][0]], edge);

                // 删除点
                Vert& deleteVert = verts_[sonCell0.cellToVert[3]];
                verts_[deleteVert.neiVertices[1]].neiVertices[3] = deleteVert.neiVertices[3];
                verts_[deleteVert.neiVertices[3]].neiVertices[1] = deleteVert.neiVertices[1];
                verts_[deleteVert.neiVertices[1]].vertToVEdge[1] = newVEdgeL;
                verts_[deleteVert.neiVertices[3]].vertToVEdge[0] = newVEdgeL;
                verts_[deleteVert.neiVertices[1]].vertToCell[2] = sonCell0.parent;
                verts_[deleteVert.neiVertices[3]].vertToCell[1] = sonCell0.parent;

                Vert vert;
                assignResetEntity(deleteVert, vert);
            }

            if (sonCell1.neiCell[1][0] != -1)// 下
            {
                Cell& neiCell1 = cells_[sonCell1.neiCell[1][0]];
                if (neiCell1.neiCell[3][1] != -1)
                {
                    neiCell1.neiCell[3][1] = -1;
                    neiCell1.neiCell[3][0] = sonCell1.parent;
                    parentCell.neiCell[1][0] = sonCell1.neiCell[1][0];
                    parentCell.neiCell[1][1] = -1;

                    // 合并边
                    int newHEdgeL = hEdges_[sonCell1.cellToHEdge[0][0]].parent;
                    parentCell.cellToHEdge[0][0] = newHEdgeL;
                    parentCell.cellToHEdge[0][1] = -1;
                    neiCell1.cellToHEdge[1][0] = newHEdgeL;
                    neiCell1.cellToHEdge[1][1] = -1;
                    hEdges_[newHEdgeL].masterCell = sonCell1.parent;
                    hEdges_[newHEdgeL].slaveCell = sonCell1.neiCell[1][0];
                    hEdges_[newHEdgeL].son[0] = -1;
                    hEdges_[newHEdgeL].son[1] = -1;

                    hEdges_[newHEdgeL].setLabel(newHEdgeL);

                    findSubmergeRelatedHEdges(hEdges_[newHEdgeL], sonCell1, sonCell0,workspace);

                    Edge edge;
                    assignResetEntity(hEdges_[sonCell1.cellToHEdge[0][0]], edge);
                    assignResetEntity(hEdges_[sonCell0.cellToHEdge[0][0]], edge);

                    // 删除点
                    Vert& deleteVert = verts_[sonCell1.cellToVert[0]];
                    verts_[deleteVert.neiVertices[2]].neiVertices[0] = deleteVert.neiVertices[0];
                    verts_[deleteVert.neiVertices[0]].neiVertices[2] = deleteVert.neiVertices[2];
                    verts_[deleteVert.neiVertices[2]].vertToHEdge[0] = newHEdgeL;
                    verts_[deleteVert.neiVertices[0]].vertToHEdge[1] = newHEdgeL;
                    verts_[deleteVert.neiVertices[2]].vertToCell[3] = sonCell1.parent;
                    verts_[deleteVert.neiVertices[0]].vertToCell[2] = sonCell1.parent;

                    Vert vert;
                    assignResetEntity(deleteVert, vert);
                }
                else if (neiCell1.neiCell[3][0] != -1)
                {
                    neiCell1.neiCell[3][0] = sonCell1.parent;
                    cells_[sonCell0.neiCell[1][0]].neiCell[3][0] = sonCell1.parent;

                    parentCell.neiCell[1][0] = sonCell0.neiCell[1][0];
                    parentCell.neiCell[1][1] = sonCell1.neiCell[1][0];

                    parentCell.cellToHEdge[0][0] = sonCell0.cellToHEdge[0][0];
                    parentCell.cellToHEdge[0][1] = sonCell1.cellToHEdge[0][0];

                    hEdges_[sonCell0.cellToHEdge[0][0]].masterCell = sonCell1.parent;
                    hEdges_[sonCell1.cellToHEdge[0][0]].masterCell = sonCell1.parent;

                    verts_[sonCell1.cellToVert[0]].neiVertices[3] = -1;
                    verts_[sonCell1.cellToVert[0]].vertToVEdge[1] = -1;
                    verts_[sonCell1.cellToVert[0]].vertToCell[2] = sonCell1.parent;
                    verts_[sonCell1.cellToVert[0]].vertToCell[3] = sonCell1.parent;
                    verts_[sonCell1.cellToVert[1]].vertToCell[3] = sonCell1.parent;
                    verts_[sonCell0.cellToVert[0]].vertToCell[2] = sonCell1.parent;
                }
                else
                {
                    cerr << "[QuadtreeMesh] ERROR: neighbor mapping error. Line: " << __LINE__ << endl;
                    exit(-1);
                }
            }
            else
            {
                parentCell.neiCell[1][0] = -1;

                // 合并边
                int newHEdgeL = hEdges_[sonCell1.cellToHEdge[0][0]].parent;
                parentCell.cellToHEdge[0][0] = newHEdgeL;
                parentCell.cellToHEdge[0][1] = -1;
                hEdges_[newHEdgeL].masterCell = sonCell1.parent;
                hEdges_[newHEdgeL].son[0] = -1;
                hEdges_[newHEdgeL].son[1] = -1;

                hEdges_[newHEdgeL].setLabel(newHEdgeL);

                findSubmergeRelatedHEdges(hEdges_[newHEdgeL], sonCell1, sonCell0,workspace);

                Edge edge;
                assignResetEntity(hEdges_[sonCell1.cellToHEdge[0][0]], edge);
                assignResetEntity(hEdges_[sonCell0.cellToHEdge[0][0]], edge);

                // 删除点
                Vert& deleteVert = verts_[sonCell1.cellToVert[0]];
                verts_[deleteVert.neiVertices[2]].neiVertices[0] = deleteVert.neiVertices[0];
                verts_[deleteVert.neiVertices[0]].neiVertices[2] = deleteVert.neiVertices[2];
                verts_[deleteVert.neiVertices[2]].vertToHEdge[0] = newHEdgeL;
                verts_[deleteVert.neiVertices[0]].vertToHEdge[1] = newHEdgeL;
                verts_[deleteVert.neiVertices[2]].vertToCell[3] = sonCell1.parent;
                verts_[deleteVert.neiVertices[0]].vertToCell[2] = sonCell1.parent;

                Vert vert;
                assignResetEntity(deleteVert, vert);
            }

            if (sonCell2.neiCell[2][0] != -1)// 右
            {
                Cell& neiCell2 = cells_[sonCell2.neiCell[2][0]];
                if (neiCell2.neiCell[0][1] != -1)
                {
                    neiCell2.neiCell[0][1] = -1;
                    neiCell2.neiCell[0][0] = sonCell2.parent;
                    parentCell.neiCell[2][0] = sonCell2.neiCell[2][0];
                    parentCell.neiCell[2][1] = -1;

                    // 合并边
                    int newVEdgeL = vEdges_[sonCell2.cellToVEdge[1][0]].parent;
                    parentCell.cellToVEdge[1][0] = newVEdgeL;
                    parentCell.cellToVEdge[1][1] = -1;
                    neiCell2.cellToVEdge[0][0] = newVEdgeL;
                    neiCell2.cellToVEdge[0][1] = -1;
                    vEdges_[newVEdgeL].slaveCell = sonCell2.parent;
                    vEdges_[newVEdgeL].masterCell = sonCell2.neiCell[2][0];
                    vEdges_[newVEdgeL].son[0] = -1;
                    vEdges_[newVEdgeL].son[1] = -1;

                    vEdges_[newVEdgeL].setLabel(newVEdgeL);

                    findSubmergeRelatedVEdges(vEdges_[newVEdgeL], sonCell2, sonCell1,workspace);

                    Edge edge;
                    assignResetEntity(vEdges_[sonCell2.cellToVEdge[1][0]], edge);
                    assignResetEntity(vEdges_[sonCell1.cellToVEdge[1][0]], edge);

                    // 删除点
                    Vert& deleteVert = verts_[sonCell2.cellToVert[1]];
                    verts_[deleteVert.neiVertices[3]].neiVertices[1] = deleteVert.neiVertices[1];
                    verts_[deleteVert.neiVertices[1]].neiVertices[3] = deleteVert.neiVertices[3];
                    verts_[deleteVert.neiVertices[3]].vertToVEdge[0] = newVEdgeL;
                    verts_[deleteVert.neiVertices[1]].vertToVEdge[1] = newVEdgeL;
                    verts_[deleteVert.neiVertices[3]].vertToCell[0] = sonCell2.parent;
                    verts_[deleteVert.neiVertices[1]].vertToCell[3] = sonCell2.parent;

                    Vert vert;
                    assignResetEntity(deleteVert, vert);
                }
                else if (neiCell2.neiCell[0][0] != -1)
                {
                    neiCell2.neiCell[0][0] = sonCell2.parent;
                    cells_[sonCell1.neiCell[2][0]].neiCell[0][0] = sonCell2.parent;

                    parentCell.neiCell[2][0] = sonCell1.neiCell[2][0];
                    parentCell.neiCell[2][1] = sonCell2.neiCell[2][0];

                    parentCell.cellToVEdge[1][0] = sonCell1.cellToVEdge[1][0];
                    parentCell.cellToVEdge[1][1] = sonCell2.cellToVEdge[1][0];

                    vEdges_[sonCell1.cellToVEdge[1][0]].slaveCell = sonCell2.parent;
                    vEdges_[sonCell2.cellToVEdge[1][0]].slaveCell = sonCell2.parent;

                    verts_[sonCell2.cellToVert[1]].neiVertices[0] = -1;
                    verts_[sonCell2.cellToVert[1]].vertToHEdge[0] = -1;
                    verts_[sonCell2.cellToVert[1]].vertToCell[0] = sonCell2.parent;
                    verts_[sonCell2.cellToVert[1]].vertToCell[3] = sonCell2.parent;
                    verts_[sonCell2.cellToVert[2]].vertToCell[0] = sonCell2.parent;
                    verts_[sonCell1.cellToVert[1]].vertToCell[3] = sonCell2.parent;
                }
                else
                {
                    cerr << "[QuadtreeMesh] ERROR: neighbor mapping error. Line: " << __LINE__ << endl;
                    exit(-1);
                }
            }
            else
            {
                parentCell.neiCell[2][0] = -1;

                // 合并边
                int newVEdgeL = vEdges_[sonCell2.cellToVEdge[1][0]].parent;
                parentCell.cellToVEdge[1][0] = newVEdgeL;
                parentCell.cellToVEdge[1][1] = -1;
                vEdges_[newVEdgeL].slaveCell = sonCell2.parent;
                vEdges_[newVEdgeL].son[0] = -1;
                vEdges_[newVEdgeL].son[1] = -1;

                vEdges_[newVEdgeL].setLabel(newVEdgeL);

                findSubmergeRelatedVEdges(vEdges_[newVEdgeL], sonCell2, sonCell1,workspace);

                Edge edge;
                assignResetEntity(vEdges_[sonCell2.cellToVEdge[1][0]], edge);
                assignResetEntity(vEdges_[sonCell1.cellToVEdge[1][0]], edge);

                // 删除点
                Vert& deleteVert = verts_[sonCell2.cellToVert[1]];
                verts_[deleteVert.neiVertices[3]].neiVertices[1] = deleteVert.neiVertices[1];
                verts_[deleteVert.neiVertices[1]].neiVertices[3] = deleteVert.neiVertices[3];
                verts_[deleteVert.neiVertices[3]].vertToVEdge[0] = newVEdgeL;
                verts_[deleteVert.neiVertices[1]].vertToVEdge[1] = newVEdgeL;
                verts_[deleteVert.neiVertices[3]].vertToCell[0] = sonCell2.parent;
                verts_[deleteVert.neiVertices[1]].vertToCell[3] = sonCell2.parent;

                Vert vert;
                assignResetEntity(deleteVert, vert);
            }

            if (sonCell3.neiCell[3][0] != -1)// 上
            {
                Cell& neiCell3 = cells_[sonCell3.neiCell[3][0]];
                if (neiCell3.neiCell[1][1] != -1)
                {
                    neiCell3.neiCell[1][1] = -1;
                    neiCell3.neiCell[1][0] = sonCell3.parent;
                    parentCell.neiCell[3][0] = sonCell3.neiCell[3][0];
                    parentCell.neiCell[3][1] = -1;

                    // 合并边
                    int newHEdgeL = hEdges_[sonCell3.cellToHEdge[1][0]].parent;
                    parentCell.cellToHEdge[1][0] = newHEdgeL;
                    parentCell.cellToHEdge[1][1] = -1;
                    neiCell3.cellToHEdge[0][0] = newHEdgeL;
                    neiCell3.cellToHEdge[0][1] = -1;
                    hEdges_[newHEdgeL].slaveCell = sonCell3.parent;
                    hEdges_[newHEdgeL].masterCell = sonCell3.neiCell[3][0];
                    hEdges_[newHEdgeL].son[0] = -1;
                    hEdges_[newHEdgeL].son[1] = -1;

                    hEdges_[newHEdgeL].setLabel(newHEdgeL);

                    findSubmergeRelatedHEdges(hEdges_[newHEdgeL], sonCell3, sonCell2,workspace);

                    Edge edge;
                    assignResetEntity(hEdges_[sonCell3.cellToHEdge[1][0]], edge);
                    assignResetEntity(hEdges_[sonCell2.cellToHEdge[1][0]], edge);

                    // 删除点
                    Vert& deleteVert = verts_[sonCell3.cellToVert[2]];
                    verts_[deleteVert.neiVertices[0]].neiVertices[2] = deleteVert.neiVertices[2];
                    verts_[deleteVert.neiVertices[2]].neiVertices[0] = deleteVert.neiVertices[0];
                    verts_[deleteVert.neiVertices[0]].vertToHEdge[1] = newHEdgeL;
                    verts_[deleteVert.neiVertices[2]].vertToHEdge[0] = newHEdgeL;
                    verts_[deleteVert.neiVertices[0]].vertToCell[1] = sonCell3.parent;
                    verts_[deleteVert.neiVertices[2]].vertToCell[0] = sonCell3.parent;

                    Vert vert;
                    assignResetEntity(deleteVert, vert);
                }
                else if (neiCell3.neiCell[1][0] != -1)
                {
                    neiCell3.neiCell[1][0] = sonCell3.parent;
                    cells_[sonCell2.neiCell[3][0]].neiCell[1][0] = sonCell3.parent;

                    parentCell.neiCell[3][0] = sonCell3.neiCell[3][0];
                    parentCell.neiCell[3][1] = sonCell2.neiCell[3][0];

                    parentCell.cellToHEdge[1][0] = sonCell3.cellToHEdge[1][0];
                    parentCell.cellToHEdge[1][1] = sonCell2.cellToHEdge[1][0];

                    hEdges_[sonCell2.cellToHEdge[1][0]].slaveCell = sonCell3.parent;
                    hEdges_[sonCell3.cellToHEdge[1][0]].slaveCell = sonCell3.parent;

                    verts_[sonCell3.cellToVert[2]].neiVertices[1] = -1;
                    verts_[sonCell3.cellToVert[2]].vertToVEdge[0] = -1;
                    verts_[sonCell3.cellToVert[2]].vertToCell[0] = sonCell3.parent;
                    verts_[sonCell3.cellToVert[2]].vertToCell[1] = sonCell3.parent;
                    verts_[sonCell3.cellToVert[3]].vertToCell[1] = sonCell3.parent;
                    verts_[sonCell2.cellToVert[2]].vertToCell[0] = sonCell3.parent;
                }
                else
                {
                    cerr << "[QuadtreeMesh] ERROR: neighbor mapping error. Line: " << __LINE__ << endl;
                    exit(-1);
                }
            }
            else
            {
                parentCell.neiCell[3][0] = -1;

                // 合并边
                int newHEdgeL = hEdges_[sonCell3.cellToHEdge[1][0]].parent;
                parentCell.cellToHEdge[1][0] = newHEdgeL;
                parentCell.cellToHEdge[1][1] = -1;
                hEdges_[newHEdgeL].slaveCell = sonCell3.parent;
                hEdges_[newHEdgeL].son[0] = -1;
                hEdges_[newHEdgeL].son[1] = -1;

                hEdges_[newHEdgeL].setLabel(newHEdgeL);

                findSubmergeRelatedHEdges(hEdges_[newHEdgeL], sonCell3, sonCell2,workspace);

                Edge edge;
                assignResetEntity(hEdges_[sonCell3.cellToHEdge[1][0]], edge);
                assignResetEntity(hEdges_[sonCell2.cellToHEdge[1][0]], edge);

                // 删除点
                Vert& deleteVert = verts_[sonCell3.cellToVert[2]];
                verts_[deleteVert.neiVertices[0]].neiVertices[2] = deleteVert.neiVertices[2];
                verts_[deleteVert.neiVertices[2]].neiVertices[0] = deleteVert.neiVertices[0];
                verts_[deleteVert.neiVertices[0]].vertToHEdge[1] = newHEdgeL;
                verts_[deleteVert.neiVertices[2]].vertToHEdge[0] = newHEdgeL;
                verts_[deleteVert.neiVertices[0]].vertToCell[1] = sonCell3.parent;
                verts_[deleteVert.neiVertices[2]].vertToCell[0] = sonCell3.parent;

                Vert vert;
                assignResetEntity(deleteVert, vert);
            }

            // 中心点
            Vert vert;
            assignResetEntity(verts_[sonCell0.cellToVert[2]], vert);

            // 中心边
            Edge edge;
            assignResetEntity(hEdges_[sonCell0.cellToHEdge[1][0]], edge);
            assignResetEntity(hEdges_[sonCell1.cellToHEdge[1][0]], edge);
            assignResetEntity(vEdges_[sonCell0.cellToVEdge[1][0]], edge);
            assignResetEntity(vEdges_[sonCell3.cellToVEdge[1][0]], edge);

            cells_[cell.parent].setLabel(cell.parent);

            // 子单元
            for (int i = 0; i < 4; i++)
            {
                Cell cell;
                assignResetEntity(cells_[parentCell.son[i]], cell);
                parentCell.son[i] = -1;
            }
        }
    }
    statistics(workspace);
}

void QuadtreeMesh::findSubdivideRelatedCells(Cell& cell,NativeAdaptWorkspace* workspace)
{
    for(int k=0;k<4;++k)recordNativeWrite(workspace,&cells_[cell.son[k]]);
    for (int i = 0; i < 4; i++)
    {
        cells_[cell.son[i]].relatedPreIndex.push_back(cell.preIndex());
        cells_[cell.son[i]].relatedPreCoord.push_back(cell.centerCoord());
    }

    // Thin background lattices can have no donor on either side. Each child
    // already has its genuine parent above; never sample an absent (-1) cell.
    if (cell.neiCell[0][0] != -1)
    {
        if (cell.neiCell[0][1] == -1)
        {
            cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[0][0]].preIndex());
            cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[0][0]].centerCoord());

            cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[0][0]].preIndex());
            cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[0][0]].centerCoord());
        }
        else
        {
            if (cells_[cell.neiCell[0][1]].preIndex() == -1)
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[0][1]].parent].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[0][1]].parent].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[0][1]].parent].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[0][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[0][0]].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[0][0]].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[0][1]].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[0][1]].centerCoord());
            }
        }
    }
    else if (cell.neiCell[2][0] != -1)
    {
        if (cell.neiCell[2][1] == -1)
        {
            cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[2][0]].preIndex());
            cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[2][0]].centerCoord());

            cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[2][0]].preIndex());
            cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[2][0]].centerCoord());
        }
        else
        {
            if (cells_[cell.neiCell[2][1]].preIndex() == -1)
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[2][1]].parent].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[2][1]].parent].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[2][1]].parent].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[2][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[2][0]].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[2][0]].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[2][1]].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[2][1]].centerCoord());
            }
        }
    }

    if (cell.neiCell[1][0] != -1)
    {
        if (cell.neiCell[1][1] == -1)
        {
            cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[1][0]].preIndex());
            cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[1][0]].centerCoord());

            cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[1][0]].preIndex());
            cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[1][0]].centerCoord());
        }
        else
        {
            if (cells_[cell.neiCell[1][1]].preIndex() == -1)
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[1][1]].parent].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[1][1]].parent].centerCoord());

                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[1][1]].parent].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[1][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[1][0]].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[1][0]].centerCoord());

                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[1][1]].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[1][1]].centerCoord());
            }
        }
    }
    else if (cell.neiCell[3][0] != -1)
    {
        if (cell.neiCell[3][1] == -1)
        {
            cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[3][0]].preIndex());
            cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[3][0]].centerCoord());

            cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[3][0]].preIndex());
            cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[3][0]].centerCoord());
        }
        else
        {
            if (cells_[cell.neiCell[3][1]].preIndex() == -1)
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[3][1]].parent].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[3][1]].parent].centerCoord());

                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[3][1]].parent].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[3][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[0]].relatedPreIndex.push_back(cells_[cell.neiCell[3][0]].preIndex());
                cells_[cell.son[0]].relatedPreCoord.push_back(cells_[cell.neiCell[3][0]].centerCoord());

                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[3][1]].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[3][1]].centerCoord());
            }
        }
    }

    if (cell.neiCell[2][0] != -1)
    {
        if (cell.neiCell[2][1] == -1)
        {
            cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[2][0]].preIndex());
            cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[2][0]].centerCoord());

            cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[2][0]].preIndex());
            cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[2][0]].centerCoord());
        }
        else 
        {
            if (cells_[cell.neiCell[2][1]].preIndex() == -1)
            {
                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[2][1]].parent].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[2][1]].parent].centerCoord());

                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[2][1]].parent].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[2][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[2][0]].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[2][0]].centerCoord());

                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[2][1]].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[2][1]].centerCoord());
            }
        }
    }
    else if (cell.neiCell[0][0] != -1)
    {
        if (cell.neiCell[0][1] == -1)
        {
            cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[0][0]].preIndex());
            cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[0][0]].centerCoord());

            cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[0][0]].preIndex());
            cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[0][0]].centerCoord());
        }
        else
        {
            if (cells_[cell.neiCell[0][1]].preIndex() == -1)
            {
                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[0][1]].parent].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[0][1]].parent].centerCoord());

                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[0][1]].parent].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[0][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[1]].relatedPreIndex.push_back(cells_[cell.neiCell[0][0]].preIndex());
                cells_[cell.son[1]].relatedPreCoord.push_back(cells_[cell.neiCell[0][0]].centerCoord());

                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[0][1]].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[0][1]].centerCoord());
            }
        }
    }

    if (cell.neiCell[3][0] != -1)
    {
        if (cell.neiCell[3][1] == -1)
        {
            cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[3][0]].preIndex());
            cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[3][0]].centerCoord());

            cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[3][0]].preIndex());
            cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[3][0]].centerCoord());
        }
        else
        {
            if (cells_[cell.neiCell[3][1]].preIndex() == -1)
            {
                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[3][1]].parent].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[3][1]].parent].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[3][1]].parent].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[3][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[3][0]].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[3][0]].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[3][1]].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[3][1]].centerCoord());
            }
        }
    }
    else if (cell.neiCell[1][0] != -1)
    {
        if (cell.neiCell[1][1] == -1)
        {
            cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[1][0]].preIndex());
            cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[1][0]].centerCoord());

            cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[1][0]].preIndex());
            cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[1][0]].centerCoord());
        }
        else
        {
            if (cells_[cell.neiCell[1][1]].preIndex() == -1)
            {
                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[1][1]].parent].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[1][1]].parent].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cells_[cell.neiCell[1][1]].parent].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cells_[cell.neiCell[1][1]].parent].centerCoord());
            }
            else
            {
                cells_[cell.son[2]].relatedPreIndex.push_back(cells_[cell.neiCell[1][0]].preIndex());
                cells_[cell.son[2]].relatedPreCoord.push_back(cells_[cell.neiCell[1][0]].centerCoord());

                cells_[cell.son[3]].relatedPreIndex.push_back(cells_[cell.neiCell[1][1]].preIndex());
                cells_[cell.son[3]].relatedPreCoord.push_back(cells_[cell.neiCell[1][1]].centerCoord());
            }
        }
    }
}

void QuadtreeMesh::findSubmergeRelatedCells(Cell& cell,NativeAdaptWorkspace* workspace)
{
    recordNativeWrite(workspace,&cell);
    cell.relatedPreIndex.clear();
    cell.relatedPreCoord.clear();
    for (int i = 0; i < 4; i++)
    {
        cell.relatedPreIndex.push_back(cells_[cell.son[i]].preIndex());
        cell.relatedPreCoord.push_back(cells_[cell.son[i]].centerCoord());
    }
}

void QuadtreeMesh::findSubdivideRelatedHEdges(const Cell& cell, Edge& hEdge0, Edge& hEdge1, Edge& hEdge2, const int direction, const bool isRightSide,NativeAdaptWorkspace* workspace)
{
    recordNativeWrite(workspace,&hEdge0);recordNativeWrite(workspace,&hEdge1);recordNativeWrite(workspace,&hEdge2);
    // A real one-column/one-row background cell can have NO donor on
    // either opposite side. Never interpret the exterior -1 as a cell label.
    // Retain only actual old edge samples (or their one-level old parent).
    if(cell.neiCell[0][0]==-1&&cell.neiCell[2][0]==-1)
    {
        if(cell.neiCell[0][1]!=-1||cell.neiCell[2][1]!=-1)
            throw std::invalid_argument("QuadtreeMesh: incomplete thin-H neighbor pair");
        const auto appendOwn=[&](Edge& target,int side){
            for(int slot=0;slot<2;++slot){
                const int label=cell.cellToHEdge[side][slot];if(label==-1)continue;
                if(label<0||label>=hEdgeLabelCapacity_||hEdges_[label].label()!=label)
                    throw std::invalid_argument("QuadtreeMesh: invalid thin-H own edge label");
                const Edge* old=&hEdges_[label];
                if(old->preIndex()==-1){
                    const int parent=old->parent;
                    if(parent<0||parent>=hEdgeLabelCapacity_||hEdges_[parent].label()!=parent)
                        throw std::invalid_argument("QuadtreeMesh: missing thin-H old edge parent");
                    old=&hEdges_[parent];
                }
                if(old->preIndex()<0)
                    throw std::invalid_argument("QuadtreeMesh: thin-H edge has no actual old source");
                const auto found=std::find(target.relatedPreIndex.begin(),target.relatedPreIndex.end(),old->preIndex());
                if(found==target.relatedPreIndex.end()){target.relatedPreIndex.push_back(old->preIndex());target.relatedPreCoord.push_back(old->midpointCoord());}
            }
        };
        appendOwn(hEdge1,0);appendOwn(hEdge1,1);
        // Legacy edge fields average <3 samples. The internal midpoint uses
        // both normal-side points; boundary child edges use only their own
        // side, preserving constants and the identifiable normal linear field.
        if(hEdge0.relatedPreIndex.empty())appendOwn(hEdge0,0);
        if(hEdge2.relatedPreIndex.empty())appendOwn(hEdge2,1);
        if(hEdge0.relatedPreIndex.empty()||hEdge1.relatedPreIndex.empty()||hEdge2.relatedPreIndex.empty())
            throw std::invalid_argument("QuadtreeMesh: empty thin-H actual edge stencil");
        return;
    }
    if (hEdges_[cell.cellToHEdge[0][0]].preIndex() == -1)
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[hEdges_[cell.cellToHEdge[0][0]].parent].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[hEdges_[cell.cellToHEdge[0][0]].parent].midpointCoord());
    }
    else if (isRightSide && cell.cellToHEdge[0][1] != -1)
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[cell.cellToHEdge[0][1]].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[cell.cellToHEdge[0][1]].midpointCoord());
    }
    else
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[cell.cellToHEdge[0][0]].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[cell.cellToHEdge[0][0]].midpointCoord());
    }

    if (hEdges_[cell.cellToHEdge[1][0]].preIndex() == -1)
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[hEdges_[cell.cellToHEdge[1][0]].parent].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[hEdges_[cell.cellToHEdge[1][0]].parent].midpointCoord());
    }
    else if (isRightSide && cell.cellToHEdge[1][1] != -1)
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[cell.cellToHEdge[1][1]].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[cell.cellToHEdge[1][1]].midpointCoord());
    }
    else
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[cell.cellToHEdge[1][0]].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[cell.cellToHEdge[1][0]].midpointCoord());
    }

    Cell& neiCell0 = cells_[cell.neiCell[direction][0]];
    if (cell.neiCell[direction][1] != -1)
    {
        Cell& neiCell1 = cells_[cell.neiCell[direction][1]];
        if (hEdges_[neiCell0.cellToHEdge[1][0]].preIndex() == -1)
        {
            if (hEdges_[neiCell1.cellToHEdge[1][0]].preIndex() == -1)
            {
                hEdge1.relatedPreIndex.push_back(hEdges_[hEdges_[neiCell1.cellToHEdge[1][0]].parent].preIndex());
                hEdge1.relatedPreCoord.push_back(hEdges_[hEdges_[neiCell1.cellToHEdge[1][0]].parent].midpointCoord());
            }
            else if (!isRightSide && neiCell1.cellToHEdge[1][1] != -1)
            {
                hEdge1.relatedPreIndex.push_back(hEdges_[neiCell1.cellToHEdge[1][1]].preIndex());
                hEdge1.relatedPreCoord.push_back(hEdges_[neiCell1.cellToHEdge[1][1]].midpointCoord());
            }
            else
            {
                hEdge1.relatedPreIndex.push_back(hEdges_[neiCell1.cellToHEdge[1][0]].preIndex());
                hEdge1.relatedPreCoord.push_back(hEdges_[neiCell1.cellToHEdge[1][0]].midpointCoord());
            }
        }
        else
        {
            hEdge1.relatedPreIndex.push_back(hEdges_[neiCell0.cellToHEdge[1][0]].preIndex());
            hEdge1.relatedPreCoord.push_back(hEdges_[neiCell0.cellToHEdge[1][0]].midpointCoord());
        }
    }
    else
    {
        if (hEdges_[neiCell0.cellToHEdge[1][0]].preIndex() == -1)
        {
            hEdge1.relatedPreIndex.push_back(hEdges_[hEdges_[neiCell0.cellToHEdge[1][0]].parent].preIndex());
            hEdge1.relatedPreCoord.push_back(hEdges_[hEdges_[neiCell0.cellToHEdge[1][0]].parent].midpointCoord());
        }
        else if (!isRightSide && neiCell0.cellToHEdge[1][1] != -1)
        {
            hEdge1.relatedPreIndex.push_back(hEdges_[neiCell0.cellToHEdge[1][1]].preIndex());
            hEdge1.relatedPreCoord.push_back(hEdges_[neiCell0.cellToHEdge[1][1]].midpointCoord());
        }
        else
        {
            hEdge1.relatedPreIndex.push_back(hEdges_[neiCell0.cellToHEdge[1][0]].preIndex());
            hEdge1.relatedPreCoord.push_back(hEdges_[neiCell0.cellToHEdge[1][0]].midpointCoord());
        }
    }
    
    if (hEdges_[neiCell0.cellToHEdge[0][0]].preIndex() == -1)
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[hEdges_[neiCell0.cellToHEdge[0][0]].parent].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[hEdges_[neiCell0.cellToHEdge[0][0]].parent].midpointCoord());
    }
    else if (!isRightSide && neiCell0.cellToHEdge[0][1] != -1)
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[neiCell0.cellToHEdge[0][1]].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[neiCell0.cellToHEdge[0][1]].midpointCoord());
    }
    else
    {
        hEdge1.relatedPreIndex.push_back(hEdges_[neiCell0.cellToHEdge[0][0]].preIndex());
        hEdge1.relatedPreCoord.push_back(hEdges_[neiCell0.cellToHEdge[0][0]].midpointCoord());
    }

    if (hEdge0.relatedPreIndex.empty())
    {
        hEdge0.relatedPreIndex = hEdge1.relatedPreIndex;
        hEdge0.relatedPreCoord = hEdge1.relatedPreCoord;
    }
    
    if (hEdge2.relatedPreIndex.empty())
    {
        hEdge2.relatedPreIndex = hEdge1.relatedPreIndex;
        hEdge2.relatedPreCoord = hEdge1.relatedPreCoord;
        if (cell.neiCell[direction][1] != -1)
        {
            Cell& neiCell1 = cells_[cell.neiCell[direction][1]];
            if (hEdges_[neiCell1.cellToHEdge[1][0]].preIndex() == -1)
            {
                hEdge2.relatedPreIndex.back() = hEdges_[hEdges_[neiCell1.cellToHEdge[1][0]].parent].preIndex();
                hEdge2.relatedPreCoord.back() = hEdges_[hEdges_[neiCell1.cellToHEdge[1][0]].parent].midpointCoord();
            }
            else if (!isRightSide && neiCell1.cellToHEdge[1][1] != -1)
            {
                if (hEdges_[neiCell1.cellToHEdge[1][1]].preIndex() != -1)
                {
                    hEdge2.relatedPreIndex.back() = hEdges_[neiCell1.cellToHEdge[1][1]].preIndex();
                    hEdge2.relatedPreCoord.back() = hEdges_[neiCell1.cellToHEdge[1][1]].midpointCoord();
                }
                else
                {
                    hEdge2.relatedPreIndex.back() = hEdges_[neiCell1.cellToHEdge[1][0]].preIndex();
                    hEdge2.relatedPreCoord.back() = hEdges_[neiCell1.cellToHEdge[1][0]].midpointCoord();
                }
            }
            else
            {
                hEdge2.relatedPreIndex.back() = hEdges_[neiCell1.cellToHEdge[1][0]].preIndex();
                hEdge2.relatedPreCoord.back() = hEdges_[neiCell1.cellToHEdge[1][0]].midpointCoord();
            }
        }
    }
}

void QuadtreeMesh::findSubdivideRelatedVEdges(const Cell& cell, Edge& vEdge0, Edge& vEdge1, Edge& vEdge2, const int direction, const bool isUpSide,NativeAdaptWorkspace* workspace)
{
    recordNativeWrite(workspace,&vEdge0);recordNativeWrite(workspace,&vEdge1);recordNativeWrite(workspace,&vEdge2);
    // A real one-column/one-row background cell can have NO donor on
    // either opposite side. Never interpret the exterior -1 as a cell label.
    // Retain only actual old edge samples (or their one-level old parent).
    if(cell.neiCell[1][0]==-1&&cell.neiCell[3][0]==-1)
    {
        if(cell.neiCell[1][1]!=-1||cell.neiCell[3][1]!=-1)
            throw std::invalid_argument("QuadtreeMesh: incomplete thin-V neighbor pair");
        const auto appendOwn=[&](Edge& target,int side){
            for(int slot=0;slot<2;++slot){
                const int label=cell.cellToVEdge[side][slot];if(label==-1)continue;
                if(label<0||label>=vEdgeLabelCapacity_||vEdges_[label].label()!=label)
                    throw std::invalid_argument("QuadtreeMesh: invalid thin-V own edge label");
                const Edge* old=&vEdges_[label];
                if(old->preIndex()==-1){
                    const int parent=old->parent;
                    if(parent<0||parent>=vEdgeLabelCapacity_||vEdges_[parent].label()!=parent)
                        throw std::invalid_argument("QuadtreeMesh: missing thin-V old edge parent");
                    old=&vEdges_[parent];
                }
                if(old->preIndex()<0)
                    throw std::invalid_argument("QuadtreeMesh: thin-V edge has no actual old source");
                const auto found=std::find(target.relatedPreIndex.begin(),target.relatedPreIndex.end(),old->preIndex());
                if(found==target.relatedPreIndex.end()){target.relatedPreIndex.push_back(old->preIndex());target.relatedPreCoord.push_back(old->midpointCoord());}
            }
        };
        appendOwn(vEdge1,0);appendOwn(vEdge1,1);
        // Legacy edge fields average <3 samples. The internal midpoint uses
        // both normal-side points; boundary child edges use only their own
        // side, preserving constants and the identifiable normal linear field.
        if(vEdge0.relatedPreIndex.empty())appendOwn(vEdge0,0);
        if(vEdge2.relatedPreIndex.empty())appendOwn(vEdge2,1);
        if(vEdge0.relatedPreIndex.empty()||vEdge1.relatedPreIndex.empty()||vEdge2.relatedPreIndex.empty())
            throw std::invalid_argument("QuadtreeMesh: empty thin-V actual edge stencil");
        return;
    }
    if (vEdges_[cell.cellToVEdge[0][0]].preIndex() == -1)
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[vEdges_[cell.cellToVEdge[0][0]].parent].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[vEdges_[cell.cellToVEdge[0][0]].parent].midpointCoord());
    }
    else if (isUpSide && cell.cellToVEdge[0][1] != -1)
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[cell.cellToVEdge[0][1]].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[cell.cellToVEdge[0][1]].midpointCoord());
    }
    else
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[cell.cellToVEdge[0][0]].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[cell.cellToVEdge[0][0]].midpointCoord());
    }

    if (vEdges_[cell.cellToVEdge[1][0]].preIndex() == -1)
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[vEdges_[cell.cellToVEdge[1][0]].parent].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[vEdges_[cell.cellToVEdge[1][0]].parent].midpointCoord());
    }
    else if (isUpSide && cell.cellToVEdge[1][1] != -1)
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[cell.cellToVEdge[1][1]].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[cell.cellToVEdge[1][1]].midpointCoord());
    }
    else
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[cell.cellToVEdge[1][0]].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[cell.cellToVEdge[1][0]].midpointCoord());
    }

    Cell& neiCell0 = cells_[cell.neiCell[direction][0]];
    if (cell.neiCell[direction][1] != -1)
    {
        Cell& neiCell1 = cells_[cell.neiCell[direction][1]];
        if (vEdges_[neiCell0.cellToVEdge[1][0]].preIndex() == -1)
        {
            if (vEdges_[neiCell1.cellToVEdge[1][0]].preIndex() == -1)
            {
                vEdge1.relatedPreIndex.push_back(vEdges_[vEdges_[neiCell1.cellToVEdge[1][0]].parent].preIndex());
                vEdge1.relatedPreCoord.push_back(vEdges_[vEdges_[neiCell1.cellToVEdge[1][0]].parent].midpointCoord());
            }
            else if (!isUpSide && neiCell1.cellToVEdge[1][1] != -1)
            {
                vEdge1.relatedPreIndex.push_back(vEdges_[neiCell1.cellToVEdge[1][1]].preIndex());
                vEdge1.relatedPreCoord.push_back(vEdges_[neiCell1.cellToVEdge[1][1]].midpointCoord());
            }
            else
            {
                vEdge1.relatedPreIndex.push_back(vEdges_[neiCell1.cellToVEdge[1][0]].preIndex());
                vEdge1.relatedPreCoord.push_back(vEdges_[neiCell1.cellToVEdge[1][0]].midpointCoord());
            }
        }
        else
        {
            vEdge1.relatedPreIndex.push_back(vEdges_[neiCell0.cellToVEdge[1][0]].preIndex());
            vEdge1.relatedPreCoord.push_back(vEdges_[neiCell0.cellToVEdge[1][0]].midpointCoord());
        }
    }
    else
    {
        if (vEdges_[neiCell0.cellToVEdge[1][0]].preIndex() == -1)
        {
            vEdge1.relatedPreIndex.push_back(vEdges_[vEdges_[neiCell0.cellToVEdge[1][0]].parent].preIndex());
            vEdge1.relatedPreCoord.push_back(vEdges_[vEdges_[neiCell0.cellToVEdge[1][0]].parent].midpointCoord());
        }
        else if (!isUpSide && neiCell0.cellToVEdge[1][1] != -1)
        {
            vEdge1.relatedPreIndex.push_back(vEdges_[neiCell0.cellToVEdge[1][1]].preIndex());
            vEdge1.relatedPreCoord.push_back(vEdges_[neiCell0.cellToVEdge[1][1]].midpointCoord());
        }
        else
        {
            vEdge1.relatedPreIndex.push_back(vEdges_[neiCell0.cellToVEdge[1][0]].preIndex());
            vEdge1.relatedPreCoord.push_back(vEdges_[neiCell0.cellToVEdge[1][0]].midpointCoord());
        }
    }

    if (vEdges_[neiCell0.cellToVEdge[0][0]].preIndex() == -1)
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[vEdges_[neiCell0.cellToVEdge[0][0]].parent].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[vEdges_[neiCell0.cellToVEdge[0][0]].parent].midpointCoord());
    }
    else if (!isUpSide && neiCell0.cellToVEdge[0][1] != -1)
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[neiCell0.cellToVEdge[0][1]].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[neiCell0.cellToVEdge[0][1]].midpointCoord());
    }
    else
    {
        vEdge1.relatedPreIndex.push_back(vEdges_[neiCell0.cellToVEdge[0][0]].preIndex());
        vEdge1.relatedPreCoord.push_back(vEdges_[neiCell0.cellToVEdge[0][0]].midpointCoord());
    }

    if (vEdge0.relatedPreIndex.empty())
    {
        vEdge0.relatedPreIndex = vEdge1.relatedPreIndex;
        vEdge0.relatedPreCoord = vEdge1.relatedPreCoord;
    }

    if (vEdge2.relatedPreIndex.empty())
    {
        vEdge2.relatedPreIndex = vEdge1.relatedPreIndex;
        vEdge2.relatedPreCoord = vEdge1.relatedPreCoord;
        if (cell.neiCell[direction][1] != -1)
        {
            Cell& neiCell1 = cells_[cell.neiCell[direction][1]];
            if (vEdges_[neiCell1.cellToVEdge[1][0]].preIndex() == -1)
            {
                vEdge2.relatedPreIndex.back() = vEdges_[vEdges_[neiCell1.cellToVEdge[1][0]].parent].preIndex();
                vEdge2.relatedPreCoord.back() = vEdges_[vEdges_[neiCell1.cellToVEdge[1][0]].parent].midpointCoord();
            }
            else if (!isUpSide && neiCell1.cellToVEdge[1][1] != -1)
            {
                if (vEdges_[neiCell1.cellToVEdge[1][1]].preIndex() != -1)
                {
                    vEdge2.relatedPreIndex.back() = vEdges_[neiCell1.cellToVEdge[1][1]].preIndex();
                    vEdge2.relatedPreCoord.back() = vEdges_[neiCell1.cellToVEdge[1][1]].midpointCoord();
                }
                else
                {
                    vEdge2.relatedPreIndex.back() = vEdges_[neiCell1.cellToVEdge[1][0]].preIndex();
                    vEdge2.relatedPreCoord.back() = vEdges_[neiCell1.cellToVEdge[1][0]].midpointCoord();
                }
            }
            else
            {
                vEdge2.relatedPreIndex.back() = vEdges_[neiCell1.cellToVEdge[1][0]].preIndex();
                vEdge2.relatedPreCoord.back() = vEdges_[neiCell1.cellToVEdge[1][0]].midpointCoord();
            }
        }
    }
}

void QuadtreeMesh::findSubmergeRelatedHEdges(Edge& hEdge, const Cell& sonCellA, const Cell& sonCellB,NativeAdaptWorkspace* workspace)
{
    recordNativeWrite(workspace,&hEdge);
    hEdge.relatedPreIndex.clear();
    hEdge.relatedPreCoord.clear();
    hEdge.relatedPreIndex.push_back(hEdges_[sonCellA.cellToHEdge[0][0]].preIndex());
    hEdge.relatedPreCoord.push_back(hEdges_[sonCellA.cellToHEdge[0][0]].midpointCoord());
    hEdge.relatedPreIndex.push_back(hEdges_[sonCellA.cellToHEdge[1][0]].preIndex());
    hEdge.relatedPreCoord.push_back(hEdges_[sonCellA.cellToHEdge[1][0]].midpointCoord());
    hEdge.relatedPreIndex.push_back(hEdges_[sonCellB.cellToHEdge[0][0]].preIndex());
    hEdge.relatedPreCoord.push_back(hEdges_[sonCellB.cellToHEdge[0][0]].midpointCoord());
    hEdge.relatedPreIndex.push_back(hEdges_[sonCellB.cellToHEdge[1][0]].preIndex());
    hEdge.relatedPreCoord.push_back(hEdges_[sonCellB.cellToHEdge[1][0]].midpointCoord());
}

void QuadtreeMesh::findSubmergeRelatedVEdges(Edge& vEdge, const Cell& sonCellA, const Cell& sonCellB,NativeAdaptWorkspace* workspace)
{
    recordNativeWrite(workspace,&vEdge);
    vEdge.relatedPreIndex.clear();
    vEdge.relatedPreCoord.clear();
    vEdge.relatedPreIndex.push_back(vEdges_[sonCellA.cellToVEdge[0][0]].preIndex());
    vEdge.relatedPreCoord.push_back(vEdges_[sonCellA.cellToVEdge[0][0]].midpointCoord());
    vEdge.relatedPreIndex.push_back(vEdges_[sonCellA.cellToVEdge[1][0]].preIndex());
    vEdge.relatedPreCoord.push_back(vEdges_[sonCellA.cellToVEdge[1][0]].midpointCoord());
    vEdge.relatedPreIndex.push_back(vEdges_[sonCellB.cellToVEdge[0][0]].preIndex());
    vEdge.relatedPreCoord.push_back(vEdges_[sonCellB.cellToVEdge[0][0]].midpointCoord());
    vEdge.relatedPreIndex.push_back(vEdges_[sonCellB.cellToVEdge[1][0]].preIndex());
    vEdge.relatedPreCoord.push_back(vEdges_[sonCellB.cellToVEdge[1][0]].midpointCoord());
}

void QuadtreeMesh::statistics()
{statistics(nullptr);}
void QuadtreeMesh::statistics(NativeAdaptWorkspace* workspace)
{
    bool completed=false;
    if(refreshNativeMembership(workspace)){
        // On optional-pressure allocation failure, free the workspace first
        // and finish only original maintenance from this ALREADY mutated tree.
        NativeMaintenanceEnvironment attempt;
        try{
            workspace->vertices.exportIndices(vertIndexNum_,vertIndexToLabel,nativeAdaptationReport_);
            workspace->horizontal.exportIndices(hEdgeIndexNum_,hEdgeIndexToLabel,nativeAdaptationReport_);
            workspace->vertical.exportIndices(vEdgeIndexNum_,vEdgeIndexToLabel,nativeAdaptationReport_);
            workspace->cells.exportIndices(cellIndexNum_,cellIndexToLabel,nativeAdaptationReport_);
            completed=true;attempt.dismiss();
        }catch(const std::bad_alloc&){workspace->fail();}
        catch(const std::length_error&){workspace->fail();}
    }
    if(!completed){
        count(verts_, vertPtr_, vertIndexNum_, vertIndexToLabel);
        count(hEdges_, hEdgePtr_, hEdgeIndexNum_, hEdgeIndexToLabel);
        count(vEdges_, vEdgePtr_, vEdgeIndexNum_, vEdgeIndexToLabel);
        count(cells_, cellPtr_, cellIndexNum_, cellIndexToLabel);
    }

    nCellsPerLevel.clear();
    nCellsPerLevel.resize(maxLevel_ + 1, 0);
    for (int i = 0; i < cellIndexNum_; i++)
    {
        const int cellLabel = cellIndexToLabel[i];

        nCellsPerLevel[cells_[cellLabel].level()]++;
    }
}

void QuadtreeMesh::output(const string& name) const
{
    cout << "Writing mesh..." << endl;

    ostringstream fileName;
	fileName << "mesh_" << name;
	const auto filePath = (OutputLayout::fields() / (fileName.str() + ".dat")).string();

	ofstream fp(filePath, ofstream::out);
	if (!fp.is_open())
	{
        cerr << filePath << " write mesh file failed!" << endl;
		exit(-1);
	}

	fp << "Title = \"Quadtree Mesh" << "\"" << endl;
	fp << "Variables = \"X\",\"Y\"" << endl;
	fp << "ZONE N = " << vertIndexNum_ << ", E = " << cellIndexNum_
		<< ", F = FEPOINT, ET = QUADRILATERAL" << endl;

	for (int i = 0; i < vertIndexNum_; i++)
    {
        const int vertLabel = vertIndexToLabel[i];
        const Vert& vert = verts_[vertLabel];
        fp << verts_[vertLabel].coord().x << " " << verts_[vertLabel].coord().y << endl;
    }

	for (int i = 0; i < cellIndexNum_; i++)
    {
        const int cellLabel = cellIndexToLabel[i];
        const Cell& cell = cells_[cellLabel];
        for (int j = 0; j < 4; j++)
        {
            const Vert& vert = verts_[cell.cellToVert[j]];
            fp << vert.index() + 1 << " ";
        }
        fp << endl;
    }

	fp.close();
    cout << "Mesh output done." << endl << endl;
}

void QuadtreeMesh::writeCache(const string& filePath) const
{
    cout << "Writing mesh cache..." << endl;

    ostringstream fileName;
	fileName << filePath << "mesh_cache.bin";

    // QMC2 stores a fixed number of LS entries per allocated slot (including
    // inactive/free slots), without vector length tags for these arrays. Check
    // every source before opening/truncating the destination file.
    const auto levelSets = static_cast<std::size_t>(levelSetCount_);
    auto requireStorage = [&](bool valid, const char* kind, int label)
    {
        if (!valid)
            throw std::runtime_error("Quadtree mesh cache '" + fileName.str() +
                "' cannot be written: " + kind + " " + std::to_string(label) +
                " has inconsistent level-set storage; no cache file was created or changed.");
    };
    requireStorage(levelSetCount_ >= 0, "mesh", -1);
    for (int label = 0; label < maxVertNum_; ++label)
    {
        const auto& vertex = verts_.peek(label);
        requireStorage(vertex.onCutEdge.size() == levelSets && vertex.onCutCell.size() == levelSets &&
                       vertex.nearestNodeCoord.size() == levelSets && vertex.nearestNodeLabel.size() == levelSets,
                       "vertex", label);
    }
    for (int label = 0; label < maxHEdgeNum_; ++label)
        requireStorage(hEdges_.peek(label).isCut.size() == levelSets, "horizontal edge", label);
    for (int label = 0; label < maxVEdgeNum_; ++label)
        requireStorage(vEdges_.peek(label).isCut.size() == levelSets, "vertical edge", label);
    for (int label = 0; label < maxCellNum_; ++label)
        requireStorage(cells_.peek(label).isCut.size() == levelSets, "cell", label);

    auto cache = FileIO::openOutputFile(fileName.str(), ios::binary | ios::out);

    const char magic[4] = {'Q', 'M', 'C', '2'};
    const std::uint32_t version = 2;
    MeshCacheIO::writeRawArray(cache, magic, 4);
    MeshCacheIO::writeRaw(cache, version);

    MeshCacheIO::writeRaw(cache, levelSetCount_);
    MeshCacheIO::writeRaw(cache, vertIndexNum_);
    MeshCacheIO::writeRaw(cache, hEdgeIndexNum_);
    MeshCacheIO::writeRaw(cache, vEdgeIndexNum_);
    MeshCacheIO::writeRaw(cache, cellIndexNum_);

    MeshCacheIO::writeRawArray(cache, vertIndexToLabel, static_cast<std::size_t>(vertIndexNum_));
    MeshCacheIO::writeRawArray(cache, hEdgeIndexToLabel, static_cast<std::size_t>(hEdgeIndexNum_));
    MeshCacheIO::writeRawArray(cache, vEdgeIndexToLabel, static_cast<std::size_t>(vEdgeIndexNum_));
    MeshCacheIO::writeRawArray(cache, cellIndexToLabel, static_cast<std::size_t>(cellIndexNum_));

    MeshCacheIO::writeRaw(cache, maxLevel_);
    MeshCacheIO::writeRaw(cache, maxVertNum_);
    MeshCacheIO::writeRaw(cache, maxHEdgeNum_);
    MeshCacheIO::writeRaw(cache, maxVEdgeNum_);
    MeshCacheIO::writeRaw(cache, maxCellNum_);

    MeshCacheIO::writeRaw(cache, vertPtr_);
    MeshCacheIO::writeRaw(cache, hEdgePtr_);
    MeshCacheIO::writeRaw(cache, vEdgePtr_);
    MeshCacheIO::writeRaw(cache, cellPtr_);

    // QMC2 retains full-capacity arrays. Only the free prefix is meaningful;
    // never read old/stale suffix bytes or mutate the live allocation order.
    // Count occupied parents too, and include holes below the high-water mark.
    const auto writeLabels = [&cache](const auto& objects, const int* labels, int capacity)
    {
        int freeCount = 0;
        for (int i = 0; i < capacity; ++i)
            if (objects.peek(i).label() == -1) ++freeCount;
        MeshCacheIO::writeRawArray(cache, labels, static_cast<std::size_t>(freeCount));
        const std::array<int, 256> zeros{};
        for (int remaining = capacity - freeCount; remaining > 0;)
        {
            const int count = std::min(remaining, static_cast<int>(zeros.size()));
            MeshCacheIO::writeRawArray(cache, zeros.data(), static_cast<std::size_t>(count));
            remaining -= count;
        }
    };
    writeLabels(verts_, vertLabelSpace_, maxVertNum_);
    writeLabels(hEdges_, hEdgeLabelSpace_, maxHEdgeNum_);
    writeLabels(vEdges_, vEdgeLabelSpace_, maxVEdgeNum_);
    writeLabels(cells_, cellLabelSpace_, maxCellNum_);

    for (int i = 0; i < maxVertNum_; ++i)
        MeshCacheIO::writeVert(cache, verts_.peek(i), levelSetCount_);
    for (int i = 0; i < maxHEdgeNum_; ++i)
        MeshCacheIO::writeEdge(cache, hEdges_.peek(i), levelSetCount_);
    for (int i = 0; i < maxVEdgeNum_; ++i)
        MeshCacheIO::writeEdge(cache, vEdges_.peek(i), levelSetCount_);
    for (int i = 0; i < maxCellNum_; ++i)
        MeshCacheIO::writeCell(cache, cells_.peek(i), levelSetCount_);

    MeshCacheIO::writeVector(cache, minLevel_neg);
    MeshCacheIO::writeVector(cache, minLevel_pos);
    MeshCacheIO::writeVector(cache, maxLevel_local);
    MeshCacheIO::writeVector(cache, nCellsPerLevel);

    FileIO::finishOutputFile(cache, fileName.str());
    cout << "Mesh cache written." << endl;
}

void QuadtreeMesh::readCache(const std::string& filePath)
{
    const std::string name = filePath + "mesh_cache.bin";
    try
    {
        std::ifstream input(name, std::ios::binary | std::ios::in | std::ios::ate);
        if (!input) throw std::runtime_error("file is missing or cannot be opened");
        const auto end = input.tellg();
        if (end < std::streampos(0)) throw std::runtime_error("cannot determine cache file size");
        input.seekg(0, std::ios::beg);
        if (!input) throw std::runtime_error("cannot seek cache file");
        MeshCacheIO::ReadBudget budget(input, static_cast<std::uint64_t>(end));
        QuadtreeCacheDetail::Snapshot staged;
        staged.read(input, budget);
        staged.validate(*this);

        // Prepare all entity owners before changing any receiving storage.
        auto stagedVerts=decltype(verts_)::adoptDense(staged.vertices,staged.capacity[0],static_cast<std::size_t>(staged.levelSets));
        auto stagedHEdges=decltype(hEdges_)::adoptDense(staged.hEdges,staged.capacity[1],static_cast<std::size_t>(staged.levelSets));
        auto stagedVEdges=decltype(vEdges_)::adoptDense(staged.vEdges,staged.capacity[2],static_cast<std::size_t>(staged.levelSets));
        auto stagedCells=decltype(cells_)::adoptDense(staged.cells,staged.capacity[3],static_cast<std::size_t>(staged.levelSets));
        // Everything that can allocate or reject has completed. Domain/Nx/Ny
        // are const and have been checked against all recovered root geometry.
        // Transfer ownership without copying/recomputing geometry or cut data.
        auto commitArray = [](auto*& destination, auto& source) noexcept
        {
            auto* previous = destination;
            destination = source.release();
            source.reset(previous);
        };
        verts_.swap(stagedVerts);
        hEdges_.swap(stagedHEdges);
        vEdges_.swap(stagedVEdges);
        cells_.swap(stagedCells);
        commitArray(vertIndexToLabel, staged.indexToLabel[0]);
        commitArray(hEdgeIndexToLabel, staged.indexToLabel[1]);
        commitArray(vEdgeIndexToLabel, staged.indexToLabel[2]);
        commitArray(cellIndexToLabel, staged.indexToLabel[3]);
        commitArray(vertLabelSpace_, staged.labelSpace[0]);
        commitArray(hEdgeLabelSpace_, staged.labelSpace[1]);
        commitArray(vEdgeLabelSpace_, staged.labelSpace[2]);
        commitArray(cellLabelSpace_, staged.labelSpace[3]);
        levelSetCount_ = staged.levelSets;
        maxLevel_ = staged.maxLevel;
        vertIndexNum_ = staged.count[0]; hEdgeIndexNum_ = staged.count[1];
        vEdgeIndexNum_ = staged.count[2]; cellIndexNum_ = staged.count[3];
        maxVertNum_ = vertLabelCapacity_ = staged.capacity[0];
        maxHEdgeNum_ = hEdgeLabelCapacity_ = staged.capacity[1];
        maxVEdgeNum_ = vEdgeLabelCapacity_ = staged.capacity[2];
        maxCellNum_ = cellLabelCapacity_ = staged.capacity[3];
        vertPtr_ = staged.ptr[0]; hEdgePtr_ = staged.ptr[1];
        vEdgePtr_ = staged.ptr[2]; cellPtr_ = staged.ptr[3];
        minLevel_neg.swap(staged.minNegative);
        minLevel_pos.swap(staged.minPositive);
        maxLevel_local.swap(staged.maxLocal);
        nCellsPerLevel.swap(staged.cellsPerLevel);
        // maxInitialIterations_, C_Lipschitz and other nonserialized settings
        // deliberately retain the receiving mesh's existing values.
    }
    catch (const std::exception& error)
    {
        throw std::runtime_error("Quadtree mesh cache '" + name + "' rejected: " + error.what() +
            ". Rebuild this cache from geometry with the current executable; no cache file was changed.");
    }
}
