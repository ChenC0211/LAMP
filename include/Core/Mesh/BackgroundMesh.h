//                       网格编号信息
// ---------------------------------------------------------------------------
//      *-----*-----*             3
//      |  3  |  2  |           *---*
//      *-----*-----*          0|   |2
//      |  0  |  1  |           *---*
//      *-----*-----*             1
//    父网格与子网格关系    单元与相邻单元的关系
// ---------------------------------------------------------------------------
//      3-----------2        *-----1-----*
//      |           |        |           |
//      |           |        0           1
//      |           |        |           |
//      0-----------1        *-----0-----*
//      单元到点的关系         单元到边的关系
// ---------------------------------------------------------------------------
//            3                    *
//            |                    1
//      0-----*-----2        *--0--*--1--*
//            |                    0
//            1                    *
//     点与相邻点的关系         点到边的关系
// ---------------------------------------------------------------------------

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>
#include "Core/Base/Vector2D.h"
#include "Core/Base/Base.h"
#include "Core/Mesh/Vert.h"
#include "Core/Mesh/Edge.h"
#include "Core/Mesh/Cell.h"
#include "Core/IO/DictionaryMesh.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <memory>
#include <mutex>
#include <limits>
#include <new>
#include <stdexcept>
#include <cstdint>
#include <cstddef>
#include <cfenv>
#include <cerrno>
#include <exception>
#if defined(_M_X64) || defined(__SSE2__)
#include <xmmintrin.h>
#endif

// Label-addressable unique storage. Dense paths retain contiguous arrays;
// LS-free clones leave never-exposed default slots unconstructed.
namespace LAMP { namespace Mesh { namespace Detail {
class EntityPoolEnvironment
{
    std::fenv_t native_{};
    int error_=errno;
    bool captured_=std::fegetenv(&native_)==0;
#if defined(_M_X64) || defined(__SSE2__)
    unsigned csr_=_mm_getcsr();
#endif
public:
    bool captured() const noexcept {return captured_;}
    ~EntityPoolEnvironment() noexcept
    {
        if(captured_&&std::fesetenv(&native_)!=0)std::terminate();
#if defined(_M_X64) || defined(__SSE2__)
        _mm_setcsr(csr_);
#endif
        errno=error_;
    }
};
template<class T> class EntityPool
{
    static constexpr std::size_t chunkSize=256;
    struct Chunk
    {
        std::atomic<T*> published{nullptr};
        std::unique_ptr<T[]> objects;
    };
    struct Range
    {
        std::uintptr_t begin=0,end=0;
        int first=0,count=0;
    };
    struct Storage
    {
        int capacity=0;
        std::size_t levelSets=0,chunkCount=0,rangeCount=0;
        T defaultRecord;
        std::unique_ptr<T[]> dense;
        std::unique_ptr<Chunk[]> chunks;
        std::unique_ptr<Range[]> ranges;
        std::atomic<std::size_t> constructed{0},publishedChunks{0};
        mutable std::mutex mutex;
        Storage(int n,std::size_t ls,bool lazy):capacity(n),levelSets(ls)
        {
            if(n<0||static_cast<std::size_t>(n)>std::numeric_limits<std::size_t>::max()/sizeof(T))
                throw std::bad_array_new_length();
            if(lazy){
                defaultRecord.ensureLevelSetSize(ls);
                chunkCount=(static_cast<std::size_t>(n)+chunkSize-1)/chunkSize;
                if(chunkCount>std::numeric_limits<std::size_t>::max()/(sizeof(Chunk)+sizeof(Range)))
                    throw std::bad_array_new_length();
                chunks=std::make_unique<Chunk[]>(chunkCount);
                ranges=std::make_unique<Range[]>(chunkCount);
            }
        }
        T* materialized(int label) const noexcept
        {
            if(dense)return dense.get()+label;
            const auto index=static_cast<std::size_t>(label);
            T* p=chunks[index/chunkSize].published.load(std::memory_order_acquire);
            return p?p+index%chunkSize:nullptr;
        }
        T& get(int label)
        {
            if(T* p=materialized(label))return *p;
            EntityPoolEnvironment environment;
            if(!environment.captured())throw std::runtime_error("entity pool cannot capture floating-point environment");
            std::lock_guard<std::mutex> lock(mutex);
            if(T* p=materialized(label))return *p;
            const auto index=static_cast<std::size_t>(label),block=index/chunkSize;
            const auto first=block*chunkSize;
            const auto count=std::min(chunkSize,static_cast<std::size_t>(capacity)-first);
            auto objects=std::make_unique<T[]>(count);
            for(std::size_t i=0;i<count;++i)objects[i].ensureLevelSetSize(levelSets);
            const auto begin=reinterpret_cast<std::uintptr_t>(objects.get());
            const auto bytes=count*sizeof(T);
            if(bytes>std::numeric_limits<std::uintptr_t>::max()-begin)throw std::bad_array_new_length();
            const auto end=begin+bytes;
            std::size_t pos=0,last=rangeCount;
            while(pos<last){const auto mid=pos+(last-pos)/2;
                if(ranges[mid].begin<begin)pos=mid+1;else last=mid;}
            for(auto k=rangeCount;k>pos;--k)ranges[k]=ranges[k-1];
            ranges[pos]={begin,end,static_cast<int>(first),static_cast<int>(count)};
            ++rangeCount;
            T* result=objects.get();
            chunks[block].objects=std::move(objects);
            constructed.fetch_add(count,std::memory_order_relaxed);
            publishedChunks.fetch_add(1,std::memory_order_relaxed);
            chunks[block].published.store(result,std::memory_order_release);
            return result[index%chunkSize];
        }
        int labelOf(const void* object) const noexcept
        {
            struct RestoreErrno {int value=errno;~RestoreErrno() noexcept {errno=value;}} restoreErrno;
            const auto p=reinterpret_cast<std::uintptr_t>(object);
            if(dense){const auto base=reinterpret_cast<std::uintptr_t>(dense.get());
                if(p<base)return -1;const auto delta=p-base;
                if(delta%sizeof(T)||delta/sizeof(T)>=static_cast<std::size_t>(capacity))return -1;
                return static_cast<int>(delta/sizeof(T));}
            try{
                std::lock_guard<std::mutex> lock(mutex);
                std::size_t first=0,last=rangeCount;
                while(first<last){const auto mid=first+(last-first)/2;
                    if(ranges[mid].begin<=p)first=mid+1;else last=mid;}
                if(first==0)return -1;const auto& range=ranges[first-1];
                if(p>=range.end||(p-range.begin)%sizeof(T))return -1;
                const auto offset=(p-range.begin)/sizeof(T);
                if(offset>=static_cast<std::size_t>(range.count))return -1;
                return range.first+static_cast<int>(offset);
            }catch(...){return -1;}
        }
    };
    std::unique_ptr<Storage> storage_;
public:
    using value_type=T;
    EntityPool() noexcept=default;
    ~EntityPool() noexcept {clear();}
    EntityPool(const EntityPool&)=delete;
    EntityPool& operator=(const EntityPool&)=delete;
    EntityPool(EntityPool&&) noexcept=default;
    EntityPool& operator=(EntityPool&&) noexcept=default;
    void clear() noexcept {EntityPoolEnvironment environment;storage_.reset();}
    void swap(EntityPool& other) noexcept {storage_.swap(other.storage_);}
    explicit operator bool() const noexcept {return static_cast<bool>(storage_);}
    const void* identity() const noexcept {return storage_.get();}
    int capacity() const noexcept {return storage_?storage_->capacity:0;}
    bool denseMode() const noexcept {return storage_&&static_cast<bool>(storage_->dense);}
    T& operator[](int label) const {return storage_->get(label);}
    const T& peek(int label) const noexcept
    {
        const T* p=storage_->materialized(label);
        return p?*p:storage_->defaultRecord;
    }
    T* materialized(int label) const noexcept {return storage_->materialized(label);}
    std::size_t materializedCount() const noexcept
    {return storage_?storage_->constructed.load(std::memory_order_relaxed):0;}
    std::size_t publishedChunks() const noexcept
    {return storage_?storage_->publishedChunks.load(std::memory_order_relaxed):0;}
    std::array<std::size_t,3> storageStatistics() const noexcept
    {
        if(!storage_)return {};
        EntityPoolEnvironment environment;
        std::lock_guard<std::mutex> lock(storage_->mutex);
        return {static_cast<std::size_t>(storage_->capacity),
            storage_->constructed.load(std::memory_order_relaxed),
            storage_->publishedChunks.load(std::memory_order_relaxed)};
    }
    // Actual fixed address-table capacity, not just the currently published prefix.
    std::size_t addressMetadataBytes() const noexcept
    {return storage_?storage_->chunkCount*(sizeof(Chunk)+sizeof(Range)):0;}
    int labelOf(const void* object) const noexcept {return storage_?storage_->labelOf(object):-1;}
    void resetDense(int n)
    {
        auto next=std::make_unique<Storage>(n,0,false);
        next->dense=std::make_unique<T[]>(static_cast<std::size_t>(n));
        next->constructed.store(static_cast<std::size_t>(n),std::memory_order_relaxed);
        storage_.swap(next);
    }
    void resetLazy(int n,std::size_t levelSets)
    {
        EntityPoolEnvironment environment;
        if(!environment.captured()){
            resetDense(n);
            for(int i=0;i<n;++i)(*this)[i].ensureLevelSetSize(levelSets);
            return;
        }
        auto next=std::make_unique<Storage>(n,levelSets,true);
        storage_.swap(next);
    }
    // All four preparations precede the cache receiver's no-throw commit.
    static EntityPool adoptDense(std::unique_ptr<T[]>& objects,int n,std::size_t levelSets)
    {
        EntityPool result;
        auto next=std::make_unique<Storage>(n,levelSets,false);
        next->dense=std::move(objects);
        next->constructed.store(static_cast<std::size_t>(n),std::memory_order_relaxed);
        result.storage_=std::move(next);return result;
    }
};
}}}

class BackgroundMesh
{
public:
    const Vector2D minCoord;
    const Vector2D maxCoord;
    const int Nx;
    const int Ny;
    const Vector2D cellSize;

    int* vertIndexToLabel = nullptr;
    int* hEdgeIndexToLabel = nullptr;
    int* vEdgeIndexToLabel = nullptr;
    int* cellIndexToLabel = nullptr;

    std::vector<std::vector<int>> negFaceVerts;
    std::vector<std::vector<int>> posFaceVerts;

    int vertIndexNum() const { return vertIndexNum_; }
    int hEdgeIndexNum() const { return hEdgeIndexNum_; }
    int vEdgeIndexNum() const { return vEdgeIndexNum_; }
    int cellIndexNum() const { return cellIndexNum_; }

    int vertLabelCapacity() const { return vertLabelCapacity_; }
    int hEdgeLabelCapacity() const { return hEdgeLabelCapacity_; }
    int vEdgeLabelCapacity() const { return vEdgeLabelCapacity_; }

    // Logical capacity of label-addressable storage (cells_/verts_/edges).
    // For uniform meshes this equals cellIndexNum(); for quadtree meshes it can
    // be larger after internal buffer growth.
    int cellLabelCapacity() const { return cellLabelCapacity_; }

    Vert& verts(int label) const { return verts_[label]; }
    Edge& hEdges(int label) const { return hEdges_[label]; }
    Edge& vEdges(int label) const { return vEdges_[label]; }
    Cell& cells(int label) const { return cells_[label]; }

    int levelSetCount() const { return levelSetCount_; }

    BackgroundMesh(Vector2D _minCoord, Vector2D _maxCoord, int _Nx, int _Ny);
    // Explicit dictMesh.levelSetCount == 0 constructs a mesh without LS storage.
    // The default count and setLevelSetCount retain their legacy minimum of one.
    BackgroundMesh(const DictionaryMesh& dictMesh);
    ~BackgroundMesh();

    int findCellLabel(const Vector2D& coord) const;
    // Explicit physical boundary tolerance; zero performs a strict domain
    // query without the legacy absolute MINVAL length. Descent is shared.
    int findCellLabel(const Vector2D& coord, double boundaryTolerance) const;
    int findBackgroundCellLabel(const Vector2D& coord) const;

    // Fast in-cell check and hint-based cell lookup (stay-in-cell + neighbor ring + fallback).
    // eps is a small tolerance for boundary inclusion.
    bool isInsideCell(int cellLabel, const Vector2D& coord, double eps = MINVAL) const;
    int findCellLabelFromHint(int hintCellLabel, const Vector2D& coord, double eps = MINVAL) const;

    void output(const std::string& name) const;
    void writeCache(const std::string& filePath) const;
    void readCache();

protected:
    // For an owning derived clone only: copy lattice metadata, leave storage
    // empty. The clone fills native records directly, without buildMesh().
    struct CopyMetadataTag {};
    BackgroundMesh(const BackgroundMesh& source, CopyMetadataTag);

    void setLevelSetCount(int n) { levelSetCount_ = (n > 0) ? n : 1; }
    void ensureLevelSetStorage();

protected:
    LAMP::Mesh::Detail::EntityPool<Vert> verts_;
    LAMP::Mesh::Detail::EntityPool<Edge> hEdges_;
    LAMP::Mesh::Detail::EntityPool<Edge> vEdges_;
    LAMP::Mesh::Detail::EntityPool<Cell> cells_;

    int vertIndexNum_ = 0;
    int hEdgeIndexNum_ = 0;
    int vEdgeIndexNum_ = 0;
    int cellIndexNum_ = 0;

    int vertLabelCapacity_ = 0;
    int hEdgeLabelCapacity_ = 0;
    int vEdgeLabelCapacity_ = 0;

    // Logical capacity of cells_ in label space.
    int cellLabelCapacity_ = 0;

    int levelSetCount_ = 1;
    
    void buildMesh();

    void createVertices();
    void createEdges();
    void createCells();

    void findNeighbourCells();

    void defineTypes();

    template<typename T>
    void statistics(LAMP::Mesh::Detail::EntityPool<T>& t, int num);

};
