#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <cstdint>
#include <vector>
#include <utility>

class QuadtreeMesh;

// Exact native geometry/connectivity identity extracted from Hydro's fitted
// domain. No physical boundary/material data or mutable interpolation caches.
// Capture once per mesh state; validate before reusing fields or operators.
class NativeMeshIdentity
{
public:
    // Audited, lexical read-only mesh segment only. First validation remains
    // exact; no proof survives an opaque callback, exception or mesh mutation.
    class ValidationBatch
    {
    public:
        explicit ValidationBatch(const QuadtreeMesh& mesh) noexcept;
        ~ValidationBatch() noexcept;
        ValidationBatch(const ValidationBatch&)=delete;
        ValidationBatch& operator=(const ValidationBatch&)=delete;
        class Pause
        {
        public:
            Pause() noexcept;
            ~Pause() noexcept;
            Pause(const Pause&)=delete;
            Pause& operator=(const Pause&)=delete;
        private:
            ValidationBatch* suspended_;
        };
    private:
        friend class NativeMeshIdentity;
        static thread_local ValidationBatch* active_;
        const QuadtreeMesh* mesh_;
        ValidationBatch* previous_;
        int uncaught_;
        std::array<std::shared_ptr<const unsigned char>,16> proofs_{};
        std::size_t size_=0;
        void invalidateChain() noexcept;
        bool contains(const std::shared_ptr<const unsigned char>& token) const noexcept;
        void remember(const std::shared_ptr<const unsigned char>& token) noexcept;
    };
    explicit NativeMeshIdentity(const QuadtreeMesh& mesh);
    NativeMeshIdentity(const NativeMeshIdentity&)=default;
    NativeMeshIdentity& operator=(const NativeMeshIdentity&)=default;
    NativeMeshIdentity(NativeMeshIdentity&&) noexcept=default;
    NativeMeshIdentity& operator=(NativeMeshIdentity&& other) noexcept
    {
        if(this!=&other){words_=std::move(other.words_);token_=std::move(other.token_);}
        return *this;
    }
    void validate(const QuadtreeMesh& mesh) const;
    // Exact read-only persistence fingerprint; restoring an identity still
    // requires constructing it from an actual mesh and comparing these words.
    const std::vector<std::uint64_t>& checkpointWords() const noexcept{return words_;}
private:
    std::vector<std::uint64_t> words_;
    // Separate lifetime object prevents address reuse after an identity dies.
    std::shared_ptr<const unsigned char> token_;
};

