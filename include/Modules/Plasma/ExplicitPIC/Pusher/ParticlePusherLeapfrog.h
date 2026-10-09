#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Core/Base/Vector2D.h"
#include "Modules/Plasma/Particles/ParticleType.h"
#include "Modules/Plasma/Particles/Particle.h"
#include "Core/Mesh/Cell.h"
#include "Core/Mesh/Vert.h"
#include "Core/Geometry/ImmersedGeometry.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ExplicitPIC
{
namespace Pusher
{

// Domain side classification for particles absorbed at the domain boundary.
// The domain is the mesh bounding box [minCoord, maxCoord].
enum class DomainSide
{
    kLeft,
    kRight,
    kBottom,
    kTop,
    kUnknown,
};

struct DomainImpact
{
    Particle particle;
    DomainSide side = DomainSide::kUnknown;
};

// A lightweight buffer to record absorbed particles and where they exited.
// This enables wall-interaction models (e.g. SEE) to be implemented at the module level.
class DomainImpactBuffer
{
public:
    void Clear() { impacts_.clear(); }
    bool empty() const { return impacts_.empty(); }
    std::size_t size() const { return impacts_.size(); }

    void Push(const Particle& p, const DomainSide side)
    {
        impacts_.push_back(DomainImpact{p, side});
    }

    const std::vector<DomainImpact>& impacts() const { return impacts_; }

private:
    std::vector<DomainImpact> impacts_;
};

struct EmbeddedImpact
{
    EmbeddedImpact() = default;
    EmbeddedImpact(const Particle& p, const Vector2D& position, int cell, bool didCross)
        : particle(p), pos(position), cellLabel(cell), crossing(didCross) {}

    Particle particle;
    Vector2D pos; // representative position on/near the EB
    int cellLabel = -1;
    int boundaryId = 0;
    int levelSetIndex = -1;
    int lineIndex = -1;
    Vector2D normal;
    bool crossing = false;
};

// Buffer to record EB-absorbed particles (inside or crossing).
class EmbeddedImpactBuffer
{
public:
    void Clear() { impacts_.clear(); }
    bool empty() const { return impacts_.empty(); }
    std::size_t size() const { return impacts_.size(); }

    void Push(const Particle& p, const Vector2D& pos, int cellLabel, bool crossing)
    {
        EmbeddedImpact imp;
        imp.particle = p;
        imp.pos = pos;
        imp.cellLabel = cellLabel;
        imp.crossing = crossing;
        impacts_.push_back(std::move(imp));
    }

    const std::vector<EmbeddedImpact>& impacts() const { return impacts_; }

private:
    std::vector<EmbeddedImpact> impacts_;
};

// ============================================================================
// ParticlePusherLeapfrog
// ----------------------------------------------------------------------------
// Kick-drift leapfrog pusher:
//   v^{n+1/2} = v^{n-1/2} + kappa_s E(x^n) dt
//   x^{n+1}   = x^n + v^{n+1/2} dt
//
// After advancing, the particle's host cell label is updated via
//   mesh.findCellLabel(position).
//
// Absorption (deletion) policy:
//   - Domain: if x is outside [minCoord, maxCoord] -> delete
//   - EB: if any queried level set psi(x) < 0 -> delete
//     (negative region is treated as "inside")
//
// Notes:
//   * This pusher operates on ParticleType::particle (vector<Particle>) and
//     compacts survivors in-place.
//   * The input Ep must correspond 1-to-1 with the particle ordering *before*
//     deletion.
// ============================================================================

class ParticlePusherLeapfrog
{
public:
    struct Options
    {
        bool absorbAtDomain = true;
        bool absorbAtEmbedded = true;
        bool absorbEmbeddedCrossing = false;

        // Periodic wrap for domain coordinates. When enabled, particles that
        // drift beyond the mesh bounding box in the corresponding direction
        // are wrapped back into the domain before any domain absorption logic.
        //
        // Convention: coordinates are wrapped into [minCoord, maxCoord).
        // This half-open interval avoids corner cases where positions land
        // exactly on maxCoord and later fail cell lookup.
        bool periodicX = false;
        bool periodicY = false;

        // If non-empty, only these level-set indices are used for EB absorption.
        // Otherwise, all iGeometry.levelSets are checked.
        std::vector<int> embeddedLevelSetIds;

        // Level-set sign convention for EB absorption:
        //   true  -> inside (solid) is psi < 0
        //   false -> inside (solid) is psi > 0
        bool embeddedInsideIsNegative = true;

        // Treat the domain box as inclusive (<=, >=). Particles strictly outside
        // are absorbed.
        double eps = 0.0;
    };

    // Workaround for GCC: avoid Options{} as a default argument inside the enclosing class.
    static Options defaultOptions() { return Options{}; }

    // Initialize leapfrog by applying a half-kick:
    //   v^{1/2} = v^{0} + kappa_s E(x^0) (dt/2)
    // This does NOT drift the position.
    template<typename Mesh>
    static void initializeHalfStep(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dt,
        const Options& opt,
        ImmersedGeometry<Mesh>* iGeometry)
    {
        kickImpl(mesh, species, Ep, 0.5 * dt, opt, iGeometry);
    }

    // Convenience overloads without Options/impact buffers
    template<typename Mesh>
    static void initializeHalfStep(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dt)
    {
        const Options opt = defaultOptions();
        kickImpl(mesh, species, Ep, 0.5 * dt, opt, /*iGeometry=*/nullptr);
    }

    template<typename Mesh>
    static void initializeHalfStep(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dt,
        ImmersedGeometry<Mesh>* iGeometry)
    {
        const Options opt = defaultOptions();
        kickImpl(mesh, species, Ep, 0.5 * dt, opt, iGeometry);
    }

    // Advance one full step (kick + drift) and absorb particles.
    template<typename Mesh>
    static void advance(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dt,
        const Options& opt,
        ImmersedGeometry<Mesh>* iGeometry,
        DomainImpactBuffer* impactBuffer,
        EmbeddedImpactBuffer* embeddedBuffer = nullptr)
    {
        if (dt <= 0.0)
            return;

        const int n = static_cast<int>(species.particle.size());
        if (static_cast<int>(Ep.size()) != n)
            throw std::runtime_error("[ParticlePusherLeapfrog] Ep size mismatch with particle list.");

        // aliveMask + in-place stable compaction.
        // Use a reusable scratch buffer to avoid per-step allocations.
        Workspace& ws = workspace_();
        ws.alive.resize(static_cast<size_t>(n));
        std::fill(ws.alive.begin(), ws.alive.end(), static_cast<std::uint8_t>(0));

        // Thread-local impact buffers (merged after the parallel region).
        if (impactBuffer != nullptr)
        {
            #ifdef _OPENMP
            ws.impactsTLS.resize(static_cast<size_t>(omp_get_max_threads()));
            for (auto& v : ws.impactsTLS)
                v.clear();
            #else
            ws.impactsSerial.clear();
            #endif
        }
        if (embeddedBuffer != nullptr)
        {
            #ifdef _OPENMP
            ws.embeddedImpactsTLS.resize(static_cast<size_t>(omp_get_max_threads()));
            for (auto& v : ws.embeddedImpactsTLS)
                v.clear();
            #else
            ws.embeddedImpactsSerial.clear();
            #endif
        }

        const double kappa = species.kappa_s();
        #ifdef _OPENMP
        const bool useOmp = (n >= 20000);
        #pragma omp parallel for if(useOmp) schedule(static)
        #endif
        for (int p = 0; p < n; ++p)
        {
            Particle particle = species.particle[p];
            const Vector2D pos0 = particle.position;
            const int cell0 = particle.inCellL;

            // ---- Kick: v^{n+1/2} ----
            particle.velocity.x += kappa * Ep[static_cast<size_t>(p)].x * dt;
            particle.velocity.y += kappa * Ep[static_cast<size_t>(p)].y * dt;

            // ---- Drift: x^{n+1} ----
            particle.position.x += particle.velocity.x * dt;
            particle.position.y += particle.velocity.y * dt;

            // ---- Periodic wrap (optional) ----
            if (opt.periodicX)
                wrapPeriodicCoord(particle.position.x, mesh.minCoord.x, mesh.maxCoord.x, opt.eps);
            if (opt.periodicY)
                wrapPeriodicCoord(particle.position.y, mesh.minCoord.y, mesh.maxCoord.y, opt.eps);

            // ---- Domain absorption ----
            if (opt.absorbAtDomain && isOutsideDomain(mesh, particle.position, opt.eps))
            {
                if (impactBuffer != nullptr)
                {
                    const DomainSide side = classifyDomainSide(mesh, particle.position, opt.eps);
                    #ifdef _OPENMP
                    const int tid = omp_get_thread_num();
                    ws.impactsTLS[static_cast<size_t>(tid)].push_back(DomainImpact{particle, side});
                    #else
                    ws.impactsSerial.push_back(DomainImpact{particle, side});
                    #endif
                }
                continue;
            }

            // ---- Update host cell label ----
            particle.inCellL = mesh.findCellLabelFromHint(particle.inCellL, particle.position);
            if (particle.inCellL < 0)
            {
                if (impactBuffer != nullptr)
                {
                    const DomainSide side = classifyDomainSide(mesh, particle.position, opt.eps);
                    #ifdef _OPENMP
                    const int tid = omp_get_thread_num();
                    ws.impactsTLS[static_cast<size_t>(tid)].push_back(DomainImpact{particle, side});
                    #else
                    ws.impactsSerial.push_back(DomainImpact{particle, side});
                    #endif
                }
                continue;
            }

            // ---- EB absorption (optional) ----
            if (opt.absorbAtEmbedded && iGeometry != nullptr)
            {
                if (isInsideEmbedded(mesh, *iGeometry, particle, opt))
                {
                    if (embeddedBuffer != nullptr)
                    {
                        #ifdef _OPENMP
                        const int tid = omp_get_thread_num();
                        ws.embeddedImpactsTLS[static_cast<size_t>(tid)].push_back(
                            EmbeddedImpact{particle, particle.position, particle.inCellL, /*crossing=*/false});
                        #else
                        ws.embeddedImpactsSerial.push_back(
                            EmbeddedImpact{particle, particle.position, particle.inCellL, /*crossing=*/false});
                        #endif
                    }
                    continue;
                }
                if (opt.absorbEmbeddedCrossing)
                {
                    int cellStart = cell0;
                    if (cellStart < 0)
                        cellStart = mesh.findCellLabel(pos0);
                    if (cellStart >= 0)
                    {
                        Vector2D hit = (pos0 + particle.position) * 0.5;
                        int hitCell = cellStart;
                        if (crossesEmbedded(mesh, *iGeometry, pos0, cellStart,
                                            particle.position, particle.inCellL, opt, &hit, &hitCell))
                        {
                            if (embeddedBuffer != nullptr)
                            {
                                #ifdef _OPENMP
                                const int tid = omp_get_thread_num();
                                ws.embeddedImpactsTLS[static_cast<size_t>(tid)].push_back(
                                    EmbeddedImpact{particle, hit, hitCell, /*crossing=*/true});
                                #else
                                ws.embeddedImpactsSerial.push_back(
                                    EmbeddedImpact{particle, hit, hitCell, /*crossing=*/true});
                                #endif
                            }
                            continue;
                        }
                    }
                }
            }

            // ---- Background cell label (best-effort) ----
            particle.inBackgroundCellL = mesh.findBackgroundCellLabel(particle.position);

            ws.alive[static_cast<size_t>(p)] = 1;
            // Each thread writes a distinct index; this is thread-safe.
            species.particle[static_cast<size_t>(p)] = particle;
        }

        // Merge impacts into the user buffer in deterministic thread-order.
        if (impactBuffer != nullptr)
        {
            #ifdef _OPENMP
            for (const auto& local : ws.impactsTLS)
            {
                for (const auto& imp : local)
                    impactBuffer->Push(imp.particle, imp.side);
            }
            #else
            for (const auto& imp : ws.impactsSerial)
                impactBuffer->Push(imp.particle, imp.side);
            #endif
        }
        if (embeddedBuffer != nullptr)
        {
            #ifdef _OPENMP
            for (const auto& local : ws.embeddedImpactsTLS)
            {
                for (const auto& imp : local)
                    embeddedBuffer->Push(imp.particle, imp.pos, imp.cellLabel, imp.crossing);
            }
            #else
            for (const auto& imp : ws.embeddedImpactsSerial)
                embeddedBuffer->Push(imp.particle, imp.pos, imp.cellLabel, imp.crossing);
            #endif
        }

        // Stable in-place compaction (write-pointer). No prefix-sum, no extra vectors.
        int w = 0;
        for (int i = 0; i < n; ++i)
        {
            if (ws.alive[static_cast<size_t>(i)] == 0) continue;
            if (w != i)
                species.particle[static_cast<size_t>(w)] = species.particle[static_cast<size_t>(i)];
            ++w;
        }
        species.particle.resize(static_cast<size_t>(w));
        species.number = static_cast<int>(w);
    }

    template<typename Mesh>
    static void advance(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dt)
    {
        const Options opt = defaultOptions();
        advance(mesh, species, Ep, dt, opt, /*iGeometry=*/nullptr, /*impactBuffer=*/nullptr);
    }

    template<typename Mesh>
    static void advance(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dt,
        ImmersedGeometry<Mesh>* iGeometry)
    {
        const Options opt = defaultOptions();
        advance(mesh, species, Ep, dt, opt, iGeometry, /*impactBuffer=*/nullptr);
    }

    template<typename Mesh>
    static void advance(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dt,
        DomainImpactBuffer* impactBuffer)
    {
        const Options opt = defaultOptions();
        advance(mesh, species, Ep, dt, opt, /*iGeometry=*/nullptr, impactBuffer);
    }

private:
    // Reusable scratch buffers to avoid per-step allocations in the hot loop.
    // This is per-thread (C++ thread_local) to remain safe if multiple modules
    // ever advance in parallel.
    struct Workspace
    {
        std::vector<std::uint8_t> alive;

        #ifdef _OPENMP
        std::vector<std::vector<DomainImpact>> impactsTLS;
        std::vector<std::vector<EmbeddedImpact>> embeddedImpactsTLS;
        #else
        std::vector<DomainImpact> impactsSerial;
        std::vector<EmbeddedImpact> embeddedImpactsSerial;
        #endif
    };

    static Workspace& workspace_()
    {
        static thread_local Workspace ws;
        return ws;
    }

    // Wrap coordinate into [xmin, xmax). Robust for large excursions.
    static void wrapPeriodicCoord(double& x, const double xmin, const double xmax, const double eps)
    {
        const double L = xmax - xmin;
        if (!(L > 0.0) || !std::isfinite(L) || !std::isfinite(x))
            return;

        // Fast path: already inside the half-open interval.
        if (x >= xmin && x < xmax)
            return;

        // Map to [0, L) then shift by xmin.
        double t = std::fmod(x - xmin, L);
        if (t < 0.0)
            t += L;
        x = xmin + t;

        // Guard against round-off that could place x at xmax.
        const double tol = (eps > 0.0) ? eps : 0.0;
        if (x >= xmax - tol)
            x = xmin;
    }
    template<typename Mesh>
    static bool isOutsideDomain(const Mesh& mesh, const Vector2D& x, const double eps)
    {
        return (x.x < mesh.minCoord.x - eps) || (x.x > mesh.maxCoord.x + eps) ||
               (x.y < mesh.minCoord.y - eps) || (x.y > mesh.maxCoord.y + eps);
    }

    template<typename Mesh>
    static DomainSide classifyDomainSide(const Mesh& mesh, const Vector2D& x, const double eps)
    {
        if (x.x < mesh.minCoord.x - eps) return DomainSide::kLeft;
        if (x.x > mesh.maxCoord.x + eps) return DomainSide::kRight;
        if (x.y < mesh.minCoord.y - eps) return DomainSide::kBottom;
        if (x.y > mesh.maxCoord.y + eps) return DomainSide::kTop;
        return DomainSide::kUnknown;
    }

    // Evaluate psi(x) in the particle's current cell using bilinear weights.
    template<typename Mesh>
    static double evalPsiInCell(
        Mesh& mesh,
        ImmersedGeometry<Mesh>& iGeometry,
        const int levelSetI,
        Particle& particle)
    {
        Cell& cell = mesh.cells(particle.inCellL);
        particle.updateWeight(mesh, cell);

        double psi = 0.0;
        for (int i = 0; i < 4; ++i)
        {
            Vert& v = mesh.verts(cell.cellToVert[i]);
            psi += iGeometry.levelSets[levelSetI].value[v.index()] * particle.weight[i];
        }
        return psi;
    }

    template<typename Mesh>
    static double evalPsiAt(
        Mesh& mesh,
        ImmersedGeometry<Mesh>& iGeometry,
        const int levelSetI,
        const Vector2D& pos,
        const int cellLabel)
    {
        if (cellLabel < 0)
            return std::numeric_limits<double>::quiet_NaN();
        Particle tmp;
        tmp.position = pos;
        tmp.inCellL = cellLabel;
        return evalPsiInCell(mesh, iGeometry, levelSetI, tmp);
    }

    template<typename Mesh>
    static bool posInCell(const Cell& c, const Vector2D& pos)
    {
        const Vector2D cc = c.centerCoord();
        const Vector2D s = c.size();
        const double h = std::min(s.x, s.y);
        const double tol = 1e-12 + 1e-9 * h;
        if (pos.x < cc.x - 0.5 * s.x - tol || pos.x > cc.x + 0.5 * s.x + tol)
            return false;
        if (pos.y < cc.y - 0.5 * s.y - tol || pos.y > cc.y + 0.5 * s.y + tol)
            return false;
        return true;
    }

    template<typename Mesh>
    static int resolveCellLabel(const Mesh& mesh, const Vector2D& pos, const int hint)
    {
        if (hint >= 0)
        {
            const Cell& c = mesh.cells(hint);
            if (posInCell<Mesh>(c, pos))
                return hint;
            int cell = mesh.findCellLabelFromHint(hint, pos);
            if (cell >= 0)
            {
                const Cell& c2 = mesh.cells(cell);
                if (posInCell<Mesh>(c2, pos))
                    return cell;
            }
        }
        int cell = mesh.findCellLabel(pos);
        if (cell >= 0)
        {
            const Cell& c = mesh.cells(cell);
            if (posInCell<Mesh>(c, pos))
                return cell;
        }
        return cell;
    }

    static double psiTolFromCell(const Cell& c)
    {
        const Vector2D s = c.size();
        const double h = std::min(s.x, s.y);
        return std::max(1e-12, 1e-3 * h);
    }

    static int signWithEps(const double psi, const double eps)
    {
        if (psi > eps)
            return 1;
        if (psi < -eps)
            return -1;
        return 0;
    }

    template<typename Mesh>
    static bool crossesEmbedded(
        Mesh& mesh,
        ImmersedGeometry<Mesh>& iGeometry,
        const Vector2D& x0,
        const int cell0,
        const Vector2D& x1,
        const int cell1,
        const Options& opt,
        Vector2D* hitPos,
        int* hitCellLabel)
    {
        const int nLS = static_cast<int>(iGeometry.levelSets.size());
        if (nLS <= 0)
            return false;
        std::vector<int> lsList;
        if (opt.embeddedLevelSetIds.empty())
        {
            lsList.reserve(static_cast<std::size_t>(nLS));
            for (int i = 0; i < nLS; ++i)
                lsList.push_back(i);
        }
        else
        {
            lsList = opt.embeddedLevelSetIds;
        }

        const Vector2D dx = x1 - x0;
        const double fluidSign = opt.embeddedInsideIsNegative ? 1.0 : -1.0;

        auto samplePsi = [&](int lsI,
                             const Vector2D& pos,
                             int hint,
                             double* psiOut,
                             int* cellOut,
                             double* epsOut,
                             int* signOut) -> bool {
            int cell = resolveCellLabel(mesh, pos, hint);
            if (cell < 0)
                return false;
            const Cell& c = mesh.cells(cell);
            const double psi = evalPsiAt(mesh, iGeometry, lsI, pos, cell);
            if (!std::isfinite(psi))
                return false;
            const double epsPsi = psiTolFromCell(c);
            if (psiOut) *psiOut = psi;
            if (cellOut) *cellOut = cell;
            if (epsOut) *epsOut = epsPsi;
            if (signOut) *signOut = signWithEps(psi, epsPsi);
            return true;
        };

        struct Sample
        {
            double t = 0.0;
            Vector2D pos;
            int cell = -1;
            double psi = 0.0;
            double eps = 0.0;
            int sign = 0;
            bool ok = false;
        };

        const double tSamples[5] = {0.0, 0.25, 0.5, 0.75, 1.0};
        for (const int lsI : lsList)
        {
            if (lsI < 0 || lsI >= nLS)
                continue;

            Sample s[5];
            int hint = cell0;
            for (int i = 0; i < 5; ++i)
            {
                s[i].t = tSamples[i];
                s[i].pos = x0 + dx * s[i].t;
                s[i].ok = samplePsi(lsI, s[i].pos, hint, &s[i].psi, &s[i].cell, &s[i].eps, &s[i].sign);
                if (s[i].ok)
                    hint = s[i].cell;
            }

            for (int i = 0; i < 4; ++i)
            {
                if (!s[i].ok || !s[i + 1].ok)
                    continue;

                const int s0 = s[i].sign;
                const int s1 = s[i + 1].sign;
                if (s0 == 0 || s1 == 0 || (s0 * s1 < 0))
                {
                    double tL = s[i].t;
                    double tR = s[i + 1].t;
                    int signL = s0;
                    int signR = s1;
                    double psiL = s[i].psi;
                    double psiR = s[i + 1].psi;

                    double tHit = (s0 == 0) ? tL : ((s1 == 0) ? tR : 0.5 * (tL + tR));
                    if (s0 != 0 && s1 != 0)
                    {
                        int hintCell = s[i].cell;
                        for (int it = 0; it < 20; ++it)
                        {
                            const double tM = 0.5 * (tL + tR);
                            Vector2D xM = x0 + dx * tM;
                            double psiM = 0.0;
                            double epsM = 0.0;
                            int signM = 0;
                            int cellM = -1;
                            if (!samplePsi(lsI, xM, hintCell, &psiM, &cellM, &epsM, &signM))
                                break;
                            hintCell = cellM;
                            if (signM == 0)
                            {
                                tHit = tM;
                                break;
                            }
                            if (signM == signL)
                            {
                                tL = tM;
                                psiL = psiM;
                                signL = signM;
                            }
                            else
                            {
                                tR = tM;
                                psiR = psiM;
                                signR = signM;
                            }
                            tHit = 0.5 * (tL + tR);
                        }
                    }

                    Vector2D xHit = x0 + dx * tHit;
                    int cellHit = resolveCellLabel(mesh, xHit, s[i].cell);
                    if (cellHit >= 0)
                    {
                        const Cell& cHit = mesh.cells(cellHit);
                        const double h = std::min(cHit.size().x, cHit.size().y);
                        const double epsShift = std::max(1e-12, 1e-6 * h);
                        const double len = std::max(1e-30, dx.mod());
                        Vector2D dir(dx.x / len, dx.y / len);
                        if (signL == static_cast<int>(fluidSign))
                            dir = Vector2D(-dir.x, -dir.y);
                        else if (signR == static_cast<int>(fluidSign))
                            dir = Vector2D(dir.x, dir.y);
                        xHit = xHit + dir * epsShift;
                        cellHit = resolveCellLabel(mesh, xHit, cellHit);
                    }

                    if (hitPos)
                        *hitPos = xHit;
                    if (hitCellLabel)
                        *hitCellLabel = cellHit;
                    return true;
                }
            }
        }
        return false;
    }

    template<typename Mesh>
    static bool isInsideEmbedded(
        Mesh& mesh,
        ImmersedGeometry<Mesh>& iGeometry,
        Particle& particle,
        const Options& opt)
    {
        const int nLS = static_cast<int>(iGeometry.levelSets.size());
        if (nLS <= 0)
            return false;

        if (!opt.embeddedLevelSetIds.empty())
        {
            for (const int lsI : opt.embeddedLevelSetIds)
            {
                if (lsI < 0 || lsI >= nLS) continue;
                const double psi = evalPsiInCell(mesh, iGeometry, lsI, particle);
                const Cell& c = mesh.cells(particle.inCellL);
                const double epsPsi = psiTolFromCell(c);
                const bool inside = opt.embeddedInsideIsNegative ? (psi <= epsPsi) : (psi >= -epsPsi);
                if (inside)
                    return true;
            }
            return false;
        }

        for (int lsI = 0; lsI < nLS; ++lsI)
        {
            const double psi = evalPsiInCell(mesh, iGeometry, lsI, particle);
            const Cell& c = mesh.cells(particle.inCellL);
            const double epsPsi = psiTolFromCell(c);
            const bool inside = opt.embeddedInsideIsNegative ? (psi <= epsPsi) : (psi >= -epsPsi);
            if (inside)
                return true;
        }
        return false;
    }

    template<typename Mesh>
    static void kickImpl(
        Mesh& mesh,
        ParticleType& species,
        const std::vector<Vector2D>& Ep,
        const double dtKick,
        const Options& opt,
        ImmersedGeometry<Mesh>* iGeometry)
    {
        if (dtKick == 0.0)
            return;

        const int n = static_cast<int>(species.particle.size());
        if (static_cast<int>(Ep.size()) != n)
            throw std::runtime_error("[ParticlePusherLeapfrog] Ep size mismatch with particle list.");

        Workspace& ws = workspace_();
        ws.alive.resize(static_cast<size_t>(n));
        std::fill(ws.alive.begin(), ws.alive.end(), static_cast<std::uint8_t>(0));

        const double kappa = species.kappa_s();
        #ifdef _OPENMP
        const bool useOmp = (n >= 20000);
        #pragma omp parallel for if(useOmp) schedule(static)
        #endif
        for (int p = 0; p < n; ++p)
        {
            Particle particle = species.particle[p];

            if (opt.absorbAtDomain && isOutsideDomain(mesh, particle.position, opt.eps))
                continue;

            particle.inCellL = mesh.findCellLabelFromHint(particle.inCellL, particle.position);
            if (particle.inCellL < 0)
                continue;

            if (opt.absorbAtEmbedded && iGeometry != nullptr)
            {
                if (isInsideEmbedded(mesh, *iGeometry, particle, opt))
                    continue;
            }

            particle.velocity.x += kappa * Ep[static_cast<size_t>(p)].x * dtKick;
            particle.velocity.y += kappa * Ep[static_cast<size_t>(p)].y * dtKick;

            particle.inBackgroundCellL = mesh.findBackgroundCellLabel(particle.position);

            ws.alive[static_cast<size_t>(p)] = 1;
            species.particle[static_cast<size_t>(p)] = particle;
        }

        int w = 0;
        for (int i = 0; i < n; ++i)
        {
            if (ws.alive[static_cast<size_t>(i)] == 0) continue;
            if (w != i)
                species.particle[static_cast<size_t>(w)] = species.particle[static_cast<size_t>(i)];
            ++w;
        }
        species.particle.resize(static_cast<size_t>(w));
        species.number = static_cast<int>(w);
    }
};

} // namespace Pusher
} // namespace ExplicitPIC
