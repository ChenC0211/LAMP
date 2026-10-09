#include "Modules/Plasma/Boundary/PicBoundarySetup.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/ConsoleControl.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
std::string lower(std::string s)
{
    for (char& c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string inferSide(const DictionaryGeometry& dictGeometry,
    const DictionaryGeometry::BoundaryEntity& ent)
{
    const double dxDom = dictGeometry.domainMax.x - dictGeometry.domainMin.x;
    const double dyDom = dictGeometry.domainMax.y - dictGeometry.domainMin.y;
    const double tol = 1e-10 * std::max(dxDom, dyDom) + 1e-14;

    auto approxEq = [&](double a, double b) -> bool {
        return std::abs(a - b) <= tol;
    };

    std::string side;

    auto classifyLine = [&](int lineId) -> std::string {
        if (lineId < 0 || lineId >= static_cast<int>(dictGeometry.lines.size()))
            return std::string();

        const Line& L = dictGeometry.lines[static_cast<size_t>(lineId)];
        const double x0 = L.points[0].x;
        const double y0 = L.points[0].y;
        const double x1 = L.points[1].x;
        const double y1 = L.points[1].y;

        const bool onLeft   = approxEq(x0, dictGeometry.domainMin.x) && approxEq(x1, dictGeometry.domainMin.x);
        const bool onRight  = approxEq(x0, dictGeometry.domainMax.x) && approxEq(x1, dictGeometry.domainMax.x);
        const bool onBottom = approxEq(y0, dictGeometry.domainMin.y) && approxEq(y1, dictGeometry.domainMin.y);
        const bool onTop    = approxEq(y0, dictGeometry.domainMax.y) && approxEq(y1, dictGeometry.domainMax.y);

        if (onLeft)   return "left";
        if (onRight)  return "right";
        if (onBottom) return "bottom";
        if (onTop)    return "top";
        return std::string();
    };

    for (int lid : ent.lineIds)
    {
        const std::string s = classifyLine(lid);
        if (s.empty())
            continue;

        if (side.empty())
            side = s;
        else if (side != s)
            return std::string();
    }

    if (!side.empty())
        return side;

    const std::string n = lower(ent.name);
    if (n == "left" || n == "right" || n == "bottom" || n == "top")
        return n;

    return std::string();
}
}

namespace PicBoundarySetup
{
Fields mapFields(const PicConfig& cfg, const DictionaryGeometry& dictGeometry,
    bool haveGeometry, Mode mode)
{
    Fields result;
    result.domain = cfg.domainBC;
    result.embedded = cfg.embeddedBC;
    if (!haveGeometry) return result;
    auto& domainBCMapped = result.domain;
    auto& embeddedBCMapped = result.embedded;
    auto& embeddedBCByIdMapped = result.embeddedById;
    auto& lineToBoundaryId = result.lineToBoundaryId;
    const char* prefix = mode == Mode::Explicit ? "[ExplicitPIC]" : "[ImplicitPIC]";
    // --------------------------------------------------------------
    // Build geometry line -> boundaryId mapping (for EB cut-cell tagging)
    // --------------------------------------------------------------
    lineToBoundaryId.assign(dictGeometry.lines.size(), 0);
    for (const auto& ent : dictGeometry.boundaries)
    {
        if (ent.ignored)
            continue;
        if (ent.kind == GeometricBoundary::domain)
            continue;
        for (int lid : ent.lineIds)
        {
            if (lid < 0 || lid >= static_cast<int>(lineToBoundaryId.size()))
                continue;
            const int prev = lineToBoundaryId[static_cast<size_t>(lid)];
            if (prev != 0 && prev != ent.id)
            {
                if (mode == Mode::Explicit)
                    std::cerr << prefix << "[Geometry] WARNING: geometry lineId " << (lid + 1)
                              << " is mapped to multiple boundaryIds (" << prev << " and " << ent.id
                              << "). Using the first." << std::endl;
                continue;
            }
            lineToBoundaryId[static_cast<size_t>(lid)] = ent.id;
        }
    }

    auto relEq = [&](double x, double y) -> bool {
        const double s = 1.0 + std::max(std::abs(x), std::abs(y));
        return std::abs(x - y) <= 1e-12 * s;
    };

    auto timeFnEqual = [&](const PicBoundarySide::TimeFunction& a,
                           const PicBoundarySide::TimeFunction& b) -> bool {
        if (a.enabled != b.enabled)
            return false;
        if (!a.enabled)
            return true;
        if (lower(a.type) != lower(b.type))
            return false;
        return relEq(a.dc, b.dc) && relEq(a.ac, b.ac) && relEq(a.freq_hz, b.freq_hz) &&
               relEq(a.phase_rad, b.phase_rad) && relEq(a.t0_s, b.t0_s);
    };

    auto sideEqual = [&](const PicBoundarySide& a, const PicBoundarySide& b) -> bool {
        if (a.type != b.type)
            return false;
        if (!relEq(a.value, b.value))
            return false;
        return timeFnEqual(a.timeDirichlet, b.timeDirichlet);
    };

    auto embeddedEqual = [&](const PicEmbeddedBC& a, const PicEmbeddedBC& b) -> bool {
        if (a.enabled != b.enabled)
            return false;
        if (!a.enabled)
            return true;
        if (a.type != b.type)
            return false;
        if (!relEq(a.a, b.a) || !relEq(a.b, b.b) || !relEq(a.g, b.g))
            return false;
        return timeFnEqual(a.timeDirichlet, b.timeDirichlet);
    };

    auto applyDomainSide = [&](const std::string& side, const PicBoundarySide& bc) {
        if (side == "left")
            domainBCMapped.left = bc;
        else if (side == "right")
            domainBCMapped.right = bc;
        else if (side == "bottom")
            domainBCMapped.bottom = bc;
        else if (side == "top")
            domainBCMapped.top = bc;
    };

    auto sideToEmbedded = [&](const PicBoundarySide& bc, int boundaryId) -> PicEmbeddedBC {
        PicEmbeddedBC e;
        e.enabled = true;
        e.timeDirichlet = bc.timeDirichlet;

        if (bc.type == PicBCType::Dirichlet)
        {
            e.type = PicBCType::Dirichlet;
            e.a = 1.0;
            e.b = 0.0;
            e.g = bc.value;
        }
        else if (bc.type == PicBCType::Neumann)
        {
            e.type = PicBCType::Neumann;
            e.a = 0.0;
            e.b = 1.0;
            e.g = bc.value;
        }
        else
        {
            std::ostringstream oss;
            oss << prefix << "[BC] ERROR: boundaryId " << boundaryId
                << " requests Robin BC on an embedded boundary via BC/Field, but BoundarySide has no (a,b). "
                << "Use the global EmbeddedBC block in pic.dict.";
            throw std::runtime_error(oss.str());
        }

        return e;
    };

    // Track per-side assignments (in case geometry splits a side into multiple boundaryIds).
    struct SideAssign
    {
        bool set = false;
        PicBoundarySide bc;
        int id = -1;
    };
    std::unordered_map<std::string, SideAssign> sideAssign;

    auto assignSide = [&](const std::string& side, const PicBoundarySide& bc, int boundaryId) {
        if (side.empty())
        {
            std::ostringstream oss;
            oss << prefix << "[BC] ERROR: cannot infer domain side for domainBoundary id="
                << boundaryId
                << ". Ensure its line segments lie exactly on a domain edge (x=min/max or y=min/max).";
            throw std::runtime_error(oss.str());
        }

        auto& slot = sideAssign[side];
        if (!slot.set)
        {
            slot.set = true;
            slot.bc = bc;
            slot.id = boundaryId;
            return;
        }

        if (!sideEqual(slot.bc, bc))
        {
            std::ostringstream oss;
            oss << prefix << "[BC] ERROR: multiple boundaryId Field BCs mapped to domain side '"
                << side << "' but conflict (ids " << slot.id << " and " << boundaryId << ").";
            throw std::runtime_error(oss.str());
        }
    };

    bool embeddedAssigned = false;
    bool embeddedAllSame = true;
    PicEmbeddedBC embeddedFromIds = embeddedBCMapped;

    // -----------------------------------------------------------------
    // Field BC mapping (boundaryId -> domain side or embedded)
    // -----------------------------------------------------------------
    for (const auto& bc : cfg.fieldBCById)
    {
        auto it = dictGeometry.boundaryById.find(bc.id);
        if (it == dictGeometry.boundaryById.end())
        {
            if (mode == Mode::Explicit)
                std::cerr << prefix << "[BC] WARNING: boundaryId " << bc.id
                          << " not found in geometry.dict; ignoring Field BC.\n";
            continue;
        }
        const auto& ent = it->second;
        if (ent.ignored)
            continue;

        if (ent.kind == GeometricBoundary::domain)
        {
            const std::string side = inferSide(dictGeometry, ent);
            assignSide(side, bc.bc, bc.id);
            applyDomainSide(side, bc.bc);
        }
        else
        {
            const PicEmbeddedBC cand = sideToEmbedded(bc.bc, bc.id);
            auto itMap = embeddedBCByIdMapped.find(bc.id);
            if (itMap != embeddedBCByIdMapped.end() && !embeddedEqual(itMap->second, cand))
            {
                std::ostringstream oss;
                oss << prefix << "[BC] ERROR: boundaryId " << bc.id
                    << " appears multiple times with conflicting embedded BC definitions.";
                throw std::runtime_error(oss.str());
            }
            embeddedBCByIdMapped[bc.id] = cand;
            if (!embeddedAssigned)
            {
                embeddedFromIds = cand;
                embeddedAssigned = true;
            }
            else if (!embeddedEqual(embeddedFromIds, cand))
            {
                embeddedAllSame = false;
            }
        }
    }

    if (embeddedAssigned && embeddedAllSame)
    {
        embeddedBCMapped = embeddedFromIds;
    }
    else if (embeddedAssigned && !embeddedAllSame && mode == Mode::Explicit)
    {
        std::cerr << prefix << "[BC] WARNING: multiple embedded boundaryId Field BCs differ; "
                  << "using embeddedBCById map and leaving global EmbeddedBC as fallback only.\n";
    }

    auto getSideBoundaryId = [&](const std::string& side) -> int {
        auto itS = sideAssign.find(side);
        if (itS == sideAssign.end() || !itS->second.set)
            return -1;
        return itS->second.id;
    };
    result.domainIds = {{getSideBoundaryId("left"), getSideBoundaryId("right"),
                         getSideBoundaryId("bottom"), getSideBoundaryId("top")}};
    return result;
}

Explicit prepareExplicit(const PicConfig& cfg, const DictionaryGeometry& dictGeometry, bool haveGeometry)
{
    Explicit result;
    result.fields = mapFields(cfg, dictGeometry, haveGeometry, Mode::Explicit);
    result.particles.see = cfg.see;
    const auto& domainBCMapped = result.fields.domain;
    const auto& embeddedBCMapped = result.fields.embedded;
    const auto& embeddedBCByIdMapped = result.fields.embeddedById;
    const int idLeft = result.fields.domainIds[0], idRight = result.fields.domainIds[1];
    const int idBottom = result.fields.domainIds[2], idTop = result.fields.domainIds[3];
    auto& periodicX = result.periodicX;
    auto& periodicY = result.periodicY;
    auto& seeCfg = result.particles.see;
    auto& seeBoundaryIds = result.particles.embeddedIds;
    if (haveGeometry)
    {
        if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
        {
            auto bcTypeName = [](PicBCType t) -> const char* {
                switch (t)
                {
                case PicBCType::Dirichlet: return "Dirichlet";
                case PicBCType::Neumann:   return "Neumann";
                case PicBCType::Robin:     return "Robin";
                case PicBCType::Periodic:  return "Periodic";
                default:                   return "Unknown";
                }
            };
            auto timeBrief = [](const PicBoundarySide::TimeFunction& tf) -> std::string {
                if (!tf.enabled)
                    return "time=off";
                std::ostringstream oss;
                oss << "time=" << tf.type << " dc=" << tf.dc << " ac=" << tf.ac
                    << " f=" << tf.freq_hz << "Hz phase=" << tf.phase_rad << " t0=" << tf.t0_s;
                return oss.str();
            };
            auto sideDesc = [&](const PicBoundarySide& s) -> std::string {
                std::ostringstream oss;
                oss << bcTypeName(s.type) << " value=" << s.value;
                if (s.timeDirichlet.enabled && s.type == PicBCType::Dirichlet)
                    oss << " " << timeBrief(s.timeDirichlet);
                return oss.str();
            };
            auto embeddedDesc = [&](const PicEmbeddedBC& e) -> std::string {
                std::ostringstream oss;
                oss << bcTypeName(e.type) << " a=" << e.a << " b=" << e.b << " g=" << e.g;
                if (e.timeDirichlet.enabled && e.type == PicBCType::Dirichlet)
                    oss << " " << timeBrief(e.timeDirichlet);
                return oss.str();
            };

            std::cout << "[ExplicitPIC][BC] domainBC mapped: left=" << sideDesc(domainBCMapped.left)
                      << ", right=" << sideDesc(domainBCMapped.right)
                      << ", bottom=" << sideDesc(domainBCMapped.bottom)
                      << ", top=" << sideDesc(domainBCMapped.top) << "\n";

            if (!embeddedBCByIdMapped.empty())
            {
                std::cout << "[ExplicitPIC][BC] embeddedBCById size=" << embeddedBCByIdMapped.size() << "\n";
                for (const auto& kv : embeddedBCByIdMapped)
                {
                    std::cout << "  [ExplicitPIC][BC] embedded id=" << kv.first
                              << " -> " << embeddedDesc(kv.second) << "\n";
                }
            }
            else if (embeddedBCMapped.enabled)
            {
                std::cout << "[ExplicitPIC][BC] embeddedBC global: " << embeddedDesc(embeddedBCMapped) << "\n";
            }
        }

        const bool bcPeriodicX = (domainBCMapped.left.type == PicBCType::Periodic) &&
                                 (domainBCMapped.right.type == PicBCType::Periodic);
        const bool bcPeriodicY = (domainBCMapped.bottom.type == PicBCType::Periodic) &&
                                 (domainBCMapped.top.type == PicBCType::Periodic);

        auto isMate = [&](int a, int b) -> bool {
            if (a <= 0 || b <= 0)
                return false;
            if (!dictGeometry.isPeriodicBoundary(a))
                return false;
            return dictGeometry.periodicMate(a) == b;
        };

        const bool geomPeriodicX = isMate(idLeft, idRight) || isMate(idRight, idLeft);
        const bool geomPeriodicY = isMate(idBottom, idTop) || isMate(idTop, idBottom);

        periodicX = bcPeriodicX && geomPeriodicX;
        periodicY = bcPeriodicY && geomPeriodicY;

        if ((bcPeriodicX && !geomPeriodicX) || (bcPeriodicY && !geomPeriodicY))
        {
            std::cerr << "[ExplicitPIC][periodic] WARNING: Field BC requests Periodic on a side pair, "
                         "but geometry.dict PeriodicPairs does not define a matching mate. "
                         "Periodic wrap/ghost will remain disabled for safety.\n";
        }

        if (ConsoleControl::enabled(ConsoleControl::Level::kInfo))
        {
            std::cout << "[ExplicitPIC][periodic] enabled: periodicX=" << (periodicX ? "true" : "false")
                      << " periodicY=" << (periodicY ? "true" : "false")
                      << " (FieldBC: x=" << (bcPeriodicX ? "Periodic" : "non")
                      << ", y=" << (bcPeriodicY ? "Periodic" : "non")
                      << "; GeometryPairs: x=" << (geomPeriodicX ? "ok" : "missing")
                      << ", y=" << (geomPeriodicY ? "ok" : "missing")
                      << ")\n";
        }

        // -----------------------------------------------------------------
        // Particle BC mapping (SEE + placeholder for future wall models)
        // -----------------------------------------------------------------
        for (const auto& bc : cfg.particleBCById)
        {
            auto it = dictGeometry.boundaryById.find(bc.id);
            if (it == dictGeometry.boundaryById.end())
            {
                std::cerr << "[ExplicitPIC][BC] WARNING: boundaryId " << bc.id
                          << " not found in geometry.dict; ignoring Particle BC.\n";
                continue;
            }
            const auto& ent = it->second;
            if (ent.ignored)
                continue;

            // Merge SEE configs: per-boundary overrides global enable, but keep the global parameters.
            if (bc.seeSpecified)
            {
                if (!seeCfg.enabled)
                    seeCfg = bc.see;
                else if (bc.see.enabled)
                    seeCfg.enabled = true;
            }

            if (bc.seeSpecified && bc.see.enabled)
            {
                if (ent.kind == GeometricBoundary::domain)
                {
                    const std::string side = inferSide(dictGeometry, ent);
                    if (side.empty())
                    {
                        std::cerr << "[ExplicitPIC][SEE] WARNING: cannot infer domain side for boundaryId "
                                  << bc.id << "; SEE not applied.\n";
                    }
                    else
                    {
                        if (std::find(seeCfg.cathodeSides.begin(), seeCfg.cathodeSides.end(), side) == seeCfg.cathodeSides.end())
                            seeCfg.cathodeSides.push_back(side);
                    }
                }
                else
                {
                    seeBoundaryIds.insert(bc.id);
                }
            }

            if (!bc.model.empty())
            {
                const std::string model = lower(bc.model);
                if (model != "absorb")
                {
                    std::cerr << "[ExplicitPIC][BC] WARNING: Particle boundary model '" << bc.model
                              << "' is not implemented (only absorb). boundaryId=" << bc.id << "\n";
                }
            }
        }
    }

    // -------------------------------------------------------------
    // Startup logs: boundaryId -> BC (Field/Particle)
    // -------------------------------------------------------------
    auto bcTypeName = [](PicBCType t) -> const char* {
        switch (t)
        {
            case PicBCType::Dirichlet: return "Dirichlet";
            case PicBCType::Neumann:   return "Neumann";
            case PicBCType::Robin:     return "Robin";
            case PicBCType::Periodic:  return "Periodic";
        }
        return "Unknown";
    };
    auto geomKindName = [](GeometricBoundary k) -> const char* {
        switch (k)
        {
            case GeometricBoundary::domain: return "domainBoundary";
            case GeometricBoundary::fixed:  return "fixedBoundary";
            case GeometricBoundary::moving: return "movingBoundary";
        }
        return "unknown";
    };

    std::cerr << "[ExplicitPIC][BC] Field BCs by boundaryId (from pic.dict):\n";
    for (const auto& bc : cfg.fieldBCById)
    {
        std::string gname;
        std::string gkind;
        if (haveGeometry)
        {
            auto it = dictGeometry.boundaryById.find(bc.id);
            if (it != dictGeometry.boundaryById.end())
            {
                gname = it->second.name;
                gkind = geomKindName(it->second.kind);
            }
        }
        std::cerr << "  - id=" << bc.id;
        if (!gkind.empty()) std::cerr << " (" << gkind << ")";
        if (!gname.empty()) std::cerr << " name='" << gname << "'";
        std::cerr << ": " << bcTypeName(bc.bc.type) << " value=" << bc.bc.value;
        if (bc.bc.timeDirichlet.enabled)
        {
            const auto& tf = bc.bc.timeDirichlet;
            std::cerr << " waveform=" << tf.type << " dc=" << tf.dc << " ac=" << tf.ac
                      << " f=" << tf.freq_hz << " phase=" << tf.phase_rad << " t0=" << tf.t0_s;
        }
        std::cerr << "\n";
    }

    std::cerr << "[ExplicitPIC][BC] Particle BCs by boundaryId (from pic.dict):\n";
    for (const auto& bc : cfg.particleBCById)
    {
        std::string gname;
        std::string gkind;
        if (haveGeometry)
        {
            auto it = dictGeometry.boundaryById.find(bc.id);
            if (it != dictGeometry.boundaryById.end())
            {
                gname = it->second.name;
                gkind = geomKindName(it->second.kind);
            }
        }
        std::cerr << "  - id=" << bc.id;
        if (!gkind.empty()) std::cerr << " (" << gkind << ")";
        if (!gname.empty()) std::cerr << " name='" << gname << "'";
        std::cerr << ": model=" << (bc.model.empty() ? "absorb" : bc.model);
        if (bc.seeSpecified)
            std::cerr << " SEE.enabled=" << (bc.see.enabled ? "true" : "false")
                      << " gamma=" << bc.see.gamma;
        std::cerr << "\n";
    }

    return result;
}

Particles prepareImplicitParticles(const PicConfig& cfg, const DictionaryGeometry& dictGeometry)
{
    Particles result;
    result.see = cfg.see;
    auto& seeCfg = result.see;
    bool seeTemplateSelected = seeCfg.enabled;
    auto& seeBoundaryIds = result.embeddedIds;
    seeBoundaryIds.reserve(cfg.particleBCById.size());
    for (const auto& bc : cfg.particleBCById)
    {
        auto it = dictGeometry.boundaryById.find(bc.id);
        if (it == dictGeometry.boundaryById.end())
            continue;
        const auto& ent = it->second;
        if (ent.ignored)
            continue;

        if (bc.seeSpecified)
        {
            if (bc.see.enabled && !seeTemplateSelected)
            {
                seeCfg = bc.see;
                seeTemplateSelected = true;
            }
            else if (bc.see.enabled)
                seeCfg.enabled = true;
        }

        if (bc.seeSpecified && bc.see.enabled && ent.kind != GeometricBoundary::domain)
            seeBoundaryIds.insert(bc.id);
    }

    return result;
}
}
