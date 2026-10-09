#include <iomanip>
#include <sstream>
#include "Modules/Plasma/IO/PicFieldOutput.h"
#include "Modules/Plasma/IO/PicFieldFrame.h"
#include "Modules/Plasma/IO/PicProfiles.h"
#include "Modules/Plasma/Boundary/PicEmbeddedGeometry.h"
#include "Core/Units/UnitsContext.h"
#include "Core/Units/NormalizationContext.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/VertField.h"
#include "Core/Field/CellField.h"
#include "Core/Geometry/LevelSet.h"
#include "Core/Mesh/FaceSegmentGraph.h"
#include "Solvers/Poisson/CutCellMetrics.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
int compareStableOutputCoord(double a, double b)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    const double tol = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    if (a < b - tol)
        return -1;
    if (a > b + tol)
        return 1;
    return 0;
}

std::vector<int> activeCells(const QuadtreeMesh& mesh, const std::vector<unsigned char>* cellMask)
{
    std::vector<int> activeCellIndices;
    activeCellIndices.reserve(static_cast<std::size_t>(mesh.cellIndexNum()));
    for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
    {
        if (cellMask != nullptr &&
            cellMask->size() == static_cast<std::size_t>(mesh.cellIndexNum()) &&
            ((*cellMask)[static_cast<std::size_t>(ci)] == 0u))
        {
            continue;
        }
        activeCellIndices.push_back(ci);
    }


    std::stable_sort(activeCellIndices.begin(), activeCellIndices.end(), [&](int a, int b) {
        const Cell& ca = mesh.cells(mesh.cellIndexToLabel[a]);
        const Cell& cb = mesh.cells(mesh.cellIndexToLabel[b]);
        const Vector2D alo = ca.centerCoord() - ca.size() * 0.5;
        const Vector2D blo = cb.centerCoord() - cb.size() * 0.5;

        int cmp = compareStableOutputCoord(alo.y, blo.y);
        if (cmp != 0)
            return cmp < 0;
        cmp = compareStableOutputCoord(alo.x, blo.x);
        if (cmp != 0)
            return cmp < 0;

        const Vector2D ahi = ca.centerCoord() + ca.size() * 0.5;
        const Vector2D bhi = cb.centerCoord() + cb.size() * 0.5;
        cmp = compareStableOutputCoord(ahi.y, bhi.y);
        if (cmp != 0)
            return cmp < 0;
        cmp = compareStableOutputCoord(ahi.x, bhi.x);
        if (cmp != 0)
            return cmp < 0;

        if (ca.level() != cb.level())
            return ca.level() < cb.level();
        return a < b;
    });
    return activeCellIndices;
}
}

void PicFieldOutput::writeExplicit(int step, double dt, const ExplicitInput& input,
    PicProfiles& profiles, bool writeSIProfile)
{
    std::ostringstream tag;
    tag << std::setw(6) << std::setfill('0') << step;
    if (!input.units || !input.normalization)
    {
        std::cerr << "[ExplicitPIC][output] ERROR: SI outputs require units.dict and normalization.\n";
        throw std::runtime_error("ExplicitPIC SI outputs require UnitsContext and NormalizationContext");
    }
    const auto& mesh = input.mesh;
    const auto& embedded = input.embedded;
    const auto& alpha = embedded.alphaByCell();
    std::vector<double> x(mesh.vertIndexNum(), 0.0), y(mesh.vertIndexNum(), 0.0);
    for (int i = 0; i < mesh.vertIndexNum(); ++i)
    {
        const auto& v = mesh.verts(mesh.vertIndexToLabel[i]);
        x[i] = input.units->LengthCodeToSI(v.coord().x);
        y[i] = input.units->LengthCodeToSI(v.coord().y);
    }
    const PicFieldSelection selection(PicFieldMode::Explicit, input.fields);
    PicFieldFrame frame(mesh, selection,
        PicOutputReconstruction::buildActiveCellMask(mesh, input.graph,
            embedded.metricsReady() ? &embedded.metrics() : nullptr,
            embedded.metricsReady() ? &alpha : nullptr), input.potential);
    auto buildCache = [&](std::size_t si, CouplingWeightsCache& cache) {
        PicCurrentWeights::buildDiagnostic(input.coupling, embedded, input.species[si], cache);
    };
    if (profiles.instantaneousDue(step))
    {
        const double timeCode = static_cast<double>(step) * dt;
        const double timeSI = input.units ? input.units->TimeCodeToSI(timeCode) : timeCode;
        const auto* profileAlpha = embedded.metricsReady() && alpha.size() == static_cast<std::size_t>(mesh.cellIndexNum())
            ? &alpha : nullptr;
        profiles.writeInstantaneous(step, timeCode, timeSI,
            {mesh, input.units, input.normalization, input.potential, input.electric, input.species, input.weights, profileAlpha},
            input.charge, buildCache, writeSIProfile);
    }
    frame.addElectrostatic(input.potential, input.charge, *input.units, *input.normalization);
    frame.addLevelSets(*input.units, embedded.ready() ? embedded.levelSet() : nullptr);
    frame.addParticles(input.species, input.weights, *input.normalization,
        embedded.metricsReady() ? &alpha : nullptr, buildCache);
    frame.addElectric(input.electric, *input.units);
    writeExplicitVertices(mesh, "fields_" + tag.str(), x, y, frame.names, frame.vertices, &frame.validCells);
}

void PicFieldOutput::writeExplicitVertices(const QuadtreeMesh& mesh,
                                    const std::string& nameNoExt,
                                    const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    const std::vector<std::string>& varNames,
                                    const std::vector<VertField<double>>& vars,
                                    const std::vector<unsigned char>* cellMask)
{
    if (static_cast<int>(x.size()) != mesh.vertIndexNum() ||
        static_cast<int>(y.size()) != mesh.vertIndexNum())
    {
        std::cerr << "[ExplicitPIC][output] ERROR: coordinate array size mismatch for '" << nameNoExt << "'\n";
        std::cerr << "  x.size=" << x.size() << " y.size=" << y.size()
                  << " mesh.vertIndexNum=" << mesh.vertIndexNum() << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuad: coordinate size mismatch");
    }
    if (varNames.size() != vars.size())
    {
        std::cerr << "[ExplicitPIC][output] ERROR: varNames/vars size mismatch for '" << nameNoExt << "'\n";
        std::cerr << "  varNames=" << varNames.size() << " vars=" << vars.size() << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuad: varNames/vars size mismatch");
    }

    const std::filesystem::path fieldsDir("./output/fields");
    std::filesystem::create_directories(fieldsDir);
    const std::string filePath = (fieldsDir / (nameNoExt + ".dat")).string();
    std::ofstream fp(filePath, std::ofstream::out);
    if (!fp.is_open())
    {
        std::cerr << "[ExplicitPIC][output] ERROR: failed to open output file: " << filePath << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuad: failed to open file");
    }

    fp << "Title= \"" << nameNoExt << "\"\n";
    fp << "Variables= \"X\" \"Y\" ";
    for (const auto& n : varNames)
        fp << "\"" << n << "\" ";
    fp << "\n";
    const auto activeCellIndices = activeCells(mesh, cellMask);
    std::vector<unsigned char> usedVerts(static_cast<std::size_t>(mesh.vertIndexNum()), 0u);


    for (const int ci : activeCellIndices)
    {
        const int cellLabel = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cellLabel);
        for (int j = 0; j < 4; ++j)
        {
            const Vert& v = mesh.verts(cell.cellToVert[j]);
            if (v.index() >= 0 && v.index() < mesh.vertIndexNum())
                usedVerts[static_cast<std::size_t>(v.index())] = 1u;
        }
    }

    std::vector<int> activeVertIndices;
    activeVertIndices.reserve(static_cast<std::size_t>(mesh.vertIndexNum()));
    std::vector<int> vertMap(static_cast<std::size_t>(mesh.vertIndexNum()), -1);
    for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
    {
        if (usedVerts[static_cast<std::size_t>(vi)] == 0u)
            continue;
        vertMap[static_cast<std::size_t>(vi)] = static_cast<int>(activeVertIndices.size()) + 1;
        activeVertIndices.push_back(vi);
    }

    fp << "ZONE N = " << activeVertIndices.size()
       << ", E = " << activeCellIndices.size()
       << ", F = FEPOINT, ET = QUADRILATERAL\n";

    for (const int i : activeVertIndices)
    {
        fp << x[i] << " " << y[i] << " ";
        for (std::size_t j = 0; j < vars.size(); ++j)
            fp << vars[j].value[i] << " ";
        fp << "\n";
    }

    for (const int ci : activeCellIndices)
    {
        const int cellLabel = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cellLabel);
        for (int j = 0; j < 4; ++j)
        {
            const Vert& v = mesh.verts(cell.cellToVert[j]);
            const int mapped = (v.index() >= 0 && v.index() < mesh.vertIndexNum())
                                 ? vertMap[static_cast<std::size_t>(v.index())]
                                 : -1;
            if (mapped <= 0)
            {
                std::cerr << "[ExplicitPIC][output] ERROR: invalid masked vertex connectivity for '" << nameNoExt
                          << "' cellIndex=" << ci << std::endl;
                throw std::runtime_error("writeTecplotFEPointQuad: invalid masked connectivity");
            }
            fp << mapped << " ";
        }
        fp << "\n";
    }
}

void PicFieldOutput::writeImplicitVertices(const QuadtreeMesh& mesh,
                                    const std::string& nameNoExt,
                                    const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    const std::vector<std::string>& varNames,
                                    const std::vector<VertField<double>>& vars,
                                    const std::vector<unsigned char>* cellMask)
{
    if (static_cast<int>(x.size()) != mesh.vertIndexNum() ||
        static_cast<int>(y.size()) != mesh.vertIndexNum())
    {
        std::cerr << "[ImplicitPIC][output] ERROR: coordinate array size mismatch for '" << nameNoExt << "'\n";
        std::cerr << "  x.size=" << x.size() << " y.size=" << y.size()
                  << " mesh.vertIndexNum=" << mesh.vertIndexNum() << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuad: coordinate size mismatch");
    }
    if (varNames.size() != vars.size())
    {
        std::cerr << "[ImplicitPIC][output] ERROR: varNames/vars size mismatch for '" << nameNoExt << "'\n";
        std::cerr << "  varNames=" << varNames.size() << " vars=" << vars.size() << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuad: varNames/vars size mismatch");
    }

    const std::filesystem::path filePath = std::filesystem::path("./output") / (nameNoExt + ".dat");
    const std::filesystem::path parent = filePath.parent_path();
    if (!parent.empty())
        std::filesystem::create_directories(parent);
    std::ofstream fp(filePath.string(), std::ofstream::out);
    if (!fp.is_open())
    {
        std::cerr << "[ImplicitPIC][output] ERROR: failed to open output file: " << filePath << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuad: failed to open file");
    }

    fp << "Title= \"" << nameNoExt << "\"\n";
    fp << "Variables= \"X\" \"Y\" ";
    for (const auto& n : varNames)
        fp << "\"" << n << "\" ";
    fp << "\n";
    const auto activeCellIndices = activeCells(mesh, cellMask);
    std::vector<unsigned char> usedVerts(static_cast<std::size_t>(mesh.vertIndexNum()), 0u);


    for (const int ci : activeCellIndices)
    {
        const int cellLabel = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cellLabel);
        for (int j = 0; j < 4; ++j)
        {
            const Vert& v = mesh.verts(cell.cellToVert[j]);
            usedVerts[static_cast<std::size_t>(v.index())] = 1u;
        }
    }

    std::vector<int> vertOldToNew(static_cast<std::size_t>(mesh.vertIndexNum()), -1);
    std::vector<int> activeVertIndices;
    activeVertIndices.reserve(static_cast<std::size_t>(mesh.vertIndexNum()));
    for (int vi = 0; vi < mesh.vertIndexNum(); ++vi)
    {
        if (usedVerts[static_cast<std::size_t>(vi)] == 0u)
            continue;
        activeVertIndices.push_back(vi);
    }

    std::stable_sort(activeVertIndices.begin(), activeVertIndices.end(), [&](int a, int b) {
        int cmp = compareStableOutputCoord(y[a], y[b]);
        if (cmp != 0)
            return cmp < 0;
        cmp = compareStableOutputCoord(x[a], x[b]);
        if (cmp != 0)
            return cmp < 0;
        return a < b;
    });
    for (std::size_t newIndex = 0; newIndex < activeVertIndices.size(); ++newIndex)
        vertOldToNew[static_cast<std::size_t>(activeVertIndices[newIndex])] = static_cast<int>(newIndex);

    fp << "ZONE N = " << activeVertIndices.size()
       << ", E = " << activeCellIndices.size()
       << ", F = FEPOINT, ET = QUADRILATERAL\n";

    for (const int i : activeVertIndices)
    {
        fp << x[i] << " " << y[i] << " ";
        for (std::size_t j = 0; j < vars.size(); ++j)
            fp << vars[j].value[i] << " ";
        fp << "\n";
    }

    for (const int ci : activeCellIndices)
    {
        const int cellLabel = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cellLabel);
        for (int j = 0; j < 4; ++j)
        {
            const Vert& v = mesh.verts(cell.cellToVert[j]);
            const int newIndex = vertOldToNew[static_cast<std::size_t>(v.index())];
            if (newIndex < 0)
                throw std::runtime_error("writeTecplotFEPointQuad: active cell references filtered vertex");
            fp << newIndex + 1 << " ";
        }
        fp << "\n";
    }
}

void PicFieldOutput::writeImplicitCells(const QuadtreeMesh& mesh,
                                                const std::string& nameNoExt,
                                                const std::vector<double>& x,
                                                const std::vector<double>& y,
                                                const std::vector<std::string>& varNames,
                                                const std::vector<CellField<double>>& vars,
                                                const std::vector<unsigned char>* cellMask)
{
    if (static_cast<int>(x.size()) != mesh.vertIndexNum() ||
        static_cast<int>(y.size()) != mesh.vertIndexNum())
    {
        std::cerr << "[ImplicitPIC][output] ERROR: coordinate array size mismatch for '" << nameNoExt << "'\n";
        std::cerr << "  x.size=" << x.size() << " y.size=" << y.size()
                  << " mesh.vertIndexNum=" << mesh.vertIndexNum() << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuadCellConstant: coordinate size mismatch");
    }
    if (varNames.size() != vars.size())
    {
        std::cerr << "[ImplicitPIC][output] ERROR: varNames/vars size mismatch for '" << nameNoExt << "'\n";
        std::cerr << "  varNames=" << varNames.size() << " vars=" << vars.size() << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuadCellConstant: varNames/vars size mismatch");
    }

    const std::filesystem::path filePath = std::filesystem::path("./output") / (nameNoExt + ".dat");
    const std::filesystem::path parent = filePath.parent_path();
    if (!parent.empty())
        std::filesystem::create_directories(parent);
    std::ofstream fp(filePath.string(), std::ofstream::out);
    if (!fp.is_open())
    {
        std::cerr << "[ImplicitPIC][output] ERROR: failed to open output file: " << filePath << std::endl;
        throw std::runtime_error("writeTecplotFEPointQuadCellConstant: failed to open file");
    }

    const auto activeCellIndices = activeCells(mesh, cellMask);


    fp << "Title= \"" << nameNoExt << "\"\n";
    fp << "Variables= \"X\" \"Y\" ";
    for (const auto& n : varNames)
        fp << "\"" << n << "\" ";
    fp << "\n";
    fp << "ZONE N = " << (4 * activeCellIndices.size())
       << ", E = " << activeCellIndices.size()
       << ", F = FEPOINT, ET = QUADRILATERAL\n";

    for (const int ci : activeCellIndices)
    {
        const int cellLabel = mesh.cellIndexToLabel[ci];
        const Cell& cell = mesh.cells(cellLabel);
        for (int j = 0; j < 4; ++j)
        {
            const Vert& v = mesh.verts(cell.cellToVert[j]);
            const int vi = v.index();
            fp << x[vi] << " " << y[vi] << " ";
            for (std::size_t k = 0; k < vars.size(); ++k)
                fp << vars[k].value[ci] << " ";
            fp << "\n";
        }
    }

    for (std::size_t e = 0; e < activeCellIndices.size(); ++e)
    {
        const int base = static_cast<int>(4 * e);
        fp << (base + 1) << " " << (base + 2) << " " << (base + 3) << " " << (base + 4) << "\n";
    }
}

void PicFieldOutput::writeMaskDiagnostics(const QuadtreeMesh& mesh, int step,
    const std::vector<unsigned char>& validCells, const CellField<double>& potential,
    const FaceSegmentGraph& faceGraph, const Solvers::Poisson::CutCellMetrics* cutMetricsPtr,
    const std::vector<double>* alphaByCell, const std::vector<LevelSet>* levelSets,
    const std::vector<int>& levelSetIndices, bool insideIsNegative)
{
    const bool outputHasAlphaMask = alphaByCell && alphaByCell->size() == static_cast<std::size_t>(mesh.cellIndexNum());
    const std::filesystem::path diagPath =
        std::filesystem::path("./output") / "meta" /
        (std::string("packed_mask_diag_") + std::to_string(step) + ".csv");
    std::filesystem::create_directories(diagPath.parent_path());
    std::ofstream diag(diagPath.string(), std::ios::out | std::ios::trunc);
    if (diag)
    {
        diag << "cellIndex,x_m,y_m,isLeaf,outputValid,cutActive,faceSegCount,poissonActive,alpha,allInsideAnyLS,phiCell,isCut,A_gamma,d_gamma,gammaBoundaryId,gammaLineIndex\n";
        for (int ci = 0; ci < mesh.cellIndexNum(); ++ci)
        {
            const int cellLabel = mesh.cellIndexToLabel[ci];
            const Cell& cell = mesh.cells(cellLabel);
            const Vector2D cc = cell.centerCoord();
            const bool isLeaf =
                (cell.son[0] == -1) && (cell.son[1] == -1) && (cell.son[2] == -1) && (cell.son[3] == -1);
            const bool outputValid = (validCells[static_cast<std::size_t>(ci)] != 0u);
            const bool cutActive = (cutMetricsPtr != nullptr) ? cutMetricsPtr->isActiveCellIndex(ci) : false;
            double alpha = 1.0;
            bool isCut = false;
            double A_gamma = 0.0;
            double d_gamma = 0.0;
            int gammaBoundaryId = 0;
            int gammaLineIndex = -1;
            const auto segRange = faceGraph.segmentRangeForCellIndex(ci);
            const int faceSegCount = static_cast<int>(segRange.second - segRange.first);
            const bool poissonActive = (faceSegCount > 0);
            if (outputHasAlphaMask)
            {
                alpha = (*alphaByCell)[static_cast<std::size_t>(ci)];
                if (!std::isfinite(alpha))
                    alpha = 1.0;
                alpha = std::clamp(alpha, 0.0, 1.0);
            }
            if (cutMetricsPtr != nullptr)
            {
                const auto& cm = cutMetricsPtr->cellMetricsByIndex();
                if (static_cast<std::size_t>(ci) < cm.size())
                {
                    isCut = cm[static_cast<std::size_t>(ci)].isCut;
                    A_gamma = cm[static_cast<std::size_t>(ci)].A_gamma;
                    d_gamma = cm[static_cast<std::size_t>(ci)].d_gamma;
                    gammaBoundaryId = cm[static_cast<std::size_t>(ci)].gammaBoundaryId;
                    gammaLineIndex = cm[static_cast<std::size_t>(ci)].gammaLineIndex;
                }
            }

            bool allInsideAnyLS = false;
            if (levelSets && !levelSetIndices.empty())
            {
                for (int lsI : levelSetIndices)
                {
                    if (lsI < 0 || lsI >= static_cast<int>((*levelSets).size()))
                        continue;
                    bool allInside = true;
                    for (int j = 0; j < 4; ++j)
                    {
                        const Vert& v = mesh.verts(cell.cellToVert[j]);
                        const double psi = (*levelSets)[lsI].value[v.index()];
                        const bool inside = insideIsNegative ? (psi < 0.0) : (psi > 0.0);
                        if (!inside)
                        {
                            allInside = false;
                            break;
                        }
                    }
                    if (allInside)
                    {
                        allInsideAnyLS = true;
                        break;
                    }
                }
            }

            diag << ci << ','
                 << cc.x << ','
                 << cc.y << ','
                 << (isLeaf ? 1 : 0) << ','
                 << (outputValid ? 1 : 0) << ','
                 << (cutActive ? 1 : 0) << ','
                 << faceSegCount << ','
                 << (poissonActive ? 1 : 0) << ','
                 << alpha << ','
                 << (allInsideAnyLS ? 1 : 0) << ','
                 << potential.value[ci] << ','
                 << (isCut ? 1 : 0) << ','
                 << A_gamma << ','
                 << d_gamma << ','
                 << gammaBoundaryId << ','
                 << gammaLineIndex
                 << "\n";
        }
    }
}

void PicFieldOutput::writeCutcellState(const std::filesystem::path& outputRoot,
    const std::filesystem::path& configuredPath, int stepOut, double time_s,
    bool onlyActive, const PicCutcellStateView::Input& input)
{
    std::filesystem::path dumpPath = configuredPath.empty()
        ? std::filesystem::path()
        : std::filesystem::path(configuredPath);
    if (dumpPath.empty())
    {
        std::ostringstream name;
        name << "cutcell_state_step_" << std::setw(6) << std::setfill('0') << stepOut << ".csv";
        dumpPath = outputRoot / name.str();
    }
    else if (dumpPath.is_relative())
    {
        const auto parent = dumpPath.parent_path();
        if (parent.empty() || parent == ".")
            dumpPath = outputRoot / dumpPath;
        else
            dumpPath = std::filesystem::path(".") / dumpPath;
    }

    const auto parent = dumpPath.parent_path();
    if (!parent.empty())
        std::filesystem::create_directories(parent);

    const PicCutcellStateView view(input);
    auto bcTypeName = [](const PicBCType type) -> const char* {
        switch (type)
        {
        case PicBCType::Dirichlet: return "Dirichlet";
        case PicBCType::Neumann:   return "Neumann";
        case PicBCType::Robin:     return "Robin";
        case PicBCType::Periodic:  return "Periodic";
        default:                   return "Unknown";
        }
    };
    std::ofstream os(dumpPath.string(), std::ios::out | std::ios::trunc);
    if (!os)
    {
        std::cerr << "[ImplicitPIC] WARNING: failed to create cutcell state dump: "
                  << dumpPath.string() << std::endl;
        return;
    }

    os << "step,time_s,cellIndex,cellLabel,level,boundary,center_x_code,center_y_code,"
          "center_x_m,center_y_m,area_code2,area_m2,active,isCut,alpha,A_gamma_code,"
          "A_gamma_m,d_gamma_code,d_gamma_m,n_gamma_x,n_gamma_y,gammaBoundaryId,"
          "gammaLineIndex,rho_code,rho_Cpm3,kappa_cell,phi_code,phi_V,Ex_code,Ey_code,"
          "Ex_Vpm,Ey_Vpm,EnGamma_code,EnGamma_Vpm,eb_bc_source,eb_bc_type,eb_bc_a,"
          "eb_bc_b,eb_bc_g";
    if (view.levelSetCount() > 0)
    {
        for (int lsI = 0; lsI < view.levelSetCount(); ++lsI)
        {
            os << ",psi_ls" << std::setw(2) << std::setfill('0') << (lsI + 1) << "_code";
        }
        os << std::setfill(' ');
    }
    os << "\n";
    os << std::setprecision(16) << std::fixed;

    for (int ci = 0; ci < view.size(); ++ci)
    {
        if (onlyActive && !view.active(ci)) continue;
        const auto r = view.sample(ci);
        os << stepOut << ',' << time_s << ','
           << r.ci << ',' << r.cLabel << ',' << r.level << ',' << r.boundary << ','
           << r.xc.x << ',' << r.xc.y << ',' << r.x_m << ',' << r.y_m << ','
           << r.areaCode << ',' << r.area_m2 << ','
           << (r.active ? 1 : 0) << ',' << (r.cm.isCut ? 1 : 0) << ','
           << r.cm.alpha << ',' << r.cm.A_gamma << ',' << r.aGamma_m << ','
           << r.cm.d_gamma << ',' << r.dGamma_m << ','
           << r.cm.n_gamma.x << ',' << r.cm.n_gamma.y << ','
           << r.cm.gammaBoundaryId << ',' << r.cm.gammaLineIndex << ','
           << r.rhoCode << ',' << r.rhoSI << ',' << r.kappa << ','
           << r.phiCode << ',' << r.phiV << ','
           << r.eCode.x << ',' << r.eCode.y << ','
           << r.exVpm << ',' << r.eyVpm << ','
           << r.enGammaCode << ',' << r.enGammaVpm << ','
           << r.ebSource << ',' << (r.eb ? bcTypeName(r.eb->type) : "none") << ','
           << (r.eb ? r.eb->a : 0.0) << ',' << (r.eb ? r.eb->b : 0.0) << ',' << (r.eb ? r.eb->g : 0.0)
           ;
        for (double psiCode : r.psi) os << ',' << psiCode;
        os << '\n';
    }
    os.flush();
    std::cout << "[ImplicitPIC] cutcell state dump -> " << dumpPath.string()
              << " at step=" << stepOut << std::endl;
}
