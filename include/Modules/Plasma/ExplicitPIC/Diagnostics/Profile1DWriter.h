#pragma once

#include <cstddef>
#include <functional>
#include <fstream>
#include <string>
#include <vector>

#include "Core/Base/Vector2D.h"
#include "Core/Field/CellField.h"
#include "Core/IO/ConsoleControl.h" // toLowerCopy
#include "Modules/Plasma/Particles/ParticleType.h"

class QuadtreeMesh;
class UnitsContext;
class NormalizationContext;
class CouplingWeightsCache;

namespace ExplicitPIC
{
namespace Diagnostics
{

struct Profile1DWriterConfig
{
    bool enabled = false;
    int every = 50;
    int nbinsX = 128;
    std::string path = "./output/profile1d.csv";

    bool includePhi = true;
    bool includeE = true;
    bool includeRho = true;
};

class Profile1DWriter
{
public:
    using BuildCacheFn = std::function<void(std::size_t speciesIndex, CouplingWeightsCache& cache)>;

    Profile1DWriter() = default;
    explicit Profile1DWriter(const Profile1DWriterConfig& cfg) : cfg_(cfg) {}

    void setConfig(const Profile1DWriterConfig& cfg) { cfg_ = cfg; }
    const Profile1DWriterConfig& config() const { return cfg_; }

    bool enabled() const { return cfg_.enabled; }
    bool isOpen() const { return out_.is_open(); }

    void open(const std::string& filePath, int startStep = 0);
    void close();

    void write(int step,
               double time_code,
               double time_SI_s,
               const QuadtreeMesh& mesh,
               const UnitsContext* units,
               const NormalizationContext* norm,
               const CellField<double>& phi,
               const CellField<Vector2D>& Ecell,
               const CellField<double>& rho_tilde,
               const std::vector<ParticleType>& species,
               const std::vector<double>& macroWeight,
               const BuildCacheFn& buildCacheForSpecies,
               const std::vector<double>* alphaByCellIndex,
               bool writeSIOutputs);

private:
    Profile1DWriterConfig cfg_;
    std::ofstream out_;
    bool headerWritten_ = false;

    static std::string safeKey(const std::string& name);
    void writeHeader(const std::vector<ParticleType>& species, bool writeSIOutputs);
};

} // namespace Diagnostics
} // namespace ExplicitPIC
