#include "Modules/Plasma/ExplicitPIC/Diagnostics/CurrentWaveformWriter.h"

#include <filesystem>
#include <system_error>
#include <iomanip>

#include "Modules/Plasma/Diagnostics/TimeSeriesIO.h"

namespace ExplicitPIC
{
namespace Diagnostics
{

CurrentWaveformWriter::~CurrentWaveformWriter()
{
    close();
}

void CurrentWaveformWriter::open(const std::string& filePath, int startStep)
{
    close();

    std::filesystem::path p(filePath);
    const std::filesystem::path parent = p.parent_path();
    if (!parent.empty())
    {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
        // Best-effort: even if directory creation fails, try opening the file.
        (void)ec;
    }

    headerWritten_ = TimeSeriesIO::openPreparedStepCsv(out_, p, startStep);
    if (!headerWritten_)
        writeHeader();
}

void CurrentWaveformWriter::close()
{
    if (out_.is_open())
    {
        out_.flush();
        out_.close();
    }
    headerWritten_ = false;
}

void CurrentWaveformWriter::write(std::int64_t step,
                                  double t_si_s,
                                  double V_drive_V,
                                  double J_cond_Am2,
                                  double J_disp_Am2)
{
    if (!out_.is_open())
        return;

    if (!headerWritten_)
        writeHeader();

    const double J_total_Am2 = J_cond_Am2 + J_disp_Am2;

    out_ << step << ','
         << std::setprecision(17) << std::scientific << t_si_s << ','
         << std::setprecision(17) << std::scientific << V_drive_V << ','
         << std::setprecision(17) << std::scientific << J_cond_Am2 << ','
         << std::setprecision(17) << std::scientific << J_disp_Am2 << ','
         << std::setprecision(17) << std::scientific << J_total_Am2
         << '\n';
}

void CurrentWaveformWriter::writeHeader()
{
    if (!out_.is_open() || headerWritten_)
        return;

    out_ << "step,t_si_s,V_drive_V,J_cond_Am2,J_disp_Am2,J_total_Am2\n";
    headerWritten_ = true;
}

} // namespace Diagnostics
} // namespace ExplicitPIC
