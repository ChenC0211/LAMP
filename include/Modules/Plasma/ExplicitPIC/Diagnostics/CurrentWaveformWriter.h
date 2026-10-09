#pragma once

#include <cstdint>
#include <fstream>
#include <string>

namespace ExplicitPIC
{
namespace Diagnostics
{

// Lightweight CSV writer for RF current waveform diagnostics.
//
// CSV header:
//   step,t_si_s,V_drive_V,J_cond_Am2,J_disp_Am2,J_total_Am2
//
// where J_total = J_cond + J_disp.
class CurrentWaveformWriter
{
public:
    CurrentWaveformWriter() = default;
    ~CurrentWaveformWriter();

    bool isOpen() const { return out_.is_open(); }

    void open(const std::string& filePath, int startStep = 0);
    void close();

    void write(std::int64_t step,
               double t_si_s,
               double V_drive_V,
               double J_cond_Am2,
               double J_disp_Am2);

private:
    std::ofstream out_;
    bool headerWritten_ = false;

    void writeHeader();
};

} // namespace Diagnostics
} // namespace ExplicitPIC
