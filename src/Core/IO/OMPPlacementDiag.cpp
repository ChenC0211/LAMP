#include "Core/IO/OMPPlacementDiag.h"

#if defined(_OPENMP) && defined(_WIN32)

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string>

#include <omp.h>
#include <windows.h>

namespace
{
std::mutex g_diag_mtx;
bool g_printed_once = false;
}

void PrintOMPPlacementOnce(const char* tag)
{
    const char* safeTag = (tag && tag[0] != '\0') ? tag : "omp";
    char* enabled = nullptr;
    std::size_t enabledLen = 0;
    if (_dupenv_s(&enabled, &enabledLen, "QUANTUM_OMP_PLACEMENT_DIAG") != 0 || enabled == nullptr)
        return;
    std::string flag(enabled);
    std::free(enabled);
    std::transform(flag.begin(), flag.end(), flag.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (!(flag == "1" || flag == "true" || flag == "yes" || flag == "on"))
        return;

    {
        std::lock_guard<std::mutex> lock(g_diag_mtx);
        if (g_printed_once)
            return;
        g_printed_once = true;
    }

#pragma omp parallel
    {
        PROCESSOR_NUMBER pn{};
        GetCurrentProcessorNumberEx(&pn);
        const int tid = omp_get_thread_num();

        std::lock_guard<std::mutex> lock(g_diag_mtx);
        std::cout << "[" << safeTag << "] tid=" << tid
                  << " group=" << static_cast<unsigned int>(pn.Group)
                  << " cpu=" << static_cast<unsigned int>(pn.Number)
                  << "\n";
    }
}

#else

void PrintOMPPlacementOnce(const char* tag)
{
    (void)tag;
}

#endif
