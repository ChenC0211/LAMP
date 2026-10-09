#pragma once

#include <algorithm>
#include <cctype>
#include <string>

// -----------------------------------------------------------------------------
// ConsoleControl (O4-lite)
//
// A tiny, header-only console verbosity gate.
//
// - ExplicitPICModule sets the global config from pic.dict:
//     Output.console.level / Output.console.every
// - Other components (DepositQK / Poisson / etc.) consult this gate before
//   printing to std::cout.
//
// This is intentionally minimal and does *not* attempt to be a full logger.
// -----------------------------------------------------------------------------

namespace ConsoleControl
{

enum class Level : int
{
    kError = 0,
    kWarn = 1,
    kInfo = 2,
    kDebug = 3,
    kTrace = 4,
};

struct Config
{
    Level level = Level::kInfo;
    int every = 10; // for info-level periodic prints
};

// C++17 inline variable: one instance per program.
inline Config g_config{};

inline std::string toLowerCopy(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

inline Level parseLevel(const std::string& s, Level fallback = Level::kInfo)
{
    const std::string ls = toLowerCopy(s);
    if (ls == "error" || ls == "err")
        return Level::kError;
    if (ls == "warn" || ls == "warning")
        return Level::kWarn;
    if (ls == "info")
        return Level::kInfo;
    if (ls == "debug")
        return Level::kDebug;
    if (ls == "trace")
        return Level::kTrace;
    return fallback;
}

inline void setConfig(const std::string& levelStr, int every)
{
    g_config.level = parseLevel(levelStr, Level::kInfo);
    g_config.every = (every > 0) ? every : 1;
}

inline bool enabled(Level msgLevel)
{
    return static_cast<int>(g_config.level) >= static_cast<int>(msgLevel);
}

inline bool shouldPrintPeriodic(int step, int every)
{
    if (every <= 0)
        return false;
    if (step <= 0)
        return true;
    return (step % every) == 0;
}

// For physical summaries at info-level (and above).
inline bool shouldPrintInfoSummary(int step)
{
    return enabled(Level::kInfo) && shouldPrintPeriodic(step, g_config.every);
}

} // namespace ConsoleControl
