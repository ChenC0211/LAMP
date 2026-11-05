#include <algorithm>
#include <cctype>
#include <clocale>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#include <cstdio>
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>

namespace {

LONG WINAPI quantumUnhandledExceptionFilter(EXCEPTION_POINTERS* ep)
{
    std::ofstream ofs("./output/lamp_crash_unhandled.log", std::ios::app);
    if (!ofs)
        return EXCEPTION_CONTINUE_SEARCH;

    const DWORD tid = GetCurrentThreadId();
    const DWORD pid = GetCurrentProcessId();
    const DWORD code = (ep && ep->ExceptionRecord) ? ep->ExceptionRecord->ExceptionCode : 0;
    const void* faultAddr = (ep && ep->ExceptionRecord) ? ep->ExceptionRecord->ExceptionAddress : nullptr;
    const uintptr_t base = reinterpret_cast<uintptr_t>(GetModuleHandleA(nullptr));
    const uintptr_t fault = reinterpret_cast<uintptr_t>(faultAddr);

    ofs << "[LAMP][crash] pid=" << pid
        << " tid=" << tid
        << " code=0x" << std::hex << code
        << " fault=0x" << fault;
    if (base != 0 && fault >= base)
        ofs << " rva=0x" << (fault - base);
    ofs << std::dec << "\n";

    void* frames[64] = {};
    USHORT n = CaptureStackBackTrace(0, 64, frames, nullptr);
    ofs << "[LAMP][crash] stack_frames=" << n << "\n";
    for (USHORT i = 0; i < n; ++i)
    {
        const uintptr_t a = reinterpret_cast<uintptr_t>(frames[i]);
        ofs << "  #" << i << " addr=0x" << std::hex << a;
        if (base != 0 && a >= base)
            ofs << " rva=0x" << (a - base);
        ofs << std::dec << "\n";
    }
    ofs.flush();
    return EXCEPTION_CONTINUE_SEARCH;
}

} // namespace
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Core/Base/Base.h"
#include "RunCase.h"

using namespace std;

static void configure_console_utf8()
{
#ifdef _WIN32
    // 让 Windows 控制台按 UTF-8 解释/显示字节流，避免中文输出乱码
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");
#endif
}

namespace {

static bool stdinIsTTY()
{
#ifdef _WIN32
    return _isatty(_fileno(stdin)) != 0;
#else
    return isatty(fileno(stdin)) != 0;
#endif
}

static std::string trimCopy(std::string s)
{
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

void applyDefaultOpenMPEnv()
{
#ifdef _OPENMP
    auto setIfMissing = [](const char* name, const char* value) {
        const char* cur = std::getenv(name);
        if (cur && cur[0] != '\0')
            return;
#ifdef _WIN32
        // Keep both the Windows environment and the CRT getenv cache in sync;
        // the OpenMP runtime reads the latter when it initializes.
        if (_putenv_s(name, value) != 0)
            throw std::runtime_error(std::string("Cannot set OpenMP default: ") + name);
#else
        setenv(name, value, 0);
#endif
    };

    // Let the OS place independent simulations by default. Close binding can
    // pin every process's serial master to the same first core. An explicit
    // launcher/user binding policy is retained; cores then defines its places.
    setIfMissing("OMP_PLACES", "cores");
    setIfMissing("OMP_PROC_BIND", "false");
    setIfMissing("OMP_MAX_ACTIVE_LEVELS", "1");
    setIfMissing("OMP_DYNAMIC", "false");


#endif
}

void printUsage(const char* exe)
{
    cout
        << "Usage: " << exe << " [--case <caseDir>] [--check-case] [--threads <N>] [--help]\n\n"
        << "Options:\n"
        << "  --case <caseDir>   Run from the specified case directory.\n"
        << "                    A case directory must contain a 'system/' subdirectory.\n"
        << "  --check-case       Validate formal Hydrodynamics, Interface or Electric inputs without writing.\n"
        << "  --threads <N>     Override run.dict Parallel.threads (1..256; default 2).\n"
        << "  --help, -h         Show this help and exit.\n\n"
        << "Notes:\n"
        << "  - If --case is not provided and './system' is missing:\n"
        << "      * In an interactive terminal, the program will prompt you to select\n"
        << "        a case under './cases/'. Press ENTER to accept the default.\n"
        << "      * In a non-interactive run, it will try './cases/default' (searching\n"
        << "        up a few parents).\n";
}

// Resolve the case directory from CLI args. Returns empty string if not set.
string parseCaseArg(int argc, char** argv, bool& showHelp, bool& checkCase, std::optional<int>& threads)
{
    showHelp = false;
    checkCase = false;
    string caseDir;
    for (int i = 1; i < argc; ++i)
    {
        const string a = argv[i];
        if (a == "--threads" || a.rfind("--threads=",0)==0)
        {
            if(threads)throw std::invalid_argument("Duplicate --threads option");
            const std::string value=a=="--threads"
                ?(i+1<argc?argv[++i]:std::string{}):a.substr(10);
            if(value.empty() || value.find_first_not_of("0123456789")!=std::string::npos
                || value.size()>3)throw std::invalid_argument("--threads requires an integer in [1,256]");
            const int count=std::stoi(value);
            if(count<1 || count>256)throw std::invalid_argument("--threads requires an integer in [1,256]");
            threads=count;continue;
        }
        if (a == "--check-case")
        {
            checkCase = true;
            continue;
        }
        if (a == "--help" || a == "-h")
        {
            showHelp = true;
            return {};
        }

        if (a == "--case")
        {
            if (i + 1 < argc)
                caseDir = argv[++i];
            else
                showHelp = true;
            continue;
        }

        // Support --case=path
        const string prefix = "--case=";
        if (a.rfind(prefix, 0) == 0)
        {
            caseDir = a.substr(prefix.size());
            continue;
        }
        throw std::invalid_argument("Unknown option: "+a);
    }
    return caseDir;
}

bool maybeEnterCaseDirectory(const string& caseDir)
{
    namespace fs = std::filesystem;

    const fs::path cwd = fs::current_path();
    fs::path target;

    auto findCasesDirUpwards = [](fs::path start, int maxUp = 4) -> fs::path {
        for (int k = 0; k <= maxUp; ++k)
        {
            fs::path cand = start / "cases";
            if (fs::exists(cand) && fs::is_directory(cand))
                return cand;
            if (!start.has_parent_path())
                break;
            start = start.parent_path();
        }
        return {};
    };

    auto findDefaultCaseUpwards = [](fs::path start, int maxUp = 4) -> fs::path {
        for (int k = 0; k <= maxUp; ++k)
        {
            fs::path cand = start / "cases" / "default";
            if (fs::exists(cand / "system"))
                return cand;
            if (!start.has_parent_path())
                break;
            start = start.parent_path();
        }
        return {};
    };

    auto listCases = [](const fs::path& casesDir) -> std::vector<fs::path> {
        std::vector<fs::path> out;
        if (!fs::exists(casesDir) || !fs::is_directory(casesDir))
            return out;
        for (const auto& e : fs::directory_iterator(casesDir))
        {
            if (!e.is_directory())
                continue;
            const fs::path p = e.path();
            if (fs::exists(p / "system") && fs::exists(p / "system" / "run.dict"))
                out.push_back(p);
        }
        std::sort(out.begin(), out.end(), [](const fs::path& a, const fs::path& b) {
            return a.filename().string() < b.filename().string();
        });
        return out;
    };

    auto promptSelectCase = [&](const fs::path& casesDir) -> fs::path {
        const auto cases = listCases(casesDir);
        if (cases.empty())
            return {};

        // Default: prefer 'default' if present, else use the first case.
        int defaultIdx = 0;
        for (size_t i = 0; i < cases.size(); ++i)
        {
            if (cases[i].filename().string() == "default")
            {
                defaultIdx = static_cast<int>(i);
                break;
            }
        }

        cout << "[LAMP] Select a case under: " << casesDir.string() << "\n";
        for (size_t i = 0; i < cases.size(); ++i)
        {
            cout << "  [" << i << "] " << cases[i].filename().string();
            if (static_cast<int>(i) == defaultIdx)
                cout << "  (default)";
            cout << "\n";
        }

        cout << "Enter case index or name (ENTER = default): ";
        std::string input;
        std::getline(std::cin, input);
        input = trimCopy(input);
        if (input.empty())
            return cases[defaultIdx];

        // If digits: treat as index.
        bool allDigits = !input.empty() && std::all_of(input.begin(), input.end(), [](unsigned char c) {
            return std::isdigit(c) != 0;
        });
        if (allDigits)
        {
            const int idx = std::stoi(input);
            if (idx >= 0 && idx < static_cast<int>(cases.size()))
                return cases[idx];
            cerr << "[LAMP] WARNING: invalid index, using default.\n";
            return cases[defaultIdx];
        }

        // Otherwise: treat as case directory name (relative to casesDir), or a path.
        fs::path cand = fs::path(input);
        if (cand.is_relative())
        {
            // Allow users to type either 'electrospray' or 'cases/electrospray'.
            if (cand.begin() != cand.end() && (*cand.begin()).string() == "cases")
                cand = fs::path(input); // keep as-is
            else
                cand = casesDir / cand;
        }

        if (fs::exists(cand / "system") && fs::exists(cand / "system" / "run.dict"))
            return fs::canonical(cand);

        cerr << "[LAMP] WARNING: case '" << input << "' not found or missing system/run.dict, using default.\n";
        return cases[defaultIdx];
    };

    if (!caseDir.empty())
    {
        target = fs::path(caseDir);
        if (target.is_relative())
            target = cwd / target;
    }
    else
    {
        // No --case: if ./system exists, keep legacy behavior.
        if (fs::exists(cwd / "system"))
            return true;

        // Interactive: prompt for case selection.
        if (stdinIsTTY())
        {
            fs::path casesDir = findCasesDirUpwards(cwd);
            if (!casesDir.empty())
            {
                fs::path chosen = promptSelectCase(casesDir);
                if (!chosen.empty())
                    target = chosen;
            }
        }

        // Non-interactive (or prompt failed): fall back to default case.
        if (target.empty())
        {
            fs::path defaultCase = findDefaultCaseUpwards(cwd);
            if (!defaultCase.empty())
                target = defaultCase;
            else
                return true; // caller will fail later when reading dicts.
        }
    }

    if (!target.empty())
    {
        if (!fs::exists(target / "system"))
        {
            cerr << "[LAMP] ERROR: case directory '" << target.string()
                 << "' does not contain a 'system' subdirectory.\n";
            return false;
        }
        fs::current_path(target);
        cout << "[LAMP] Case: " << fs::current_path().string() << endl;
    }

    return true;
}

} // namespace

int main(int argc, char** argv) try
{
    applyDefaultOpenMPEnv();
    configure_console_utf8();

    bool showHelp = false;
    bool checkCase = false;
    std::optional<int> threadsOverride;
    const string caseDir = parseCaseArg(argc, argv, showHelp, checkCase, threadsOverride);
    if (showHelp)
    {
        printUsage(argv[0]);
        return 0;
    }

    if (!maybeEnterCaseDirectory(caseDir))
        return 2;

    // Public configurations and small-test fixtures are immutable inputs.
    // The generic launcher clones them into runs/ before invoking this entry.
    if (!checkCase) {
        for (auto p=std::filesystem::current_path(); p.has_parent_path(); p=p.parent_path()) {
            if ((p.filename()=="cases" && std::filesystem::exists(p.parent_path()/"CMakeLists.txt")) ||
                (p.filename()=="fixtures" && p.parent_path().filename()=="tests")) {
                throw std::invalid_argument("Use tools\\run_case.cmd <case-name>; physical runs require a copied run directory.");
            }
            if(p==p.parent_path()) break;
        }
    }

#ifdef _WIN32
    SetUnhandledExceptionFilter(quantumUnhandledExceptionFilter);
#endif

    // Thin entry: enter a case directory and hand over to the App dispatcher.
    return runCase(checkCase, threadsOverride);
}
catch (const std::exception& e)
{
    std::cerr << "[LAMP] ERROR: " << e.what() << std::endl;
    return 2;
}
