#include <algorithm>
#include <cctype>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Edges = std::set<std::pair<std::string, std::string>>;
static void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
static bool starts(const std::string& s, const std::string& prefix) {
    return s.compare(0, prefix.size(), prefix) == 0;
}
static std::string read(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    require(bool(f), "Cannot read " + p.generic_string());
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}
static bool sourceFile(const fs::path& p) {
    static const std::set<std::string> extensions{ ".h", ".hpp", ".hxx", ".inl", ".tpp", ".c", ".cc", ".cpp", ".cxx" };
    return extensions.count(p.extension().string()) != 0;
}
static bool translationUnit(const fs::path& p) {
    static const std::set<std::string> extensions{ ".c", ".cc", ".cpp", ".cxx" };
    return extensions.count(p.extension().string()) != 0;
}
static const std::set<std::string> legacyFiles{
    "include/Modules/Electricity/Legacy/Electricity.h",
    "include/Modules/Electricity/Legacy/Electricity_impl.h",
    "src/Modules/Electricity/Legacy/Electricity.cpp"
};
// Exact live compatibility edges, never a directory-wide exemption. Removing
// a consumer must also remove its allowance; new consumers fail by default.
static const Edges compatibility{
    {"src/Apps/Electricity/Compatibility/LegacyElectricityModule.cpp", "include/Modules/Electricity/Legacy/Electricity.h"},
    {"include/Modules/Electricity/Legacy/Electricity.h", "include/Modules/Hydrodynamics/LegacyNodalState.h"},
    {"include/Modules/Electricity/Legacy/Electricity.h", "include/Modules/Plasma/Particles/ParticleType.h"},
    {"include/Modules/Electricity/Legacy/Electricity_impl.h", "include/Modules/Hydrodynamics/LegacyNodalState.h"},
    {"include/Modules/Electricity/Legacy/Electricity_impl.h", "include/Modules/Electricity/LegacyNodalMaxwell.h"},
    {"include/Apps/Electrospray/Coupling/ElectrosprayInterfaceMotion.h", "include/Modules/Interface/Compatibility/InterfaceMotion.h"},
    {"include/Apps/Electrospray/IO/ElectrosprayHydroDiagnostics.h", "include/Modules/Interface/Compatibility/InterfacePhaseFlux.h"},
    {"include/Apps/Electrospray/Setup/ElectrosprayInterface.h", "include/Modules/Interface/Compatibility/InterfaceInitialization.h"},
    {"include/Apps/Electrospray/Setup/ElectrosprayModuleInputs.h", "include/Modules/Interface/Compatibility/InterfaceCaseConfig.h"},
    {"include/Apps/Hydrodynamics/Coupling/HydroCaseCoupling.h", "include/Modules/Interface/Compatibility/InterfaceNeckCap.h"},
    {"include/Apps/Hydrodynamics/HydroConservativeInterface.h", "include/Modules/Interface/Compatibility/InterfaceColumnMode.h"},
    {"include/Apps/Hydrodynamics/HydroConservativeInterface.h", "include/Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"},
    {"include/Apps/Hydrodynamics/HydroConservativeInterface.h", "include/Modules/Interface/Compatibility/InterfaceMeshTransfer.h"},
    {"include/Apps/Hydrodynamics/HydroConservativeInterface.h", "include/Modules/Interface/Compatibility/InterfaceNeckCap.h"},
    {"include/Apps/Hydrodynamics/HydroConservativeInterface.h", "include/Modules/Interface/Compatibility/InterfaceState.h"},
    {"include/Apps/Hydrodynamics/HydroConservativeInterface.h", "include/Modules/Interface/Compatibility/IO/DictionaryInterface.h"},
    {"include/Apps/Hydrodynamics/IO/HydroCaseInput.h", "include/Modules/Interface/Compatibility/InterfaceNodalConfig.h"},
    {"include/Modules/Hydrodynamics/HydroCompatibleMomentFlow.h", "include/Modules/Interface/Compatibility/InterfaceParametricNeck.h"},
    {"include/Modules/Hydrodynamics/HydroCompatibleMomentFlow.h", "include/Modules/Interface/Compatibility/InterfaceParametricProjection.h"},
    {"include/Modules/Hydrodynamics/HydroCompatibleMomentSpline.h", "include/Modules/Interface/Compatibility/InterfaceParametricTransport.h"},
    {"include/Modules/Hydrodynamics/HydroCompatibleMomentSpline.h", "include/Modules/Interface/Compatibility/InterfaceTransportInput.h"},
    {"include/Modules/Hydrodynamics/HydroCompatibleMomentTransport.h", "include/Modules/Interface/Compatibility/InterfaceParametricTransport.h"},
    {"include/Modules/Hydrodynamics/HydroInterfaceBoundary.h", "include/Modules/Interface/Compatibility/InterfaceNodalConfig.h"},
    {"include/Modules/Hydrodynamics/HydroInterfaceTransport.h", "include/Modules/Interface/Compatibility/InterfaceParametricTransport.h"},
    {"include/Modules/Hydrodynamics/HydroInterfaceTransport.h", "include/Modules/Interface/Compatibility/InterfaceTransportInput.h"},
    {"include/Modules/Hydrodynamics/HydroParametricStep.h", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"include/Modules/Hydrodynamics/IO/DictionaryHydrodynamics.h", "include/Modules/Interface/Compatibility/InterfaceNodalConfig.h"},
    {"include/Modules/Hydrodynamics/IO/HydroGeometryLoader.h", "include/Modules/Interface/Compatibility/IO/DictionaryInterface.h"},
    {"src/Apps/Electrohydrodynamics/Coupling/ConservativeColumnCoupling.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"},
    {"src/Apps/Electrohydrodynamics/Coupling/ConservativeColumnCoupling.cpp", "include/Modules/Interface/Compatibility/InterfaceSurfaceTransfer.h"},
    {"src/Apps/Electrohydrodynamics/Coupling/ConservativeColumnCoupling.cpp", "include/Modules/Interface/Compatibility/InterfaceTopology.h"},
    {"src/Apps/Electrohydrodynamics/Coupling/DielectricCoupling.cpp", "include/Modules/Interface/Compatibility/InterfaceNodalConfig.h"},
    {"src/Apps/Electrohydrodynamics/Coupling/DielectricCoupling.cpp", "include/Modules/Interface/Compatibility/InterfaceNodalState.h"},
    {"src/Apps/Electrospray/Coupling/ElectrosprayCapillarity.cpp", "include/Modules/Interface/Compatibility/InterfaceState.h"},
    {"src/Apps/Electrospray/Coupling/ElectrosprayMeshAdaptation.cpp", "include/Modules/Interface/Compatibility/InterfaceMeshTransfer.h"},
    {"src/Apps/Electrospray/ElectrosprayTimeLoop.cpp", "include/Modules/Interface/Compatibility/InterfacePhaseFlux.h"},
    {"src/Apps/Electrospray/IO/ElectrosprayInitialInterfaceDiagnostics.cpp", "include/Modules/Interface/Compatibility/InterfaceCaseConfig.h"},
    {"src/Apps/Electrospray/IO/ElectrosprayVersionedInput.cpp", "include/Modules/Interface/Compatibility/IO/DictionaryInterface.h"},
    {"src/Apps/Electrospray/IO/ElectrosprayVolumeDiagnostics.cpp", "include/Modules/Interface/Compatibility/InterfaceVolumeCorrection.h"},
    {"src/Apps/Electrospray/Setup/ElectrosprayCaseBindings.cpp", "include/Modules/Interface/Compatibility/InterfaceCaseConfig.h"},
    {"src/Apps/Electrospray/Setup/ElectrosprayElectrostatics.cpp", "include/Modules/Interface/Compatibility/InterfaceCaseConfig.h"},
    {"src/Apps/Electrospray/Setup/ElectrosprayInterface.cpp", "include/Modules/Interface/Compatibility/InterfaceCaseConfig.h"},
    {"src/Apps/Hydrodynamics/HydroCase.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"},
    {"src/Apps/Hydrodynamics/HydroCase.cpp", "include/Modules/Interface/Compatibility/InterfaceNodalState.h"},
    {"src/Apps/Hydrodynamics/HydroCase.cpp", "include/Modules/Interface/Compatibility/InterfaceTopology.h"},
    {"src/Apps/Hydrodynamics/HydroConservativeInterface.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativeGeometry.h"},
    {"src/Apps/Hydrodynamics/HydroConservativeInterface.cpp", "include/Modules/Interface/Compatibility/InterfaceMotion.h"},
    {"src/Apps/Hydrodynamics/HydroConservativeInterface.cpp", "include/Modules/Interface/Compatibility/InterfaceNeckCap.h"},
    {"src/Apps/Hydrodynamics/HydroConservativeInterface.cpp", "include/Modules/Interface/Compatibility/InterfaceNodalConfig.h"},
    {"src/Apps/Hydrodynamics/HydroConservativeInterface.cpp", "include/Modules/Interface/Compatibility/InterfaceTopology.h"},
    {"src/Apps/Hydrodynamics/HydroMomentCase.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Apps/Hydrodynamics/HydroMomentCase.cpp", "include/Modules/Interface/Compatibility/InterfacePrescribedData.h"},
    {"src/Apps/Hydrodynamics/HydroMomentCase.cpp", "include/Modules/Interface/Compatibility/IO/DictionaryInterface.h"},
    {"src/Apps/Hydrodynamics/IO/HydroCaseInput.cpp", "include/Modules/Interface/Compatibility/IO/DictionaryInterface.h"},
    {"src/Modules/Hydrodynamics/DictionaryHydrodynamics.cpp", "include/Modules/Interface/Compatibility/IO/DictionaryInterface.h"},
    {"src/Modules/Hydrodynamics/HydroCompatibleMomentFlow.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydroCompatibleMomentMassCurrent.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydroCompatibleMomentTransport.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydroConservativeStress.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydroConservativeTransport.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativePhase.h"},
    {"src/Modules/Hydrodynamics/HydroConservativeTransport.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydrodynamicsConservative.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativePhase.h"},
    {"src/Modules/Hydrodynamics/HydrodynamicsConservativeTransfer.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativePhase.h"},
    {"src/Modules/Hydrodynamics/HydrodynamicsConservativeTransfer.cpp", "include/Modules/Interface/Compatibility/InterfaceState.h"},
    {"src/Modules/Hydrodynamics/HydrodynamicsState.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativePhase.h"},
    {"src/Modules/Hydrodynamics/HydrodynamicsState.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydrodynamicsStep.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydroInterfaceBoundary.cpp", "include/Modules/Interface/Compatibility/InterfaceBoundaryGeometry.h"},
    {"src/Modules/Hydrodynamics/HydroInterfaceForce.cpp", "include/Modules/Interface/Compatibility/UniformInterfaceContour.h"},
    {"src/Modules/Hydrodynamics/HydroNodalTraction.cpp", "include/Modules/Interface/Compatibility/UniformInterfaceContour.h"},
    {"src/Modules/Hydrodynamics/HydroUniformLoad.cpp", "include/Modules/Interface/Compatibility/InterfaceConservativePhase.h"},
    {"src/Modules/Hydrodynamics/HydroUniformLoad.cpp", "include/Modules/Interface/Compatibility/InterfaceCurvature.h"},
    {"src/Modules/Hydrodynamics/HydroUniformLoad.cpp", "include/Modules/Interface/Compatibility/InterfaceParametricMotion.h"},
    {"src/Modules/Hydrodynamics/HydroUniformLoad.cpp", "include/Modules/Interface/Compatibility/InterfaceTopology.h"},
    {"src/Modules/Hydrodynamics/IO/HydroGeometryLoader.cpp", "include/Modules/Interface/Compatibility/InterfaceNodalConfig.h"},
    {"src/Modules/Hydrodynamics/IO/HydroGeometryLoader.cpp", "include/Modules/Interface/Compatibility/InterfacePrescribedData.h"},
};
static std::string layer(const std::string& file) {
    const auto slash = file.find('/');
    if (slash == std::string::npos || (file.substr(0, slash) != "src" && file.substr(0, slash) != "include")) return "Unknown";
    const auto p = file.substr(slash + 1);
    if (starts(p, "Core/")) return "Core";
    if (starts(p, "Solvers/")) return "Solvers";
    if (starts(p, "Apps/Runtime/")) return "Runtime";
    if (starts(p, "Apps/")) return "Apps";
    if (starts(p, "Modules/Electricity/Legacy/")) return legacyFiles.count(file) ? "Legacy" : "Unknown";
    if (starts(p, "Modules/Interface/Compatibility/")) return "InterfaceCompatibility";
    if (starts(p, "Modules/Hydrodynamics/Cch") || starts(p, "Modules/Hydrodynamics/Detail/Cch"))
        return "HydrodynamicsCch";
    for (const auto* m : {"Hydrodynamics", "Electricity", "Interface", "Plasma"})
        if (starts(p, std::string("Modules/") + m + '/')) return m;
    return "Unknown";
}
static bool allowed(const std::string& from, const std::string& to) {
    if (from == "Unknown" || to == "Unknown") return false;
    if (from == to || to == "Core") return true;
    if (to == "Legacy" || to == "InterfaceCompatibility") return false;
    if (from == "InterfaceCompatibility") return to == "Interface";
    if (from == "Apps") return true;
    if (from == "Runtime") return to == "Solvers";
    if (from == "Hydrodynamics") return to == "Interface" || to == "Solvers";
    if (from == "HydrodynamicsCch") return to == "Interface" || to == "Solvers";
    if (from == "Electricity" || from == "Plasma") return to == "Solvers";
    return false;
}
static bool allowedInclude(const std::string& from, const std::string& to) {
    return compatibility.count({from, to}) || allowed(layer(from), layer(to));
}

// Splice physical lines first, then mask comments and raw strings. Preserve
// ordinary quoted header names and newlines, including inactive #if branches.
// This is a source dependency guard, not a replacement C++ preprocessor.
static std::string directives(std::string s) {
    s = std::regex_replace(s, std::regex("\\\\\r?\n"), "");
    if (starts(s, "\xEF\xBB\xBF")) s.erase(0, 3);
    for (size_t i = 0; i < s.size();) {
        size_t end = i;
        if (s.compare(i, 2, "//") == 0) {
            end = s.find('\n', i + 2);
            if (end == std::string::npos) end = s.size();
        } else if (s.compare(i, 2, "/*") == 0) {
            end = s.find("*/", i + 2);
            require(end != std::string::npos, "Unterminated block comment");
            end += 2;
        } else if (s.compare(i, 2, "R\"") == 0) {
            const auto open = s.find('(', i + 2);
            require(open != std::string::npos && open - i <= 18, "Invalid raw string delimiter");
            const auto closing = ")" + s.substr(i + 2, open - i - 2) + "\"";
            end = s.find(closing, open + 1);
            require(end != std::string::npos, "Unterminated raw string");
            end += closing.size();
        } else if (std::isdigit(static_cast<unsigned char>(s[i])) && (i == 0 ||
            (!std::isalnum(static_cast<unsigned char>(s[i - 1])) && s[i - 1] != '_'))) {
            // Apostrophes in C++14 numeric separators are not character quotes.
            while (i < s.size() && (std::isalnum(static_cast<unsigned char>(s[i])) ||
                s[i] == '\'' || s[i] == '_' || s[i] == '.')) ++i;
            continue;
        } else if (s[i] == '\"' || s[i] == '\'') {
            const char quote = s[i++];
            while (i < s.size()) {
                if (s[i] == '\\') i += std::min(size_t(2), s.size() - i);
                else if (s[i++] == quote) break;
            }
            continue;
        } else { ++i; continue; }
        while (i < end) { if (s[i] != '\n') s[i] = ' '; ++i; }
    }
    return s;
}
struct Include { std::string name; bool quoted; };
static std::vector<Include> includes(const std::string& text) {
    std::vector<Include> result;
    std::istringstream lines(directives(text));
    static const std::regex directive(R"(^\s*(?:#|%:)\s*(include|include_next|import)\b(.*)$)");
    static const std::regex header(R"re(^\s*("([^"\r\n]+)"|<([^>\r\n]+)>)\s*$)re");
    for (std::string line; std::getline(lines, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::smatch m, h;
        if (!std::regex_match(line, m, directive)) continue;
        require(m[1] == "include", "Unsupported dependency directive: " + line);
        const std::string tail = m[2];
        require(std::regex_match(tail, h, header), "Include must name a literal header: " + line);
        result.push_back({h[2].matched ? h[2].str() : h[3].str(), h[2].matched});
    }
    return result;
}
static std::string projectRelative(const fs::path& p, const fs::path& root) {
    return fs::weakly_canonical(p).lexically_relative(root).generic_string();
}
static std::string resolve(const fs::path& root, const std::string& file, const Include& inc) {
    require(!fs::path(inc.name).is_absolute(), file + ": absolute include " + inc.name);
    // include/Apps aliases are resolved even for modules, so aliases cannot
    // hide reverse dependencies. Actual missing search paths fail compilation.
    std::vector<fs::path> candidates;
    if (inc.quoted) candidates.push_back((root / file).parent_path() / inc.name);
    candidates.push_back(root / "include" / inc.name);
    candidates.push_back(root / "include/Apps" / inc.name);
    for (const auto& p : candidates) {
        if (!fs::is_regular_file(p)) continue;
        const auto dest = projectRelative(p, root);
        require(!starts(dest, "../"), file + ": include escapes project: " + inc.name);
        return dest;
    }
    require(!inc.quoted && !starts(inc.name, "Core/") && !starts(inc.name, "Modules/") &&
        !starts(inc.name, "Solvers/") && !starts(inc.name, "Apps/") && !starts(inc.name, "Runtime/"),
        file + ": unresolved project header " + inc.name);
    return {}; // External angle-bracket headers are the compiler's responsibility.
}
static std::string targetLayer(const std::string& target) {
    static const std::map<std::string, std::string> layers{
        {"quantum_compile_options", "Build"},
        {"quantum_core", "Core"}, {"quantum_hermite_geometry", "Core"},
        {"quantum_linalg", "Solvers"}, {"quantum_solvers", "Solvers"},
        {"quantum_interface", "Interface"}, {"lamp_interface_compatibility", "InterfaceCompatibility"}, {"quantum_electricity", "Electricity"},
        {"quantum_hydrodynamics", "Hydrodynamics"}, {"lamp_hydrodynamics_cch", "HydrodynamicsCch"}, {"quantum_plasma", "Plasma"},
        {"quantum_legacy_electricity", "Legacy"}, {"quantum_electrospray", "Apps"}, {"LAMP", "Apps"}
    };
    const auto it = layers.find(target);
    return it == layers.end() ? "Unknown" : it->second;
}
static bool allowedLink(const std::string& from, const std::string& to) {
    // Compiler settings cannot become a back door to a production library.
    if (from == "quantum_compile_options") return false;
    if (to == "quantum_compile_options") return targetLayer(from) != "Unknown";
    if (from == "LAMP" && to == "quantum_legacy_electricity") return true;
    if (from == "quantum_hydrodynamics" && to == "lamp_interface_compatibility") return true;
    if (from == "lamp_interface_compatibility" && to == "quantum_linalg") return true;
    if (from == "quantum_legacy_electricity")
        return to == "quantum_electricity" || to == "quantum_plasma" || to == "quantum_hydrodynamics";
    // Interface uses linalg for geometry, without admitting generic Solvers.
    if (from == "quantum_interface" && to == "quantum_linalg") return true;
    return allowed(targetLayer(from), targetLayer(to));
}
static bool owns(const std::string& target, const std::string& file) {
    const auto a = targetLayer(target), b = layer(file);
    return a != "Unknown" && (a == b || (a == "Apps" && b == "Runtime"));
}
static void checkGraph(const fs::path& root, const fs::path& graph, const std::set<std::string>& units) {
    std::istringstream lines(read(graph));
    std::set<std::string> targets;
    std::map<std::string, int> owners;
    for (std::string line; std::getline(lines, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto first = line.find('|'), second = line.find('|', first + 1);
        require(first != std::string::npos, "Malformed build graph: " + line);
        const auto kind = line.substr(0, first), target = line.substr(first + 1, second - first - 1);
        require(targetLayer(target) != "Unknown", "Unclassified production target: " + target);
        if (kind == "target") { targets.insert(target); continue; }
        require(second != std::string::npos, "Missing graph value: " + line);
        const auto value = line.substr(second + 1);
        if (value.empty()) continue;
        require(value.find("$<") == std::string::npos, "Unevaluated build dependency: " + line);
        if (kind == "link") {
            if (targetLayer(value) != "Unknown" || starts(value, "quantum_"))
                require(allowedLink(target, value), "Forbidden target link: " + target + " -> " + value);
        } else if (kind == "source") {
            const auto file = projectRelative(value, root);
            require(starts(file, "src/") && owns(target, file), "Wrong source owner: " + target + " -> " + file);
            if (translationUnit(file)) ++owners[file];
        } else if (kind == "include") {
            // Windows cannot express a relative path across drive roots.
            // Classify that canonical path as external before applying layers.
            const auto absolute = fs::weakly_canonical(value);
#ifdef _WIN32
            auto includeDrive = absolute.root_name().generic_string();
            auto projectDrive = root.root_name().generic_string();
            const auto lower = [](unsigned char c) { return static_cast<char>(std::tolower(c)); };
            std::transform(includeDrive.begin(), includeDrive.end(), includeDrive.begin(), lower);
            std::transform(projectDrive.begin(), projectDrive.end(), projectDrive.begin(), lower);
            if (includeDrive != projectDrive) continue;
#endif
            const auto dir = absolute.lexically_relative(root).generic_string();
            if (starts(dir, "../")) continue;
            require(targetLayer(target) != "Build", "Build settings contain a project include directory: " + dir);
            if (dir == "include") continue;
            const auto to = layer(dir + "/__header.h");
            require(allowed(targetLayer(target), to) || (targetLayer(target) == "Legacy" && to == "Legacy"),
                "Forbidden include directory: " + target + " -> " + dir);
        } else require(false, "Unknown graph record: " + line);
    }
    require(targets.size() == 14, "Expected the 14 classified production/build targets");
    for (const auto& file : units) require(owners[file] == 1, "Source needs exactly one production owner: " + file);
    require(owners.size() == units.size(), "Build graph contains missing or unscanned translation units");
    std::cout << "Build graph: " << targets.size() << " targets, " << units.size() << " translation units\n";
}
static void check(const fs::path& root, const fs::path& graph) {
    std::set<std::string> files, units;
    Edges usedCompatibility;
    size_t edgeCount = 0;
    for (const auto* base : {"src", "include"}) {
        require(fs::is_directory(root / base), std::string("Missing source root: ") + base);
        for (const auto& entry : fs::recursive_directory_iterator(root / base)) {
            require(!entry.is_symlink(), "Source symlinks are not supported: " + entry.path().string());
            if (!entry.is_regular_file() || !sourceFile(entry.path())) continue;
            const auto file = projectRelative(entry.path(), root);
            require(starts(file, "src/") ? translationUnit(file) : !translationUnit(file),
                "Headers belong in include/ and implementations in src/: " + file);
            require(layer(file) != "Unknown", "Unclassified source layout: " + file);
            for (const auto* directory : {"/Utils/", "/Common/", "/Manager/"})
                require(file.find(directory) == std::string::npos, "Unclassified utility directory: " + file);
            files.insert(file);
            if (starts(file, "src/") && translationUnit(file)) units.insert(file);
        }
    }
    for (const auto& file : files) {
        try {
            for (const auto& inc : includes(read(root / file))) {
                const auto to = resolve(root, file, inc);
                if (to.empty()) continue;
                require(allowedInclude(file, to), "Forbidden include: " + file + " -> " + to);
                if (compatibility.count({file, to})) usedCompatibility.insert({file, to});
                ++edgeCount;
            }
        } catch (const std::exception& e) { throw std::runtime_error(file + ": " + e.what()); }
    }
    for (const auto& edge : compatibility)
        require(usedCompatibility.count(edge), "Stale compatibility allowance: " + edge.first + " -> " + edge.second);
    for (const auto& file : legacyFiles) require(files.count(file), "Stale legacy file allowance: " + file);
    checkGraph(root, graph, units);
    std::cout << "Source layout/dependencies: " << files.size() << " files, " << edgeCount
        << " includes, " << usedCompatibility.size()
        << " exact legacy edges; no module -> App or Runtime -> concrete App allowance\n";
}
static void selfTest() {
    size_t cases = 0;
    const auto test = [&](bool ok, const char* name) { ++cases; require(ok, std::string("Self-test failed: ") + name); };
    const auto scan = includes("/*\n#include \"bad.h\"\n*/\n// ignored\\\n#include \"also_bad.h\"\n"
        "const char* s = R\"tag(\n#include \"raw_bad.h\"\n)tag\";\n"
        "const char* c = \"/* not a comment */\";\n"
        "# /* gap */ include \\\n \"Core/Base/Base.h\" // tail\n"
        "#if 0\n#include <Modules/Plasma/PicStepSession.h>\n#endif\n");
    test(scan.size() == 2 && scan[0].quoted && !scan[1].quoted && scan[0].name == "Core/Base/Base.h", "comments/raw strings/splices/inactive branch");
    const auto windows = includes("\xEF\xBB\xBF#include \"Core/Base/Base.h\"\r\nconst auto n = 1'000;\r\n#include <vector>\r\n");
    test(windows.size() == 2 && windows[1].name == "vector", "BOM/CRLF/numeric separators");
    test(includes("const auto s = u8R\"(\n#include \"fake.h\"\n)\";\n%:include <vector>\n").size() == 1,
        "prefixed raw string and preprocessor digraph");
    for (const auto* bad : {"#include HEADER\n", "#include_next <x>\n", "#import <x>\n"}) {
        bool rejected = false;
        try { includes(bad); } catch (const std::exception&) { rejected = true; }
        test(rejected, "unsupported directive rejected");
    }
    for (const auto* from : {"Core", "Solvers", "Hydrodynamics", "HydrodynamicsCch", "Electricity", "Interface", "Plasma", "Legacy"}) {
        test(!allowed(from, "Apps"), "reverse App dependency");
        test(!allowed(from, "Runtime"), "Runtime alias dependency");
    }
    test(allowed("Hydrodynamics", "Interface"), "Hydro public Interface dependency");
    test(!allowed("Interface", "Hydrodynamics"), "Interface independent of Hydro");
    test(!allowed("Electricity", "Plasma") && !allowed("Electricity", "Hydrodynamics"), "formal Electricity independence");
    test(!allowed("Solvers", "Plasma") && !allowed("Core", "Solvers"), "lower-layer direction");
    test(!allowed("Runtime", "Apps"), "Runtime cannot include concrete App");
    for (const auto* module : {"Hydrodynamics", "Electricity", "Interface", "Plasma"})
        test(!allowed("Runtime", module), "Runtime cannot include physics module");
    for (const auto& edge : compatibility) test(allowedInclude(edge.first, edge.second), "live exact compatibility edge");
    for (const auto* app : {"Electrospray/ElectrosprayApp.h", "Electrohydrodynamics/ElectrohydrodynamicsApp.h",
        "SingleModule/SingleModuleApp.h", "Interface/InterfaceApp.h", "Electricity/ElectricityApp.h", "Verification/Pic/PicVerificationApp.h"})
        test(!allowedInclude("src/Apps/Runtime/RunController.cpp", std::string("include/Apps/") + app), "retired registration exemption stays forbidden");
    test(!allowedInclude("src/Apps/Other.cpp", "include/Modules/Electricity/Legacy/Electricity.h"), "no new legacy consumer");
    test(!allowedInclude("include/Modules/Electricity/Legacy/Electricity.h", "include/Modules/Plasma/PicStepSession.h"), "no broad legacy exemption");
    test(layer("include/Modules/Electricity/Legacy/New.h") == "Unknown", "no new legacy file");
    test(layer("src/Unknown/New.cpp") == "Unknown" && layer("src/Modules/New/New.cpp") == "Unknown", "unknown source layout");
    test(!owns("quantum_plasma", "src/Apps/Runtime/RunController.cpp"), "no App source in module target");
    test(!owns("quantum_electricity", "src/Modules/Electricity/Legacy/Electricity.cpp"), "formal target excludes legacy source");
    test(owns("quantum_legacy_electricity", "src/Modules/Electricity/Legacy/Electricity.cpp"), "legacy source has explicit owner");
    test(!allowedLink("quantum_plasma", "LAMP") && !allowedLink("quantum_electricity", "quantum_hydrodynamics"), "forbidden build edges");
    test(allowedLink("quantum_interface", "quantum_linalg"), "Interface linear algebra dependency");
    test(!allowedLink("quantum_interface", "lamp_interface_compatibility"), "GALS cannot link compatibility");
    test(!allowedInclude("src/Modules/Hydrodynamics/New.cpp", "include/Modules/Interface/Compatibility/InterfaceState.h"), "new compatibility consumer rejected");
    test(!allowedInclude("src/Modules/Interface/HermiteGalsFirstJetStep.cpp", "include/Modules/Interface/Compatibility/InterfaceState.h"), "FirstJet compatibility dependency rejected");
    test(allowedLink("lamp_hydrodynamics_cch", "quantum_core"), "CCH Core dependency");
    test(!allowedLink("lamp_hydrodynamics_cch", "quantum_hydrodynamics"), "CCH cannot link frozen Hydro");
    test(!allowedLink("lamp_hydrodynamics_cch", "lamp_interface_compatibility"), "CCH cannot link compatibility");
    test(!allowedLink("lamp_hydrodynamics_cch", "LAMP"), "CCH cannot link App");
    test(owns("lamp_hydrodynamics_cch", "src/Modules/Hydrodynamics/CchProbe.cpp"), "native CCH source owner");
    test(!owns("quantum_hydrodynamics", "src/Modules/Hydrodynamics/CchProbe.cpp"), "frozen Hydro cannot compile CCH");
    test(!owns("lamp_hydrodynamics_cch", "src/Modules/Hydrodynamics/HydroProbe.cpp"), "CCH cannot compile frozen Hydro");
    test(!allowedInclude("include/Modules/Hydrodynamics/CchProbe.h", "include/Modules/Hydrodynamics/HydroInterfaceBoundary.h"), "CCH cannot transitively consume legacy Hydro headers");
    std::cout << "Architecture scanner/policy self-tests: " << cases << " passed\n";
}
static void fixtureTests(const fs::path& scratch) {
    // A fresh, tiny synthetic repository exercises resolution and the actual
    // graph reader without editing, copying or compiling production sources.
    const auto root = fs::absolute(scratch / std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::map<std::string, std::string> files;
    for (const auto& file : legacyFiles) files[file] = "";
    for (const auto& edge : compatibility) {
        files[edge.first] += "#include \"" + edge.second.substr(8) + "\"\n";
        files.try_emplace(edge.second, "");
    }
    files["include/Apps/Electricity/ElectricityApp.h"] = "";
    files["src/Apps/Runtime/AppRegistry.cpp"] = "";
    files["include/Core/Base/Base.h"] = "";
    files["include/Apps/Runtime/RunConfig.h"] = "";
    files["src/Core/Base/Probe.cpp"] = "#include \"Core/Base/Base.h\"\n";
    files["src/Modules/Plasma/Probe.cpp"] = "#include <Core/Base/Base.h>\n";
    files["src/Modules/Hydrodynamics/CchProbe.cpp"] = "#include \"Core/Base/Base.h\"\n";
    const auto write = [&](const std::string& path, const std::string& text) {
        fs::create_directories((root / path).parent_path());
        std::ofstream out(root / path, std::ios::binary);
        out << text;
        require(bool(out), "Cannot write fixture " + path);
    };
    for (const auto& file : files) write(file.first, file.second);
    std::string graph;
    for (const auto* target : {"quantum_compile_options", "quantum_core", "quantum_hermite_geometry", "quantum_linalg", "quantum_solvers",
        "quantum_interface", "lamp_interface_compatibility", "quantum_electricity", "quantum_hydrodynamics", "lamp_hydrodynamics_cch", "quantum_plasma",
        "quantum_legacy_electricity", "quantum_electrospray", "LAMP"}) graph += std::string("target|") + target + '\n';
    for (const auto& file : files) {
        if (!starts(file.first, "src/") || !translationUnit(file.first)) continue;
        const auto group = layer(file.first);
        const std::string owner = group == "Core" ? "quantum_core" : group == "Plasma" ? "quantum_plasma" :
            group == "Legacy" ? "quantum_legacy_electricity" : group == "InterfaceCompatibility" ? "lamp_interface_compatibility" : group == "HydrodynamicsCch" ? "lamp_hydrodynamics_cch" : group == "Hydrodynamics" ? "quantum_hydrodynamics" : "LAMP";
        graph += "source|" + owner + '|' + (root / file.first).generic_string() + '\n';
    }
    write("graph.txt", graph);
    check(fs::canonical(root), root / "graph.txt");
#ifdef _WIN32
    const auto externalDrive = root.root_name().generic_string() == "C:" ? "D:" : "C:";
    write("graph.txt", graph + "include|quantum_core|" + externalDrive + "/LAMP_external_include\n");
    check(fs::canonical(root), root / "graph.txt");
    write("graph.txt", graph);
    std::cout << "Architecture cross-drive external include fixture: passed\n";
#endif
    size_t rejected = 0;
    const auto reject = [&](const std::string& path, const std::string& value, const std::string& expected) {
        const bool existed = fs::exists(root / path);
        const auto before = existed ? read(root / path) : std::string();
        write(path, value);
        std::string error;
        try { check(fs::canonical(root), root / "graph.txt"); } catch (const std::exception& e) { error = e.what(); }
        if (existed) write(path, before); else fs::remove(root / path);
        require(error.find(expected) != std::string::npos, "Mutation was not rejected as expected: " + path + ": " + error);
        ++rejected;
    };
    reject("src/Modules/Plasma/Probe.cpp", "#include \"Runtime/RunConfig.h\"\r\n", "Forbidden include");
    reject("src/Modules/Plasma/Probe.cpp", "#include \"../../../include/Apps/Runtime/RunConfig.h\"\n", "Forbidden include");
    reject("src/Core/Base/Probe.cpp", "#if 0\n#include <Modules/Plasma/Particles/ParticleType.h>\n#endif\n", "Forbidden include");
    reject("src/Core/Base/Probe.cpp", "#include HEADER\n", "literal header");
    reject("src/Core/Base/Probe.cpp", "#include \"Core/Missing.h\"\n", "unresolved project header");
    reject("src/Apps/New.cpp", "#include \"Modules/Electricity/Legacy/Electricity.h\"\n", "Forbidden include");
    reject("src/Apps/Runtime/New.cpp", "#include \"Electricity/ElectricityApp.h\"\n", "Forbidden include");
    reject("src/Apps/Runtime/New.cpp", "#include \"Modules/Plasma/Particles/ParticleType.h\"\n", "Forbidden include");
    reject("src/Modules/Electricity/Legacy/New.cpp", "", "Unclassified source layout");
    reject("src/Modules/New/New.cpp", "", "Unclassified source layout");
    reject("src/Core/Base/Unbuilt.cpp", "", "exactly one production owner");
    reject("include/Core/Utils/New.h", "", "Unclassified utility directory");
    reject("src/Core/Base/New.h", "", "Headers belong in include/");
    reject("include/Core/Base/New.cpp", "", "Headers belong in include/");
    reject("src/Apps/Electricity/Compatibility/LegacyElectricityModule.cpp", "", "Stale compatibility allowance");
    reject("graph.txt", graph + "link|quantum_electricity|quantum_hydrodynamics\n", "Forbidden target link");
    reject("graph.txt", graph + "link|quantum_compile_options|LAMP\n", "Forbidden target link");
    reject("graph.txt", graph + "link|quantum_compile_options|quantum_core\n", "Forbidden target link");
    reject("graph.txt", graph + "include|quantum_compile_options|" + (root / "include").generic_string() + '\n', "Build settings contain");
    reject("graph.txt", graph + "source|quantum_compile_options|" + (root / "src/Core/Base/Probe.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "include|quantum_plasma|" + (root / "include/Apps").generic_string() + '\n', "Forbidden include directory");
    reject("graph.txt", graph + "source|quantum_plasma|" + (root / "src/Apps/Runtime/AppRegistry.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "source|quantum_core|" + (root / "src/Core/Base/Probe.cpp").generic_string() + '\n', "exactly one production owner");
    reject("graph.txt", graph + "link|lamp_hydrodynamics_cch|quantum_hydrodynamics\n", "Forbidden target link");
    reject("graph.txt", graph + "link|lamp_hydrodynamics_cch|lamp_interface_compatibility\n", "Forbidden target link");
    reject("graph.txt", graph + "link|lamp_hydrodynamics_cch|LAMP\n", "Forbidden target link");
    reject("graph.txt", graph + "source|quantum_hydrodynamics|" + (root / "src/Modules/Hydrodynamics/CchProbe.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "source|lamp_hydrodynamics_cch|" + (root / "src/Core/Base/Probe.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "source|lamp_hydrodynamics_cch|" + (root / "src/Modules/Hydrodynamics/HydroInterfaceBoundary.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "source|lamp_hydrodynamics_cch|" + (root / "src/Modules/Hydrodynamics/CchProbe.cpp").generic_string() + '\n', "exactly one production owner");
    reject("src/Modules/Hydrodynamics/CchProbe.cpp", "#include \"Modules/Hydrodynamics/HydroInterfaceBoundary.h\"\n", "Forbidden include");
    std::cout << "Architecture end-to-end fixtures: baseline passed, " << rejected << " mutations rejected\n";
    // Remove only the exact files created above. Empty fixture directories stay
    // in the disposable check tree; no recursive deletion is needed.
    for (const auto& file : files) fs::remove(root / file.first);
    fs::remove(root / "graph.txt");
}
int main(int argc, char** argv) {
    try {
        require(argc == 3, "Usage: architecture_check <project root> <generated graph> | --self-test <fixture directory>");
        if (std::string(argv[1]) == "--self-test") { selfTest(); fixtureTests(argv[2]); }
        else check(fs::canonical(argv[1]), argv[2]);
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
