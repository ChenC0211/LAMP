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
// The removed SEVOF route cannot return as unbuilt reference sources.
static bool retiredSeVof(const std::string& file) {
    auto name = file.substr(file.find_last_of('/') + 1);
    std::transform(name.begin(), name.end(), name.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return name.find("sevof") != std::string::npos;
}
// Retained HG sources are reference assets, never production owners.
// Match the same frozen retirement scope as the production CMake source lists.
static bool referenceOnly(const std::string& file) {
    const auto p = file.substr(file.find('/') + 1);
    if (starts(p, "Apps/Interface/")) return true;
    // Exact new VOF reader, not a directory-wide revival of the old HG App.
    if (starts(p, "Apps/Electrohydrodynamics/"))
        return p != "Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"
            && p != "Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.cpp"
            && p != "Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"
            && p != "Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.cpp"
            && p != "Apps/Electrohydrodynamics/VofElectrohydrodynamicsApp.h"
            && p != "Apps/Electrohydrodynamics/VofElectrohydrodynamicsApp.cpp"
            && p != "Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.h"
            && p != "Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp";
    return false;
}
static bool translationUnit(const fs::path& p) {
    static const std::set<std::string> extensions{ ".c", ".cc", ".cpp", ".cxx" };
    return extensions.count(p.extension().string()) != 0;
}
// P1 consumes the authoritative immutable publication and same-source phase
// cubature only through these exact adapter files. These are current typed data
// dependencies, not legacy or directory-wide exemptions. Removed consumers
// must remove their allowance.
static const Edges publicationViews{
    {"include/Modules/Electric/ElectricVofSharpField.h", "include/Modules/Interface/VofSmoothGeometry.h"},
    {"include/Modules/Electric/ElectricSurfaceChargeTransport.h", "include/Modules/Interface/GeometricVofTransport.h"},
    {"include/Modules/Hydrodynamics/HydroMaxwellLoad.h", "include/Modules/Electric/ElectricVofSharpMaxwellStress.h"},
    {"include/Modules/Electric/ElectricPhaseGeometry.h", "include/Modules/Interface/HermiteGalsPublication.h"},
    {"src/Modules/Electric/ElectricPhaseGeometry.cpp", "include/Modules/Interface/HermiteGalsPhaseQuadrature.h"},
    {"include/Modules/Electric/ElectricVofField.h", "include/Modules/Interface/GeometricVofPublication.h"},
    {"include/Modules/Hydrodynamics/HydroMaxwellLoad.h", "include/Modules/Electric/ElectricVofMaxwellStress.h"},
    {"include/Modules/Electric/ElectricVofChargeTransport.h", "include/Modules/Interface/GeometricVofTransport.h"},
    {"include/Modules/Electric/ElectricVofChargeTransport.h", "include/Modules/Interface/GeometricVofMeshTransfer.h"},
};
static std::string layer(const std::string& file) {
    const auto slash = file.find('/');
    if (slash == std::string::npos || (file.substr(0, slash) != "src" && file.substr(0, slash) != "include")) return "Unknown";
    const auto p = file.substr(slash + 1);
    if (starts(p, "Core/")) return "Core";
    if (starts(p, "Solvers/")) return "Solvers";
    // Retired solver/adapter directories are not classified production sources.
    if (starts(p, "Modules/Electricity/") || starts(p, "Apps/Electricity/")) return "Unknown";
    if (starts(p, "Modules/Electric/Legacy/") || starts(p, "Apps/Electric/Compatibility/")
        || starts(p, "Modules/Hydrodynamics/LegacyNodal")
        || starts(p, "Modules/Hydrodynamics/IO/LegacyNodal")) return "Unknown";
    if (starts(p, "Apps/Runtime/") || starts(p, "Apps/Coupling/") || starts(p, "Apps/Verification/")
        || starts(p, "Apps/SingleModule/") || starts(p, "Apps/Electrospray/")) return "Unknown";
    if (starts(p, "Apps/IO/TwoPhaseInputComponents.") || starts(p, "Apps/IO/HydroOutputComponents.")) return "AppFluidIO";
    for (const auto* app : {"Electric", "Electrohydrodynamics", "Hydrodynamics", "Interface", "Pic"})
        if (starts(p, std::string("Apps/") + app + '/')) return std::string(app) + "App";
    if (p == "Apps/RunCase.cpp" || p == "Apps/RunCase.h" || p == "Apps/lamp_main.cpp") return "Apps";
    // Internal policy label retained for module-neutral entry/IO facilities;
    // it does not authorize a Runtime directory or concrete App dependency.
    if (starts(p, "Apps/")) return "Runtime";
    if (starts(p, "Modules/Interface/Compatibility/")) return "InterfaceCompatibility";
    if (starts(p, "Modules/Hydrodynamics/Hydro") || starts(p, "Modules/Hydrodynamics/Detail/Hydro"))
        return "HydrodynamicsHydro";
    for (const auto* m : {"Hydrodynamics", "Electric", "Interface", "Plasma"})
        if (starts(p, std::string("Modules/") + m + '/')) return m;
    return "Unknown";
}
static bool allowed(const std::string& from, const std::string& to) {
    if (from == "Unknown" || to == "Unknown") return false;
    if (from == to || to == "Core") return true;
    if (to == "InterfaceCompatibility") return false;
    if (from == "InterfaceCompatibility") return to == "Interface";
    if (from == "Apps") return true;
    if (from == "Runtime") return to == "Solvers";
    if (from == "AppFluidIO") return to == "Runtime" || to == "HydrodynamicsHydro" || to == "Interface" || to == "Solvers";
    if (from == "HydrodynamicsApp") return to == "Runtime" || to == "AppFluidIO" || to == "HydrodynamicsHydro" || to == "Interface" || to == "Solvers";
    if (from == "ElectrohydrodynamicsApp") return to == "Runtime" || to == "AppFluidIO" || to == "HydrodynamicsApp"
        || to == "HydrodynamicsHydro" || to == "Interface" || to == "Electric" || to == "Solvers";
    if (from == "ElectricApp") return to == "Runtime" || to == "Electric" || to == "Interface" || to == "Solvers";
    if (from == "InterfaceApp") return to == "Runtime" || to == "Interface" || to == "Solvers";
    if (from == "PicApp") return to == "Runtime" || to == "Plasma" || to == "Solvers";
    if (from == "Hydrodynamics") return to == "Interface" || to == "Solvers";
    if (from == "HydrodynamicsHydro") return to == "Interface" || to == "Solvers";
    if (from == "Electric" || from == "Plasma") return to == "Solvers";
    return false;
}
static bool allowedInclude(const std::string& from, const std::string& to) {
    if (retiredSeVof(from) || retiredSeVof(to)) return false;
    if (referenceOnly(from) || referenceOnly(to)) return false;
    // The two shared fluid serializers borrow the pure typed input, never an App controller.
    const bool fluidInput = (from == "include/Apps/IO/TwoPhaseInputComponents.h" || from == "include/Apps/IO/HydroOutputComponents.h")
        && to == "include/Apps/Hydrodynamics/IO/HydrodynamicsInput.h";
    return publicationViews.count({from, to}) || fluidInput || allowed(layer(from), layer(to));
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
        {"lamp_compile_options", "Build"},
        {"lamp_core", "Core"}, {"lamp_hermite_geometry", "Core"},
        {"lamp_linalg", "Solvers"}, {"lamp_solvers", "Solvers"},
        {"lamp_interface", "Interface"}, {"lamp_electric", "Electric"},
        {"lamp_hydrodynamics", "HydrodynamicsHydro"}, {"lamp_plasma", "Plasma"},
        {"LAMP", "Apps"}
    };
    const auto it = layers.find(target);
    return it == layers.end() ? "Unknown" : it->second;
}
static bool allowedLink(const std::string& from, const std::string& to) {
    // Compiler settings cannot become a back door to a production library.
    if (from == "lamp_compile_options") return false;
    if (to == "lamp_compile_options") return targetLayer(from) != "Unknown";
    // Interface uses linalg for geometry, without admitting generic Solvers.
    if (from == "lamp_interface" && to == "lamp_linalg") return true;
    // PUBLIC linkage supplies the exact publication API and its implementation;
    // include-directory policy still rejects Electric -> Interface generally.
    if (from == "lamp_electric" && to == "lamp_interface") return true;
    // This one native load consumes the exact readonly Maxwell stress owner.
    if (from == "lamp_hydrodynamics" && to == "lamp_electric") return true;
    return allowed(targetLayer(from), targetLayer(to));
}
static bool owns(const std::string& target, const std::string& file) {
    const auto a = targetLayer(target), b = layer(file);
    return a != "Unknown" && (a == b || (a == "Apps" && (b == "Runtime" || b == "AppFluidIO" || (b.size() >= 3 && b.compare(b.size()-3,3,"App") == 0))));
}
static void checkGraph(const fs::path& root, const fs::path& graph, const std::set<std::string>& units) {
    std::istringstream lines(read(graph));
    std::set<std::string> targets;
    std::map<std::string, int> owners;
    bool maxwellImplementationLink=false;
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
            if (target=="lamp_hydrodynamics" && value=="lamp_electric") maxwellImplementationLink=true;
            if (targetLayer(value) != "Unknown" || starts(value, "lamp_") || starts(value, "lamp_"))
                require(allowedLink(target, value), "Forbidden target link: " + target + " -> " + value);
        } else if (kind == "source") {
            const auto file = projectRelative(value, root);
            require(!retiredSeVof(file), "Retired SEVOF source in production graph: " + file);
            require(!referenceOnly(file), "Reference-only source in production graph: " + file);
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
            require(allowed(targetLayer(target), to),
                "Forbidden include directory: " + target + " -> " + dir);
        } else require(false, "Unknown graph record: " + line);
    }
    require(targets.size() == 10, "Expected the 10 classified production/build targets");
    require(maxwellImplementationLink,"Missing exact Maxwell implementation link");
    for (const auto& file : units) require(owners[file] == 1, "Source needs exactly one production owner: " + file);
    require(owners.size() == units.size(), "Build graph contains missing or unscanned translation units");
    std::cout << "Build graph: " << targets.size() << " targets, " << units.size() << " translation units\n";
}
static void check(const fs::path& root, const fs::path& graph) {
    std::set<std::string> files, units;
    Edges usedPublicationViews;
    size_t edgeCount = 0;
    for (const auto* base : {"src", "include"}) {
        require(fs::is_directory(root / base), std::string("Missing source root: ") + base);
        for (const auto& entry : fs::recursive_directory_iterator(root / base)) {
            require(!entry.is_symlink(), "Source symlinks are not supported: " + entry.path().string());
            if (!entry.is_regular_file() || !sourceFile(entry.path())) continue;
            const auto file = projectRelative(entry.path(), root);
            require(!retiredSeVof(file), "Retired SEVOF source must remain removed: " + file);
            require(starts(file, "src/") ? translationUnit(file) : !translationUnit(file),
                "Headers belong in include/ and implementations in src/: " + file);
            if (referenceOnly(file)) continue;
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
                if (publicationViews.count({file, to})) usedPublicationViews.insert({file, to});
                ++edgeCount;
            }
        } catch (const std::exception& e) { throw std::runtime_error(file + ": " + e.what()); }
    }
    for (const auto& edge : publicationViews)
        require(usedPublicationViews.count(edge), "Stale publication view allowance: " + edge.first + " -> " + edge.second);
    checkGraph(root, graph, units);
    std::cout << "Source layout/dependencies: " << files.size() << " files, " << edgeCount
        << " includes, " << usedPublicationViews.size()
        << " exact immutable publication edges; no module -> App or Runtime -> concrete App allowance\n";
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
    for (const auto* from : {"Core", "Solvers", "Hydrodynamics", "HydrodynamicsHydro", "Electric", "Interface", "Plasma"}) {
        test(!allowed(from, "Apps"), "reverse App dependency");
        test(!allowed(from, "Runtime"), "Runtime alias dependency");
    }
    test(allowed("Hydrodynamics", "Interface"), "Hydro public Interface dependency");
    test(!allowed("Interface", "Hydrodynamics"), "Interface independent of Hydro");
    test(!allowed("Electric", "Plasma") && !allowed("Electric", "Hydrodynamics"), "formal Electric independence");
    test(!allowed("Solvers", "Plasma") && !allowed("Core", "Solvers"), "lower-layer direction");
    test(!allowed("Runtime", "Apps"), "Runtime cannot include concrete App");
    for (const auto* module : {"Hydrodynamics", "Electric", "Interface", "Plasma"})
        test(!allowed("Runtime", module), "Runtime cannot include physics module");
    test(allowedInclude("include/Modules/Electric/ElectricPhaseGeometry.h",
        "include/Modules/Interface/HermiteGalsPublication.h"), "exact electric immutable publication dependency");
    test(allowedInclude("src/Modules/Electric/ElectricPhaseGeometry.cpp",
        "include/Modules/Interface/HermiteGalsPhaseQuadrature.h"), "exact electric same-source phase cubature dependency");
    test(!allowedInclude("src/Modules/Electric/Other.cpp",
        "include/Modules/Interface/HermiteGalsPhaseQuadrature.h"), "other Electric phase cubature consumer rejected");
    test(!allowedInclude("src/Modules/Electric/ElectricPhaseGeometry.cpp",
        "include/Modules/Interface/HermiteGalsState.h"), "phase cubature adapter cannot directly include primary state");
    // E1/E2 borrow only these actual immutable VOF data contracts.
    const Edges vofViews{
    {"include/Modules/Electric/ElectricVofField.h", "include/Modules/Interface/GeometricVofPublication.h"},
    {"include/Modules/Electric/ElectricVofChargeTransport.h", "include/Modules/Interface/GeometricVofTransport.h"},
    {"include/Modules/Electric/ElectricVofChargeTransport.h", "include/Modules/Interface/GeometricVofMeshTransfer.h"},
    };
    for (const auto& edge : vofViews) {
        test(allowedInclude(edge.first, edge.second), "exact VOF Electric publication/transport dependency");
        test(!allowedInclude("include/Modules/Electric/Other.h", edge.second),
            "VOF contracts do not authorize another Electric consumer");
        test(!allowedInclude(edge.second, edge.first), "no reverse VOF Interface Electric dependency");
        test(!allowedInclude(edge.first, "include/Modules/Interface/Detail/GeometricVofPrivate.h"),
            "exact VOF consumers cannot include private Interface helpers");
    }
    test(!allowedInclude("include/Modules/Electric/ElectricVofField.h",
        "include/Modules/Interface/GeometricVofTransport.h"), "E1 field cannot inherit E2 transport allowance");
    test(!allowedInclude("include/Modules/Electric/ElectricVofChargeTransport.h",
        "include/Modules/Interface/GeometricVofPublication.h"), "E2 transport cannot inherit E1 publication allowance");
    const Edges b1Views{
        {"include/Modules/Electric/ElectricVofSharpField.h", "include/Modules/Interface/VofSmoothGeometry.h"},
        {"include/Modules/Electric/ElectricSurfaceChargeTransport.h", "include/Modules/Interface/GeometricVofTransport.h"},
        {"include/Modules/Hydrodynamics/HydroMaxwellLoad.h", "include/Modules/Electric/ElectricVofSharpMaxwellStress.h"},
    };
    for (const auto& edge : b1Views) {
        const bool hydro = edge.first == "include/Modules/Hydrodynamics/HydroMaxwellLoad.h";
        test(allowedInclude(edge.first, edge.second), "exact immutable B1 publication dependency");
        test(!allowedInclude(hydro ? "include/Modules/Hydrodynamics/HydroOther.h" : "include/Modules/Electric/Other.h", edge.second),
            "B1 allowance cannot authorize another consumer");
        test(!allowedInclude(edge.second, edge.first), "B1 owner cannot depend on its consumer");
        test(!allowedInclude(edge.first, hydro ? "include/Modules/Electric/Detail/ElectricVofPrivate.h"
            : "include/Modules/Interface/Detail/GeometricVofPrivate.h"), "B1 allowance cannot authorize private helpers");
    }
    test(!allowedInclude("include/Modules/Electric/ElectricVofSharpField.h",
        "include/Modules/Interface/GeometricVofTransport.h"), "B1 field cannot inherit material transport allowance");
    test(!allowedInclude("include/Modules/Electric/ElectricSurfaceChargeTransport.h",
        "include/Modules/Interface/VofSmoothGeometry.h"), "B1 material transport cannot inherit field geometry allowance");
    test(!allowedInclude("include/Modules/Electric/ElectricVofField.h",
        "include/Modules/Interface/VofSmoothGeometry.h"), "A field cannot inherit B1 geometry allowance");
    test(!allowed("Electric", "Interface"), "no Electric Interface directory allowance");
    test(!allowedInclude("include/Modules/Electric/Other.h",
        "include/Modules/Interface/HermiteGalsPublication.h"), "other Electric publication consumer rejected");
    test(!allowedInclude("include/Modules/Electric/ElectricPhaseGeometry.h",
        "include/Modules/Interface/HermiteGalsState.h"), "electric cannot directly include primary state");
    test(!allowedInclude("include/Modules/Interface/HermiteGalsPublication.h",
        "include/Modules/Electric/ElectricPhaseGeometry.h"), "Interface cannot depend on electric adapter");
    test(allowedLink("lamp_electric", "lamp_interface"), "exact electric publication implementation link");
    test(!allowedLink("lamp_interface", "lamp_electric"), "no reverse Interface Electric link");
    test(allowedInclude("include/Modules/Hydrodynamics/HydroMaxwellLoad.h",
        "include/Modules/Electric/ElectricVofMaxwellStress.h"), "exact readonly Hydro Maxwell stress dependency");
    test(!allowed("HydrodynamicsHydro","Electric"),"no generic Hydro Electric directory allowance");
    test(!allowedInclude("include/Modules/Hydrodynamics/HydroOther.h",
        "include/Modules/Electric/ElectricVofMaxwellStress.h"),"other Hydro Maxwell consumer rejected");
    test(!allowedInclude("include/Modules/Electric/ElectricVofMaxwellStress.h",
        "include/Modules/Hydrodynamics/HydroMaxwellLoad.h"),"no reverse Maxwell Hydro dependency");
    test(!allowedInclude("include/Modules/Hydrodynamics/HydroMaxwellLoad.h",
        "include/Modules/Electric/ElectricVofField.h"),"Maxwell load cannot include a different electric owner directly");
    test(!allowedInclude("include/Modules/Hydrodynamics/HydroMaxwellLoad.h",
        "include/Modules/Electric/Detail/ElectricVofPrivate.h"),"Maxwell load cannot include private Electric helper");
    test(allowedLink("lamp_hydrodynamics","lamp_electric"),"exact readonly Maxwell implementation link");
    test(!allowedLink("lamp_electric","lamp_hydrodynamics"),"Electric remains independent of Hydro");

    for (const auto* app : {"Electrospray/ElectrosprayApp.h", "Electrohydrodynamics/ElectrohydrodynamicsApp.h",
        "Pic/PicApp.h", "Interface/InterfaceApp.h", "Electric/ElectricApp.h", "Electrohydrodynamics/ElectrohydrodynamicsApp.h"})
        test(!allowedInclude("src/Apps/RunController.cpp", std::string("include/Apps/") + app), "Runtime cannot depend on concrete Apps");
    test(!allowedInclude("src/Apps/Other.cpp", "include/Modules/Electric/Legacy/Electric.h"), "no new legacy consumer");
    test(!allowedInclude("include/Modules/Electric/Legacy/Electric.h", "include/Modules/Plasma/PicStepSession.h"), "no broad legacy exemption");
    test(layer("include/Modules/Electric/Legacy/New.h") == "Unknown", "no new legacy file");
    test(layer("src/Unknown/New.cpp") == "Unknown" && layer("src/Modules/New/New.cpp") == "Unknown", "unknown source layout");
    test(!owns("lamp_plasma", "src/Apps/RunController.cpp"), "no App source in module target");
    test(!owns("lamp_electric", "src/Modules/Electric/Legacy/Electric.cpp"), "formal target excludes legacy source");
    test(!owns("lamp_legacy_electric", "src/Modules/Electric/Legacy/Electric.cpp"), "retired source has no production owner");
    test(!allowedLink("lamp_plasma", "LAMP") && !allowedLink("lamp_electric", "lamp_legacy_nodal_adapter"), "forbidden build edges");
    test(allowedLink("lamp_interface", "lamp_linalg"), "Interface linear algebra dependency");
    test(!allowedLink("lamp_interface", "lamp_interface_compatibility"), "GALS cannot link compatibility");
    test(!allowedInclude("src/Modules/Hydrodynamics/New.cpp", "include/Modules/Interface/Compatibility/InterfaceState.h"), "new compatibility consumer rejected");
    test(!allowedInclude("src/Modules/Interface/HermiteGalsFirstJetStep.cpp", "include/Modules/Interface/Compatibility/InterfaceState.h"), "FirstJet compatibility dependency rejected");
    test(allowedLink("lamp_hydrodynamics", "lamp_core"), "HYDRO Core dependency");
    test(!allowedLink("lamp_hydrodynamics", "lamp_legacy_nodal_adapter"), "HYDRO cannot link legacy electric array adapter");
    test(!allowedLink("lamp_hydrodynamics", "lamp_interface_compatibility"), "HYDRO cannot link compatibility");
    test(!allowedLink("lamp_hydrodynamics", "LAMP"), "HYDRO cannot link App");
    test(owns("lamp_hydrodynamics", "src/Modules/Hydrodynamics/HydroProbe.cpp"), "native HYDRO source owner");
    test(!owns("lamp_legacy_nodal_adapter", "src/Modules/Hydrodynamics/HydroProbe.cpp"), "legacy electric array adapter cannot compile HYDRO");
    test(!owns("lamp_hydrodynamics", "src/Modules/Hydrodynamics/LegacyNodalProbe.cpp"), "HYDRO cannot compile legacy electric array adapter");
    test(!allowedInclude("include/Modules/Hydrodynamics/HydroProbe.h", "include/Modules/Hydrodynamics/LegacyNodalState.h"), "HYDRO cannot transitively consume legacy electric array headers");
    // AppNormalization structural guards; original module policy tests above remain intact.
    test(!allowedInclude("src/Apps/Hydrodynamics/HydrodynamicsRun.cpp", "include/Apps/Electrohydrodynamics/ElectrohydrodynamicStep.h"), "pure Hydro rejects joint owner");
    test(!allowedInclude("src/Apps/Hydrodynamics/IO/HydrodynamicsOutput.cpp", "include/Modules/Electric/ElectricPotential.h"), "pure Hydro rejects electric module");
    test(!allowedInclude("src/Apps/Electrohydrodynamics/ElectrohydrodynamicsRun.cpp", "include/Apps/Hydrodynamics/HydrodynamicsRun.h"), "retired joint owner is not production");
    test(!allowedInclude("src/Apps/Electrohydrodynamics/IO/ElectrohydrodynamicsOutput.cpp", "include/Modules/Electric/ElectricPotential.h"), "retired joint output is not production");
    test(allowedInclude("include/Apps/IO/HydroOutputComponents.h", "include/Apps/Hydrodynamics/IO/HydrodynamicsInput.h"), "shared serializer borrows pure typed input");
    test(!allowedInclude("src/Apps/IO/HydroOutputComponents.cpp", "include/Modules/Electric/ElectricPotential.h"), "shared serializer rejects electric receipt");
    test(!allowedInclude("src/Apps/IO/TwoPhaseInputComponents.cpp", "include/Apps/Electrohydrodynamics/IO/ElectrohydrodynamicsInput.h"), "shared parser rejects joint owner");
    test(layer("src/Apps/Verification/VerificationApp.cpp") == "Unknown" && layer("src/Apps/Coupling/Step.cpp") == "Unknown"
        && layer("src/Apps/Runtime/RunCase.cpp") == "Unknown", "retired directories rejected");
    test(!owns("lamp_hydrodynamics", "src/Apps/IO/HydroOutputComponents.cpp"), "shared App serialization has no numerical module owner");
    test(!referenceOnly("src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsApp.cpp"), "exact new VOF EHD entry is production");
    test(!referenceOnly("src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp"), "exact joint VOF EHD Run is production");
    test(!referenceOnly("include/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.h"), "exact joint VOF EHD Run header is production");
    test(!referenceOnly("include/Apps/Electrohydrodynamics/VofElectrohydrodynamicsApp.h"), "exact new VOF EHD entry header is production");
    test(referenceOnly("src/Apps/Electrohydrodynamics/ElectrohydrodynamicsApp.cpp"), "old HG entry remains reference");
    test(referenceOnly("src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsFuture.cpp"), "new prefix is not a directory allowance");
    test(allowedInclude("src/Apps/RunCase.cpp", "include/Apps/Electrohydrodynamics/VofElectrohydrodynamicsApp.h"), "dispatcher consumes exact new entry");
    test(!referenceOnly("src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.cpp"), "exact VOF EHD reader is production");
    test(!referenceOnly("include/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"), "exact VOF EHD reader header is production");
    test(!referenceOnly("src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.cpp"), "exact passive VOF EHD output is production");
    test(!referenceOnly("include/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"), "exact passive VOF EHD output header is production");
    test(referenceOnly("src/Apps/Electrohydrodynamics/IO/ElectrohydrodynamicsOutput.cpp"), "retained HG lifecycle output stays reference only");
    test(referenceOnly("src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsFutureOutput.cpp"), "new output prefix is not a directory allowance");
    test(allowedInclude("src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.cpp", "include/Apps/IO/HydroOutputComponents.h"), "new output borrows shared passive fluid serializer");
    test(allowedInclude("src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.cpp", "include/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"), "new output consumes its exact header");
    test(!allowedInclude("src/Modules/Electric/ElectricVofField.cpp", "include/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.h"), "electric physical kernel cannot depend on App output");
    test(referenceOnly("src/Apps/Electrohydrodynamics/IO/ElectrohydrodynamicsInput.cpp"), "old HG EHD reader remains reference");
    test(referenceOnly("src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsFuture.cpp"), "VOF name prefix does not grant production scope");
    test(allowedInclude("src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.cpp", "include/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"), "new reader consumes its actual header");
    test(allowedInclude("include/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h", "include/Modules/Electric/ElectricVofField.h"), "new reader consumes typed VOF electric config");
    test(retiredSeVof("src/Modules/Interface/SeVofPhaseFraction.cpp"), "removed SEVOF source scope");
    test(retiredSeVof("include/Modules/Hydrodynamics/HydroSeVofStep.h"), "removed Hydro SEVOF header scope");
    test(retiredSeVof("src/Apps/Electrohydrodynamics/IO/SeVofInput.cpp"), "removed SEVOF reference App scope");
    test(!referenceOnly("src/Modules/Interface/SeVofPhaseFraction.cpp"), "SEVOF has no reference exemption");
    test(!retiredSeVof("src/Modules/Interface/GeometricVofTransport.cpp"), "VOF source is retained");
    test(!referenceOnly("src/Modules/Interface/GeometricVofTransport.cpp"), "VOF remains production");
    test(!allowedInclude("src/Apps/Hydrodynamics/HydrodynamicsApp.cpp", "include/Modules/Interface/SeVofState.h"), "no removed SEVOF include");
    std::cout << "Architecture scanner/policy self-tests: " << cases << " passed\n";
}
static void fixtureTests(const fs::path& scratch) {
    // A fresh, tiny synthetic repository exercises resolution and the actual
    // graph reader without editing, copying or compiling production sources.
    const auto root = fs::absolute(scratch / std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::map<std::string, std::string> files;
    for (const auto& edge : publicationViews) {
        files[edge.first] += "#include \"" + edge.second.substr(8) + "\"\n";
        files.try_emplace(edge.second, "");
    }
    files["include/Modules/Interface/HermiteGalsState.h"] = "";
    files["include/Modules/Interface/Detail/GeometricVofPrivate.h"] = "";
    files["include/Modules/Electric/Detail/ElectricVofPrivate.h"] = "";
    files["include/Modules/Plasma/Particles/ParticleType.h"] = "";
    files["include/Apps/Electric/ElectricApp.h"] = "";
    files["src/Apps/AppRegistry.cpp"] = "";
    files["include/Core/Base/Base.h"] = "";
    files["include/Apps/RunConfig.h"] = "";
    files["src/Core/Base/Probe.cpp"] = "#include \"Core/Base/Base.h\"\n";
    files["src/Modules/Plasma/Probe.cpp"] = "#include <Core/Base/Base.h>\n";
    files["src/Modules/Hydrodynamics/HydroProbe.cpp"] = "#include \"Core/Base/Base.h\"\n";
    const auto write = [&](const std::string& path, const std::string& text) {
        fs::create_directories((root / path).parent_path());
        std::ofstream out(root / path, std::ios::binary);
        out << text;
        require(bool(out), "Cannot write fixture " + path);
    };
    for (const auto& file : files) write(file.first, file.second);
    std::string graph;
    for (const auto* target : {"lamp_compile_options", "lamp_core", "lamp_hermite_geometry", "lamp_linalg", "lamp_solvers",
        "lamp_interface", "lamp_electric", "lamp_hydrodynamics", "lamp_plasma", "LAMP"}) graph += std::string("target|") + target + '\n';
    for (const auto& file : files) {
        if (!starts(file.first, "src/") || !translationUnit(file.first)) continue;
        const auto group = layer(file.first);
        const std::string owner = group == "Core" ? "lamp_core" : group == "Plasma" ? "lamp_plasma" :
            group == "Electric" ? "lamp_electric" : group == "InterfaceCompatibility" ? "lamp_interface_compatibility" : group == "HydrodynamicsHydro" ? "lamp_hydrodynamics" : "LAMP";
        graph += "source|" + owner + '|' + (root / file.first).generic_string() + '\n';
    }
    graph += "link|lamp_electric|lamp_interface\n";
    graph += "link|lamp_hydrodynamics|lamp_electric\n";
    write("graph.txt", graph);
    check(fs::canonical(root), root / "graph.txt");
#ifdef _WIN32
    const auto externalDrive = root.root_name().generic_string() == "C:" ? "D:" : "C:";
    write("graph.txt", graph + "include|lamp_core|" + externalDrive + "/LAMP_external_include\n");
    check(fs::canonical(root), root / "graph.txt");
    write("graph.txt", graph);
    std::cout << "Architecture cross-drive external include fixture: passed\n";
#endif
    const std::string referenceUnit = "src/Apps/Interface/RetainedReference.cpp";
    const std::string referenceHeader = "include/Apps/Interface/RetainedReference.h";
    files[referenceUnit] = ""; files[referenceHeader] = "";
    write(referenceUnit, ""); write(referenceHeader, "");
    check(fs::canonical(root), root / "graph.txt");
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
    reject("src/Modules/Plasma/Probe.cpp", "#include \"RunConfig.h\"\r\n", "Forbidden include");
    reject("src/Modules/Plasma/Probe.cpp", "#include \"../../../include/Apps/RunConfig.h\"\n", "Forbidden include");
    reject("src/Core/Base/Probe.cpp", "#if 0\n#include <Modules/Plasma/Particles/ParticleType.h>\n#endif\n", "Forbidden include");
    reject("src/Core/Base/Probe.cpp", "#include HEADER\n", "literal header");
    reject("graph.txt", graph + "source|lamp_interface|" + (root / referenceUnit).generic_string() + '\n', "Reference-only source in production graph");
    reject("src/Modules/Interface/SeVofReference.cpp", "", "Retired SEVOF source must remain removed");
    reject("include/Modules/Hydrodynamics/HydroSeVofReference.h", "", "Retired SEVOF source must remain removed");
    reject("src/Apps/Interface/SeVofReference.cpp", "", "Retired SEVOF source must remain removed");
    reject("graph.txt", graph + "source|lamp_interface|" + (root / "src/Modules/Interface/SeVofReference.cpp").generic_string() + '\n', "Retired SEVOF source in production graph");
    reject("src/Modules/Hydrodynamics/HydroProbe.cpp", "#include \"Modules/Interface/SeVofReference.h\"\n", "unresolved project header");
    reject("src/Core/Base/Probe.cpp", "#include \"Core/Missing.h\"\n", "unresolved project header");
    reject("src/Apps/New.cpp", "#include \"Modules/Electric/Legacy/Electric.h\"\n", "unresolved project header");
    reject("src/Apps/New.cpp", "#include \"Electric/ElectricApp.h\"\n", "Forbidden include");
    reject("src/Apps/New.cpp", "#include \"Modules/Plasma/Particles/ParticleType.h\"\n", "Forbidden include");
    reject("src/Modules/Electric/Legacy/New.cpp", "", "Unclassified source layout");
    reject("src/Modules/New/New.cpp", "", "Unclassified source layout");
    reject("src/Core/Base/Unbuilt.cpp", "", "exactly one production owner");
    reject("include/Core/Utils/New.h", "", "Unclassified utility directory");
    reject("src/Core/Base/New.h", "", "Headers belong in include/");
    reject("include/Core/Base/New.cpp", "", "Headers belong in include/");
    reject("src/Apps/Electric/Compatibility/LegacyElectricModule.cpp", "", "Unclassified source layout");
    reject("include/Modules/Electric/ElectricPhaseGeometry.h", "", "Stale publication view allowance");
    reject("src/Modules/Electric/ElectricPhaseGeometry.cpp", "", "Stale publication view allowance");
    reject("include/Modules/Electric/Other.h",
        "#include \"Modules/Interface/HermiteGalsPhaseQuadrature.h\"\n", "Forbidden include");
    reject("src/Modules/Electric/ElectricPhaseGeometry.cpp",
        "#include \"Modules/Interface/HermiteGalsState.h\"\n", "Forbidden include");
    reject("include/Modules/Electric/Other.h",
        "#include \"Modules/Interface/HermiteGalsPublication.h\"\n", "Forbidden include");
    reject("include/Modules/Electric/ElectricPhaseGeometry.h",
        "#include \"Modules/Interface/HermiteGalsState.h\"\n", "Forbidden include");
    reject("include/Modules/Interface/HermiteGalsPublication.h",
        "#include \"Modules/Electric/ElectricPhaseGeometry.h\"\n", "Forbidden include");
    // Each exact VOF allowance must be used; removing one E2 include cannot
    // be hidden by the other E2 include or by the older electric publication.
    for (const auto& edge : publicationViews) {
        if (edge.first != "include/Modules/Electric/ElectricVofField.h"
            && edge.first != "include/Modules/Electric/ElectricVofChargeTransport.h") continue;
        const std::string directive = "#include \"" + edge.second.substr(8) + "\"\n";
        auto without = files.at(edge.first);
        const auto position = without.find(directive);
        require(position != std::string::npos, "VOF exact-edge fixture is missing its actual include");
        without.erase(position, directive.size());
        reject(edge.first, without, "Stale publication view allowance");
        reject("include/Modules/Electric/Other.h", directive, "Forbidden include");
        reject(edge.second, "#include \"" + edge.first.substr(8) + "\"\n", "Forbidden include");
        reject(edge.first, files.at(edge.first)
            + "#include \"Modules/Interface/Detail/GeometricVofPrivate.h\"\n", "Forbidden include");
    }
    reject("include/Modules/Electric/ElectricVofField.h",
        files.at("include/Modules/Electric/ElectricVofField.h")
            + "#include \"Modules/Interface/GeometricVofTransport.h\"\n", "Forbidden include");
    reject("include/Modules/Electric/ElectricVofChargeTransport.h",
        files.at("include/Modules/Electric/ElectricVofChargeTransport.h")
            + "#include \"Modules/Interface/GeometricVofPublication.h\"\n", "Forbidden include");
    for (const auto& edge : publicationViews) {
        const bool hydro = edge.first == "include/Modules/Hydrodynamics/HydroMaxwellLoad.h";
        if (edge.first != "include/Modules/Electric/ElectricVofSharpField.h"
            && edge.first != "include/Modules/Electric/ElectricSurfaceChargeTransport.h"
            && !(hydro && edge.second == "include/Modules/Electric/ElectricVofSharpMaxwellStress.h")) continue;
        const std::string directive = "#include \"" + edge.second.substr(8) + "\"\n";
        auto without = files.at(edge.first);const auto position = without.find(directive);
        require(position != std::string::npos, "B1 exact-edge fixture is missing its actual include");
        without.erase(position, directive.size());
        reject(edge.first, without, "Stale publication view allowance");
        reject(hydro ? "include/Modules/Hydrodynamics/HydroOther.h" : "include/Modules/Electric/Other.h", directive, "Forbidden include");
        reject(edge.second, "#include \"" + edge.first.substr(8) + "\"\n", "Forbidden include");
        reject(edge.first, files.at(edge.first) + (hydro
            ? "#include \"Modules/Electric/Detail/ElectricVofPrivate.h\"\n"
            : "#include \"Modules/Interface/Detail/GeometricVofPrivate.h\"\n"), "Forbidden include");
    }
    reject("include/Modules/Electric/ElectricVofSharpField.h", files.at("include/Modules/Electric/ElectricVofSharpField.h")
        + "#include \"Modules/Interface/GeometricVofTransport.h\"\n", "Forbidden include");
    reject("include/Modules/Electric/ElectricSurfaceChargeTransport.h", files.at("include/Modules/Electric/ElectricSurfaceChargeTransport.h")
        + "#include \"Modules/Interface/VofSmoothGeometry.h\"\n", "Forbidden include");
    reject("include/Modules/Electric/ElectricVofField.h", files.at("include/Modules/Electric/ElectricVofField.h")
        + "#include \"Modules/Interface/VofSmoothGeometry.h\"\n", "Forbidden include");
    reject("graph.txt", graph + "link|lamp_interface|lamp_electric\n", "Forbidden target link");
    reject("include/Modules/Hydrodynamics/HydroMaxwellLoad.h","","Stale publication view allowance");
    reject("include/Modules/Hydrodynamics/HydroOther.h",
        "#include \"Modules/Electric/ElectricVofMaxwellStress.h\"\n","Forbidden include");
    reject("include/Modules/Electric/ElectricVofMaxwellStress.h",
        "#include \"Modules/Hydrodynamics/HydroMaxwellLoad.h\"\n","Forbidden include");
    reject("include/Modules/Hydrodynamics/HydroMaxwellLoad.h",files.at("include/Modules/Hydrodynamics/HydroMaxwellLoad.h")+
        "#include \"Modules/Electric/ElectricVofField.h\"\n","Forbidden include");
    reject("include/Modules/Hydrodynamics/HydroMaxwellLoad.h",files.at("include/Modules/Hydrodynamics/HydroMaxwellLoad.h")+
        "#include \"Modules/Electric/Detail/ElectricVofPrivate.h\"\n","Forbidden include");
    reject("graph.txt",graph+"link|lamp_electric|lamp_hydrodynamics\n","Forbidden target link");
    reject("graph.txt",graph+"include|lamp_hydrodynamics|"+(root/"include/Modules/Electric").generic_string()+'\n',"Forbidden include directory");
    auto missingMaxwellLink=graph;const std::string maxwellLink="link|lamp_hydrodynamics|lamp_electric\n";
    missingMaxwellLink.erase(missingMaxwellLink.find(maxwellLink),maxwellLink.size());
    reject("graph.txt",missingMaxwellLink,"Missing exact Maxwell implementation link");
    reject("graph.txt", graph + "include|lamp_electric|" + (root / "include/Modules/Interface").generic_string() + '\n',
        "Forbidden include directory");
    reject("graph.txt", graph + "link|lamp_electric|lamp_legacy_nodal_adapter\n", "Forbidden target link");
    reject("graph.txt", graph + "link|lamp_compile_options|LAMP\n", "Forbidden target link");
    reject("graph.txt", graph + "link|lamp_compile_options|lamp_core\n", "Forbidden target link");
    reject("graph.txt", graph + "include|lamp_compile_options|" + (root / "include").generic_string() + '\n', "Build settings contain");
    reject("graph.txt", graph + "source|lamp_compile_options|" + (root / "src/Core/Base/Probe.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "include|lamp_plasma|" + (root / "include/Apps").generic_string() + '\n', "Forbidden include directory");
    reject("graph.txt", graph + "source|lamp_plasma|" + (root / "src/Apps/AppRegistry.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "source|lamp_core|" + (root / "src/Core/Base/Probe.cpp").generic_string() + '\n', "exactly one production owner");
    reject("graph.txt", graph + "link|lamp_hydrodynamics|lamp_legacy_nodal_adapter\n", "Forbidden target link");
    reject("graph.txt", graph + "link|lamp_hydrodynamics|lamp_interface_compatibility\n", "Forbidden target link");
    reject("graph.txt", graph + "link|lamp_hydrodynamics|LAMP\n", "Forbidden target link");
    reject("graph.txt", graph + "source|lamp_legacy_nodal_adapter|" + (root / "src/Modules/Hydrodynamics/HydroProbe.cpp").generic_string() + '\n', "Unclassified production target");
    reject("graph.txt", graph + "source|lamp_hydrodynamics|" + (root / "src/Core/Base/Probe.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "source|lamp_hydrodynamics|" + (root / "src/Modules/Hydrodynamics/LegacyNodalState.cpp").generic_string() + '\n', "Wrong source owner");
    reject("graph.txt", graph + "source|lamp_hydrodynamics|" + (root / "src/Modules/Hydrodynamics/HydroProbe.cpp").generic_string() + '\n', "exactly one production owner");
    reject("src/Modules/Hydrodynamics/HydroProbe.cpp", "#include \"Modules/Hydrodynamics/LegacyNodalState.h\"\n", "unresolved project header");
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
