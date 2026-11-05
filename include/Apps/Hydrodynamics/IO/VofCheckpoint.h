#pragma once
#include "Hydrodynamics/IO/HydrodynamicsInput.h"
#include "Core/Base/TimeLayer.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

struct VofCheckpointState {
    LAMP::TimeLayer clock;
    double nextDt=0;
    int withoutRetry=0;
    std::string transaction;
    std::filesystem::path directory;
};
// The App owns persistence, compatibility and publication. Modules provide
// sealed state payloads; no field is inferred from plotting output.
class VofCheckpointStore {
public:
    VofCheckpointStore(const HydrodynamicsInput&,std::string app);
    std::optional<VofCheckpointState> read() const;
    void restoreMesh(QuadtreeMesh&,const VofCheckpointState&) const;
    void write(const QuadtreeMesh&,const VofCheckpointState&);
private:
    const HydrodynamicsInput& input_;
    std::string app_,compatibility_;
    std::filesystem::path root_;
    std::vector<std::filesystem::path> generations_;
};
