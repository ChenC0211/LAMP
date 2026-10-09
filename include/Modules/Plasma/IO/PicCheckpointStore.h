#pragma once
#include <filesystem>
#include <random>
#include <string>
#include <vector>

struct PicCaseInput;
struct CollisionsConfig;
enum class PicIntegrator;

// Legacy PIC checkpoint directory/metadata policy, shared by both drivers.
// Component caches remain owned by their physical services.
class PicCheckpointStore
{
public:
    struct Metadata
    {
        int version = 1;
        int nextStep = 0;
        double t_code = 0.0;
        double dt_code = 0.0;
    };
    struct Requirements
    {
        std::vector<std::string> species;
        bool geometry = false;
        bool implicitFields = false;
        bool movingBoundary = false;
        std::filesystem::path movingBoundaryCache;
        bool collisionState = false;
    };
    inline static constexpr const char* metaName = "pic_checkpoint.txt";
    inline static constexpr const char* rngName = "rng_state.txt";
    inline static constexpr const char* fieldName = "implicit_pic_fields_cache.bin";
    inline static constexpr const char* collisionName = "pic_collisions_state.txt";
    inline static constexpr const char* incompleteName = ".pic_checkpoint_incomplete";

    PicCheckpointStore(std::filesystem::path base, Requirements requirements);
    // Snapshot the prepared run's component requirements without accessing disk.
    // Explicit geometry also requires named level sets; implicit always requires
    // its field history. Keep these legacy selection/completion rules distinct.
    static PicCheckpointStore forRun(PicIntegrator integrator, const PicCaseInput& input,
        const CollisionsConfig* collisions);
    std::filesystem::path directory(int step) const;
    std::filesystem::path latest() const;
    // Legacy explicit paths remain supported; interrupted new writes are rejected.
    std::filesystem::path select(int step) const;
    bool usable(const std::filesystem::path& directory) const;
    bool createDirectory(int step, std::filesystem::path& directory) const;
    void beginWrite(const std::filesystem::path& directory) const;
    void completeWrite(const std::filesystem::path& directory, bool requireRng) const;
    void keepLatest(int count) const;
    static std::string pathPrefix(const std::filesystem::path& directory);
    static bool readMetadata(const std::filesystem::path& directory, Metadata& metadata);
    // The explicit output preflight historically skips the version line.
    static bool readOutputMetadata(const std::filesystem::path& directory, Metadata& metadata);
    static bool writeMetadata(const std::filesystem::path& directory, int nextStep, double time, double dt);
    static bool readRng(const std::filesystem::path& directory, std::mt19937_64& rng);
    static bool writeRng(const std::filesystem::path& directory, const std::mt19937_64& rng);

private:
    const std::filesystem::path base_;
    const Requirements requirements_;
    bool componentsPresent(const std::filesystem::path& directory) const;
    static bool parseStep(const std::string& name, int& step);
    static bool readMetadataImpl(const std::filesystem::path& directory, Metadata& metadata, bool readVersion);
};
