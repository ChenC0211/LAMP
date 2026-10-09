#pragma once
#include "Modules/Plasma/Boundary/LegacyGridMotion.h"

// Owns the undeformed reference and the last accepted prescribed motion.
// The caller rebuilds mesh-dependent geometry before accepting a candidate.
class LegacyGridMotionEvolution
{
public:
    struct Reference
    {
        std::vector<Line> lines;
        std::vector<Vector2D> innerPoints;
    };
    struct InnerPointBindings
    {
        std::vector<int> accel;
        std::vector<int> decel;
    };

    LegacyGridMotionEvolution(LegacyGridMotionConfig config, LegacyGridMotion::Bindings lines,
                              Reference reference, InnerPointBindings innerPoints);
    LegacyGridMotionState prepare(int step, bool force = false) const;
    void deformReference(std::vector<Line>& lines, const LegacyGridMotionState& state) const;
    std::size_t innerPointCount() const { return reference_.innerPoints.size(); }
    Vector2D innerPoint(std::size_t index, const LegacyGridMotionState& state) const;
    void accept(LegacyGridMotionState state, const std::vector<Line>& lines);
    const LegacyGridMotionState& current() const { return applied_; }

    // Checkpoint nextStep counts completed steps; it is not the last motion index.
    // Same input/reference required. Loaded geometry is checked, never moved here.
    void restore(int nextStep, const std::vector<Line>& checkpointLines);

private:
    const LegacyGridMotionConfig config_;
    const LegacyGridMotion motion_;
    const Reference reference_;
    const InnerPointBindings innerPoints_;
    LegacyGridMotionState applied_;
    bool haveApplied_ = false;
};
