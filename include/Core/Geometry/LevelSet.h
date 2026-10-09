#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <type_traits>
#include "Core/Base/Vector2D.h"
#include "Core/Field/VertField.h"
#include "Core/Base/Structures.h"

class LevelSet : public VertField<double>
{
public:
    LevelSet();
    LevelSet(const LevelSet&) = default;
    LevelSet& operator=(const LevelSet&) = default;

    template<typename Mesh, typename = std::enable_if_t<!std::is_same_v<std::decay_t<Mesh>, LevelSet>>>
    LevelSet(Mesh& mesh);
    ~LevelSet();

    std::unique_ptr<double[]> convertToHyperTan(double thinckness) const;

    inline Vector2D innerPoint() const { return innerPoint_; }
    inline std::vector<SignInt> lines() const { return lines_; }
    inline bool moving() const { return moving_; }
    inline double reinitFineEps() const { return reinitFineEps_; }
    inline double reinitCoarseEps() const { return reinitCoarseEps_; }
    inline double narrowbandScale() const { return narrowbandScale_; }

    inline void setInnerPoint(const Vector2D& innerPoint) { innerPoint_ = innerPoint; }
    inline void setLines(const std::vector<SignInt>& lines) { lines_ = lines; }
    inline void setMoving(const bool moving) { moving_ = moving; }
    inline void setReinitFineEps(const double reinitFineEps) { reinitFineEps_ = reinitFineEps; }
    inline void setReinitCoarseEps(const double reinitCoarseEps) { reinitCoarseEps_ = reinitCoarseEps; }
    inline void setNarrowbandScale(const double narrowbandScale) { narrowbandScale_ = narrowbandScale; }

private:
    Vector2D innerPoint_;
    std::vector<SignInt> lines_;

    bool moving_ = false;

    double reinitFineEps_ = 1e-4;
    double reinitCoarseEps_ = 1e-3;

    double narrowbandScale_ = 5.0;
};

#include "Core/Geometry/LevelSet_impl.h"
