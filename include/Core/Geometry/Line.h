#pragma once

#include <iostream>
#include <vector>
#include "Core/Base/Enumeration.h"
#include "Core/Base/Vector2D.h"

struct Parameters
{
    Vector2D arcCenter;
    double radius;
    double radian;

    std::vector<Vector2D> featurePoints;
};

class Line
{
public:
    Vector2D points[2];

    int levelSetIndex = -1;
    int phaseIndex = -1;
    int opposedPhaseIndex = -1;

    Line();
    virtual ~Line() = default;

    void reverse();

    void setFeaturePoint(int i, Vector2D _featurePoint);

    inline LineType type() const { return type_; }
    inline Vector2D arcCenter() const { return parameters_.arcCenter; }
    inline double radius() const { return parameters_.radius; }
    inline double radian() const { return parameters_.radian; }
    inline Vector2D featurePoint(int i) const { return parameters_.featurePoints[i]; }
    inline std::vector<Vector2D> featurePoints() const { return parameters_.featurePoints; }
    inline GeometricBoundary geoBoundary() const { return geoBoundary_; }
    inline bool ignored() const { return ignored_; }
    inline bool immersed() const { return immersed_; }

    inline void setType(LineType _type) { type_ = _type; }
    inline void setArcCenter(Vector2D _arcCenter) { parameters_.arcCenter = _arcCenter; }
    inline void setRadius(double _radius) { parameters_.radius = _radius; }
    inline void setRadian(double _radian) { parameters_.radian = _radian; }
    inline void addFeaturePoint(Vector2D _featurePoint) { parameters_.featurePoints.push_back(_featurePoint); }
    inline void setGeoBoundary(GeometricBoundary _geoBoundary) { geoBoundary_ = _geoBoundary; }
    inline void setIgnored(bool _ignored) { ignored_ = _ignored; }
    inline void setImmersed(bool _immersed) { immersed_ = _immersed; }

private:
    LineType type_ = LineType::straight;
    Parameters parameters_;

    GeometricBoundary geoBoundary_;

    bool ignored_ = false;
    bool immersed_ = false;
};