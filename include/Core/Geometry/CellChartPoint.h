#pragma once

#include "Core/Base/Vector2D.h"

namespace Quantum::Geometry
{
// A dimensionless cell coordinate retained relative to either endpoint.
// Its mathematical value is offset (lower) or 1-offset (upper). Never form
// the latter in floating point to recover topology: a positive tail may be
// smaller than one ULP at one. Both orientations admit finite offsets [0,1].
struct CellChartCoordinate
{
    double offset = 0.0;
    bool fromUpper = false;
};

struct CellChartPoint
{
    CellChartCoordinate x, y;
};

// Out-of-line strict-FP helpers; no expansion arithmetic crosses the ABI.
CellChartPoint cellChartFromLocal(Vector2D local);
void validateCellChart(const CellChartPoint& point);
bool cellChartInterior(const CellChartPoint& point);
bool cellChartEqual(const CellChartPoint& a, const CellChartPoint& b);
bool cellChartEqual(const CellChartPoint& a, Vector2D local);
int compareCellChartCoordinate(CellChartCoordinate a, CellChartCoordinate b);
// Conversion is permitted only if the mathematical coordinate is EXACTLY
// representable as one double. In particular, 1-tiny is not rounded to one.
double cellChartCoordinateValue(CellChartCoordinate coordinate);
// Terminal integration parameter on an original represented interval.
// Bounds/endpoint identity are exact; a true interior fraction must remain
// strictly interior after conversion. Reversed intervals are supported.
double cellChartIntervalFraction(CellChartCoordinate coordinate, double begin, double end);
}
