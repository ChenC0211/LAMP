#pragma once

enum class LineType { straight, arc, spline };

enum class GeometricBoundary { domain, fixed, moving };

enum class BoundaryCondition { dirichlet, neumann };

enum class LevelSetOperation { non, min, max };

enum class Boundary { non, west, south, east, north, southWest, southEast, northEast, northWest };

enum class TVertType { non, west, south, east, north };

enum class DifferenceOrder { non, dx, dy, dxx, dyy, dxy };

enum class Time { current, previous, prePrevious };