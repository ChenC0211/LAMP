#pragma once

#include <cstddef>
#include <vector>
#include "Core/Base/Vector2D.h"

// Caller-owned optional observations of the ACTUAL sharp global assembly.
// Pure data only: no solver/Detail engine, shape reference or reconstruction.
// One record per corrected physical interior face, with the graph's original
// owner-outward orientation (not necessarily the positive Cartesian normal).
struct DielectricSharpFaceDiagnostic
{
    struct Coefficient { int cellIndex = -1; double value = 0.0; };
    struct Piece { Vector2D begin, end; bool negative = false; };

    std::size_t graphSegmentIndex = 0;
    int edgeIndex = -1, ownerCellIndex = -1, neighborCellIndex = -1;
    bool horizontal = false, ownerNegative = false, neighborNegative = false;
    Vector2D faceCenter, outwardNormal, interfacePoint, interfaceNormal;
    double faceArea = 0.0, centerDistance = 0.0, epsilonFace = 0.0;
    double linkNegativeFraction = 0.0, faceNegativeFraction = 0.0;
    double curvature = 0.0, curvatureTangentialDerivative = 0.0;
    double condition = 0.0, reproduction = 0.0, correctionConstant = 0.0;
    std::vector<int> sampleCellIndices;
    std::vector<unsigned char> sampleNegative; // 0/1; same actual sample order.
    std::vector<Piece> facePieces; // Actual ordered phase partition used in Q.
    std::vector<Coefficient> correction; // Final sorted nonlocal Q-TP weights.
    double physicalFaceMeasure = 0.0; // faceArea in planar, 2*pi*r*faceArea in RZ.
};

struct DielectricSharpAssemblyDiagnostics
{
    // Replaced only after a successful sharp solve. Empty after a successful
    // homogeneous-material solve. On failure all values/storage are retained.
    std::vector<DielectricSharpFaceDiagnostic> correctedFaces;
};
