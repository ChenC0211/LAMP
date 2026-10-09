#include "Core/Mesh/Vert.h"

using namespace std;

Vert::Vert() {}

Vert::Vert
(
    int label,
    Vector2D coord
) :
    label_(label),
    coord_(coord)
{}

Vert::~Vert() {}