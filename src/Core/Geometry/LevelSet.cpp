#include "Core/Geometry/LevelSet.h"

using namespace std;

LevelSet::LevelSet() : VertField<double>() {}

LevelSet::~LevelSet() {}

unique_ptr<double[]> LevelSet::convertToHyperTan(double thinckness) const
{
    unique_ptr<double[]> hyperTanValue(new double[size_]);
    for (int i = 0; i < size_; ++i)
    {
        hyperTanValue[i] = 0.5 * (1.0 + tanh(value[i] / (2.0 * thinckness)));
    }
    return hyperTanValue;
}
