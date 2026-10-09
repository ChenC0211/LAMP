#pragma once

#include <vector>

struct SignInt
{
    int value = 0;
    bool isNegative = false;
};

struct IIMCoeff
{
    double k_vac = 0.0;
    double k_liq = 0.0;
    double f_vac = 0.0;
    double f_liq = 0.0;
    double beta_vac = 0.0;
    double beta_liq = 0.0;
    double omega = 0.0;
    double lambda = 0.0;
};

struct DiscCoeff
{
    int discNum = 0;
    double discOpr[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double rhs = 0.0;
};

struct supCoeff
{
    std::vector<int> ind;
    std::vector<double> coe;
};

struct virtualSupCell
{
    supCoeff center;
    supCoeff vEdges[2];
    supCoeff hEdges[2];
};
