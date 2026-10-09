#pragma once

#define thread_num 12
#define MINVAL 1e-8
#define MAXVAL 1e8
#define Pi 3.14159265358979323846

#include <iostream>
#include <cmath>

const double epsilon0 = 8.854187817e-12;// 真空介电常数(F/m) 
const double epsilonLiq = 8.854187817e-11;// 液体介电常数(F/m) 
const double temperature = 298.0; // K

const double character_rho = 1240.0; // kg/m^3
const double character_mu = 0.038; // Pa.s
const double character_sigma = 0.0452; // N/m
const double character_r = 1e-6; // m
// const double character_u = 0.25; // m/s
const double character_u = sqrt(character_sigma / (character_rho * character_r)); // m/s
const double character_E = 5e7;// 特征电场强度(V/m)
const double character_m = character_rho * pow(character_r, 3.0);// 质量量纲
const double character_e = epsilon0 * character_E * character_r * character_r;// 电荷量纲
const double character_t = character_r / character_u; // s

const double conductivity = 1.4;// 电导率(S/m)

const double number_Boltzmann = 1.38064852e-23; // J/K
const double number_Planck = 6.62607004e-34; // J.s
const double number_Oh = character_mu / sqrt(character_rho * character_sigma * character_r); // 无量纲Ohnesorge数
const double number_Re = character_rho * character_u * character_r / character_mu; // 无量纲Re数
const double number_Re_inv = 1.0 / number_Re; // 无量纲Re数的倒数
const double number_We = character_rho * character_u * character_u * character_r / character_sigma; // 无量纲We数
const double number_We_inv = 1.0 / number_We; // 无量纲We数的倒数
const double number_Bond = epsilon0 * character_E * character_E * character_r / character_sigma; // 无量纲Bond数

const double particleCharge = 1.602e-19;// 元电荷(C)
const double activeEnergy_anion = 1.602e-19;// 阴离子活化能(J),1eV
const double activeEnergy_cation = 1.00 * activeEnergy_anion;// 阳离子活化能(J),1eV
const double amu = 1.661e-27;// 原子质量单位(kg)
const double cationMass = (111.16 + (111.16 + 86.81)*0.0) * amu;// 阳离子质量(kg)
const double anionMass = (86.81 + (111.16 + 86.81)*0.0) * amu;// 阴离子质量(kg)

const double nonDim_epsilon = epsilonLiq / epsilon0;// 介电常数比
const double number_lambda = epsilonLiq / (conductivity * character_t);// 电荷弛豫时间与毛细时间比

const double CHI_cation = activeEnergy_cation / (number_Boltzmann * temperature);
const double CHI_anion = activeEnergy_anion / (number_Boltzmann * temperature);

const double chi = number_Boltzmann * temperature / number_Planck * character_t; // 毛细时间与场蒸发时间之比
// const double chi = 551984.0;

const double xi_cation = sqrt(pow(particleCharge, 3.0) * character_E / (4.0 * Pi * epsilon0)) / activeEnergy_cation;// 阳离子电能与活化能之比
const double xi_anion = sqrt(pow(particleCharge, 3.0) * character_E / (4.0 * Pi * epsilon0)) / activeEnergy_anion;// 阴离子电能与活化能之比

const double nonDim_rho_l = 1.0;
const double nonDim_rho_g = 1.184 / character_rho;
const double nonDim_mu_l = 1.0;
const double nonDim_mu_g = 1.849e-5 / character_mu;

const double nonDim_anodePotential = 8.0;// 阳极电位
const double nonDim_cathodePotential = -8.0;// 阴极电位

const double g = 0.0 / character_r; // m/s^2
const double nonDim_g = g * character_r / (character_u * character_u); // 无量纲重力加速度
