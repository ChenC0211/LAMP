#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include "Modules/Plasma/Particles/Particle.h"

class Particle;

class ParticleType
{
public:
    ParticleType(const string& name, const double& mass, const double& charge);
    ~ParticleType();

    int number = 0;

    std::vector<Particle> particle;

    void output(const string& name);

    void writeCache(const std::string& filePath) const;
    void readCache(const std::string& filePath);

    inline const string& name() const { return name_; }
    inline const double& mass() const { return mass_; }
    inline const double& charge() const { return charge_; }
    inline const double& q_m() const { return q_m_; }

    // Normalized acceleration coefficient for this species (dimensionless):
    //   v_code^{n+1/2} = v_code^{n-1/2} + kappa_s * E_code(x^n) * dt_code.

    // kappa_s is computed in the module from physical q/m and normalization scales:
    //   With v_code = v/v0, E_code = E/E0, v0 = L0/T0, E0 = Phi0/L0,
    //   dv_code = (q/m) * (E0/v0) * E_code * dt_SI.
    //   Since the pusher uses dt_code = dt_SI / T0,
    //   dv_code = (q/m) * (E0/v0) * E_code * (dt_code * T0)
    //           = (q/m) * (Phi0*T0^2/L0^2) * E_code * dt_code.
    inline double kappa_s() const { return kappa_s_; }
    inline void set_kappa_s(const double kappa) { kappa_s_ = kappa; }

private:
    string name_;
        
    double mass_;
    double charge_;

    // 荷质比
    double q_m_;

    // Normalized acceleration coefficient (dimensionless).
    double kappa_s_ = 0.0;
};