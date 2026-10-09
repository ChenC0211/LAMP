#include "Modules/Plasma/Physics/BuiltInCrossSections.h"

#include <cmath>
#include <memory>
#include <string>

#include "Modules/Plasma/Physics/FunctionalCrossSection.h"

namespace
{
constexpr double kPi = 3.141592653589793238462643383279502884;

inline double log10_safe(double x)
{
    return std::log10((x > 0.0) ? x : 1e-300);
}

// -----------------------------------------------------------------------------
// Built-in Ar cross-section set v1
//   Convention: energy in eV, sigma in m^2
// -----------------------------------------------------------------------------

double eAr_elastic(double E)
{
    if (E < 0.0) return 0.0;

    if (E < 0.21)
        return 2.6082e-13 * std::pow(E, 8) - 3.8935e-13 * std::pow(E, 7) + 2.3977e-13 * std::pow(E, 6) -
               7.873e-14 * std::pow(E, 5) + 1.4878e-14 * std::pow(E, 4) - 1.6321e-15 * std::pow(E, 3) +
               1.008e-16 * E * E - 3.4053e-18 * E + 6.6976e-20;
    if (E <= 0.9953)
        return -1.0254e-19 * std::pow(E, 5) + 3.665e-19 * std::pow(E, 4) - 5.0768e-19 * std::pow(E, 3) +
               3.3939e-19 * E * E - 9.09e-20 * E + 9.1308e-21;
    if (E <= 25.6073)
        return -1.0991e-28 * std::pow(E, 8) + 1.5096e-26 * std::pow(E, 7) - 8.3753e-25 * std::pow(E, 6) +
               2.3939e-23 * std::pow(E, 5) - 3.6958e-22 * std::pow(E, 4) + 2.9261e-21 * std::pow(E, 3) -
               1.0744e-20 * E * E + 3.1043e-20 * E - 9.688e-21;
    if (E <= 97.8553)
        return 6.4428e-33 * std::pow(E, 8) - 3.0209e-30 * std::pow(E, 7) + 6.0335e-28 * std::pow(E, 6) -
               6.7123e-26 * std::pow(E, 5) + 4.5663e-24 * std::pow(E, 4) - 1.961e-22 * std::pow(E, 3) +
               5.2804e-21 * E * E - 8.4786e-20 * E + 7.0372e-19;
    if (E <= 575.0)
        return 1.8699e-42 * std::pow(E, 8) - 1.242e-38 * std::pow(E, 7) + 3.2053e-35 * std::pow(E, 6) -
               4.321e-32 * std::pow(E, 5) + 3.3718e-29 * std::pow(E, 4) - 1.5735e-26 * std::pow(E, 3) +
               4.3523e-24 * E * E - 6.7724e-22 * E + 5.1844e-20;

    return std::pow(10.0, (-0.9296 * log10_safe(E) - 18.03));
}

double eAr_ionization(double E)
{
    if (E < 15.8)
        return 0.0;
    if (E < 49.78)
        return 2.1173e-28 * std::pow(E, 6) - 4.1353e-26 * std::pow(E, 5) + 3.2899e-24 * std::pow(E, 4) -
               1.3614e-22 * std::pow(E, 3) + 3.0549e-21 * E * E - 3.3567e-20 * E + 1.3716e-19;
    if (E < 154.3)
        return 2.8278e-32 * std::pow(E, 6) - 1.955e-29 * std::pow(E, 5) + 5.405e-27 * std::pow(E, 4) -
               7.5343e-25 * std::pow(E, 3) + 5.3962e-23 * E * E - 1.7472e-21 * E + 4.3804e-20;
    if (E < 950.0)
        return -7.4976e-38 * std::pow(E, 6) + 2.1918e-34 * std::pow(E, 5) - 2.0094e-31 * std::pow(E, 4) +
               4.1294e-31 * std::pow(E, 3) + 1.2192e-25 * E * E - 9.4819e-23 * E + 3.838e-20;

    return std::pow(10.0, (-0.6099 * log10_safe(E) - 18.2));
}

double eAr_excitation1(double E)
{
    if (E < 11.59)
        return 0.0;
    if (E < 16.0)
        return -2.089e-22 * (std::sin(E - kPi)) + 5.503e-23 * std::pow(E - 10.0, 2) + 9.626e-23;
    if (E < 48.0)
        return -1.539e-23 * (std::sin(E - kPi)) + 8.882e-25 * std::pow(E - 10.0, 2) + 2.11e-21;
    if (E < 950.0)
        return 2.112e-20 * std::pow(E, -0.4831);

    return std::pow(10.0, (-0.7155 * log10_safe(E) - 19.16));
}

double eAr_excitation2(double E)
{
    if (E < 13.04)
        return 0.0;
    if (E < 20.0)
        return -1.857e-23 * E * E + 1.202e-21 * E - 1.247e-20;
    if (E < 30.2)
        return -4.56e-25 * E * E + 1.425e-22 * E + 1.526e-21;
    if (E < 40.0)
        return -1.559e-26 * std::pow(E, 4.0) + 2.152e-24 * std::pow(E, 3.0) - 1.111e-22 * std::pow(E, 2.0) +
               2.53e-21 * E - 1.602e-20;
    if (E < 950.0)
        return -1.287e-43 * std::pow(E, 8.0) + 2.918e-40 * std::pow(E, 7.0) + 1.471e-37 * std::pow(E, 6.0) -
               1.103e-33 * std::pow(E, 5.0) + 1.464e-30 * std::pow(E, 4.0) - 9.49e-28 * std::pow(E, 3.0) +
               3.377e-25 * std::pow(E, 2.0) - 6.663e-23 * E + 7.546e-21;

    return std::pow(10.0, (-0.8218 * log10_safe(E) - 18.74));
}

double eAr_excitation3(double E)
{
    if (E < 14.01)
        return 0.0;
    if (E <= 77.0)
        return -6.415e-28 * std::pow(E, 4) + 1.56e-25 * std::pow(E, 3) - 1.406e-23 * E * E + 5.861e-22 * E -
               5.782e-21;
    if (E <= 600.0)
        return -2.462e-42 * std::pow(E, 8) + 1.117e-38 * std::pow(E, 7) - 2.129e-35 * std::pow(E, 6) +
               2.213e-32 * std::pow(E, 5) - 1.361e-29 * std::pow(E, 4) + 5.016e-27 * std::pow(E, 3) -
               1.055e-24 * E * E + 1.034e-22 * E + 1.05e-21;
    if (E < 950.0)
        return 9.948e-28 * E * E - 3.026e-24 * E + 3.146e-21;

    return std::pow(10.0, (-0.7558 * log10_safe(E) - 18.67));
}

double arIon_ar_elastic(double E)
{
    if (E <= 0.0) return 0.0;
    return 1e-20 * (20.0 / std::pow(2.0 * E, 0.5) / (1.0 + 2.0 * E) + 30.0 * 2.0 * E / std::pow(1.0 + 2.0 * E / 3.0, 2.3));
}

double arIon_ar_charge_exchange(double E)
{
    if (E <= 0.0) return 0.0;
    return 0.5e-20 * (115.0 * std::pow(1.0 + 0.015 / 2.0 / E, 0.6) / std::pow(2.0 * E, 0.1) -
                      20.0 / std::pow(2.0 * E, 0.5) / (1.0 + 2.0 * E) -
                      30.0 * 2.0 * E / std::pow(1.0 + 2.0 * E / 3.0, 2.3));
}

double arIon_ar_ionization(double E)
{
    if (E < 39.0)
        return 0.0;
    if (E < 1000.0)
        return 1e-20 * (2.309 * std::sin(0.001579 * E - 0.04039) + 0.05629 * std::sin(0.007407 * E - 1.268));

    return 1e-20 * (-574.8 * std::pow(E, -0.9485) + 3.109);
}

double arIon_ar_excitation1(double E)
{
    if (E < 24.1)
        return 0.0;
    if (E < 56.23)
        return 1e-20 * (0.236 * std::pow(E, 0.3704) - 0.7651);
    if (E < 749.89)
        return 1e-20 * (-9.642e-17 * std::pow(E, 6) + 3.045e-13 * std::pow(E, 5) - 3.787e-10 * std::pow(E, 4) +
                        2.376e-07 * std::pow(E, 3) - 7.977e-05 * std::pow(E, 2) + 0.01392 * E - 0.2754);
    if (E < 2371.4)
        return 1e-20 * (-8.913e-05 * E + 0.8399);

    return 1e-20 * (-0.1791 * std::pow(E, 0.1765) + 1.343);
}

double arIon_ar_excitation2(double E)
{
    if (E < 74.0)
        return 0.0;
    if (E < 1000.0)
        return 1e-20 * (-6.342e-14 * std::pow(E, 4) + 2.063e-10 * std::pow(E, 3) - 2.545e-07 * std::pow(E, 2) +
                        0.0001445 * E - 0.00882);
    if (E < 2400.0)
        return 1e-20 * (3.116e-16 * std::pow(E, 4) - 1.505e-13 * std::pow(E, 3) - 4.963e-09 * std::pow(E, 2) +
                        1.018e-05 * E + 0.01862);

    return 1e-20 * (1.024e-18 * std::pow(E, 4) - 4.52e-14 * std::pow(E, 3) + 6.45e-10 * std::pow(E, 2) -
                    3.959e-06 * E + 0.02905);
}

double arIon_ar_excitation3(double E)
{
    if (E < 125.0)
        return 0.0;
    if (E < 2327.0)
        return 1e-20 * (3.19e-16 * std::pow(E, 4) - 1.674e-12 * std::pow(E, 3) + 6.114e-10 * std::pow(E, 2) +
                        9.334e-06 * E - 0.001169);

    return 1e-20 * (-1.584e-18 * std::pow(E, 4) + 5.281e-14 * std::pow(E, 3) - 6.72e-10 * std::pow(E, 2) +
                    3.999e-06 * E + 0.005948);
}

double arFast_ar_elastic(double E)
{
    // In this built-in set, "fast-neutral" Ar collisions reuse the same elastic fit as Ar+-Ar.
    return arIon_ar_elastic(E);
}

double arFast_ar_ionization(double E)
{
    if (E < 55.2)
        return 0.0;
    return 1e-20 * (-11.6 * std::pow(E, -0.1171) + 7.255);
}

double arFast_ar_excitation1(double E)
{
    if (E < 55.2)
        return 0.0;
    return 1e-20 * (-11.6 * std::pow(E, -0.1171) + 7.255);
}

double arFast_ar_excitation2(double E)
{
    if (E < 88.0)
        return 0.0;
    return 1e-20 * (-0.794 * std::pow(E, -0.4493) + 0.1068);
}

double arFast_ar_excitation3(double E)
{
    if (E < 47.0)
        return 0.0;
    return 1e-20 * (-5.022 * std::pow(E, -0.6449) + 0.4203);
}

std::shared_ptr<ICrossSectionModel> make_xs(double (*f)(double), const std::string& name)
{
    return std::make_shared<FunctionalCrossSection>([f](double E) { return f(E); }, name);
}

bool registerOne(CrossSectionDB* db, const std::string& key, std::shared_ptr<ICrossSectionModel> model, bool overwrite, std::string* err)
{
    if (!db)
    {
        if (err) *err = "CrossSectionDB pointer is null";
        return false;
    }
    return db->Register(key, std::move(model), overwrite, err);
}
}

namespace BuiltInCrossSections
{

bool RegisterArV1(CrossSectionDB* db, const std::string& prefix, bool overwrite, std::string* err)
{
    std::string local_err;

    const std::string p = prefix;

    struct Item { const char* suffix; double (*f)(double); const char* name; };

    const Item items[] = {
        {"e-Ar:ELASTIC", eAr_elastic, "Ar v1 e-Ar elastic"},
        {"e-Ar:IONIZATION", eAr_ionization, "Ar v1 e-Ar ionization"},
        {"e-Ar:EXCITATION_1", eAr_excitation1, "Ar v1 e-Ar excitation 1"},
        {"e-Ar:EXCITATION_2", eAr_excitation2, "Ar v1 e-Ar excitation 2"},
        {"e-Ar:EXCITATION_3", eAr_excitation3, "Ar v1 e-Ar excitation 3"},

        {"Ar+-Ar:ELASTIC", arIon_ar_elastic, "Ar v1 Ar+-Ar elastic"},
        {"Ar+-Ar:CHARGE_EXCHANGE", arIon_ar_charge_exchange, "Ar v1 Ar+-Ar charge exchange"},
        {"Ar+-Ar:IONIZATION", arIon_ar_ionization, "Ar v1 Ar+-Ar ionization"},
        {"Ar+-Ar:EXCITATION_1", arIon_ar_excitation1, "Ar v1 Ar+-Ar excitation 1"},
        {"Ar+-Ar:EXCITATION_2", arIon_ar_excitation2, "Ar v1 Ar+-Ar excitation 2"},
        {"Ar+-Ar:EXCITATION_3", arIon_ar_excitation3, "Ar v1 Ar+-Ar excitation 3"},

        {"ArFast-Ar:ELASTIC", arFast_ar_elastic, "Ar v1 ArFast-Ar elastic"},
        {"ArFast-Ar:IONIZATION", arFast_ar_ionization, "Ar v1 ArFast-Ar ionization"},
        {"ArFast-Ar:EXCITATION_1", arFast_ar_excitation1, "Ar v1 ArFast-Ar excitation 1"},
        {"ArFast-Ar:EXCITATION_2", arFast_ar_excitation2, "Ar v1 ArFast-Ar excitation 2"},
        {"ArFast-Ar:EXCITATION_3", arFast_ar_excitation3, "Ar v1 ArFast-Ar excitation 3"},
    };

    for (const auto& it : items)
    {
        const std::string key = p + it.suffix;
        if (!registerOne(db, key, make_xs(it.f, it.name), overwrite, &local_err))
        {
            if (err) *err = local_err;
            return false;
        }
    }

    return true;
}

}  // namespace BuiltInCrossSections
