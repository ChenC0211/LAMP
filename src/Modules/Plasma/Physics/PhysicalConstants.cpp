#include "Modules/Plasma/Physics/PhysicalConstants.h"

#include <iomanip>
#include <ostream>
#include <sstream>

void PhysicalConstants::PrintTable(std::ostream& os)
{
    os << "[LAMP] PhysicalConstants (SI)\n";
    os << std::scientific << std::setprecision(16);

    os << "  epsilon0 [F/m]   = " << epsilon0_Fpm << "\n";
    os << "  mu0      [H/m]   = " << mu0_Hpm << "\n";
    os << "  c0       [m/s]   = " << c0_mps << "\n";
    os << "  kB       [J/K]   = " << kB_JpK << "\n";
    os << "  qe       [C]     = " << qe_C << "\n";
    os << "  me       [kg]    = " << me_kg << "\n";
    os << "  mp       [kg]    = " << mp_kg << "\n";
    os << "  amu      [kg]    = " << amu_kg << "\n";
    os << "  eV       [J]     = " << eV_J << "\n";
    os << "  NA       [1/mol] = " << NA_1pmol << "\n";
    os << "  R        [J/(mol*K)] = " << R_JpMolK << "\n";
}

std::string PhysicalConstants::TableString()
{
    std::ostringstream oss;
    PrintTable(oss);
    return oss.str();
}
