#include "Modules/Electricity/Legacy/Electricity.h"

using namespace std;

Electricity::~Electricity() {}

void Electricity::computeInterfaceRhoe(LegacyNodalState& sharpInterface, const int levelSetI, const double& dt, const double& preDt, const int _loop)
{
    if (_loop == 1)
    {
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(); i++)
        {    
            double a = 0.0;
			if (EnVac[levelSetI].value[i] > 0.0)
			{
                a = 1.0 / dt - expansivity[levelSetI].value[i] + chi * exp(activeEnergy_cation / (number_Boltzmann * temperature) * (xi_cation * sqrt(fabs(EnVac[levelSetI].value[i])) - 1.0));
            }
			else
			{
                a = 1.0 / dt - expansivity[levelSetI].value[i] + chi * exp(activeEnergy_anion / (number_Boltzmann * temperature) * (xi_anion * sqrt(fabs(EnVac[levelSetI].value[i])) - 1.0));
            }

			double b = 1.0 / dt * preRhoe[levelSetI].value[i] + nonDim_epsilon / number_lambda * EnLiq[levelSetI].value[i];

			nextRhoe[levelSetI].value[i] = b / a;
        }
    }
    else
    {
        double _alpha = (2.0 * dt + preDt) / (dt + preDt);
        double _beta = -dt / (dt + preDt);
#pragma omp parallel for num_threads(thread_num)
        for (int i = 0; i < sharpInterface.geometry().nodeIndexToLabel[levelSetI].size(); i++)
        {     
            double a = 0.0;
            if (EnVac[levelSetI].value[i] > 0.0)
			{
                a = _alpha / dt - expansivity[levelSetI].value[i] + chi * exp(activeEnergy_cation / (number_Boltzmann * temperature) * (xi_cation * sqrt(fabs(EnVac[levelSetI].value[i])) - 1.0));
            }
            else
			{
                a = _alpha / dt - expansivity[levelSetI].value[i] + chi * exp(activeEnergy_anion / (number_Boltzmann * temperature) * (xi_anion * sqrt(fabs(EnVac[levelSetI].value[i])) - 1.0));
            }

			double b = _alpha / dt * preRhoe[levelSetI].value[i] - _beta * (preRhoe[levelSetI].value[i] - prePreRhoe[levelSetI].value[i]) / preDt + nonDim_epsilon / number_lambda * EnLiq[levelSetI].value[i];

			nextRhoe[levelSetI].value[i] = b / a;
        }
    }
}