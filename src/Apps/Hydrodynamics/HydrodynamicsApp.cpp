#include "Hydrodynamics/HydrodynamicsApp.h"
#include "Hydrodynamics/HydroCase.h"
#include <exception>
#include <iostream>

int HydrodynamicsApp::run()
{
    try
    {
        HydroCase flow;
        if (cfg_.checkCase)
        {
            std::cout << "[Hydrodynamics] Case check passed\n";
            return 0;
        }
        int result = flow.initialize();
        if (result != 0) return result;
        while(flow.hasNextStep())
        {
            result = flow.advanceNext();
            if (result != 0) return result;
        }
        flow.finish();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "[Hydrodynamics] ERROR: " << error.what() << '\n';
        return 2;
    }
}
