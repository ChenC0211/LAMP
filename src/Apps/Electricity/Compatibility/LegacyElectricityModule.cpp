#include "Electricity/Compatibility/LegacyElectricityModule.h"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

#include "Electricity/Compatibility/DictionaryElectricity.h"
#include "Core/IO/DictionaryGeometry.h"
#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/OutputLayout.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Geometry/ImmersedGeometry.h"
#include "Modules/Hydrodynamics/LegacyNodalState.h"
#include "Modules/Electricity/Legacy/Electricity.h"
#include "Modules/Hydrodynamics/LegacyNodalForcing.h"

using namespace std;

ElectricityModule::ElectricityModule(const RunConfig& cfg)
    : cfg_(cfg)
{
}

int ElectricityModule::run()
{
    cout << "[Electricity] Init-only run (deterministic, non-interactive)" << endl;

    DictionaryElectricity dictEle;
    if (!dictEle.readDict())
    {
        cerr << "[Electricity] ERROR: cannot open ./system/electricity.dict" << endl;
        return 2;
    }
    const ElectricityConfig& eleCfg = dictEle.config();

    std::filesystem::create_directories("./output");
    std::filesystem::create_directories("./cache");

    DictionaryMesh dictMesh;
    dictMesh.readDict(cfg_.units.get());

    DictionaryGeometry dictGeometry;
    dictGeometry.readDict(cfg_.units.get());
    dictMesh.ensureLevelSetCount(static_cast<int>(dictGeometry.levelSetName.size()));

    QuadtreeMesh mesh(dictMesh);
    ImmersedGeometry<QuadtreeMesh> iGeometry(dictGeometry, mesh);

    // Bring geometry to a consistent state (similar to the legacy init path).
    for (int ls = 0; ls < static_cast<int>(iGeometry.levelSets.size()); ++ls)
    {
        iGeometry.checkLevelSet(mesh, ls);
    }

    iGeometry.markPhases(mesh);
    const double t0 = 0.0;
    iGeometry.reinitializeLevelSets(mesh, t0);
    iGeometry.checkCutInfo(mesh, true);

    LegacyNodalState sharpInterface;
    sharpInterface.initializeGeometry(mesh, iGeometry);
    for (int ls = 0; ls < static_cast<int>(iGeometry.levelSets.size()); ++ls)
        if (iGeometry.levelSets[ls].moving())
            LegacyNodalForcing::initializeCapillary(sharpInterface.geometry().nodeIndexToLabel[ls], sharpInterface.geometry().allNodes(ls),
                number_We_inv, sharpInterface.Fxs[ls], sharpInterface.Fys[ls]);

    Electricity electricity(mesh, sharpInterface, iGeometry);
    electricity.initialize(mesh, sharpInterface, iGeometry);

    if (eleCfg.writeFields)
    {
        // Output vertex fields (Tecplot .dat) for quick verification.
        vector<string> names;
        vector<VertField<double>> fields;

        for (int ls = 0; ls < static_cast<int>(iGeometry.levelSets.size()); ++ls)
        {
            names.push_back("levelSet" + to_string(ls));
            fields.push_back(iGeometry.levelSets[ls]);
        }

        names.push_back("phi");
        names.push_back("Ex");
        names.push_back("Ey");
        names.push_back("res");

        fields.push_back(electricity.phi);
        fields.push_back(electricity.Ex);
        fields.push_back(electricity.Ey);
        fields.push_back(electricity.res);

        OutputLayout::fields();
        mesh.outputFields(eleCfg.outputName, names, fields);
        cout << "[Electricity] Done. Wrote: ./output/" << eleCfg.outputName << ".dat" << endl;
    }
    else
    {
        cout << "[Electricity] Done. (writeFields=false)" << endl;
    }
    return 0;
}
