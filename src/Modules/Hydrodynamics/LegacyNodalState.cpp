#include "Modules/Hydrodynamics/LegacyNodalState.h"

#include "Modules/Hydrodynamics/IO/LegacyNodalCache.h"
#include <cstdint>
#include <cstring>

using namespace std;

namespace
{
constexpr char kSharpInterfaceCacheMagic[4] = {'S', 'I', 'C', '2'};
constexpr std::uint32_t kSharpInterfaceCacheVersion = 2;

inline bool MagicEquals(const char a[4], const char b[4])
{
    return std::memcmp(a, b, 4) == 0;
}
} // namespace

LegacyNodalState::LegacyNodalState() {}

LegacyNodalState::~LegacyNodalState() {}

const std::vector<double>& LegacyNodalState::evaporationVelocity(int levelSetI, Time layer) const
{
    if (layer == Time::current) return evaporationVelocity_.at(levelSetI);
    if (layer == Time::previous) return preEvaporationVelocity_.at(levelSetI);
    return prePreEvaporationVelocity_.at(levelSetI);
}

void LegacyNodalState::setEvaporationVelocity(int levelSetI, int nodeLabel, double value)
{
    evaporationVelocity_.at(levelSetI).at(nodeLabel) = value;
}

double LegacyNodalState::noSlipEvaporationVelocity(int nodeLabel) const
{
    return noSlipEvaporationVelocity_.at(nodeLabel);
}

void LegacyNodalState::setNoSlipEvaporationVelocity(int nodeLabel, double value)
{
    noSlipEvaporationVelocity_.resize(geometry_.noSlipNodes().size(), 0.0);
    noSlipEvaporationVelocity_.at(nodeLabel) = value;
}

void LegacyNodalState::output(const std::string& _name, const int levelSetI, const std::vector<string>& valNames, const std::vector<Field<double>>& fields) const
{
    cout << endl << "[SharpInterface] Writing interface data for " << _name << " at level set " << levelSetI << "..." << endl;

    ostringstream fileName;
	fileName << _name;
	string filePath = "./output/" + fileName.str() + ".dat";

	ofstream fp(filePath, ofstream::out);
	if (!fp.is_open())
	{
        cerr << "[SharpInterface] ERROR: failed to write field file: " << filePath << endl;
		exit(-1);
	}

	fp << "Title= \"" << _name << "\"" << endl;
	fp << "Variables= \"X\" \"Y\" \"curvature\" \"normalX\" \"normalY\" \"Ux\" \"Uy\" \"evapVel\"";
    for (int i = 0; i < fields.size(); i++)
    {
        fp << " \"" << valNames[i] << "\"";
    }
    fp << endl;

    for (int i = 0; i < geometry_.nodeIndexToLabel[levelSetI].size(); i++)
    {
        const int nodeLabel = geometry_.nodeIndexToLabel[levelSetI][i];

        fp << geometry_.allNodes(levelSetI)[nodeLabel].coord().x << " " << geometry_.allNodes(levelSetI)[nodeLabel].coord().y << " " << geometry_.allNodes(levelSetI)[nodeLabel].curvature() << " " << geometry_.allNodes(levelSetI)[nodeLabel].normal().x << " " << geometry_.allNodes(levelSetI)[nodeLabel].normal().y << " " << Uxs[levelSetI].value[i] << " " << Uys[levelSetI].value[i] << " " << evaporationVelocity_[levelSetI][nodeLabel] << " ";

        for (int j = 0; j < fields.size(); j++)
        {
            fp << fields[j].value[i] << " ";
        }

        fp << endl;
    }
    fp.close();

    cout << "[SharpInterface] Interface data written." << endl;
}

void LegacyNodalState::temporalOutput(const std::string& _name, const std::vector<std::string>& valNames, const std::vector<Field<double>>& fields, const int levelSetI, const InterfaceNode& node, double nodeEvaporationVelocity, const double& t, const int loop, const bool needInitialize) const
{
    cout << endl << "[SharpInterface] Writing temporal history for " << _name << "..." << endl;

    ostringstream fileName;
	fileName << "time_" << _name;
	string filePath = "./output/" + fileName.str() + ".dat";

	if (needInitialize == true)
	{
		ofstream fp(filePath, ofstream::out);
		if (!fp.is_open())
		{
            cerr << "[SharpInterface] ERROR: failed to write file: " << filePath << endl;
			exit(-1);
		}

		fp << "Title= \"Temporal Node Results" << "\"" << endl;
		fp << "Variables= \"Time\" \"Loop\" \"X\" \"Y\" \"curvature\" \"normalX\" \"normalY\" \"Ux\" \"Uy\" \"evapVel\"";
        for (int i = 0; i < fields.size(); i++)
        {
            fp << " \"" << valNames[i] << "\"";
        }
        fp << endl;

		fp.close();
	}

	ofstream fp(filePath, ofstream::app);
	if (!fp.is_open())
	{
        cerr << "[SharpInterface] ERROR: failed to write field file: " << filePath << endl;
		exit(-1);
	}

    fp << t << " " << loop << " " << node.coord().x << " " << node.coord().y << " " << node.curvature() << " " << node.normal().x << " " << node.normal().y << " " << Uxs[1].value[node.index()] << " " << Uys[1].value[node.index()] << " " << nodeEvaporationVelocity << " ";
    for (int i = 0; i < fields.size(); i++)
    {
        fp << fields[i].value[node.index()] << " ";
    }
    fp << endl;

    fp.close();

    cout << "[SharpInterface] Temporal interface-node history written." << endl;
}

void LegacyNodalState::writeCache(const std::string& filePath) const
{
    cout << "[SharpInterface] Writing interface cache..." << endl;

    ostringstream fileName;
    fileName << filePath << "interface_cache.bin";

    ofstream cache(fileName.str(), ios::binary | ios::out);
    if (cache.fail())
    {
        cout << "Writing sharp interface cache \"" << fileName.str() << "\" failed!" << endl;
        exit(-1);
    }

    cache.write(reinterpret_cast<const char*>(kSharpInterfaceCacheMagic), 4);
    const std::uint32_t version = kSharpInterfaceCacheVersion;
    cache.write(reinterpret_cast<const char*>(&version), sizeof(std::uint32_t));
    const int levelSetCount = static_cast<int>(geometry_.nodeIndexToLabel.size());
    cache.write(reinterpret_cast<const char*>(&levelSetCount), sizeof(int));

    for (int i = 0; i < levelSetCount; i++)
    {
        int nodeIndexNum = static_cast<int>(geometry_.nodeIndexToLabel[i].size());
        cache.write(reinterpret_cast<const char*>(&nodeIndexNum), sizeof(int));
        cache.write(reinterpret_cast<const char*>(geometry_.nodeIndexToLabel[i].data()), nodeIndexNum * sizeof(int));
        int preNodeIndexNum = static_cast<int>(geometry_.preNodeIndexToLabel[i].size());
        cache.write(reinterpret_cast<const char*>(&preNodeIndexNum), sizeof(int));
        cache.write(reinterpret_cast<const char*>(geometry_.preNodeIndexToLabel[i].data()), preNodeIndexNum * sizeof(int));
        int prePreNodeIndexNum = static_cast<int>(geometry_.prePreNodeIndexToLabel[i].size());
        cache.write(reinterpret_cast<const char*>(&prePreNodeIndexNum), sizeof(int));
        cache.write(reinterpret_cast<const char*>(geometry_.prePreNodeIndexToLabel[i].data()), prePreNodeIndexNum * sizeof(int));
        int field_size = static_cast<int>(Uxs[i].size());
        cache.write(reinterpret_cast<const char*>(&field_size), sizeof(int));
        for (int j = 0; j < field_size; j++)
        {
            int label = Uxs[i].label(j);
            cache.write(reinterpret_cast<const char*>(&label), sizeof(int));
            cache.write(reinterpret_cast<const char*>(&Uxs[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&Uys[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&Fxs[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&Fys[i].value[j]), sizeof(double));
        }
        int preField_size = static_cast<int>(preUxs[i].size());
        cache.write(reinterpret_cast<const char*>(&preField_size), sizeof(int));
        for (int j = 0; j < preField_size; j++)
        {
            int label = preUxs[i].label(j);
            cache.write(reinterpret_cast<const char*>(&label), sizeof(int));
            cache.write(reinterpret_cast<const char*>(&preUxs[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&preUys[i].value[j]), sizeof(double));
        }
        int prePreField_size = static_cast<int>(prePreUxs[i].size());
        cache.write(reinterpret_cast<const char*>(&prePreField_size), sizeof(int));
        for (int j = 0; j < prePreField_size; j++)
        {
            int label = prePreUxs[i].label(j);
            cache.write(reinterpret_cast<const char*>(&label), sizeof(int));
            cache.write(reinterpret_cast<const char*>(&prePreUxs[i].value[j]), sizeof(double));
            cache.write(reinterpret_cast<const char*>(&prePreUys[i].value[j]), sizeof(double));
        }

        int nodeNum = static_cast<int>(geometry_.allNodes(i).size());
        cache.write(reinterpret_cast<const char*>(&nodeNum), sizeof(int));
        LegacyNodalCache::writeNodes(cache, geometry_.allNodes(i), evaporationVelocity_[i], nodeNum);
        int preNodeNum = static_cast<int>(geometry_.allPreNodes(i).size());
        cache.write(reinterpret_cast<const char*>(&preNodeNum), sizeof(int));
        LegacyNodalCache::writeNodes(cache, geometry_.allPreNodes(i), preEvaporationVelocity_[i], preNodeNum);
        int prePreNodeNum = static_cast<int>(geometry_.allPrePreNodes(i).size());
        cache.write(reinterpret_cast<const char*>(&prePreNodeNum), sizeof(int));
        LegacyNodalCache::writeNodes(cache, geometry_.allPrePreNodes(i), prePreEvaporationVelocity_[i], prePreNodeNum);
    }

    cache.write(reinterpret_cast<const char*>(&geometry_.noSlipNodeNum), sizeof(int));
    int noSlipFieldSize = static_cast<int>(noSlipUx.size());
    cache.write(reinterpret_cast<const char*>(&noSlipFieldSize), sizeof(int));
    for (int i = 0; i < noSlipFieldSize; i++)
    {
        int label = noSlipUx.label(i);
        cache.write(reinterpret_cast<const char*>(&label), sizeof(int));
        cache.write(reinterpret_cast<const char*>(&noSlipUx.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&noSlipUy.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&noSlipFx.value[i]), sizeof(double));
        cache.write(reinterpret_cast<const char*>(&noSlipFy.value[i]), sizeof(double));
    }
    LegacyNodalCache::writeNodes(cache, geometry_.noSlipNodes(), noSlipEvaporationVelocity_, geometry_.noSlipNodeNum);

    cache.close();
    cout << "[SharpInterface] Interface cache written." << endl;
}

void LegacyNodalState::readCache(const std::string& filePath)
{
    cout << "[SharpInterface] Reading interface cache..." << endl;

    ostringstream fileName;
    fileName << filePath << "interface_cache.bin";

    ifstream cache(fileName.str(), ios::binary | ios::in);
    if (cache.fail())
    {
        cout << "Reading sharp interface cache \"" << fileName.str() << "\" failed!" << endl;
        exit(-1);
    }

    char magic[4] = {0, 0, 0, 0};
    cache.read(magic, 4);
    if (!MagicEquals(magic, kSharpInterfaceCacheMagic))
    {
        cout << "[SharpInterface] ERROR: incompatible cache format in \"" << fileName.str()
             << "\" (delete cache)." << endl;
        exit(-1);
    }
    std::uint32_t version = 0;
    cache.read(reinterpret_cast<char*>(&version), sizeof(std::uint32_t));
    if (version != kSharpInterfaceCacheVersion)
    {
        cout << "[SharpInterface] ERROR: unsupported cache version " << version
             << " in \"" << fileName.str() << "\"." << endl;
        exit(-1);
    }
    int levelSetCount = 0;
    cache.read(reinterpret_cast<char*>(&levelSetCount), sizeof(int));
    if (levelSetCount < 0)
        levelSetCount = 0;
    ensureLevelSetCount(levelSetCount);

    for (int i = 0; i < levelSetCount; i++)
    {
        int nodeIndexNum;
        cache.read(reinterpret_cast<char*>(&nodeIndexNum), sizeof(int));
        geometry_.nodeIndexToLabel[i].resize(nodeIndexNum);
        cache.read(reinterpret_cast<char*>(geometry_.nodeIndexToLabel[i].data()), nodeIndexNum * sizeof(int));
        int preNodeIndexNum;
        cache.read(reinterpret_cast<char*>(&preNodeIndexNum), sizeof(int));
        geometry_.preNodeIndexToLabel[i].resize(preNodeIndexNum);
        cache.read(reinterpret_cast<char*>(geometry_.preNodeIndexToLabel[i].data()), preNodeIndexNum * sizeof(int));
        int prePreNodeIndexNum;
        cache.read(reinterpret_cast<char*>(&prePreNodeIndexNum), sizeof(int));
        geometry_.prePreNodeIndexToLabel[i].resize(prePreNodeIndexNum);
        cache.read(reinterpret_cast<char*>(geometry_.prePreNodeIndexToLabel[i].data()), prePreNodeIndexNum * sizeof(int));

        int field_size;
        cache.read(reinterpret_cast<char*>(&field_size), sizeof(int));
        Uxs[i].initialize(field_size, 0.0);
        Uys[i].initialize(field_size, 0.0);
        Fxs[i].initialize(field_size, 0.0);
        Fys[i].initialize(field_size, 0.0);
        for (int j = 0; j < field_size; j++)
        {
            int label;
            cache.read(reinterpret_cast<char*>(&label), sizeof(int));
            Uxs[i].setLabel(j, label);
            Uys[i].setLabel(j, label);
            Fxs[i].setLabel(j, label);
            Fys[i].setLabel(j, label);
            cache.read(reinterpret_cast<char*>(&Uxs[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&Uys[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&Fxs[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&Fys[i].value[j]), sizeof(double));
        }
        
        int preField_size;
        cache.read(reinterpret_cast<char*>(&preField_size), sizeof(int));
        preUxs[i].initialize(preField_size, 0.0);
        preUys[i].initialize(preField_size, 0.0);
        for (int j = 0; j < preField_size; j++)
        {
            int label;
            cache.read(reinterpret_cast<char*>(&label), sizeof(int));
            preUxs[i].setLabel(j, label);
            preUys[i].setLabel(j, label);
            cache.read(reinterpret_cast<char*>(&preUxs[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&preUys[i].value[j]), sizeof(double));
        }
        int prePreField_size;
        cache.read(reinterpret_cast<char*>(&prePreField_size), sizeof(int));
        prePreUxs[i].initialize(prePreField_size, 0.0);
        prePreUys[i].initialize(prePreField_size, 0.0);
        for (int j = 0; j < prePreField_size; j++)
        {
            int label;
            cache.read(reinterpret_cast<char*>(&label), sizeof(int));
            prePreUxs[i].setLabel(j, label);
            prePreUys[i].setLabel(j, label);
            cache.read(reinterpret_cast<char*>(&prePreUxs[i].value[j]), sizeof(double));
            cache.read(reinterpret_cast<char*>(&prePreUys[i].value[j]), sizeof(double));
        }

        int nodeNum;
        cache.read(reinterpret_cast<char*>(&nodeNum), sizeof(int));
        geometry_.allNodes(i).resize(nodeNum);
        LegacyNodalCache::readNodes(cache, geometry_.allNodes(i), evaporationVelocity_[i], nodeNum);
        int preNodeNum;
        cache.read(reinterpret_cast<char*>(&preNodeNum), sizeof(int));
        geometry_.allPreNodes(i).resize(preNodeNum);
        LegacyNodalCache::readNodes(cache, geometry_.allPreNodes(i), preEvaporationVelocity_[i], preNodeNum);
        int prePreNodeNum;
        cache.read(reinterpret_cast<char*>(&prePreNodeNum), sizeof(int));
        geometry_.allPrePreNodes(i).resize(prePreNodeNum);
        LegacyNodalCache::readNodes(cache, geometry_.allPrePreNodes(i), prePreEvaporationVelocity_[i], prePreNodeNum);
    }

    cache.read(reinterpret_cast<char*>(&geometry_.noSlipNodeNum), sizeof(int));
    int noSlipFieldSize;
    cache.read(reinterpret_cast<char*>(&noSlipFieldSize), sizeof(int));
    noSlipUx.initialize(noSlipFieldSize, 0.0);
    noSlipUy.initialize(noSlipFieldSize, 0.0);
    noSlipFx.initialize(noSlipFieldSize, 0.0);
    noSlipFy.initialize(noSlipFieldSize, 0.0);
    for (int i = 0; i < noSlipFieldSize; i++)
    {
        int label;
        cache.read(reinterpret_cast<char*>(&label), sizeof(int));
        noSlipUx.setLabel(i, label);
        noSlipUy.setLabel(i, label);
        noSlipFx.setLabel(i, label);
        noSlipFy.setLabel(i, label);
        cache.read(reinterpret_cast<char*>(&noSlipUx.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&noSlipUy.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&noSlipFx.value[i]), sizeof(double));
        cache.read(reinterpret_cast<char*>(&noSlipFy.value[i]), sizeof(double));
    }
    geometry_.noSlipNodes().resize(geometry_.noSlipNodeNum);
    LegacyNodalCache::readNodes(cache, geometry_.noSlipNodes(), noSlipEvaporationVelocity_, geometry_.noSlipNodeNum);

    cache.close();
    cout << "[SharpInterface] Interface cache read." << endl;
}
