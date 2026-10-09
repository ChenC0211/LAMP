#include "Core/IO/FileManagement.h"
#include "Modules/Plasma/Particles/ParticleType.h"

#include <algorithm>
#include <cstdint>
#include <cstring>

using namespace std;

namespace
{
constexpr char kParticleTypeCacheMagic[4] = {'P', 'T', 'C', '2'};
constexpr std::uint32_t kParticleTypeCacheVersion = 4;

struct ParticleCacheV1Fixed
{
    Vector2D position;
    Vector2D velocity;
    int inCellL;
    int inBackgroundCellL;
    int discreteVertL;
    double weight[4];
    double emitChance;
};

struct ParticleCacheV2Fixed
{
    Vector2D position;
    Vector2D velocity;
    double vz;
    int inCellL;
    int inBackgroundCellL;
    int discreteVertL;
    double weight[4];
    double emitChance;
};

struct ParticleCacheV3Fixed
{
    Vector2D position;
    Vector2D velocity;
    double vz;
    int inCellL;
    int inBackgroundCellL;
    int discreteVertL;
    double w;
    double weight[4];
    double emitChance;
};

struct ParticleCacheV4Fixed
{
    Vector2D position;
    Vector2D velocity;
    double vz;
    int inCellL;
    int inBackgroundCellL;
    int discreteVertL;
    double w;
    double weight[4];
    double emitChance;
    std::uint64_t rngState;
};

inline std::uint32_t levelSetCountFromBytes(std::uint32_t recordBytes, std::uint32_t fixedBytes)
{
    if (recordBytes <= fixedBytes)
        return 0u;
    const std::uint32_t payload = recordBytes - fixedBytes;
    return payload / static_cast<std::uint32_t>(sizeof(double));
}

inline bool MagicEquals(const char a[4], const char b[4])
{
    return std::memcmp(a, b, 4) == 0;
}
} // namespace

ParticleType::ParticleType(const string& name, 
    const double& mass, const double& charge)
:
name_(name),
mass_(mass),
charge_(charge),
q_m_(charge / mass)
{}

ParticleType::~ParticleType() {}

void ParticleType::output(const string& name)
{
	cout << "[ParticleType] Writing particle data: " << this->name_ << "..." << endl;

	ostringstream fileName;
	fileName << "Particle_" << this->name_ << "_" << name;
	string filePath = "./output/" + fileName.str() + ".dat";

	ofstream fp(filePath, ofstream::out);
	if (!fp.is_open())
	{
		cout << "[ParticleType] ERROR: failed to write file: " << filePath << endl;
		exit(-1);
	}

	fp << "Title = \"" << name_ << "\"" << endl;
	fp << "Variables = \"X\" \"Y\" \"Vx\" \"Vy\" \"Vz\""
		<< endl;

	for(Particle& p : particle)
	{
		fp << p.position.x << " " << p.position.y << " " << p.velocity.x << " " << p.velocity.y << " " << p.vz << endl;
	}

	fp.close();
	cout << "[ParticleType] Particle data written." << endl;
}


void ParticleType::writeCache(const std::string& filePath) const
{
	cout << "[ParticleType] Writing particle cache: " << name_ << "..." << endl;

    ostringstream fileName;
	fileName << filePath << "particle_" << name_ << "_cache.bin";

    auto cache = FileIO::openOutputFile(fileName.str(), ios::binary | ios::out);

	// Header (v4): magic + version + levelSetCount.
	cache.write(reinterpret_cast<const char*>(kParticleTypeCacheMagic), 4);
	const std::uint32_t version = kParticleTypeCacheVersion;
	cache.write(reinterpret_cast<const char*>(&version), sizeof(std::uint32_t));
    std::uint32_t levelSetCount = 0;
    for (const Particle& p : particle)
    {
        if (p.levelSetValues.size() > levelSetCount)
            levelSetCount = static_cast<std::uint32_t>(p.levelSetValues.size());
    }
    cache.write(reinterpret_cast<const char*>(&levelSetCount), sizeof(std::uint32_t));

	cache.write(reinterpret_cast<const char*>(&number), sizeof(int));
	
	int nameLen = static_cast<int>(name_.size());
	cache.write(reinterpret_cast<const char*>(&nameLen), sizeof(int));
	cache.write(name_.data(), nameLen);

	cache.write(reinterpret_cast<const char*>(&mass_), sizeof(double));
	cache.write(reinterpret_cast<const char*>(&charge_), sizeof(double));
	cache.write(reinterpret_cast<const char*>(&q_m_), sizeof(double));

	int particleSize = static_cast<int>(particle.size());
	cache.write(reinterpret_cast<const char*>(&particleSize), sizeof(int));

	// Store particle payload as a stable POD layout.
	const std::uint32_t fixedBytes = static_cast<std::uint32_t>(sizeof(ParticleCacheV4Fixed));
	const std::uint32_t particleBytes = fixedBytes + levelSetCount * static_cast<std::uint32_t>(sizeof(double));
	cache.write(reinterpret_cast<const char*>(&particleBytes), sizeof(std::uint32_t));
	for (int i = 0; i < particleSize; ++i)
	{
		const Particle& p = particle[static_cast<std::size_t>(i)];
		ParticleCacheV4Fixed q{};
		q.position = p.position;
		q.velocity = p.velocity;
		q.vz = p.vz;
		q.inCellL = p.inCellL;
		q.inBackgroundCellL = p.inBackgroundCellL;
		q.discreteVertL = p.discreteVertL;
		q.w = p.w;
		std::memcpy(q.weight, p.weight, sizeof(q.weight));
		q.emitChance = p.emitChance;
		q.rngState = p.rngState;
		cache.write(reinterpret_cast<const char*>(&q), sizeof(ParticleCacheV4Fixed));

		if (levelSetCount > 0)
		{
			std::vector<double> values(levelSetCount, 0.0);
			const std::size_t copyCount = std::min(values.size(), p.levelSetValues.size());
			if (copyCount > 0)
				std::copy_n(p.levelSetValues.data(), copyCount, values.data());
			cache.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(double));
		}
	}

	// Optional tail payload (kappa_s).
	cache.write(reinterpret_cast<const char*>(&kappa_s_), sizeof(double));

	FileIO::finishOutputFile(cache, fileName.str());
	cout << "[ParticleType] Particle cache written: " << name_ << "." << endl;
}

void ParticleType::readCache(const std::string& filePath)
{
	cout << "[ParticleType] Reading particle cache: " << name_ << "..." << endl;

	ostringstream fileName;
	fileName << filePath << "particle_" << name_ << "_cache.bin";

	ifstream cache(fileName.str(), ios::binary | ios::in);
	if (cache.fail())
	{
		cout << "Reading cache \"" << fileName.str() << "\" failed!" << endl;
		exit(-1);
	}

	// Versioned cache format (v2/v3/v4) starts with a magic.
	char magic[4] = {0, 0, 0, 0};
	cache.read(magic, 4);
	if (cache.fail())
	{
		cout << "Reading cache \"" << fileName.str() << "\" failed (empty?)!" << endl;
		exit(-1);
	}

	if (!MagicEquals(magic, kParticleTypeCacheMagic))
	{
		cout << "[ParticleType] ERROR: unsupported legacy cache format for \"" << fileName.str()
			 << "\" (delete cache)." << endl;
		exit(-1);
	}

	std::uint32_t version = 0;
	cache.read(reinterpret_cast<char*>(&version), sizeof(std::uint32_t));
	if (!(version == 4u || version == 3u || version == 2u))
	{
		cout << "[ParticleType] ERROR: unsupported cache version " << version
			 << " for \"" << fileName.str() << "\"" << endl;
		exit(-1);
	}

	std::uint32_t levelSetCount = 0;
	if (version >= 4u)
	{
		cache.read(reinterpret_cast<char*>(&levelSetCount), sizeof(std::uint32_t));
	}

	cache.read(reinterpret_cast<char*>(&number), sizeof(int));
	int nameLen = 0;
	cache.read(reinterpret_cast<char*>(&nameLen), sizeof(int));
	name_.resize(nameLen);
	cache.read(&name_[0], nameLen);

	cache.read(reinterpret_cast<char*>(&mass_), sizeof(double));
	cache.read(reinterpret_cast<char*>(&charge_), sizeof(double));
	cache.read(reinterpret_cast<char*>(&q_m_), sizeof(double));

	int particleSize = 0;
	cache.read(reinterpret_cast<char*>(&particleSize), sizeof(int));
	std::uint32_t particleBytes = 0;
	cache.read(reinterpret_cast<char*>(&particleBytes), sizeof(std::uint32_t));

	if (particleSize < 0)
		particleSize = 0;
	particle.clear();
	particle.resize(static_cast<std::size_t>(particleSize));

	if (version >= 4u)
	{
		const std::uint32_t fixedBytes = static_cast<std::uint32_t>(sizeof(ParticleCacheV4Fixed));
		if (particleBytes < fixedBytes)
		{
			cout << "[ParticleType] ERROR: unexpected particle record size " << particleBytes
				 << " in cache \"" << fileName.str() << "\"" << endl;
			exit(-1);
		}
		const std::uint32_t valuesInRecord = levelSetCountFromBytes(particleBytes, fixedBytes);
		const std::uint32_t readValues = (levelSetCount > 0) ? std::min(valuesInRecord, levelSetCount) : valuesInRecord;
		const std::uint32_t skipBytes = particleBytes - fixedBytes - readValues * static_cast<std::uint32_t>(sizeof(double));

		for (int i = 0; i < particleSize; ++i)
		{
			ParticleCacheV4Fixed q{};
			cache.read(reinterpret_cast<char*>(&q), sizeof(ParticleCacheV4Fixed));

			Particle& p = particle[static_cast<std::size_t>(i)];
			p.position = q.position;
			p.velocity = q.velocity;
			p.vz = q.vz;
			p.inCellL = q.inCellL;
			p.inBackgroundCellL = q.inBackgroundCellL;
			p.discreteVertL = q.discreteVertL;
			p.w = q.w;
			std::memcpy(p.weight, q.weight, sizeof(q.weight));
			p.emitChance = q.emitChance;
			p.rngState = q.rngState;
			p.ensureLevelSetSize(levelSetCount);

			if (readValues > 0)
			{
				std::vector<double> values(readValues, 0.0);
				cache.read(reinterpret_cast<char*>(values.data()), values.size() * sizeof(double));
				const std::size_t copyCount = std::min<std::size_t>(p.levelSetValues.size(), values.size());
				if (copyCount > 0)
					std::copy_n(values.data(), copyCount, p.levelSetValues.data());
			}
			if (skipBytes > 0)
				cache.ignore(static_cast<std::streamsize>(skipBytes));
		}
	}
	else
	{
		const std::uint32_t fixedV3 = static_cast<std::uint32_t>(sizeof(ParticleCacheV3Fixed));
		const std::uint32_t fixedV2 = static_cast<std::uint32_t>(sizeof(ParticleCacheV2Fixed));
		const std::uint32_t fixedV1 = static_cast<std::uint32_t>(sizeof(ParticleCacheV1Fixed));

		enum class LegacyLayout { V3, V2, V1 };
		LegacyLayout layout = LegacyLayout::V1;
		if (version == 3u && particleBytes >= fixedV3)
			layout = LegacyLayout::V3;
		else if (particleBytes >= fixedV2)
			layout = LegacyLayout::V2;
		else if (particleBytes >= fixedV1)
			layout = LegacyLayout::V1;
		else
		{
			cout << "[ParticleType] ERROR: unexpected particle record size " << particleBytes
				 << " in cache \"" << fileName.str() << "\"" << endl;
			exit(-1);
		}

		std::uint32_t fixedBytes = 0;
		switch (layout)
		{
		case LegacyLayout::V3: fixedBytes = fixedV3; break;
		case LegacyLayout::V2: fixedBytes = fixedV2; break;
		case LegacyLayout::V1: fixedBytes = fixedV1; break;
		}
		const std::uint32_t legacyCount = levelSetCountFromBytes(particleBytes, fixedBytes);
		levelSetCount = legacyCount;

		for (int i = 0; i < particleSize; ++i)
		{
			Particle& p = particle[static_cast<std::size_t>(i)];
			p.ensureLevelSetSize(levelSetCount);
			p.rngState = 0;

			if (layout == LegacyLayout::V3)
			{
				ParticleCacheV3Fixed q{};
				cache.read(reinterpret_cast<char*>(&q), sizeof(ParticleCacheV3Fixed));
				p.position = q.position;
				p.velocity = q.velocity;
				p.vz = q.vz;
				p.inCellL = q.inCellL;
				p.inBackgroundCellL = q.inBackgroundCellL;
				p.discreteVertL = q.discreteVertL;
				p.w = q.w;
				std::memcpy(p.weight, q.weight, sizeof(q.weight));
				p.emitChance = q.emitChance;
			}
			else if (layout == LegacyLayout::V2)
			{
				ParticleCacheV2Fixed q{};
				cache.read(reinterpret_cast<char*>(&q), sizeof(ParticleCacheV2Fixed));
				p.position = q.position;
				p.velocity = q.velocity;
				p.vz = q.vz;
				p.inCellL = q.inCellL;
				p.inBackgroundCellL = q.inBackgroundCellL;
				p.discreteVertL = q.discreteVertL;
				p.w = 1.0;
				std::memcpy(p.weight, q.weight, sizeof(q.weight));
				p.emitChance = q.emitChance;
			}
			else
			{
				ParticleCacheV1Fixed q{};
				cache.read(reinterpret_cast<char*>(&q), sizeof(ParticleCacheV1Fixed));
				p.position = q.position;
				p.velocity = q.velocity;
				p.vz = 0.0;
				p.inCellL = q.inCellL;
				p.inBackgroundCellL = q.inBackgroundCellL;
				p.discreteVertL = q.discreteVertL;
				p.w = 1.0;
				std::memcpy(p.weight, q.weight, sizeof(q.weight));
				p.emitChance = q.emitChance;
			}

			if (legacyCount > 0)
			{
				std::vector<double> values(legacyCount, 0.0);
				cache.read(reinterpret_cast<char*>(values.data()), values.size() * sizeof(double));
				const std::size_t copyCount = std::min<std::size_t>(p.levelSetValues.size(), values.size());
				if (copyCount > 0)
					std::copy_n(values.data(), copyCount, p.levelSetValues.data());
			}
		}
	}

	// Optional tail payload (kappa_s).
	if (cache.peek() != EOF)
	{
		cache.read(reinterpret_cast<char*>(&kappa_s_), sizeof(double));
		if (cache.fail())
			kappa_s_ = 0.0;
	}

	cache.close();
	cout << "[ParticleType] Particle cache read: " << name_ << "." << endl;
}
