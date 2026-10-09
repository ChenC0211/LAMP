#include "Core/IO/FileManagement.h"
#include <stdexcept>
#include "Modules/Plasma/IO/ImplicitPicFieldCheckpoint.h"
#include "Modules/Plasma/IO/PicCheckpointStore.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Field/CellField.h"
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
namespace ImplicitPicFieldCheckpoint {
namespace {
    struct FieldCheckpointHeader
    {
        char magic[4];
        std::uint32_t version;
        std::uint32_t cellCount;
    };

}
void write(const QuadtreeMesh& checkpointMesh,
    const std::filesystem::path& dir,
    const CellField<double>& phiField,
    const CellField<double>& rhsField,
    const CellField<double>& rhoField,
    const CellField<double>& qField,
    const CellField<Vector2D>& ecellField,
    const CellField<Vector2D>& ecellNewField,
    const CellField<Vector2D>& ehalfField) {
    const std::filesystem::path filePath = dir / PicCheckpointStore::fieldName;
    const int cellCount = checkpointMesh.cellIndexNum();
    const auto requireSize = [&](const auto& field, const char* name) {
        if (field.size() != cellCount)
            throw std::invalid_argument(std::string("[ImplicitPIC][Checkpoint] Field size mismatch: ") + name);
    };
    requireSize(phiField, "phi");
    requireSize(rhsField, "rhs");
    requireSize(rhoField, "rho");
    requireSize(qField, "Q");
    requireSize(ecellField, "Ecell");
    requireSize(ecellNewField, "EcellNew");
    requireSize(ehalfField, "Ehalf");
    auto out = FileIO::openOutputFile(filePath, std::ios::binary | std::ios::out | std::ios::trunc);

    FieldCheckpointHeader header{{'I', 'P', 'F', '1'}, 1u, static_cast<std::uint32_t>(cellCount)};
    out.write(reinterpret_cast<const char*>(&header), sizeof(header));

    auto writeScalarField = [&](const CellField<double>& field) {
        out.write(reinterpret_cast<const char*>(field.value.get()), sizeof(double) * static_cast<std::size_t>(cellCount));
    };
    auto writeVectorField = [&](const CellField<Vector2D>& field) {
        out.write(reinterpret_cast<const char*>(field.value.get()), sizeof(Vector2D) * static_cast<std::size_t>(cellCount));
    };

    writeScalarField(phiField);
    writeScalarField(rhsField);
    writeScalarField(rhoField);
    writeScalarField(qField);
    writeVectorField(ecellField);
    writeVectorField(ecellNewField);
    writeVectorField(ehalfField);
    FileIO::finishOutputFile(out, filePath);
}

bool read(const QuadtreeMesh& checkpointMesh,
    const std::filesystem::path& dir,
    CellField<double>& phiField,
    CellField<double>& rhsField,
    CellField<double>& rhoField,
    CellField<double>& qField,
    CellField<Vector2D>& ecellField,
    CellField<Vector2D>& ecellNewField,
    CellField<Vector2D>& ehalfField) {
    const std::filesystem::path filePath = dir / PicCheckpointStore::fieldName;
    std::ifstream in(filePath, std::ios::binary | std::ios::in);
    if (!in)
        return false;

    FieldCheckpointHeader header{};
    in.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (!in ||
        std::memcmp(header.magic, "IPF1", 4) != 0 ||
        header.version != 1u ||
        header.cellCount != static_cast<std::uint32_t>(checkpointMesh.cellIndexNum()))
    {
        std::cerr << "[ImplicitPIC][Checkpoint] WARNING: incompatible field state checkpoint: "
                  << filePath.string() << std::endl;
        return false;
    }

    phiField.initialize(checkpointMesh, 0.0);
    rhsField.initialize(checkpointMesh, 0.0);
    rhoField.initialize(checkpointMesh, 0.0);
    qField.initialize(checkpointMesh, 0.0);
    ecellField.initialize(checkpointMesh, Vector2D(0.0, 0.0));
    ecellNewField.initialize(checkpointMesh, Vector2D(0.0, 0.0));
    ehalfField.initialize(checkpointMesh, Vector2D(0.0, 0.0));

    const std::size_t fieldSize = static_cast<std::size_t>(checkpointMesh.cellIndexNum());
    auto readScalarField = [&](CellField<double>& field) {
        in.read(reinterpret_cast<char*>(field.value.get()), sizeof(double) * fieldSize);
        return static_cast<bool>(in);
    };
    auto readVectorField = [&](CellField<Vector2D>& field) {
        in.read(reinterpret_cast<char*>(field.value.get()), sizeof(Vector2D) * fieldSize);
        return static_cast<bool>(in);
    };

    if (!readScalarField(phiField) ||
        !readScalarField(rhsField) ||
        !readScalarField(rhoField) ||
        !readScalarField(qField) ||
        !readVectorField(ecellField) ||
        !readVectorField(ecellNewField) ||
        !readVectorField(ehalfField))
    {
        std::cerr << "[ImplicitPIC][Checkpoint] WARNING: failed to read complete field state: "
                  << filePath.string() << std::endl;
        return false;
    }
    return true;
}
}
