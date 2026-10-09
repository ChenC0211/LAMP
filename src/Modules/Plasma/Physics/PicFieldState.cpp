#include "Modules/Plasma/Physics/PicFieldState.h"
#include "Modules/Plasma/IO/ImplicitPicFieldCheckpoint.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include <stdexcept>

PicFieldState::PicFieldState(QuadtreeMesh& mesh)
    : mesh_(mesh), phi_(mesh, 0.0), rhs_(mesh, 0.0), rho_(mesh, 0.0),
      q_(mesh, 0.0), electric_(mesh, Vector2D(0.0, 0.0))
{}

void PicFieldState::copyChargeDensityToRhs()
{
    if (rhs_.size() != mesh_.cellIndexNum() || rho_.size() != mesh_.cellIndexNum())
        throw std::logic_error("PIC charge/RHS fields must be reset after mesh adaptation");
    for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
        rhs_.value[ci] = rho_.value[ci];
}

void PicFieldState::resetAfterMeshChange(const CellField<double>* remappedPotential)
{
    if (remappedPotential)
    {
        if (remappedPotential->size() != mesh_.cellIndexNum())
            throw std::invalid_argument("PIC remapped potential has the wrong cell count");
        for (int ci = 0; ci < mesh_.cellIndexNum(); ++ci)
            if (remappedPotential->label(ci) != mesh_.cellIndexToLabel[ci])
                throw std::invalid_argument("PIC remapped potential has the wrong cell labels");
        phi_ = *remappedPotential;
    }
    else
        phi_.initialize(mesh_, 0.0);
    rhs_.initialize(mesh_, 0.0);
    rho_.initialize(mesh_, 0.0);
    q_.initialize(mesh_, 0.0);
    electric_.initialize(mesh_, Vector2D(0.0, 0.0));
}

ImplicitPicFieldHistory::ImplicitPicFieldHistory(PicFieldState& fields) : fields_(fields)
{
    resetAfterMeshChange();
}

void ImplicitPicFieldHistory::resetAfterMeshChange()
{
    next_.initialize(fields_.mesh_, Vector2D(0.0, 0.0));
    half_.initialize(fields_.mesh_, Vector2D(0.0, 0.0));
}

void ImplicitPicFieldHistory::requireCurrentSizes() const
{
    const int count = fields_.mesh_.cellIndexNum();
    if (fields_.electric_.size() != count || next_.size() != count || half_.size() != count)
        throw std::logic_error("PIC electric history must be rebuilt after mesh adaptation");
}

void ImplicitPicFieldHistory::seedFromCurrent()
{
    requireCurrentSizes();
    for (int ci = 0; ci < fields_.mesh_.cellIndexNum(); ++ci)
    {
        next_.value[ci] = fields_.electric_.value[ci];
        half_.value[ci] = fields_.electric_.value[ci];
    }
}

void ImplicitPicFieldHistory::formHalfStep()
{
    requireCurrentSizes();
    for (int ci = 0; ci < fields_.mesh_.cellIndexNum(); ++ci)
        half_.value[ci] = (fields_.electric_.value[ci] + next_.value[ci]) * 0.5;
}

void ImplicitPicFieldHistory::acceptNext()
{
    requireCurrentSizes();
    for (int ci = 0; ci < fields_.mesh_.cellIndexNum(); ++ci)
        fields_.electric_.value[ci] = next_.value[ci];
}

bool ImplicitPicFieldHistory::readCheckpoint(const std::filesystem::path& directory)
{
    return ImplicitPicFieldCheckpoint::read(fields_.mesh_, directory,
        fields_.phi_, fields_.rhs_, fields_.rho_, fields_.q_, fields_.electric_, next_, half_);
}

void ImplicitPicFieldHistory::writeCheckpoint(const std::filesystem::path& directory) const
{
    ImplicitPicFieldCheckpoint::write(fields_.mesh_, directory,
        fields_.phi_, fields_.rhs_, fields_.rho_, fields_.q_, fields_.electric_, next_, half_);
}
