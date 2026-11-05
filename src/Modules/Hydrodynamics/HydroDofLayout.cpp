#include "Modules/Hydrodynamics/HydroDofLayout.h"

namespace LAMP::Hydrodynamics
{
HydroDofLayout HydroDofLayout::build(const QuadtreeMesh& mesh, const NativeLeafDomain& domain)
{
    domain.validateMeshIdentity(mesh);
    auto data = std::make_shared<Data>();
    data->domain = domain; // Retains Core identity; copies compact selected-cell metadata.
    data->velocityByNativeFace.assign(static_cast<std::size_t>(domain.nativeFaceCount()), -1);
    for (int row = 0; row < domain.nativeFaceCount(); ++row)
    {
        const auto face = domain.face(mesh, row);
        if (face.negativeCell < 0 && face.positiveCell < 0) continue;
        data->velocityByNativeFace[row] = static_cast<int>(data->nativeFaceByVelocity.size());
        data->nativeFaceByVelocity.push_back(row);
        data->normalAxis.push_back(face.horizontal ? NormalAxis::Second : NormalAxis::First);
    }
    return HydroDofLayout(std::move(data));
}
}
