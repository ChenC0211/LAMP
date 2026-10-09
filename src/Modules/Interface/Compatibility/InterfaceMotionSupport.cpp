#include "Modules/Interface/Compatibility/Detail/InterfaceMotionSupport.h"
#include "Modules/Interface/Compatibility/InterfaceTransportInput.h"
#include "Core/Geometry/LevelSetInterfacePanels.h"
#include "Core/Mesh/QuadtreeMesh.h"
#include "Core/Mesh/UniformLeafGrid.h"
#include "Core/Mesh/NativeLeafDomain.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace {
void require(bool valid, const char* reason)
{
    if (!valid) throw std::runtime_error(std::string("[InterfaceMotion] ") + reason);
}
void checkBand(const QuadtreeMesh &mesh, const LevelSetHermiteGeometry::Snapshot &geometry,
               const std::set<int> &cells, double width, bool inside, const char *stage)
{
    const std::vector<int> ordered(cells.begin(),cells.end());
    std::size_t first=ordered.size();
    if (geometry.usesNativeLeaves())
    {
        using Requirement=LevelSetHermiteGeometry::Snapshot::BandRequirement;
        first=geometry.firstNativeBandViolation(ordered,width,
            inside?Requirement::EntireCellInside:Requirement::EntireCellOutside);
    }
    else for (std::size_t i=0;i<ordered.size();++i)
    {
        const auto point=mesh.cells(mesh.cellIndexToLabel[ordered[i]]).centerCoord();
        if (inside?geometry.closedCellInsideBand(point,width):geometry.closedCellOutsideBand(point,width)) continue;
        first=i;break;
    }
    if (first==ordered.size()) return;
    const int k=ordered[first];
    const auto point=mesh.cells(mesh.cellIndexToLabel[k]).centerCoord();
    std::ostringstream message;
    message << std::setprecision(17) << "[InterfaceMotion] " << stage << ": cell=" << k
            << ", center=(" << point.x << ',' << point.y << "), phiWidth=" << width << ", required="
            << (inside ? "entire closed cell inside full-speed band" : "entire closed cell outside cutoff halo");
    throw std::runtime_error(message.str());
}
} // namespace

InterfaceMotionSupport::Band::Band(double requiredHalfWidthCells, bool quartic, double maxh)
{
    // Recursive-extension data need up to seven cells on a one-sided chart.
    // Enlarge only the numerical full-speed/support band, not force thickness.
    const double supportHalo = quartic ? 7.0 : 4.0;
    fullWidth = (requiredHalfWidthCells + supportHalo) * maxh;
    zeroWidth = fullWidth + 4.0 * maxh;
    portWidth = zeroWidth + 2.0 * maxh;
    require(std::isfinite(portWidth) && portWidth > zeroWidth, "transport band is unresolved");
}

InterfaceMotionSupport::InterfaceMotionSupport(const QuadtreeMesh& mesh,
    const InterfaceTransportDomain& transport, const LevelSetHermiteGeometry::Snapshot& source,
    const UniformLeafGrid* grid, const std::vector<std::array<int, 2>>& logical, double maxh, Band band)
    : mesh_(mesh), domain_(transport.geometry()), source_(source), grid_(grid), logical_(logical),
      maxh_(maxh), supportHalo_(grid ? 4.0 : 7.0), band_(band)
{
    for (const auto &face : domain_.faces(mesh_))
        if (transport.boundary(face.nativeRow) == InterfaceTransportDomain::Boundary::SinglePhasePort)
        {
            const int k = face.negativeCell >= 0 ? face.negativeCell : face.positiveCell;
            require(k >= 0 && k < mesh_.cellIndexNum(), "pressure face has no active owner");
            pressureCells_.insert(k);
        }
    if (!grid_)
        locator_.emplace(mesh_);
}

void InterfaceMotionSupport::append(const LevelSetInterfacePanels& panels)
{
    for (const auto &panel : panels.panels())
    {
        require(panel.cellIndex >= 0 && panel.cellIndex < mesh_.cellIndexNum() &&
                    domain_.cellRowByIndex()[panel.cellIndex] >= 0,
                "free panel lacks an active cell");
        if (!grid_)
        {
            const auto bounds = source_.nativeCellBounds(panel.cellIndex);
            const auto h = bounds[1] - bounds[0];
            const Vector2D lo{std::max(mesh_.minCoord.x, bounds[0].x - supportHalo_ * h.x),
                              std::max(mesh_.minCoord.y, bounds[0].y - supportHalo_ * h.y)};
            const Vector2D hi{std::min(mesh_.maxCoord.x, bounds[1].x + supportHalo_ * h.x),
                              std::min(mesh_.maxCoord.y, bounds[1].y + supportHalo_ * h.y)};
            for (int k : locator_->intersectingBox(lo, hi))
            {
                const auto candidate = source_.nativeCellBounds(k);
                if (domain_.cellRowByIndex()[k] >= 0 &&
                    std::max(lo.x, candidate[0].x) < std::min(hi.x, candidate[1].x) &&
                    std::max(lo.y, candidate[0].y) < std::min(hi.y, candidate[1].y))
                    supportCells_.insert(k);
            }
            continue;
        }
        const auto ij = logical_[panel.cellIndex];
        // Conservative four-cell Cartesian halo: covers both central and
        // one-sided 5x5 supports. Missing solid cells are not data/extension.
        for (int j = std::max(0, ij[1] - 4); j <= std::min(grid_->ny() - 1, ij[1] + 4); ++j)
            for (int i = std::max(0, ij[0] - 4); i <= std::min(grid_->nx() - 1, ij[0] + 4); ++i)
            {
                const int k =
                    mesh_.cells(grid_->cells()[static_cast<std::size_t>(j) * grid_->nx() + i]).index();
                if (domain_.cellRowByIndex()[k] >= 0)
                    supportCells_.insert(k);
            }
    }
}

void InterfaceMotionSupport::prepare(const LevelSetInterfacePanels& previous)
{
    append(previous);
    if (!grid_)
    {
        // phi is transported data, not necessarily a signed distance. A
        // geometric halo's diagonal reach cannot be converted to a scalar
        // amplitude by its cell count alone. Select one OLD-scalar band by
        // the existing full-cell certificate before constructing any trace.
        // Preserve successful powers-of-two choices. If a growth step reaches
        // a pressure port, search its bracket using these SAME certificates;
        // doubling alone can skip a feasible band. No transport trial or
        // residual controls this selection, and no port margin is reduced.
        using Requirement=LevelSetHermiteGeometry::Snapshot::BandRequirement;
        const std::vector<int> orderedSupport(supportCells_.begin(),supportCells_.end());
        const std::vector<int> orderedPressure(pressureCells_.begin(),pressureCells_.end());
        const auto covered=[&](double width) {
            return source_.firstNativeBandViolation(orderedSupport,width,
                Requirement::EntireCellInside)==orderedSupport.size();
        };
        const auto separated=[&](const Band& band) {
            return source_.firstNativeBandViolation(orderedPressure,band.portWidth,
                Requirement::EntireCellOutside)==orderedPressure.size();
        };
        const auto atWidth=[&](double width) {
            Band band=band_;
            band.fullWidth=width;
            band.zeroWidth=band.fullWidth+4.0*maxh_;
            band.portWidth=band.zeroWidth+2.0*maxh_;
            require(std::isfinite(band.portWidth) && band.portWidth>band.zeroWidth
                && band.zeroWidth>band.fullWidth,"required scalar transport band is unresolved");
            return band;
        };
        double lower=0.0;
        for (;;)
        {
            if (!separated(band_))
            {
                // lower is port-safe but support-incomplete; upper is not
                // port-safe. Monotone full-cell certificates preserve this
                // bracket. Stop at adjacent representable widths, never at
                // an arbitrary tolerance that could hide a narrow interval.
                double upper=band_.fullWidth;
                while (lower>0.0)
                {
                    const double middle=0.5*lower+0.5*upper;
                    if (!(middle>lower && middle<upper)) break;
                    const auto trial=atWidth(middle);
                    if (!separated(trial)) upper=middle;
                    else if (covered(middle)) { band_=trial; break; }
                    else lower=middle;
                }
                // An initial conflict cannot shrink below the caller's
                // minimum band. An unresolved bracket retains the original
                // failed growth band for the existing precise diagnostic.
                checkBand(mesh_, source_, pressureCells_, band_.portWidth, false,
                    "required support/pressure-port separation");
                break;
            }
            if (covered(band_.fullWidth)) break;
            lower=band_.fullWidth;
            band_=atWidth(2.0*band_.fullWidth);
        }
    }
    checkBand(mesh_, source_, pressureCells_, band_.portWidth, false, "old pressure-port separation");
    checkBand(mesh_, source_, supportCells_, band_.fullWidth, true, "old interface support");
}

void InterfaceMotionSupport::checkCandidate(const LevelSetInterfacePanels& panels,
    const LevelSetHermiteGeometry::Snapshot& geometry, Stage stage)
{
    append(panels);
    const bool corrected = stage == Stage::Corrected;
    checkBand(mesh_, source_, supportCells_, band_.fullWidth, true,
        corrected ? "corrected support in old full-speed band" : "new support in old full-speed band");
    checkBand(mesh_, geometry, supportCells_, band_.fullWidth, true,
        corrected ? "corrected interface support" : "new interface support");
    checkBand(mesh_, geometry, pressureCells_, band_.portWidth, false,
        corrected ? "corrected pressure-port separation" : "new pressure-port separation");
}
