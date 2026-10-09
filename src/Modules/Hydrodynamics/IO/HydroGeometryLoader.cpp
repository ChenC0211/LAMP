#include "Modules/Hydrodynamics/IO/HydroGeometryLoader.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <vector>

#include "Core/IO/DictionaryGeometry.h"
#include "Modules/Hydrodynamics/HydrodynamicsConfig.h"
#include "Modules/Interface/Compatibility/InterfaceNodalConfig.h"
#include "Modules/Interface/Compatibility/InterfacePrescribedData.h"
#include "Modules/Hydrodynamics/IO/HydroInputParser.h"

#include "Core/IO/DictionaryMesh.h"
#include "Core/IO/ConfigurationTrace.h"

namespace
{
void uniqueKey(HydroInput::Reader& input, std::set<std::string>& seen, const std::string& key)
{
    if (!seen.insert(key).second) input.fail("duplicate key/block '" + key + "'");
}

Vector2D coordinate(HydroInput::Reader& input, const std::string& key)
{
    input.expect("(");
    const double x = input.number(key);
    if (input.peek() == ",") input.next();
    const double y = input.number(key);
    input.expect(")");
    return Vector2D(x, y);
}

bool insideBox(const DictionaryGeometry& geometry, const Vector2D& point)
{
    return point.x > geometry.domainMin.x && point.x < geometry.domainMax.x &&
           point.y > geometry.domainMin.y && point.y < geometry.domainMax.y;
}

bool insidePolygon(const std::vector<Vector2D>& polygon, const Vector2D& point)
{
    bool inside = false;
    for (std::size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
    {
        const auto& a = polygon[i];
        const auto& b = polygon[j];
        if ((a.y > point.y) != (b.y > point.y) &&
            point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y) + a.x)
            inside = !inside;
    }
    return inside;
}

bool insideReflectedBox(const DictionaryGeometry& geometry, const Vector2D& point)
{
    return point.x > -geometry.domainMax.x && point.x < geometry.domainMax.x &&
           point.y > geometry.domainMin.y && point.y < geometry.domainMax.y;
}

struct BoundaryExtension
{
    bool west = false, south = false, east = false, north = false;
};

BoundaryExtension validateInterfaceContacts(HydroInput::Reader& input,
    const HydrodynamicsConfig& config, const InterfaceNodalConfig& interface, const DictionaryGeometry& geometry)
{
    BoundaryExtension extension;
    if (interface.boundary.model.empty()) return extension;
    if (interface.boundary.model != "prescribedScalar" || config.materialModel != "twoPhase"
        || config.interfaceForceGeometry != "hermite" || config.interfaceForceLayer != "regularized"
        || interface.reinitializeEvery != 0 || config.boundaryPatches.empty()
        || interface.boundary.fixedContacts.empty())
        input.fail("invalid prescribedScalar interface boundary configuration");
    const auto& contacts = interface.boundary.fixedContacts;
    for (std::size_t i = 0; i < contacts.size(); ++i)
    {
        const auto point = contacts[i];
        const bool west = point.x == geometry.domainMin.x, east = point.x == geometry.domainMax.x;
        const bool south = point.y == geometry.domainMin.y, north = point.y == geometry.domainMax.y;
        if (!std::isfinite(point.x) || !std::isfinite(point.y)
            || point.x < geometry.domainMin.x || point.x > geometry.domainMax.x
            || point.y < geometry.domainMin.y || point.y > geometry.domainMax.y
            || static_cast<int>(west)+east+south+north != 1)
            input.fail("fixedContacts must be finite noncorner points on exactly one physical Domain side");
        if (west && config.coordinateSystem == "axisymmetricRZ")
            input.fail("the radial symmetry axis is not a prescribedScalar fixed contact wall");
        for (std::size_t j = 0; j < i; ++j)
            if (point.x == contacts[j].x && point.y == contacts[j].y) input.fail("duplicate fixed contact");
        const auto side = west ? HydroBoundarySide::West : east ? HydroBoundarySide::East
            : south ? HydroBoundarySide::South : HydroBoundarySide::North;
        const double s = west || east ? point.y : point.x;
        bool constrained = false;
        for (const auto& patch : config.boundaryPatches)
            if (patch.side == side && s >= patch.start && s <= patch.end)
            {
                if (patch.type != HydroBoundaryType::NormalVelocity)
                    input.fail("fixed contact cannot lie on or at the endpoint of a pressure outlet");
                constrained = true;
            }
        if (!constrained) input.fail("every fixed contact must lie on an explicitly declared normalVelocity patch");
        extension.west = extension.west || west; extension.east = extension.east || east;
        extension.south = extension.south || south; extension.north = extension.north || north;
    }
    // Exact grid identity, same-scalar zeros and the complete contact/root
    // inventory are checked by the runtime public geometry, not approximated
    // here from the auxiliary polygon's sampled coordinates.
    return extension;
}

bool insideAuxiliaryBox(const DictionaryGeometry& geometry, const Vector2D& point,
    const BoundaryExtension& extension, bool reflected)
{
    return std::isfinite(point.x) && std::isfinite(point.y)
        && (reflected ? point.x > -geometry.domainMax.x && point.x < geometry.domainMax.x
            : (extension.west || point.x > geometry.domainMin.x)
                && (extension.east || point.x < geometry.domainMax.x))
        && (extension.south || point.y > geometry.domainMin.y)
        && (extension.north || point.y < geometry.domainMax.y);
}

// These analytic checks apply only to the explicitly permitted RZ extension.
// The original planar/off-axis polygon sampling path remains unchanged.
bool arcContainsAngle(const Line& line, double angle)
{
    const double pi = std::acos(-1.0);
    const auto center = line.arcCenter();
    const double begin = std::atan2(line.points[0].y - center.y, line.points[0].x - center.x);
    const double end = begin + line.radian() * pi;
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(begin), std::abs(end)});
    // Conservatively include an angle uncertain only at a sweep endpoint.
    // This can reject an unresolved contact; it never silently removes one.
    for (int turn = -2; turn <= 2; ++turn)
    {
        const double candidate = angle + turn * (2.0 * pi);
        if (candidate >= begin - tolerance && candidate <= end + tolerance) return true;
    }
    return false;
}

struct ContourBounds { double minX, maxX, minY, maxY; };

ContourBounds contourBounds(HydroInput::Reader& input, const Line& line)
{
    ContourBounds bounds{line.points[0].x, line.points[0].x, line.points[0].y, line.points[0].y};
    const auto add = [&](const Vector2D& point)
    {
        if (!std::isfinite(point.x) || !std::isfinite(point.y)) input.fail("unrepresentable reflected contour bounds");
        bounds.minX = std::min(bounds.minX, point.x); bounds.maxX = std::max(bounds.maxX, point.x);
        bounds.minY = std::min(bounds.minY, point.y); bounds.maxY = std::max(bounds.maxY, point.y);
    };
    add(line.points[1]);
    if (line.type() == LineType::arc)
    {
        const auto center = line.arcCenter(); const double radius = line.radius();
        const double pi = std::acos(-1.0);
        const double begin = std::atan2(line.points[0].y - center.y, line.points[0].x - center.x);
        const double end = begin + line.radian() * pi;
        for (double angle : {begin, end})
            add(center + Vector2D(radius * std::cos(angle), radius * std::sin(angle)));
        if (arcContainsAngle(line, 0.0)) add({center.x + radius, center.y});
        if (arcContainsAngle(line, .5 * pi)) add({center.x, center.y + radius});
        if (arcContainsAngle(line, pi)) add({center.x - radius, center.y});
        if (arcContainsAngle(line, 1.5 * pi)) add({center.x, center.y - radius});
    }
    return bounds;
}

bool validateRadialExtension(HydroInput::Reader& input, const HydrodynamicsConfig& config, const InterfaceNodalConfig& interface,
    const DictionaryGeometry& geometry, const std::vector<SignInt>& lineIds,
    const BoundaryExtension& boundaryExtension = {})
{
    if (config.coordinateSystem != "axisymmetricRZ" || geometry.domainMin.x != 0.0 ||
        !((config.interfaceForceGeometry == "hermite" && config.interfaceForceLayer == "regularized")
          || (config.interfaceForceGeometry == "conservative" && config.interfaceForceLayer == "inventory"))) return false;
    std::vector<ContourBounds> bounds;
    bool extension = false;
    for (const auto& point : geometry.points) extension = extension || point.x <= 0.0;
    for (const auto& id : lineIds)
    {
        bounds.push_back(contourBounds(input, geometry.lines[id.value]));
        extension = extension || bounds.back().minX <= 0.0;
    }
    if (!extension) return false;
    if (interface.reinitializeEvery > 0)
        input.fail("axis-reflected auxiliary contours do not support level-set reinitialization");
    if (config.interfaceNormalTraction != 0.0 || config.interfaceTangentialTraction != 0.0)
        input.fail("axis-reflected auxiliary contours do not support legacy prescribed interface traction; use the shared face-load coupling");
    bool negative = false, positive = false;
    for (std::size_t i = 0; i < lineIds.size(); ++i)
    {
        const auto& line = geometry.lines[lineIds[i].value]; const auto& box = bounds[i];
        if (!(box.minX > -geometry.domainMax.x && box.maxX < geometry.domainMax.x &&
              (boundaryExtension.south || box.minY > geometry.domainMin.y)
              && (boundaryExtension.north || box.maxY < geometry.domainMax.y)))
            input.fail("axis-reflected contour must remain strictly inside the reflected radial box and every non-axis wall");
        negative = negative || box.minX < 0.0; positive = positive || box.maxX > 0.0;
        if (line.type() == LineType::straight)
        {
            if (line.points[0].x == 0.0 && line.points[1].x == 0.0)
                input.fail("a straight interface segment along the radial axis is unsupported");
            if (box.minX <= 0.0 && box.maxX >= 0.0 && line.points[0].y != line.points[1].y)
                input.fail("straight axis crossings must be horizontal (smooth orthogonal contact), not conical tips");
        }
        else
        {
            const auto center = line.arcCenter(); const double radius = line.radius();
            bool crosses = line.points[0].x == 0.0 || line.points[1].x == 0.0;
            if (std::abs(center.x) <= radius)
            {
                const double angle = std::acos(-center.x / radius);
                crosses = crosses || arcContainsAngle(line, angle) || arcContainsAngle(line, -angle);
            }
            if (crosses && center.x != 0.0)
                input.fail("arc axis crossings require center.r=0 for smooth orthogonal contact");
        }
    }
    if (!negative || !positive) input.fail("axis contact requires a complete positive/negative-r reflected auxiliary contour");
    const auto mirroredPoint = [](const Vector2D& a, const Vector2D& b)
    { return a.x == -b.x && a.y == b.y; };
    std::vector<std::size_t> mate(lineIds.size(), lineIds.size());
    for (std::size_t i = 0; i < lineIds.size(); ++i)
    {
        const auto& a = geometry.lines[lineIds[i].value];
        for (std::size_t j = 0; j < lineIds.size(); ++j)
        {
            const auto& b = geometry.lines[lineIds[j].value];
            if (a.type() != b.type() || !mirroredPoint(a.points[0], b.points[1]) ||
                !mirroredPoint(a.points[1], b.points[0])) continue;
            if (a.type() == LineType::arc && (!mirroredPoint(a.arcCenter(), b.arcCenter()) ||
                a.radius() != b.radius() || a.radian() != b.radian())) continue;
            if (mate[i] != lineIds.size()) input.fail("reflected contour has duplicate/ambiguous mirror segments");
            mate[i] = j;
        }
        if (mate[i] == lineIds.size()) input.fail("every auxiliary line/arc must have an exact reversed radial mirror mate");
    }
    for (std::size_t i = 0; i < mate.size(); ++i)
        if (mate[mate[i]] != i || mate[(i + 1) % mate.size()] != (mate[i] + mate.size() - 1) % mate.size())
            input.fail("the complete ordered contour is not invariant under reversed radial reflection");
    return true;
}

int loadGeometry(const HydrodynamicsConfig& config, const InterfaceNodalConfig& interface, DictionaryGeometry& geometry,
    const std::string& geometryPath,std::optional<HydroGeometryInput::AxisymmetricColumn>* column,std::optional<HydroGeometryInput::ParametricLoop>* parametric)
{
    HydroInput::Reader input(geometryPath);
    input.expect("Domain");
    input.expect("{");
    std::set<std::string> seen;
    while (input.peek() != "}")
    {
        const auto key = input.next();
        uniqueKey(input, seen, key);
        if (key == "min") geometry.domainMin = coordinate(input, "Domain.min");
        else if (key == "max") geometry.domainMax = coordinate(input, "Domain.max");
        else if (key == "dimension")
        {
            geometry.dimension = input.integer("Domain.dimension");
            geometry.hasDimension = true;
        }
        else input.fail("unsupported Domain key '" + key + "'");
    }
    input.expect("}");
    if (!seen.count("min") || !seen.count("max")) input.fail("Domain requires both min and max");
    if (geometry.dimension != 2) input.fail("Domain.dimension must be 2");
    const double width = geometry.domainMax.x - geometry.domainMin.x;
    const double height = geometry.domainMax.y - geometry.domainMin.y;
    if (!(width > 0.0) || !(height > 0.0) || !std::isfinite(width) || !std::isfinite(height))
        input.fail("Domain bounds must be finite and max greater than min on both axes");
    if (config.coordinateSystem == "axisymmetricRZ" && geometry.domainMin.x != 0.0)
        input.fail("axisymmetricRZ requires Domain.min.x = 0 (the radial axis)");
    geometry.hasDomain = true;
    const auto boundaryExtension = validateInterfaceContacts(input, config, interface, geometry);
    const bool prescribedBoundary = !interface.boundary.model.empty();
    if (config.materialModel == "constant")
    {
        if (!input.done()) input.fail("unsupported geometry block '" + input.peek() +
            "'; constant materialModel requires Domain-only geometry, without Phase, Boundary or LevelSet");
        return -1;
    }

    if(interface.transport=="parametricGals"){
        if(!parametric||config.timeScheme!="cchMomentMidpoint")
            input.fail("parametricGals requires a moment consumer on the physical XY/RZ rectangle");
        input.expect("ParametricInterface");input.expect("{");seen.clear();HydroGeometryInput::ParametricLoop result;std::string name,type;
        const bool rz=config.coordinateSystem=="axisymmetricRZ";
        if(rz&&geometry.domainMin.x!=0)input.fail("parametric spheroid/supply initial data require an RZ domain starting at r=0");
        const auto validateShape=[&](const HydroGeometryInput::ParametricShape& part,const std::string& kind,const std::set<std::string>& keys){
            std::set<std::string> expected{"type","segments"};
            if(keys.count("name"))expected.insert("name");
            if(kind=="supplyCap"){expected.insert("radius");expected.insert("height");}
            else if(kind=="analytic"){
                expected.insert(rz?"r":"x");expected.insert(rz?"z":"y");
                if(keys.count("constants"))expected.insert("constants");
            }
            else if(kind=="implicit"){
                expected.insert({"F","seed"});
                if(keys.count("constants"))expected.insert("constants");
                if(keys.count("traceTolerance"))expected.insert("traceTolerance");
            }
            else {expected.insert({"center","radii"});if(kind=="neckedSpheroid")expected.insert("neckRatio");}
            if(keys!=expected)input.fail("shape keys do not match parametric type '"+kind+"'");
            if(kind=="analytic"||kind=="implicit"){
                if(part.segments<8||part.segments>512)input.fail("analytic shape requires 8..512 segments");
                if(keys.count("traceTolerance")&&!(part.traceTolerance>0))input.fail("traceTolerance must be positive");
                try{
                    const Quantum::Geometry::Parametric::DomainRectangle box{geometry.domainMin.x,geometry.domainMax.x,geometry.domainMin.y,geometry.domainMax.y};
                    const double tolerance=part.traceTolerance>0?part.traceTolerance:1e-9*std::min(box.right-box.left,box.top-box.bottom);
                    const auto curve=kind=="analytic"?InterfacePrescribedData::parametricAnalytic(
                        Quantum::ScalarExpression(part.firstExpression,part.constants),Quantum::ScalarExpression(part.secondExpression,part.constants),part.segments,rz):
                        InterfacePrescribedData::parametricImplicit(Quantum::ScalarExpression(part.implicitExpression,part.constants,
                            rz?Quantum::ScalarExpression::Coordinates::Meridional:Quantum::ScalarExpression::Coordinates::Planar),box,part.seed,part.segments,rz,tolerance);
                    Quantum::Geometry::Parametric::validateContainedSurface(curve,
                        {geometry.domainMin.x,geometry.domainMax.x,geometry.domainMin.y,geometry.domainMax.y},rz);
                }catch(const std::exception& e){input.fail(std::string("invalid analytic shape: ")+e.what());}
                return;
            }
            if(part.segments<8||part.segments>512||!(part.radii.x>0&&part.radii.y>0))input.fail("component requires positive radii and 8..512 segments");
            if(kind=="supplyCap"){
                if(!rz||!keys.count("radius")||!keys.count("height")||part.radii.x>=geometry.domainMax.x||part.radii.y>=geometry.domainMax.y-geometry.domainMin.y)
                    input.fail("supplyCap requires RZ and positive radius/height smaller than the physical domain");
                bool supplied=false;for(const auto& patch:config.boundaryPatches)if(patch.side==HydroBoundarySide::South){
                    if(patch.start==0&&patch.end==part.radii.x&&patch.type==HydroBoundaryType::NormalVelocity&&patch.value<0&&patch.incomingPhase==-1)supplied=true;
                    if(patch.start<=part.radii.x&&patch.end>part.radii.x)input.fail("supplyCap rim must adjoin the unlisted no-slip wall");}
                if(!supplied)input.fail("supplyCap requires a matching south liquid inlet");
            }else if((kind!=(rz?"spheroid":"ellipse")&&!(rz&&kind=="neckedSpheroid"))||!keys.count("center")||!keys.count("radii")||!(rz?(part.center.x==0&&part.radii.x<geometry.domainMax.x&&part.center.y-part.radii.y>geometry.domainMin.y&&part.center.y+part.radii.y<geometry.domainMax.y)
                :(insideBox(geometry,part.center-part.radii)&&insideBox(geometry,part.center+part.radii))))input.fail("component must be an interior ellipse or axis-centered RZ spheroid");
            if(kind=="neckedSpheroid"&&(!keys.count("neckRatio")||!(part.neckRatio>0&&part.neckRatio<2./3)))input.fail("neckedSpheroid requires neckRatio in (0,2/3)");
        };
        while(input.peek()!="}"){const auto key=input.next();uniqueKey(input,seen,key);
            if(key=="name")name=input.next();else if(key=="type")type=input.next();
            else if(key=="center")result.center=coordinate(input,key);else if(key=="radii")result.radii=coordinate(input,key);
            else if(key=="radius")result.radii.x=input.number(key);else if(key=="height")result.radii.y=input.number(key);
            else if(key=="segments")result.segments=input.integer(key);
            else if(key=="neckRatio")result.neckRatio=input.number(key);
            else if(key==(rz?"r":"x"))result.firstExpression=input.next();
            else if(key==(rz?"z":"y"))result.secondExpression=input.next();
            else if(key=="F")result.implicitExpression=input.next();
            else if(key=="seed")result.seed=coordinate(input,key);
            else if(key=="traceTolerance")result.traceTolerance=input.number(key);
            else if(key=="constants"){
                input.expect("{");while(input.peek()!="}"){const auto label=input.next();const double value=input.number(label);
                    if(!result.constants.emplace(label,value).second)input.fail("duplicate analytic constant '"+label+"'");}input.expect("}");
            }
            else if(key=="components"){
                input.expect("{");std::set<std::string> labels;
                while(input.peek()!="}"){const auto label=input.next();uniqueKey(input,labels,label);input.expect("{");
                    HydroGeometryInput::ParametricShape part;std::string kind;std::set<std::string> keys;
                    while(input.peek()!="}"){const auto field=input.next();uniqueKey(input,keys,field);
                        if(field=="type")kind=input.next();else if(field=="center")part.center=coordinate(input,field);else if(field=="radii")part.radii=coordinate(input,field);
                        else if(field=="radius")part.radii.x=input.number(field);else if(field=="height")part.radii.y=input.number(field);
                        else if(field=="neckRatio")part.neckRatio=input.number(field);
                        else if(field==(rz?"r":"x"))part.firstExpression=input.next();
                        else if(field==(rz?"z":"y"))part.secondExpression=input.next();
                        else if(field=="F")part.implicitExpression=input.next();
                        else if(field=="seed")part.seed=coordinate(input,field);
                        else if(field=="traceTolerance")part.traceTolerance=input.number(field);
                        else if(field=="constants"){
                            input.expect("{");while(input.peek()!="}"){const auto constant=input.next();const double value=input.number(constant);
                                if(!part.constants.emplace(constant,value).second)input.fail("duplicate analytic constant '"+constant+"'");}input.expect("}");
                        }
                        else if(field=="segments")part.segments=input.integer(field);else input.fail("unknown parametric component key '"+field+"'");}
                    input.expect("}");
                    validateShape(part,kind,keys);part.suppliedCap=kind=="supplyCap";part.necked=kind=="neckedSpheroid";
                    part.analytic=kind=="analytic";part.implicit=kind=="implicit";result.components.push_back(part);
                }input.expect("}");
            }else input.fail("unknown ParametricInterface key '"+key+"'");}
        input.expect("}");
        if(type=="components"){
            if(!input.done()||seen.size()!=3||name!=config.interfaceLevelSet||result.components.empty())input.fail("components requires name, type components, and a nonempty components block");
            result.segments=0;for(const auto& part:result.components)result.segments=std::max(result.segments,part.segments);*parametric=result;return -1;
        }
        result.suppliedCap=type=="supplyCap";result.necked=type=="neckedSpheroid";
        result.analytic=type=="analytic";
        result.implicit=type=="implicit";
        if(result.analytic||result.implicit){
            if(!input.done()||name!=config.interfaceLevelSet||!seen.count("name"))input.fail("analytic name must match interfaceLevelSet, without trailing geometry blocks");
            validateShape(result,type,seen);*parametric=result;return -1;
        }
        if(result.necked){
            if(!input.done()||seen.size()!=6||name!=config.interfaceLevelSet)input.fail("neckedSpheroid requires name, type, center, radii, neckRatio and segments");
            validateShape(result,type,seen);*parametric=result;return -1;
        }
        if(result.suppliedCap){
            if(!input.done()||seen.size()!=5||name!=config.interfaceLevelSet)
                input.fail("supplyCap requires name, type, radius, height and segments");
            validateShape(result,type,seen);result.center={0,geometry.domainMin.y};
            *parametric=result;return -1;
        }
        if(!input.done()||seen.size()!=5||name!=config.interfaceLevelSet)
            input.fail("ParametricInterface requires name, type, center, radii and segments");
        validateShape(result,type,seen);
        *parametric=result;return -1;
    }
    if(input.peek()=="AxisymmetricColumn") {
        if(!column || config.coordinateSystem!="axisymmetricRZ" || interface.transport!="conservative"
            || !interface.advect || interface.reinitializeEvery!=0 || !interface.boundary.model.empty()
            || !interface.portPhases.empty() || !config.boundaryPatches.empty()
            || config.interfaceForceGeometry!="conservative" || config.interfaceForceLayer!="inventory")
            input.fail("AxisymmetricColumn requires a conservative RZ consumer, advected C, and closed slip boundaries without contacts/ports");
        input.expect("AxisymmetricColumn");input.expect("{");seen.clear();
        HydroGeometryInput::AxisymmetricColumn result;std::string name;
        while(input.peek()!="}") {
            const auto key=input.next();uniqueKey(input,seen,key);
            if(key=="levelSet")name=input.next();
            else if(key=="radius")result.radius=input.number(key);
            else if(key=="relativeAmplitude")result.relativeAmplitude=input.number(key);
            else if(key=="halfWaves")result.halfWaves=input.integer(key);
            else input.fail("unsupported AxisymmetricColumn key '"+key+"'");
        }
        input.expect("}");
        if(!input.done() || !seen.count("levelSet") || !seen.count("radius") || !seen.count("relativeAmplitude")
            || name!=config.interfaceLevelSet || !(result.radius>0) || std::abs(result.relativeAmplitude)>=1
            || !std::isfinite(result.radius) || !std::isfinite(result.relativeAmplitude)
            || (result.halfWaves!=1 && result.halfWaves!=2)
            || !(result.radius*(1+std::abs(result.relativeAmplitude))/std::sqrt(1+result.relativeAmplitude*result.relativeAmplitude/2)<geometry.domainMax.x))
            input.fail("AxisymmetricColumn requires the named interface, finite radius>0, |relativeAmplitude|<1, and radial separation; no auxiliary contour");
        *column=result;return -1;
    }

    const bool mayReflectAxis = config.coordinateSystem == "axisymmetricRZ" && geometry.domainMin.x == 0.0 &&
        ((config.interfaceForceGeometry == "hermite" && config.interfaceForceLayer == "regularized")
        || (config.interfaceForceGeometry == "conservative" && config.interfaceForceLayer == "inventory"));
    input.expect("Point");
    input.expect("{");
    while (input.peek() != "}")
    {
        const auto point = coordinate(input, "Point");
        const bool allowed = prescribedBoundary
            ? insideAuxiliaryBox(geometry, point, boundaryExtension, false)
                || (mayReflectAxis && insideAuxiliaryBox(geometry, point, boundaryExtension, true))
            : insideBox(geometry, point) || (mayReflectAxis && insideReflectedBox(geometry, point));
        if (!allowed && prescribedBoundary)
            input.fail("auxiliary interface points may leave Domain only through a declared fixed-contact side or the checked radial reflection box");
        if (!allowed)
            input.fail(mayReflectAxis
                ? "fluid interface points must be strictly inside Domain or the checked RZ radial reflection box"
                : "fluid interface points must be strictly inside Domain; wall/contact-line interfaces are unsupported");
        geometry.points.push_back(point);
        if (geometry.points.size() > 4096) input.fail("Point count exceeds the supported 4096-point interface limit");
    }
    input.expect("}");
    if (geometry.points.size() < 3) input.fail("a closed fluid interface requires at least three points");

    input.expect("Line");
    input.expect("{");
    std::vector<std::pair<int, int>> endpoints;
    const double pi = std::acos(-1.0);
    const double geometryTolerance = 1e-9 * std::max(width, height);
    while (input.peek() != "}")
    {
        input.expect("(");
        const int first = input.integer("Line point id");
        if (input.peek() == ",") input.next();
        const int second = input.integer("Line point id");
        input.expect(")");
        if (first < 1 || second < 1 || first > static_cast<int>(geometry.points.size()) ||
            second > static_cast<int>(geometry.points.size()) || first == second)
            input.fail("Line requires two distinct valid 1-based Point ids");
        Line line;
        line.points[0] = geometry.points[first - 1];
        line.points[1] = geometry.points[second - 1];
        line.setGeoBoundary(GeometricBoundary::moving);
        line.setImmersed(true);
        if ((line.points[1] - line.points[0]).mod() <= geometryTolerance)
            input.fail("Line endpoints must have nonzero separation");
        input.expect("{");
        std::set<std::string> keys;
        while (input.peek() != "}")
        {
            const auto key = input.next();
            uniqueKey(input, keys, key);
            if (key == "type")
            {
                const auto type = input.next();
                if (type == "straight") line.setType(LineType::straight);
                else if (type == "arc") line.setType(LineType::arc);
                else input.fail("supported fluid Line types are straight and arc (positive sweep only)");
            }
            else if (key == "center") line.setArcCenter(coordinate(input, "Line.center"));
            else if (key == "radius") line.setRadius(input.number("Line.radius"));
            else if (key == "radian") line.setRadian(input.number("Line.radian"));
            else if (key == "name") input.next(); // geometry metadata; binding is by LevelSet name.
            else input.fail("unsupported Line key '" + key + "'; solid boundary assignments are unsupported");
        }
        input.expect("}");
        if (!keys.count("type")) input.fail("Line.type is required");
        if (line.type() == LineType::arc)
        {
            if (!keys.count("center") || !keys.count("radius") || !keys.count("radian"))
                input.fail("arc Line requires center, radius and radian (in multiples of pi)");
            if (!(line.radius() > 0.0) || !(line.radian() > 0.0) || line.radian() >= 2.0)
                input.fail("arc radius must be positive and radian must be in (0,2)");
            const auto center = line.arcCenter();
            const auto start = line.points[0] - center;
            const double angle = std::atan2(start.y, start.x) + line.radian() * pi;
            const auto expected = center + Vector2D(line.radius() * std::cos(angle), line.radius() * std::sin(angle));
            if (std::abs(start.mod() - line.radius()) > geometryTolerance ||
                (expected - line.points[1]).mod() > geometryTolerance)
                input.fail("arc endpoints do not match center/radius/positive radian sweep");
        }
        else if (keys.count("center") || keys.count("radius") || keys.count("radian"))
            input.fail("straight Line must not specify arc parameters");
        endpoints.push_back({first - 1, second - 1});
        geometry.lines.push_back(line);
        if (geometry.lines.size() > 4096) input.fail("Line count exceeds the supported 4096-line interface limit");
    }
    input.expect("}");
    if (geometry.lines.size() < 3) input.fail("a closed fluid interface requires at least three Lines");

    input.expect("LevelSet");
    input.expect("{");
    const std::string name = input.next();
    if (name.empty() || name == "}" || name == "{") input.fail("a named LevelSet entry is required");
    input.expect("{");
    seen.clear();
    std::vector<SignInt> lineIds;
    Vector2D inner;
    bool moving = false;
    while (input.peek() != "}")
    {
        const auto key = input.next();
        uniqueKey(input, seen, key);
        if (key == "lines")
        {
            input.expect("(");
            std::set<int> ids;
            while (input.peek() != ")")
            {
                const int id = input.integer("LevelSet.lines");
                if (id < 1 || id > static_cast<int>(geometry.lines.size()))
                    input.fail("LevelSet.lines requires valid positive 1-based Line ids; legacy zero/signed ids are unsupported");
                if (!ids.insert(id).second) input.fail("LevelSet.lines contains a duplicate Line id");
                lineIds.push_back(SignInt{id - 1, false});
                if (input.peek() == ",") input.next();
            }
            input.expect(")");
        }
        else if (key == "innerPoint") inner = coordinate(input, "LevelSet.innerPoint");
        else if (key == "moving") moving = input.boolean("LevelSet.moving");
        else input.fail("unsupported LevelSet key '" + key + "'");
    }
    input.expect("}");
    if (input.peek() != "}") input.fail("twoPhase currently requires exactly one named LevelSet");
    input.expect("}");
    if (!input.done()) input.fail("unsupported geometry block '" + input.peek() +
        "'; twoPhase supports only Domain/Point/Line/LevelSet, without Phase, Boundary, solid bodies or prescribed motion");
    if (!seen.count("lines") || !seen.count("innerPoint") || !seen.count("moving"))
        input.fail("LevelSet requires lines, innerPoint and moving");
    if (lineIds.size() != geometry.lines.size()) input.fail("the single fluid LevelSet must consume every declared Line");
    if (config.interfaceLevelSet != name) input.fail("unknown interfaceLevelSet '" + config.interfaceLevelSet + "'; available LevelSet: '" + name + "'");
    if (interface.advect && !moving) input.fail("advectInterface requires the selected LevelSet moving true");
    if (interface.reinitializeEvery > 0 && !moving) input.fail("reinitializeEvery requires the selected LevelSet moving true");

    const bool axisExtension = validateRadialExtension(input, config, interface, geometry, lineIds, boundaryExtension);
    std::vector<Vector2D> polygon;
    for (std::size_t i = 0; i < lineIds.size(); ++i)
    {
        const int id = lineIds[i].value;
        const int nextId = lineIds[(i + 1) % lineIds.size()].value;
        if (endpoints[id].second != endpoints[nextId].first)
            input.fail("LevelSet.lines must form one ordered closed contour");
        const auto& line = geometry.lines[id];
        if (prescribedBoundary)
        {
            // Check analytic arc extrema as well as all vertices: permitting
            // auxiliary continuation through a declared wall must not waive
            // any other wall just because polygon samples missed its crossing.
            const auto box = contourBounds(input, line);
            if (!insideAuxiliaryBox(geometry, {box.minX, box.minY}, boundaryExtension, axisExtension)
                || !insideAuxiliaryBox(geometry, {box.maxX, box.maxY}, boundaryExtension, axisExtension))
                input.fail("auxiliary fluid contour crosses an undeclared wall or leaves the allowed radial reflection box");
        }
        const int samples = line.type() == LineType::arc ? 32 : 1;
        const auto center = line.type() == LineType::arc ? line.arcCenter() : Vector2D();
        const double start = line.type() == LineType::arc ? std::atan2(line.points[0].y - center.y, line.points[0].x - center.x) : 0.0;
        for (int sample = 0; sample < samples; ++sample)
        {
            const double t = static_cast<double>(sample) / samples;
            const auto point = line.type() == LineType::arc
                ? center + Vector2D(line.radius() * std::cos(start + t * line.radian() * pi),
                                    line.radius() * std::sin(start + t * line.radian() * pi))
                : line.points[0];
            if (!(prescribedBoundary ? insideAuxiliaryBox(geometry, point, boundaryExtension, axisExtension)
                : axisExtension ? insideReflectedBox(geometry, point) : insideBox(geometry, point)))
                input.fail("fluid contour must stay strictly inside Domain");
            polygon.push_back(point);
        }
    }
    if (!insideBox(geometry, inner) || !insidePolygon(polygon, inner))
        input.fail("LevelSet.innerPoint must lie inside the closed fluid contour");
    geometry.levelSetName.push_back(name);
    geometry.linesOnInterface.push_back(lineIds);
    geometry.innerPoints.push_back(inner);
    geometry.moving.push_back(moving);
    return static_cast<int>(geometry.levelSetName.size()) - 1;
}

// Uniform-mesh preflight remains independent of the geometry parser. Both
// constant and two-phase runs use exactly the same mesh owner and solver path.
void prepareMesh(const DictionaryGeometry& domain, DictionaryMesh& mesh, const std::string& meshPath,
    std::optional<InterfaceAdaptationControls>* adaptation,bool momentRefinement)
{
    const double width = domain.domainMax.x - domain.domainMin.x;
    const double height = domain.domainMax.y - domain.domainMin.y;
    mesh.minCoord = domain.domainMin;
    mesh.maxCoord = domain.domainMax;

    HydroInput::Reader input(meshPath);
    std::set<std::string> blocks;
    bool hasDivide = false;
    bool hasCapacity = false;
    bool adaptive=false;
    mesh.maxInitialIterations = 0;
    while (!input.done())
    {
        const auto block = input.next();
        uniqueKey(input, blocks, block);
        if(block=="RefineBoxes"){
            if(!momentRefinement)input.fail("RefineBoxes currently requires cchMomentMidpoint");
            input.expect("{");std::set<std::string> names;
            while(input.peek()!="}"){
                const auto name=input.next();uniqueKey(input,names,name);input.expect("{");std::set<std::string> keys;DictionaryMesh::RefineBox box;
                const auto point=[&](){input.expect("(");const double x=input.number("coordinate");if(input.peek()==",")input.next();const double y=input.number("coordinate");input.expect(")");return Vector2D{x,y};};
                while(input.peek()!="}"){const auto key=input.next();uniqueKey(input,keys,key);
                    if(key=="min")box.min=point();else if(key=="max")box.max=point();else if(key=="level")box.level=input.integer(key);else input.fail("unknown RefineBoxes key '"+key+"'");}
                input.expect("}");
                if(keys.size()!=3 || box.level<1 || box.level>8 || !(box.min.x<box.max.x&&box.min.y<box.max.y)
                    || box.min.x<domain.domainMin.x || box.min.y<domain.domainMin.y || box.max.x>domain.domainMax.x || box.max.y>domain.domainMax.y)
                    input.fail("RefineBoxes needs min/max inside Domain and level in [1,8]");
                mesh.refineBoxes.push_back(box);if(mesh.refineBoxes.size()>32)input.fail("at most 32 RefineBoxes are supported");
            }
            input.expect("}");continue;
        }
        if (block == "adaptiveMesh")
        {
            adaptive=input.boolean(block);
            if(adaptive && !adaptation) input.fail("adaptiveMesh must be false for this Hydro consumer");
            continue;
        }
        if(block=="Adaptation" && adaptation) {
            // The shared structured reader below owns all entries/duplicates.
            input.expect("{");int depth=1;
            while(depth) {const auto token=input.next();if(token=="{")++depth;else if(token=="}")--depth;}
            continue;
        }
        if (block != "BackgroundMesh" && block != "QuadtreeMesh")
            input.fail("unsupported mesh block '" + block + "'; only a fixed uniform mesh is supported");
        input.expect("{");
        std::set<std::string> keys;
        while (input.peek() != "}")
        {
            const auto key = input.next();
            uniqueKey(input, keys, key);
            if (block == "BackgroundMesh" && key == "divide")
            {
                input.expect("(");
                mesh.Nx = input.integer("divide.Nx");
                if (input.peek() == ",") input.next();
                mesh.Ny = input.integer("divide.Ny");
                input.expect(")");
                hasDivide = true;
            }
            else if (block == "QuadtreeMesh" && key == "maxCellNum")
            {
                mesh.maxCellNum = input.integer(key);
                hasCapacity = true;
            }
            else if (block == "QuadtreeMesh" && key == "maxInitialIterations")
                mesh.maxInitialIterations = input.integer(key);
            else input.fail("unsupported " + block + " key '" + key +
                            "'; legacy domains, level-set rules and refinement are unsupported");
        }
        input.expect("}");
    }
    if (!hasDivide) input.fail("BackgroundMesh.divide is required");
    if (!blocks.count("adaptiveMesh")) input.fail("adaptiveMesh must be explicit");
    if(adaptive!=static_cast<bool>(blocks.count("Adaptation")))input.fail("adaptiveMesh true requires Adaptation; false must omit it");
    if(adaptive) {
        if(hasCapacity)input.fail("adaptive Hydro uses Adaptation.cellBudget; omit redundant QuadtreeMesh.maxCellNum");
        DictionaryInput::Document document(meshPath);DictionaryInput::ConfigurationTrace trace;
        auto controls=readInterfaceAdaptation(document.root().block("Adaptation"),4,trace,momentRefinement);
        if(controls.cellBudget<=static_cast<long long>(mesh.Nx)*mesh.Ny || controls.cellBudget>250000)
            input.fail("Hydro Adaptation.cellBudget must exceed root cells and be at most 250000");
        mesh.maxCellNum=4*controls.cellBudget;mesh.maxLevel_local={controls.policy.maximumLevel};
        *adaptation=std::move(controls);
    } else if(!hasCapacity)input.fail("QuadtreeMesh.maxCellNum is required");
    if (mesh.Nx < 4 || mesh.Ny < 4 || mesh.Nx > 256 || mesh.Ny > 256)
        input.fail("divide dimensions must be integers in [4, 256] for the first-batch solver stencils");
    const long long requiredCapacity = (static_cast<long long>(mesh.Nx) + 1) * (mesh.Ny + 1);
    if (mesh.maxCellNum < requiredCapacity || mesh.maxCellNum > 1000000)
        input.fail("maxCellNum must cover (Nx+1)*(Ny+1) vertices and must not exceed 1000000");
    if (mesh.maxInitialIterations != 0) input.fail("maxInitialIterations must be 0; initialization AMR is unsupported");
    const double dx = width / mesh.Nx;
    const double dy = height / mesh.Ny;
    if (!(dx > 0.0) || !(dy > 0.0) || !std::isfinite(1.0 / dx) || !std::isfinite(1.0 / dy))
        input.fail("Domain/divide produces an invalid mesh spacing");
    mesh.ensureLevelSetCount(static_cast<int>(domain.levelSetName.size()));
}

}

int HydroGeometryInput::load(const HydrodynamicsConfig& config, const InterfaceNodalConfig& interface, DictionaryGeometry& geometry,
    DictionaryMesh& mesh, const std::string& geometryPath, const std::string& meshPath,std::optional<AxisymmetricColumn>* column,
    std::optional<InterfaceAdaptationControls>* adaptation,std::optional<ParametricLoop>* parametric)
{
    if(column)column->reset();
    if(adaptation)adaptation->reset();
    if(parametric)parametric->reset();
    const int interfaceIndex = loadGeometry(config, interface, geometry, geometryPath,column,parametric);
    prepareMesh(geometry, mesh, meshPath,adaptation,config.timeScheme=="cchMomentMidpoint");
    const bool parametricMoment=config.timeScheme=="cchMomentMidpoint"&&parametric&&*parametric;
    if(parametricMoment){
        const double hx=(mesh.maxCoord.x-mesh.minCoord.x)/mesh.Nx,hy=(mesh.maxCoord.y-mesh.minCoord.y)/mesh.Ny;
        if(std::abs(hx-hy)>64*std::numeric_limits<double>::epsilon()*std::max(hx,hy))
            throw std::invalid_argument("parametric moment motion requires square physical leaves, not equal divide counts");
    }
    if(adaptation && *adaptation && ((interface.transport!="conservative" && !(parametric&&*parametric)) || !interface.advect
        || (!parametricMoment&&(!config.boundaryPatches.empty() || !interface.boundary.model.empty()))))
        throw std::invalid_argument("Hydro Adaptation requires an advected conservative interface; boundary patches require parametric CCH");
    if(adaptation&&*adaptation&&(*adaptation)->coarsening&&!parametricMoment)
        throw std::invalid_argument("Hydro Coarsening requires the complete parametric CCH state");
    return interfaceIndex;
}
