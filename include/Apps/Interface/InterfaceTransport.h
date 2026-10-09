#pragma once

#include "Core/Base/ScalarExpression.h"
#include "Core/IO/ConfigurationTrace.h"
#include "Core/Geometry/PolynomialVelocitySlab.h"
#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

// Standalone prescribed transport input. Problem names, exact references and
// benchmark acceptance thresholds never belong in this production contract.
struct InterfaceTransportConfig
{
    enum class Field { Resolution, Steps, Dt, EndTime, WriteEverySteps, RegridSteps, Domain, Refinement };
    enum class Velocity { Constant, Polynomial };
    bool radial=false, initialAmr=false, writeFields=false, prescribedBoundary=false;
    int resolution=0, steps=0, writeEverySteps=1, refinementLevel=1;
    double dt=0, endTime=0;
    Vector2D domainMin{0,0}, domainMax{1,1}, refineMin{}, refineMax{};
    std::vector<int> regridSteps;
    Quantum::ScalarExpression::Constants constants;
    std::string initialPhi, boundaryPhi;
    Velocity velocity=Velocity::Constant;
    Vector2D constantVelocity{}, velocityOrigin{}, velocityScale{1,1};
    struct Coefficient { std::string x="0", y="0"; bool specified=false; };
    std::array<Coefficient,16> coefficients; // [4*i+j] * xi^i eta^j
    void validate() const;
};
struct InterfaceTransportConfigError : std::invalid_argument
{
    InterfaceTransportConfigError(InterfaceTransportConfig::Field,const char*);
    const InterfaceTransportConfig::Field field;
};
// Expressions own immutable compiled functions. Time coefficients feed the
// existing polynomial endpoint provider, including its gradient/CFL contract.
class InterfacePrescribedData
{
public:
    explicit InterfacePrescribedData(const InterfaceTransportConfig&);
    Quantum::ScalarExpression::Jet initial(Vector2D) const;
    Quantum::ScalarExpression::Jet boundary(Vector2D,double time) const;
    Quantum::Geometry::PolynomialVelocityData polynomial(double previous,double arrival) const;
private:
    Quantum::ScalarExpression initial_;
    std::shared_ptr<const Quantum::ScalarExpression> boundary_;
    Vector2D origin_,scale_;
    std::array<std::array<std::shared_ptr<const Quantum::ScalarExpression>,2>,16> coefficients_;
};
struct InterfaceTransportInput
{
    InterfaceTransportConfig config;
    DictionaryInput::ConfigurationTrace trace;
};
// The standalone workflow consumes exactly these files. IO uses the same list
// for input provenance; numerical modules never see paths or dictionaries.
inline constexpr std::array<const char*,5> interfaceTransportInputFiles{
    "system/run.dict","system/geometry.dict","system/mesh.dict","system/fields.dict","system/interface.dict"};
InterfaceTransportInput readInterfaceTransportInput();
int runInterfaceTransport(const InterfaceTransportInput&,bool checkOnly);
