#pragma once

#include <chrono>
#include <cstdint>
#include <memory>
#include <string>

class QuadtreeMesh;
struct InterfaceTransportInput;
namespace LAMP::Interface { class HermiteGalsPublication; struct HermiteGalsMeshTransferReport; struct HermiteGalsFirstJetStepReport; }

// Read-only serialization of an App-accepted publication. No geometry owner,
// reconstruction, physical acceptance decision or restart state lives here.
class InterfaceFirstJetOutput
{
public:
    using Clock=std::chrono::steady_clock;
    explicit InterfaceFirstJetOutput(const InterfaceTransportInput&);
    ~InterfaceFirstJetOutput();
    void begin();
    void accepted(const QuadtreeMesh&,const LAMP::Interface::HermiteGalsPublication&,
        Clock::time_point start,double stepSeconds,double transferSeconds,
        const LAMP::Interface::HermiteGalsMeshTransferReport* = nullptr,
        const LAMP::Interface::HermiteGalsFirstJetStepReport* = nullptr);
    void complete(long double peakDrift);
    void rejected(int attempted,std::int64_t acceptedStep,double acceptedTime,const std::string& stage,const std::string& reason);
private:
    struct Data;
    std::unique_ptr<Data> data_;
};
