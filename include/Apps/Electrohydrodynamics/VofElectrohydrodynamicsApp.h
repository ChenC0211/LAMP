#pragma once
#include "IApp.h"
#include "RunConfig.h"
#include "Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.h"
// Entry only: the original Hydro run loop owns clock/retry/output lifecycle;
// the VOF-EHD transaction owns the single accepted fluid/charge/field bundle.
class VofElectrohydrodynamicsApp final:public IApp
{
public:
    using Application=VofElectrohydrodynamicsInput::Application;
    explicit VofElectrohydrodynamicsApp(const RunConfig& config,Application application=Application::General)
        :config_(config),application_(application){}
    std::string name() const override{return vofElectrohydrodynamicsApplicationName(application_);}
    int run() override;
private:
    RunConfig config_;
    Application application_;
};
