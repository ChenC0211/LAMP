#pragma once

#include <string>

// Every registered App owns a complete user workflow: validated input,
// accepted state/clock, module coordination, and output lifecycle. Single
// module and coupled workflows obey this same entry contract. RunCase owns
// selection and process resources; numerical algorithms remain in modules.
class IApp
{
public:
    virtual ~IApp() = default;
    virtual std::string name() const = 0;
    // Return a main()-style status. Preflight support/rejection is explicit.
    virtual int run() = 0;
};
