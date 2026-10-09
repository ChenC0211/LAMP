#pragma once

#include <string>

// ----------------------------------------------------------------------------
// IModule
//
// A lightweight module interface for QUANTUM.
// - Each module owns its config and internal objects.
// - The entrypoint is run() for now.
//
// In later tasks, we can extend to initialize()/step()/finalize(), or pass a
// shared context that holds mesh/geometry/output services.
// ----------------------------------------------------------------------------

class IModule
{
public:
    virtual ~IModule() = default;

    virtual std::string name() const = 0;

    // Return code like main().
    virtual int run() = 0;
};
