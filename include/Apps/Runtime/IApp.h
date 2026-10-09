#pragma once

#include <string>

// ----------------------------------------------------------------------------
// IApp
//
// Application-level entrypoint.
//
// Compared to IModule (which represents a reusable physics feature entry), an
// App owns a higher-level workflow: coupling, multi-module orchestration,
// restart strategy, I/O policy, etc.
//
// The interface is intentionally minimal and similar to IModule.
// Later we can extend this to initialize()/step()/finalize() or pass a shared
// SimContext.
// ----------------------------------------------------------------------------

class IApp
{
public:
    virtual ~IApp() = default;

    virtual std::string name() const = 0;

    // Return code like main().
    virtual int run() = 0;
};
