#pragma once

#include "RunConfig.h"

// Application composition root for an already selected case directory.
// Binds built-in workflows and their existing input/preflight policies before
// dispatch through the generic AppRegistry. Does not own a physics time loop.
int runCase(bool checkCase = false, std::optional<int> threadsOverride = {});
