#pragma once

struct PicRunSetup;

// Original unit/normalization console and metadata serialization; no solve,
// runtime context construction, accepted state or main clock is owned here.
void writePicRunReports(const PicRunSetup&);
