LAMP portable (Windows 10 version 1903 or later / Windows 11, x64)

Extract the complete ZIP to a writable local directory. Keep bin/, tools/ and
cases/ beside Start.cmd. Double-click Start.cmd, choose a case, then choose:
  1 - New run
  2 - Resume from a completed checkpoint
  3 - Validate dictionaries without advancing
Enter q or close the input to cancel. Errors are displayed before the window
waits for a key. Visual Studio, CMake and Python are not required.

Edit cases/<case>/system/*.dict before a new run. Results go into a fresh
runs/<case>/<timestamp>/ directory. The executable and its DLLs are copied
there with the dictionaries, so each run retains its runtime snapshot.

For restart, keep the complete old run under this package's runs/<case>/,
including run.json, system/, input/ if present, and output/checkpoints/.
Resume uses frozen dictionaries from that run and writes to a new directory.
Hydro/A VOF checkpoints are supported; B1 restart is not supported. A run with
no completed checkpoint cannot be resumed. No automatic fallback is performed
if the selected newest checkpoint is corrupt.

Command line (PowerShell or cmd):
  .\Start.cmd -Case vof_a_drop_deformation_amr -Action new -NoPause
  .\Start.cmd -Case vof_a_drop_deformation_amr -Action check -NoPause
  .\Start.cmd -Case vof_a_drop_deformation_amr -Action resume -NoPause
  .\Start.cmd -Case vof_a_drop_deformation_amr -Action resume -Restart "runs\vof_a_drop_deformation_amr\<old-run>" -EndTime 0.2 -NoPause
Add -Threads 8 to override the dictionary's thread count.

The package manifest records distributed files. Startup verifies the executable,
DLLs and run launcher. Case dictionary edits are allowed; changing executable
assets requires a new package. Do not copy individual files out of the package.
Only run packages from a trusted source. See LICENSES/ for runtime notices.
