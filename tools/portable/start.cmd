@echo off
setlocal
set "LAMP_PORTABLE_PAUSE=1"
:scan
if "%~1"=="" goto launch
if /i "%~1"=="-NoPause" set "LAMP_PORTABLE_PAUSE=0"
shift /1
goto scan
:launch
if not exist "%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" (
  echo [PORTABLE] ERROR: Windows PowerShell 5.1 is required.
  set "LAMP_PORTABLE_EXIT=2"
  goto finish
)
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy Bypass -File "%~dp0Start.ps1" %*
set "LAMP_PORTABLE_EXIT=%errorlevel%"
:finish
if "%LAMP_PORTABLE_PAUSE%"=="1" pause
exit /b %LAMP_PORTABLE_EXIT%
