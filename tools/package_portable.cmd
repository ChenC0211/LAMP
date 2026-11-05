@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0package_portable.ps1" %*
exit /b %errorlevel%
