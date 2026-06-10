@echo off
setlocal
cd /d "%~dp0"

echo ============================================================
echo  SpheroidPA -- build standalone Windows installer
echo ============================================================

:: Step 1: install / upgrade build deps
echo.
echo [1/3] Installing build dependencies ...
pip install pyinstaller scikit-image matplotlib nd2 numpy scipy --quiet
if errorlevel 1 (echo ERROR: pip install failed & pause & exit /b 1)

:: Step 2: PyInstaller -- produces dist\SpheroidPA\
echo.
echo [2/3] Running PyInstaller ...
python -m PyInstaller SpheroidPA.spec --clean --noconfirm
if errorlevel 1 (echo ERROR: PyInstaller failed & pause & exit /b 1)

:: Step 3: Inno Setup (optional -- skipped if iscc.exe not found)
echo.
echo [3/3] Checking for Inno Setup ...
set ISCC=
for %%P in (
    "C:\Program Files (x86)\Inno Setup 6\ISCC.exe"
    "C:\Program Files\Inno Setup 6\ISCC.exe"
) do if exist %%P set ISCC=%%P

if defined ISCC (
    echo Found Inno Setup: %ISCC%
    %ISCC% installer.iss
    if errorlevel 1 (echo ERROR: Inno Setup failed & pause & exit /b 1)
    echo.
    echo Installer written to: installer_out\SpheroidPA_Setup.exe
    echo Copy SpheroidPA_Setup.exe to the NIS-E PC and double-click to install.
) else (
    echo Inno Setup not found -- skipping installer packaging.
    echo.
    echo Distribute manually: copy the entire dist\SpheroidPA\ folder to the NIS-E PC.
    echo Run SpheroidPA\SpheroidPA.exe to launch.
    echo Download Inno Setup from: https://jrsoftware.org/isdl.php
)

echo.
echo ============================================================
echo  Build complete.
echo ============================================================
pause
