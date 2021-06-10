::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: CHAMP (CHerednik Algebra Magma Package)
:: Copyright (C) 2013–2021 Ulrich Thiel
:: Licensed under GNU GPLv3, see COPYING.
:: http://ulthiel.com/math/champ
:: thiel@mathematik.uni-kl.de
::
:: Start CHAMP (and automatically set the necessary environment variables).
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

::prevent output
@echo off

::otherwise variables will exist globally
setlocal

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: You can define the Magma installation directory here if Magma
:: cannot be found or you want to use some specific installation.
:: Don't put a backslash at the end of the name, and no spaces around the
:: equal sign.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set MAGMA_DIR="..."


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Find Magma (first user defined path, then environment path, then trying
:: some directories, then give up).
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

::first, use path provided by user
if not %MAGMA_DIR% == "..." (
    if not exist %MAGMA_DIR%\magma.exe (
        goto startchamp
    ) else (
    	set MAGMA_EXEC=%MAGMA_DIR%\magma.exe
        goto startchamp
    )
)

::now, check if magma is in path
where magma >NUL 2>&1
if %ERRORLEVEL% equ 0 (
    set MAGMA_EXEC="magma"
    goto startchamp
)

::now, check the program files paths
if exist "%ProgramFiles%\Magma\magma.exe" (
    set MAGMA_EXEC="%ProgramFiles%\Magma\magma.exe"
    goto startchamp
)

if exist "%ProgramFiles(x86)%\Magma\magma.exe" (
    set MAGMA_EXEC="%ProgramFiles(x86)%\Magma\magma.exe"
    goto startchamp
)

if exist "%ProgramW6432%\Magma\magma.exe" (
    set MAGMA_EXEC="%ProgramW6432%\Magma\magma.exe"
    goto startchamp
)

::at last, try to read Magma installation directory from registry
::(i'm not sure if this is always the correct key though)
setlocal ENABLEEXTENSIONS
set KEY_NAME="HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\Windows\CurrentVersion\Uninstall\Magma_is1"
set VALUE_NAME=InstallLocation

FOR /F "usebackq skip=2 tokens=2,*" %%A IN (`REG QUERY %KEY_NAME% /v %VALUE_NAME% 2^>nul`) DO (
    set MAGMA_DIR="%%B"
)

if defined MAGMA_DIR (
	if exist %MAGMA_DIR%magma.exe (
	    set MAGMA_EXEC=%MAGMA_DIR%magma.exe
	    goto startchamp
	)

)

:nomagma
echo Error: cannot find Magma.
echo Either add Magma installation directory in champ.bat or in PATH environment variable.
pause
exit /b


:startchamp
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:: Set up all the environment variables for CHAMP
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

:: current directory. in contrast to the unix version this comes with a
:: trailing backslash. i therefore add dot at the end.
set CHAMP_DIR=%~dp0.

:: no quotation marks here
set CHAMP_OS_TYPE=Windows

:: the value of %OS% is Windows_NT since Windows XP i think.
:: uname.exe -s outputs WindowsNT and I'll just go for this.
:: i don't care about any older Windows any more but i think this will also not
:: get relevant.
set CHAMP_OS=WindowsNT

:: Get CHAMP version. I'll cd into the CHAMP directory to avoid git -C since
:: this didn't exist before Git 1.8.5 (Nov 2013).
:: Setting the command output to a variable is a bit ugly in Windows.
set curdir=%cd%
cd %CHAMP_DIR%
git describe >NUL 2>NUL
if %ERRORLEVEL% equ 0 (
    FOR /F "tokens=* USEBACKQ" %%F IN (`git describe`) DO (
        SET CHAMP_VER=%%F
    )
) else (
    if exist "version.txt" (
        set /p CHAMP_VER= < version.txt
    )
)
cd %curdir%

:: Get the OS version (more detailed than CHAMP_OS)
for /f "delims== tokens=2-3" %%f in ('wmic os get Version /value ^| find "="') do (
    set CHAMP_OS_VER=Microsoft Windows %%f
)

:: Hostname
FOR /F "tokens=* USEBACKQ" %%F IN (`hostname`) DO (
    SET CHAMP_HOSTNAME=%%F
)

:: CPU
for /f "delims== tokens=2-3" %%f in ('wmic cpu get Name /value ^| find "="') do (
    set CHAMP_CPU=%%f
)

:: Architecture
for /f "delims== tokens=2-3" %%f in ('wmic os get OSArchitecture /value ^| find "="') do (
    set CHAMP_OS_ARCH=%%f
)

:: Add the CHAMP spec file to the Magma startup spec variable
set MAGMA_USER_SPEC=%CHAMP_DIR%\Source\CHAMP\CHAMP.s.m;%MAGMA_USER_SPEC%

::Set the CHAMP startup file as Magma startup file
set MAGMA_STARTUP_FILE=%CHAMP_DIR%\Source\CHAMP\CHAMP.m

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
::now, start magma with the Startup script from the Config directory
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%MAGMA_EXEC% -b %*
