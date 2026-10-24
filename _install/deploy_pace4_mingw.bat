@echo off

set "APPDIR=c:\Users\taras\Documents\Codex\PACE4\_install\"
set "WDEP=C:\Qt\6.7.0\mingw_64\bin\windeployqt.exe"

if not exist "%APPDIR%\PACE4.exe" (
  echo ERROR: Cannot find "%APPDIR%\PACE4.exe"
  echo Check APPDIR path.
  pause
  exit /b 1
)

if not exist "%WDEP%" (
  echo ERROR: Cannot find "%WDEP%"
  echo Check Qt path / version / kit folder.
  pause
  exit /b 1
)

cd /d "%APPDIR%" || (
  echo ERROR: cd failed to "%APPDIR%"
  pause
  exit /b 1
)

"%WDEP%" --release --compiler-runtime --dir . "PACE4.exe"

if errorlevel 1 (
  echo.
  echo ERROR: windeployqt failed.
  pause
  exit /b 1
)

echo.
echo DONE: PACE4 MSVC deployment completed.
pause
