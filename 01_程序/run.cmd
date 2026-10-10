@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo Python environment is not ready.
  echo Run setup.cmd first, then run this file again.
  exit /b 1
)
".venv\Scripts\python.exe" "check_environment.py"
if errorlevel 1 exit /b 1
".venv\Scripts\python.exe" "example_vector.py" --dim 2 --n 4 --p 1 --steps 4 --output "output\demo" %*
exit /b %errorlevel%
