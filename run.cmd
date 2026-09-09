@echo off
setlocal
"%~dp0.venv\Scripts\python.exe" "%~dp0check_environment.py"
if errorlevel 1 exit /b 1
"%~dp0.venv\Scripts\python.exe" "%~dp0cdr_solver.py" %*
if errorlevel 1 exit /b 1
"%~dp0.venv\Scripts\python.exe" "%~dp0verify_integrators.py"
if errorlevel 1 exit /b 1
"%~dp0.venv\Scripts\python.exe" "%~dp0verify_formulation.py"
exit /b %errorlevel%
