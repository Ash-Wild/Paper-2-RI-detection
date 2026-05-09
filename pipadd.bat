@echo off

IF NOT EXIST .venv\Scripts\python.exe (
    echo venv not found. Create it first with:
    echo python -m venv .venv
    exit /b 1
)

IF EXIST .venv\Scripts\uv.exe (
    .venv\Scripts\uv.exe pip install %*
) ELSE (
    .venv\Scripts\python -m pip install %*
)

IF ERRORLEVEL 1 exit /b 1

.venv\Scripts\python -m pip freeze > requirements.txt

echo.
echo requirements.txt updated.
