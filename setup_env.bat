@echo off

IF NOT EXIST .venv (
    echo Creating virtual environment...
    python -m venv .venv
)

echo Activating environment...
call .venv\Scripts\activate

echo Installing requirements...
pip install --upgrade pip
pip install uv
uv pip install -r requirements.txt

echo.
echo Environment ready.
