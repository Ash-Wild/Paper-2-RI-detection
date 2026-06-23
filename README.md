# Paper-2-RI-detection

Repository of Python scripts used for cyclone rapid intensification (RI) detection, GNSS/CYGNSS data handling, and related meteorological analyses. This collection contains data downloaders, event-finders, analysis pipelines and visualisation helpers used in research workflows.

## Status
Research/analysis scripts — intended for interactive use and reproducible research, not a packaged application.

## Quick Start
1. Create and activate a Python virtual environment (PowerShell):

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

2. (Optional) Run helper batch scripts on Windows:

```powershell
./setup_env.bat   # configure environment shortcuts
./pipadd.bat     # helper to install additional packages
```

## Usage Examples
- Run the main RI analysis pipeline:

```powershell
python RI_analysis.py
```

- Find RI events:

```powershell
python RI_finder.py --input Merged_Measurements_RCG35.csv --out ri_events.csv
```

- Run a case study finder:

```powershell
python RI_casestudy_insitu_finder.py --case 2020-01
```

- Download CYGNSS L1 data (example):

```powershell
python cyg_L1_downloader.py --start 2023-01-01 --end 2023-01-07
```

Adjust script flags by opening the script headers or running `--help` when implemented.

## Files Overview
- `RI_analysis.py`: Main analysis pipeline and plotting routines for RI studies.
- `RI_finder.py`: Event detection utilities to generate RI candidate lists.
- `RI_casestudy.py` / `RI_casestudy_insitu_finder.py`: Case-study helpers and in-situ matching tools.
- `cyg_L1_downloader.py`: CYGNSS L1 data downloader and archive helpers.
- `CMS_CDS_api.py`: API wrappers for dataset access.
- `RI_plotter.py`, `Storm_timeline_visualiser.py`: Visualization and plotting utilities.
- `TCA_analysis.py`, `NOAA_Measurements_RCG35.csv`, `Merged_Measurements_RCG35.csv`: analysis inputs and example data.

See the `Old/` folder for legacy scripts and experiments; some are exploratory and may require modifications.

## Data
- Small example CSVs are included (`Merged_Measurements_RCG35.csv`, `NOAA_Measurements_RCG35.csv`).
- For larger/raw datasets, place files in a data directory and update script paths or pass an `--input` argument where supported.

## Dependencies
- See [requirements.txt](requirements.txt) for exact package versions.
- Recommended Python: 3.10+ (adjust if your environment differs).

## Configuration
- API keys and credentials: store required keys (NOAA, CYGNSS, other providers) as environment variables or in a local config file that is excluded from version control.
- Common env var examples (replace with actual variable names used by scripts):

```powershell
setx NOAA_API_KEY "your_key_here"
```

## Examples Workflow

```powershell
.venv\Scripts\Activate.ps1
python cyg_L1_downloader.py --start 2023-01-01 --end 2023-01-07
python RI_finder.py --input Merged_Measurements_RCG35.csv --out ri_events.csv
python RI_analysis.py --events ri_events.csv
```

## Contributing
- Open issues for bugs or feature requests.
- Submit changes via pull requests; keep commits focused and add tests where practical.

## License & Citation
- Add a `LICENSE` file to indicate the project's license.
- If used for academic work, include a short citation description here (author, year, DOI or arXiv).

## Contact
- Maintainer: Ashley Siedlecki (Wild), RMIT University Australia, ashley.siedlecki@rmit.edu.au

