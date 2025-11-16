# Miedema Enthalpy Calculation Program

Program for calculating enthalpy of formation using Miedema's model.

## Project Structure

```
miedema/
├── src/miedema/          # Source code
│   ├── core/             # Core calculation modules
│   ├── ui/               # User interface modules
│   └── utils/            # Utility functions
├── assets/               # UI assets (KV files)
├── data/                 # Data files (coefficients)
├── output/               # Output files
├── run.py                # Application entry point
└── pyproject.toml        # Project configuration
```

## Installation

### Quick setup

```bash
./setup.sh
```

### Manual setup

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

## Run

```bash
source .venv/bin/activate
python run.py
```

## Development

All dependencies are managed via `pyproject.toml` using `uv`. The project uses Python 3.10+.

## Usage

The application provides a graphical interface for:

- Defining crystal structures
- Calculating formation enthalpies
- Visualizing results
