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

Quick start (auto setup if needed):

```bash
./start.sh
```

Or manually:

```bash
source .venv/bin/activate
python run.py
```

## Development

All dependencies are managed via `pyproject.toml` using `uv`. The project uses Python 3.10+.

## Usage

The application provides a graphical interface for:

- Defining crystal structures
- Calculating formation enthalpies using Miedema model
- Visualizing results

### Cutoff Parameter

- Use `-1` for automatic cutoff calculation (recommended)
- Automatic: `cutoff = min(a, b, c) × 0.8`
- Manual: specify custom value (e.g., `2.0`, `3.5`)

## Known Issues

**Matplotlib visualization**: Garden matplotlib has compatibility issues with newer matplotlib versions. Interactive features (zoom, pan) are disabled but static visualization works. See `docs/FIXES.md` for details.

Other known issues and solutions are documented in `docs/FIXES.md`.
