#!/bin/bash

set -e

echo "Creating virtual environment..."
uv venv

echo "Activating virtual environment..."
source .venv/bin/activate

echo "Installing dependencies..."
uv pip install -e .

echo "Setup complete! Run the application with:"
echo "  source .venv/bin/activate"
echo "  python run.py"

