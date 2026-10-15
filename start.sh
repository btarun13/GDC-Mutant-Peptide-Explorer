#!/usr/bin/env bash
# First run: create .venv, install packages, fetch MHCflurry models. Every run: start the app.
set -euo pipefail
cd "$(dirname "$0")"

PYTHON="${PYTHON:-python3.11}"

if [ ! -x .venv/bin/python ]; then
  if ! command -v "$PYTHON" >/dev/null 2>&1; then
    echo "Could not find '$PYTHON'. Install Python 3.11, or run: PYTHON=/path/to/python ./start.sh" >&2
    exit 1
  fi
  "$PYTHON" -c 'import sys; sys.exit(0 if (3, 10) <= sys.version_info[:2] <= (3, 12) else 1)' || {
    echo "Python 3.10-3.12 is required; '$PYTHON' is $("$PYTHON" --version 2>&1)." >&2
    exit 1
  }
  echo "Creating virtual environment in .venv ..."
  "$PYTHON" -m venv .venv
  .venv/bin/python -m pip install --upgrade pip
  echo "Installing packages (about 1 GB, mostly PyTorch) ..."
  if ! .venv/bin/python -m pip install -r requirements.txt; then
    rm -rf .venv  # do not leave a half-installed environment behind
    echo "Package installation failed; see the messages above." >&2
    exit 1
  fi
fi

# One-time model download. MHCflurry keeps models in its own user data folder, outside .venv.
if ! .venv/bin/python -c 'from mhcflurry import Class1PresentationPredictor as P; P.load()' \
    >/dev/null 2>&1; then
  echo "Downloading MHCflurry models (one-time) ..."
  .venv/bin/mhcflurry-downloads fetch
fi

echo "Starting the app at http://localhost:8501 (Ctrl+C to stop) ..."
# Headless skips Streamlit's first-run email prompt (it blocks the terminal, and crashes
# without one); open the browser ourselves instead.
if [ -z "${NO_BROWSER:-}" ]; then
  opener="$(command -v open || command -v xdg-open || true)"
  [ -n "$opener" ] && (sleep 4 && "$opener" http://localhost:8501 >/dev/null 2>&1) &
fi
exec .venv/bin/streamlit run app.py --server.headless true --browser.gatherUsageStats false
