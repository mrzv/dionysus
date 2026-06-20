#!/bin/bash
set -euo pipefail

uv lock --check

rm -rf .venv
for PYTHON_VERSION in 3.10 3.11 3.12 3.13 3.14; do
  uv run --python "$PYTHON_VERSION" --with build python -m build
done
./repair-macos-wheels.sh

rm -rf .venv
for PYTHON_VERSION in 3.10 3.11 3.12 3.13 3.14; do
  orb -m ubuntu-plucky uv run --python "$PYTHON_VERSION" --with build python -m build
done

orb -m ubuntu-plucky ./repair-wheels.sh
