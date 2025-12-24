#!/bin/sh

rm -r .venv
for PYTHON_VERSION in 3.9 3.10 3.11 3.12 3.13 3.14; do
  uv run --python $PYTHON_VERSION --with build python -m build || exit 1
done

rm -r .venv
for PYTHON_VERSION in 3.9 3.10 3.11 3.12 3.13 3.14; do
  orb -m ubuntu-plucky uv run --python $PYTHON_VERSION --with build python -m build || exit 1
done
