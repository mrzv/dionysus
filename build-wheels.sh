#!/bin/sh

PROJECT_VERSION=$(uv run python -c 'import tomllib; print(tomllib.load(open("pyproject.toml", "rb"))["project"]["version"])') || exit 1

rm -r .venv
for PYTHON_VERSION in 3.10 3.11 3.12 3.13 3.14; do
  uv run --python $PYTHON_VERSION --with build python -m build || exit 1
done

rm -r .venv
for PYTHON_VERSION in 3.10 3.11 3.12 3.13 3.14; do
  orb -m ubuntu-plucky uv run --python $PYTHON_VERSION --with build python -m build || exit 1
done

orb -m ubuntu-plucky sh -c 'for fn in "dist/dionysus-$1"-*-linux*.whl; do uvx auditwheel repair --wheel-dir dist/wheelhouse --plat manylinux_2_39_x86_64 "$fn"; done' sh "$PROJECT_VERSION" || exit 1
