#!/bin/bash
set -euo pipefail
shopt -s nullglob

DIONYSUS_VERSION=$(uv run --python 3.11 python -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])")

WHEELS=(dist/dionysus-$DIONYSUS_VERSION-*-linux*.whl)
if (( ${#WHEELS[@]} == 0 )); then
  echo "No Linux wheels found for dionysus $DIONYSUS_VERSION" >&2
  exit 1
fi

REPAIRED_DIR=$(mktemp -d "${TMPDIR:-/tmp}/dionysus-auditwheel.XXXXXX")
trap 'rm -rf "$REPAIRED_DIR"' EXIT
UNREPAIRED_DIR=dist/unrepaired-linux
mkdir -p "$UNREPAIRED_DIR"

RAW_WHEELS=()
for fn in "${WHEELS[@]}"; do
  raw_fn="$UNREPAIRED_DIR/$(basename "$fn")"
  mv "$fn" "$raw_fn"
  RAW_WHEELS+=("$raw_fn")
done

for fn in "${RAW_WHEELS[@]}"; do
  uvx auditwheel repair --plat manylinux_2_39_x86_64 --wheel-dir "$REPAIRED_DIR" "$fn"
done

REPAIRED_WHEELS=("$REPAIRED_DIR"/*.whl)
if (( ${#REPAIRED_WHEELS[@]} == 0 )); then
  echo "auditwheel did not produce any repaired wheels" >&2
  exit 1
fi

for fn in "${REPAIRED_WHEELS[@]}"; do
  mv "$fn" dist/
done
