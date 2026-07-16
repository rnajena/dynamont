#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="build"
BIN_DIR="src/bin"

if [[ "${1:-}" == "--clean" ]]; then
    echo "[INFO] Cleaning build..."
    rm -rf "${BUILD_DIR}" "${BIN_DIR}"
fi

mkdir -p "${BIN_DIR}"

for cmd in cmake hatch python; do
    if ! command -v "${cmd}" >/dev/null; then
        echo "[ERROR] ${cmd} not found"
        exit 1
    fi
done

echo "[INFO] Getting project version..."
PROJECT_VERSION=$(hatch version)
echo " -> Version: ${PROJECT_VERSION}"
export PROJECT_VERSION

if command -v nproc >/dev/null; then
    JOBS=$(nproc)
elif command -v sysctl >/dev/null; then
    JOBS=$(sysctl -n hw.ncpu)
else
    JOBS=1
fi

echo "[INFO] Configuring CMake..."

cmake -S . \
      -B "${BUILD_DIR}" \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${BUILD_DIR}" \
      -DPROJECT_VERSION="${PROJECT_VERSION}" \
      -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc"


echo "[INFO] Building..."
cmake --build "${BUILD_DIR}" --parallel "${JOBS}"
echo "[INFO] Installing binaries..."
cmake --install "${BUILD_DIR}"
echo "[INFO] Copying binaries..."

cp "${BUILD_DIR}/bin/dynamont-NTC" \
   "${BUILD_DIR}/bin/dynamont-NT-banded" \
   "${BUILD_DIR}/bin/dynamont-NT" \
   "${BIN_DIR}/"

echo "[INFO] Stripping binaries..."

strip \
    "${BIN_DIR}/dynamont-NTC" \
    "${BIN_DIR}/dynamont-NT-banded" \
    "${BIN_DIR}/dynamont-NT"

echo "[INFO] Building Python package..."
python -m build
echo "[INFO] Build completed successfully."