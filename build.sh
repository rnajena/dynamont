#!/usr/bin/env bash
set -euo pipefail

BUILD_DIR="build"
PYTHON="python3"

if [[ "${1:-}" == "--clean" ]]; then
    echo "[INFO] Cleaning build..."
    rm -rf "${BUILD_DIR}"
fi

PYBIND11_DIR=$("${PYTHON}" -c 'import pybind11; print(pybind11.get_cmake_dir())')

echo "[INFO] Getting project version..."
PROJECT_VERSION=$("${PYTHON}" -c 'from scikit_build_core.metadata.setuptools_scm import dynamic_metadata; print(dynamic_metadata("version"))')

echo " -> Version: ${PROJECT_VERSION}"

JOBS=$(nproc)

echo "[INFO] Configuring CMake..."

cmake -S . \
    -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${BUILD_DIR}" \
    -DBUILD_PYTHON_BINDING=ON \
    -Dpybind11_DIR="${PYBIND11_DIR}" \
    -DPROJECT_VERSION="${PROJECT_VERSION}" \
    -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc"


echo "[INFO] Building..."
cmake --build "${BUILD_DIR}" --parallel "${JOBS}"
echo "[INFO] Installing CMake targets..."
cmake --install "${BUILD_DIR}"

echo "[INFO] Building Python package..."
"${PYTHON}" -m build
echo "[INFO] Build completed successfully."