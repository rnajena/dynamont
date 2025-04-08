#!/bin/bash

BUILD_DIR="build"
BIN_DIR="src/bin/"

# Ensure clean build dirs
rm -rf "$BUILD_DIR" "$BIN_DIR"
mkdir -p "$BUILD_DIR" "$BIN_DIR"
cd $BUILD_DIR

# Extract version from VCS (via hatch)
echo "[INFO] Getting project version..."
PROJECT_VERSION=$(hatch version)
export PROJECT_VERSION
echo " -> Version: $PROJECT_VERSION"

# Build C++ binaries with static libstdc++/libgcc
echo "[INFO] Building C++ binaries with CMake..."
cmake .. \
  -DPROJECT_VERSION="$PROJECT_VERSION" \
  -DCMAKE_INSTALL_PREFIX=. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXE_LINKER_FLAGS="-static-libstdc++ -static-libgcc" \
  -DCMAKE_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"
make -j$(nproc)
make install
cd ..

mkdir -p src/bin
cp $BUILD_DIR/bin/dynamont-NTC $BIN_DIR
cp $BUILD_DIR/bin/dynamont-NT-banded $BIN_DIR
strip "$BIN_DIR/dynamont-NTC" "$BIN_DIR/dynamont-NT-banded"

# Build Python package (wheel + sdist)
echo "[INFO] Building Python package..."
# python -m pip install --user --upgrade setuptools wheel build --ignore-installed -v .
python -m build
# twine upload dist/*
