#!/bin/bash

# Get the version from Hatch (Python)
VERSION=$(hatch version)  # Extract the short version string
export PROJECT_VERSION=$VERSION  # Set environment variable

mkdir -p _build
cd _build
# cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local -DCMAKE_BUILD_TYPE=Release
cmake .. -DPROJECT_VERSION="$PROJECT_VERSION" -DCMAKE_INSTALL_PREFIX=. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make install
cd ..

mkdir -p src/bin
cp _build/bin/dynamont-NTC src/bin/
cp _build/bin/dynamont-NT-banded src/bin/

# python -m pip install --user --upgrade setuptools wheel build --ignore-installed -v .
python -m build
# twine upload dist/*
