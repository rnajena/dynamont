#!/bin/bash
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
make install
cd ..
python -m pip install --user --upgrade setuptools wheel build --no-deps --ignore-installed -v .
python -m build
# twine upload dist/*