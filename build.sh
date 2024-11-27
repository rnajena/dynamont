#!/bin/bash
mkdir -p build
cd build
cmake ..
make
cd ..
python3 -m pip install --user --upgrade setuptools wheel build --no-deps --ignore-installed -v .
python3 -m build
twine upload dist/*