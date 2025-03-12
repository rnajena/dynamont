# ![Dynamont](figures/logo.png)

A **Dynam**ic Programming Approach to Segment **ONT** Signals
Dynamont is a segmentation/resquiggling tool for ONT signals.

![PyPI - Python Version](https://img.shields.io/pypi/pyversions/dynamont)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI](https://img.shields.io/pypi/v/dynamont) ![PyPI - Downloads](https://img.shields.io/pypi/dm/dynamont)](https://pypi.org/project/dynamont/)
[![Anaconda-Server Badge](https://anaconda.org/jannessp/dynamont/badges/version.svg)](https://anaconda.org/jannessp/dynamont) ![Conda](https://img.shields.io/conda/dn/jannessp/dynamont) [![Conda package](https://anaconda.org/jannessp/dynamont/badges/latest_release_date.svg)](https://anaconda.org/jannessp/dynamont) [![Conda package](https://anaconda.org/jannessp/dynamont/badges/platforms.svg)](https://anaconda.org/jannessp/dynamont)

<!-- [![DOI](https://zenodo.org/badge/633012569.svg)](https://zenodo.org/badge/latestdoi/633012569) -->

---

- [Installation](#installation)
  - [Pypi/pip](#pypipip)
  - [Conda](#conda)
- [Usage](#usage)
- [Exit-Codes](#exit-codes)

---

## Installation

### Pypi/pip

```bash
pip install dynamont
```

### Conda

```bash
conda create -n dynamont jannessp::dynamont
conda activate dynamont
```

## Usage

```bash
# segment a dataset
dynamont-resquiggle -r <path/to/pod5/dataset/> -b <basecalls.bam> --mode basic --model_path <path/to/model> -o <output.csv> -p <pore>

# train model
dynamont-train -r <path/to/pod5/dataset/> -b <basecalls.bam> --mode basic --model_path <path/to/init/model> -o <output/path> -p <pore>
```

## Exit-Codes

- -11: Segmentation fault
- -9: Out of Memory error. Decrease the number of processes or move to a system with more memory.
- -6: std::bad_alloc
- 1: `resquiggle mode` specific: alignment score (Z) does not match between forward and backward run in preprocessing on signal (T) and read (N).
- 2: `resquiggle mode` specific: alignment score (Z) does not match between forward and backward run in preprocessing on signal (T) and error correction (C).
- 3: Alignment score (Z) does not match between forward and backward pass or is -Infinity
- 4: Input signal is missing or not found in stdin stream
- 5: Input read is missing or not found in stdin stream
- 6: 
- 7: Invalid model path was provided
- 8: Provided ONT signal is too short
- 9: Read is too short
- 10: Signal is smaller than read
- 11: Read is smaller than `kmerSize` of provided pore model
