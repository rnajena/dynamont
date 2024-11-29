# ![Alt text](figures/logo.png)

- [**Dynamont**](https://github.com/JannesSP/dynamont) is a segmentation/resquiggling tool for ONT signals. ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/dynamont)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-teal.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![PyPI](https://img.shields.io/pypi/v/dynamont) ![PyPI - Downloads](https://img.shields.io/pypi/dm/dynamont) ![PyPI - Status](https://img.shields.io/pypi/status/dynamont)](https://pypi.org/project/dynamont/)

[![Anaconda-Server Badge](https://anaconda.org/jannessp/dynamont/badges/version.svg)](https://anaconda.org/jannessp/dynamont) ![Conda](https://img.shields.io/conda/dn/jannessp/dynamont) [![Conda package](https://anaconda.org/jannessp/dynamont/badges/latest_release_date.svg)](https://anaconda.org/jannessp/dynamont) [![Conda package](https://anaconda.org/jannessp/dynamont/badges/platforms.svg)](https://anaconda.org/jannessp/dynamont)

<!-- [![DOI](https://zenodo.org/badge/633012569.svg)](https://zenodo.org/badge/latestdoi/633012569) -->

---

- [](#)
  - [Dynamont](#dynamont)
  - [Installation](#installation)
    - [Pypi/pip](#pypipip)
    - [Conda](#conda)
  - [Usage](#usage)
  - [Exit-Codes](#exit-codes)

---

## Dynamont

TODO: change title
A **dynam**ic programming algorithm to resquiggle and segment the raw **ONT** signal.

TODO: add readme
TODO: create conda package of tool
TODO: build tool with subtools

## Installation

### Pypi/pip

```bash
pip install dynamont
```

### Conda

```bash
conda install mamba
mamba create -n dynamont -c jannessp dynamont
conda activate dynamont
```

## Usage

TODO

## Exit-Codes

- -9: Out of Memory error. Deacrease the number of processes or move to a system with more memory.
- 1: `dynamont-NTK` specific: alignment score (Z) does not match between forward and backward run in preprocessing on signal (T) and read (N) in function: `preProcTN`
- 2: `dynamont-NTK` specific: alignment score (Z) does not match between forward and backward run in preprocessing on signal (T) and de novo calling (K) in function: `preProcTK`
- 3: Alignment score (Z) does not match between forward and backward run in `main` function or is -Infinity
- 4: Input signal is missing or not found in stdin stream
- 5: Input read is missing or not found in stdin stream
- 6: Model file contains a k-mer that does not match the k-mer size of the pore
- 7: Invalid model path was provided
- 8: Provided ONT signal too short
- 9: Provided reads too short
- 10: Signal is smaller than read
- 11: Read is smaller than `kmerSize` of pore type
