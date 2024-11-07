# Dynamont

TODO: change title
A **dynam**ic programming algorithm to resquiggle and segment the raw **ONT** signal.

TODO: add readme
TODO: create conda package of tool
TODO: build tool with subtools

## Exit Codes

- 1: dynamont_NTK specific: alignment score (Z) does not match between forward and backward run in preprocessing on signal (T) and read (N) in function: `preProcTN`
- 2: dynamont_NTK specific: alignment score (Z) does not match between forward and backward run in preprocessing on signal (T) and de novo calling (K) in function: `preProcTK`
- 3: alignment score (Z) does not match between forward and backward run in `main` function
- 4: Input signal is missing or not found in stdin stream
- 5: Input read is missing or not found in stdin stream
- 6: Model file contains a k-mer that does not match the k-mer size of the pore
