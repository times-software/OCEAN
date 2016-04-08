## v. 1.1.5

#### Bug Fix
Changed mpi_double to mpi_double_precision in several files. 
Thanks to Liang Li at ANL for the bug report.

## v. 1.1.4

#### Bug Fix
The array ordering of calls to FFTW was switched (C/Fort). Only compilations using 
FFTW (-D__FFTW3 in Makefile.arch) would have been affected.

## v. 1.1.3

#### Minor feature updates:
 1. Option to disable spin-orbit splitting from the input file
 2. Code can generate its own cnbse.spect_range
 3. Code detects when number of xpoints is too small compared to bands


#### Under the hood:
 1. Some checks for NaN cropping up in BSE
 2. Some support for FFTW instead of legacy FFT in OBF pathway
 3. Some support for RIXS workflows
