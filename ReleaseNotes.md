## v. 2.0.3

#### Bugfixes
 1. A call to OpenMP wasn't behind a sentinel
 2. Fixed a crashing bug with FFTs in the BSE

## v. 2.0.2

#### Minor Features
 1. Cleaned up the openmp directives in the BSE section

## v. 2.0.1

#### Bugfixes
 1. Fixed major bug in the valence BSE regarding FFTs. 

## v. 2.0

#### Major Features:
1. RIXS! The valence BSE code (formerly AI2NBSE by Hadley Lawler, et al) is now included, allowing for direct RIXS calculations with a relatively simple workflow. See the diamond example. 
2. The optimal basis functions (PAW-style reconstruction) has been improved. Thanks to Eric Shirley for these changes.

#### Minor Features:
 1. Changes under the hood to MPI calls in the BSE.
 2. More consistent determination of the Fermi level for metallic systems.
 3. Better switching between (optional) FFTW3 and legacy FFT

#### Bug fixes
 1. Some bugs that surfaced with MPI calls, see Makefile.arch.example for new -D__OLD_MPI flag.

## v. 1.2.1

#### Minor Feature: Averaging spectra
At the end of the CNBSE section the code will automatically create averages: 
 1. By site over all photon files, e.g. for multiple polarizations
 2. By photon file over all sites
 3. Everything for the given edge

## v. 1.2.0

#### Major Feature: LDA+U
Fixed up LDA+U support (QuantumEspresso only). Thanks to Yufeng Liang for 
spearheading this effort.

#### Minor Features:
 1. Improvements to the parsing of QE wavefunctions
 2. More consistent determination of the Fermi level for metallic systems. 

#### Bug fixes
 1. Some bugs that surfaced with MPICH
 2. Cleaned up the OCEAN2/Makefile

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
