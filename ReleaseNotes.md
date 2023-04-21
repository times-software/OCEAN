## v. 3.0.4

#### Minor features
 1. The valence BSE can now self-consistently determine the static dielectric 
    constant by setting bse.val.epsilon_threshold to be > 0. 
 2. When using QuantumESPRESSO, the static dielctric constant of insulating 
    systems can be determined using density-functional perturbation theory. 
    This is attempted automatically if there is no dielectric constant set in 
    the input file. Systems run as metals will use a default value. 

#### Bug fixes
 1. Improved compatibility with older versions of QE. 
 2. Problems in the valence BSE with spin systems for both valence and RIXS.

## v. 3.0.3

#### Bug fixes
 1. A small, odd number for the first k-point dimension combined with a 
    different number of x-points in the first two dimensions could lead to 
    problems. 
 2. Some x-mesh grids could lead to failures in processing the density for
    valence calculations

## v. 3.0.2

#### Minor features
 1. MPSE support added back. The many-pole self-energy model from AI2NBSE has 
    been re-enabled. The scripts and source are in POST/MPSE. This is a post-
    processing step for after valence UV/optical calculations. 

#### Bug fixes
 1. Screening had small bug. Mostly didn't seem to change spectra, but 
    induced potential could change slightly from run to run. 

## v. 3.0.1

#### Major Features
 1. Completely new formatting for input files (the old format or a mix of old 
    and new are supported). The new style maps to a json file. See the file 
    Common/postDefaultsOceanDatafile for a complete picture of the input 
    following parsing the user-supplied inputs and filling in with defaults.

#### Bug fixes
 1. The Umklapp wasn't functioning correctly for large-q valence calculations
    (thanks Ishiaka Mansaray for the bug report).
 2. Potential problem in parsing the 'nelec' element in the QE xml output 
    (thanks Max Radin for the bug report).

## v. 2.9.7

### This is the pre-release/test run for 3.0

#### Major Features
 1. Support for D. R. Hamann's Optimized Norm-conserving Vanderbilt pseudopotentials
    1. Create your own
    2. Use built in database from PseudoDoJo collection
 2. Improvements to parallelism
    1. Screening calculation is faster
    2. Prep stage is faster too
 3. The DFT calculations for the BSE are split in two whenever there is finite-q

#### Minor Features
 1. OPF calculation should be slightly more robust, but noticeable differences are unlikely
 2. Large finite-q valence calculations should be working again

## v. 2.5.2

#### Major Features
 1. Impoved screening calculations
    1. Augmentation of pseudo-wavefunctions to restore all-electron character
    2. Much faster
 2. Better core-level shift calculations

#### Bugfixes
 * Fixed bugs in new screening method in 2.5.0 and 2.5.1

## v. 2.1.1

#### Bugfixes
 1. Fixed bug in valence (UV/optical) calculations
 2. Failed to initialize a variable in mpi_avg.x
 3. Rare bug gives bad densities in mpi_avg.x

## v. 2.1.0

#### Major Features
 1. No longer need iotk.a for QE support. Faster PREP stage for QE
 2. Cleaned up inputs and defaults

#### Minor Features
 1. Improvements for spin with QE
 2. Database of core hole lifetime widths
 3. Changes to alignment! By default, OCEAN sets the DFT LUMO = 0. However! 
    When core_offset is specified the DFT energies will **not** be adjusted. 

#### Bugfixes
 1. Fixed bug in calculating the core-hole screening for some metals

## v. 2.0.4

#### Minor Features
 1. Improvements to OPF generation and error checking
 2. Consistency between ABINIT and QuantumESPRESSO occupation numbers (occopt/smearing/etc)
 3. Updates to core-level shifts to improve feedback to user. Updated documentation.

## v. 2.0.3

#### Bugfixes
 1. A call to OpenMP wasn't behind a sentinel
 2. Fixed a crashing bug with FFTs in the BSE
 3. Fixed DFT+U support for QE-6.0+

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
