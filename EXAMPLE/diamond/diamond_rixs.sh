#!/bin/sh

OCEAN_BIN=/home/jtv1/bin/intel/rixs_test2/

## Run for XAS
echo "Running bse for spectra"
$OCEAN_BIN/ocean.pl diamond.in

mv CNBSE CNBSE_spect
##  


mkdir CNBSE

## Run XES
echo xes > Common/cnbse.mode
cd CNBSE
$OCEAN_BIN/cnbse_mpi.pl > xes.log
cd ..
##


## Calculate excitons
echo xas > Common/cnbse.mode
echo 1 > Common/nphoton
# switch to gmres
echo gmres > Common/cnbse.solver
# set the energy points to calculate
echo 0.58 > Common/cnbse.gmres.elist
echo 15 >> Common/cnbse.gmres.elist
echo .true. > CNBSE/echamp.inp

echo "Running bse for excitons"
cd CNBSE
$OCEAN_BIN/cnbse_mpi.pl > cnbse.log
cd ..
##


echo "Running RIXS w/ ai2nbse"
mkdir RIXS
cd RIXS
$OCEAN_BIN/RIXS.pl > rixs.log
cd ..

echo "Finished"
