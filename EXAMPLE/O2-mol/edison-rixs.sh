#!/usr/bin/env bash
#SBATCH --job-name=O2         
#SBATCH --partition=debug            
#SBATCH --time=00:30:00
#SBATCH -N 2   
##SBATCH -C haswell            


export OCEAN_BIN=/global/homes/p/pemmaraj/Builds/Edison/OCEAN/BIN-213               

$OCEAN_BIN/ocean.pl O2.in > out.VAL

## Switch from using Haydock to GMRES
#echo gmres > Common/cnbse.solver
#
## Run for every energy in a range (same energy scale as the absorption spectra)
## start stop step, ie, run for 0, 2, 4, ... 30
#echo -5 25 0.5 > Common/cnbse.gmres.erange
#
## Which photon file(s) is/are incoming (can be a list)
#echo 1 > Common/photon_in
## Which photon file(s) is/are outgoing (can be a list)
#echo 2 > Common/photon_out
#
#mkdir RIXS
#cd RIXS
#$OCEAN_BIN/rixs.pl
