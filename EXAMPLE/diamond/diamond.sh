export OCEAN_BIN=~/cluster/bin/ocean/develop

$OCEAN_BIN/ocean.pl diamond.in

# Switch from using Haydock to GMRES
echo gmres > Common/cnbse.solver

# Run for every energy in a range (same energy scale as the absorption spectra)
# start stop step, ie, run for 0, 2, 4, ... 30
echo 0 30 2 > Common/cnbse.gmres.erange

# Which photon file(s) is/are incoming (can be a list)
echo 1 > Common/photon_in
# Which photon file(s) is/are outgoing (can be a list)
echo 2 > Common/photon_out

mkdir RIXS
cd RIXS
$OCEAN_BIN/rixs.pl
