export OCEAN_BIN=/path/to/ocean/

$OCEAN_BIN/ocean.pl Fe.in

# Switch from using Haydock to GMRES
echo gmres > Common/cnbse.solver

# Run for a set of energies 
echo 0 > Common/cnbse.gmres.elist
echo 4.2 >> Common/cnbse.gmres.elist
echo 7.4 >> Common/cnbse.gmres.elist
echo 17.2 >> Common/cnbse.gmres.elist
echo 42.7 >> Common/cnbse.gmres.elist
echo 45.4 >> Common/cnbse.gmres.elist

# fix emission spectra range
echo 1500 0 150 > Common/spect.h

mkdir RIXS
# Which photon file(s) is/are incoming (can be a list)
echo 1 > RIXS/photon_in
# Which photon file(s) is/are outgoing (can be a list)
echo 2 > RIXS/photon_out

cd RIXS
$OCEAN_BIN/rixs.pl
