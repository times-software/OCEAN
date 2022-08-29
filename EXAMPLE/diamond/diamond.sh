export OCEAN_BIN=~/cluster/bin/ocean/develop

$OCEAN_BIN/ocean.pl diamond.in
mv CNBSE CNBSE_XAS
$OCEAN_BIN/ocean.pl diamond.rixs.in
