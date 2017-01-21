

ocean.pl diamond.in
echo gmres > Common/cnbse.solver
echo 0 30 2 > Common/cnbse.gmres.erange

mkdir RIXS
cd RIXS
rixs.pl
