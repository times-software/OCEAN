calc
{
xas
} calc
dft
{
abinit
} dft
nkpt
{
old
} nkpt
ngkpt
{
old
} ngkpt
ngkpt.auto
{
F
} ngkpt.auto
photon_q
{
0.00 0.00 0.00
} qinunitsofbvectors.ipt
dft.split
{
.true.
} dft.split
dft.qe_redirect
{
.false.
} dft.qe_redirect
nbands
{
range
} nbands
dft_energy_range
{
25
} dft_energy_range.ipt
obf_energy_range
{
-1
} obf_energy_range
obkpt
{
-1
} obkpt.ipt
obf.nbands
{
-1
} obf.nbands
trace_tol
{
5.0d-14
} trace_tol
acc_level
{
1
} acc_level.ipt
k0
{
0.125 0.250 0.375
} k0.ipt
fband
{
0.125
} fband
occopt
{
1
} occopt
mixing
{
0.7
} mixing
acell
{
} rscale
rprim
{
} rprim
ntypat
{
-1
} ntype
typat
{
} typat
znucl
{
} znucl
zsymb
{
 
} zsymb
pp_list
{
NULL
} pplist
pp_database
{
ONCVPSP-PBE-PDv0.4-stringent
} ppdatabase
ecut.quality
{
normal
} ecut.qualtiy
natom
{
-1
} natoms
coord
{
xred
} coord
xred
{
} taulist
ecut
{
} ecut
diemac
{
} epsilon
toldfe
{
1.0d-10
} etol
tolwfr
{
1.0d-12
} wftol
nstep
{
50
} nrun
dft.startingwfc
{
atomic+random
} dft.startingwfc
dft.diagonalization
{
david
} dft.diagonalization
dft.ndiag
{
-1
} dft.ndiag
dft.functional
{
none
} dft.functional
dft.exx.qmesh
{
1 1 1 
} dft.exx.qmesh
dft.nscf.poolsize
{
-1
} dft.nscf.poolsize
verbatim
{
#
} verbatim
para_prefix
{
# 
} para_prefix
ser_prefix
{
!
} ser_prefix
abpad
{
4
} abpad
scfac
{
1.00
} scfac
screen.shells
{
3.5
} screen.shells
opf.hfkgrid
{
2000  32
} opf.hfkgrid
opf.fill
{
!
} opf.fill
opf.opts
{
!
} opf.opts
opf.program
{
shirley
} opf.program
screen.nkpt
{
old
} screen.nkpt
screen.k0
{
0.125 0.250 0.375
} screen.k0
screen.nbands
{
0
} screen.nbands
caution
{
.false.
} caution
nedges
{
-1
} nedges.ipt
edges
{
#
} edges.ipt
cnbse.nbuse
{
0
} nbuse.ipt
cnbse.xmesh
{
-1
} xmesh.ipt
cnbse.rad
{
3.5
} cnbse.rad
metal
{
.false.
} metal
cksshift
{
0.00
} cksshift
cksstretch
{
1.00
} cksstretch
cnbse.niter
{
100
} cnbse.niter
haydock_convergence
{
-0.001 5 
} haydockconv.ipt
cnbse.spect_range
{
1200 -40.817 68.028
} cnbse.spect_range
cnbse.broaden
{
-1
} cnbse.broaden
cnbse.strength
{
1.0
} cnbse.strength
cnbse.solver
{
haydock
} cnbse.solver
cnbse.gmres.elist
{
false
} cnbse.gmres.elist
cnbse.gmres.erange
{
false
} cnbse.gmres.erange
cnbse.gmres.nloop
{
80
} cnbse.gmres.nloop
cnbse.gmres.gprc
{
13.605
} cnbse.gmres.gprc
cnbse.gmres.ffff
{
0.00000005
} cnbse.gmres.ffff
cnbse.write_rhs
{
.false.
} cnbse.write_rhs
cnbse.gw.control
{
none
} gw_control
bse.gw.cstr
{
0.0
} gwcstr
bse.gw.vstr
{
0.0
} gwvstr
bse.gw.gap
{
0.0 .false.
} gwgap
degauss
{
0.02
} degauss
ibrav
{
0
} ibrav
isolated
{
'mp'
} isolated
noncolin
{
.false.
} noncolin
prefix
{
system
} prefix
ppdir
{
../
} ppdir
dft.calc_stress
{
.false.
} dft.calc_stress
dft.calc_force
{
.true.
} dft.calc_force
spinorb
{
.false.
} spinorb
work_dir
{
'./Out'
} work_dir
tmp_dir
{
undefined
} tmp_dir
den.kshift
{
1 1 1
} den.kshift
core_offset
{
.false.
} core_offset
ham_kpoints
{
4 4 4 
} ham_kpoints
nbse.niter
{
100
} niter
nbse.backf
{
0
} backf
nbse.aldaf
{
0
} aldaf
nbse.qpflg
{
1
} qpflg
nbse.bwflg
{
0
} bwflg
nbse.bande
{
1
} bande
nbse.bflag
{
1
} bflag
nbse.lflag
{
1
} lflag
nbse.convergence
{
0
} convergence
nbse.decut
{
10 2
} decut
nbse.se_rs
{
-1
} se_rs
nbse.se_metal
{
0
} se_metal
nbse.se_niter
{
1
} se_niter
nbse.spect_range
{
1000 0 100 
} spect.h
tot_charge
{
0.0
} tot_charge
nspin
{
1
} nspin
smag
{

} smag
ldau
{

}ldau
qe_scissor
{

}qe_scissor
nphoton
{
-1
} nphoton
ser_bse
{
0
} serbse
spin_orbit
{
-1
} spin_orbit
screen_energy_range
{
100
} screen_energy_range.ipt
screen.grid.scheme
{
central
} screen.grid.scheme
screen.grid.rmode
{
uniform
} screen.grid.rmode
screen.grid.ninter
{
2
} screen.grid.ninter
screen.grid.shells
{
-1
8
} screen.grid.shells
screen.grid.xyz
{
2 2 2 
} screen.grid.xyz
screen.grid.rmax
{
8
} screen.grid.rmax
screen.grid.nr
{
25
} screen.grid.nr
screen.grid.ang
{
5
} screen.grid.ang
screen.grid.deltar
{
0.1
0.2
} screen.grid.deltar
screen.lmax
{
0
} screen.lmax
screen.grid.nb
{
24
} screen.grid.nb
screen.final.rmax
{
10
} screen.final.rmax
screen.final.dr
{
0.1
} screen.final.dr
screen.model.dq
{
0.01
} screen.model.dq
screen.model.qmax
{
10
} screen.model.qmax
screen.legacy
{
0
} screen.legacy
screen.augment
{
.true.
} screen.augment
screen.wvfn
{
new
} screen.wvfn
screen.convertstyle
{
intp
} screen.convertstyle
screen.inversionstyle
{
direct
} screen.inversionstyle
screen.mode
{
core
} screen.mode
bse.wvfn
{
new
} bse.wvfn
hamnum
{
1 1 1 1 1 1
} hamnum
echamp
{
.true.
} echamp.inp
bshift
{
0
} bshift
