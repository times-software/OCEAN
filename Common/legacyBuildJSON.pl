#!/usr/bin/env perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
use strict;
require JSON::PP;
JSON::PP->import;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/legacyBuildJSON\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }


my $json = JSON::PP->new;
my $input;

if( open( my $in, "oparse.json" ))
{
  local $/ = undef;
  $input = $json->decode(<$in>);
  close($in);
} 
else
{
  die "Failed to open oparse.json which should have been copied in from the OCEAN install directory!\m$!";
}

my %decoder = (
  'calc' => 'calc.mode',
  'dft' => 'dft.program',
  'nkpt' => 'bse.kmesh',
  'ngkpt' => 'dft.den.kmesh',
  'ngkpt.auto' => 'dft.den.auto',
  'qinunitsofbvectors.ipt' => 'calc.photon_q',
  'dft.split' => 'dft.bse.split',
  'dft.qe_redirect' => 'dft.redirect',
  'nbands' => 'bse.nbands',
  'dft_energy_range.ipt' => 'bse.dft_energy_range',
  'obf_energy_range' => 'nope.obf_energy_range',
  'obkpt.ipt' => 'nope.obkpt',
  'obf.nbands' => 'nope.obf_nbands',
  'trace_tol' => 'nope.trace_tol',
  'acc_level.ipt' => 'nope.acc_level',
  'k0.ipt' => 'bse.kshift',
  'fband' => 'dft.fband',
  'occopt' => 'dft.occopt',
  'mixing' => 'dft.mixing',
  'rscale' => 'structure.rscale',
  'rprim' => 'structure.rprim',
  'ntype' => 'nope.ntype',
  'typat' => 'structure.typat',
  'znucl' => 'structure.znucl',
  'zsymb' => 'structure.zsymb',
  'pplist' => 'psp.pp_list',
  'ppdatabase' => 'psp.pp_database',
  'ecut.qualtiy' => 'psp.ecut_quality',
  'natoms' => 'nope.natoms',
  'coord' => 'structure.coord',
  'taulist' => 'structure.taulist',
  'ecut' => 'dft.ecut', 
  'epsilon' => 'structure.epsilon',
  'etol' => 'dft.toldfe',
  'wftol' => 'dft.tolwfr',
  'nrun' => 'dft.nstep',
  'dft.startingwfc' => 'dft.startingwfc',
  'dft.diagonalization' => 'dft.diagonalization',
  'dft.ndiag' => 'dft.ndiag',
  'dft.functional' => 'dft.functional',
  'dft.exx.qmesh' => 'dft.exx.qmesh',
  'dft.nscf.poolsize' => 'dft.bse.poolsize',
  'verbatim' => 'dft.verbatim',
  'para_prefix' => 'computer.para_prefix',
  'ser_prefix' => 'computer.ser_prefix',
  'abpad' => 'dft.abpad',
  'scfac' => 'bse.core.scfac',
  'screen.shells' => 'screen.shells',
  'opf.hfkgrid' => 'opf.shirley.hfkgrid',
  'opf.fill' => 'opf.shirley.fill',
  'opf.opts' => 'opf.shirley.opts',
  'opf.program' => 'opf.program',
  'screen.nkpt' => 'screen.kmesh',
  'screen.k0' => 'screen.kshift',
  'screen.nbands' => 'screen.nbands',
  'caution' => 'opf.shirley.caution',
  'nedges.ipt' => 'nope.nedges',
  'edges.ipt' => 'calc.edges',
  'nbuse.ipt' => 'nope.nbuse',
  'xmesh.ipt' => 'bse.xmesh',
  'cnbse.rad' => 'bse.core.screen_radius',
  'metal' => 'structure.metal',
  'cksshift' => 'nope.cksshift',
  'cksstretch' => 'nope.cksstretch',
  'cnbse.niter' => 'bse.core.haydock.niter',
  'haydockconv.ipt' => 'bse.core.haydock.converge.thresh',
  'cnbse.spect_range' => 'bse.core.plot.range',
  'cnbse.broaden' => 'bse.core.broaden',
  'cnbse.strength' => 'bse.core.strength',
  'cnbse.solver' => 'bse.core.solver',
  'cnbse.gmres.elist' => 'bse.core.gmres.elist',
  'cnbse.gmres.erange' => 'bse.core.gmres.erange',
  'cnbse.gmres.nloop' => 'bse.core.gmres.nloop',
  'cnbse.gmres.gprc' => 'bse.core.gmres.gprc',
  'cnbse.gmres.ffff' => 'bse.core.gmres.ffff',
  'cnbse.write_rhs' => 'bse.core.write_rhs',
  'gw_control' => 'bse.core.gw.control',
  'gwcstr' => 'nope.bse_gw_cstr',
  'gwvstr' => 'nope.bse_gw_vstr',
  'gwgap' => 'nope.gwgap',
  'degauss' => 'dft.degauss',
  'ibrav' => 'nope.ibrav',
  'isolated' => 'nope.isolated',
  'noncolin' => 'dft.isolated',
  'prefix' => 'nope.prefix',
  'ppdir' => 'psp.ppdir',
  'dft.calc_stress' => 'dft.calc_stress',
  'dft.calc_force' => 'dft.calc_force',
  'spinorb' => 'dft.spinorb',
  'work_dir' => 'nope.wordir',
  'tmp_dir' => 'dft.tmp_dir',
  'den.kshift' => 'dft.den.kshift',
  'core_offset' => 'screen.core_offset',
  'ham_kpoints' => 'nope.ham_kpoints',
  'niter' => 'bse.val.haydock.niter',
  'backf' => 'bse.val.backf',
  'aldaf' => 'bse.val.aldaf',
  'qpflg' => 'bse.val.qpflg',
  'bwflg' => 'bse.val.bwflg',
  'bande' => 'bse.val.bande',
  'bflag' => 'bse.val.bflag',
  'lflag' => 'bse.val.lflag',
  'convergence' => 'nope.convergence',
  'decut' => 'bse.val.decut',
  'se_rs' => 'nope.se.rs',
  'se_metal' => 'nope.se.metal',
  'se_niter' => 'nope.se.niter',
  'spect.h' => 'bse.val.plot.range',
  'tot_charge' => 'dft.tot_charge',
  'nspin' => 'dft.nspin',
  'smag' => 'dft.smag',
  'ldau' => 'dft.ldau.Hubbard_U',
  'qe_scissor' => 'nope.qe_scissor',
  'nphoton' => 'nope.nphoton',
  'serbse' => 'nope.ser_bse',
  'spin_orbit' => 'bse.core.spin_orbit',
  'screen_energy_range.ipt' => 'screen.dft_energy_range',
  'screen.grid.scheme' => 'screen.grid.scheme',
  'screen.grid.rmode' => 'screen.grid.rmode',
  'screen.grid.ninter' => 'nope.screen_grid_ninter',
  'screen.grid.shells' => 'screen.grid.shells',
  'screen.grid.xyz' => 'nope.screen_grid_xyz',
  'screen.grid.rmax' => 'screen.grid.rmax', 
  'screen.grid.nr' => 'nope.screen_grid_nr', 
  'screen.grid.ang' => 'screen.grid.ang',
  'screen.grid.deltar' => 'screen.grid.deltar', 
  'screen.lmax' => 'screen.grid.lmax', 
  'screen.grid.nb' => 'nope.screen_grid_nb',
  'screen.final.rmax' => 'screen.final.rmax', 
  'screen.final.dr' => 'screen.final.dr', 
  'screen.model.dq' => 'screen.model.SLL.dq',
  'screen.model.qmax' => 'screen.model.SLL.qmax',
  'screen.legacy' => 'nope.screen_legacy',
  'screen.augment' => 'screen.augment',
  'screen.wvfn' => 'nope.screen_wvfn',
  'screen.convertstyle' => 'screen.convertstyle',
  'screen.inversionstyle' => 'screen.inversionstyle',
  'screen.mode' => 'screen.mode', 
  'bse.wvfn' => 'nope.bse_wvfn',
  'hamnum' => 'nope.hamnum',
  'echamp.inp' => 'bse.core.gmres.echamp',
  'bshift' => 'nope.bshift' );
  

# TEST MODE -- make sure all the input flags are correctly here
if( 0 )
{
  open IN, "oparse.h" or die;
  while( my $line = <IN> )
  { 
    if( $line =~ m/[\{\}]/ ) { die "Mix up in oparse.h: expected input\n$line";}
    $line = <IN>;
    die "$line. Mix up in oparse.h\n" unless( $line =~ m/\{/ );
    $line = <IN>;
    print "$line";
    $line = <IN> unless( $line =~ m/\}/ );
    $line =~ m/\}\s*(\S+)/ or die "$line";
    my $tag = $1;
    chomp($tag);
    die "Failed to find tag: $tag\n" unless( exists $decoder{ $tag } );
  }
}

my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);

open OUT, ">", "json.test" or die;
print OUT $json->encode($input);
close OUT;
