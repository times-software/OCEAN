#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;


sub QErunDensity
{
  my $hashRef = $_[0];
  my $errorCode = 1;

  # open input file

  # Modify QEprintInput to take two different hashes and a string 
  # QEprintInput( $fh, $hashRef->{'general'}, $hashRef->{'scf'}, "scf" ) ;


  # test run for k-points

  # parse test results

  # full run

  return $errorCode;
}


# subroutine to handle running each QE step
# make one of these for pp and ph too!
sub QErunPW
{

}


# Subroutine to print a QE-style input file
#   Pass in a file handle and a hashtable
sub QEprintInput
{
  my ($fh, %inputs ) = @_;


  # Array of QE names for smearing by occopt
  my @QE_smear;
  $QE_smear[1] = "'gaussian'";     # Still need to fix to be insulator
  $QE_smear[3] = "'fermi-dirac'";  # ABINIT = fermi-dirac
  $QE_smear[4] = "'marzari-vanderbilt'";  # ABINIT = Marzari cold smearing a = -0.5634
  $QE_smear[5] = "'marzari-vanderbilt'";  # ABINIT = Marzari a = -0.8165
  $QE_smear[6] = "'methfessel-paxton'";  # ABINIT = Methfessel and Paxton PRB 40, 3616 (1989)
  $QE_smear[7] = "'gaussian'";     # ABINIT = Gaussian

  print $fh "&control\n"
        .  "  calculation = \'$inputs{'calctype'}\'\n"
        .  "  prefix = \'$inputs{'prefix'}\'\n"
        .  "  pseudo_dir = \'$inputs{'ppdir'}\'\n"
        .  "  outdir = \'$inputs{'work_dir'}\'\n"
        .  "  wfcdir = \'$inputs{'tmp_dir'}\'\n"
        .  "  tstress = $inputs{'dft.calc_stress'}\n"
        .  "  tprnfor = $inputs{'dft.calc_force'}\n"
        .  "  wf_collect = .true.\n"
        .  "  disk_io = 'low'\n"
        .  "/\n";
  print $fh "&system\n"
        .  "  ibrav = $inputs{'ibrav'}\n"
        .  "  nat = $inputs{'natoms'}\n"
        .  "  ntyp = $inputs{'ntype'}\n"
        .  "  noncolin = $inputs{'noncolin'}\n"
        .  "  lspinorb = $inputs{'spinorb'}\n"
        .  "  ecutwfc = $inputs{'ecut'}\n"
        .  "  occupations = '$inputs{'occtype'}'\n"
        .  "  smearing = $QE_smear[$inputs{'occopt'}]\n"
        .  "  degauss = $inputs{'degauss'}\n"
        .  "  nspin  = $inputs{'nspin'}\n"
        .  "  tot_charge  = $inputs{'tot_charge'}\n"
        .  "  nosym = $inputs{'nosym'}\n"
        .  "  noinv = $inputs{'noinv'}\n";
  unless( $inputs{'dft.functional'} =~ m/default/ )
  {
    print $fh "  input_dft = \'$inputs{'dft.functional'}\'\n";
#    foreach( @exx )
#    {
#      if( $inputs{'dft.functional'} =~ m/$_/i )
#      {
        print $fh "  nqx1 = $inputs{'nqx1'}, nqx2 = $inputs{'nqx2'}, nqx3 = $inputs{'nqx3'}\n";
 #       last;
 #     }
 #   }
  }

  if( $inputs{'print nbands'} > 0 ) # for scf no nbnd is set. 
                                    # Therefore -1 is passed in and nothing is written to the input file
  {
    print $fh "  nbnd = $inputs{'print nbands'}\n";
  }
  if( $inputs{'smag'}  ne "" )
  {
    print $fh "$inputs{'smag'}\n";
  }
  if( $inputs{'ldau'}  ne "" )
  {
    print $fh "$inputs{'ldau'}\n";
  }
  if( $inputs{'qe_scissor'}  ne "" )
  {
    print $fh "$inputs{'qe_scissor'}\n";
  }
  if( $inputs{'ibrav'} != 0 )
  {
    print $fh "  celldim(1) = $inputs{'celldm1'}\n";
  }
  print $fh "/\n"
        .  "&electrons\n"
        .  "  conv_thr = $inputs{'etol'}\n"
        .  "  mixing_beta = $inputs{'mixing'}\n"
        .  "  electron_maxstep = $inputs{'nrun'}\n"
        .  "  startingwfc = \'$inputs{'dft.startingwfc'}\'\n"
        .  "  startingpot = \'$inputs{'dft.startingpot'}\'\n"
        .  "  diagonalization = \'$inputs{'dft.diagonalization'}\'\n";
  if( $inputs{'print nbands'} > 100 && $inputs{'calctype'} =~ m/nscf/i )
  {
    print $fh "  diago_david_ndim = 2\n";
  }
  if( $inputs{'nscfEXX'} == 1 )
  {
    # Since (at the moment) we are loading the SCF density
    #  don't converge the density for the first iteration w/o EXX
    print $fh "  adaptive_thr = .true., conv_thr_init = 1\n";
  }
    print $fh "ATOMIC_SPECIES\n" . $inputs{'atompp'} . "\n";

  if ($inputs{'ibrav'} == 0) {
    print $fh "CELL_PARAMETERS cubic\n" . $inputs{'acell'} . "\n";
  }

#  if( $coord_type =~ m/angst/ )
#  {
#    print $fh "ATOMIC_POSITIONS angstrom\n";
#  }
#  elsif( $coord_type =~ m/bohr/ || $coord_type =~ m/cart/ )
#  {
#    print $fh "ATOMIC_POSITIONS bohr\n";
#  }
#  else
#  {
    print $fh "ATOMIC_POSITIONS crystal\n";
#  }

  print $fh $inputs{'coords'} . "\n";

  print $fh $inputs{'print kpts'};

}


1;
