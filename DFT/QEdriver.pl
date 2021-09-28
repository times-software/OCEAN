#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use POSIX;

sub QErunDensity
{
  my $hashRef = $_[0];

  # Make SCF k-point string
  my $kptString;
  if( $hashRef->{'scf'}->{'kmesh'}[0] * $hashRef->{'scf'}->{'kmesh'}[1] * $hashRef->{'scf'}->{'kmesh'}[2] == 1 
      && $hashRef->{'scf'}->{'kshift'}[0] == 0 && $hashRef->{'scf'}->{'kshift'}[1] == 0 
      && $hashRef->{'scf'}->{'kshift'}[2] == 0 ) {
   $kptString = "K_POINTS gamma\n";
  } else { 
    $kptString = "K_POINTS automatic\n" 
               . sprintf "%i %i %i %i %i %i\n", $hashRef->{'scf'}->{'kmesh'}[0], $hashRef->{'scf'}->{'kmesh'}[1],
                                                $hashRef->{'scf'}->{'kmesh'}[2], $hashRef->{'scf'}->{'kshift'}[0],
                                                $hashRef->{'scf'}->{'kshift'}[1], $hashRef->{'scf'}->{'kshift'}[2];
  }
  

  # open input file
  open my $fh, ">scf.in" or die "Failed to open scf.in.\n$!";

  my $calcFlag = 1;

  # Modify QEprintInput to take two different hashes and a string 
  QEprintInput( $fh, $hashRef, $hashRef->{'scf'}, $calcFlag, $kptString ) ;

  close $fh;

  # test run for k-points
  print "Testing parallel QE execution\n";
  open TMP, '>', "system.EXIT" or die "Failed to open file system.EXIT\n$!";
  close TMP;

  QErunPW( $hashRef->{'general'}->{'redirect'}, $hashRef->{'computer'}->{'ser_prefix'}, "" , "scf.in", "test.out" );
  
  # parse test results
  my $npool = 1;
  my $ncpus = $hashRef->{'computer'}->{'ncpus'};
  if( open TMP, "test.out" )
  {
    my $actualKpts = -1;
    my $numKS = -1;
    my $mem = -1;
    while (<TMP>)
    {
      if( $_ =~ m/number of Kohn-Sham states=\s+(\d+)/ )
      {
        $numKS = $1;
      }
      elsif( $_ =~ m/number of k points=\s+(\d+)/ )
      {
        $actualKpts = $1;
      }
      elsif( $_ =~ m/Estimated max dynamical RAM per process\s+>\s+(\d+\.\d+)/ )
      {
        $mem = $1;
      }
      last if( $actualKpts > 0 && $numKS > 0 && $mem > 0 );
    }
    close TMP;

    if( $actualKpts == -1 )
    {
      print "Had trouble parsing test.out\nDidn't find number of k points\n";
    }
    else
    {
      ($ncpus, $npool) = QEPoolControl( $actualKpts, $numKS, $mem, $hashRef->{'computer'} );
    }
  }
  else
  {
    print "Had trouble parsing test.out\n. Will attempt to continue.\n";
  }

  

  # full run

  print "$ncpus  $npool\n";

  my $prefix = $hashRef->{'computer'}->{'para_prefix'};
  $prefix =~ s/$hashRef->{'computer'}->{'ncpus'}/$ncpus/ if( $hashRef->{'computer'}->{'ncpus'} != $ncpus );

  my $cmdLine = " -npool $npool ";

  QErunPW( $hashRef->{'general'}->{'redirect'}, $prefix, $cmdLine, "scf.in", "scf.out" );


  my ($errorCode, $energy ) = QEparseOut( "scf.out" );

  return ( $errorCode, $energy );
}


sub QEparseOut
{
  my $outFile = $_[0];

  my $errorCode = 1;
  my $totEnergy = "";

  open OUTFILE, "<", $outFile or die "Failed to open $outFile\n";
  while( my $line=<OUTFILE> )
  {
    if( $line =~ m/!\s+total energy\s+=\s+(-?\d+\.\d+)/ )
    {
      $totEnergy = $1;
#      print ">>>>> $totEnergy\n";
    }
    elsif( $line =~ m/convergence has been achieved/ )
    {
      $errorCode = 0;
    }
    elsif( $line =~ m/End of band/ )
    {
      $errorCode = 0;
    }
  }
  close OUTFILE;
  return ( $errorCode, $totEnergy );
}


# Figure out how many pools to run with
sub QEPoolControl
{
  my ( $actualKpts, $numKS, $mem, $hashRef ) = @_;

  my $maxMem = 2000;
  my $minPool = 0.9;
  $minPool = $mem / $maxMem if( $mem > $maxMem );

#  print "$hashRef->{'ncpus'}\n";

  my %cpusAndPools;
  # j is number of processors
  for( my $j = $hashRef->{'ncpus'}; $j > 0; $j-- )
  {
    # i is number of pools
    for( my $i = 1; $i <= sqrt( $j+.1); $i++ )
    {
      # gotta be an even factor
      next if ( $j % $i );
      # pool size (j/i) must be larger than minimum size
      next if ( $j / $i < $minPool );

      next if ( $numKS > 0 && $j / $i > $numKS );

      my $cpuPerPool = $j / $i;
      my $kPerPool = ceil( $actualKpts / $i );
      my $cost = $kPerPool / ( $cpuPerPool * ( 0.98**$cpuPerPool ) );
#      print "$j  $i  $kPerPool  $cost\n";
      $cpusAndPools{ $cost } = [ $j, $i ];
    }
  }

  my $ncpus; my $npool;
  foreach my $i (sort {$b <=> $a} keys %cpusAndPools )
  {
#    print "$i  $cpusAndPools{$i}[0]  $cpusAndPools{$i}[1]\n";
    $ncpus = $cpusAndPools{$i}[0];
    $npool = $cpusAndPools{$i}[1];
  }

  return ( $ncpus, $npool );
}

# subroutine to handle running each QE step
# make one of these for pp and ph too!
sub QErunPW
{
  my( $redirect, $prefix, $cmdLine, $in, $out ) = @_;

  if( $redirect )
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine < $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine < $in > $out 2>&1");
  }
  else
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine -inp $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine -inp $in > $out 2>&1");
  }
}


# Subroutine to print a QE-style input file
#   Pass in a file handle and a hashtable
sub QEprintInput
{
  my ($fh, $generalRef, $specificRef, $calcFlag, $kptString ) = @_;

  my $calc;
  my $tstress = '.false';
  my $tprnfor = '.false.';
  my $nosyminv;
  my $startingPot;
  if( $calcFlag )
  {
    $calc = 'scf';
    $tstress = 'true' if( $generalRef->{'general'}->{'calc_stress'} );
    $tprnfor = 'true' if( $generalRef->{'general'}->{'calc_force'} );
    $nosyminv = 'false';
    $startingPot = 'atomic';
  }
  else
  {
    $calc = 'nscf';
    $tstress = 'false';
    $tprnfor = 'false';
    $nosyminv = 'true';
    $startingPot = 'file';
  }

  my $noncolin = '.false.';
  $noncolin = '.true.' if( $generalRef->{'general'}->{'noncolin'} );
  my $spinorb = '.false.';
  $spinorb = '.true.' if( $generalRef->{'general'}->{'spinorb'} );

  my $nbnd;
  if( exists $specificRef->{'nbands'} )
  {
    $nbnd = $specificRef->{'nbands'};
  } else {
    $nbnd = $generalRef->{'structure'}->{'valence_electrons'} / 2
          + $generalRef->{'general'}->{'fband'} * scalar @{$generalRef->{'structure'}->{'xred'}};
    $nbnd = ceil( $nbnd );
    $nbnd++ if( $nbnd % 2 );
  }
    

  # Array of QE names for smearing by occopt
  my @QE_smear;
  $QE_smear[1] = "'gaussian'";     # Still need to fix to be insulator
  $QE_smear[3] = "'fermi-dirac'";  # ABINIT = fermi-dirac
  $QE_smear[4] = "'marzari-vanderbilt'";  # ABINIT = Marzari cold smearing a = -0.5634
  $QE_smear[5] = "'marzari-vanderbilt'";  # ABINIT = Marzari a = -0.8165
  $QE_smear[6] = "'methfessel-paxton'";  # ABINIT = Methfessel and Paxton PRB 40, 3616 (1989)
  $QE_smear[7] = "'gaussian'";     # ABINIT = Gaussian

  my $occopt = "fixed";
  $occopt = "smearing" if( $generalRef->{'general'}->{'occopt'} != 1 );

  print $fh "&control\n"
        .  "  calculation = \'$calc\'\n"
        .  "  prefix = \'system\'\n"
        .  "  pseudo_dir = \'$generalRef->{'psp'}->{'ppdir'}\'\n"
        .  "  outdir = \'Out\'\n"
        .  "  wfcdir = \'$generalRef->{'general'}->{'tmp_dir'}\'\n"
        .  "  tstress = $tstress\n"
        .  "  tprnfor = $tprnfor\n"
        .  "  wf_collect = .true.\n"
        .  "  disk_io = 'low'\n"
        .  "/\n";
  print $fh "&system\n"
        .  "  ibrav = 0\n"
        .  "  nat = " . scalar @{$generalRef->{'structure'}->{'xred'}} . "\n"
        .  "  ntyp = " . scalar @{$generalRef->{'structure'}->{'typat'}} . "\n"
        .  "  noncolin = $noncolin\n" 
        .  "  lspinorb = $spinorb\n"
        .  "  ecutwfc = $generalRef->{'general'}->{'ecut'}\n"
        .  "  occupations = \'$occopt\'\n"
        .  "  smearing = $QE_smear[$generalRef->{'general'}->{'occopt'}]\n"
        .  "  degauss = $generalRef->{'general'}->{'degauss'}\n"
        .  "  nspin  = $generalRef->{'general'}->{'nspin'}\n"
        .  "  tot_charge  = $generalRef->{'general'}->{'tot_charge'}\n"
        .  "  nosym = $nosyminv\n"
        .  "  noinv = $nosyminv\n";
  unless( $generalRef->{'general'}->{'functional'} =~ m/default/ )
  {
    print $fh "  input_dft = \'$generalRef->{'general'}->{'functional'}\'\n";
    print $fh "  nqx1 = $generalRef->{'general'}->{'exx'}->{'qmesh'}[0],"
             . " nqx2 = $generalRef->{'general'}->{'exx'}->{'qmesh'}[1],"
             . " nqx3 = $generalRef->{'general'}->{'exx'}->{'qmesh'}[2]\n";
  }

  print $fh "  nbnd = $nbnd\n";

  if( $generalRef->{'general'}->{'smag'}  ne "" )
  {
    print $fh "$generalRef->{'general'}->{'smag'}\n";
  }
  if( $generalRef->{'general'}->{'ldau'}->{'enable'} )
  {
    print $fh "lda_plus_u = true\n" 
            . "lda_plus_u_kind = $generalRef->{'general'}->{'ldau'}->{'lda_plus_u_kind'}\n"
            . "U_projection_type = $generalRef->{'general'}->{'ldau'}->{'U_projection_type'}\n";
    print $fh "Hubbard_U = $generalRef->{'general'}->{'ldau'}->{'Hubbard_U'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_U'} ne "" );
    print $fh "Hubbard_V = $generalRef->{'general'}->{'ldau'}->{'Hubbard_V'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_V'} ne "" );
    print $fh "Hubbard_J = $generalRef->{'general'}->{'ldau'}->{'Hubbard_J'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_J'} ne "" );
    print $fh "Hubbard_J0 = $generalRef->{'general'}->{'ldau'}->{'Hubbard_J0'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_J0'} ne "" );
  }
#  if( $inputs{'qe_scissor'}  ne "" )
#  {
#    print $fh "$inputs{'qe_scissor'}\n";
#  }
#  if( $inputs{'ibrav'} != 0 )
#  {
#    print $fh "  celldim(1) = $inputs{'celldm1'}\n";
#  }
  print $fh "/\n"
        .  "&electrons\n"
        .  "  conv_thr = $specificRef->{'toldfe'}\n"
        .  "  mixing_beta = $generalRef->{'general'}->{'mixing'}\n"
        .  "  electron_maxstep = $generalRef->{'general'}->{'nstep'}\n"
        .  "  startingwfc = \'$generalRef->{'general'}->{'startingwfc'}\'\n"
        .  "  startingpot = \'$startingPot\'\n"
        .  "  diagonalization = \'$generalRef->{'general'}->{'diagonalization'}\'\n"
        .  "/\n";
#  if( $inputs{'print nbands'} > 100 && $inputs{'calctype'} =~ m/nscf/i )
#  {
#    print $fh "  diago_david_ndim = 2\n";
#  }
#  if( $inputs{'nscfEXX'} == 1 )
#  {
#    # Since (at the moment) we are loading the SCF density
#    #  don't converge the density for the first iteration w/o EXX
#    print $fh "  adaptive_thr = .true., conv_thr_init = 1\n";
#  }

  print $fh "ATOMIC_SPECIES\n";# . $inputs{'atompp'} . "\n";
  for( my $a=0; $a < scalar @{$generalRef->{'structure'}->{'zsymb'}}; $a++ )
  {
    printf $fh "%4s  0.0  %s\n", $generalRef->{'structure'}->{'zsymb'}[$a], $generalRef->{'psp'}->{'pp_list'}[$a];
  }

#  if ($inputs{'ibrav'} == 0) {
    print $fh "CELL_PARAMETERS cubic\n";  # . $inputs{'acell'} . "\n";
    for( my $i = 0; $i<3; $i++ )
    {
      printf $fh " %.16f    %.16f   %.16f \n", $generalRef->{'structure'}->{'avecs'}[$i][0], 
          $generalRef->{'structure'}->{'avecs'}[$i][1], $generalRef->{'structure'}->{'avecs'}[$i][2];
    }
#  }

  print $fh "ATOMIC_POSITIONS crystal\n";

  for( my $a=0; $a < scalar @{$generalRef->{'structure'}->{'xred'}}; $a++ )
  {
    printf $fh "%4s  %.16f    %.16f    %.16f \n", 
          $generalRef->{'structure'}->{'zsymb'}[$generalRef->{'structure'}->{'typat'}[$a]-1], 
          $generalRef->{'structure'}->{'xred'}[$a][0],
          $generalRef->{'structure'}->{'xred'}[$a][1],
          $generalRef->{'structure'}->{'xred'}[$a][2];
  }

  print $fh $kptString;

}


1;
