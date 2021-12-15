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
use File::Spec;
use File::Copy;

#sub QErunNSCF {}
#sub QEsetupNSCF {}

sub ABIrunDensity 
{
  my $hashRef = $_[0];

  open my $input, ">", "scf.in" or die "Failed to open scf.in.\n$!";
  open my $files, ">", "scf.files" or die "Failed to open scf.files\n$!";

  my $calcFlag = 1;

  # Modify QEprintInput to take two different hashes and a string 
  ABIprintInput( $input, $files, $hashRef, $hashRef->{'scf'}, $calcFlag, "scf" ) ;

  close $input;
  close $files;

  ABIrunABINIT( $hashRef->{'computer'}->{'para_prefix'}, "scf.files", "scf.log" );

  my $errorCode = ABIparseOut( "scf.log", $hashRef->{'scf'} );
  return $errorCode;
#  return 1;
}


#sub QErunTest {}

sub ABIparseOut 
{
  my( $outFile, $hashRef ) = @_;

  my $errorCode = 1;

  my $fermi = 'no';

  open SCF, "<", $outFile or die "Failed to open $outFile\n$!";

  while( my $scf_line = <SCF> )
  { 
    if( $scf_line  =~  m/Fermi \(or HOMO\) energy \(hartree\) =\s+([+-]?\d+\.?\d+)/ )
    {
      $fermi = $1 * 2;
      $hashRef->{'fermi'} = $fermi;
    }
    elsif( $scf_line =~ m/Version\s+(\d+\.\d+\.\d+)\s+of ABINIT/ ) {
      $hashRef->{'version'} = $1;
    }
    elsif( $scf_line =~ m/wall=\s+(\d+\.\d+[eE]?[+-]?\d?)/ ) {
      $hashRef->{'time'} = $1*1;
    }
    elsif( $scf_line =~ m/etotal\s+(-?\d+\.\d+[eE]?[-+]?\d*)/ )
    { 
      $hashRef->{'etot'} = $1*1;
    }
    elsif( $scf_line =~ m/is converged/ )
    {
      $errorCode = 0;
    }
    elsif( $scf_line =~ m/nelect=\s+(\d+\.?\d*)/ ) 
    {
      $hashRef->{'nelec'} = $1*1;
    }

  }

  close SCF;

  return $errorCode;
}

sub ABIparseDensityPotential
{
  my ( $hashRef, $type ) = @_;

  my $infile;
  my $datafile;
  my $outfile;
  my $logfile;
  if( $type eq 'density' ) {
    $infile = "cut3d.in";
    $datafile = "SCFx_DEN";
    $outfile = "rhoofr";
    $logfile = 'cut3d.log'
  } elsif( $type eq 'potential' ) {
    $infile = "cut3d2.in";
    $datafile = "SCFx_POT";
    $outfile = "potofr";
    $logfile = 'cut3d2.log'
  } else {
    return 1;
  }
  
  
  my $AbiVersion = 9;
  if( $hashRef->{'scf'}->{'version'} =~ m/(\d+)\.(\d+)\.(\d+)/ ) {
    $AbiVersion = $1;
  }

  open CUTIN, ">", $infile or die "Failed to open $infile for writing.\n$!\n";
  if( $AbiVersion <= 7 )
  {
    if( $hashRef->{'general'}->{'nspin'} == 2 )
    {
      print CUTIN "$datafile\n1\n0\n6\n$outfile\n0\n";
    }
    else
    {
      print CUTIN "$datafile\n1\n6\n$outfile\n0\n";
    }
  }
  else
  {
    if(  $hashRef->{'general'}->{'nspin'}  == 2 )
    {
      print CUTIN "$datafile\n0\n6\n$outfile\n0\n";
    }
    else
    {
      print CUTIN "$datafile\n6\n$outfile\n0\n";
    }
  }
  close CUTIN;
  

  print "$hashRef->{'computer'}->{'ser_prefix'} $ENV{'OCEAN_CUT3D'} < $infile > $logfile 2>&1\n";
  unless (system("$hashRef->{'computer'}->{'ser_prefix'} $ENV{'OCEAN_CUT3D'} < $infile > $logfile 2>&1") == 0)
  {
    print "Failed to run cut3d\n";
    return 2;
  }
      
  return 0;
}


#sub QEPoolControl {}

sub ABIrunABINIT
{
  my( $prefix, $in, $out ) = @_;

  print  "$prefix $ENV{'OCEAN_ABINIT'} < $in > $out 2>&1\n";
  system("$prefix $ENV{'OCEAN_ABINIT'} < $in > $out 2>&1");
  
}

#sub QErunPP {} 

sub ABIprintInput 
{
  my ($input, $files, $generalRef, $specificRef, $calcFlag, $baseString ) = @_;


  print $files $baseString . ".in\n" . $baseString . ".out\n" 
          . uc($baseString) . "\n" . uc($baseString) . "x\n" 
          . uc($baseString) . "xx\n";
  foreach my $p ( @{$generalRef->{'psp'}->{'pp_list'}}) {
    print $files catfile( $generalRef->{'psp'}->{'ppdir'}, $p ) . "\n";
  }
  

  print $input "symmorphi 0\nautoparal 1\nchksymbreak 0\n"
             . "rprim\n";
  printf $input "% .16g  % .16g  % .16g\n% .16g  % .16g  % .16g\n% .16g  % .16g  % .16g\n", 
        $generalRef->{'structure'}->{'avecs'}[0][0], $generalRef->{'structure'}->{'avecs'}[0][1],
        $generalRef->{'structure'}->{'avecs'}[0][2], $generalRef->{'structure'}->{'avecs'}[1][0],
        $generalRef->{'structure'}->{'avecs'}[1][1], $generalRef->{'structure'}->{'avecs'}[1][2],
        $generalRef->{'structure'}->{'avecs'}[2][0], $generalRef->{'structure'}->{'avecs'}[2][1],
        $generalRef->{'structure'}->{'avecs'}[2][2];

  printf $input "ntypat %i\nznucl", scalar @{$generalRef->{'structure'}->{'znucl'}};
  foreach my $z (@{$generalRef->{'structure'}->{'znucl'}}) {
    printf $input " %i", $z;
  }
  printf $input "\nnatom %i\ntypat", scalar @{$generalRef->{'structure'}->{'typat'}};
  foreach my $z (@{$generalRef->{'structure'}->{'typat'}}) {
    printf $input " %i", $z;
  }
  print $input "\nxred\n";
  for( my $i = 0; $i < scalar @{$generalRef->{'structure'}->{'typat'}}; $i++ ) {
    printf $input " % .16g  % .16g  % .16g\n", $generalRef->{'structure'}->{'xred'}[$i][0], 
        $generalRef->{'structure'}->{'xred'}[$i][1], $generalRef->{'structure'}->{'xred'}[$i][2];
  }

  printf $input "ecut %g Ry\n", $generalRef->{'general'}->{'ecut'};
  printf $input "nstep %i\n", $generalRef->{'general'}->{'nstep'};
  if( $generalRef->{'epsilon'}->{'method'} eq "input" ) {
    printf $input "diemac %g\n", $generalRef->{'structure'}->{'epsilon'};
  } else {  # If no input dielectric then need to guess
    print $input "diemac 10\n";
  }

  printf $input "occopt %i\ntsmear %g Ry\n", 
      $generalRef->{'general'}->{'occopt'}, $generalRef->{'general'}->{'degauss'};

  print $input "npfft 1\n";
  printf $input "charge %g\n", $generalRef->{'general'}->{'tot_charge'};
  
  printf $input "nsppol %i\n", $generalRef->{'general'}->{'nspin'};
  if( $generalRef->{'general'}->{'nspin'} != 1 ) {
    printf $input "spinat\n%s", $generalRef->{'general'}->{'smag'};
  }

  if( $calcFlag == 1 ) {
    printf $input "fband %g\n", $generalRef->{'general'}->{'fband'};
    print $input "prtden 1\nprtpot 1\nkptopt 1\n";
    printf $input "ngkpt %i %i %i\nnshiftk 1\nshiftk", $specificRef->{'kmesh'}[0], 
          $specificRef->{'kmesh'}[1], $specificRef->{'kmesh'}[2];
    for( my $i = 0; $i < 3; $i++ ) {
      if( $specificRef->{'kshift'}[$i]==1 ) {
        print $input " 0.5";  
      } else {
        print $input " 0.0";
      }
    }
    print $input "\n";

    printf $input "toldfe  %g\n", $specificRef->{'toldfe'};
    if( $generalRef->{'calc_stress'} ) {
      print $input "optstress 1\n";
    } else {
      print $input "optstress 0\n";
    }

    if( $generalRef->{'calc_force'} ) {
      print $input "optforces 1\n";
    } else {
      print $input "optforces 0\n";
    }

  } else {

  }

}

#sub kptGen {}
#sub QEparseEnergies {}
#sub QEparseEnergies62  {}


1;
