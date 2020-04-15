#!/usr/bin/perl

# Copyright (C) 2015 - 2017, 2019 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
# Adapted from cnbse_mpi.pl

use strict;
use File::Copy;
use File::Spec::Functions;
use POSIX;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/nbse\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }

my @CommonFiles = ("epsilon", "xmesh.ipt", "k0.ipt", "nbuse.ipt", 
  "metal", "cksshift", "cksstretch", 
  "cnbse.niter", "cnbse.spect_range", "cnbse.broaden", "dft", 
  "para_prefix", "cnbse.strength", "serbse", "core_offset", "avecsinbohr.ipt", 
  "cnbse.solver", "cnbse.gmres.elist", "cnbse.gmres.erange", "cnbse.gmres.nloop", 
  "cnbse.gmres.gprc", "cnbse.gmres.ffff", "cnbse.write_rhs", "spin_orbit", "nspin", 
  "niter", "backf", "aldaf", "bwflg", "bande", "bflag", "lflag", "decut", "spect.h", 
  "gw_control", "gwcstr", "gwvstr", "gwgap", "bse.wvfn" );

my @DFTFiles = ("nelectron", "rhoofr");

my @DenDipFiles = ("kmesh.ipt", "efermiinrydberg.ipt", 
                   "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "nelectron", 
                   "eshift.ipt" );

my @LegacyPrepFiles = ( "masterwfile", "listwfile", "tmels", "enkfile" );

my @OceanPrepFiles = ( "ptmels.dat", "enkfile", "tmels.info" );


my @WFNFiles = ("kmesh.ipt",  "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", 
                "wvfcninfo", "wvfvainfo", "obf_control", "ibeg.h", "q.out", "tmels.info", 
                "val_energies.dat", "con_energies.dat" );

my @ScreenFiles = ( "screen.mode", "screen.lmax", "screen.final.rmax", "screen.final.dr", 
                    "cnbse.rad" );

foreach (@CommonFiles) {
  copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}

foreach (@ScreenFiles)
{
  copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}

open DFT, "dft" or die "Failed to open dft\n";
<DFT> =~ m/(\w+)/ or die;
my $dft_type = $1;
close DTF;
my $obf;
if( $dft_type =~ m/obf/i )
{
  $obf = 1;
}
else
{
  $obf = 0;
}

open IN, "bse.wvfn" or die "Failed to open bse.wvfn\n$!";
my $bseWvfn = <IN>;
close IN;


if( $obf == 1 ) 
{
  foreach (@DFTFiles) {
    copy( "../DFT/$_", $_) or die "Failed to get DFT/$_\n$!";
  }
  foreach (@WFNFiles) {
    copy( "../zWFN/$_", $_ ) or die "Failed to get zWFN/$_\n$!";
  }
  open OUT, ">tmel_selector" or die;
  print OUT "1\n";
  close OUT;
  open OUT, ">enk_selector" or die;
  print OUT "1\n";
  close OUT;
  open OUT, ">bloch_selector" or die;
  print OUT "1\n";
  close OUT;
}
else
{
  foreach (@DFTFiles) {
    copy( "../PREP/$_", $_) or die "Failed to get DFT/$_\n$!";
  }

  foreach (@DenDipFiles) {
    copy( "../PREP/BSE/$_", $_ ) or die "Failed to get PREP/BSE/$_\n$!" ;
  }
  if( $bseWvfn =~ m/legacy/ )
  {
    foreach (@LegacyPrepFiles)
    {
      copy( "../PREP/BSE/$_", $_ ) or die "Failed to get PREP/BSE/$_\n$!" ;
    }
    open OUT, ">tmel_selector" or die;
    print OUT "0\n";
    close OUT;
    open OUT, ">enk_selector" or die;
    print OUT "0\n";
    close OUT;
    open OUT, ">bloch_selector" or die;
    print OUT "0\n";
    close OUT;
    if (-e "../PREP/BSE/u2.dat")
    {
      `ln -s ../PREP/BSE/u2.dat`;
    }
    else
    {
      die "Required file ../PREP/BSE/u2.dat is missing!\n";
    }
  }
  else  #### THIS SECTION STILL NEEDS WORK
        #### MAKE SURE LEGACY AND NEW ARE PULLING IN CORRECT FILES
  {
    foreach (@OceanPrepFiles)
    {
      copy( "../PREP/BSE/$_", $_ ) or die "Failed to get PREP/BSE/$_\n$!" ;
    }
    foreach( "val.u2.dat", "con.u2.dat" )
    {
      if (-e "../PREP/BSE/$_" )
      {
        symlink( "../PREP/BSE/$_", $_ );
      }
      else
      {
        die "Required file ../PREP/BSE/$_ is missing!\n";
      }
    }
    open OUT, ">tmel_selector" or die;
    print OUT "1\n";
    close OUT;
    open OUT, ">enk_selector" or die;
    print OUT "0\n";
    close OUT;
    open OUT, ">bloch_selector" or die;
    print OUT "3\n";
    close OUT;
  }
}


##### Trigger serial bse fallback here
# later we should remove this and fold these two perl scripts together
my $run_serial = 0;
if( -e "serbse" )
{
  open IN, "serbse" or die "$!";
  <IN> =~ m/(\d)/;
  if( $1 == 1 )
  {
    print "Serial BSE requested!!\n";
    print "Parallel required\nCannot continue\n";
    exit 1;
  }
  close IN;
}
unless( -e "$ENV{'OCEAN_BIN'}/ocean.x" )
{
  print "Parallel BSE executable not present in $ENV{'OCEAN_BIN'}\n";
  print "Parallel required\nCannot continue\n";
  exit 1;
}



##### misc other setup
#`echo gmanner > format65`;
copy( "kmesh.ipt", "kgrid" ) or die "$!";
copy( "k0.ipt", "scaledkzero.ipt" ) or die "$!";
copy( "qinunitsofbvectors.ipt", "cksdq" ) or die "$!";

##### qin must be non-zero
open QIN, "qinunitsofbvectors.ipt" or die "Failed to open qinunitsofbvectors.ipt\n$!";
<QIN> =~ m/(\d*\.?\d+)\s+(\d*\.?\d+)\s+(\d*\.?\d+)/ or die "Failed to parse qinunitsofbvectors.ipt\n";
close QIN;
if( abs( $1 ) + abs( $2 ) + abs( $3 ) < 0.0000001 ) 
{
  print "Non-zero q required for VAL\n";
  exit 1;
}

my $para_prefix = "";
if( open PARA_PREFIX, "para_prefix" )
{
  $para_prefix = <PARA_PREFIX>;
  chomp($para_prefix);
  close( PARA_PREFIX);
} else
{
  print "Failed to open para_prefix. Error: $!\nRunning serially\n";
}


system("$ENV{'OCEAN_BIN'}/getnval.x") == 0 or die "Failed to get nval\n";


#####################
if( $obf == 1 )
{
  if( -e "../zWFN/u2par.dat" )
  {
    `ln -s ../zWFN/u2par.dat`;
    `ln -s ../zWFN/ptmels.dat`;
    `cp ../zWFN/tmels.info .`;
    open OUT, ">bloch_type" or die;
    print OUT "new\n";
    close OUT;
  }
  else
  {
    `ln -s ../zWFN/u2.dat`;
  }
}
else  # We are using abi/qe path w/o obfs
{
#  # grab .Psi
#  `touch .Psi`;
#  system("rm .Psi*");
#  open LISTW, "listwfile" or die "Failed to open listwfile\n";
#  while (<LISTW>) 
#  {
#    $_ =~ m/(\d+)\s+(\S+)/ or die "Failed to parse listwfile\n";
#    system("ln -sf ../PREP/BSE/$2 .") == 0 or die "Failed to link $2\n";
#  }  
#
#  print "Running setup\n";
#  system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";
#
#  if (-e "../PREP/BSE/u2.dat")
#  {
#    `ln -s ../PREP/BSE/u2.dat`;
#  }
#  else
#  {
#    print "conugtoux\n";
#    system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
#    print "orthog\n";
#    system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;
#  }
}

open MODE, "screen.mode" or die "Failed to open screen.mode\n$!";
if( <MODE> =~ m/grid/i )
{
  open XMESH, "xmesh.ipt" or die "Failed to open xmesh.ipt\n$!";
  <XMESH> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse xmesh.ipt\n$_";
  my $nx = $1*$2*$3;
  close XMESH;

  open RAD, "cnbse.rad" or die "Failed to open cnbse.rad\n$!";
  <RAD> =~ m/(\d+\.?\d*)/ or die "Failed to parse cnbse.rad\n$_";
  my $rad = sprintf "zR%03.2f", $1;
  close RAD;

  for( my $i=1; $i <= $nx; $i++ )
  {
    my $xdir = sprintf "x%06i", $i;
    my $targfile = catfile( updir(), "SCREEN", $xdir, $rad, "rpot" );
    my $destfile = sprintf "rpotx%06i", $i;
    copy( $targfile, $destfile ) or die "Failed to copy $targfile to $destfile\n";
  }

}
close MODE;

my $gamma0 = `cat cnbse.broaden`;
chomp($gamma0);

open TMPFILE, "cnbse.niter" or die "Failed to open cnbse.niter\n$!";
<TMPFILE> =~ m/(\d+)/ or die "Failed to parse cnbse.niter";
my $num_haydock_iterations = $1;
close TMPFILE;

open TMPFILE, "cnbse.strength" or die "Failed to open cnbse.strength\n$!";
<TMPFILE> =~ m/([0-9]*\.?[0-9]+)/ or die "Failed to parse cnbse.strength\n";
my $interaction_strength = $1;
close TMPFILE;

#write mode file
open TMPFILE, ">mode" or die "Failed to open mode for writing\n$!";
print TMPFILE "$interaction_strength    $num_haydock_iterations\n";
close TMPFILE;

open CKSNORM, ">cks.normal" or die "$!\n";
print CKSNORM ".true.\n";
close CKSNORM;


open RUNLIST, ">runlist";
print RUNLIST "1\n0  0  0  __  __  0  0  VAL\n";
close RUNLIST;

# Set up for running the valence pathway for RIXS

print "\nSetting up valence section\n";


`tail -n 1 rhoofr > nfft`;

open INFILE, ">bse.in" or die "Failed to open bse.in\n";
print INFILE "0 0 0\n";
print INFILE "0 0 0\n";
my $spectrange = `cat spect.h`;
chomp($spectrange);
my $num_haydock_iterations = `cat niter`;
chomp($num_haydock_iterations);


print INFILE "hay\n";
print INFILE "$num_haydock_iterations  $spectrange  $gamma0  0.000\n";
close INFILE;

`cat bse.in`;

`echo 0 0 0 > ZNL`;

print "Running valence\n";
print "$para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log";
system("$para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log") == 0 or die "Failed to finish\n";
print "\n";


exit 0;

