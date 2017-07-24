#!/usr/bin/perl

# Copyright (C) 2015 - 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
# Adapted from cnbse_mpi.pl

use strict;
use File::Copy;
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
  "gw_control", "gwcstr", "gwvstr", "gwgap" );

my @DFTFiles = ("nelectron", "rhoofr");

my @DenDipFiles = ("kmesh.ipt", "masterwfile", "listwfile", "efermiinrydberg.ipt", 
                   "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "tmels", "nelectron", 
                   "eshift.ipt" );

my @WFNFiles = ("kmesh.ipt",  "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", 
                "wvfcninfo", "wvfvainfo", "obf_control", "ibeg.h", "q.out", "tmels.info", 
                "val_energies.dat", "con_energies.dat" );

foreach (@CommonFiles) {
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
  open OUT, ">tmel_selector" or die;
  print OUT "0\n";
  close OUT;
  open OUT, ">enk_selector" or die;
  print OUT "0\n";
  close OUT;
  open OUT, ">bloch_selector" or die;
  print OUT "0\n";
  close OUT;

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
  # grab .Psi
  `touch .Psi`;
  system("rm .Psi*");
  open LISTW, "listwfile" or die "Failed to open listwfile\n";
  while (<LISTW>) 
  {
    $_ =~ m/(\d+)\s+(\S+)/ or die "Failed to parse listwfile\n";
    system("ln -sf ../PREP/BSE/$2 .") == 0 or die "Failed to link $2\n";
  }  

  print "Running setup\n";
  system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";

  if (-e "../PREP/BSE/u2.dat")
  {
    `ln -s ../PREP/BSE/u2.dat`;
  }
  else
  {
    print "conugtoux\n";
    system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
    print "orthog\n";
    system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;
  }
}

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
print "time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log";
system("time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log") == 0 or die "Failed to finish\n";


exit 0;

