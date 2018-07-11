#!/usr/bin/perl
# Copyright (C) 2014, 2016, 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

# Adapted from qe_dendip.pl this is the jdftx only version
use strict;
use File::Copy;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/jdftx_dendip\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################
my $RunDenDip = 0;


my $stat = 0;
$stat = 1 if (-e "done");
`rm -f done`;
my $oldden = 0;
$oldden = 1 if (-e "../DFT/old");


my @DFTFiles = ( "scf.out", "rhoofr", "nfft", "brange.ipt" );
my @CommonFiles = (  "screen.nkpt", "nkpt", "qinunitsofbvectors.ipt", "avecsinbohr.ipt", "dft",
                    "nspin", "xmesh.ipt", "dft.split", "prefix", "calc" );
                    

foreach (@DFTFiles) {
  system("cp ../DFT/$_ .") == 0 or die "Failed to copy $_\n";
}
foreach (@CommonFiles)
{
  system("cp ../Common/$_ .") == 0 or die "Failed to copy $_\n";
}
system("mv nkpt bse.nkpt") == 0 or die "Failed to rename nkpt $_\n";
print "$stat  $oldden\n";

#JTV at the moment assume jdftx is always split
my $split_dft = 1;
open IN, "dft.split" or die "Failed to open dft.split\n$!";
if( <IN> =~ m/t/i )
{
  $split_dft = 1;
  print "DFT has been split into to runs\n";
}

# PREP density
system("$ENV{'OCEAN_BIN'}/nelectron.x") == 0 
  or die "Failed to run nelectron.x\n";

system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0
  or die "Failed to run bvecs.pl\n";

system("$ENV{'OCEAN_BIN'}/gvecs2.pl") == 0
  or die "Failed to run gvecs2.pl\n";


open IN, "calc" or die "Failed to open calc\n";
<IN> =~ m/(\w+)/ or die "Failed to parse calc\n";
my $calc = $1;
close IN;
my $run_screen = 1;
if( $calc =~ m/val/i )
{
  $run_screen = 0;
}

my $rundir;
my @nkpt;

if( $run_screen == 1 )
{
  die "John hasn't programed run screen yet for jdftx\n";
}


open NKPT, "bse.nkpt" or die "Failed to open nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse nkpt\n";
@nkpt = ($1, $2, $3);
close NKPT;
my $nkpt = $nkpt[0]*$nkpt[1]*$nkpt[2];


my $run_bse = 1;
if( $run_bse == 1 )
{
  
  `rm -r BSE` if (-e "BSE");
  mkdir "BSE";
  chdir "BSE";

  open OUT, ">nkpts" or die "Failed to open nkpts for writing\n$!";
  print OUT "$nkpt\n";
  close OUT;

  my @fileList = ( "dft.split", "brange.ipt", "nspin", "nelectron", "xmesh.ipt", "avecsinbohr.ipt", 
                   "qinunitsofbvectors.ipt", "bvecs", "dft" );
  foreach (@fileList)
  {
    copy( "../$_", "$_" ) or die "Failed to copy $_\n$!";
  }
  copy( "../bse.nkpt", "kmesh.ipt" ) or die "Failed to copy bse.nkpt\n";

  # do energies here
  $rundir = sprintf("../../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);
  print "Looking for $rundir/nscf.eigenvals\n";
  copy( "$rundir/nscf.eigenvals", "nscf.eigenvals" ) or die "$!";
  `ln -s $rundir unshifted`;
  if( $split_dft == 1 )
  {
    $rundir = sprintf("../../DFT/%03u%03u%03u.shift", $nkpt[0], $nkpt[1], $nkpt[2]);
    print "Looking for $rundir/nscf.eigenvals\n";
    copy( "$rundir/nscf.eigenvals", "nscf.eigenvals.shift" ) or die "$!";
    `ln -s $rundir shifted`;
  }

  print "Converting energies\n";
  system("$ENV{'OCEAN_BIN'}/jdftx_energy.x") == 0 
    or die "Failed to run jdftx_energy.x\n";

  system("$ENV{'OCEAN_BIN'}/ofermi.pl") == 0
    or die "Failed to run ofermi.pl\n";

  #FOR NOW!!!
  copy( "efermiinrydberg2.ipt", "efermiinrydberg.ipt" ) or die "$!";

  # end energies

  open OUT, ">masterwfile" or die "Failed to open masterwfile for writing\n$!";
  print OUT "$nkpt\n";
  close OUT;

  open OUT, ">listwfile" or die "Failed to open listwfile for writing\n$!";
  for( my $i = 0; $i < $nkpt; $i++ )
  {
    printf OUT "%6i kpoint%i\n", ($i+1), $i;
  }
  close OUT;

  print "Running setup\n";
  system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0
    or die "Failed to run setup\n";

  print "split_conugtoux\n";
  system("$ENV{'OCEAN_BIN'}/shifted_conugtoux.x > conugtoux.log");# == 0 or die;
  print "orthog\n";
  system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;

  chdir "../";
}

unless ($stat && $oldden) 
{
  system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0
    or die "Failed to run rhoofg.x\n";
  `wc -l rhoG2 > rhoofg`;
  `sort -n -k 6 rhoG2 >> rhoofg`;
}

`touch done`;

exit 0;
