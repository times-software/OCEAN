#!/usr/bin/perl
# Copyright (C) 2014, 2016 - 2019 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Path qw( rmtree );
use Cwd 'abs_path';
use File::Compare;
use File::Copy;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/qe_dendip\.pl/;
#  $ENV{"OCEAN_BIN"} = $1;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################
my $RunDenDip = 0;

my $stat = 0;
$stat = 1 if (-e "done");
`rm -f done`;
my $oldden = 0;
$oldden = 1 if (-e "../DFT/old");


my @QEFiles     = ( "rhoofr", "efermiinrydberg.ipt" );
my @CommonFiles = ( "screen.nkpt", "nkpt", "qinunitsofbvectors.ipt", "avecsinbohr.ipt", "dft", 
                    "nspin", "xmesh.ipt", "dft.split", "prefix", "calc", "screen.wvfn", "screen.legacy", 
                    "screen.mode");

foreach (@QEFiles) {
  system("cp ../DFT/$_ .") == 0 or die "Failed to copy $_\n";
}
foreach (@CommonFiles) {
  system("cp ../Common/$_ .") == 0 or die "Failed to copy $_\n";
}
system("mv nkpt bse.nkpt") == 0 or die "Failed to rename nkpt $_\n";
print "$stat  $oldden\n";
unless ($stat && $oldden) {


  `tail -n 1 rhoofr > nfft`;

  system("$ENV{'OCEAN_BIN'}/nelectron.x") == 0 
    or die "Failed to run nelectron.x\n";

  system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0
    or die "Failed to run bvecs.pl\n";

  system("$ENV{'OCEAN_BIN'}/gvecs2.pl") == 0
    or die "Failed to run gvecs2.pl\n";
}


open IN, "calc" or die "Failed to open calc\n";
<IN> =~ m/(\w+)/ or die "Failed to parse calc\n";
my $calc = $1;
close IN;

open IN, "screen.mode" or die "Failed to open screen.mode\n";
<IN> =~ m/(\w+)/ or die "Failed to parse screen.mode\n";
my $screen_mode = $1;
close IN;
my $run_screen = 1;
if( $calc =~ m/val/i )
{
  $run_screen = 0 unless( $screen_mode =~ m/grid/i );
}

open IN, "screen.wvfn" or die "Failed to open screen.wvfn\n$!";
my $screenWvfn = <IN>;
close IN;

open IN, "screen.legacy" or die "Failed to open screen.legacy\n$!";
my $screenLegacy = <IN>;
chomp $screenLegacy;
close IN;

my $rundir;
my @nkpt;

## process screen wf files ##
if( $run_screen == 1 )
{
  open NKPT, "screen.nkpt" or die "Failed to open screen.nkpt";
  <NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse screen.nkpt\n";
  @nkpt = ($1, $2, $3);
  close NKPT;
  $rundir = "../DFT/SCREEN";

  unless( -e "PAW/done" && -e "${rundir}/old" ) {
  `rm -r PAW` if (-e "PAW");
  mkdir "PAW"; 
  chdir "PAW";

  open NKPTS, ">nkpts" or die "Failed to open nkpts for writing\n";
  print NKPTS $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
  close NKPTS;


  foreach ("kmesh.ipt", "brange.ipt", "qinunitsofbvectors.ipt" ) {
    system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
  }
  #`cp ../qinunitsofbvectors.ipt .`;
  `cp ../bvecs .`;
  `cp ../dft .`;
  `cp ../nspin .`;
  `cp ../${rundir}/umklapp .`;
  `cp ../prefix .`;

  my $prefix;
  open PREFIX, "prefix";
  $prefix = <PREFIX>;
  close (PREFIX);
  chomp( $prefix );


  #`cp -r ../${rundir}/Out .`;
  if( -l "Out" )  # Out is an existing link
  {
    unlink "Out" or die "Problem cleaning old 'Out' link\n$!";
  }
  elsif(  -d "Out" ) #or Out is existing directory
  {
    rmtree( "Out" );
  }
  elsif( -e "Out" ) #or Out is some other file
  {
    unlink "Out";
  }
  print "../$rundir/Out\n";
  symlink ("../$rundir/Out", "Out") == 1 or die "Failed to link Out\n$!";

  # New methods for skipping wfconvert
  if( $screenWvfn =~ m/qe54/ && $screenLegacy == 0 )
  {
    print "Don't convert DFT. Using new method for screening wavefunctions\n";
    `touch listwfile masterwfile enkfile`;
  }
  else  # old method, run wfconvert
  {
    print "$screenWvfn $screenLegacy\n";
    print "$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml\n";
    system("$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml") == 0 
      or die "Failed to run qe_data_file.pl\n$!";

    system("$ENV{'OCEAN_BIN'}/wfconvert.x") == 0 
      or die "Failed to run wfconvert.x\n$!";
  }


  `touch done`;
  chdir "../";
  }
  else {
    `touch PAW/old`;
    print  "Nothing needed for SCREEN wfns\n";
  }


  print "Done with SCREEN files\n";
}
else
{
  print "Nothing needed for SCREEN wvfn\n";
}


## process bse wf files ##
my $split_dft = 0;
open IN, "dft.split" or die "$!\n";
if( <IN> =~ m/t/i )
{
  $split_dft = 1;
}
close IN;



open NKPT, "bse.nkpt" or die "Failed to open nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse nkpt\n";
@nkpt = ($1, $2, $3);
close NKPT;

$rundir = sprintf("../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

my @BSECommonFiles = ( "qinunitsofbvectors.ipt", "bvecs", "dft", "nelectron", "avecsinbohr.ipt", 
                       "xmesh.ipt", "nspin", "dft.split", "prefix" );
my @rundirFiles = ( "kmesh.ipt", "brange.ipt", "umklapp" );
#Checks for BSE prep
my $runBSE = 1;
if( -e "BSE/done" && -e "${rundir}/old" )
{
  $runBSE = 0;
  foreach (@BSECommonFiles)
  {
    if( compare( "$_", "BSE/$_") != 0 )
    {
      $runBSE = 1;
      print "Difference found in $_\n";
      last;
    }
  }
  if( $runBSE == 0 )
  {
    foreach( @rundirFiles )
    {
      if( compare( "${rundir}/$_", "BSE/$_") != 0 )
      {
        $runBSE = 1;
        print "Difference found in $_\n";
        last;
      }
    }
  }
}

if( $runBSE != 0 )
{
  `rm -r BSE` if (-e "BSE");
  mkdir "BSE";
  chdir "BSE";

  open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
  print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
  close NKPT;

  foreach( @rundirFiles )
  {
    copy( "../${rundir}/$_", $_ );
  }
  foreach( @BSECommonFiles )
  {
    copy ( "../$_", $_ );
  }
#  foreach ("kmesh.ipt", "brange.ipt") {
#    system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
#  }
#  `cp ../qinunitsofbvectors.ipt .`;
#  `cp ../bvecs .`;
#  `cp ../dft .`;
#  `cp ../nelectron .`;
#  `cp ../avecsinbohr.ipt .`;
#  `cp ../xmesh.ipt .`;
#  `cp ../nspin .`;
#  `cp ../${rundir}/umklapp .`;
#  `cp ../dft.split .`;
#  `cp ../prefix .`;
##  my $Nfiles = `cat Nfiles`;

  my $prefix;
  open PREFIX, "prefix" or die "$!\n";
  $prefix = <PREFIX>;
  chomp($prefix);
  close (PREFIX);



#  `cp -r ../${rundir}/Out .`;
  if( -l "Out" )  # Out is an existing link
  {
    unlink "Out" or die "Problem cleaning old 'Out' link\n$!";
  } 
  elsif(  -d "Out" ) #or Out is existing directory
  {
    rmtree( "Out" );
  } 
  elsif( -e "Out" ) #or Out is some other file
  {
    unlink "Out";
  } 
  print "../$rundir/Out\n";
  symlink ("../$rundir/Out", "Out") == 1 or die "Failed to link Out\n$!";

  if( $split_dft ) 
  {
    print "$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml Out/${prefix}_shift.save/data-file.xml\n";
    system("$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml Out/${prefix}_shift.save/data-file.xml") == 0
      or die "Failed to run qe_data_file.pl\n$!";
  }
  else
  {
    print "$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml\n";
    system("$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml") == 0
      or die "Failed to run qe_data_file.pl\n$!";
  }

  system("$ENV{'OCEAN_BIN'}/wfconvert.x system") == 0 
    or die "Failed to run wfconvert.x\n";

  system("$ENV{'OCEAN_BIN'}/ofermi.pl") == 0
    or die "Failed to run ofermi.pl\n";

  `cp eshift.ipt ../`;
  system("cp ../efermiinrydberg.ipt ./") == 0 
    or die "Failed to copy efermiinrydberg.ipt\n";

  print "Running setup\n";
  system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0
    or die "Failed to run setup\n";

  print "conugtoux\n";
  system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
  print "orthog\n";
  system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;

  `touch done`;
  chdir "../";
}


else {
  `touch BSE/old`;
  print "Nothing needed for bse wfns\n";
}

unless ($stat && $oldden) {
system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0
  or die "Failed to run rhoofg.x\n";
`wc -l rhoG2 > rhoofg`;
`sort -n -k 6 rhoG2 >> rhoofg`;
}

`touch done`;
exit 0;

