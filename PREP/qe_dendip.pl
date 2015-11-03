#!/usr/bin/perl
# Copyright (C) 2014 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Path qw( rmtree );

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/qe_dendip\.pl/;
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


my @QEFiles     = ( "rhoofr" );
my @CommonFiles = ( "paw.nkpt", "nkpt", "qinunitsofbvectors.ipt", "avecsinbohr.ipt" );

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

my $rundir;
## process paw wf files ##

open NKPT, "paw.nkpt" or die "Failed to open paw.nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse paw.nkpt\n";
my @nkpt = ($1, $2, $3);
close NKPT;
$rundir = sprintf("../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

unless( -e "PAW/done" && -e "${rundir}/old" ) {
`rm -r PAW` if (-e "PAW");
mkdir "PAW"; 
chdir "PAW";

open NKPTS, ">nkpts" or die "Failed to open nkpts for writing\n";
print NKPTS $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
close NKPTS;


foreach ("Nfiles", "kmesh.ipt", "brange.ipt", "qinunitsofbvectors.ipt" ) {
  system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
}
`cp ../qinunitsofbvectors.ipt .`;
`cp ../bvecs .`;
open BRANGE, "brange.ipt";
my @brange;
<BRANGE> =~ m/(\d+)\s+(\d+)/;
$brange[0] = $1;
$brange[1] = $2;
<BRANGE> =~ m/(\d+)\s+(\d+)/;
$brange[2] = $1;
$brange[3] = $2;
close BRANGE;
my $nelectron = `cat ../nelectron`;
open BRANGE, ">brange.ipt";
print BRANGE "1  " . $nelectron/2 . "\n";
print BRANGE $nelectron/2+1 . "    $brange[3]\n";
close BRANGE;

my $prefix;
open PREFIX, "../../Common/prefix";
while (<PREFIX>) {
  chomp;
  $prefix = $_;
}
close (PREFIX);


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

# make prep.in

#open  PREP, ">prep.in";
#print PREP  "&input\n";
#print PREP  "  prefix = '";
#print PREP  ${prefix};
#print PREP  "'\n";
#print PREP  "  work_dir = './Out'\n";
#print PREP  "/\n";
#close PREP;


#system("$ENV{'OCEAN_BIN'}/qe_wfconvert.x ${prefix}") == 0 
system("$ENV{'OCEAN_BIN'}/qe_wfconvert.x") == 0 
  or die "Failed to run wfconvert.x\n$!";


`touch done`;
chdir "../";
}
else {
  `touch PAW/old`;
  print  "Nothing needed for PAW wfns\n";
}


print "Done with PAW files";


## process bse wf files ##

open NKPT, "bse.nkpt" or die "Failed to open nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse nkpt\n";
@nkpt = ($1, $2, $3);
close NKPT;

$rundir = sprintf("../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

unless( -e "BSE/done" && -e "${rundir}/old" ) {

  `rm -r BSE` if (-e "BSE");
  mkdir "BSE";
  chdir "BSE";

  open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
  print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
  close NKPT;

  foreach ("Nfiles", "kmesh.ipt", "brange.ipt") {
    system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
  }
  `cp ../qinunitsofbvectors.ipt .`;
  `cp ../bvecs .`;
  `cp ../nelectron .`;
#  `cp ../${rundir}/umklapp .`;
  my $Nfiles = `cat Nfiles`;

  open BRANGE, "brange.ipt";
  my @brange;
  <BRANGE> =~ m/(\d+)\s+(\d+)/;
  $brange[0] = $1;
  $brange[1] = $2;
  <BRANGE> =~ m/(\d+)\s+(\d+)/;
  $brange[2] = $1;
  $brange[3] = $2;
  close BRANGE;
  my $nelectron = `cat ../nelectron`;
  open BRANGE, ">brange.ipt";
  print BRANGE "1  " . $nelectron/2 . "\n";
  print BRANGE $nelectron/2+1 . "    $brange[3]\n";
  close BRANGE;

  my $prefix;
  open PREFIX, "../../Common/prefix";
  while (<PREFIX>) {
    chomp;
    $prefix = $_;
  }
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


  system("$ENV{'OCEAN_BIN'}/qeband.x") == 0
    or die "Failed to run qeband.x\n";

  system("$ENV{'OCEAN_BIN'}/qe_wfconvert.x system") == 0 
    or die "Failed to run wfconvert.x\n";

#KG#  system("$ENV{'OCEAN_BIN'}/ofermi.pl") == 0
#KG#    or die "Failed to run ofermi.pl\n";
  `cp eshift.ipt ../`;
  system("cp efermiinrydberg.ipt ../") == 0 
    or die "Failed to copy efermiinrydberg.ipt\n";

  `cp ../avecsinbohr.ipt .`;
  `cp ../../Common/xmesh.ipt .`;
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

