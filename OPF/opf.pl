#!/usr/bin/perl
# Copyright (C) 2014 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Copy;
use Cwd;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/opf\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}

my $dir = getcwd;
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################


my @CommonFiles = ("znucl", "opf.hfkgrid", "opf.fill", "opf.opts", "pplist", "screen.shells", "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution", "epsilon", "k0.ipt", "ibase", "scfac" );


my @DenDipFiles = ("rhoofg", "bvecs", "efermiinrydberg.ipt");
my @DenDipFiles2 = ( "masterwfile", "listwfile", "enkfile", "kmesh.ipt", "brange.ipt" );

my @ExtraFiles = ("specpnt", "Pquadrature" );

my $runPAW = 1;
if (-e "done" ) {
  $runPAW = 0;
  foreach (@CommonFiles) {
    if (`diff -q $_ ../Common/$_`) {
      $runPAW = 1;
      last;
    }
  }
}

if ($runPAW == 0 ) {
  print "Nothing new needed for PAW stage\n";
  exit 0;
}
`rm -f done`;


foreach (@CommonFiles) {
  `cp ../Common/$_ .` == 0 or die "Failed to get $_ from Common/\n";
}

foreach (@ExtraFiles) {
  `cp $ENV{'OCEAN_BIN'}/$_ .` == 0 or die;
}

open PPLIST, "pplist" or die "Failed to open pplist\n";
while (<PPLIST>) {
  chomp;
  s/\s+//;
  `cp "../$_" .`;
}

###################################

# Setup
###################################
print "Running PAW Setup\n";
system("$ENV{'OCEAN_BIN'}/pawsetup.x") == 0 or die "Failed to run pawsetup.x\n";

#`mkdir -p zdiag/ zpawinfo/`;
unless( -d "zpawinfo" )
{
  mkdir "zpawinfo" or die "$!";
}
###################################


# PP mod
###################################
print "Converting psp format\n";
open SITELIST, "sitelist" or die "Failed to open sitelist\n";
<SITELIST> =~ m/\s*(\d+)/ or die "Malformaed sitelst\n";
my $numsites = $1;
my $ppfilename;
my $elname;
my $elnum;
my $ppmodname;
for (my $i=0; $i < $numsites; $i++ ) {
  <SITELIST> =~ m/\s*(\S+)\s+(\d+)\s+(\S+)/ or die "Malformaed sitelst\n";
  $elname = $1;
  $elnum = $2;
  $ppfilename = $3;
  $ppmodname = $ppfilename . ".mod";
  system("echo '$ppfilename\n$ppmodname' | $ENV{'OCEAN_BIN'}/fhi2eric.x") == 0 
      or die "Failed to convert psp file $ppfilename\n";
}
close SITELIST;
###################################

# Load up pspoptions
###################################
open PSPO, "pspoptions" or die "Failed to open pspoptions\n";
my %pspopts;
my %pspfill;
while (<PSPO>) {
  $_ =~ m/\s*(\d+)\s+(\w+)\s+\S+\s+(\S+)\s+(\S+)/ or die "Malformed pspoptions\n";
  $pspopts{"$1"} = $3;
  $pspfill{"$1"} = $4;
  `cp "../$3" .`;
  `cp "../$4" . `;
}
close PSPO;


# shells
##################################
open SHELLS, "screen.shells" or die "Failed to open screen.shells\n";
my $numshells = 0;
my $allshells = '';
while (<SHELLS>) {
  chomp;
  $allshells .= $_ ." ";
}
close SHELLS;
my @rads = split(/\s+/, $allshells);
$numshells = $#rads + 1;
open SHELLS, ">shells" or die "Failed to open shells for writing\n";
print SHELLS "$numshells\n";
print SHELLS "$allshells\n";
close SHELLS;

# HFK
###################################
print "Starting HFK section\n";

open GRID, "opf.hfkgrid" or die;
my $grid = <GRID>;
chomp $grid;
close GRID;

my $optionfilename;

open HFINLIST, "hfinlist" or die "Failed top open hfinlist\n";
my $hfinline;
my $znucl;
my $nnum;
my $lnum;

while ( $hfinline = <HFINLIST> ) {
  ($hfinline =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)/) or die "Malformed hfinlist\n";
  $ppfilename = $1;
  $znucl = $2;
  $nnum = $3;
  $lnum = $4;
  $elname = $5;
  $elnum = $6;  

  $optionfilename = $pspopts{"$znucl"};
#  print "$optionfilename\n";
  die "No option file for $ppfilename\n" unless (-e $optionfilename );
  `cp $optionfilename atomoptions`;
  `cp ${ppfilename}.mod ppot`;
  open HFIN, ">HFIN" or die "Failed to open HFIN for writing\n";
  print HFIN "initgrid\n";
  print HFIN "$znucl $grid\n";
  print HFIN "ppload\nmkcorcon\nscreencore\ncalcso\nspartanfip\n";
#  print HFIN "$nnum  $lnum\n";
  print HFIN `cat "$pspfill{"$znucl"}"`;
#  print HFIN "calcso\n";
  print HFIN "quit\n";
  close HFIN;

  open HFIN, ">hfin1" or die;
  print HFIN "initgrid\n";
  print HFIN "$znucl $grid\n";
  print HFIN "ppload\nmkcorcon\nscreencore\nquit\n";
  close HFIN;

  open HFIN, ">hfin2" or die;
  print HFIN "initgrid\n";
  print HFIN "$znucl $grid\n";
  print HFIN "ppload\nmkcorcon\ncalcso\nquit\n";
  close HFIN;

  open HFIN, ">hfin3" or die;
  print HFIN "initgrid\n";
  print HFIN "$znucl $grid\n";
  print HFIN "ppload\nmkcorcon\nspartanfip\n";
  print HFIN `cat "$pspfill{"$znucl"}"`;
  print HFIN "quit\n";
  close HFIN;


  print "Running hfk.x\n";


  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin1 > hfk.${znucl}_${nnum}_${lnum}.1.log") == 0 or die;
  # Check the end of the log to see if we are ok
  my $hfk_status = `tail -n 1 hfk.${znucl}_${nnum}_${lnum}.1.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly.\nExiting ...\n";
  }

  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin2 > hfk.${znucl}_${nnum}_${lnum}.2.log") == 0 or die;
  $hfk_status = `tail -n 1 hfk.${znucl}_${nnum}_${lnum}.2.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly.\nExiting ...\n";
  }
  my $corezfile = sprintf("corezetaz%03i",$znucl);
  move("xifile","$corezfile");

  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin3 > hfk.${znucl}_${nnum}_${lnum}.3.log") == 0 or die;
  $hfk_status = `tail -n 1 hfk.${znucl}_${nnum}_${lnum}.3.log`;
  print "Done running hfk.x\n";


# Clean-up time
#  `mv xifile $corezfile`;
  unless( -d "zdiag_${znucl}_${nnum}_${lnum}" )
  {
    ( mkdir "zdiag_${znucl}_${nnum}_${lnum}" ) or die "$!";
  }

# Wander through the directory making some changes
  opendir DIR, $dir;
  while( my $file = readdir(DIR) )
  {
    if( ( $file =~ m/^[fg]k/ ) or ( $file =~ m/^(ae|ft|ps).z/ ) or
        ( $file =~ m/^(deflin|melfile|corez)/ ) or ( $file =~ m/^(rad|prj)filez/ ) or
#        ( $file =~ m/^(vcallel|vvallel|vc_bare|vcxxxxx|vvpseud)/ ) or 
        ( $file =~ m/^(vcallel|vvallel|vc_bare|vpseud|valence)/ ) or 
        ( $file =~ m/^coreorb/ ) or ($file =~ m/^phr/ ) )
    {
      move( $file, "zpawinfo" );
    }
    elsif( ( $file =~ m/^melfilez\w$/ ) or ( $file =~ m/^(sm|am|di|pr|psft|aeft)\w$/ ) or
           ( $file =~ m/^(mt|dif)\w\w$/ ) or ( $file =~ m/^(map|ex)/ ) )
    {
      move( $file, "zdiag_${znucl}_${nnum}_${lnum}" );
    }
  }

  my @paw_files_to_kill = ( 'radpot', 'hapot', 'hfpot', 'config', 'corcon', 'valcon', 'chgsum',
                            'atmf99', 'skip', 'psld', 'rslt', 'tmp', 'aprog', 'vxcofr' );
  foreach (@paw_files_to_kill)
  {
    unlink $_;
  }



#  `mv {f,g}k* shellR* {ae,ft,ps}?z* {deflin,melfile}z*n*l* {rad,prj}filez* corez* zpawinfo/`;
#  # vcsz*
#  `mkdir -p zdiag_${znucl}_${nnum}_${lnum}`;
#  `mv melfilez??? mt?? {sm,am,di,pr}? dif?? map* ex* {ps,ae}ft? zdiag_${znucl}_${nnum}_${lnum}/`;
#  `mv vcallel* vvallel* vc_bare* vcxxxxx* vvpseud* zpawinfo/`;
#  `rm -f radpot hapot hfpot config corcon valcon chgsum atmf99 skip psld rslt tmp aprog vxcofr`;
  last;
}
close HFINLIST;


######################################
print "PAW section done\n";
exit 0;
