#!/usr/bin/perl
# Copyright (C) 2014, 2016-2017 OCEAN collaboration
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


my @CommonFiles = ("znucl", "opf.hfkgrid", "opf.fill", "opf.opts", "pplist", "screen.shells", 
                   "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution", 
                   "scfac", "calc" );



my $runOBF = 1;
if (-e "done" ) {
  $runOBF = 0;
  foreach (@CommonFiles) {
    if (`diff -q $_ ../Common/$_`) {
      $runOBF = 1;
      last;
    }
  }
}

if ($runOBF == 0 ) {
  print "Nothing new needed for OBF stage\n";
  exit 0;
}

unlink "done";


foreach (@CommonFiles) {
  copy( "../Common/$_", "$_") == 1 or die "Failed to get $_ from Common/\n";
}


if( open CALC, "calc" )
{
  if( <CALC> =~ m/VAL/i )
  {
    print "No OPF calc for valence run\n";
    close CALC;
    exit 0;
  }
  close CALC;
}

###################################
# Quickcheck
open IN, "opf.opts" or die "Failed to open opf.opts: $!\n";
if( <IN> =~ m/^!/ )
{
  close IN;
  die "The input parameter opf.opts was not specified!\n";
}
close IN;

open IN, "opf.fill" or die "Failed to open opf.fill: $!\n";
if( <IN> =~ m/^!/ )
{
  close IN;
  die "The input parameter opf.fill was not specified!\n";
}
close IN;


# Setup
###################################
print "Running OBF Setup\n";
system("$ENV{'OCEAN_BIN'}/pawsetup.x") == 0 or die "Failed to run pawsetup.x\n";

unless( -d "zpawinfo" )
{
  mkdir "zpawinfo" or die "$!";
}
###################################


my $ppfilename;
my $ppmodname;

# Load up pspoptions
###################################
open PSPO, "pspoptions" or die "Failed to open pspoptions\n";
my %pspopts;
my %pspfill;
my %psplist;
while (<PSPO>) {
  $_ =~ m/\s*(\d+)\s+(\w+)\s+(\S+)\s+(\S+)\s+(\S+)/ or die "Malformed pspoptions\n";
  $psplist{"$1"} = $3;
  $pspopts{"$1"} = $4;
  $pspfill{"$1"} = $5;
  copy( "../$3", "$3" ) == 1 or die "Failed to copy $3\n$1";
  copy( "../$4", "$4" ) == 1 or die "Failed to copy $4\n$!";
  copy( "../$5", "$5" ) == 1 or die "Failed to copy $5\n$!";
}
close PSPO;

# PP mod
###################################
foreach  my $znucl (keys %psplist )
{
  $ppfilename = $psplist{$znucl};
  $ppmodname = $ppfilename . ".mod";
  if( -e "../$ppmodname" )
  {
    print "Using found modified psps\n";
    copy( "../$ppmodname", $ppmodname ) == 1 or die "Failed to copy $ppmodname\n$!";
  }
  else
  {
    print "Converting $ppfilename\n";
    system("echo '$ppfilename\n$ppmodname' | $ENV{'OCEAN_BIN'}/fhi2eric.x") == 0
        or die "Failed to convert psp file $ppfilename\n";
  }
}

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

foreach my $znucl (keys %psplist )
{
  $ppfilename = $psplist{$znucl};
  $optionfilename = $pspopts{"$znucl"};
#  print "$optionfilename\n";
  die "No option file for $ppfilename\n" unless (-e $optionfilename );
  copy( $optionfilename, "atomoptions" ) == 1 or die "Failed to copy $optionfilename to atomoptions\n$!";
  copy( "${ppfilename}.mod", "ppot" ) == 1 or die "Failed to copy ${ppfilename}.mod to ppot\n$!";

  system( "$ENV{'OCEAN_BIN'}/validate_opts.pl ${ppfilename} $optionfilename" ) == 0 
    or die "Failed to validate options file\nCheck $optionfilename\n";

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


  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin1 > hfk.${znucl}.1.log") == 0 or die;
  # Check the end of the log to see if we are ok
  my $hfk_status = `tail -n 1 hfk.${znucl}.1.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly for hfin1.\nExiting ...\n";
  }

  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin2 > hfk.${znucl}.2.log") == 0 or die;
  $hfk_status = `tail -n 1 hfk.${znucl}.2.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly for hfin2.\nExiting ...\n";
  }
  my $corezfile = sprintf("corezetaz%03i",$znucl);
  move("xifile","$corezfile");

  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin3 > hfk.${znucl}.3.log") == 0 or die;
  $hfk_status = `tail -n 1 hfk.${znucl}.3.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly for hfin3.\nExiting ...\n";
  }
  print "Done running hfk.x\n";


  unless( -d "zdiag_${znucl}" )
  {
    ( mkdir "zdiag_${znucl}" ) or die "$!";
  }

# Wander through the directory making some changes
  opendir DIR, $dir;
  while( my $file = readdir(DIR) )
  {
    if( ( $file =~ m/^[fg]k/ ) or ( $file =~ m/^(ae|ft|ps).z/ ) or
        ( $file =~ m/^(deflin|melfile|corez)/ ) or ( $file =~ m/^(rad|prj)filez/ ) or
        ( $file =~ m/^(vcallel|vvallel|vc_bare|vpseud|valence)/ ) or 
        ( $file =~ m/^coreorb/ ) )
    {
      move( $file, "zpawinfo" );
    }
    elsif( $file =~ m/^phr/ )
    {
      my $dest = sprintf( "${file}z%03i", $znucl );
      move( $file, "zpawinfo/$dest" );
    }
    elsif( ( $file =~ m/^melfilez\w$/ ) or ( $file =~ m/^(sm|am|di|pr|psft|aeft)\w$/ ) or
           ( $file =~ m/^(mt|dif)\w\w$/ ) or ( $file =~ m/^(map|ex)/ ) or ( $file =~ /hfin/ ) or
           ( $file =~ m/hfk.+log/ ) or ( $file =~ m/aetotal/ ) or ( $file =~ m/radf/ ) )
    {
      move( $file, "zdiag_${znucl}" );
    }
    elsif( ( $file =~ m/^angs\d$/ ) or ( $file =~ m/^ldep\d$/ ) or ( $file =~ m/^shellR\d\.\d\d$/ ) )
    {
      unlink $file;
    }
  }

  my @paw_files_to_kill = ( 'radpot', 'hapot', 'hfpot', 'config', 'corcon', 'valcon', 'chgsum',
                            'atmf99', 'skip', 'psld', 'rslt', 'tmp', 'aprog', 'vxcofr', 'xifile',
                            'actual', 'ppot' );
  foreach (@paw_files_to_kill)
  {
    unlink $_;
  }


}


######################################
print "OBF section done\n";

open DONE, ">done" or exit 0;
print DONE "1\n";
close DONE;

exit 0;
