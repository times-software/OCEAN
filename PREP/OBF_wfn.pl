#!/usr/bin/perl
# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/OBF_wfn\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################

#FAKE INPUTS FOR NOW

#my $para_prefix = "mpirun -n 16 ";
# These set the spline for the V_xc in the interpolation
my $ham_kpoints = "4 4 4 ";


# Step 1: Create support files
my @CommonFiles = ("nkpt", "k0.ipt", "qinunitsofbvectors.ipt", "nbands", "xmesh.ipt", "para_prefix", 
                   "pool_control", "ham_kpoints", "core_offset", "avecsinbohr.ipt", "tmp_dir" );
my @ExtraFiles = ("specpnt", "Pquadrature", "sphpts" );
my @DFTFiles = ("rhoofr", "efermiinrydberg.ipt", "nelectron");
my @PawFiles = ("hfinlist", "xyz.wyck");


foreach(@ExtraFiles)
{
  `cp $ENV{'OCEAN_BIN'}/$_ .` == 0 or die;
}

foreach (@DFTFiles)
{
  `cp ../DFT/$_ .` == 0 or die "Failed to get $_ from DFT/\n";
}

foreach (@CommonFiles) {
  `cp ../Common/$_ .` == 0 or die "Failed to get $_ from Common/\n";
}

foreach(@PawFiles)
{
  `cp ../OPF/$_ .` == 0 or die "Failed to get $_ from OPF/\n";
}


# TERRIBLE THINGS
`cp nkpt kmesh.ipt`;
`cp k0.ipt scaledkzero.ipt`;


### load up the para_prefix

my $pool_size = 1;
open INPUT, "pool_control" or die;
while (<INPUT>)
{
  if( $_ =~ m/interpolate kpt\s+(\d+)/ )
  {
    $pool_size = $1;
    last;
  }
}
close INPUT;

$ham_kpoints = `cat ham_kpoints`;
chomp($ham_kpoints);


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

my $tmpdir = "undefined";
if( open TMPDIR, "tmp_dir" )
{
  $tmpdir = <TMPDIR>;
  chomp( $tmpdir );
  close TMPDIR;
}

system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0 or die "$!\nBVECS.PL Failed\n";


`ln -sf ../DFT/Out .`;
#`cp -r ../DFT/Out .`;

`echo 1 > core`;
system("$ENV{'OCEAN_BIN'}/kgen_qe.x") == 0 or die "KGEN.X Failed\n";


# Prep input file
open BOFX, ">bofx.in" or die "$!\nFailed to open bofx.in for writing\n";
print BOFX "&input\n" .
           "  prefix = 'system_opt'\n" .
           "  outdir = './Out'\n";
unless( $tmpdir =~ m/undefined/ )
{
  print BOFX "  wfcdir = '$tmpdir'\n";
}
print BOFX "  updatepp = .false.\n" .
           "  ncpp = .true.\n" .
           "  calculation = 'ocean_bofx'\n" .
           "/\n" .
           " K_POINTS\n" .
           "$ham_kpoints 0 0 0\n";
close BOFX;

open OBF2, ">obf2loc.in" or die "$!\nFailed to open obf2loc.in for writing\n";
print OBF2 "&input\n" .
          "  prefix = 'system_opt'\n" .
          "  outdir = './Out'\n" .
          "  updatepp = .false.\n" .
          "  ncpp = .true.\n" .
          "  calculation = 'obf2local'\n" .
          "/\n" .
          " K_POINTS\n" .
          "$ham_kpoints 0 0 0\n";
close OBF2;


my $nbands = `cat nbands`;
chomp( $nbands );

open QIN, ">q.in" or die "$!\nFailed to open q.in for writing\n";
print QIN "&input\n"
        . "  prefix = 'system_opt'\n"
        . "  outdir = './Out/'\n"
        . "  band_subset = 1 $nbands\n"
        . "/\n"
        . "K_POINTS\n";
open KPTS, "kpts.0001" or die "Failed to open kpts.0001\n";
while (<KPTS>)
{
  print QIN $_;
}
close QIN;
close KPTS;

#`echo "&input" > q.in`;
#`echo "  prefix = 'system_opt'" >> q.in`;
#`echo "  outdir = './Out/'" >> q.in`;
#`echo -n "  band_subset = 1 " >> q.in`;
#`cat nbands >> q.in`;
#`echo "/" >> q.in`;
#`echo "K_POINTS" >> q.in`;
##`echo " crystal" >> q.in`;
##`wc -l kpts.0001 | awk '{print $1}' >> q.in`;
#`cat kpts.0001 >> q.in`;

my $cksout = "";
my $ckscount = 0;
open EDGE, "hfinlist" or die "Failed to open hfinlist\n$!";
while (<EDGE>) {
  $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
  my $ppname = $1;
  my $znum = $2;
  my $nnum = $3;
  my $lnum = $4;
  my $elname = $5;
  my $elnum = $6;

  my $zstring = sprintf("z%03i", $znum);
  print $zstring ."\n";
#  `ln -sf ../PAW/zpawinfo/*${zstring}* .`;
  `cp ../OPF/zpawinfo/prjfile$zstring . `;
  `cp ../OPF/zpawinfo/ft?$zstring .`;
  my $templine = `ls ../OPF/zpawinfo/*$zstring`;
  chomp($templine);
  my @pslist = split(/\s+/, $templine);
  foreach (@pslist)
  {
    $_ =~ m/ae(\S+)/;
    `cp ../OPF/zpawinfo/ae$1 .`;
    `ln -sf ae$1 ps$1`;
  }

  open ZNL, ">ZNL" or die;
  print ZNL "$znum  $nnum  $lnum\n";
  close ZNL;

  $cksout .= "$elname  $elnum  cbinf\n";
  $ckscount++;
}
  open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
#  print CKSIN "1\n$elname  $elnum  cbinf\n";
  print CKSIN "$ckscount\n";
  print CKSIN "$cksout";
  close CKSIN;

print "\nRunning QDIAG\n$para_prefix $ENV{'OCEAN_BIN'}/ocean_qdiagp.x $pool_size < q.in > q.out";
system("$para_prefix $ENV{'OCEAN_BIN'}/ocean_qdiagp.x $pool_size < q.in > q.out") == 0
  or die "Failed to run qdiag\n$!";

print "\nRunning BOFX";
system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_ham_o.x < bofx.in > bofx.out") == 0
  or die "Failed to run bofx\n$!";

print "\nSkipping OBF2LOC";
#print "\nRunning OBF2LOC";
#system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_ham_o.x < obf2loc.in > obf2loc.out") == 0
#  or die "Failed to run obf2loc\n$!";


print "\n";
`touch done`;


exit 0;
