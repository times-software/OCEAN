#!/usr/bin/perl
# Copyright (C) 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
use strict;

if( scalar @ARGV < 2 ) 
{
  die "Usage validate_fill psp fill\n";
}

my @shellsize;
$shellsize[0] = 2;
$shellsize[1] = 6;
$shellsize[2] = 10;
$shellsize[3] = 14;



open PSP, $ARGV[0] or die "Failed to open $ARGV[0]\n$!";
open FILL, $ARGV[1] or die "Failed to open $ARGV[1]\n$!";

# First line in fhi is comments
<PSP>;
#second line has zatom and zion which can be floats? but must be at least ints
<PSP> =~ m/^\s+(\d+(\.\d*)?)\s+(\d+(\.\d*)?)/ or die "Failed to parse $ARGV[0]\n";
my $z = $1;
my $z_eff = $3;

close PSP;


my $line = <FILL>;
#<FILL> =~ m/0*(\d+)/ or die "Failed to parse Z from $ARGV[1]\n";
$line =~ m/0*(\d+)/ or die "Failed to parse Z from $ARGV[1]\n";
my $z_fill = $1;
if( abs( $z - $z_fill ) > 0.01 ) 
{
  die "Mismatch between psp and fill: $z  $z_fill\n";
}

my @occs;
<FILL> =~ m/(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse occs from $ARGV[1]\n";
$occs[0] = $1;
$occs[1] = $2;
$occs[2] = $3;
$occs[3] = $4;

my $core_occ = 0;
for( my $i = 0; $i < 4; $i++ )
{
  $core_occ += $occs[$i] * $shellsize[$i];
}
if( abs( $z - $core_occ - $z_eff )  > 0.01 ) 
{
  die "Mismatch between psp and fill for Z effective:\n\t$z  $core_occ  $z_eff\n";
}

print "Z = $z.  Zeff = $z_eff. Core occ. = $core_occ \n";

exit 0;
