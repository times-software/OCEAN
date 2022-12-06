#!/usr/bin/perl
# Copyright (C) 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#


use strict;
use File::Copy;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/eplot\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }


my @cnbseFiles = ( "nbuse.ipt", "brange.ipt", "kmesh.ipt", "qinunitsofbvectors.ipt", "k0.ipt", 
                   "xmesh.ipt", "avecsinbohr.ipt", "xyz.wyck", "bloch_selector" );


print "Welcome to EPLOT.pl\n";
print "Input path to CNBSE/RIXS directory (blank for local)\n";
my $path = <STDIN>;
chomp $path;

if( length $path > 0 )
{
  print "Looking for necessary files in $path ...\n";
  foreach (@cnbseFiles) {
    copy( "$path/$_", $_ ) or die "Failed to get $path/$_\n$!\n";
  }
}
else
{
  $path = ".";
}

print "What is the name of the exciton file?\n";
my $echamp_file = <STDIN>;
chomp $echamp_file;
unless( -e "$path/$echamp_file" )
{
  die "Failed to find $path/$echamp_file";
}

copy( "$path/$echamp_file", $echamp_file ) or die "Failed to copy $path/$echamp_file\n$!\n";

print "What is the name of the output (cube) file?\n(Default is $echamp_file.cube)\n";

my $cube_file = <STDIN>;
chomp $cube_file;
if( length( $cube_file ) == 0 )
{
  $cube_file = $echamp_file . ".cube";
}


print "What supercell size?\n";
<STDIN> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die;
my @rmesh;
$rmesh[0] = $1; $rmesh[1] = $2; $rmesh[2] = $3;

print "What offset for the supercell?\n";
<STDIN> =~ m/(-?\d+)\s+(-?\d+)\s+(-?\d+)/ or die;
my @rstart;
$rstart[0] = $1; $rstart[1] = $2; $rstart[2] = $3;


open OUT, ">exciton_plot.ipt" or die;
print OUT "$echamp_file\n$cube_file\n";
print OUT "$rmesh[0]  $rmesh[1]  $rmesh[2]\n";
print OUT "$rstart[0]  $rstart[1]  $rstart[2]\n";
close OUT;

open IN, "bloch_selector" or die;
my $b_sel = <IN>;
chomp $b_sel;
if( $b_sel =~ m/1/ )
{
  symlink "$path/u2par.dat", "u2par.dat";
}
else
{
  symlink "$path/u2.dat", "u2.dat";
}
system( "$ENV{'OCEAN_BIN'}/exciton_plot.x") == 0 or die "Failed to run exciton_plot.x correctly\n";

print "All done\n";
exit 0;


