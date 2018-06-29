#!/usr/bin/perl
# Copyright (C) 2018 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/structure\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }

open NTYPE, "ntype" or die "Failed to open ntype\n$!";
my $ntype = <NTYPE>;
close NTYPE;
chomp $ntype;

my @znucl;
open ZNUCL, "znucl" or die "Failed to open znucl\n$!";
while (<ZNUCL>)
{
  chomp;
  push @znucl, split / /;
}
close ZNUCL;

open NATOMS, "natoms" or die "Failed to open natoms\n$!";
my $natoms = <NATOMS>;
close NATOMS;
chomp $natoms;

my @typat;
open TYPAT, "typat" or die "Failed to open typat\n$!";
while (<TYPAT>)
{
  chomp;
  push @typat, split / /;
}
close TYPAT;

# assume no more than 110 elements
if( $ntype < 1 || $ntype > 110 )
{
  $ntype = scalar @znucl;
  open NTYPE, ">", "ntype" or die "Failed to open ntype\n$!";
  print NTYPE "$ntype\n";
  close NTYPE;
  print "Detected $ntype unique atom types\n";
}

# test for positive integer, allow some dumb options
unless( $natoms =~ m/\A\+?[0-9]*[1-9][0-9]*\z/ )
{
  $natoms = scalar @typat;
  open NATOMS, ">", "natoms", or die "Failed to open natoms\n$!";
  print NATOMS "$natoms\n";
  close NATOMS;
  print "Detected $natoms atoms\n";
}

# Later add some detection to make sure that the indices in typat are within the number of elements
# specfied by znucl

exit 0;
