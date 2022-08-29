#!/usr/bin/perl
# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#


######################
# Edges is used to expand the user shorthand for which edges into the full OCEAN-compatible listing
#
#
# In the input file 
# * nedges.ipt determines the number of edges listed in the input file
# * edges.ipt lists the edges, but these can be in shorthand
#
# The standard format is
#  AtomNumber  Principle  Angular
#
# The shorthand is 
#  -Z          Principle  Angular
#     where for we will calcuate for every Z in the input file


use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/edges\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }

if( -e "postDefaultsOceanDatafile" )
{
  print "New JSON detected. Not running edges.pl\n";
  exit 0;
}

open IN, "calc" or die "Failed to open calc: $!\n";
if( <IN> =~ m/val/i )
{
  close IN;
  print "No edges for a valence calculation\n";

  my $overwriteEdges = 0;
  if( -e "screen.mode" )
  {
    open MODE, "screen.mode" or die "Failed to open screen.mode\n$!";
    if( <MODE> =~ m/grid/i )
    {
      $overwriteEdges = 1;
    }
    close MODE;
  }
  
  if( $overwriteEdges == 1 || not -e "edges" )
  {
    print "Over-writing edges\n";
    open OUT, ">nedges" or die "Failed to open edges for writing\n$!";
    print OUT "1\n";
    close OUT;

    open OUT, ">edges" or die "Failed to open edges for writing\n$!";
    print OUT "0 0 0\n";
    close OUT;
  }
  exit 0;
}
close IN;


# fill $line with edges.ipt
my $line;
open IN, "edges.ipt" or die "Failed to open edges.ipt for reading\n$!";
while (<IN>)
{
  $line .= $_;
}
close IN;

# get rid of line breaks and extra spaces
$line =~ s/[\r]/ /g;
$line =~ s/[\n]/ /g;
$line =~ s/\s+/ /g;

my @orig_edges = split( " ", $line );
  
if( scalar @orig_edges < 3 )
{
  print "The input edges was not long enough!\nCannot continue!!\n";
  exit 1;
}


# bring in all the atom info
open IN, "typat" or die "Failed to open typat\n$!";
$line = '';
while( <IN> )
{
  $line .= $_;
}
close IN;
$line =~ s/[\r]/ /g;
$line =~ s/[\n]/ /g;
$line =~ s/\s+/ /g;
my @typat = split( " ", $line );

open IN, "znucl" or die "Failed to open znucl\n$!";
$line = '';
while( <IN> ) 
{
  $line .= $_;
}
close IN;
$line =~ s/[\r]/ /g;
$line =~ s/[\n]/ /g;
$line =~ s/\s+/ /g;
my @znucl = split( " ", $line );

# Check that all of typat is 1) an integer and 2) within the range of the znucl lookup
my $z_count = scalar @znucl;
foreach my $i (@typat)
{
  $i =~ m/^(\d+)$/ or die "Unexpected string in typat: $_\n";
  if( $i < 1 || $i > $z_count )
  {
    die "Value in typat is too small or too large!\n\t$i\n";
  }
}

# Check that every znucl is 1) and integer and 2) an allowable element > 0 < 110?
my $j = 0;
foreach my $i (@znucl)
{
  $j++;
  $i =~ m/^(\d+)$/ or die "Unexpected string in znucl: $_\n";
  if( $i < 0 || $i > 110 ) 
  {
    die "Element $j in znucl is negative or too large!\t$i\n";
  }
}


# Build a hash with each unique Z
# Key = Z, Value = list of sites
my %SitesByZ;
for( $j = 0; $j <  scalar @typat; $j++ )
{
  # zero ordered arrays!
  my $i = $typat[$j] - 1;
  my $z = $znucl[ $i ];
  push @{ $SitesByZ{ $z } }, $j+1;
}

foreach my $key ( keys %SitesByZ )
{
  print $key . ":\n";
  foreach my $site ( @{ $SitesByZ{ $key } } )
  {
    print "  $site";
  }
  print "\n";
}


# Start in on the edges
# is edges divisible by 3?
my $orig_nedges = scalar @orig_edges;
unless( $orig_nedges % 3 == 0 )
{
  die "Malformed edges in input\n";
}

# Create a hash of Z
#   which stores hash NL
#     which has array of sites

# When writing out edges we want to group by this, and sort sites into numerical order
my %FullEdges;

for( my $i = 0; $i < $orig_nedges; $i+=3 )
{
  my $Z = $orig_edges[$i];
  my $N = $orig_edges[$i+1];
  my $L = $orig_edges[$i+2];

  if( $Z == 0 ) 
  {
    die "Illegal Z value in edges:\t$Z\n";
  }
  if( $N < 1 || $N > 7 )
  {
    die "Illegal principle qunatum number:\t$N\n";
  }
  if( $L < 0 || $L > 4 )
  {
    die "Illegal angular quantum number:\t$L\n";
  }

# Should be able to set value of hash of hash of hash
  # Are we selecting all Z?
  if( $Z < 0 )
  {
    my $ZEE = abs( $Z );
    unless( exists $SitesByZ{ $ZEE } )
    {
      die "Requested Z in edges that is not present in znucl/typat:\t$Z\t$ZEE\n";
    }
    foreach my $site ( @{ $SitesByZ{ $ZEE } } )
    {
      $FullEdges{ $ZEE }{ "$N $L" }{ $site } = 1;
    }
  }
  else
  {
    # Z is the site number
    if ( $Z > scalar @typat )
    {
      die "Site requested in edges is larger than the number of site:\t$Z\n";
    }
    my $ZEE = $typat[ $Z - 1 ];
    my $site = $Z;
    $FullEdges{ $ZEE }{ "$N $L" }{ $site } = 1;
  }
}


my @edges;
foreach my $ZEE ( keys %FullEdges )
{
  foreach my $NL ( keys %{ $FullEdges{ $ZEE } } )
  {
    my @tempsite;
    foreach my $site ( keys %{ $FullEdges{ $ZEE }{ $NL } } )
    {
      push @tempsite, $site;
    }

    my @sortsite = sort { $a <=> $b } @tempsite;
    foreach my $site ( @sortsite )
    {
      push @edges, "$site $NL";
    }
  }
}

open OUT, ">nedges" or die "Failed to open edges for writing\n$!";
print OUT scalar @edges . "\n";
close OUT;

open OUT, ">edges" or die "Failed to open edges for writing\n$!";
foreach my $i ( @edges )
{
  print $i . "\n";
  print OUT $i . "\n";
}
close OUT;
