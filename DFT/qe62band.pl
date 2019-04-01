#!/usr/bin/perl
# Copyright (C) 2015, 2017 - 2019 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
# qeband.pl looks at the SCF run and pulls up the occupation numbers for 
# all of the k-pts. The tolerance parameter below sets what it considers
# to be an occupied (for 1 - tolerance ) unoccupied band. This means 
# that qeband should get a good reading for min and max ranges for metallic 
# systems as well.
#
# QE now puts all the eigenvalues in xml. This script pulls them out and drops 
# them into enkfile
use strict;

my $tolerance = 0.000001;
open IN, "../work_dir" or die "Failed to open work_dir\n$!";
my $work_dir = <IN>;
close IN;
chomp $work_dir;
$work_dir =~ s/\'//g;
$work_dir =~ s/^\.//;
$work_dir =~ s/^\///;

open IN, "../prefix" or die "Failed to open prefix\n$!";
my $prefix = <IN>;
close IN;
chomp $prefix;

# Different behavior for metals/nonmetals
open IN, "../metal" or die "Failed to open metal\n$!";
my $metal_line = <IN>;
close IN;
my $metal = 0;
if( $metal_line =~ m/true/i )
{
  $metal = 1;
}

# Spin=2 needs to use the metals version
if( $metal == 0 )
{
  open IN, "../nspin" or die "Failed to open nspin\n$!\n";
  my $spin = <IN>;
  close IN;
  $metal = 2 if( $spin =~ m/2/ );
}


my $band_max = 0;
my $band_min = -1;


my $datafile = $work_dir . '/' . $prefix . ".save/data-file-schema.xml";
open IN, $datafile or die "Failed to open $datafile\n$!";

my $num_e;
while( 1 )
{
  die "Parsing $datafile failed\n" unless( my $line = <IN> );
# <nelec>8.000000000000000e0</nelec>
  if( $line =~ m/\<nelec\>(\d+\.\d+[Ee][+-]?\d+)/ )
  {
    $num_e = $1;
    last;
  }
}


if( $metal == 0 ) 
{
  $band_max = $num_e / 2;
#  print "Number of electrons $num_e\n";
  if( $band_max * 2 - $num_e > 0.01 )
  {
    die "Odd number of electrons, but metal = false\n$band_max\t$num_e\n";
  }
  $band_min = $band_max + 1;
#  close IN;
#  print $band_max . "\t" . $band_min . "\n";
}

my @energies;

my $nkpt = 0;

while( my $line = <IN> )
{
  if( $metal != 0 && $line =~ m/<occupations size=\"\d+\">/ ) #([\s\d\.eE+-]+)<\/occupations>/ )
  {
    my $occs = '';
    $line = <IN>;

    until( $line =~ m/occupations/ )
    {
      chomp $line;
      $occs .= $line . ' ';
      $line = <IN>;
    }
    my @occs = split( ' ', $occs );

    my $count = 0;
    my $min_count = 0;
#    foreach my $occ (@occs )
#    {
#      $count++ if( $occ > $tolerance );
#      $min_count++ if( $occ > 1-$tolerance );
#    }
    for( my $i = 0; $i < scalar @occs; $i++ )
    {
      $count = $i + 1 if( $occs[$i] > $tolerance );
    }
    for( my $i = scalar @occs - 1; $i >= 0; $i-- )
    {
      $min_count = $i + 1 if( $occs[$i] < $tolerance );
#      if( $nkpt == 0 ) { print "$occs[$i]\t$min_count\t$i\n"; }
    }
    $nkpt ++;
    $band_max = $count if( $count > $band_max );
    if( $band_min > 0 )
    {
      $band_min = $min_count if( $min_count < $band_min );
    }
    else
    {
      $band_min = $min_count;
    }
  }

  if( $line =~ m/<eigenvalues size=\"\d+\">/ )
  {
    my $eigs = '';
    $line = <IN>;
    until( $line =~ m/eigenvalues/ )
    {
      chomp $line;
      $eigs .= $line . ' ';
      $line = <IN>;
    }

    my @eigs = split( /\s+/, $eigs );
    push @energies, \@eigs;

  }
}



close IN;
print "Found $nkpt k-points\n";
print $band_max . "\t" . $band_min . "\n";

open OUT, ">brange.stub" or die "Failed to open brange.stub for writing\n$!";
print OUT "1    $band_max\n$band_min    ";
close OUT;


open ENK, ">", "enkfile" or die "Failed to open enkfile for writing\n";
for( my $k = 0; $k < $nkpt; $k++ )
{
  my $n = 3;
  my $start = 0;
  my $stop = $band_max ;
  my $delim = " ";
  my @eslice = @{ $energies[$k] }[ $start .. $stop ];
#    while (my @x = splice @{ $energies[$k] }, 1, $n) {
  while (my @x = splice @eslice, 1, $n) 
  {
     print ENK join($delim, @x), "\n";
  }   
#  print "\n";
  
  
  my $start = $band_min - 1;
  my $stop = scalar @{ $energies[$k] };
  @eslice = @{ $energies[$k] }[ $start .. $stop ];
  while (my @x = splice @eslice, 1, $n) 
  {
    print ENK join($delim, @x), "\n";
  }   
#  print "\n";
}
close ENK;

open ENK, ">", "enk_un" or die "Failed top open enk_un\n$!";
for( my $k = 0; $k < $nkpt; $k++ )
{
  my $n = 3;
  my $delim = " ";
#  my $start = 0;
#  my $stop = scalar @{ $energies[$k] };
#  my @eslice = @{ $energies[$k] }[ $start .. $stop ];

  while (my @x = splice @{ $energies[$k] }, 1, $n)
  {
     print ENK join($delim, @x), "\n";
  }
} 
close ENK;



exit 0;

