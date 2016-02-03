#!/usr/bin/perl
# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
#
# qeband.pl looks at the SCF run and pulls up the occupation numbers for 
# all of the k-pts. The tolerance parameter below sets what it considers
# to be an occupied (for 1 - tolerance ) unoccupied band. This means 
# that qeband should get a good reading for min and max ranges for metallic 
# systems as well.

# Output: creates the file brange.stub to which the total number of bands
#         can be appended to get brange.ipt

use strict;

my $tolerance = 0.000001;

open IN, "work_dir" or die "Failed to open work_dir\n$!";
my $work_dir = <IN>;
close IN;
chomp $work_dir;
$work_dir =~ s/\'//g;
$work_dir =~ s/^\.//;
$work_dir =~ s/^\///;

open IN, "prefix" or die "Failed to open prefix\n$!";
my $prefix = <IN>;
close IN;
chomp $prefix;

my $datafile = $work_dir . '/' . $prefix . ".save/data-file.xml";
open IN, $datafile or die "Failed to open $datafile\n$!";

while( 1 )
{
  die "Parsing $datafile failed\n" unless( my $line = <IN> );
  last if( $line =~ m/<EIGENVALUES>/ );
}

my @filelist;
while( 1 )
{
  die "Parsing $datafile failed\n" unless( my $line = <IN> );
  if( $line =~ m/<DATAFILE\.?\d? iotk_link=\"(\S+)\">/ )
  {
    push @filelist, $1;
  }
  last if( $line =~ m/<\/EIGENVALUES>/ );
}
close IN;

print "Counted " . scalar @filelist . " kpts\n";

my $band_max = 0;
my $band_min = -1;
foreach my $file_stub (@filelist)
{
  $file_stub =~ s/^\.//;
  $file_stub =~ s/^\///;
  my $file = $work_dir . '/' . $prefix . ".save/"  . $file_stub . "\n";
#  print $file;

  my @occ;
  open IN, $file or die "Failed to open $file\n$!";
  while( 1 )
  {
    die "Parsing $file failed\n" unless( my $line = <IN> );
    last if( $line =~ m/<OCC/ );
  }
  while( 1 )
  {
    die "Parsing $file failed\n" unless( my $line = <IN> );
    if( $line =~ m/(\d+\.\d+[Ee]?[+-]?\d*)/ )
    {
      push @occ, $1;
    }
    last if( $line =~ m/<\/OCC/ );
  }
  close IN;
  my $count = 0;
  my $min_count = 0;
  foreach my $i (@occ)
  {
    $count++ if( $i > $tolerance );
    $min_count++ if( $i > 1-$tolerance );
  }
#  print $file_stub . "\t" . $count . "\n";
  $band_max = $count if( $count > $band_max );
  if( $band_min > 0 )
  {
    $band_min = $min_count if( $min_count < $band_min );
  }
  else
  {
    $band_min = $min_count;
  }
#  print "$file_stub:\t$count\t $min_count\n";
}
$band_min++;
# band_min is the lowest count whereas band_max is the highest count
#   Therefore brange will be 
#   1         band_max
#   band_min  total bands
print $band_max . "\t" . $band_min . "\n";


open OUT, ">brange.stub" or die "Failed to open brange.stub for writing\n$!";
print OUT "1    $band_max\n$band_min    ";
close OUT;
