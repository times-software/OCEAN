#!/usr/bin/perl
# Copyright (C) 2015, 2017 OCEAN collaboration
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

# If the run is an insulator (metal = .false.) then instead grab the number 
# of electrons

# Output: creates the file brange.stub to which the total number of bands
#         can be appended to get brange.ipt

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

my $datafile = $work_dir . '/' . $prefix . ".save/data-file.xml";
open IN, $datafile or die "Failed to open $datafile\n$!";

if( $metal == 0 )
{
  while( 1 )
  {
    die "Parsing $datafile failed\n" unless( my $line = <IN> );
    last if( $line =~ m/<NUMBER_OF_ELECTRONS[\s=\"\d\w]*>/ );
  }
  my $line = <IN>;
  $line =~ m/(\d+\.\d+E\+\d+)/ or die "Failed to parse number of electrons\n$line";
  my $num_e = $1;
  $band_max = $num_e / 2;
#  print "Number of electrons $num_e\n";
  if( $band_max * 2 - $num_e > 0.01 )
  {
    die "Odd number of electrons, but metal = false\n$band_max\t$num_e\n";
  }
  $band_min = $band_max + 1;
  close IN;
  print $band_max . "\t" . $band_min . "\n";
}
elsif( $metal == 1 || $metal == 2 )
{
  my $fermi = 'no';
  my $units;
  while( 1 )
  {
    die "Parsing $datafile failed\n" unless( my $line = <IN> );
    if( $line =~ m/\<UNITS_FOR_ENERGIES UNITS=\"(\w+)/ )
    { 
      $units = $1;
    }
    if( $line =~ m/\<FERMI_ENERGY/ )
    {
      $line = <IN>;
      $line =~ m/([-]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ or die "$line";
      if( $fermi eq 'no' ) 
      {
        $fermi = $1;
      }
      else
      {
        $fermi = $1 if( $fermi < $1 );
      }
    }

    last if( $line =~ m/<EIGENVALUES>/ );
  }

  unless( $fermi =~ m/no/ )
  {
    if( $units =~ m/ryd/i )
    {
      $fermi /= 2;
    }
    elsif( $units =~ m/eV/i )
    {
      $fermi /= 13.60569253 / 2;
    }
    print "$fermi\n";
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

  my $loud = 1;
  foreach my $file_stub (@filelist)
  {
    $file_stub =~ s/^\.//;
    $file_stub =~ s/^\///;
    my $file = $work_dir . '/' . $prefix . ".save/"  . $file_stub . "\n";
#    print $file;

    if( $fermi =~ m/no/ ) 
    {
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
        if( $line =~ m/([-]?\d+\.\d+[Ee]?[+-]?\d*)/ )
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
    }
    else  # Was able to load the fermi level
    {
      my @energy;
      open IN, $file or die "Failed to open $file\n$!";
      while( 1 )
      {
        die "Parsing $file failed\n" unless( my $line = <IN> );
        last if( $line =~ m/<EIGENVAL/ );
      }
      while( 1 )
      {
        die "Parsing $file failed\n" unless( my $line = <IN> );
        if( $line =~ m/([-]?\d+\.\d+[Ee]?[+-]?\d*)/ )
        {
          push @energy, $1;
        }
        last if( $line =~ m/<\/EIGENVAL/ );
      }
      close IN;
      my $count = 0;
      my $min_count = 0;
      foreach my $i (@energy)
      {
        if( $loud == 1 ) 
        {
          print "0\t$i\n" if ( $i < $fermi );
          print "1\t$i\n" if ( $i > $fermi );
        }
        $count++ if( $i < $fermi );
#        $min_count++ if( $i > $fermi );
      }
      # for compatibility with occ verison
      $min_count = $count;
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

    }
    $loud = 0;
  #  print "$file_stub:\t$count\t $min_count\n";
  }
  $band_min++;
  # band_min is the lowest count whereas band_max is the highest count
  #   Therefore brange will be 
  #   1         band_max
  #   band_min  total bands

  # pad band_max by 1
  $band_max;
  print $band_max . "\t" . $band_min . "\n";
}
else
{
  die "Error in metal flag. Programmer is to blame\n";
}


open OUT, ">brange.stub" or die "Failed to open brange.stub for writing\n$!";
print OUT "1    $band_max\n$band_min    ";
close OUT;

exit 0;

