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
# The above statement is no good
# New plan for brange. Will find either fermi_energy or highest/lowest energy 
# in the data-file-scema.xml. The latter will be converted to Fermi. Then a simple
# higher/lower test for occupation. 
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
my $spin;
#if( $metal == 0 )
#{
  open IN, "../nspin" or die "Failed to open nspin\n$!\n";
  $spin = <IN>;
  close IN;
  $metal = 2 if( $spin =~ m/2/ );
#}


my $dft_split = 0;
if( -e "dft.split" )
{
  open IN, "dft.split" or die "Failed to open dft.split\n$!";
  $dft_split = 1 if( <IN> =~ m/t/i );
  close IN;
}

my $dft_shift = 0;
if( -e "qinunitsofbvectors.ipt" )
{
  open IN, "qinunitsofbvectors.ipt" or die "Failed to open qinunitsofbvectors.ipt\n$!";
  <IN> =~ m/(\S+)\s+(\S+)\s+(\S+)/ or die "Failed to parse qinunitsofbvectors.ipt\n$_";
  $dft_shift = 1 if( abs($1)+abs($2)+abs($3) > $tolerance );
  close IN;
}

$dft_split = 0 if( $dft_shift == 0 );

print "Split = " . $dft_split . ". Shift = " . $dft_shift . ".\n";


my $band_max = 0;
my $band_min;


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
  print $band_max . "\t" . $band_min . "\n";
}

my @energies;

my $nkpt = 0;
my $fermi = 'no';
my $highest = 'no';
my $lowest = 'no';
my $restart = 0;
while( my $line =<IN> )
{
  # in case things move out of order? At some point should just
  # load the xml into a perl hash or something
  $restart = 1 if( $line =~ m/<occupations size=\"\d+\">/ );
  $restart = 1 if( $line =~ m/<eigenvalues size=\"\d+\">/ );

  if( $line =~ m/\<highestOccupiedLevel\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
  {
    $highest = $1;
    last if( $dft_split == 1 );
    unless( $lowest eq 'no' )
    {
      $fermi = ($highest + $lowest ) / 2;
      last;
    }
  }
  elsif( $line =~ m/\<lowestUnoccupiedLevel\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
  {
    next if( $dft_split == 1 );
    $lowest = $1;
    unless( $highest eq 'no' )
    {
      $fermi = ($highest + $lowest ) / 2;
      last;
    }
  }
  elsif( $line =~ m/\<fermi_energy\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
  {
    $fermi = $1;
    last;
  }
}

#print "Fermi energy: $fermi\n";

if( $restart == 1 )
{
  close IN;
  open IN, $datafile or die "Failed to open $datafile\n$!";
}

if( $dft_split == 1 )
{
  
  my $datafileSplit = $work_dir . '/' . $prefix . "_shift.save/data-file-schema.xml";
  open FH, $datafileSplit or die "Failed to open $datafileSplit\n$!";

  if( $fermi eq 'no' )
  {
    while( my $line =<FH> )
    {
      if( $line =~ m/\<lowestUnoccupiedLevel\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $lowest = $1;
        $fermi = ($highest + $lowest ) / 2;
        last; 
      }
    }
  }
  close FH;
}
else
{ print "NO SPLIT\n"; }

print "Fermi energy: $fermi\n";

for( my $isplit = 0; $isplit <= $dft_split; $isplit++ )
{
  while( my $line = <IN> )
  {

    if( $line =~ m/<eigenvalues size=\"\d+\">/ )
    {
      $nkpt ++;
      my $eigs = '';
      $line = <IN>;
      until( $line =~ m/eigenvalues/ )
      {
        chomp $line;
        $eigs .= $line . ' ';
        $line = <IN>;
      }

      my @eigs = split( ' ', $eigs );
      push @energies, \@eigs;

      #occupations from Fermi
      if( $spin == 2 )
      {
        my $start = 0;
        my $stop = (scalar @eigs) / 2;
        for( my $j = 0; $j <= 1; $j++ )
        {
          my $count = 0;
          for( my $i = $start; $i < $stop; $i++ )
          {
            if( $eigs[$i] < $fermi )
            {
              $count++;
            }
            else
            {
              last;
            }
          }
          
  #        print "$count $start $stop\n";
          # when does count matter
          # 
          if( ( $dft_shift == 1 && $dft_split == 0 && $nkpt%2 == 1 ) ||
              ( $dft_shift == 1 && $dft_split == 1 && $isplit == 0 ) || $dft_shift == 0 )
          {
            $band_max = $count if( $count > $band_max );
          }
   
          if( ( $dft_shift == 1 && $dft_split == 0 && $nkpt%2 == 0 ) ||
              ( $dft_shift == 1 && $dft_split == 1 && $isplit == 1 ) || $dft_shift == 0 )
          {
            $band_min = ($count+1) if( !defined $band_min );
            $band_min = ($count+1) if( ( $count + 1 ) < $band_min );
          }
          $start = $stop;
          $stop = scalar @eigs;
        }
      }
      else
      {
        my $count = 0;
        for( my $i = 0; $i < scalar @eigs; $i++ )
        {
          if( $eigs[$i] < $fermi )
          {
            $count++;
          }
          else
          {
            last;
          }
        }
        
        if( ( $dft_shift == 1 && $dft_split == 0 && $nkpt%2 == 1 ) ||
            ( $dft_shift == 1 && $dft_split == 1 && $isplit == 0 ) || $dft_shift == 0 )
        {
          $band_max = $count if( $count > $band_max ); 
        }
        if( ( $dft_shift == 1 && $dft_split == 0 && $nkpt%2 == 0 ) ||
            ( $dft_shift == 1 && $dft_split == 1 && $isplit == 1 ) || $dft_shift == 0 )
        {
          $band_min = ($count+1) if( !defined $band_min );
          $band_min = ($count+1) if( ( $count + 1 ) < $band_min );
        }
      }
    }
  }
  close IN;

  if( $isplit < $dft_split )
  {
    my $datafileSplit = $work_dir . '/' . $prefix . "_shift.save/data-file-schema.xml";
    open IN, $datafileSplit or die "Failed to open $datafileSplit\n$!";
    print "$datafileSplit\n";
  }
}


$nkpt /= 2 if( $dft_split == 1 );
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
  my $stop = $band_max - 1;
  my $delim = " ";
  my @eslice = @{ $energies[$k] }[ $start .. $stop ];
  # move to Ryd
  foreach my $x (@eslice) { $x = $x * 2; }
  while (my @x = splice @eslice, 0, $n) 
  {
     print ENK join($delim, @x), "\n";
  }   
#  print "\n";
#    print ENK "\n";
  
  
  my $kk = $k;
  if( $dft_shift == 1 )
  {
    if( $dft_split == 0 )
    {
      $k++;
      $kk = $k;
    }
    else
    {
      $kk = $k + $nkpt;
    }
  }
  $start = $band_min - 1;
  if( $spin == 2 )
  {
    $stop = scalar @{ $energies[$kk] }/2 - 1;
  }
  else
  {
    $stop = scalar @{ $energies[$kk] } - 1;
  }
  @eslice = @{ $energies[$kk] }[ $start .. $stop ];
  # move to Ryd
  foreach my $x (@eslice) { $x = $x * 2; }
  while (my @x = splice @eslice, 0, $n) 
  {
    print ENK join($delim, @x), "\n";
  }   
}


#NEED TO REPLICATE FROM ABOVE

if( $spin == 2 )
{
  for( my $k = 0; $k < $nkpt; $k++ )
  {
    my $n = 3;
    my $start = scalar @{ $energies[$k] }/2 ;
    my $stop = $band_max + $start - 1;
    my $delim = " ";
    my @eslice = @{ $energies[$k] }[ $start .. $stop ];
    # move to Ryd
    foreach my $x (@eslice) { $x = $x * 2; }
    while (my @x = splice @eslice, 0, $n) 
    {   
       print ENK join($delim, @x), "\n";
    }   
    print "$k $energies[$k][$start] $start ";
        
    my $kk = $k;
    $start += $band_min-1 ;
    if( $dft_shift == 1 )
    {
      if( $dft_split == 0 )
      {
        $k++;
        $kk = $k;
      }
      else
      { 
        $kk = $k + $nkpt;
        $start = $band_min-1+ (scalar @{ $energies[$kk] }/2) ;
      }
    }
        
    $stop = scalar @{ $energies[$kk] } - 1; 
    @eslice = @{ $energies[$kk] }[ $start .. $stop ];
    # move to Ryd
    foreach my $x (@eslice) { $x = $x * 2; }
    while (my @x = splice @eslice, 0, $n)
    {   
      print ENK join($delim, @x), "\n";
    }     
    print "  $energies[$kk][$start] $start\n";
  } 
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

