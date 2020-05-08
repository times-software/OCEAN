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
  push @znucl, split ' ';
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
  push @typat, split ' ';
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
if( 0 ) {
#This section is just the beginning for trying to add in support for 
# specifying coordinates in real units instead of crystal/reduced coordinates.
# It won't really work yet, but it also won't break anything (since it wouldn't have worked anyway)
open COORD, "coord" or die "Failed to open coord\n$!";
my $coord = <COORD>;
close COORD;
if( $coord =~ m/bohr/i ) #later match angstrom too, but would need to correct other parts of input?
{
  my @avecs;
  open AVECS, "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    my $line = <AVECS>;
    chomp $line;
    my @temp = split ' ', $line;
    push @avecs, \@temp;
  }
  close AVECS;
  
  my @invA;
  my $det = 0;
  for( my $i = 0; $i < 3; $i++ )
  {
    my $j = 0;
    $det += $avecs[$i][0] *
         ( $avecs[($i+1)%3][1]*$avecs[($i+2)%3][2]
          - $avecs[($i+1)%3][2]*$avecs[($i+2)%3][1] );
  }
  print $det . "\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    for( my $j = 0; $j < 3; $j++ )
    {
      $invA[$j][$i] = ( $avecs[($i+1)%3][($j+1)%3]*$avecs[($i+2)%3][($j+2)%3]
                    - $avecs[($i+1)%3][($j+2)%3]*$avecs[($i+2)%3][($j+1)%3] ) / $det;
    }
  }

  for( my $i = 0; $i < 3; $i++ )
  {
    my @line;
    for( my $j = 0; $j < 3; $j++ )
    {
      my $l = 0;
      for( my $k=0;$k<3;$k++ )
      {
        $l+= $avecs[$i][$k]*$invA[$k][$j];
      }
      push @line, $l;
    }
    printf "%f %f %f\n", $line[0], $line[1], $line[2];
  }

  for ( my $i = 0; $i < 3; $i++ )
  {
    printf "%f %f %f\n", $avecs[$i][0], $avecs[$i][1], $avecs[$i][2];
  }
  for ( my $i = 0; $i < 3; $i++ )
  {
    printf "%f %f %f\n", $invA[$i][0], $invA[$i][1], $invA[$i][2];
  }


  # In the future use exact coordinates instead of crystal later in the program
  open IN, "taulist" or die "Failed to open taulist\n$!";
  my @tau;
  while( my $line = <IN> )
  {
    chomp $line;
    my @temp = split ' ', $line;
    push @tau, \@temp;
  }
  close IN;
  open OUT, ">", "taulist" or die;
  for( my $i = 0; $i < scalar @tau; $i++ )
  {

    my @xyz;
    for( my $j = 0; $j<3; $j ++ )
    {
      for( my $k = 0; $k< 3; $k++ )
      {
          $xyz[$j] += $invA[$k][$j] * $tau[$i][$k]
      }
    }
    my $tol = 0.0; #00001; # no -0.00000 crap
    printf OUT "%24.16f   %24.16f   %24.16f\n", $xyz[0]+$tol, $xyz[1]+$tol, $xyz[2]+$tol;
  }
  close OUT;

}
}

exit 0;
