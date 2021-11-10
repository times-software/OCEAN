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
open OPTS, $ARGV[1] or die "Failed to open $ARGV[1]\n$!";

my $z;
my $z_eff;

if( $ARGV[0] =~ m/upf$/i ) {
  while( my $line = <PSP> ) {
#    print $line;
    if( $line =~ m/(\w+)\s+Element/ ) {
      my $element = $1;
      $z = &symb2z( $element );
      print "$element $z\n";
      last if( defined( $z_eff ) );
    } elsif( $line =~ m/(\d+\.?\d*)\s+Z valence/i ) {
      $z_eff = $1;
      print "$z_eff\n";
      last if( defined( $z ) );
    }
  }
  die "Failed to parse $ARGV[0]" unless( defined( $z) && defined( $z_eff ) );
} else {
  # First line in fhi is comments
  <PSP>;
  #second line has zatom and zion which can be floats? but must be at least ints
  <PSP> =~ m/^\s+(\d+(\.\d*)?)\s+(\d+(\.\d*)?)/ or die "Failed to parse $ARGV[0]\n";
  $z = $1;
  $z_eff = $3;
}

close PSP;


my $line = <OPTS>;
#<OPTS> =~ m/0*(\d+)/ or die "Failed to parse Z from $ARGV[1]\n";
$line =~ m/0*(\d+)/ or die "Failed to parse Z from $ARGV[1]\n";
my $z_fill = $1;
if( abs( $z - $z_fill ) > 0.01 ) 
{
  die "Mismatch between psp and fill: $z  $z_fill\n";
}

my @occs;
<OPTS> =~ m/(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse occs from $ARGV[1]\n";
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


sub symb2z
{
  my ($symb) = @_;
  my @z2symb =          ( '', 'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' ,
      'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
      'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
      'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
      'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
      'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
      'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
      'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
      'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu',
      'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',
      'Sg', 'Bh', 'Hs', 'Mt' );

  for( my $i=1; $i<scalar @z2symb; $i++ )
  {
    return $i  if( lc( $symb ) eq lc( $z2symb[$i] ) );
  }
  return -1;
}


