#!/usr/bin/perl
# Copyright (C) 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
use strict;

#die "Expected usage qe2rhoofr.pl infile outfile\n" unless( $ARGV > 1 );

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

my $ngx = -1; my $ngy; my $ngz;
my $ntypat; my $natom;
open IN, "$infile" or die "Failed to open infile\n$!";
while( my $line=<IN> )
{
  # density fft grid, other fft grid, natom, ntypat
  if( $line =~ m/^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+/ )
  {
    $ngx = $1;
    $ngy = $2;
    $ngz = $3;
    $natom = $7;
    $ntypat = $8;
    last;
  }
}

if( $ngx < 0 ) 
{
  die "Failed to parse $infile\n";
}

# lattice info
<IN> =~ m/^\s*(\d+)/ or die "Failed to parse lattice line in $infile\n";
if( $1 == 0 ) {
  <IN>; # unit vec x
  <IN>; # unit vec y
  <IN>; # unit vec z
}
<IN>; # volume, nelectron, nbands? something?

for( my $i = 0; $i < $ntypat; $i++ )
{
  <IN>;
}

for( my $i = 0; $i < $natom; $i++ )
{
  <IN>;
}

my $rhoofr;

while( my $line=<IN> )
{
  chomp( $line );
  $rhoofr .= $line . " ";
}
close IN;

my @rhoofr = split /\s+/, $rhoofr;
# First line/element is blank
#print $rhoofr[0] . "\n";
print $rhoofr[1] . "\n";

if( scalar @rhoofr != $ngz*$ngy*$ngx+1 ) 
{
  print scalar @rhoofr - 1 . "\t". $ngz*$ngy*$ngx . "\n";
  die "Mismatch between read in rhoofr and FFT dims!\n";
}

open OUT, ">$outfile" or die "Failed to open $outfile for writing.\n$!";
printf OUT  "%8s %8s %8s %s\n", "i1", "i2", "i3", "data";

my $iter = 1;
for( my $k = 1; $k <= $ngz; $k++ )
{
  for( my $j = 1; $j <= $ngy; $j++ )
  {
    for( my $i = 1; $i <= $ngx; $i++ )
    {
      printf OUT "%8i %8i %8i %.9e\n", $i, $j, $k, $rhoofr[$iter];
      $iter++;
    }
  }
}

close OUT;

my $fftFile = "nfft";
$fftFile = "nfft.pot" if( $outfile =~ m/potofr/ );
open OUT, ">", $fftFile or die $!;
printf OUT "%8i %8i %8i\n", $ngx, $ngy, $ngz;
close OUT;
exit 0;
