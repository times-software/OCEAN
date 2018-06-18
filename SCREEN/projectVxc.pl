#!/usr/bin/perl
# Copyright (C) 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

my @hfin;
open IN, "hfinlist" or die "Failed to open hfinlist\n$!";
while( my $line = <IN> )
{
  chomp $line;
  push @hfin, [split ' ', $line];
}
close IN;


open OUT, ">pot.txt" or die "Failed to open pot.txt for writing\n$!";

for( my $h = 0; $h <= $#hfin; $h++ )
{
  my $filename = sprintf("../zpawinfo/coreorbz%03in%02il%02i", $hfin[$h][1], $hfin[$h][2], $hfin[$h][3]);
#  print $filename . "\n";
#  exit 0;
#  open IN, "../zpawinfo/coreorbz009n01l00" or die;
  open IN, $filename or die "Failed to open $filename\n$!";

  <IN>;
  my @rad;
  my @wvfn;

  while ( my $line = <IN> )
  {
    $line =~ m/(\S+)\s+(\S+)/;
    push @rad, $1;
    push @wvfn, $2;
  }

  close IN;

  my $rat = $rad[-1]/$rad[0];
  my $dl = log( $rat )/ ($#rad);
  my $xrat = exp( $dl );

  my $xr1 = sqrt( $xrat ) - sqrt( 1.0/$xrat);
#  print "  $rat\t$dl\t$xrat\t$xr1\n";

  my $rmin = $rad[0] / $xrat;

  ### TEST CODE
  my $sum = 0;
  my $sum2 = 0;
  for( my $i = 0; $i <= $#rad; $i++ )
  {
    my $temp_rad = $rmin * $xrat**($i+1);
  #  print "$rad[$i]\t$temp_rad\t" . ($rad[$i]-$temp_rad) . "\n";
    $sum += $rad[$i] * $rad[$i] * $rad[$i] * $xr1 * $wvfn[$i]**2;
    $sum2 += $rad[$i] * $xr1 * $wvfn[$i]**2;
  }


  #print "$sum\t" . 4*3.14159*$sum . "\n";
  #print "$sum2\t" . 4*3.14159*$sum2 . "\n";
  print "$sum2\n";
  #######


  $filename = sprintf("avg%2s%04i", $hfin[$h][4], $hfin[$h][5]);
  open IN, $filename or die "Failed to open $filename\n$!";
#  open IN, "avgF_0001" or die;
  my @prad;
  my @pot;
  while( my $line = <IN> )
  {
    $line =~ m/\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+/;
    push @prad, $1;
    push @pot, $2;
  }
  close IN;

  $sum = 0;
  my $j = 0;
  for( my $i = 0; $i <= $#rad; $i++ )
  {
    while( $rad[$i] > $prad[$j] )
    {
      die if( $j > $#prad );
      $j++;
    }
    my $interp = @pot[$j];
    if( $j > 0 && $j < $#prad )
    {
      my $run = $prad[$j+1]-$prad[$j];
      my $slope = ($pot[$j+1]-$pot[$j]) / $run;
      $interp += $slope * ( $rad[$i] - $prad[$j] );   
    }
    $sum += $rad[$i] * $xr1 * $wvfn[$i]**2 * $interp;
  }

#  print "   $sum  \n";
  print OUT "$sum    $hfin[$h][1]  $hfin[$h][2]  $hfin[$h][3]  $hfin[$h][4]  $hfin[$h][5]\n"
}

close OUT;

exit 0;
