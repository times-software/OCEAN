#!/usr/bin/perl
# Copyright (C) 2017, 2021 OCEAN collaboration
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


# Load up all the radii we need
open RAD, "screen.shells" or die "Failed to open screen.shells\n";
my $line;
while( <RAD> )
{
  chomp;
  $line .= $_ . " ";
}
close RAD;
my @rads = split( /\s+/, $line );

my @Wshift;
my @Wsum;

open OUT, ">W.txt" or die "Failed to open pot.txt for writing\n$!";

for( my $h = 0; $h <= $#hfin; $h++ )
{
  my $filename = sprintf("zpawinfo/coreorbz%03in%02il%02i", $hfin[$h][1], $hfin[$h][2], $hfin[$h][3]);
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


  my $nn = $hfin[$h][2];
  my $ll = $hfin[$h][3];
  my $el = $hfin[$h][4];
  my $el_rank = $hfin[$h][5];
  my $string = sprintf("z%s%04d/n%02dl%02d",$el, $el_rank,$nn,$ll);
#  if( $debug != 0 )
#  {
#    print "$string\n";
#  }
  # W shift is in Ha., but we want to multiple by 1/2 anyway, so the units work out

  for( my $j = 0; $j < scalar @rads; $j++ )
  {
    my $rad_dir = sprintf("zR%03.2f", $rads[$j] );

    $filename = "$string/$rad_dir/ropt";
    open IN, $filename or die "Failed to open $filename\n$!";


    my @prad;
    my @pot;
    while( my $line = <IN> )
    {
      $line =~ m/(\S+)\s+\S+\s+\S+\s+(\S+)/;
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
      if( $j > 1 && $j < $#prad - 1 )
      {
        $interp = $pot[$j-2] * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j+1])
                            / (($prad[$j-2]-$prad[$j-1])*($prad[$j-2]-$prad[$j])*($prad[$j-2]-$prad[$j+1]) )
                + $pot[$j-1] * ($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j+1])*($rad[$i]-$prad[$j-2]) 
                            / (($prad[$j-1]-$prad[$j])*($prad[$j-1]-$prad[$j+1])*($prad[$j-1]-$prad[$j-2]) )
                + $pot[$j]   * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j+1])*($rad[$i]-$prad[$j-2]) 
                            / (($prad[$j]-$prad[$j-1])*($prad[$j]-$prad[$j+1])*($prad[$j]-$prad[$j-2]) )
                + $pot[$j+1] * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j-2])
                            / (($prad[$j+1]-$prad[$j-1])*($prad[$j+1]-$prad[$j])*($prad[$j+1]-$prad[$j-2]) );
      } 
      elsif( $j == $#prad - 1 || $j == 1 || $j == 0 ) 
      {
        my $run = $prad[$j+1]-$prad[$j];
        my $slope = ($pot[$j+1]-$pot[$j]) / $run;
        $interp += $slope * ( $rad[$i] - $prad[$j] );   
      }
      $sum += $rad[$i] * $xr1 * $wvfn[$i]**2 * $interp;
    }
    $Wshift[$h][$j] = $sum;
    $Wsum[$j] += $sum;
    print OUT "$sum    $rad_dir   $hfin[$h][1]  $hfin[$h][2]  $hfin[$h][3]  $hfin[$h][4]  $hfin[$h][5]\n"
  }

#  print "   $sum  \n";
}

close OUT;

exit 0;

