#!/usr/bin/perl

use strict;

# bvecs are defined without the traditional 2pi.

my @avecs;
open AVECS, "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n";
for (my $i = 0; $i < 3; $i++ ) {
 <AVECS> =~ m/([+-]?\d*\.?\d+([Ee][+-]\d+)?)\s+([+-]?\d*\.?\d+([Ee][+-]\d+)?)\s+([+-]?\d*\.?\d+([Ee][+-]\d+)?)/ or die "Failed to parse AVECS\n";
 $avecs[$i][0] = $1/(2*4*atan2(1,1));
 $avecs[$i][1] = $3/(2*4*atan2(1,1));
 $avecs[$i][2] = $5/(2*4*atan2(1,1));
}
close AVECS;


#print "$avecs[0][0]\t$avecs[0][1]\t$avecs[0][2]\n";
my $vol = 0;
open BVEC, ">bvecs" or die "Failed to open bvecs\n";
for( my $i = 0; $i < 3; $i++ ) {
  $vol += $avecs[0][$i] * ( $avecs[1][$i-2] * $avecs[2][$i-1] - $avecs[2][$i-2] * $avecs[1][$i-1]);
}

my @bvec;
for (my $i = 0; $i < 3; $i++ ) {
  for (my $j = 0; $j < 3; $j++ ) {
    $bvec[$i][$j] = ($avecs[$i-2][$j-2] * $avecs[$i-1][$j-1] - $avecs[$i-2][$j-1] * $avecs[$i-1][$j-2]) / $vol;
    printf BVEC  " % .16f ", $bvec[$i][$j];
  }
  print BVEC "\n";
}
close BVEC;
