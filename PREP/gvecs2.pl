#!/usr/bin/perl
# Copyright (C) 2010 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

#         modG = (bv1(1)*b(1) + bv2(1)*b(2) + bv3(1)*b(3))**2 +
#     &          (bv1(2)*b(1) + bv2(2)*b(2) + bv3(2)*b(3))**2 +
#     &          (bv1(3)*b(1) + bv2(3)*b(2) + bv3(3)*b(3))**2 

use strict;

open BVEC, "bvecs";
my @bvecs;
for (my $i = 0; $i < 3; $i++ ) {
 <BVEC> =~ m/(\S+)\s+(\S+)\s+(\S+)/ or die; #at somepoint should match floats
 $bvecs[$i][0] = $1;
 $bvecs[$i][1] = $2;
 $bvecs[$i][2] = $3;
}
close BVEC;
#print "$bvecs[2][0]\t$bvecs[2][1]\t$bvecs[2][2]\n";


open NFFT, "nfft";
my @nfft;
<NFFT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die;
@nfft = ($1,$2,$3);
close NFFT;

my @bounds;
for (my $i = 0; $i < 3; $i++ ) {
 $bounds[$i] = int($nfft[$i]/3+.5);
}
#print "$bounds[0]\t$bounds[1]\t$bounds[2]\n";

my $mod;
open TMP, ">gvtemp";
for (my $i = -$bounds[0]; $i <= $bounds[0]; $i++ ) {
  for (my $j = -$bounds[1]; $j <= $bounds[1]; $j++ ) {
    for (my $k = -$bounds[2]; $k <= $bounds[2]; $k++ ) {
      $mod = (($bvecs[0][0]*$i + $bvecs[1][0]*$j + $bvecs[2][0]*$k)**2 +
              ($bvecs[0][1]*$i + $bvecs[1][1]*$j + $bvecs[2][1]*$k)**2 +
              ($bvecs[0][2]*$i + $bvecs[1][2]*$j + $bvecs[2][2]*$k)**2 );
      printf TMP "%i\t%i\t%i\t%.8f\n", $i,$j,$k,$mod; 
    }
  }
}
close TMP;

system("sort -n -k 4 gvtemp > gvtemp2");

open IN, "gvtemp2" or die;
my @gvec = <IN>;
close IN;

my $skip = 0;
my $mag = 0;
my @vecArray;
my $vec;
my $outsream = "";
my $GROUPcounter = 0;
my $GROUPnum = 1;
my $CONTINOUScounter = 1;
open GVEC, ">gvectors";
my $i;
foreach (@gvec) {
        if ($_ =~ /\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ ) {
                if ($mag != $4) {
                        print GVEC "\t$GROUPnum\t$GROUPcounter\t$mag\n";
                        $i = $CONTINOUScounter - $GROUPcounter;
                        foreach $vec (@vecArray) {
                                print GVEC "\t$i\t$GROUPcounter\t$vec\n";
                                $i++
                        }
                        print GVEC "\n";
                        @vecArray = ();
                        $mag = $4;
                        $GROUPnum++;
                        $GROUPcounter = 0;
                }
                push(@vecArray, "$1\t$2\t$3");
                $CONTINOUScounter++;
                $GROUPcounter++;
        }
        else {
                if ($skip) {
                        print "skipped line, could be bad.\n";
                }
                $skip++;
        }
}
print GVEC "\t$GROUPnum\t$GROUPcounter\t$mag\n";
$i = $CONTINOUScounter - $GROUPcounter;
foreach $vec (@vecArray) {
        print GVEC "\t$i\t$GROUPcounter\t$vec\n";
        $i++;
}
print GVEC "\n";

close GVEC;

exit 0;
