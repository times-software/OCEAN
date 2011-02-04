#!/usr/bin/perl
#
# Written by J. Vinson  Dec. 08
# This reads through enkfile, grabs all the occ_band energies and sorts
# counts through based on the number of kpts and the number of electrons per unit cell
# 
# If enkfile writeout changes this will break
# The number of kpts should always be even if the number of electrons is odd.

use strict;

open BRANGE, "brange.ipt" or die;
<BRANGE> =~ m/(\d+)\s+(\d+)/ ;
my @brange = ($1,$2);
<BRANGE> =~ m/(\d+)\s+(\d+)/ ;
$brange[2] = $1;
$brange[3] = $2;
close BRANGE;
print "$brange[3]\t$brange[2]\t$brange[1]\t$brange[0]\n";
open NE, "nelectron" or die;
<NE> =~ m/(\d+)/;
my $ne = $1 / 2;
print "$ne\n";
close NE;

open NKPT, "nkpts" or die;
<NKPT> =~ m/(\d+)/;
my $nkpt = $1 ;#/ 2;
print "$nkpt\n";
close NKPT;

$ne *= $nkpt;

#open ENK, "enkfile" or die;
open ENK, "enk_un" or die;
my $line;
my @energy;
my @Cenergy;
while ($line = <ENK> ) {
#for (my $kptcount = 0; $kptcount < $nkpt; $kptcount ++) {
# for (my $i = 0; $i <($brange[1])/3; $i++ ) {
#  $line = <ENK> or die; 
  chomp($line);
  $line =~ m/^\s*(.+)/ ;
  push @energy, split (/ +/, $1);
#  print "occ\n";
# } 
# for (my $i = 0; $i <($brange[3]-$brange[2]+1)/3; $i++ ) {
#  $line = <ENK> or die;
#  print "unocc\n";
# }
}
close ENK;
#my $bah = ($brange[3]-$brange[2]+1)/3;
#print "$bah\n";

my @sorted = sort { $a <=> $b } @energy;
open TEST,">test";
foreach (@sorted) {
 print TEST "$_\n";
}


my $ldagap = ($sorted[$ne]-$sorted[$ne-1]);
my $efermi = $ldagap / 2 + $sorted[$ne-1];

#print "$nkpt\t$ne\n";
print "$sorted[$ne-2]\t$sorted[$ne-1]\t$sorted[$ne]\t$sorted[$ne+1]\n";
open FERMI, ">efermiinrydberg.ipt" or die;
print FERMI $efermi . "\n";
close FERMI;

open LDA, ">ldagap" or die;
printf LDA "%.6f\n", (($sorted[$ne]-$sorted[$ne-1])*13.605698) ;
close LDA;

open ESHIFT, ">eshift.ipt" or die;
printf ESHIFT "%.6f\n", (-$sorted[$ne]*13.605698) ;
close ESHIFT;

exit 0;
