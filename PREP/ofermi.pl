#!/usr/bin/perl
# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
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

my $nspin = 1;
if( open NSPIN, "nspin" )
{
  <NSPIN> =~ m/(\d)/;
  $nspin = $1;
  close NSPIN;
}

$ne *= $nkpt;
$ne *= $nspin;

open ENK, "enk_un" or die;
my $line;
my @energy;
while ($line = <ENK> ) {
  chomp($line);
  push @energy, split( ' ', $line );
#  $line =~ m/^\s*(.+)/ ;
#  push @energy, split (/ +/, $1);
}
close ENK;

open ENK, "enk_sh" or die;
my @sh_energy;
while ($line = <ENK> ) {
  chomp($line);
  push @sh_energy, split( ' ', $line );
}
close ENK;


my @sorted = sort { $a <=> $b } @energy;
open TEST,">test";
foreach (@sorted) {
 print TEST "$_\n";
}

my @sh_sorted = sort { $a <=> $b } @sh_energy;
open TEST,">test_shifted";
foreach (@sh_sorted) {
 print TEST "$_\n";
}



my $ldagap = ($sh_sorted[$ne]-$sorted[$ne-1]);
my $efermi = $ldagap / 2 + $sorted[$ne-1];

#print "$nkpt\t$ne\n";
print "$sorted[$ne-2]\t$sorted[$ne-1]\t$sh_sorted[$ne]\t$sh_sorted[$ne+1]\n";
open FERMI, ">efermiinrydberg2.ipt" or die;
print FERMI $efermi . "\n";
close FERMI;

open LDA, ">ldagap" or die;
printf LDA "%.6f\n", $ldagap *13.605698;
#printf LDA "%.6f\n", (($sorted[$ne]-$sorted[$ne-1])*13.605698) ;
close LDA;

open ESHIFT, ">eshift.ipt" or die;
printf ESHIFT "%.6f\n", (-$sh_sorted[$ne]*13.605698) ;
close ESHIFT;

exit 0;
