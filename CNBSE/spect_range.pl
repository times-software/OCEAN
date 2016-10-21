#!/usr/bin/perl
#
## Copyright (C) 2015 OCEAN collaboration
##
## This file is part of the OCEAN project and distributed under the terms 
## of the University of Illinois/NCSA Open Source License. See the file 
## `License' in the root directory of the present distribution.
##
##
#
use strict;


# Figure out DFT energy range
my $dft_min;
my $dft_max;
open IN, "dft" or die "Failed to open dft\n$!";
my $dft = <IN>;
close IN;

if( $dft =~ m/obf/i )
{
  open IN, "q.out" or die "Failed to open q.out\n$!";
  while( my $line = <IN> )
  {
    last if( $line =~ m/done building Hamiltonian/ );
  }
  my $bonus_line = <IN>;
  $bonus_line .= <IN>;
  $bonus_line .= <IN>;
  $bonus_line =~ m/\d+\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)/ or die "Failed to parse q.out\n$_\n";
  $dft_min = $1;
  $dft_max = $2;
  close IN;
}
else
{
  open IN, "enkfile" or die "Failed to open enkfile\n$!";
  my $enk;
  while( my $line = <IN> )
  {
    chomp $line;
    $enk .= " " . $line;
  }
  close IN;
  my @energies = split /\s+/, $enk;
  my @sorted = sort{ $a <=> $b } @energies;
  $dft_min = $sorted[0] * 13.605;
  $dft_max = $sorted[-1] * 13.605;
}

open IN, "efermiinrydberg.ipt" or die "Failed to open efermiinrydberg.ipt\n$!";
my $efermi = <IN>;
close IN;
chomp $efermi;
$efermi *= 13.605;


print "$dft_min\t$efermi\t$dft_max\n";

my $is_xas;
open TMPFILE, "cnbse.mode" or die "Failed to open cnbse.mode\n";
my $mode = <TMPFILE>;
close TMPFILE;
chomp($mode);
if( lc($mode) =~ m/xes/ )
{
  $is_xas = 0;
}
elsif( lc($mode) =~ m/xas/ )
{
  $is_xas = 1;
}
else
{
  $is_xas = 1;
}

open IN, "cls_average" or die "$!";
my $cls = <IN>;
chomp $cls;
close IN;

open IN, "so.ipt" or die "$!";
my $so = <IN>;
chomp $so;
close IN;

open IN, "cnbse.broaden" or die "$!";
my $broaden = <IN>;
chomp $broaden;
close IN;

my $start;
my $stop;
my $niter;
if( $is_xas == 1 )
{
  $start = 0;
  $stop = $dft_max - $efermi;
}
else
{
  $start = $dft_min - $efermi;
  $stop = 0;
}

$start -= 0.35 * $so;
$stop += 0.70 * $so;
$start += $cls - 20;
$stop += $cls + 20;
print "$start\t$stop\n";

$start = sprintf("%.0f",($start-0.5));
my $niter = sprintf("%.0f",( $stop - $start ) / ($broaden/10) );
$stop = $start + $niter * ($broaden/10);

print "$niter\t$start\t$stop\n";
open OUT, ">cnbse.spect_range" or die "Failed to open cnbse.spect_range\n$!";
print OUT "$niter  $start  $stop";
  
exit 0;
