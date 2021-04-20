#!/usr/bin/perl
# Copyright (C) 2015-2018 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/psp_parser\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }

open OUT, ">", "nelectrons" or die "Failed to open nelectrons for writing\n$!";
print OUT "-1\n";
close OUT;

my @typat;
open TYPAT, "typat" or die "Failed to open typat\n$!";
while (<TYPAT>)
{
  chomp;
  push @typat, split ' ';
}
close TYPAT;

my @psp;
open PSP, "pplist" or die "Failed top open pplist\n$!";
while (<PSP>)
{
  chomp;
  push @psp, split ' ';
}
close PSP;

#This is to deal with the psp name difference between abinit and QE
my $sufix = '';
open DFT, "dft" or die "Failed to open file dft\n$!";
my $dft = <DFT>;
close DFT;
if( $dft =~ m/qe/i )
{
  $dft = 'qe';
  $sufix = '.UPF';
}

open PPDIR, "ppdir" or die "Failed to open file ppdir\n$!";
my $ppdir = <PPDIR>;
close PPDIR;
chomp $ppdir;
$ppdir =~ s/'//g;
$ppdir .= '/';

my @electrons;
foreach my $ppfile (@psp)
{
  $ppfile .= $sufix unless( $ppfile =~ m/$sufix$/ );
  my $pp;
#  open $pp, "../$ppfile" or die "Failed to open ../$ppfile\n$!";
  open $pp, $ppdir . $ppfile or die "Failed to open $ppdir" . "$ppfile\n$!";

  my $nelectron = -1;
  if( $sufix eq '' )
  {
    $nelectron = abinitParser( $pp );
  }
  elsif( $sufix eq '.UPF' )
  {
    $nelectron = upfParser( $pp );
  }
  else
  {
    die "Un-supported psp\n";
  }
  die "Failed to parse $ppfile: $nelectron\n" if( $nelectron < 1 );

  print "$ppfile has $nelectron electrons\n";
  push @electrons, $nelectron;

  close $pp;
}

my $totalElectrons = 0;
foreach my $typat (@typat)
{
  if( $typat > scalar @electrons ) 
  {  die "Found entry in typat that is higher than number of pseudopotentials\n";  }
  $totalElectrons += $electrons[$typat-1];
}

print "Total electrons in system: $totalElectrons\n";


open OUT, ">", "nelectrons" or die "Failed to open nelectrons for writing\n$!";
print OUT "$totalElectrons\n";
close OUT;


sub upfParser
{
  my $pp = shift;
  my $n = -1;
  while( my $line = <$pp> )
  {
#    if( $line =~ m/z_valence=\"\s*(\d+.?\d*)\s*\"/i  || 
#        $line =~ m/z\s+valence=\"\s*(\d+.?\d*)\s*\"/i )
    if( $line =~ m/z_valence/i || $line =~ m/z valence/i )
    {
      
      if( $line =~ m/(\d+.?\d*[Ee]?\+?\d*)/ )
      {
#      print $line;
        $n = $1;
      }
      else
      {
        $n = -2;
      }
      last;
    }
#    elsif( $line =~ m/z_valence/i )
#    {
#      print $line;
#      $n = -2;
#      last;
#    }
  }

  return $n;
}

sub abinitParser
{
  my $pp = shift;
  my $n = -1;
  <$pp>;
  my $line = <$pp>;
  $line =~ m/^\s*(\d+\.?\d*)\s+(\d+\.?\d*)/ or die "Failed to parse\n$line";
  $n = $2;
  return $n;
}
