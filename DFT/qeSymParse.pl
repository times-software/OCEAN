#!/usr/bin/perl
# Copyright (C) 2019 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#


use strict;

my $allowTranslate = $ARGV[0];

my $dataFile;
if( -e "data-file-schema.xml" ) 
{
  $dataFile = "data-file-schema.xml";
}
elsif( -e "data-file.xml" )
{
  $dataFile = "data-file.xml";
}
else
{  die;  }
open IN, "$dataFile" or die "Failed to open $dataFile\n$!";

my $translateTol = 0.000001;

my $nsym = -1;
my @syms;

while( my $line = <IN> )
{
# <IN>=~ m/<symmetries>/ )

#  if( $line =~ m/<nsym>[\n]?(\d+)/ )
  if( $line =~ m/\<nsym\>(\d+)/ )
  {
    $nsym = $1;
    print $nsym . "\n";
  }
  if( $line =~ m/<symmetry>/ || $line =~ m/<SYMM\.\d+>/ )
  {
    my @tmp;
    my $condense;
    until( $line =~ m/\/symmetry/ || $line =~ m/<\/SYMM\.\d+>/)
    {
#      print $line;
      chomp($line);
      $condense .= $line . " ";
      $line = <IN> or die;
    }
    my $translate = 0;
#    print "$condense\n";
#    die "Fractional translation fail\n$condense\n" 
      if( $condense =~ m/fractional\_translation>\s*(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)</ ||
          $condense =~ m/FRACTIONAL\_TRANSLATION[^>]*>\s*(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+</ ){
      print "$1 $3 $5\n";
      if( abs($1) > $translateTol || abs($3) > $translateTol || abs($5) > $translateTol )
      {
        print "Translation detected\n";
        $translate = 1;
        $translate = 0 if( $allowTranslate == 1 );
      }
    }
    if( $translate == 0 ) 
    {
      if( 0 ) 
      {
      die "Rotation fail\n"
        unless( $condense =~ m/>\s*(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)\s+(-?\d+\.\d+([eE][+-]?\d+)?)<\/rot/ );
        push @tmp, $1;
        push @tmp, $3;
        push @tmp, $5;
        push @tmp, $7;
        push @tmp, $9;
        push @tmp, $11;
        push @tmp, $13;
        push @tmp, $15;
        push @tmp, $17;
      push @syms, \@tmp;
      }
      else
      {
      die "Rotation fail\n$condense\n" unless( $condense =~ m/>([\d\s-+eE\.]*)<\/rotation>/i );
      my $rot = $1;
#      print $rot . "\n";
      @tmp = split ' ', $rot;
      push @syms, \@tmp;

#      for( my $k = 0; $k < scalar @tmp; $k++ )
#      {
#        print "$tmp[$k]\n";
#      }
#      exit;
      }
    }
  }
}
close IN;

open OUT, ">", "sym.txt" or die;
print OUT ( scalar @syms )  . "\n";

for( my $j = 0; $j < scalar @syms; $j++ )
{
  for( my $i = 0; $i < 3; $i++ )
  {
    printf OUT "%.0f  %.0f  %.0f\n", $syms[$j][$i*3], $syms[$j][$i*3+1], $syms[$j][$i*3+2];
  #  printf "%.0f  %.0f  %.0f\n", $syms[$i*3][0], $syms[$i*3+1][0], $syms[$i*3+2][0];
  }
}
close OUT;
