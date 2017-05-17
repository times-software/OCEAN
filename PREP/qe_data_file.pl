#!/usr/bin/perl


use strict;

my $MNGKV = 0;
my $NKPT = 0;
my $NBND = 0;
my $NSPN = 0;
my $NSPINOR = 0;


open IN, $ARGV[0] or die "$!";

while( my $line=<IN> )
{
  if( $line =~ m/<MAX_NUMBER_OF_GK-VECTORS/ )
  {
    <IN> =~ m/(\d+)/;
    $MNGKV = $1;
    <IN>;
  }
  elsif( $line =~ m/<NUMBER_OF_K-POINTS/ )
  {
    <IN> =~ m/(\d+)/;
    $NKPT = $1;
    <IN>;
  }
  elsif( $line =~ m/<NUMBER_OF_BANDS/ )
  {
    <IN> =~ m/(\d+)/;
    $NBND = $1;
    <IN>;
  }
  elsif( $line =~ m/<NUMBER_OF_SPIN_COMPONENTS/ )
  {
    <IN> =~ m/(\d+)/;
    $NSPN = $1;
    <IN>;
  }
  elsif( $line =~m/<NON-COLINEAR_CALCULATION/ )
  {
    if( <IN> eq 'T' )
    {
      $NSPINOR = 2;
    }
    else
    {
      $NSPINOR = 1;
    }
  }

  last if( $MNGKV * $NKPT * $NBND * $NSPN * $NSPINOR > 0 );
}
close IN;

open OUT, ">qe_data.txt" or die "$!";
print OUT "$NBND $MNGKV $NSPN $NSPINOR $NKPT\n";

if( scalar @ARGV > 1 )
{
  $MNGKV = 0;
  $NKPT = 0;
  $NBND = 0;
  $NSPN = 0;
  $NSPINOR = 0;

  open IN, $ARGV[1] or die "$!";

  while( my $line=<IN> )
  {
    if( $line =~ m/<MAX_NUMBER_OF_GK-VECTORS/ )
    {
      <IN> =~ m/(\d+)/;
      $MNGKV = $1;
      <IN>;
    }
    elsif( $line =~ m/<NUMBER_OF_K-POINTS/ )
    {
      <IN> =~ m/(\d+)/;
      $NKPT = $1;
      <IN>;
    }
    elsif( $line =~ m/<NUMBER_OF_BANDS/ )
    {
      <IN> =~ m/(\d+)/;
      $NBND = $1;
      <IN>;
    }
    elsif( $line =~ m/<NUMBER_OF_SPIN_COMPONENTS/ )
    {
      <IN> =~ m/(\d+)/;
      $NSPN = $1;
      <IN>;
    }
    elsif( $line =~m/<NON-COLINEAR_CALCULATION/ )
    {
      if( <IN> eq 'T' )
      {
        $NSPINOR = 2;
      }
      else
      {
        $NSPINOR = 1;
      }
    }

  last if( $MNGKV * $NKPT * $NBND * $NSPN * $NSPINOR > 0 );
  }
  close IN;
}

# print two lines regardless (if only one DFT then this will be a duplicate)
print OUT "$NBND $MNGKV $NSPN $NSPINOR $NKPT\n";
close OUT;

exit 0;
