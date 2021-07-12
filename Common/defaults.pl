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
  $0 =~ m/(.*)\/defaults\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }


print `pwd`;


# We can't currently default a blank in the input so fake it.
my $ncpus = 1;
open PARA_PREFIX, "para_prefix" or die "$!";
my $para_prefix = <PARA_PREFIX>;
if( $para_prefix =~ m/\#/ )
{
  close PARA_PREFIX;
  open PARA_PREFIX, ">para_prefix" or die "$!";
  print PARA_PREFIX "";
}
elsif ( $para_prefix =~ m/^\D*(\d+)/ )
{
  $ncpus = $1;
  print "GRABBED $ncpus as the number of cpus to run with!\n";
}
close PARA_PREFIX;

my @cpu_factors;
for( my $i =1; $i <= $ncpus; $i++ )
{
  push( @cpu_factors, $i ) unless ( $ncpus % $i );
}
if( scalar @cpu_factors == 1 )
{
  print "$ncpus has ", $#cpu_factors + 1, " factor\n";
}
else
{
  print "$ncpus has ", $#cpu_factors + 1, " factors\n";
}
foreach (@cpu_factors)
{ print "$_\n"; }

my @square_cpu_factors;
foreach my $cpu_count (@cpu_factors)
{
  push @square_cpu_factors, $cpu_count if( int(sqrt($cpu_count)) ** 2 == $cpu_count );
}
foreach (@square_cpu_factors)
{ print "$_\n"; }


open QE_POOL, ">pool_control" or die "Failed to open pool_control for writing\n$!\n";


my $tline;
open RSCALE, "rscale" or die;
$tline = <RSCALE>;
chomp($tline);
$tline =~  m/((\d+)?\.?\d+([eEfF][+-]?\d+)?)\s+((\d+)?\.?\d+([eEfF][+-]?\d+)?)\s+((\d+)?\.?\d+([eEfF][+-]?\d+)?)\s*$/ 
    or die "Failed to parse rscale!\n$tline\n";
my @rscale = ($1, $4, $7);
print "$1\t$4\t$7\n";
close RSCALE;

my @alength;
my @avec;
open RPRIM, "rprim" or die;
open AVECS, ">avecsinbohr.ipt" or die;
for (my $i = 0; $i < 3; $i++ ) 
{
  $tline = <RPRIM>;
  chomp( $tline );
  $tline =~  m/([+-]?(\d+)?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?(\d+)?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?(\d+)?\.?\d+([eEfF][+-]?\d+)?)\s*$/ 
      or die "Failed to parse a line of rprim!\n$tline\n";
  $avec[$i][0] = $1*$rscale[$i];
  $avec[$i][1] = $4*$rscale[$i];
  $avec[$i][2] = $7*$rscale[$i];
  $alength[$i] = sqrt( $avec[$i][0]**2 + $avec[$i][1]**2 + $avec[$i][2]**2 );
#  print AVECS $1*$rscale[0] . "  " . $4*$rscale[1] .  "  " . $7*$rscale[2] . "\n";
  printf AVECS "%.15f  %.15f  %.15f\n", $avec[$i][0], $avec[$i][1], $avec[$i][2];
  print "$1\t$4\t$7\n";
}
close RPRIM;
close AVECS;


# det(A) = det(A^T)
my $volume = $avec[0][0] * ($avec[1][1] * $avec[2][2] - $avec[2][1] * $avec[1][2] )
           - $avec[1][0] * ($avec[0][1] * $avec[2][2] - $avec[2][1] * $avec[0][2] )
           + $avec[2][0] * ($avec[0][1] * $avec[1][2] - $avec[1][1] * $avec[0][2] );
print "Volume:\t$volume\n";

my @bvec;
my @blen;

my $pref = 2*4*atan2(1,1)/$volume;

for (my $i = 0; $i < 3; $i++ ) {
  for (my $j = 0; $j < 3; $j++ ) {
    $bvec[$i][$j] = $pref * ($avec[$i-2][$j-2] * $avec[$i-1][$j-1] - $avec[$i-2][$j-1] * $avec[$i-1][$j-2]) ;
  }
  $blen[$i] = sqrt( $bvec[$i][0]**2 + $bvec[$i][1]**2 + $bvec[$i][2]**2 );
}


my %defaults;

open DEFAULTS, "$ENV{'OCEAN_BIN'}/defaults.h" or die "$!\nFailed to open defaults.h\n";
while ( my $line = <DEFAULTS>)
{
  $line =~ m/(\w+)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)/ 
      or die "Failed to parse line: $line";
#  $defaults{ $1 } = ( "$2", "$3", "$4" );
  push(@{$defaults{$1}}, $2 );
  push(@{$defaults{$1}}, $3 );
  push(@{$defaults{$1}}, $4 );
}
close DEFAULTS;
foreach my $input (keys %defaults)
{
  print "Key: $input\t";
  print "$defaults{$input}[0]\t$defaults{$input}[1]\t$defaults{$input}[2]\n";
}


open INPUT, "acc_level.ipt" or die "Failed to open acc_level.ipt\n$!";
<INPUT> =~ m/(\d)/ or die "$!";
my $acc_level = $1;
if( $acc_level < 1 || $acc_level > 3 )
{
  print "Accuracy requested (acc_level.ipt) is $acc_level.\nShould be between 1 and 3\n";
  print "Setting it to 1\n";
  $acc_level = 1;
}
$acc_level -= 1;
close INPUT;

my $input_content;
my $target;

$input_content = '';
open INPUT, "nkpt" or die "Failed to open nkpt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
my $kpt_tot;
if( $input_content =~ m/old/i )
{
  my @output;
  print "Defaults requested for kmesh.ipt\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    $target = $defaults{'kmesh'}[$acc_level];
    $target = $target / $alength[$i];
    $output[$i] = findval($target);
    if( $output[$i] == 1000 ) { die "kmesh too large\n"; }
  }
  open INPUT, ">nkpt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for kmesh.ipt:\t$output[0]\t$output[1]\t$output[2]\n";
  $kpt_tot = $output[0]*$output[1]*$output[2];
  $input_content = "$output[0]  $output[1]  $output[2]";
} 
elsif ( $input_content =~ m/^\s*(-(?=\.\d|\d)(?:0|[1-9]\d*)?(?:\.\d*)?(?:\d[eE][+\-]?\d+)?)/ ) {
  my $kden = $1;
  $kden *= -0.99999;
  my @kpt;
  printf  "Using kden=%.3f to determine k-mesh\n", $kden;
  my $klen = 3*$kden;
  for( my $i=0; $i<3; $i++ ) {
    print "$blen[$i]\n";
    $kpt[$i] = int( $kden * $blen[$i]) + 1;
    while( isAllowed( $kpt[$i] ) == 0 )
    {
      $kpt[$i]++;
    }
    $klen = $kpt[$i] / $blen[$i] if( $kpt[$i] / $blen[$i] < $klen );
  }
  printf "%3i  %3i  %3i  %16.5f\n", $kpt[0], $kpt[1], $kpt[2], $klen;

  $input_content = "$kpt[0]  $kpt[1]  $kpt[2]";
  open INPUT, ">nkpt" or die "$!\n";
  print INPUT "$kpt[0]  $kpt[1]  $kpt[2]\n";
  close INPUT;
  print "Defaults chosen for kmesh.ipt:\t$kpt[0]\t$kpt[1]\t$kpt[2]\n";
}

my @kpt;
if ( $input_content =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)/ )
{
  $kpt_tot = $1*$2*$3;
  $kpt[0] = $1;
  $kpt[1] = $2;
  $kpt[2] = $3;
}
else
{
  die "FAILED TO CORRECTLY PARSE nkpt\n";
}

if ( 0 ) {
  my @supervec;
  for( my $i =0; $i<3; $i++) {
    $supervec[$i][0] = $avec[$i][0] * $kpt[$i];
    $supervec[$i][1] = $avec[$i][1] * $kpt[$i];
    $supervec[$i][2] = $avec[$i][2] * $kpt[$i];
    print "$supervec[$i][0]  $supervec[$i][1]  $supervec[$i][2]\n";
  }
  my $rmax = sqrt( $supervec[0][0]**2 + $supervec[0][1]**2 + $supervec[0][2]**2 );
  for( my $i=0; $i<3; $i++ ) {
    my @tmp;
    $tmp[0] = $supervec[$i-2][1]*$supervec[$i-1][2] - $supervec[$i-1][1]*$supervec[$i-2][2];
    $tmp[1] = -$supervec[$i-2][0]*$supervec[$i-1][2] + $supervec[$i-1][0]*$supervec[$i-2][2];
    $tmp[2] = $supervec[$i-2][0]*$supervec[$i-1][1] - $supervec[$i-1][0]*$supervec[$i-2][1];
    print "$tmp[0] $tmp[1] $tmp[2]\n";
    my $tlen = sqrt( $tmp[0]**2 + $tmp[1]**2 + $tmp[2]**2 );
    my $r = $supervec[$i][0] * $tmp[0] + $supervec[$i][1] * $tmp[1] + $supervec[$i][2] * $tmp[2];
    $r /= $tlen;
    print "$i: $r\n";
    $rmax = $r if ( $r < $rmax );
  }
  $rmax /= 2;
  print "Max radius $rmax\n";
}


my $min_nscf_pool_size = $volume / 1400;
my $ideal_npools = 1;
if( -e "dft.nscf.poolsize" )
{
  open INPUT, "dft.nscf.poolsize";
  if( <INPUT> =~ m/(-?\d+)/ )
  {
    my $t = $1;
    $min_nscf_pool_size = $t if( $t > 0 );
  }
  $min_nscf_pool_size = $ncpus if( $min_nscf_pool_size > $ncpus );
  close INPUT;
}
print "Min pool size: $min_nscf_pool_size\n";
  
foreach (@cpu_factors)
{
  if( $ncpus / $_ >= $min_nscf_pool_size )
  {
    $ideal_npools = $_ unless( $kpt_tot % $_ );
  }
}
print "NKPTS: $kpt_tot, ideal pools: $ideal_npools\n";
if( $kpt_tot > (1.9*$ideal_npools ) && $ncpus / $ideal_npools > 4 )  # if ideal is too inefficient
{
  foreach (@cpu_factors)
  {
    if( $ncpus / $_ >= $min_nscf_pool_size )
    {
      $ideal_npools = $_ if( $_ <= $kpt_tot );
    }
  }
}
print "NKPTS: $kpt_tot, ideal pools: $ideal_npools\n";
print QE_POOL "nscf\t$ideal_npools\n";
# Not integrated cleanly right now
# Need to control pools for the interpolation routines
my $min_interp_pool_size = int( $volume / 800 );
print "Minimum size for interpolation pools $min_interp_pool_size\n";

my $pool_size = -1;
foreach( @square_cpu_factors )
{
  if( $_ >= $min_interp_pool_size && ($ncpus/$_) <= $kpt_tot )
  {
     $pool_size = $_;
     last;
  }
}
if( $pool_size == -1 )
{
  foreach( @cpu_factors )
  {
    if( $_ >= $min_interp_pool_size && ($ncpus/$_) <= $kpt_tot )
    {
       $pool_size = $_;
       last;
    }
  }
}
if( $pool_size == -1 )
{
  print "Might not be enough memory to run interpolation scheme\n";
  $pool_size = $ncpus;
}
print "KPTS: $kpt_tot, obf pool size: $pool_size\n";
print QE_POOL "interpolate kpt\t$pool_size\n";


  

$input_content = '';
open INPUT, "ngkpt" or die "Failed to open ngkpt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/old/i )
{
  my @output;
  print "Defaults requested for ngkpt\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    $target = $defaults{'ngkpt'}[$acc_level];
    $target = $target / $alength[$i];
    $output[$i] = findval($target);
    if( $output[$i] == 1000 ) { die "ngkpt too large\n"; }
  }
  open INPUT, ">ngkpt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for ngkpt:\t$output[0]\t$output[1]\t$output[2]\n";
  $kpt_tot = $output[0]*$output[1]*$output[2];
  $input_content = "$output[0]  $output[1]  $output[2]";
}
elsif ( $input_content =~ m/^\s*(-(?=\.\d|\d)(?:0|[1-9]\d*)?(?:\.\d*)?(?:\d[eE][+\-]?\d+)?)/ ) {
  my $kden = $1;
  $kden *= -0.99999;
  my @kpt;
  printf  "Using kden=%.3f to determine ngkpt\n", $kden;
  my $klen = 3*$kden;
  for( my $i=0; $i<3; $i++ ) {
    print "$blen[$i]\n";
    $kpt[$i] = int( $kden * $blen[$i]) + 1;
    while( isAllowed( $kpt[$i] ) == 0 )
    {
      $kpt[$i]++;
    }
    $klen = $kpt[$i] / $blen[$i] if( $kpt[$i] / $blen[$i] < $klen );
  }
  printf "%3i  %3i  %3i  %16.5f\n", $kpt[0], $kpt[1], $kpt[2], $klen;

  $input_content = "$kpt[0]  $kpt[1]  $kpt[2]";
  open INPUT, ">ngkpt" or die "$!\n";
  print INPUT "$kpt[0]  $kpt[1]  $kpt[2]\n";
  close INPUT;
  print "Defaults chosen for ngkpt:\t$kpt[0]\t$kpt[1]\t$kpt[2]\n";
}

if ( $input_content =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)/ )
{
  $kpt_tot = $1*$2*$3;
}
else
{
  die "FAILED TO CORRECTLY PARSE ngkpt\n";
}

my $ideal_npools = 1;
foreach (@cpu_factors)
{
  $ideal_npools = $_ unless( $kpt_tot % $_ );
}
print "NGKPTS: $kpt_tot, ideal pools: $ideal_npools\n";
print QE_POOL "scf\t$ideal_npools\n";


$input_content = '';
open INPUT, "obkpt.ipt" or die "Failed to open obkpt.ipt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/^\s*-1/ )
{
  my @output;
  print "Defaults requested for obkpt\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    $target = $defaults{'obkpt'}[$acc_level];
    $target = $target / $alength[$i];
    $output[$i] = findval($target);
    if( $output[$i] == 1000 ) { die "obkpt too large\n"; }
  }
  open INPUT, ">obkpt.ipt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for obkpt:\t$output[0]\t$output[1]\t$output[2]\n";
  $kpt_tot = $output[0]*$output[1]*$output[2];
}
elsif ( $input_content =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)/ )
{
  $kpt_tot = $1*$2*$3;
}
else
{
  die "FAILED TO CORRECTLY PARSE obkpt\n";
}
 
my $ideal_npools = 1;
foreach (@cpu_factors)
{
  $ideal_npools = $_ unless( $kpt_tot % $_ );
}
print "OBKPTS: $kpt_tot, ideal pools: $ideal_npools\n";
print QE_POOL "obf\t$ideal_npools\n";

my $nelectrons = -1;
if( -f "nelectrons" )
{
  open NE, "nelectrons" or die "Failed to open nelectrons\n$!";
  my $ne = <NE>;
  $nelectrons = sprintf "%.0f", $ne;
}

# Need to figure out the number of bands to use (for the BSE states)
# 1) User has asked for NBANDS
# 2) User has asked for and energy range (guessing time!)
open ERANGE, "dft_energy_range.ipt" or die "Failed to open dft_energy_range.ipt\n$!";
my $erange = <ERANGE>;
chomp($erange);
close ERANGE;
open NBANDS, "nbands" or die "Failed to open nbands\n$!";
my $nbands = <NBANDS>;
chomp($nbands);
close NBANDS;
die "Either dft_energy_range or nbands must be specified\n" 
    if( $erange <= 0 && $nbands <= 0);
if( $nbands <= -2 )
{
  $nbands = abs( $nbands );
  if ( $nelectrons < 1 ) 
  {
    $nbands += 0.036 * $volume;
  }
  else
  {
    $nbands += ($nelectrons/2);
  }
  open NBANDS, ">nbands" or die "Failed to open nbands for writing.\n$!";
  print NBANDS "$nbands\n";
  close NBANDS;
}
if( $nbands <= 0 || $nbands =~ m/range/i )
{
  print "Default requested for nbands. Energy range is $erange eV.\n";
  # First guess N conduction electrons
  $erange = $erange / 13.605;
  $nbands = 0.019 * $volume * ( $erange**(3/2) );
  print "      $nbands\n";
  # Then add a guess for valence
  if( $nelectrons < 1 ) 
  {
    $nbands += 0.036 * $volume;
  }
  else
  {
    $nbands += ($nelectrons/2);
  }
#  print "      $nbands\n";
  # 1.05 is a padding factor
  $nbands *= 1.05;
#  print "      $nbands\n";
  $nbands = int($nbands);
  # might as well pick an even number
  $nbands++ if( $nbands % 2 == 1 );
  print "Default chosen for nbands:\t$nbands\n";
  open NBANDS, ">nbands" or die "Failed to open nbands for writing.\n$!";
  print NBANDS "$nbands\n";
  close NBANDS;
}

# Bands for screening calculation
# 1) User has asked for a set number
# 2) Attempt to guess using erange
open ERANGE, "screen_energy_range.ipt" or die "Failed to open screen_energy_range.ipt\n$!";
my $erange = <ERANGE>;
chomp($erange);
close ERANGE;
open NBANDS, "screen.nbands" or die "Failed to open screen.nbands\n$!";
my $screen_nbands = <NBANDS>;
chomp($screen_nbands);
close NBANDS;
if( $screen_nbands <= 0 )
{
  $erange = 100 if( $erange <= 0 );
  print "Default requested for screen.nbands. Energy range is $erange eV.\n";
  $erange = $erange / 13.605;
  # First guess N conduction electrons
  $screen_nbands = 0.019 * $volume * ( $erange**(3/2) );
  # Then add a guess for valence
  if( $nelectrons < 1 )
  {
    $screen_nbands += 0.036 * $volume;
  }
  else
  {
    $screen_nbands += ($nelectrons/2);
  }
  # 1.05 is a padding factor
  $screen_nbands *= 1.05;
  $screen_nbands = int($screen_nbands);
  # might as well pick an even number
  $screen_nbands++ if( $screen_nbands % 2 == 1 );
  print "Default chosen for screen.nbands:\t$screen_nbands\n";
  open NBANDS, ">screen.nbands" or die "Failed to open nbands for writing.\n$!";
  print NBANDS "$screen_nbands\n";
  close NBANDS;
}

# OBF bands
open ERANGE, "obf_energy_range" or die "Failed to open obf_energy_range\n$!";
$erange = <ERANGE>;
chomp($erange);
close ERANGE;
open NBANDS, "obf.nbands" or die "Failed to open obf_nbands\n$!";
my $obf_nbands = <NBANDS>;
chomp($obf_nbands);
close NBANDS;
if( $obf_nbands <= 0 && $erange <= 0 )
{
  $obf_nbands = $nbands;
  print "Default chosen for obf_nbands:\t$obf_nbands\n";
} 
elsif( $obf_nbands <= 0 ) 
{
  print "Default requested for obf_nbands. Energy range is $erange eV.\n";
  $erange = $erange / 13.605;
  $obf_nbands = 0.01688686394 * $volume * ( $erange**(3/2) );
  $obf_nbands *= 1.05;
  $obf_nbands = int($obf_nbands);
  # 1.05 is a padding factor
  print "Default chosen for obf_nbands:\t$obf_nbands\n";
  open NBANDS, ">obf.nbands" or die "Failed to open nbands for writing.\n$!";
  print NBANDS "$obf_nbands\n";
  close NBANDS;
}

$input_content = '';
open INPUT, "xmesh.ipt" or die "Failed to open xmesh.ipt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/old/i )
{
  my @output;
  print "Defaults requested for xmesh.ipt\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    $target = $defaults{'xmesh'}[$acc_level];
    $target = $target * $alength[$i];
    $output[$i] = findval($target);
    if( $output[$i] == 1000 ) { die "xmesh too large\n"; }
  }
  open INPUT, ">xmesh.ipt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for xmesh.ipt:\t$output[0]\t$output[1]\t$output[2]\n";
  $input_content = "$output[0]  $output[1]  $output[2]";
}
else{
  $input_content =~ m/^\s*(-?\d+)/ or die "Failed to parse xmesh.ipt\n$input_content\n";
  my $xden = $1;
  if( $xden < 1 ) {
    $xden *= -0.99999;
    my @xpt;
    printf "Using xden=%.3f to determine x-mesh\n", $xden;
    my $xlen = 3*$xden;
    for( my $i=0; $i<3; $i++ ) {
#      my $alen = sqrt( $avec[$i][0]**2 + $avec[$i][1]**2 + $avec[$i][2]**2 );
      $xpt[$i] = int( $xden * $alength[$i]) + 1;
      while( isAllowed( $xpt[$i] ) == 0 )
      {
        $xpt[$i]++;
      }
      $xlen = $xpt[$i] / $alength[$i] if( $xpt[$i] / $alength[$i] < $xlen );
    }
    printf "%3i  %3i  %3i  %16.5f\n", $xpt[0], $xpt[1], $xpt[2], $xlen;

    $input_content = "$xpt[0]  $xpt[1]  $xpt[2]";
    open INPUT, ">xmesh.ipt" or die "$!\n";
    print INPUT "$xpt[0]  $xpt[1]  $xpt[2]\n";
    close INPUT;
    print "Defaults chosen for xmesh.ipt:\t$xpt[0]\t$xpt[1]\t$xpt[2]\n";
  }
}
# Now test xmesh to make sure it is large enough for the number of bands requested
$input_content =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Xmesh parsing failed!\n$input_content\n";
my $nxpts = $1*$2*$3;
my @orig_x;
$orig_x[0] = $1; $orig_x[1] = $2; $orig_x[2] = $3;
if( $nxpts < 1.8 * $nbands )
{
  print "Warning! Requested xmesh is too small for requested bands.\nAttempting to increase xmesh.\n";
  print "Due to orthogonalization the number of x-points must be larger than the number of bands.\n";
  print "OCEAN will increase the xmesh until the ratio exceeds 1.8.\n";

  my $scale = 1.1;
  my @output;
  while( $nxpts < 1.8 * $nbands )
  {
    for( my $i = 0; $i < 3; $i++ )
    {
      $target = $orig_x[$i] * $scale;
      $output[$i] = findval($target);
      if( $output[$i] == 1000 ) { die "xmesh too large\n"; }
    }
    $nxpts = $output[0] * $output[1] * $output[2];
    $scale += 0.1;
  }
  open INPUT, ">xmesh.ipt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "New xmesh choosen:\t$output[0]\t$output[1]\t$output[2]\n";
  $input_content = "$output[0]  $output[1]  $output[2]";
}



$input_content = '';
open INPUT, "ham_kpoints" or die "Failed to open ham_kpoints\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/^\s*-1/ )
{
  my @output;
  print "Defaults requested for ham_kpoints\n";
#  for( my $i = 0; $i < 3; $i++ )
#  {
#    $target = $defaults{'ham_kpoints'}[$acc_level];
#    $target = $target / $alength[$i];
#    $output[$i] = findval($target);
#    if( $output[$i] == 1000 ) { die "ham_kpoints too large\n"; }
#  }
  $output[0] = 4; $output[1] = 4; $output[2] = 4;
  open INPUT, ">ham_kpoints" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for ham_kpoints:\t$output[0]\t$output[1]\t$output[2]\n";
}

$input_content = '';
open INPUT, "screen.nkpt" or die "Failed to open screen.nkpt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/old/i )
{
  my @output;
  print "Defaults requested for screen.nkpt\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    $target = $defaults{'screen'}[$acc_level];
    $target = $target / $alength[$i];
    $output[$i] = findval($target);
    if( $output[$i] == 1000 ) { die "screen.nkpt too large\n"; }
  }
  open INPUT, ">screen.nkpt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for screen.nkpt:\t$output[0]\t$output[1]\t$output[2]\n";
  $kpt_tot = $output[0]*$output[1]*$output[2];
  $input_content = "$output[0]  $output[1]  $output[2]";
}
elsif ( $input_content =~ m/^\s*(-(?=\.\d|\d)(?:0|[1-9]\d*)?(?:\.\d*)?(?:\d[eE][+\-]?\d+)?)/ ) {
  my $kden = $1;
  $kden *= -0.99999;
  my @kpt;
  printf  "Using kden=%.3f to determine screen.nkpt\n", $kden;
  my $klen = 3*$kden;
  for( my $i=0; $i<3; $i++ ) {
    print "$blen[$i]\n";
    $kpt[$i] = int( $kden * $blen[$i]) + 1;
    while( isAllowed( $kpt[$i] ) == 0 )
    {
      $kpt[$i]++;
    }
    $klen = $kpt[$i] / $blen[$i] if( $kpt[$i] / $blen[$i] < $klen );
  }
  printf "%3i  %3i  %3i  %16.5f\n", $kpt[0], $kpt[1], $kpt[2], $klen;

  $input_content = "$kpt[0]  $kpt[1]  $kpt[2]";
  open INPUT, ">screen.nkpt" or die "$!\n";
  print INPUT "$kpt[0]  $kpt[1]  $kpt[2]\n";
  close INPUT;
  print "Defaults chosen for screen.nkpt:\t$kpt[0]\t$kpt[1]\t$kpt[2]\n";
}

if ( $input_content =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)/ )
{
  $kpt_tot = $1*$2*$3;
}
else
{
  die "FAILED TO CORRECTLY PARSE screen.nkpt\n";
}

my $ideal_npools = 1;
foreach (@cpu_factors)
{
  $ideal_npools = $_ unless( $kpt_tot %  $_ );
}
print "PAW: $kpt_tot, ideal pools: $ideal_npools\n";
print QE_POOL "screen\t$ideal_npools\n";

# Not integrated cleanly right now
# Need to control pools for the interpolation routines
$min_interp_pool_size = int( $volume / 800 );
if( $min_interp_pool_size < 1 )
{
  $min_interp_pool_size = 1;
}
print "Minimum size for interpolation pools $min_interp_pool_size\n";

$pool_size = -1;
foreach( @cpu_factors )
{
  if( $_ >= $min_interp_pool_size && ( $kpt_tot % ($ncpus/$_) == 0 ) ) 
#($ncpus/$_) <= $kpt_tot ) )
  {
     print "$ncpus\t$_\t$kpt_tot\n";
     $pool_size = $_;
     last;
  }
}
if( $pool_size == -1 ) 
{
  print "Might not be enough memory to run interpolation scheme\n";
  $pool_size = $ncpus;
}
print "PAW: $kpt_tot, obf pool size: $pool_size\n";
print QE_POOL "interpolate screen\t$pool_size\n";
print QE_POOL "total\t$ncpus\n";


close QE_POOL;

# qinunits
$input_content = '';
open CALC, "calc" or die "Failed to open open calc\n$!";
while (<CALC>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close CALC;

if( $input_content =~ m/VAL/i || $input_content =~ m/RXS/i )
{
  open QIN, "qinunitsofbvectors.ipt" or die "Failed to open qinunitsofbvectors.ipt\n$!";
  $input_content = '';
  while (<QIN>)
    {
      $input_content .= $_;
    }
  $input_content =~ s/\n/ /g;
  close QIN;

  $input_content =~ m/([+-]?\d+(\.\d+)?)\s+([+-]?\d+(\.\d+)?)\s+([+-]?\d+(\.\d+)?)/ 
    or die "Failed to parse qinunitsofbvectors.ipt\n\t\t$input_content\n";
  if( abs( $1 ) + abs( $2 ) + abs( 3 ) < 0.0000000001 )
  {
    print "Valence or RIXS calculation requires non-zero q-vector\n";
    print "Assigning default q-vec of ( 0.001, 0, 0 ). You should edit your input\n";
    open QIN, ">qinunitsofbvectors.ipt" or die "Failed to open qinunitsofbvectors.ipt\n$!";
    print QIN "0.001 0.0 0.0\n";
    close QIN;
  }
}


# epsilon can't be 1
open EPS, "epsilon" or die "Failed to open epsilon\n$!";
my $epsilon = <EPS>;
chomp $epsilon;
close EPS;
unless( $epsilon > 1.000001 || $epsilon =~ m/dfpt/i )
{
  print "Use of model dielectric requires that epsilon be greater than 1!\n";
  print "Input diemac too low: $epsilon\n";
  $epsilon = 1.000001;
  open EPS, ">", "epsilon" or die "Failed to open epsilon\n$!";
  print EPS "$epsilon\n";
  close EPS;
}



exit;

sub findval
{
  my $target = $_[0];
  my $out = 1000;
  my @val_list = ( 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 30, 32, 36, 40, 45, 48, 50, 
                   54, 60, 64, 72, 75, 80, 81, 90, 100, 1000 );
  foreach (@val_list)
  {
    if( $_ >= $target )
    {
      $out = $_;
      last;
    }
  }
  return $out;
}

sub isAllowed {
  my $n = $_[0];

  for( my $i=2; $i <= 7; $i++ )
  {
    while(  $n % $i == 0 )
    {
      $n /= $i;
    }
  }
  if( $n == 1 || $n == 11 || $n == 13 )
  {
    return 1;
  }
  return 0;
}

