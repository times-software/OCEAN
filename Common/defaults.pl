#!/usr/bin/perl
# Copyright (C) 2015-2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
require JSON::PP;
#JSON::PP->import;
use JSON::PP;
use POSIX;
use Math::BigFloat;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/defaults\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }


print `pwd`;


my $dataFile = "oceanDatafile";
if( -e $dataFile )
{
  my $json = JSON::PP->new;
  my $oceanData;
  if( open( my $in, "<", $dataFile ))
  {
    local $/ = undef;
    $oceanData = $json->decode(<$in>);
    close($in);
  }
  else
  {
    die "Failed to open config file $dataFile\n$!";
  }

  $oceanData->{'warnings'} = {} unless( exists $oceanData->{'warnings'} );
  $oceanData->{'warnings'}->{'defaults.pl'} = [] unless( exists $oceanData->{'warnings'}->{'defaults.pl'} );

  $oceanData->{'computer'}->{'ncpus'} = parse_para_prefix( $oceanData );
  add_cpuFactorsAndSquare( $oceanData );
  abVecsVolume( $oceanData ); 

  checkKpoints( $oceanData, $oceanData->{'dft'}->{'den'}->{'kmesh'} );
  checkKpoints( $oceanData, $oceanData->{'bse'}->{'kmesh'} );
  checkKpoints( $oceanData, $oceanData->{'screen'}->{'kmesh'} );

  checkBands( $oceanData, $oceanData->{'screen'}, $oceanData->{'screen'}->{'dft_energy_range'} );
  checkBands( $oceanData, $oceanData->{'bse'}, $oceanData->{'bse'}->{'dft_energy_range'} );

  checkXpoints( $oceanData, $oceanData->{'bse'}->{'xmesh'} );

  checkQ( $oceanData );

  checkEpsilon( $oceanData );

  checkDFTconvergence( $oceanData );

  checkDFToccopt( $oceanData );

  fixSerPrefix( $oceanData );

  makeEdges( $oceanData );
  fixCoords( $oceanData );
  makeSitesZsymb( $oceanData );
  ## xred sub to go from bohr/ang to red

  fixCNBSE( $oceanData );
  
  my $enable = 1;
  $json->canonical([$enable]);
  $json->pretty([$enable]);
  open OUT, ">", "postDefaultsOceanDatafile" or die;
  print OUT $json->encode($oceanData);
  close OUT;

  exit 0;
}
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
  $input_content =~ m/^\s*(-?\d+\.?\d*)/ or die "Failed to parse xmesh.ipt\n$input_content\n";
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

sub parse_para_prefix( ) 
{
  my $para_prefix;
  my $ncpus = 1;
  if( scalar @_ > 0 )
  {
    $para_prefix = $_[0]->{'computer'}->{'para_prefix'};
  }
  else
  {
    open PARA_PREFIX, "para_prefix" or die "$!";
    $para_prefix = <PARA_PREFIX>;
    if( $para_prefix =~ m/\#/ )
    {
      close PARA_PREFIX;
      open PARA_PREFIX, ">para_prefix" or die "$!";
      print PARA_PREFIX "";
    }
    close PARA_PREFIX;
  }
  if ( $para_prefix =~ m/^\D*(\d+)/ )
  {
    $ncpus = $1; 
    print "GRABBED $ncpus as the number of cpus to run with!\n";
  } 
  return $ncpus;
}

sub add_cpuFactorsAndSquare( )
{
  my $hashRef = $_[0];
  my $ncpus = $hashRef->{'computer'}->{'ncpus'};
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
  $hashRef->{'computer'}->{'cpu_factors'} = [@cpu_factors];

  my @square_cpu_factors;
  foreach my $cpu_count (@cpu_factors)
  {
    push @square_cpu_factors, $cpu_count if( int(sqrt($cpu_count)) ** 2 == $cpu_count );
  }
  foreach (@square_cpu_factors)
  { print "$_\n"; }
  $hashRef->{'computer'}->{'cpu_square_factors'} = [@square_cpu_factors];
}

sub abVecsVolume
{
  my $hashRef = $_[0];
  my @avec; my @alen; my @bvec; my @blen;
  for( my $i = 0; $i < 3; $i++ )
  {
    for( my $j = 0; $j< 3; $j++ )
    {
      $avec[$i][$j] = $hashRef->{'structure'}->{'rprim'}[$i*3+$j] 
                    * $hashRef->{'structure'}->{'rscale'}[$i];
    }
#  print "$avec[$i][0]  $avec[$i][1]  $avec[$i][2]\n";
  }
  my $volume = $avec[0][0] * ($avec[1][1] * $avec[2][2] - $avec[2][1] * $avec[1][2] )
             - $avec[1][0] * ($avec[0][1] * $avec[2][2] - $avec[2][1] * $avec[0][2] )
             + $avec[2][0] * ($avec[0][1] * $avec[1][2] - $avec[1][1] * $avec[0][2] );
  print "Volume:\t$volume\n";

  
  my $pref = 2*4*atan2(1,1)/$volume;

  for (my $i = 0; $i < 3; $i++ ) {
    for (my $j = 0; $j < 3; $j++ ) {
      $bvec[$i][$j] = $pref * ($avec[$i-2][$j-2] * $avec[$i-1][$j-1] - $avec[$i-2][$j-1] * $avec[$i-1][$j-2]) ;
    }
    $blen[$i] = sqrt( $bvec[$i][0]**2 + $bvec[$i][1]**2 + $bvec[$i][2]**2 );
    $alen[$i] = sqrt( $avec[$i][0]**2 + $avec[$i][1]**2 + $avec[$i][2]**2 );
  }

  $hashRef->{'structure'}->{'volume'} = $volume;
  $hashRef->{'structure'}->{'avecs'} = [@avec];
  $hashRef->{'structure'}->{'alen'} = [@alen];
  $hashRef->{'structure'}->{'bvecs'} = [@bvec];
  $hashRef->{'structure'}->{'blen'} = [@blen];
}


sub checkKpoints
{
  my $hashRef = $_[0];
  my $kpointref = $_[1];

  for( my $i = 0; $i < scalar @{ $kpointref }; $i++ )
  {
    print ${$kpointref}[$i] . "\n";
  }

  return -1 if( scalar @{ $kpointref } < 1 );
  for( my $i = scalar @{ $kpointref }; $i < 3; $i++ )
  {
    ${$kpointref}[$i] = ${$kpointref}[$i-1];
  }

  for( my $i = 0; $i < 3; $i++ )
  {
    ${$kpointref}[$i] = int(${$kpointref}[$i]) if( ${$kpointref}[$i] > 0 );
    ${$kpointref}[$i] = -3 if( ${$kpointref}[$i] == 0 );
  }

  for( my $i = 0; $i < 3; $i++ )
  {
    if( ${$kpointref}[$i] < 0 )
    {
      ${$kpointref}[$i] = int( ${$kpointref}[$i] * -0.999 * $hashRef->{'structure'}->{'blen'}[$i] ) + 1;
      while( isAllowed( ${$kpointref}[$i] ) == 0 )
      {
        ${$kpointref}[$i]++;
      }
    }
  }
}

sub checkXpoints
{
  my $hashRef = $_[0];
  my $xpointref = $_[1];

  for( my $i = 0; $i < scalar @{ $xpointref }; $i++ )
  {
    print ${$xpointref}[$i] . "\n";
  }

  return -1 if( scalar @{ $xpointref } < 1 );
  for( my $i = scalar @{ $xpointref }; $i < 3; $i++ )
  {
    ${$xpointref}[$i] = ${$xpointref}[$i-1];
  }

  for( my $i = 0; $i < 3; $i++ )
  {
    ${$xpointref}[$i] = int(${$xpointref}[$i]) if( ${$xpointref}[$i] > 0 );
    ${$xpointref}[$i] = -.82 if( ${$xpointref}[$i] == 0 );
  }

  my $totalXmesh = 0;
  my $minXtotal = 2 * $hashRef->{'bse'}->{'nbands'};
  my $warnFlag = 1;
  
#  while( $totalXmesh < $minXtotal )
  do
  {
    $totalXmesh = 1;
    for( my $i = 0; $i < 3; $i++ )
    {
      if( ${$xpointref}[$i] < 0 )
      {
        ${$xpointref}[$i] = int( ${$xpointref}[$i] * -0.999 * $hashRef->{'structure'}->{'alen'}[$i] ) + 1;
        while( isAllowed( ${$xpointref}[$i] ) == 0 )
        {
          ${$xpointref}[$i]++;
        }
      }
      $totalXmesh *= ${$xpointref}[$i];
    }
    unless( $totalXmesh > $minXtotal )
    {
      if( $warnFlag != 0 )
      {
        print "Requested xmesh density too low!\n";
        print "  Need twice as many x-points as bands\n";
        $warnFlag = 0;
      }
      print "   " . ${$xpointref}[0] . " " . ${$xpointref}[1] . " " . ${$xpointref}[2] . "\n";
      my $currentDen = (${$xpointref}[0] + 1 )/ $hashRef->{'structure'}->{'alen'}[0];
      for( my $i = 1; $i < 3; $i++ )
      {
        if( $currentDen > (${$xpointref}[$i] + 1 )/ $hashRef->{'structure'}->{'alen'}[$i] )
        {
          $currentDen =  (${$xpointref}[$i] + 1 )/ $hashRef->{'structure'}->{'alen'}[$i];
        }
      }
      for( my $i = 0; $i < 3; $i++ )
      {
        ${$xpointref}[$i] = - $currentDen;
      }
    }
  } until( $totalXmesh > $minXtotal );
  print "   " . ${$xpointref}[0] . " " . ${$xpointref}[1] . " " . ${$xpointref}[2] . "\n";
}


sub checkBands
{
  my $hashRef = $_[0];
  my $bandRef = $_[1];
  my $energyRef = $_[2];
  my $nb = $bandRef->{'nbands'};
  print "$nb\n";

  # Do nothing if we have a positive number of bands
  if( $nb =~ m/^\s*(-?\d+)/s*$/ ) {
    return if( $nb > 0); 
  }

  # Not zero means add to total valence bands, continue but WARN if no electron count
  my $flag = 0;
  if( $nb < 0 )
  {
    $flag = 1;
    $nb = abs( $nb );
  }
  else
  {
    $nb = 0.019 * $hashRef->{'structure'}->{'volume'} * ( ($energyRef/13.605)**(3/2) );
    $nb = 1 if ($nb < 1 );
    print "$nb  $hashRef->{'structure'}->{'volume'} $energyRef\n";
  }

  if( $hashRef->{'structure'}->{'valence_electrons'} < 1 )
  {
    if( $flag ) {
      print "WARNING! A specific number of conduction bands was requested, " .
            "but electrons were not parsed from the pseudopotentials.\n" .
            "  Will continue using an estimate.\n";
    }
    $nb += 0.036 * $hashRef->{'structure'}->{'volume'};
  }
  else
  {
    $nb += floor( ($hashRef->{'structure'}->{'valence_electrons'}+1)/2);
  }
  $nb = floor( ( 2 * $nb + 1 ) / 2 );
  $bandRef->{'nbands'} = $nb;
  print "$nb\n";
}


sub checkQ
{
  my $hashRef = $_[0];
  my $minQ = 0.0000000001;
  if( $hashRef->{'calc'}->{'mode'} =~ m/VAL/i || $hashRef->{'calc'}->{'mode'} =~ m/RXS/i )
  {
    my $q = 0;
    for( my $i = 0; $i < 3; $i++ )
    {
      $q += abs( $hashRef->{'calc'}->{'photon_q'}[$i] );
    }
    if( $q < $minQ )
    {
      print "Valence or RIXS calculation requires non-zero q-vector\n";
      print "Assigning default q-vec of ( 0.001, 0, 0 ). You should edit your input\n";
      $hashRef->{'calc'}->{'photon_q'}[0] = 0.001;
      $hashRef->{'calc'}->{'photon_q'}[1] = 0.0;
      $hashRef->{'calc'}->{'photon_q'}[2] = 0.0;
    }
  }
  
  my $q = 0;
  for( my $i = 0; $i < 3; $i++ )
  {
    $q += abs( $hashRef->{'calc'}->{'photon_q'}[$i] );
  }
  if( $q < $minQ ) {
    $hashRef->{'calc'}->{'nonzero_q'} = JSON::PP::false;
  } else {
    $hashRef->{'calc'}->{'nonzero_q'} = JSON::PP::true;
  }
  $hashRef->{'dft'}->{'bse'}->{'split'} = JSON::PP::false unless( $hashRef->{'calc'}->{'nonzero_q'} );
}

sub checkEpsilon
{
  my $hashRef = $_[0];
  if( $hashRef->{'dft'}->{'epsilon'}->{'method'} =~ m/input/ )
  {
    if( $hashRef->{'structure'}->{'epsilon'} < 1.000001 )
    {
      print "Use of model dielectric requires that epsilon be greater than 1!\n";
      print "Input diemac too low: " . $hashRef->{'structure'}->{'epsilon'} . "\n";
      $hashRef->{'structure'}->{'epsilon'} = 1.000001;
    }
  }
}

sub checkDFTconvergence
{
  my $hashRef = $_[0];
  
  $hashRef->{'dft'}->{'den'}->{'toldfe'} = $hashRef->{'dft'}->{'toldfe'} 
      if( $hashRef->{'dft'}->{'den'}->{'toldfe'} <= 0 );
  $hashRef->{'dft'}->{'bse'}->{'toldfe'} = $hashRef->{'dft'}->{'tolwfr'}    
      if( $hashRef->{'dft'}->{'bse'}->{'toldfe'} <= 0 );
  $hashRef->{'dft'}->{'screen'}->{'toldfe'} = $hashRef->{'dft'}->{'tolwfr'}    
      if( $hashRef->{'dft'}->{'screen'}->{'toldfe'} <= 0 );

  delete $hashRef->{'dft'}->{'toldfe'};
  delete $hashRef->{'dft'}->{'tolwfr'};
}

sub checkDFToccopt
{
  my $hashRef = $_[0];

  if( $hashRef->{'structure'}->{'valence_electrons'} > 0 ) {
    if( $hashRef->{'structure'}->{'valence_electrons'} % 2 ) {
      $hashRef->{'structure'}->{'metal'} = JSON::PP::true unless( defined( $hashRef->{'structure'}->{'metal'} ) );
      unless( $hashRef->{'structure'}->{'metal'} ) {
        print "WARN!!\n"
            . "System was not set to metallic, but it has an odd number of electrons\n";
        my $s = "Metal flag changed from false to true";
        print "  $s\n";
        $hashRef->{'warnings'} = {} unless( exists $hashRef->{'warnings'} );
        $hashRef->{'warnings'}->{'defaults.pl'} = [] unless( exists $hashRef->{'warnings'}->{'defaults.pl'} );
        push @{$hashRef->{'warnings'}->{'defaults.pl'}}, $s;
      }
    }
  }

  unless( defined( $hashRef->{'structure'}->{'metal'} ) ) {
    if( $hashRef->{'dft'}->{'occopt'} == 1 ) {
      $hashRef->{'structure'}->{'metal'} = JSON::PP::false;
    } else {
      $hashRef->{'structure'}->{'metal'} = JSON::PP::true;
    }
  }

  if( $hashRef->{'structure'}->{'metal'} ) {
    if( $hashRef->{'dft'}->{'occopt'} == 1 ) {
      print "WARN!!\n  System was set to metallic, but occopt was set to 1\n" .
                    "  Consider adjusting both occopt and degauss to match your needs.\n" .
                    "  This affects the DFT SCF calculation\n";
      $hashRef->{'warnings'} = {} unless( exists $hashRef->{'warnings'} );
      $hashRef->{'warnings'}->{'defaults.pl'} = [] unless( exists $hashRef->{'warnings'}->{'defaults.pl'} );
      my $s = "Occopt changed from 1 to 3";
      $hashRef->{'dft'}->{'occopt'} = 3;
      push @{$hashRef->{'warnings'}->{'defaults.pl'}}, $s;
      print "$s\n";
    }
  } else {
    if( $hashRef->{'dft'}->{'occopt'} != 1 ) {
      print "WARN!!\n  System was set to insulating, but occopt was set to not be 1\n";
      my $s = sprintf "Occopt changed from %i to 1", $hashRef->{'dft'}->{'occopt'};
      $hashRef->{'dft'}->{'occopt'} = 1;
      $hashRef->{'warnings'} = {} unless( exists $hashRef->{'warnings'} );
      $hashRef->{'warnings'}->{'defaults.pl'} = [] unless( exists $hashRef->{'warnings'}->{'defaults.pl'} );
      push @{$hashRef->{'warnings'}->{'defaults.pl'}}, $s;
      print "$s\n";
    }
  }

}

sub fixSerPrefix
{
  my $hashRef = $_[0];
  if( $hashRef->{'computer'}->{'ser_prefix'} eq "!" )
  {
    $hashRef->{'computer'}->{'ser_prefix'} = $hashRef->{'computer'}->{'para_prefix'};
    $hashRef->{'computer'}->{'ser_prefix'} =~ s/\d+/1/;
  }
}

sub makeEdges
{
  my $hashRef = $_[0];
  return if( $hashRef->{'calc'}->{'mode'} =~ m/val/i );


  unless( scalar @{$hashRef->{'calc'}->{'edges'}} % 3 == 0 )
  {
    die "Malformed calc.edges input\n";
  } 

  my %SitesByZ;
  for( my $j = 0; $j < scalar @{$hashRef->{'structure'}->{'typat'}}; $j++ )
  {
    my $i = $hashRef->{'structure'}->{'typat'}[$j] - 1;
    my $z = $hashRef->{'structure'}->{'znucl'}[$i];
    push @{ $SitesByZ{ $z } }, $j+1;
  }

  foreach my $key ( keys %SitesByZ )
  {
    print $key . ":\n";
    foreach my $site ( @{ $SitesByZ{ $key } } )
    {
      print "  $site";
    }
    print "\n";
  }

  foreach (@{$hashRef->{'calc'}->{'edges'}} )
  {
    print "$_\n";
  }


  # Create a hash of Z
  #   which stores hash NL
  #     which has array of sites

  # When writing out edges we want to group by this, and sort sites into numerical order

  #TODO sort on Z as well
  my %FullEdges;

  for( my $i = 0; $i < scalar @{$hashRef->{'calc'}->{'edges'}}; $i+=3 )
  {
    my $Z = $hashRef->{'calc'}->{'edges'}[$i];
    my $N = $hashRef->{'calc'}->{'edges'}[$i+1];
    my $L = $hashRef->{'calc'}->{'edges'}[$i+2];

    if( $Z == 0 )
    {
      die "Illegal Z value in edges:\t$Z\n";
    }
    if( $N < 1 || $N > 7 )
    {
      die "Illegal principle qunatum number:\t$N\n";
    }
    if( $L < 0 || $L > 4 )
    {
      die "Illegal angular quantum number:\t$L\n";
    }

    # Are we selecting all Z?
    if( $Z < 0 )
    {
      my $ZEE = abs( $Z );
      unless( exists $SitesByZ{ $ZEE } )
      {
        die "Requested Z in edges that is not present in znucl/typat:\t$Z\t$ZEE\n";
      }
      foreach my $site ( @{ $SitesByZ{ $ZEE } } )
      {
        $FullEdges{ $ZEE }{ "$N $L" }{ $site } = 1;
      }
    }
    else
    {
      # Z is the site number
      if ( $Z > scalar @{ $hashRef->{'structure'}->{'typat'}}  )
      {
        die "Site requested in edges is larger than the number of site:\t$Z\n";
      }
      my $ZEE = $hashRef->{'structure'}->{'typat'}[ $Z - 1 ];
      my $site = $Z;
      $FullEdges{ $ZEE }{ "$N $L" }{ $site } = 1;
    }
  }

  my @edges;
  foreach my $ZEE ( keys %FullEdges )
  {
    foreach my $NL ( keys %{ $FullEdges{ $ZEE } } )
    {
      my @tempsite;
      foreach my $site ( keys %{ $FullEdges{ $ZEE }{ $NL } } )
      {
        push @tempsite, $site;
      }

      my @sortsite = sort { $a <=> $b } @tempsite;
      foreach my $site ( @sortsite )
      {
        push @edges, "$site $NL";
        if( 0 ) {
          # multiply by 1 to force number instead of string context for JSON
          my $Z = $ZEE * 1;
          $NL =~ m/(\d)\s(\d)/;
          my $N = $1 * 1;
          my $L = $2 * 1;
          my @tmp = ( $site, $Z, $N, $L );
          push @edges, \@tmp;
        }
      }
    }
  }

  $hashRef->{'calc'}->{'edges'} = \@edges;
    
}


sub invertAvecs
{
  my $hashRef = $_[0];
  my @a = @{$hashRef->{'structure'}->{'avecs'}};

  my $det = $a[0][0] * ( $a[1][1] * $a[2][2] - $a[2][1] * $a[1][2] )
          - $a[0][1] * ( $a[1][0] * $a[2][2] - $a[1][2] * $a[2][0] )
          + $a[0][2] * ( $a[1][0] * $a[2][1] - $a[1][1] * $a[2][0] );


  my @invA;
  $invA[0][0] = ( $a[1][1] * $a[2][2] - $a[2][1] * $a[1][2] ) / $det;
  $invA[0][1] = ( $a[0][2] * $a[2][1] - $a[0][1] * $a[2][2] ) / $det;
  $invA[0][2] = ( $a[0][1] * $a[1][2] - $a[0][2] * $a[1][1] ) / $det;

  $invA[1][0] = ( $a[1][2] * $a[2][0] - $a[1][0] * $a[2][2] ) / $det;
  $invA[1][1] = ( $a[0][0] * $a[2][2] - $a[0][2] * $a[2][0] ) / $det;
  $invA[1][2] = ( $a[1][0] * $a[0][2] - $a[0][0] * $a[1][2] ) / $det;

  $invA[2][0] = ( $a[1][0] * $a[2][1] - $a[2][0] * $a[1][1] ) / $det;
  $invA[2][1] = ( $a[2][0] * $a[0][1] - $a[0][0] * $a[2][1] ) / $det;
  $invA[2][2] = ( $a[0][0] * $a[1][1] - $a[1][0] * $a[0][1] ) / $det;

  $hashRef->{'structure'}->{'invA'} = [ @invA ];

#  my $n = 30;
#  Math::BigFloat->accuracy($n);

#  my $x00 = Math::BigFloat->new($a[0][0]);
#  my $x01 = Math::BigFloat->new($a[0][1]);
#  my $x02 = Math::BigFloat->new($a[0][2]);
#  my $x10 = Math::BigFloat->new($a[1][0]);
#  my $x11 = Math::BigFloat->new($a[1][1]);
#  my $x12 = Math::BigFloat->new($a[1][2]);
#  my $x20 = Math::BigFloat->new($a[2][0]);
#  my $x21 = Math::BigFloat->new($a[2][1]);
#  my $x22 = Math::BigFloat->new($a[2][2]);
  
#  my $det = $x00->copy();
}

sub fixCoords
{
  my $hashRef = $_[0];

  my $bohr = 0.529177210903; #CODATA 2018

  invertAvecs( $hashRef );

  my @newReduced;
  my @newRealspace;

  my $tauLength = scalar @{$hashRef->{'structure'}->{'xred'}};
  my $realLength = scalar @{$hashRef->{'structure'}->{'xbohr'}};
  my $angLength = scalar @{$hashRef->{'structure'}->{'xangst'}};

  if( $tauLength > 0 )
  {
    print "CAUTION: Both xred and xbohr were supplied, but xred will be used!\n" if( $realLength > 0 );
    print "CAUTION: Both xred and xangst were supplied, but xred will be used!\n" if( $angLength > 0 );
    if( $tauLength % 3 != 0 )
    {
      die "Malformed xred in input\n";
    }

    $tauLength /= 3;
    for( my $i = 0; $i < $tauLength; $i++ )
    {
      $newReduced[$i][0] = $hashRef->{'structure'}->{'xred'}[3*$i];
      $newReduced[$i][1] = $hashRef->{'structure'}->{'xred'}[3*$i+1];
      $newReduced[$i][2] = $hashRef->{'structure'}->{'xred'}[3*$i+2];
    }

    for( my $i = 0; $i < $tauLength; $i++ )
    {
      $newRealspace[$i][0] = $hashRef->{'structure'}->{'avecs'}[0][0] * $newReduced[$i][0] 
                           + $hashRef->{'structure'}->{'avecs'}[0][1] * $newReduced[$i][1]
                           + $hashRef->{'structure'}->{'avecs'}[0][2] * $newReduced[$i][2];

      $newRealspace[$i][1] = $hashRef->{'structure'}->{'avecs'}[1][0] * $newReduced[$i][0]
                           + $hashRef->{'structure'}->{'avecs'}[1][1] * $newReduced[$i][1]
                           + $hashRef->{'structure'}->{'avecs'}[1][2] * $newReduced[$i][2];

      $newRealspace[$i][2] = $hashRef->{'structure'}->{'avecs'}[2][0] * $newReduced[$i][0]
                           + $hashRef->{'structure'}->{'avecs'}[2][1] * $newReduced[$i][1]
                           + $hashRef->{'structure'}->{'avecs'}[2][2] * $newReduced[$i][2];
    }
  }
  else  
  {
    if( $angLength > 0 )
    {
      if( $realLength > 0 )
      {
        print "CAUTION: Both xangst and xbohr were supplied, but xbohr will be used!\n";
      }
      else
      {
        #COPY angst over with scaling factor
        for( my $i = 0; $i < $angLength; $i++ )
        {
          $hashRef->{'structure'}->{'xbohr'}[$i] = $hashRef->{'structure'}->{'xangst'}[$1] / $bohr;
        }
      }
    }
    

    if( $realLength % 3 != 0 )
    {
      die "Malformed xbohr in input\n";
    }

    $realLength /= 3;

    for( my $i = 0; $i < $tauLength; $i++ )
    { 
      $newRealspace[$i][0] = $hashRef->{'structure'}->{'xbohr'}[3*$i];
      $newRealspace[$i][1] = $hashRef->{'structure'}->{'xbohr'}[3*$i+1];
      $newRealspace[$i][2] = $hashRef->{'structure'}->{'xbohr'}[3*$i+2];
    }
    
    for( my $i = 0; $i < $tauLength; $i++ )
    { 
      $newReduced[$i][0] = $hashRef->{'structure'}->{'invA'}[0][0] * $newRealspace[$i][0]
                         + $hashRef->{'structure'}->{'invA'}[0][1] * $newRealspace[$i][1]
                         + $hashRef->{'structure'}->{'invA'}[0][2] * $newRealspace[$i][2];
      
      $newReduced[$i][1] = $hashRef->{'structure'}->{'invA'}[1][0] * $newRealspace[$i][0]
                         + $hashRef->{'structure'}->{'invA'}[1][1] * $newRealspace[$i][1]
                         + $hashRef->{'structure'}->{'invA'}[1][2] * $newRealspace[$i][2];
      
      $newReduced[$i][2] = $hashRef->{'structure'}->{'invA'}[2][0] * $newRealspace[$i][0]
                         + $hashRef->{'structure'}->{'invA'}[2][1] * $newRealspace[$i][1]
                         + $hashRef->{'structure'}->{'invA'}[2][2] * $newRealspace[$i][2];
    }


  }

  $hashRef->{'structure'}->{'xred'} = [ @newReduced ];
  $hashRef->{'structure'}->{'xbohr'} = [ @newRealspace ];


  my @a;
  for( my $i = 0; $i < scalar @newRealspace; $i++ )
  {
    $a[$i][0] = $newRealspace[$i][0] * $bohr;
    $a[$i][1] = $newRealspace[$i][1] * $bohr;
    $a[$i][2] = $newRealspace[$i][2] * $bohr;
  }

  $hashRef->{'structure'}->{'xangst'} = [ @a ];

}


sub makeSitesZsymb
{
  my $hashRef = $_[0];

  # First make fancy site list for later

  my @z2symb =          ( '', 'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' ,                  
      'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar',
      'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
      'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' , 'Zr',
      'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
      'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',      
      'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
      'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
      'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu',
      'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',
      'Sg', 'Bh', 'Hs', 'Mt' );

#  $hashRef->{'structure'}->{'sites'} = [];
#  my %Z;
#  foreach( @{$hashRef->{'structure'}->{'typat'}} )
#  {
#    my $Z = $hashRef->{'structure'}->{'znucl'}[$_-1];
#    if( exists $Z{$Z} )
#    {
#      $Z{$Z}++;
#    } else {
#      $Z{$Z} = 1;
#    }
#    my $s = $z2symb[$Z] . '_';
#    push @{$hashRef->{'structure'}->{'sites'}}, sprintf "%.2s  %i", $s, $Z{$Z};
#  }
  

  $hashRef->{'structure'}->{'elname'} = [];
  foreach( @{$hashRef->{'structure'}->{'znucl'}} )
  {
    my $s = sprintf "%.2s", $z2symb[$_] . '_';
    push @{$hashRef->{'structure'}->{'elname'}}, $s;
  }


  # If zsymb is the correct length, assume it is fine
  return if( scalar @{$hashRef->{'structure'}->{'zsymb'}} == scalar @{$hashRef->{'structure'}->{'znucl'}} );
  # Only add number after symbol if required
  my %zsymbTracker;
  for( my $i = 0; $i < scalar @{$hashRef->{'structure'}->{'znucl'}}; $i++ )
  {
    if( exists $zsymbTracker{ $hashRef->{'structure'}->{'znucl'}[$i] } )
    {
      $zsymbTracker{ $hashRef->{'structure'}->{'znucl'}[$i] } = 1;
    }
    else
    {
      $zsymbTracker{ $hashRef->{'structure'}->{'znucl'}[$i] } = 0;
    }
  }

  for( my $i = 0; $i < scalar @{$hashRef->{'structure'}->{'znucl'}}; $i++ )
  {
    my $symb = $z2symb[ $hashRef->{'structure'}->{'znucl'}[$i] ];
    if( $zsymbTracker{ $hashRef->{'structure'}->{'znucl'}[$i] } > 0 )
    {
      $symb .= sprintf "%i", $zsymbTracker{ $hashRef->{'structure'}->{'znucl'}[$i] };
      $zsymbTracker{ $hashRef->{'structure'}->{'znucl'}[$i] }++;
    }
    $hashRef->{'structure'}->{'zsymb'}[$i] = $symb;
  }

}

sub fixCNBSE
{
  my $hashRef = $_[0];

  if( $hashRef->{'bse'}->{'core'}->{'strength'} < 0 ) {
    if( $hashRef->{'calc'}->{'mode'} eq 'xes' ) {
      $hashRef->{'bse'}->{'core'}->{'strength'} = 0;
    } else {
      $hashRef->{'bse'}->{'core'}->{'strength'} = 1;
    }
  }

  if( $hashRef->{'bse'}->{'core'}->{'screen_radius'} < 0 ) {
    if( $hashRef->{'screen'}->{'shells'}[-1] > 0 ) {
      $hashRef->{'bse'}->{'core'}->{'screen_radius'} = $hashRef->{'screen'}->{'shells'}[-1];
    } else {
      $hashRef->{'bse'}->{'core'}->{'screen_radius'} = abs( $hashRef->{'bse'}->{'core'}->{'screen_radius'} );
    }
  }

  my @tmp;
  my $found = 0;
  foreach my $r (@{$hashRef->{'screen'}->{'shells'}}) {
    if( abs( $hashRef->{'bse'}->{'core'}->{'screen_radius'} - $r ) < 0.01 ) {
      $found = 1;
    }
    push @tmp, $r if( $r > 0 );
  }

  $hashRef->{'screen'}->{'shells'} = \@tmp;

  if( $found == 0 ) {
    push @{$hashRef->{'screen'}->{'shells'}}, $hashRef->{'bse'}->{'core'}->{'screen_radius'};
  }
}
