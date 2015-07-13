# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
#!/usr/bin/perl

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
print "$ncpus has ", $#cpu_factors + 1, " factors\n";
foreach (@cpu_factors)
{ print "$_\n"; }

open QE_POOL, ">pool_control" or die "Failed to open qe_pool_control for writing\n$!\n";


open RSCALE, "rscale" or die;
open RPRIM, "rprim" or die;
<RSCALE> =~  m/(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)/ or die;
my @rscale = ($1, $3, $5);
print "$1\t$3\t$5\n";
close RSCALE;

open AVECS, ">avecsinbohr.ipt" or die;
for (my $i = 0; $i < 3; $i++ ) {
  <RPRIM> =~  m/([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)/ or die "$_";
  print AVECS $1*$rscale[0] . "  " . $3*$rscale[1] .  "  " . $5*$rscale[2] . "\n";
  print "$1\t$3\t$5\n";
}
close RPRIM;
close AVECS;

my @alength;
my @avec;
open AVECS, "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!\n";
<AVECS> =~ m/(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/ or die "$_\n";
$alength[0] = sqrt( $1*$1 + $2*$2 + $3*$3 );
$avec[0][0] = $1; $avec[1][0] = $2; $avec[2][0] = $3;
<AVECS> =~ m/(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/ or die "$_\n";
$alength[1] = sqrt( $1*$1 + $2*$2 + $3*$3 );
$avec[0][1] = $1; $avec[1][1] = $2; $avec[2][1] = $3;
<AVECS> =~ m/(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/ or die "$_\n";
$alength[2] = sqrt( $1*$1 + $2*$2 + $3*$3 );
$avec[0][2] = $1; $avec[1][2] = $2; $avec[2][2] = $3;
close AVECS;

my $volume = $avec[0][0] * ($avec[1][1] * $avec[2][2] - $avec[2][1] * $avec[1][2] )
           - $avec[1][0] * ($avec[0][1] * $avec[2][2] - $avec[2][1] * $avec[0][2] )
           + $avec[2][0] * ($avec[0][1] * $avec[1][2] - $avec[1][1] * $avec[0][2] );
print "Volume:\t$volume\n";

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
if( $input_content =~ m/^\s*-1/ )
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
  open INPUT, ">kmesh.ipt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for kmesh.ipt:\t$output[0]\t$output[1]\t$output[2]\n";
  $kpt_tot = $output[0]*$output[1]*$output[2];
} 
elsif ( $input_content =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)/ )
{
  $kpt_tot = $1*$2*$3;
}
else
{
  die "FAILED TO CORRECTLY PARSE nkpt\n";
}

my $ideal_npools = 1;
foreach (@cpu_factors)
{
  $ideal_npools = $_ unless( $kpt_tot % $_ );
}
print "NKPTS: $kpt_tot, ideal pools: $ideal_npools\n";
print QE_POOL "nscf\t$ideal_npools\n";
# Not integrated cleanly right now
# Need to control pools for the interpolation routines
my $min_interp_pool_size = int( $volume / 1200 );
print "Minimum size for interpolation pools $min_interp_pool_size\n";

my $pool_size = -1;
foreach( @cpu_factors )
{
  if( $_ >= $min_interp_pool_size && ($ncpus/$_) <= $kpt_tot )
  {
     $pool_size = $_;
     last;
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
open INPUT, "xmesh.ipt" or die "Failed to open xmesh.ipt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/^\s*-1/ )
{
  my @output;
  print "Defaults requested for xmesh.ipt\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    $target = $defaults{'xmesh'}[$acc_level];
    $target = $target * $alength[$i];
    $output[$i] = findval($target);
    if( $output[$i] == 1000 ) { die "kmesh too large\n"; }
  }
  open INPUT, ">xmesh.ipt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for xmesh.ipt:\t$output[0]\t$output[1]\t$output[2]\n";
}

$input_content = '';
open INPUT, "ngkpt" or die "Failed to open ngkpt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/^\s*-1/ )
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
}
elsif ( $input_content =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)/ )
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
die "Either dft_energy_range or nbands must be speified\n" 
    if( $erange <= 0 && $nbands <= 0);
if( $nbands <= 0 )
{
  print "Default requested for nbands. Energy range is $erange eV.\n";
  $erange = $erange / 13.605;
  $nbands = 0.01688686394 * $volume * ( $erange**(3/2) );
  $nbands *= 1.05;
  $nbands = int($nbands);
  # 1.05 is a padding factor
  print "Default chosen for nbands:\t$nbands\n";
  open NBANDS, ">nbands" or die "Failed to open nbands for writing.\n$!";
  print NBANDS "$nbands\n";
  close NBANDS;
}

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
open INPUT, "paw.nkpt" or die "Failed to open paw.nkpt\n$!\n";
while (<INPUT>)
  {
    $input_content .= $_;
  }
$input_content =~ s/\n/ /g;
close INPUT;
if( $input_content =~ m/^\s*-1/ )
{
  my @output;
  print "Defaults requested for paw.nkpt\n";
  for( my $i = 0; $i < 3; $i++ )
  {
    $target = $defaults{'paw'}[$acc_level];
    $target = $target / $alength[$i];
    $output[$i] = findval($target);
    if( $output[$i] == 1000 ) { die "paw.nkpt too large\n"; }
  }
  open INPUT, ">paw.nkpt" or die "$!\n";
  print INPUT "$output[0]  $output[1]  $output[2]\n";
  close INPUT;
  print "Defaults chosen for paw.nkpt:\t$output[0]\t$output[1]\t$output[2]\n";
  $kpt_tot = $output[0]*$output[1]*$output[2];
}
elsif ( $input_content =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)/ )
{
  $kpt_tot = $1*$2*$3;
}
else
{
  die "FAILED TO CORRECTLY PARSE paw.nkpt\n";
}

my $ideal_npools = 1;
foreach (@cpu_factors)
{
  $ideal_npools = $_ unless( $kpt_tot %  $_ );
}
print "PAW: $kpt_tot, ideal pools: $ideal_npools\n";
print QE_POOL "paw\t$ideal_npools\n";

# Not integrated cleanly right now
# Need to control pools for the interpolation routines
$min_interp_pool_size = int( $volume / 1200 );
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
print QE_POOL "interpolate paw\t$pool_size\n";


close QE_POOL;
exit;

sub findval
{
  my $target = $_[0];
  my $out = 1000;
  my @val_list = ( 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 30, 32, 36, 40, 45, 48, 50, 100, 1000 );
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
