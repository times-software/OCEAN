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
open PARA_PREFIX, "para_prefix" or die "$!";
my $para_prefix = <PARA_PREFIX>;
if( $para_prefix =~ m/\#/ )
{
  close PARA_PREFIX;
  open PARA_PREFIX, ">para_prefix" or die "$!";
  print PARA_PREFIX "";
}
close PARA_PREFIX;

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
  print INPUT "$output[0]\t$output[1]\t$output[2]\n";
  close INPUT;
  print "Defaults chosen for kmesh.ipt:\t$output[0]\t$output[1]\t$output[2]\n";
}
  
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
  print INPUT "$output[0]\t$output[1]\t$output[2]\n";
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
  print INPUT "$output[0]\t$output[1]\t$output[2]\n";
  close INPUT;
  print "Defaults chosen for ngkpt:\t$output[0]\t$output[1]\t$output[2]\n";
}

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
  print INPUT "$output[0]\t$output[1]\t$output[2]\n";
  close INPUT;
  print "Defaults chosen for obkpt:\t$output[0]\t$output[1]\t$output[2]\n";
}


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
