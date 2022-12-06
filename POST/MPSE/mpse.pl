#!/usr/bin/perl

use strict;

require JSON::PP;
use JSON::PP;
use File::Copy;
use Cwd 'abs_path';
use File::Spec::Functions;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/mpse\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
###########################


my $dataFile = catfile( updir(), updir(), "Common", "postDefaultsOceanDatafile" );

my $json = JSON::PP->new;
$json->canonical([1]);
$json->pretty([1]);

my $commonOceanData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $commonOceanData = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $dataFile\n$!";
}

my $gam = 0;
my $NIter = 1;
my $eps0;

open IN, "../opcons" or die;
open OUT, ">", "loss.dat" or die;
open OSC, ">", "osc_str.dat" or die;
open EPS2, ">", "eps2" or die;
my $line = <IN>;
while( $line =~ m/\s*#/ ) {
  $line = <IN>;
}
$line =~ m/(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)/ or die;
$eps0 = $2;
print OUT "$1  $8\n";
printf OSC " %g %g\n", $1, $1*$3;
printf EPS2 " %g %g\n", $1, $3;
while( my $line = <IN> ) {
  $line =~ m/(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(-?\d\.\d+E[+-]\d+)/ or die;
  print OUT "$1  $8\n";
  printf OSC " %g %g\n", $1, $1*$3;
  printf EPS2 " %g %g\n", $1, $3;
}
close IN;
close OUT;
close OSC;
close EPS2;

open IN, "../gap" or die;
my $ldagap = <IN>;
chomp $ldagap;
close IN;
my $egap = $ldagap;
my $efermi = $egap/2;

if( $commonOceanData->{"structure"}->{"metal"} ) {
  $eps0 = -1;
}

my $npoles = 50;

my $rs = ($commonOceanData->{"structure"}->{"valence_electrons"} / $commonOceanData->{"structure"}->{"volume"} )
       * 4.18879020478639084;
$rs = $rs**(-1.0/3.0);

open OUT, ">", "exc.inp" or die;
printf OUT "%i\n1 1\n%g\n", $npoles, $eps0;
close OUT;

  
for( my $i = 0; $i < $NIter; $i ++ ) {
  if( $i > 0 ) {
    copy "loss_SE", "loss.dat";
#    $eps0=`grep -v '#' NBSE/eps1_SE |awk '{if(NF==2) print $2;}' |head -n 1`
#    $egap=`head -n6 SE.dat |tail -n1 |awk "{print $EGap0 + \\$3}"`
#    $efermi = $egap/2;
  }

  system("$ENV{'OCEAN_BIN'}/eps2exc.x < exc.inp");

  my $excdat;
  if( open IN, "exc.dat" )
  {
    local $/ = undef;
    $excdat = <IN>;
    close IN;
  }
  open OUT, ">", "exc.dat" or die;
  print OUT "$rs 0.0\n$npoles\n$excdat";
  close OUT;

  open OUT, ">", "seconv.inp" or die;
  print OUT "# File to convolve.\nosc_str.dat\n# Output file.\nconv.dat\n# number of y columns\n1\n"
          . "# x column\n1\n# y columns to convolve\n2\n$efermi $egap $gam\n";
  close OUT;

  system("$ENV{'OCEAN_BIN'}/selfenergy.x");
  system("$ENV{'OCEAN_BIN'}/kkconv.x");

}
