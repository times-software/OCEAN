#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

use POSIX;
require JSON::PP;
use JSON::PP;

use Cwd 'abs_path';
use Cwd;
use File::Spec::Functions;
use File::Copy;

use FindBin;
use lib $FindBin::Bin;
require 'OCEANcompare.pl';

use File::Basename;
use Time::HiRes qw( gettimeofday tv_interval );

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/cnbse\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}

my @ExtraFiles = ("Pquadrature", "sphpts" );
  
my $json = JSON::PP->new;
my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);

# Load run info from Common
my $dataFile = catfile( updir(), "Common", "postDefaultsOceanDatafile");
die "Failed to find $dataFile\n" unless( -e $dataFile );

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

# Load run info from PREP
$dataFile = catfile( updir(), "PREP", "prep.json" );
die "Failed to find $dataFile\n" unless( -e $dataFile );

my $prepData; 
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $prepData = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $dataFile\n$!";
}

my $haveScreen = 1;
$haveScreen = 0 if( $commonOceanData->{'calc'}->{'mode'} eq 'val' &&
                    $commonOceanData->{'screen'}->{'mode'} ne 'grid' );
my $screenData;
if( $haveScreen ) {
  # Load run info for SCREEN
  $dataFile = catfile( updir(), "SCREEN", "screen.json" );
  die "Failed to find $dataFile\n" unless( -e $dataFile );

  if( open( my $in, "<", $dataFile ))
  {
    local $/ = undef;
    $screenData = $json->decode(<$in>);
    close($in);
  }
  else
  {
    die "Failed to open config file $dataFile\n$!";
  }
}


# Grab previous run info if it exists
my $bseDataFile = "bse.json";
my $bseData;
if( -e $bseDataFile && open( my $in, "<", $bseDataFile ) )
{ 
  local $/ = undef;
  $bseData = $json->decode(<$in>);
  close($in);
}


my $newBSEdata;
$newBSEdata->{'bse'}={};
$newBSEdata->{'structure'} = {};
$newBSEdata->{'time'} = {};

$newBSEdata->{'time'}->{'dft'} = $prepData->{'dft_time'};
$newBSEdata->{'time'}->{'prep'} = $prepData->{'time'};
if( $haveScreen ) {
  $newBSEdata->{'time'}->{'screen'} = $screenData->{'time'};
} else {
  $newBSEdata->{'time'}->{'screen'} = 0;
}

$newBSEdata->{'bse'}->{'complete'} = JSON::PP::true;

# (was previous run done)
$newBSEdata->{'bse'}->{'complete'} = JSON::PP::false
    unless( exists $bseData->{'bse'} && exists $bseData->{'bse'}->{'complete'} 
            && $bseData->{'bse'}->{'complete'} );


# Has the structure been updated? But we want this from DFT?
#copyAndCompare( $newBSEdata, $commonOceanData, $bseData,
#                $newBSEdata->{'bse'}, ["structure"] );

my @list = ( "kmesh", "kshift", "xmesh", "nspin", "photon_q", "brange", "xred", "typat", 
             "znucl", "elname", "nelec", "avecs", "bvecs", "fermi", "epsilon", "cbm", "max_occ_band", "min_unocc_band" );
copyAndCompare( $newBSEdata->{'bse'}, $prepData->{'bse'}, $bseData->{'bse'},
                $newBSEdata->{'bse'}, \@list );


copyAndCompare( $newBSEdata->{'bse'}, $commonOceanData->{'structure'}, $bseData->{'bse'},
                $newBSEdata->{'bse'}, ['metal'] );
my $fake->{'complete'} = JSON::PP::false;
$newBSEdata->{'computer'} = {};
copyAndCompare( $newBSEdata->{'computer'}, $commonOceanData->{'computer'}, $bseData->{'computer'},
                $fake, ['para_prefix'] );

CalcParams( $newBSEdata, $commonOceanData, $bseData );

BSEParams( $newBSEdata, $commonOceanData, $bseData );

if( $haveScreen ) { 
  screenParams( $newBSEdata, $screenData, $bseData );
}


#ZNL, cks.normal, mode, nelectron, epsilon, screen.mode, conveps, haydockconv.ipt, xyz.wyck, lflag, bflag

# Figure out rpot names and copy them in

# write nval.h based on avecs and nelectron
# TODO fix nval.h or remove entirely / make Shirley self-energy dampening separate

open OUT, ">", "bse.json" or die;
print OUT $json->encode($newBSEdata);
close OUT;

#TODO remove/fix
# These likely only matter for valence calcs, and should be set to those defaults
open OUT, ">tmel_selector" or die;
print OUT "1\n";
close OUT;
open OUT, ">enk_selector" or die;
print OUT "0\n";
close OUT;
open OUT, ">bloch_selector" or die;
print OUT "3\n";
close OUT;


#TODO possibly delete?
foreach (@ExtraFiles) {
  copy( "$ENV{'OCEAN_BIN'}/$_", $_ ) or die "Failed to get ../$_\n$!";
}

# Grab files
grabWaveFunctions();
my @photon_files;

unless( $newBSEdata->{'calc'}->{'mode'} eq 'val' ) {
  unless( defined $newBSEdata->{'calc'}->{'photon_in'} ) {
    print "Undef\n\n";
    $newBSEdata->{'calc'}->{'photon_in'} = [];
  }
  @photon_files = grabPhotonFiles( $newBSEdata->{'calc'}->{'photon_in'} );

  makeSiteList( $newBSEdata );
  grabOPF( $newBSEdata );
  grabCKS( $newBSEdata );

  runMels( $newBSEdata, \@photon_files );

  writeRunlistCore( $newBSEdata, \@photon_files );
  writeBSEinCore( $newBSEdata );

  writeAuxFiles( $newBSEdata );
  writeScFac( $newBSEdata );
  writeCoreAuxFiles( $newBSEdata );
} else {
  writeRunlistVAL( );
  writeBSEinVal( $newBSEdata );
  writeAuxFiles( $newBSEdata );
  writeValAuxFiles( $newBSEdata );

  open OUT, ">", "ZNL" or die "$!";
  print OUT "0 0 0\n";
  close OUT;

  prepDen( );

  grabTmels();
}


grabScreenFiles( $newBSEdata );
####

#TODO:
# If statements for core/valence
# expand writeAux for valence

open OUT, ">", "bse.json" or die;
print OUT $json->encode($newBSEdata);
close OUT;

my $t0 = [gettimeofday];

#print "Running nothing!\n";
runOCEANx( $newBSEdata, "ocean.log" );


if( $newBSEdata->{'calc'}->{'mode'} eq 'rxs' ) {
  unless( defined $newBSEdata->{'calc'}->{'photon_out'} ) {
    print "Undef\n\n";
    $newBSEdata->{'calc'}->{'photon_out'} = [];
  }
  my @photon_out_files = grabPhotonFiles( $newBSEdata->{'calc'}->{'photon_out'} );

  # Need to runMels for any photons we haven't yet
  my @uniquePhoton;
  foreach my $a (@photon_out_files) {
    my $b = 1;
    foreach (@photon_files) {
      if ($a eq $_ ) {
        $b = 0;
        last;
      }
    }
    push @uniquePhoton, $a if( $b );
  }
  runMels( $newBSEdata, \@uniquePhoton );

  move ("runlist", "runlist.xas");
  move ("bse.in", "bse.in.xas");
  writeRunlistRIXS( $newBSEdata, \@photon_files, \@photon_out_files );
  writeBSEinVal( $newBSEdata );
  writeValAuxFiles( $newBSEdata );

  prepDen( );
  
# new aux files
  
  runOCEANx( $newBSEdata, "rixs.log" );
}

$newBSEdata->{'time'}->{'bse'} = tv_interval( $t0 );
$newBSEdata->{'time'}->{'total'} = 0;
foreach my $sec ('bse', 'prep', 'screen', 'dft' ) {
  $newBSEdata->{'time'}->{'total'} += $newBSEdata->{'time'}->{$sec};
}

$newBSEdata->{'bse'}->{'complete'} = JSON::PP::true;

open OUT, ">", "bse.json" or die;
print OUT $json->encode($newBSEdata);
close OUT;

exit 0;

sub runOCEANx {

  my ($hashRef, $logFileName) = @_;

  print "$hashRef->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/ocean.x > $logFileName 2>&1\n";
  system("$hashRef->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/ocean.x > $logFileName 2>&1" );
  if ($? == -1) {
      print "failed to execute: $!\n";
      die;
  }
  elsif ($? & 127) {
      printf "ocean died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
      die;
  }
  else {
    my $errorCode = $? >> 8;
    if( $errorCode != 0 ) {
      die "CALCULATION FAILED\n  ocean exited with value $errorCode\n";
    }
    else {
      printf "ocean exited successfully with value %d\n", $errorCode;
    }
  }
}

sub grabScreenFiles {
  my ($hashRef) = @_;

  # Change this to HLL instead of core+val
  if( $hashRef->{'calc'}->{'mode'} eq 'val' ) {
    return unless( $hashRef->{'bse'}->{'screen'}->{'mode'} eq 'grid' );
    die "Grid screening not written yet\n";
  }

  if( $hashRef->{'calc'}->{'mode'} eq 'rxs' ) {
    die "Grid screening not written yet\n" if( $hashRef->{'bse'}->{'screen'}->{'mode'} eq 'grid' );
  }

  if( $hashRef->{'calc'}->{'mode'} eq 'rxs' || $hashRef->{'calc'}->{'mode'} eq 'xas' || 
      $hashRef->{'calc'}->{'mode'} eq 'xes' ) {
    grabCoreScreenFiles( $hashRef );
  }
}

sub grabCoreScreenFiles {
  my ($hashRef) = @_;
 
  my $pawrad = sprintf "%.2f", $hashRef->{'bse'}->{'core'}->{'screen_radius'};
  for( my $i = 0; $i< scalar @{$hashRef->{'calc'}->{'edges'}}; $i++ ) {
    my @edge = split ' ', $hashRef->{'calc'}->{'edges'}[$i];
    my $elname = $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][0];
    my $elnum = $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][2];
    my $nnum = $edge[1];
    my $lnum = $edge[2];


    my $zstring = sprintf("z%2s%04i/n%02il%02i", $elname, $elnum, $nnum, $lnum);
    my $compactZstring = sprintf("z%2s%04i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
    print "$zstring   $compactZstring\n";
    print "../SCREEN/${zstring}/zR${pawrad}/rpot  rpot.${compactZstring}\n";
    unless( copy( "../SCREEN/${zstring}/zR${pawrad}/rpot", "rpot.${compactZstring}" ) == 1 )
    {
      $zstring = sprintf("z%2s%04i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
      copy( "../SCREEN/${zstring}/zR${pawrad}/rpot", "rpot.${zstring}" )
        or die "Failed to grab rpot\n../SCREEN/${zstring}/zR${pawrad}/rpot ./rpot.${zstring}\n";
    }

    if( $hashRef->{'bse'}->{'screen'}->{'core_offset'}->{'enable'} ) {
      copy( "../SCREEN/${zstring}/zR${pawrad}/cls", "cls.${compactZstring}" )
      or warn "WARNING!\nCore-level shift support requested, "
            . "but could not find ../SCREEN/${zstring}/zR${pawrad}/cls\n\$!"
            . "No CLS will be done for this site!\n";

    } else {  # If we don't want CLS then make sure the file is not here
      if( -e "cls.${compactZstring}" )
      {
        unlink "cls.${compactZstring}" or die "Failed to remove cls.${compactZstring}\n$!";
      }
    }
  }
}

#TODO: option to limit these to some number?
sub grabPhotonFiles {
  my $arrayRef = $_[0];
  my @photon_files;
  my $nphoton;

  my $haveList = 0;
  if( scalar @{$arrayRef} > 0 ) {
    $haveList = 1;
    print "Photon files specified\n";
  } else {
    print "Looking for available photon files:\n";
  }
  opendir DIR, "../" or die $!;
  while( my $file = readdir( DIR ) )
  {
    if( $file =~ m/^photon(\d+)$/ )
    {
      my $p = $1;
      if( $haveList ) {
        foreach (@{$arrayRef}) {
          if( $_ == $p ) {
            push @photon_files, $file;
            last;
          }
        }
      } else {
        push @photon_files, $file;
      }
    }
  }
  closedir DIR;

  if( $#photon_files == -1 )  # no photon files, fall back to jtv for now
  {
    opendir DIR, "../" or die $!;
    while( my $file = readdir( DIR ) )
    {
      if( $file =~ m/^jtv\d+$/ )
      {
        my $p = $1;
        if( $haveList ) {
          foreach (@{$arrayRef}) {
            if( $_ == $p ) {
              push @photon_files, $file;
              last;
            }
          }
        } else {
          push @photon_files, $file;
        }
      }
    }
    closedir DIR;
  }

  if( $#photon_files == -1 )
  {
    opendir DIR, "../Common" or die $!;
    my @tmp;
    while( my $file = readdir( DIR ) )
    {
      if( $file =~ m/^default_photon\d+$/ )
      {
        $file =~ s/default_//;
        my $p = $1;
        if( $haveList ) {
          foreach (@{$arrayRef}) {
            if( $_ == $p ) {
              push @photon_files, $file;
              last;
            }
          }
        } else {
          push @tmp, $file;
        }
      }
    }
    closedir DIR;
    if( $#tmp == -1 ) { 
      print "!!!  Could not find any photon files  !!!\n  Will have to quit :(\n";
      exit 1;
    }
    print "  Using default photon files!\n";
    my @sorted_photon_files = sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } @tmp;
    foreach( @sorted_photon_files )
    {
#      print "        $_\n";
      copy( catfile( updir(), "Common", 'default_' . $_), $_ ) or die "$!";
    }
    return @sorted_photon_files;
  }
  else
  {
    $nphoton = $#photon_files+1;
    if( $nphoton > 1 )
    {
      print "    Running with $nphoton photon files\n";
    }
    else
    {
      print "    Running with $nphoton photon file\n";
    }
    # Want to sort these back to numerical order (to avoid confusing/annoying users)
    my @sorted_photon_files = sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } @photon_files;
    foreach( @sorted_photon_files )
    {
      my $targ = basename($_);
      $targ =~ s/default_//;
#      print "        ../$_ $targ\n";
      copy( "../$_", $targ ) or die "$!";
    }
    return @sorted_photon_files;
  }
}

sub grabWaveFunctions {
  my $symlink_exists = eval { symlink("",""); 1 };

  unlink( "con.u2.dat" ) if( -e "con.u2.dat" );
  unlink( "val.u2.dat" ) if( -e "val.u2.dat" );
  unlink( "u2.dat" ) if( -e "u2.dat" );

  if( -e "../PREP/BSE/con.u2.dat" )
  {
    open OUT, ">bloch_selector" or die;
    print OUT "3\n";
    close OUT;
    if( $symlink_exists == 1 )
    {
      symlink( "../PREP/BSE/con.u2.dat", "con.u2.dat" ) or die "Failed to link ../PREP/BSE/con.u2.dat\n$!";
      symlink( "../PREP/BSE/val.u2.dat", "val.u2.dat" ) or die "Failed to link ../PREP/BSE/val.u2.dat\n$!";
    }
    else
    {
      copy( "../PREP/BSE/con.u2.dat", "con.u2.dat" ) or die "Failed to copy ../PREP/BSE/con.u2.dat\n$!";
      copy( "../PREP/BSE/val.u2.dat", "val.u2.dat" ) or die "Failed to copy ../PREP/BSE/val.u2.dat\n$!";
    }
  }
  elsif (-e "../PREP/BSE/u2.dat")
  {
    if( $symlink_exists == 1 )
    {
      symlink( "../PREP/BSE/u2.dat", "u2.dat" ) or die "Failed to link ../PREP/BSE/u2.dat\n$!";
    }
    else
    {
      copy( "../PREP/BSE/u2.dat", "u2.dat" ) or die "Failed to copy ../PREP/BSE/u2.dat\n$!";
    }
  }
  else
  {
    die "Failed to get electron wave functions from PREP/BSE\n";
  }

  copy( catfile( updir(), "PREP", "BSE", "wvfcninfo" ), "wvfcninfo" ) == 1 or die "Failed to copy wvfcninfo\n$!";
  copy( catfile( updir(), "PREP", "BSE", "wvfvainfo" ), "wvfvainfo" ) == 1 or die "Failed to copy wvfvainfo\n$!";
  copy( catfile( updir(), "PREP", "BSE", "enkfile" ), "enkfile" ) == 1 or die "Failed to copy enkfile\n$!";


#  copy(  catfile( updir(), "PREP", "BSE", "eshift.ipt" ), "eshift.ipt" ) == 1 or die "Failed to copy eshift.ipt\n$!";
#  print( "FIX ESHIFT\n");
  
}

sub CalcParams {
  my ($newRef, $commonRef, $oldRef ) = @_;

  $newBSEdata->{'calc'} = {};
  copyAndCompare( $newBSEdata->{'calc'}, $commonOceanData->{'calc'}, $bseData->{'calc'},
                  $newBSEdata->{'bse'}, [ 'mode' ] );

  if( $newBSEdata->{'calc'}->{'mode'} eq 'xas' || $newBSEdata->{'calc'}->{'mode'} eq 'xes' 
   || $newBSEdata->{'calc'}->{'mode'} eq 'rxs' ) {
    copyAndCompare( $newBSEdata->{'calc'}, $commonOceanData->{'calc'}, $bseData->{'calc'},
                  $newBSEdata->{'bse'}, [ 'edges', 'photon_in', 'photon_out' ] );
  }

  if( $newBSEdata->{'calc'}->{'mode'} eq 'val' || $newBSEdata->{'calc'}->{'mode'} eq 'rxs' ) {
    copyAndCompare( $newBSEdata->{'calc'}, $commonOceanData->{'calc'}, $bseData->{'calc'},
                  $newBSEdata->{'bse'}, [ 'photon_q' ] );
  }
}

sub BSEParams {
  my ($newRef, $commonRef, $oldRef ) = @_;

  my @type;
  if( $newRef->{'calc'}->{'mode'} eq 'xas' || $newRef->{'calc'}->{'mode'} eq 'xes' ) {
    $type[0] = 'core';
    push @type, "occupation";
  } elsif ( $newRef->{'calc'}->{'mode'} eq 'val' ) {
    $type[0] = 'val';
  } elsif ( $newRef->{'calc'}->{'mode'} eq 'rxs' ) {
    $type[0] = 'core'; $type[1] = 'val';
  } else {
    die "Unrecongnized calculations type: " . $newRef->{'calc'}->{'mode'} . "\n";
  }

  push @type, "nbands";

  copyAndCompare( $newRef->{'bse'}, $commonRef->{'bse'}, $oldRef->{'bse'},
                  $newRef->{'bse'}, \@type );
}

#TODO: This assumes core-level -- it'll flag things that don't matter for valence
sub screenParams {
  my ($newRef, $commonRef, $oldRef ) = @_;


  $newRef->{'bse'}->{'screen'} = {};
  copyAndCompare( $newRef->{'bse'}->{'screen'}, $commonRef->{'screen'}, $oldRef->{'bse'}->{'screen'},
                  $newRef->{'bse'}, [ 'mode', 'core_offset' ]  );
  copyAndCompare( $newRef->{'bse'}->{'screen'}, $commonRef->{'screen'}->{'grid2'}, 
                  $oldRef->{'bse'}->{'screen'},
                  $newRef->{'bse'}, [ 'lmax' ]  );

}

sub grabOPF {
  my ($hashRef) = @_;

  return unless( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'xes' || 
                 $hashRef->{'calc'}->{'mode'} eq 'rxs' ) ;

  my %uniqueZNL;
  my %uniqueZ;
  foreach my $edge (@{$hashRef->{'calc'}->{'edges'}}) {
    my @edge = split ' ', $edge;
    my $z = sprintf "%i", $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][1];
    my $tag = sprintf "%i %i %i", $z, $edge[1], $edge[2];
    $uniqueZNL{ $tag } = [ $z, $edge[1], $edge[2] ] ;
    $uniqueZ{ $z } = 1;
  }

  foreach my $key (keys %uniqueZNL) {
    my $zstring = sprintf("z%03in%02il%02i", $uniqueZNL{$key}[0], $uniqueZNL{$key}[1], $uniqueZNL{$key}[2]);
    my @fileList = glob( catdir( updir(), "OPF", "zpawinfo", "?k???$zstring" ) );
#    push @fileList, catfile( updir(), "OPF", "zpawinfo", "melfile$zstring" );
    push @fileList, catfile( updir(), "OPF", "zpawinfo", "coreorb$zstring" );

    foreach (@fileList ) { copy( $_, basename($_ ) ); }
  }

  foreach my $zee (keys %uniqueZ ) {
    my $z = sprintf "z%03i", $zee;
    my @fileList = glob( catdir( updir(), "OPF", "zpawinfo", "phrc?$z" ) );
    push @fileList, ( catfile( updir(), "OPF", "zpawinfo", "prjfile$z" ) );
    push @fileList, glob( catdir( updir(), "OPF", "zpawinfo", "ps?$z" ) );
    push @fileList, glob( catdir( updir(), "OPF", "zpawinfo", "ae?$z" ) );
    push @fileList, ( catfile( updir(), "OPF", "zpawinfo", "corezeta$z" ) );
    push @fileList, ( catfile( updir(), "OPF", "zpawinfo", "radfile$z" ) );
#    push @fileList, glob( catdir( updir(), "OPF", "zpawinfo", "ft?$z" ) );

    push @fileList, glob( catdir( updir(), "OPF", "zdiag$z", "c2c${z}n??l??" ) );
    push @fileList, glob( catdir( updir(), "OPF", "zdiag$z", "s2c${z}n??l??" ) );
    push @fileList, glob( catdir( updir(), "OPF", "zdiag$z", "melsemi${z}n??l??" ) );
    push @fileList, glob( catdir( updir(), "OPF", "zpawinfo", "melfile${z}n??l??" ) );
    foreach (@fileList ) { copy( $_, basename($_ ) ) == 1 or die "Failed to copy $_\n$!"; }
  }

}

sub grabCKS {
  my ($hashRef) = @_;

  return unless( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'xes' ||
                 $hashRef->{'calc'}->{'mode'} eq 'rxs' ) ;

  my $symlink_exists = eval { symlink("",""); 1 };

  foreach my $edge (@{$hashRef->{'calc'}->{'edges'}}) {
    my @edge = split ' ', $edge;
    my $cksTag = sprintf "%s%04i", $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][0], 
                                   $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][2];
    foreach my $t ( "cksc.", "cksv.", "parcksc.", "parcksv." ) {
      unlink( "${t}${cksTag}" ) if( -e "${t}${cksTag}" );
    }
  
    if( -e catfile( updir(), "PREP", "BSE", "parcksc.$cksTag" ) ) {
      (copy( catfile( updir(), "PREP", "BSE", "parcksc.$cksTag" ), "parcksc.$cksTag" ) ) == 1 
        or die "Failed to grab " . catfile( updir(), "PREP", "BSE", "parcksc.$cksTag" ) ."\n$!";
    } else {
      (copy( catfile( updir(), "PREP", "BSE", "cksc.$cksTag" ), "cksc.$cksTag" ) ) == 1
        or die "Failed to grab " . catfile( updir(), "PREP", "BSE", "cksc.$cksTag" ) ."\n$!";
    }

    if( -e catfile( updir(), "PREP", "BSE", "parcksv.$cksTag" ) ) {
      (copy( catfile( updir(), "PREP", "BSE", "parcksv.$cksTag" ), "parcksv.$cksTag" ) ) == 1
        or die "Failed to grab " . catfile( updir(), "PREP", "BSE", "parcksv.$cksTag" ) ."\n$!";
    } else {
      (copy( catfile( updir(), "PREP", "BSE", "cksv.$cksTag" ), "cksv.$cksTag" ) ) == 1
        or die "Failed to grab " . catfile( updir(), "PREP", "BSE", "cksv.$cksTag" ) ."\n$!";
    }

  }
}

sub runMels {
  my ($hashRef, $photon_files ) = @_;


  return unless( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'xes' ||
                 $hashRef->{'calc'}->{'mode'} eq 'rxs' ) ;


  my %uniqueZNL;
  foreach my $edge (@{$hashRef->{'calc'}->{'edges'}}) {
    my @edge = split ' ', $edge;
    my $z = sprintf "%i", $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][1];
    my $tag = sprintf "%i %i %i", $z, $edge[1], $edge[2];
    $uniqueZNL{ $tag } = [ $z, $edge[1], $edge[2] ] ;
  }

  foreach my $key (keys %uniqueZNL) {
    open ZNL, ">ZNL" or die;
    printf ZNL "%i %i %i\n", $uniqueZNL{$key}[0], $uniqueZNL{$key}[1], $uniqueZNL{$key}[2];
    close ZNL;


    foreach my $way (@{$photon_files}) {
      copy( $way, "spectfile" ) or die "Failed to copy $way\n$!";
      system("$ENV{'OCEAN_BIN'}/meljtv.x");
      $way =~ m/(\d+)$/ or die "Misformed photon file name\n";
      my $i = $1;
      my $mel_targ = sprintf("mels.z%03un%02ul%02up%02u", $uniqueZNL{$key}[0], 
                                $uniqueZNL{$key}[1], $uniqueZNL{$key}[2], $i);
      move( "mels", $mel_targ ) or die "Failed to move mels to $mel_targ\n$!";
    }
  }

}


sub writeRunlistCore {
  my ($hashRef, $photon_files ) = @_;

  my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

  my $runLength = 1;
  if( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'xes' 
    || $hashRef->{'calc'}->{'mode'} eq 'rxs' ) {
    $runLength = scalar @{$photon_files} * scalar @{$hashRef->{'calc'}->{'edges'}};
  } else {
    die "Runlist only written for XAS or XES at the moment\n";
  }

  open RUNLIST, ">runlist";
  print RUNLIST $runLength . "\n";

  foreach my $edge ( @{$hashRef->{'calc'}->{'edges'}}) {
    my @edge = split ' ', $edge;
    my $znum = sprintf "%i", $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][1];
    my $nnum = $edge[1];
    my $lnum = $edge[2];
    my $elname = $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][0];
    my $elnum = $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][2];
    my $run_text = uc( $hashRef->{'calc'}->{'mode'} );
    $run_text = 'XAS' if( $run_text eq 'RXS' );
    foreach my $way (@{$photon_files} ) {
      $way =~ m/(\d+)$/ or die "Malformed photon file name:\t$way\n";
      my $i = $1;
      print RUNLIST "$znum  $nnum  $lnum  $elname  ${nnum}$alphal{$lnum}  $elnum  $i  $run_text\n";
    }
  }

  close RUNLIST;
}

sub writeRunlistRIXS {
  my ($hashRef, $photon_files, $photon_out_files ) = @_;

  my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

  my @edge = split ' ', $hashRef->{'calc'}->{'edges'}[0];
  my $znum = sprintf "%i", $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][1];
  my $nnum = $edge[1];
  my $lnum = $edge[2];
  my $elname = $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][0];
  my $corelevel = "${nnum}$alphal{$lnum}";
  my $calc = 'RXS';

  my $nphoton = 0;
  my @photon_combo;
  foreach my $p_in ( @{$photon_files} ) {
    $p_in =~ m/(\d+)\s*$/ or die "Malformed photon file name:\t$p_in\n"; 
    my $i = $1;
    foreach my $p_out ( @{$photon_out_files}) {
      $p_out =~ m/(\d+)\s*$/ or die "Malformed photon file name:\t$p_in\n"; 
      my $j = $1;
      next if ($p_in eq $p_out );
      $nphoton ++;
      push @photon_combo, [$i, $j ]; 
    }
  }

  open RUNLIST, ">", "runlist" or die "Failed to open file runlist\n$!";
  printf RUNLIST "%i\n", $nphoton * $hashRef->{'bse'}->{'core'}->{'gmres'}->{'count'};
  for( my $e = 1; 
       $e <= $hashRef->{'bse'}->{'core'}->{'gmres'}->{'count'};
       $e++ ) {
    for( my $i = 0; $i < scalar @photon_combo; $i++ ) { 
      print RUNLIST "$znum $nnum  $lnum  $elname  $corelevel  0  $photon_combo[$i][0]  $calc  $e  $photon_combo[$i][1]\n";
    }
  }
        
  
  close RUNLIST;
}

sub writeRunlistVAL {
  open RUNLIST, ">", "runlist" or die "Failed to open file runlist\n$!";
  print RUNLIST "1\n0  0  0  __  __  0  0  VAL\n";
  close RUNLIST;
}


sub writeBSEinVal {
  my ($hashRef) = @_;

  open BSE, ">", "bse.in" or die "$!\n";
  print BSE "0 0 0\n0 0 0\n";
  
  if( $hashRef->{'bse'}->{'val'}->{'solver'} eq 'haydock' ) {
    printf BSE "hay\n%i %i %g %g %g %g\n", 
      $hashRef->{'bse'}->{'val'}->{'haydock'}->{'niter'},
      $hashRef->{'bse'}->{'val'}->{'plot'}->{'points'},
      $hashRef->{'bse'}->{'val'}->{'plot'}->{'range'}[0],
      $hashRef->{'bse'}->{'val'}->{'plot'}->{'range'}[1],
      $hashRef->{'bse'}->{'val'}->{'broaden'}, 
      0.0;
  } elsif( $hashRef->{'bse'}->{'val'}->{'solver'} eq 'gmres' ) {
    printf BSE "inv\n%i %g %g %g %g\n", $hashRef->{'bse'}->{'val'}->{'gmres'}->{'nloop'},
        $hashRef->{'bse'}->{'val'}->{'broaden'}, $hashRef->{'bse'}->{'val'}->{'gmres'}->{'gprc'},
        $hashRef->{'bse'}->{'val'}->{'gmres'}->{'ffff'}, 0.0;
    if( $hashRef->{'bse'}->{'val'}->{'gmres'}->{'estyle'} eq 'list' ) {
      printf BSE "list\n%i\n", scalar @{$hashRef->{'bse'}->{'val'}->{'gmres'}->{'elist'}};
      $hashRef->{'bse'}->{'val'}->{'gmres'}->{'count'} =
            scalar @{$hashRef->{'bse'}->{'val'}->{'gmres'}->{'elist'}};
      foreach (@{$hashRef->{'bse'}->{'val'}->{'gmres'}->{'elist'}}) {
        printf BSE "%g\n", $_;
      }
    } elsif( $hashRef->{'bse'}->{'val'}->{'gmres'}->{'estyle'} eq 'range' ) {
      printf BSE "loop\n%g %g %g\n",
              $hashRef->{'bse'}->{'val'}->{'gmres'}->{'erange'}[0],
              $hashRef->{'bse'}->{'val'}->{'gmres'}->{'erange'}[1],
              $hashRef->{'bse'}->{'val'}->{'gmres'}->{'erange'}[2];
      $hashRef->{'bse'}->{'val'}->{'gmres'}->{'count'} =
          floor( ($hashRef->{'bse'}->{'val'}->{'gmres'}->{'erange'}[1] -
                  $hashRef->{'bse'}->{'val'}->{'gmres'}->{'erange'}[0] +
             0.9* $hashRef->{'bse'}->{'val'}->{'gmres'}->{'erange'}[2] ) /
                 $hashRef->{'bse'}->{'val'}->{'gmres'}->{'erange'}[2] );
    } else { die "Unsupported estyle in val->gmres\n"; }
      printf "%i energy steps with GMRES\n", $hashRef->{'bse'}->{'val'}->{'gmres'}->{'count'};
  } else {
    die "Only Haydock and GMRES solvers implemented for valence\n";
  }
  close BSE;

  open SPECT, ">", "spect.in" or die "Failed to open spect.in for writing\n$!";
  printf SPECT "%i %g %g %g %g %g\n", 
      $hashRef->{'bse'}->{'val'}->{'plot'}->{'points'},
      $hashRef->{'bse'}->{'val'}->{'plot'}->{'range'}[0],
      $hashRef->{'bse'}->{'val'}->{'plot'}->{'range'}[1],
      $hashRef->{'bse'}->{'val'}->{'broaden'},
      0.0;
  close SPECT;

}

sub writeBSEinCore {
  my ($hashRef) = @_;

  my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

  open BSE, ">", "bse.in" or die "$!\n";
  if ( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'xes' ||
                 $hashRef->{'calc'}->{'mode'} eq 'rxs' ) {
    #TODO only grabs the first listed, can't mix multiple edges in a single ocean.x run yet
    my @edge = split ' ', $hashRef->{'calc'}->{'edges'}[0];
    my $znum = $hashRef->{'calc'}->{'sitelist'}[$edge[0]-1][1];
    my $nnum = $edge[1];
    my $lnum = $edge[2];
    my $filename = sprintf("prjfilez%03u", $znum );
    my $filename = sprintf("prjfilez%03u", $znum );

    open TMPFILE, $filename or die "Failed to open $filename\n";
    my $line = <TMPFILE>;
    close TMPFILE;
    $line =~ m/(\d+)\s+(\d+)\s+\d+\s+\S+/ or die "Failed to match first line of $filename\n";

    print BSE "$lnum\t$1\t$2\n";
    

    if( $hashRef->{'bse'}->{'core'}->{'spin_orbit'} < 0 ) {
      my $lookup = sprintf("%1u%1s", $nnum, $alphal{$lnum}) or die;
      my $filename = sprintf("corezetaz%03u", $znum);
      print "$lookup\t$filename\n";
      $line = `grep $lookup $filename`;
    }
    else
    {
      print "Overiding spin-orbit splitting! Using: " . 
            $hashRef->{'bse'}->{'core'}->{'spin_orbit'} . " (eV)\n";
      $line = $hashRef->{'bse'}->{'core'}->{'spin_orbit'} . "\n";
    }
    print BSE $line;

    if( $hashRef->{'bse'}->{'core'}->{'solver'} eq 'haydock' ) {
      printf BSE "hay\n%i  %i  %f  %f  %f  0.000\n", $hashRef->{'bse'}->{'core'}->{'haydock'}->{'niter'}, 
          $hashRef->{'bse'}->{'core'}->{'plot'}->{'points'}, 
          $hashRef->{'bse'}->{'core'}->{'plot'}->{'range'}[0],
          $hashRef->{'bse'}->{'core'}->{'plot'}->{'range'}[1],
          $hashRef->{'bse'}->{'core'}->{'broaden'};
      open SP, ">", "spect.in" or die "$!\n";
      printf SP "%i  %f  %f  %f  0.000\n", $hashRef->{'bse'}->{'core'}->{'plot'}->{'points'},
          $hashRef->{'bse'}->{'core'}->{'plot'}->{'range'}[0],
          $hashRef->{'bse'}->{'core'}->{'plot'}->{'range'}[1],
          $hashRef->{'bse'}->{'core'}->{'broaden'};
      close SP;
    } elsif( $hashRef->{'bse'}->{'core'}->{'solver'} eq 'gmres' ) {
      printf BSE "inv\n%i %g %g %g %g\n", $hashRef->{'bse'}->{'core'}->{'gmres'}->{'nloop'}, 
          $hashRef->{'bse'}->{'core'}->{'broaden'}, $hashRef->{'bse'}->{'core'}->{'gmres'}->{'gprc'}, 
          $hashRef->{'bse'}->{'core'}->{'gmres'}->{'ffff'}, 0.0;
      if( $hashRef->{'bse'}->{'core'}->{'gmres'}->{'estyle'} eq 'list' ) {
        printf BSE "list\n%i\n", scalar @{$hashRef->{'bse'}->{'core'}->{'gmres'}->{'elist'}};
        $hashRef->{'bse'}->{'core'}->{'gmres'}->{'count'} = 
              scalar @{$hashRef->{'bse'}->{'core'}->{'gmres'}->{'elist'}};
        foreach (@{$hashRef->{'bse'}->{'core'}->{'gmres'}->{'elist'}}) {
          printf BSE "%g\n", $_;
        }
      } elsif( $hashRef->{'bse'}->{'core'}->{'gmres'}->{'estyle'} eq 'range' ) {
        printf BSE "loop\n%g %g %g\n", 
                $hashRef->{'bse'}->{'core'}->{'gmres'}->{'erange'}[0],
                $hashRef->{'bse'}->{'core'}->{'gmres'}->{'erange'}[1],
                $hashRef->{'bse'}->{'core'}->{'gmres'}->{'erange'}[2];
        $hashRef->{'bse'}->{'core'}->{'gmres'}->{'count'} =
            floor( ($hashRef->{'bse'}->{'core'}->{'gmres'}->{'erange'}[1] -
                    $hashRef->{'bse'}->{'core'}->{'gmres'}->{'erange'}[0] +
               0.9* $hashRef->{'bse'}->{'core'}->{'gmres'}->{'erange'}[2] ) /
                   $hashRef->{'bse'}->{'core'}->{'gmres'}->{'erange'}[2] );
      } else { die "Unsupported estyle in core->gmres\n"; }
      printf "%i energy steps with GMRES\n", $hashRef->{'bse'}->{'core'}->{'gmres'}->{'count'};
    } else { die "Failed to write gmres in writeBSEinCore\n"; }


    open OUT, ">", "gaussBroaden.ipt" or die "Failed to open gaussBroaden.ipt\n$!";
    printf OUT "%g\n", $hashRef->{'bse'}->{'core'}->{'gauss_broaden'};
    close OUT;

  } else { die "Unsupported calc mode in writeBSEinCore\n"; }


  close BSE;
}

sub writeAuxFiles {
  my ($hashRef) = @_;

  open OUT, ">", "xmesh.ipt" or die;
  printf OUT "%i %i %i\n", $hashRef->{'bse'}->{'xmesh'}[0], $hashRef->{'bse'}->{'xmesh'}[1], $hashRef->{'bse'}->{'xmesh'}[2];
  close OUT;

  open OUT, ">", "kmesh.ipt" or die;
  printf OUT "%i %i %i\n", $hashRef->{'bse'}->{'kmesh'}[0], $hashRef->{'bse'}->{'kmesh'}[1], $hashRef->{'bse'}->{'kmesh'}[2];
  close OUT;

  open OUT, ">", "nspin" or die;
  print OUT $hashRef->{'bse'}->{'nspin'} . "\n";
  close OUT;

  my $nb;
  if ( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'rxs' ) {
    if( $hashRef->{'bse'}->{'nbands'} < $hashRef->{'bse'}->{'brange'}[3] ) {
      $nb = $hashRef->{'bse'}->{'nbands'} - $hashRef->{'bse'}->{'brange'}[2] + 1;
    } else {
      $nb = $hashRef->{'bse'}->{'brange'}[3] - $hashRef->{'bse'}->{'brange'}[2] + 1;
    }
  } else {
    $nb = $hashRef->{'bse'}->{'brange'}[1] - $hashRef->{'bse'}->{'brange'}[0] + 1;
  }
  open OUT, ">", "nbuse.ipt" or die;
  print OUT "$nb\n";
  close OUT;

  open OUT, ">", "qinunitsofbvectors.ipt" or die;
  printf OUT "%f %f %f\n", $hashRef->{'bse'}->{'photon_q'}[0], $hashRef->{'bse'}->{'photon_q'}[1], 
                          $hashRef->{'bse'}->{'photon_q'}[2];
  close OUT;

  open OUT, ">", "brange.ipt" or die;
  printf OUT "%i %i\n%i %i\n", $hashRef->{'bse'}->{'brange'}[0], $hashRef->{'bse'}->{'brange'}[1],
                              $hashRef->{'bse'}->{'brange'}[2], $hashRef->{'bse'}->{'brange'}[3];
  close OUT;
  
  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "%.16g  %.16g  %.16g\n", $hashRef->{'bse'}->{'avecs'}[$i][0],
                                $hashRef->{'bse'}->{'avecs'}[$i][1],
                                $hashRef->{'bse'}->{'avecs'}[$i][2];

  }
  close OUT;

  open OUT, ">", "bvecs" or die "Failed to open bvecs\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "%.16g  %.16g  %.16g\n", $hashRef->{'bse'}->{'bvecs'}[$i][0],
                                $hashRef->{'bse'}->{'bvecs'}[$i][1],
                                $hashRef->{'bse'}->{'bvecs'}[$i][2];

  }
  close OUT;

  #TODO: cut out of ocean.x
  open OUT, ">", "cks.normal" or die "$!\n";
  if( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'rxs' 
                                            || $hashRef->{'calc'}->{'mode'} eq 'val' ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

  # For all options that could be either core (or valence!) 
  my $valORcore = 'core';
  $valORcore = 'val' if( $hashRef->{'calc'}->{'mode'} eq 'val' );
#  if( $hashRef->{'calc'}->{'mode'} eq 'val' ) { die "Valence not programed yet\n"; }
  #TODO: this is redundant
  open OUT, ">", "mode" or die "$!";
  printf OUT "%g %i\n", $hashRef->{'bse'}->{$valORcore}->{'strength'}, 
                       $hashRef->{'bse'}->{$valORcore}->{'haydock'}->{'niter'};
  close OUT;
    

  open OUT, ">", "nelectron" or die "$!";
  print OUT $hashRef->{'bse'}->{'nelec'} ."\n";
  close OUT;

  open OUT, ">", "nedges" or die "$!";
  if( exists $hashRef->{'calc'}->{'edges'} ) {
    printf OUT "%i\n", scalar @{$hashRef->{'calc'}->{'edges'}};
  } else {
    print OUT "0\n";
  }
  close OUT;

  open OUT, ">", "epsilon" or die "$!";
  printf OUT "%g\n", $hashRef->{'bse'}->{'epsilon'};
  close OUT;

#  if( $hashRef->{'calc'}->{'mode'} eq 'val' ) { die "Valence not programed yet\n"; }
  # TODO, enable for valence BSE, add to input options
  open OUT, ">", "cnbse.write_rhs" or die "$!";
  if( $hashRef->{'bse'}->{'core'}->{'write_rhs'} ) { 
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

#  if( $hashRef->{'calc'}->{'mode'} eq 'val' ) { die "Valence not programed yet\n"; }
  open OUT, ">", "haydockconv.ipt" or die;
  printf OUT "%g  %i\n", $hashRef->{'bse'}->{$valORcore}->{'haydock'}->{'converge'}->{'thresh'}, 
                         $hashRef->{'bse'}->{$valORcore}->{'haydock'}->{'converge'}->{'spacing'};
  close OUT;

  open OUT, ">", "xyz.wyck" or die "$!";
  printf OUT "%i\n", scalar @{$hashRef->{'bse'}->{'typat'}};
  for( my $i=0; $i< scalar @{$hashRef->{'bse'}->{'typat'}}; $i++ )
  {
    my $s = sprintf " %s %20.16f %20.16f %20.16f", $hashRef->{'bse'}->{'elname'}[$hashRef->{'bse'}->{'typat'}[$i]-1],
                  $hashRef->{'bse'}->{'xred'}[$i][0], $hashRef->{'bse'}->{'xred'}[$i][1], $hashRef->{'bse'}->{'xred'}[$i][2];
    print OUT "$s\n";
  }
  close OUT;

  my $volume = $hashRef->{'bse'}->{'avecs'}[0][0] * ($hashRef->{'bse'}->{'avecs'}[1][1] * $hashRef->{'bse'}->{'avecs'}[2][2] - $hashRef->{'bse'}->{'avecs'}[2][1] * $hashRef->{'bse'}->{'avecs'}[1][2] )
             - $hashRef->{'bse'}->{'avecs'}[1][0] * ($hashRef->{'bse'}->{'avecs'}[0][1] * $hashRef->{'bse'}->{'avecs'}[2][2] - $hashRef->{'bse'}->{'avecs'}[2][1] * $hashRef->{'bse'}->{'avecs'}[0][2] )
             + $hashRef->{'bse'}->{'avecs'}[2][0] * ($hashRef->{'bse'}->{'avecs'}[0][1] * $hashRef->{'bse'}->{'avecs'}[1][2] - $hashRef->{'bse'}->{'avecs'}[1][1] * $hashRef->{'bse'}->{'avecs'}[0][2] );
  open OUT, ">", "nval.h" or die;
  printf OUT "%.16e\n", $hashRef->{'bse'}->{'nelec'} / $volume;
  close OUT;

  open OUT, ">", "efermiinrydberg.ipt" or die "$!";
  printf OUT "%.16f\n", 2*$hashRef->{'bse'}->{'fermi'};
  close OUT;

  open OUT, ">", "metal" or die;
  if( $hashRef->{'bse'}->{'metal'} ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;
  
  open OUT, ">", "eshift.ipt" or die "Failed to open eshift.ipt\n$!";
  printf OUT "%.16f\n", $hashRef->{'bse'}->{'cbm'}*-1;
  close OUT;

  open OUT, ">", "occupation.ipt" or die "Failed to open occupation.ipt\n$!";
  if( exists $hashRef->{'bse'}->{'occupation'} ) {
    printf OUT "%s  %.10g\n", $hashRef->{'bse'}->{'occupation'}->{'type'}, $hashRef->{'bse'}->{'occupation'}->{'value'};
  } else {
    $hashRef->{'bse'}->{'occupation'}->{'type'} = 'none';
    $hashRef->{'bse'}->{'occupation'}->{'value'} = 0.0;
    print OUT "none 0.0\n";
  }
  close OUT;

  #TODO:
  # Add option for valence GMRES print out
  open OUT, ">", "echamp.inp" or die "Failed to open echamp.inp\n$!";
  if( $hashRef->{'bse'}->{'core'}->{'gmres'}->{'echamp'} || $hashRef->{'calc'}->{'mode'} eq 'rxs' ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
}

sub writeScFac {
  my ($hashRef) = @_;

  my $scfac = 1;
  if( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'rxs' ) {
    $scfac = $hashRef->{'bse'}->{'core'}->{'scfac'};
    if( $scfac < 0 ) {
      my @tmp = split( ' ', $hashRef->{'calc'}->{'edges'}[0] );
      my $z = @{$hashRef->{'bse'}->{'znucl'}}[@{$hashRef->{'bse'}->{'typat'}}[@tmp[0]-1]-1];
      print "SCFAC: $z\n";
      $scfac = 0.8;
      $scfac = 0.7 if( $z >= 55 );
    }
  }
  open OUT, ">", "scfac" or die "Failed to open scfac\n$!";
  printf OUT "%g\n", $scfac;
  close OUT;

}

#TODO fix this up
sub prepDen {
  copy ( "../DFT/rhoofr", "rhoofr" ) or die "Failed to get rhoofr\n";
  `tail -n 1 rhoofr > nfft`;
}

sub writeCoreAuxFiles {
  my ($hashRef) = @_;

  #TODO GW control
  if( exists $hashRef->{'bse'}->{'core'}->{'gw'}->{'control'} ) {
    open OUT, ">", "gw_control" or die "Failed to open gw_control\n$!";
    print OUT $hashRef->{'bse'}->{'core'}->{'gw'}->{'control'} . "\n";
    close OUT;

    if( $hashRef->{'bse'}->{'core'}->{'gw'}->{'control'} eq 'cstr' ) {
  
      open OUT, ">", "gw_core_cstr" or die "Failed to open gw_core_cstr\n$!";
      printf OUT "%g ", $hashRef->{'bse'}->{'core'}->{'gw'}->{'cstr'}->{'gap'};
      if( $hashRef->{'bse'}->{'core'}->{'gw'}->{'cstr'}->{'abs_gap'} == $JSON::PP::true ) {
        print OUT "true "; 
      } else {
        print OUT "false ";
      } 
      printf OUT "%g %g\n", $hashRef->{'bse'}->{'core'}->{'gw'}->{'cstr'}->{'vstr'},
                            $hashRef->{'bse'}->{'core'}->{'gw'}->{'cstr'}->{'cstr'};
      close OUT;
    }
  } 
  
}

sub writeValAuxFiles {
  my ($hashRef) = @_;

  my @files = ( 'aldaf', 'backf', 'bande', 'bflag', 'bwflg', 'lflag', 'qpflg', 'semitda' );
  foreach my $file (@files) {
    open OUT, ">", $file or die "Failed to open $file\n$!";
    if( $hashRef->{'bse'}->{'val'}->{ $file } == $JSON::PP::true ) {
      print OUT "1\n"; 
    } else {
      print OUT "0\n";
    }
  }

  if( $hashRef->{'bse'}->{'occupation'}->{'type'} eq 'none' ) {
    move "brange.ipt", "file_brange.ipt";
    open OUT, ">", "brange.ipt" or die;
    my $b2 = $hashRef->{'bse'}->{'brange'}[2];
    if ( exists $hashRef->{'bse'}->{'min_unocc_band'} ) {
      $b2 = $hashRef->{'bse'}->{'min_unocc_band'};
    }
    my $b1 = $hashRef->{'bse'}->{'brange'}[1];
    if ( exists $hashRef->{'bse'}->{'max_occ_band'} ) {
      $b1 = $hashRef->{'bse'}->{'max_occ_band'};
    }
    printf OUT "%i %i\n%i %i\n", $hashRef->{'bse'}->{'brange'}[0], $b1,
                                $b2, $hashRef->{'bse'}->{'brange'}[3];
    close OUT;
  }

  if( exists $hashRef->{'bse'}->{'val'}->{'epsilon_threshold'} ) {
    open OUT, ">", "conveps.ipt" or die "Failed to open conveps.ipt\n$!";
    printf OUT "%g\n", $hashRef->{'bse'}->{'val'}->{'epsilon_threshold'};
    close OUT;
  }

  #TODO GW control
  if( $hashRef->{'bse'}->{'val'}->{'gw'}->{'control'} eq 'cstr' ) {
    open OUT, ">", "gw_control" or die "Failed to open gw_control\n$!";
    print OUT "cstr\n";
    close OUT;

    open OUT, ">", "gw_val_cstr" or die "Failed to open gw_val_cstr\n$!";
    printf OUT "%g ", $hashRef->{'bse'}->{'val'}->{'gw'}->{'cstr'}->{'gap'};
    if( $hashRef->{'bse'}->{'val'}->{'gw'}->{'cstr'}->{'abs_gap'} == $JSON::PP::true ) {
      print OUT "true ";
    } else {
      print OUT "false ";
    }
    printf OUT "%g %g\n", $hashRef->{'bse'}->{'val'}->{'gw'}->{'cstr'}->{'vstr'}, 
                          $hashRef->{'bse'}->{'val'}->{'gw'}->{'cstr'}->{'cstr'};
    close OUT;
  } else {
    open OUT, ">", "gw_control" or die "Failed to open gw_control\n$!";
    print OUT "none\n";
    close OUT;
  }
    
  open OUT, ">", "decut" or die "Failed to open decut\n$!";
  print OUT $hashRef->{'bse'}->{'val'}->{'decut'} . "\n";
  close OUT;

  open OUT, ">", "haydockconv.ipt" or die;
  printf OUT "%g  %i\n", 
        $hashRef->{'bse'}->{'val'}->{'haydock'}->{'converge'}->{'thresh'},
        $hashRef->{'bse'}->{'val'}->{'haydock'}->{'converge'}->{'spacing'};
  close OUT;


}


sub makeSiteList {
  my ($hashRef) = @_;

  return unless( $hashRef->{'calc'}->{'mode'} eq 'rxs' || 
                 $hashRef->{'calc'}->{'mode'} eq 'xas' ||
                 $hashRef->{'calc'}->{'mode'} eq 'xes' );
  my %countByName;
  $hashRef->{'calc'}->{'sitelist'} = [];
        
  for( my $i=0; $i< scalar @{$hashRef->{'bse'}->{'typat'}}; $i++ ) {
    if( exists $countByName{$hashRef->{'bse'}->{'elname'}[$hashRef->{'bse'}->{'typat'}[$i]-1]} )
    {
      $countByName{$hashRef->{'bse'}->{'elname'}[$hashRef->{'bse'}->{'typat'}[$i]-1]}++;
    } else
    {       
      $countByName{$hashRef->{'bse'}->{'elname'}[$hashRef->{'bse'}->{'typat'}[$i]-1]}=1;
    }
    my @tmp;  
    $tmp[0] = $hashRef->{'bse'}->{'elname'}[$hashRef->{'bse'}->{'typat'}[$i]-1];
    $tmp[1] = $hashRef->{'bse'}->{'znucl'}[$hashRef->{'bse'}->{'typat'}[$i]-1];
    $tmp[2] = $countByName{$hashRef->{'bse'}->{'elname'}[$hashRef->{'bse'}->{'typat'}[$i]-1]};
    push @{$hashRef->{'calc'}->{'sitelist'}}, \@tmp;
    print "$tmp[0]  $tmp[1]  $tmp[2]  \n";
  }

}

sub grabTmels {
  copy( catfile( updir(), "PREP", "BSE", "tmels.info" ), "tmels.info" ) or die "Failed to copy tmels.info\n$!";
  copy( catfile( updir(), "PREP", "BSE", "ptmels.dat" ), "ptmels.dat" ) or die "Failed to copy ptmels.dat\n$!";
}
