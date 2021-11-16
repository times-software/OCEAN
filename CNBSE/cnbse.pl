#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

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

# Load run info for SCREEN
$dataFile = catfile( updir(), "SCREEN", "screen.json" );
die "Failed to find $dataFile\n" unless( -e $dataFile );

my $screenData;
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


$newBSEdata->{'bse'}->{'complete'} = JSON::PP::true;

# (was previous run done)
$newBSEdata->{'bse'}->{'complete'} = JSON::PP::false
    unless( exists $bseData->{'bse'} && exists $bseData->{'bse'}->{'complete'} 
            && $bseData->{'bse'}->{'complete'} );


# Has the structure been updated? But we want this from DFT?
#copyAndCompare( $newBSEdata, $commonOceanData, $bseData,
#                $newBSEdata->{'bse'}, ["structure"] );

my @list = ( "kmesh", "kshift", "xmesh", "nspin", "photon_q", "brange", "xred", "typat", 
             "znucl", "elname", "nelec", "avecs", "bvecs", "fermi", "epsilon" );
copyAndCompare( $newBSEdata->{'bse'}, $prepData->{'bse'}, $bseData->{'bse'},
                $newBSEdata->{'bse'}, \@list );


CalcParams( $newBSEdata, $commonOceanData, $bseData );

BSEParams( $newBSEdata, $commonOceanData, $bseData );

screenParams( $newBSEdata, $screenData, $bseData );


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
print OUT "0\n";
close OUT;
open OUT, ">enk_selector" or die;
print OUT "0\n";
close OUT;
open OUT, ">bloch_selector" or die;
print OUT "0\n";
close OUT;


#TODO possibly delete?
foreach (@ExtraFiles) {
  copy( "$ENV{'OCEAN_BIN'}/$_", $_ ) or die "Failed to get ../$_\n$!";
}

# Grab files
grabWaveFunctions();
my @photon_files = grabPhotonFiles();

grabScreenFiles( $newBSEdata );

grabOPF( $newBSEdata );
grabCKS( $newBSEdata );

runMels( $newBSEdata, \@photon_files );


####

writeRunlist( $newBSEdata, \@photon_files );
writeBSEin( $newBSEdata );

writeAuxFiles( $newBSEdata );

open OUT, ">", "bse.json" or die;
print OUT $json->encode($newBSEdata);
close OUT;

exit 0;


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
 
  # step 1, build sitelist

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
    push $hashRef->{'calc'}->{'sitelist'}, \@tmp;
    print "$tmp[0]  $tmp[1]  $tmp[2]  \n";
  }

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
  my @photon_files;
  my $nphoton;
  print "Looking for available photon files:\n";
  opendir DIR, "../" or die $!;
  while( my $file = readdir( DIR ) )
  {
    if( $file =~ m/^photon\d+$/ )
    {
      push @photon_files, $file;
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
        push @photon_files, $file;
      }
    }
    closedir DIR;
  }
  if( $#photon_files == -1 )
  {
    print "!!!  Could not find any photon files  !!!\n  Will have to quit :(\n";
    exit 1;
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
      print "        $_\n";
      copy( "../$_", $_ ) or die "$!";
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


  copy(  catfile( updir(), "PREP", "BSE", "eshift.ipt" ), "eshift.ipt" ) == 1 or die "Failed to copy eshift.ipt\n$!";
  print( "FIX ESHIFT\n");
  
}

sub CalcParams {
  my ($newRef, $commonRef, $oldRef ) = @_;

  $newBSEdata->{'calc'} = {};
  copyAndCompare( $newBSEdata->{'calc'}, $commonOceanData->{'calc'}, $bseData->{'calc'},
                  $newBSEdata->{'bse'}, [ 'mode' ] );

  if( $newBSEdata->{'calc'}->{'mode'} eq 'xas' || $newBSEdata->{'calc'}->{'mode'} eq 'xes' 
   || $newBSEdata->{'calc'}->{'mode'} eq 'rxs' ) {
    copyAndCompare( $newBSEdata->{'calc'}, $commonOceanData->{'calc'}, $bseData->{'calc'},
                  $newBSEdata->{'bse'}, [ 'edges' ] );
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
  } elsif ( $newRef->{'calc'}->{'mode'} eq 'val' ) {
    $type[0] = 'val';
  } elsif ( $newRef->{'calc'}->{'mode'} eq 'rxs' ) {
    $type[0] = 'core'; $type[1] = 'val';
  } else {
    die "Unrecongnized calculations type: " . $newRef->{'calc'}->{'mode'} . "\n";
  }

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
    push @fileList, catfile( updir(), "OPF", "zpawinfo", "melfile$zstring" );
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


sub writeRunlist {
  my ($hashRef, $photon_files ) = @_;

  my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

  my $runLength = 1;
  if( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'xes' ) {
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
    foreach my $way (@{$photon_files} ) {
      $way =~ m/(\d+)$/ or die "Malformed photon file name:\t$way\n";
      my $i = $1;
      print RUNLIST "$znum  $nnum  $lnum  $elname  ${nnum}$alphal{$lnum}  $elnum  $i  $run_text\n";
    }
  }

  close RUNLIST;
}


sub writeBSEin {
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
    } else { die "Failed to write gmres in writeBSEin\n"; }


  } else { die "Unsupported calc mode in writeBSEin\n"; }


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
    $nb = $hashRef->{'bse'}->{'brange'}[3] - $hashRef->{'bse'}->{'brange'}[2] + 1;
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
  if( $hashRef->{'calc'}->{'mode'} eq 'xas' || $hashRef->{'calc'}->{'mode'} eq 'rxs' ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

  if( $hashRef->{'calc'}->{'mode'} eq 'val' ) { die "Valence not programed yet\n"; }
  #TODO: this is redundant
  open OUT, ">", "mode" or die "$!";
  printf OUT "%g %i\n", $hashRef->{'bse'}->{'core'}->{'strength'}, 
                       $hashRef->{'bse'}->{'core'}->{'haydock'}->{'niter'};
  close OUT;
    

  open OUT, ">", "nelectron" or die "$!";
  print OUT $hashRef->{'bse'}->{'nelec'} ."\n";
  close OUT;

  open OUT, ">", "nedges" or die "$!";
  printf OUT "%i\n", scalar @{$hashRef->{'calc'}->{'edges'}};
  close OUT;

  open OUT, ">", "epsilon" or die "$!";
  printf OUT "%g\n", $hashRef->{'bse'}->{'epsilon'};
  close OUT;

  if( $hashRef->{'calc'}->{'mode'} eq 'val' ) { die "Valence not programed yet\n"; }
  open OUT, ">", "cnbse.write_rhs" or die "$!";
  if( $hashRef->{'bse'}->{'core'}->{'write_rhs'} ) { 
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

  if( $hashRef->{'calc'}->{'mode'} eq 'val' ) { die "Valence not programed yet\n"; }
  open OUT, ">", "haydockconv.ipt" or die;
  printf OUT "%g  %i\n", $hashRef->{'bse'}->{'core'}->{'haydock'}->{'converge'}->{'thresh'}, 
                         $hashRef->{'bse'}->{'core'}->{'haydock'}->{'converge'}->{'spacing'};
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
  print OUT ".false.\n";
  close OUT;
  die "Fix the metal flag!!\n";
}
