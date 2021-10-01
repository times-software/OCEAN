#!/usr/bin/perl
# Copyright (C) 2014, 2016 - 2019 OCEAN collaboration
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
use Storable qw(dclone);
use Scalar::Util qw( looks_like_number );

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/prep\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}

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

# Load run info from DFT
my $dataFile = catfile( updir(), "DFT", "dft.json" );
die "Failed to find $dataFile\n" unless( -e $dataFile );

my $dftData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $dftData = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $dataFile\n$!";
}

# Grab previous run info if it exists
my $prepDataFile = "prep.json";
my $prepData;
if( -e $prepDataFile && open( my $in, "<", $prepDataFile ) )
{
  local $/ = undef;
  $prepData = $json->decode(<$in>);
  close($in);
}

# Build to-do list
my $newPrepData;


$newPrepData->{'bse'} = {};
$newPrepData->{'bse'}->{'complete'} = JSON::PP::true;

# (was previous run done)
$newPrepData->{'bse'}->{'complete'} = JSON::PP::false
    unless( exists $prepData->{'bse'} && exists $prepData->{'bse'}->{'complete'} 
            && $prepData->{'bse'}->{'complete'});

copyAndCompare( $newPrepData->{'bse'}, $dftData->{'bse'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, [ 'hash' ] );


##### currently needed files
# TODO: make a single input file that can be read all at once by ocean_prep
# avecs
# ntype
# natoms
# znucl
# typat
# taulist
# coord = xred
my @list = ( "avecs", "typat", "znucl", "xred" );
copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'structure'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, \@list );

# kmesh.ipt
# nspin
# k0
# qinunitsofbv
# brange
# split
# xmesh
@list = ( "kmesh", "kshift", "xmesh" );
copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'bse'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, \@list );

copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'calc'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["nonzero_q", "photon_q"] );

copyAndCompare( $newPrepData->{'bse'}, $dftData->{'general'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["nspin"] );

copyAndCompare( $newPrepData->{'bse'}, $dftData->{'bse'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["split", "brange"] );


# Setting tmels and cks requires some logic
my $calc = $commonOceanData->{'calc'}->{'mode'};
my $tmels =  JSON::PP::false;
my $cks =  JSON::PP::false;

if( $calc eq 'val' ) {
  $tmels =  JSON::PP::true;
} else {
  $cks = JSON::PP::true;
  $tmels =  JSON::PP::true if ( $calc eq 'rxs' );
}


$newPrepData->{'bse'}->{'cks'} = JSON::PP::false;
$newPrepData->{'bse'}->{'tmels'} = JSON::PP::false;
# TODO, in the future make edges inclusive 
if( $cks )
{
  copyAndCompare( $newPrepData->{'bse'}, $commonOceanData->{'calc'}, $prepData->{'bse'},
                $newPrepData->{'bse'}, ["edges"] );
  $newPrepData->{'bse'}->{'complete'} = JSON::PP::false unless( $prepData->{'bse'}->{'cks'} );
  $newPrepData->{'bse'}->{'cks'} = JSON::PP::true;
}
if( $tmels )
{
  $newPrepData->{'bse'}->{'complete'} = JSON::PP::false unless( $prepData->{'bse'}->{'tmels'} );
  $newPrepData->{'bse'}->{'tmels'} = JSON::PP::true;
}

# prep.tmels
# prep.cks

# prep.u2  ## TODO: non-interacting spectra don't need u2, but do need cks or tmels

unless( $newPrepData->{'bse'}->{'complete'} )
{
######## file nonsense
  mkdir "BSE";
  chdir "BSE";

  writeOceanPrepInput( $newPrepData->{'bse'});


  chdir updir();
#######
  
  

  open OUT, ">", "prep.json" or die;
  print OUT $json->encode($newPrepData);
  close OUT;
}


exit 1;

sub writeOceanPrepInput
{
  my $hashRef = $_[0];
  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "%s  %s  %s\n", $hashRef->{'avecs'}[$i][0],
                                $hashRef->{'avecs'}[$i][1],
                                $hashRef->{'avecs'}[$i][2];

  }
  close OUT;

  open OUT, ">", "ntype" or die "$!";
  print OUT scalar @{$hashRef->{'znucl'}} . "\n";
  close OUT;

  open OUT, ">", "typat" or die "$!";
  for( my $i = 0; $i< scalar @{$hashRef->{'typat'}}; $i++ )
  {
    print OUT $hashRef->{'typat'}[$i] . "\n";
  }
  close OUT;

  open OUT, ">", "natoms" or die "$!";
  print OUT scalar @{$hashRef->{'typat'}} . "\n";
  close OUT;

  open OUT, ">", "znucl" or die "$!";
  for( my $i = 0; $i< scalar @{$hashRef->{'znucl'}}; $i++ )
  {
    print OUT $hashRef->{'znucl'}[$i] . "\n";
  }
  close OUT;

  open OUT, ">", "taulist" or die "$!";
  for( my $i = 0; $i< scalar @{$hashRef->{'xred'}}; $i++ )
  {
    printf  OUT "%s  %s  %s\n", $hashRef->{'xred'}[$i][0],
                                $hashRef->{'xred'}[$i][1],
                                $hashRef->{'xred'}[$i][2];
  }
  close OUT;

  open OUT, ">", "coord" or die $!;
  print OUT "xred\n";
  close OUT;


  open OUT, ">", "kmesh.ipt" or die $!;
  printf  OUT "%s  %s  %s\n", $hashRef->{'kmesh'}[0],
                              $hashRef->{'kmesh'}[1],
                              $hashRef->{'kmesh'}[2];
  close OUT;

  open OUT, ">", "k0.ipt" or die $!;
  printf  OUT "%s  %s  %s\n", $hashRef->{'kshift'}[0],
                              $hashRef->{'kshift'}[1],
                              $hashRef->{'kshift'}[2];
  close OUT;

  open OUT, ">", "qinunitsofbvectors.ipt" or die $!;
  printf  OUT "%s  %s  %s\n", $hashRef->{'photon_q'}[0],
                              $hashRef->{'photon_q'}[1],
                              $hashRef->{'photon_q'}[2];
  close OUT;

  open OUT, ">", "xmesh.ipt" or die $!;
  printf  OUT "%s  %s  %s\n", $hashRef->{'xmesh'}[0],
                              $hashRef->{'xmesh'}[1],
                              $hashRef->{'xmesh'}[2];
  close OUT;

  #TODO fix when nbands is less than DFT supply (
  open OUT, ">", "brange.ipt" or die $!;
  printf OUT "%i  %i  %i  %i\n", $hashRef->{'brange'}[0], $hashRef->{'brange'}[1],
                                 $hashRef->{'brange'}[2], $hashRef->{'brange'}[3];
  close OUT;

  open OUT, ">", "nspin" or die $!;
  print OUT $hashRef->{'nspin'} . "\n";
  close OUT;

  open OUT, ">", "dft.split" or die $!;
  if( $hashRef->{'split'} ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;


  open OUT, ">", "prep.tmels" or die $!;
  if( $hashRef->{'tmels'} ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

  open OUT, ">", "prep.cks" or die $!;
  if( $hashRef->{'cks'} ) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  } 
  close OUT;

  if( $hashRef->{'cks'} ) 
  {
    open IN, ">", "edges" or die $!;
#TODO: Fix edges, should be array of arrays?
#    for( my $i = 0; $i < scalar @{$hashRef->{'edges'}}; $i++ )
#    {
#      printf IN "%i  %i  %i\n", $hashRef->{'edges'}[$i][0], $hashRef->{'edges'}[$i][1], $hashRef->{'edges'}[$i][2];
#    }
    close IN;
  }
}
sub copyAndCompare
{ 
  my $newRef = $_[0];
  my $commonRef = $_[1];
  my $oldRef = $_[2];
  my $complete = $_[3];
  my @tags = @{$_[4]};
  
  my $comp;# = sub { $_[0] == $_[1] }; 
  
  foreach my $t (@tags)
  { 
    if( ref( $commonRef->{ $t } ) eq '' )
    {
#      print "$commonRef->{ $t } ---\n"; 
      $newRef->{ $t } = $commonRef->{ $t };
    }
    else
    { 
      $newRef->{ $t } = dclone $commonRef->{ $t };
    }
    
    next unless( $complete->{'complete'} );
    unless( exists $oldRef->{ $t } )
    { 
      $complete->{'complete'} = JSON::PP::false;
      next;
    }
    
    recursiveCompare( $newRef->{$t}, $oldRef->{$t}, $complete);
    print "DIFF:   $t\n" unless( $complete->{'complete'} );
  }

}


sub recursiveCompare
{
  my $newRef = $_[0];
  my $oldRef = $_[1];
  my $complete = $_[2];

  return unless( $complete->{'complete'} );


  if( ref( $newRef ) eq 'ARRAY' )
  {
    if( scalar @{ $newRef } != scalar @{ $oldRef } )
    {
      $complete->{'complete'} = JSON::PP::false;
      return;
    }
    for( my $i = 0; $i < scalar @{ $newRef }; $i++ )
    {
      recursiveCompare( @{$newRef}[$i], @{$oldRef}[$i], $complete );
      return unless( $complete->{'complete'} );
    }
  }
  else
  {
#    print "#!  $newRef  $oldRef\n";
    if( looks_like_number( $newRef ) )
    {
      $complete->{'complete'} = JSON::PP::false unless( $newRef == $oldRef );
    }
    else
    {
      $complete->{'complete'} = JSON::PP::false unless( $newRef eq $oldRef );
    }
  }
}
