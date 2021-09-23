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
#JSON::PP->import;
use JSON::PP;
use Cwd 'abs_path';
use Cwd;
use File::Spec::Functions;
use Storable qw(dclone);
use Scalar::Util qw( looks_like_number ); 


###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/dft\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}

my $dir = getcwd;
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = abs_path( catdir( updir(), $dir ) ); }

my $json = JSON::PP->new;

# Load run info from Common
my $dataFile = "../Common/postDefaultsOceanDatafile";
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


# Grab previous run info if it exists
my $dftDataFile = "dft.json";
my $dftData;
if( -e $dftDataFile && open( my $in, "<", $dftDataFile ) )
{
  local $/ = undef;
  $dftData = $json->decode(<$in>);
  close($in);
}

# Build to-do list
my $newDftData;

# First we check, using the SCF flag to store result
# 1) Was previous run?
# 2) Structure matches
# 3) PSP matches
# 
# Any failures skip future tests and SCF to not done (which in turn, cascades to all futher runs)
$newDftData->{'scf'}->{'complete'} = JSON::PP::true;

# (was previous run done)
$newDftData->{'scf'}->{'complete'} = JSON::PP::false 
    unless( exists $dftData->{'scf'} && exists $dftData->{'scf'}->{'complete'} );

# (build the structure and check )
$newDftData->{'scf'}->{'complete'} = JSON::PP::false unless( exists $dftData->{'structure'} );

my @structureNumericalList = ( "typat", "xred", "znucl" );
my @structure2DNumericalList = ( "avecs" );
my @structureStringList = ( "zsymb" );

if( 0 ) {
foreach my $tag (@structureNumericalList)
{
  $newDftData->{'structure'}->{ $tag } = dclone $commonOceanData->{'structure'}->{ $tag };

  # Don't bother checking if we've already decided to re-run
  next unless( $newDftData->{'scf'}->{'complete'} );
  if( scalar @{$commonOceanData->{'structure'}->{ $tag }} != scalar @{$dftData->{'structure'}->{ $tag }} )
  {
    $newDftData->{'scf'}->{'complete'} = JSON::PP::false;
    next;
  }
  for( my $i=0; $i < scalar @{$commonOceanData->{'structure'}->{ $tag }}; $i++ )
  {
    unless( $commonOceanData->{'structure'}->{ $tag }[$i] == $dftData->{'structure'}->{ $tag }[$i] )
    {
      $newDftData->{'scf'}->{'complete'} = JSON::PP::false;
      last;
    }
  }
}
}

$newDftData->{'structure'} = {};

my @structureList = ( "typat", "xred", "znucl", "avecs", "zsymb" );
copyAndCompare( $newDftData->{'structure'}, $commonOceanData->{'structure'}, $dftData->{'structure'}, 
                $newDftData->{'scf'}, \@structureList );


my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);
open OUT, ">", "dft.json" or die;
print OUT $json->encode($newDftData);
close OUT;

my $condition = 1;
EXIT_IF: {
  if ($condition) {

     last EXIT_IF; # break from code block

     print "never get's executed\n";
  }
}


# 4) various convergence parameters match 

# Build SCF complete list from input

# Check to-do list againt done list (including what system was done)


#NOTE
# 1. Update done list (dft.json) after each individual step

# SCF

# SCF post-processing steps

# BSE final states

# Screen


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
    $newRef->{ $t } = dclone $commonRef->{ $t };
    
    next unless( $complete->{'complete'} );
    next unless( exists $oldRef->{ $t } );
    if( scalar @{$newRef->{ $t }} != scalar @{$oldRef->{ $t }} )
    {
      $complete->{'complete'} = JSON::PP::false;
      next;
    }
  
    my $x = $newRef->{ $t }[0];
    if( ref( @{$newRef->{ $t }}[0] ) eq 'ARRAY' )
    {
      $x = $newRef->{$t}[0][0];
    }
    if( looks_like_number( $x ) )
    {
      $comp = sub { $_[0] == $_[1] }; 
    }
    else
    {
      $comp = sub { $_[0] eq $_[1] };
    }
    
    # Should change this to recursive unwrapping subroutine for any order nests?
    COMPARE_LOOP:
    for( my $i=0; $i < scalar @{$newRef->{ $t }}; $i++ )
    {
      if( ref( @{$newRef->{ $t }}[$i] ) eq 'ARRAY' ) 
      {
        for( my $j = 0; $j < scalar @{@{$newRef->{ $t }}[$i]}; $j++ )
        {
          print "###  $newRef->{$t}[$i][$j]  $oldRef->{$t}[$i][$j]  ";#  \n";
#          unless( $newRef->{$t}[$i][$j] == $oldRef->{$t}[$i][$j] )
          unless( $comp->( $newRef->{$t}[$i][$j], $oldRef->{$t}[$i][$j] ) )
          {
            $complete->{'complete'} = JSON::PP::false;
            last COMPARE_LOOP;
          }
        }
      }
      else
      {
        print "###  $newRef->{$t}[$i]  $oldRef->{$t}[$i]  ";#  \n";
        print "\n";
#        unless( $newRef->{$t}[$i] == $oldRef->{$t}[$i] )
        unless( $comp->( $newRef->{$t}[$i], $oldRef->{$t}[$i] ) )
        {
          $complete->{'complete'} = JSON::PP::false;
          last COMPARE_LOOP;
        }
      }
    }
  }

}
