#!/usr/bin/perl
# Copyright (C) 2010, 2013 - 2026 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Copy;
use File::Spec::Functions;
use File::Compare;
use Cwd 'abs_path';

use POSIX;

require JSON::PP;
use JSON::PP;
use Storable qw(dclone);
use Scalar::Util qw( looks_like_number );

use Data::Dumper;
use Time::HiRes qw( gettimeofday tv_interval );

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/screen\.pl/;
#   my $test = File::Spec->rel2abs( $0 );
#  my $test = abs_path( $1 );
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
#  print "OCEAN_BIN not set. Setting it to $1\n";
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################

my @timeSections = ( 'density', 'model', 'screen', 'combine' );


my $dataFile = catfile( updir(), "Common", "postDefaultsOceanDatafile" );
if( -e $dataFile )
{
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

  $dataFile = catfile( updir(), "DFT", "dft.json" );
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
  #TODO: add in PREP stage for legacy wave function file support here


  $dataFile = catfile( updir(), "OPF", "opf.json" );
  my $opfData;
  if( open( my $in, "<", $dataFile ))
  {
    local $/ = undef;
    $opfData = $json->decode(<$in>);
    close($in);
  } else {
    $opfData = {};
  }

  my $screenData;
  $dataFile = "screen.json";
  if( -e $dataFile && open( my $in, "<", $dataFile ) )
  {
    local $/ = undef;
    $screenData = $json->decode(<$in>);
    close($in);
  }



  # Quit if a normal valence calculation (no need for RPA screening)
  if( $commonOceanData->{'calc'}->{'mode'} eq 'val' )
  {
    unless( $commonOceanData->{'screen'}->{'mode'} eq 'grid' ) 
    {
      print "SCREEN stage not needed for valence w/ HLL screening\n";
      exit 0;
    }
  }


  #TODO:: cut this down and make it better
  my @ExtraFiles = ("specpnt", "Pquadrature", "hqp", "lqp", "EvenQuadHalf.txt" );
  foreach my $f (@ExtraFiles)
  {
    copy( catfile( $ENV{"OCEAN_BIN"}, $f ), $f );
  }





  my $newScreenData;

#  foreach my $sec (@timeSections) {
#    $newScreenData->{$sec}->{'time'} = $screenData->{$sec}->{'time'} if( exists $screenData->{$sec}->{'time'} );
#    if(  exists $screenData->{$sec}->{'time'} ) {
#      printf "%s %f\n", $sec, $screenData->{$sec}->{'time'};
#    } else {
#      printf "%s NULL\n", $sec;
#    }
#  }


  my $fake->{ 'complete' } = JSON::PP::false;
  $newScreenData->{'computer'} = {};
  copyAndCompare( $newScreenData->{'computer'}, $commonOceanData->{'computer'}, $screenData->{'computer'},
                  $fake, [ 'para_prefix' ] );

  $newScreenData->{'structure'} = {};
  my @structList = ( 'avecs', 'bvecs', 'xred', 'znucl', 'typat', 'elname', 'epsilon' );
  copyAndCompare( $newScreenData->{'structure'}, $dftData->{'structure'}, $screenData->{'structure'},
                  $fake, \@structList );

  $newScreenData->{'general'} = {};
  my $general->{'complete'} = JSON::PP::true;
  $newScreenData->{'general'}->{'complete'} = JSON::PP::true;
  copyAndCompare( $newScreenData->{'general'}, $commonOceanData->{'screen'},  $screenData->{'general'},
                  $general, [ 'mode' ] );
#  if( $newScreenData->{'general'}->{'mode'} eq 'core' ) {
    copyAndCompare( $newScreenData->{'general'}, $commonOceanData->{'calc'},  $screenData->{'general'},
                    $general, [ 'edges' ] );
    copyAndCompare( $newScreenData->{'general'}, $commonOceanData->{'bse'},  $screenData->{'general'},
                    $general, [ 'xmesh' ] );
    buildSiteEdgeWyck( $newScreenData->{'general'}, $newScreenData->{'structure'} );
#  } elsif( $newScreenData->{'general'}->{'mode'} eq 'grid' ) {

#  } else {
#    die "Malformed screening mode: $newScreenData->{'general'}->{'mode'}\n";
#  }
  
  copyAndCompare( $newScreenData->{'general'}, $dftData->{'general'}, $screenData->{'general'},
                  $general, ['program'] );

  buildAndCompareSitedir( $newScreenData->{'general'}, $screenData->{'general'} );

  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;


  #### Site/xmesh density average
  $newScreenData->{'density'} = {} unless( exists $newScreenData->{'density'});
  $newScreenData->{'density'}->{'complete'} = JSON::PP::true;
  $newScreenData->{'density'}->{'complete'} = JSON::PP::false
    unless( exists $screenData->{'density'}->{'complete'} && $screenData->{'density'}->{'complete'} );

  copyAndCompare( $newScreenData->{'density'}, $dftData->{'scf'}, $screenData->{'density'},
                  $newScreenData->{'density'}, [ 'hash' ] );

#  copyAndCompare( $newScreenData->{'density'}, $commonOceanData->{'screen'},  $screenData->{'density'},
#                  $newScreenData->{'density'}, [ 'mode' ] );
#  if( $newScreenData->{'density'}->{'mode'} eq 'core' )
#  {
#    copyAndCompare( $newScreenData->{'density'}, $commonOceanData->{'calc'},  $screenData->{'density'},
#                  $newScreenData->{'density'}, [ 'edges' ] )
#  } elsif( $newScreenData->{'density'}->{'mode'} eq 'grid' ) {
#    copyAndCompare( $newScreenData->{'density'}, $commonOceanData->{'bse'},  $screenData->{'density'},
#                  $newScreenData->{'density'}, [ 'xmesh' ] )
#
#  } else {
#    die "Malformed screening mode: $newScreenData->{'density'}->{'mode'}\n";
#  }

  $newScreenData->{'model'} = {} unless( exists $newScreenData->{'model'} );
  $newScreenData->{'model'}->{'complete'} = $newScreenData->{'density'}->{'complete'};
#  print "MODEL: " . $newScreenData->{'model'}->{'complete'} . "\n";
  $newScreenData->{'model'}->{'complete'} = JSON::PP::false
    unless( exists $screenData->{'model'}->{'complete'} && $screenData->{'model'}->{'complete'} );
#  print "MODEL: " . $newScreenData->{'model'}->{'complete'} . "\n";

  # This is necessary because we are replacing the 'model' section entirely
  $fake->{ 'complete' } = $newScreenData->{'model'}->{'complete'};
  copyAndCompare( $newScreenData, $commonOceanData->{'screen'}, $screenData, $fake, ['model'] );
  $newScreenData->{'model'}->{'complete'} = $fake->{'complete'};
  copyAndCompare( $newScreenData->{'model'}, $commonOceanData->{'screen'}, $screenData->{'model'},
                  $newScreenData->{'model'}, [ "shells" ] );
  copyAndCompare( $newScreenData->{'model'}, $newScreenData->{'structure'}, $screenData->{'model'},
                  $newScreenData->{'model'}, [ "epsilon" ] );
  $newScreenData->{'model'}->{'time'} = $screenData->{'model'}->{'time'} if( exists $screenData->{'model'}->{'time'});
  

  
#  if( $newScreenData->{'density'}->{'complete'} )
#  {
#    print "Skipping MPI_avg not enabled yet! The dev needs more coffee\n";
#    $newScreenData->{'density'}->{'complete'} = JSON::PP::false;
#  } 
  #### Site/xmesh density average

  #### SCREEN calculation
  $newScreenData->{'screen'}->{'complete'} = JSON::PP::true;
  $newScreenData->{'screen'}->{'complete'} = JSON::PP::false
    unless( exists $screenData->{'screen'}->{'complete'} && $screenData->{'screen'}->{'complete'} );
  print "No previous screening calculation detected\n" unless( $newScreenData->{'screen'}->{'complete'} );

  my @screenList = ( "all_augment", "augment", "convertstyle", "grid", "inversionstyle", "kmesh", 
                     "kshift", "mode", "nbands", "shells", "final", "vext" );
  copyAndCompare( $newScreenData->{'screen'}, $commonOceanData->{'screen'}, $screenData->{'screen'},
                  $newScreenData->{'screen'}, \@screenList );

  copyAndCompare( $newScreenData->{'screen'}, $dftData->{'screen'}, $screenData->{'screen'},
                  $newScreenData->{'screen'}, [ 'hash', 'isGamma', 'brange' ] );
  copyAndCompare( $newScreenData->{'screen'}, $commonOceanData->{'calc'},  $screenData->{'screen'},
                  $newScreenData->{'screen'}, [ 'edges' ] ) ;
  copyAndCompare( $newScreenData->{'screen'}, $dftData->{'scf'}, $screenData->{'screen'},
                  $newScreenData->{'screen'}, [ 'fermi' ] );
  copyAndCompare( $newScreenData->{'screen'}, $dftData->{'general'}, $screenData->{'screen'},
                  $newScreenData->{'screen'}, [ 'nspin' ] );

  #### check to make sure OPFs didn't re-run ####
#  print "OPF\n";
  my %Z;
  if( $newScreenData->{'general'}->{'mode'} eq 'core' ) {
    foreach my $znl (@{$newScreenData->{'general'}->{'edgelist'}}) {
      my @znl = split ' ', $znl;
      my $z = $znl[0]*1;
      $Z{$z} = 1;
    }
  } elsif( $newScreenData->{'general'}->{'mode'} eq 'grid' ) {
    foreach (@{$newScreenData->{'structure'}->{'znucl'}}) {
      $Z{$_} = 1;
    }
  }
  if( $newScreenData->{'screen'}->{'complete'} ) {
    if( exists $screenData->{'psp'} ) {
      foreach my $z ( keys %Z ) {
        unless( exists $opfData->{'completed'} ) {
          $newScreenData->{'screen'}->{'complete'} = JSON::PP::false;
          print "OPF/opf.json missing or incomplete. Will have to re-run screening. (A)\n";
          last;
        }
        unless( exists $opfData->{'completed'}->{$z} &&
                exists $opfData->{'completed'}->{$z}->{'input_hash'} &&
                exists $opfData->{'completed'}->{$z}->{'psp_hash'} ) {
          $newScreenData->{'screen'}->{'complete'} = JSON::PP::false;
          print "OPF/opf.json missing or incomplete. Will have to re-run screening. (B)\n";
          last;
        }
        if( $screenData->{'psp'}->{$z}->{'input_hash'} ne $opfData->{'completed'}->{$z}->{'input_hash'} ||
            $screenData->{'psp'}->{$z}->{'psp_hash'} ne $opfData->{'completed'}->{$z}->{'psp_hash'} ) {
          $newScreenData->{'screen'}->{'complete'} = JSON::PP::false;
          print "Psuedopotentials or OPF settings appear to have changed. Will have to re-run screening Z=$z\n";
          last;
        }
      }
    } else {
      $newScreenData->{'screen'}->{'complete'} = JSON::PP::false;
      print "SCREEN was done, but OPF info not cached. Likely a bug, but will re-run screening\n";
    }
    #Note, if someone mucks with any single psp input it'll invalidate all the screening runs regardless,
    # but also that is an edge case
      
  }
  $newScreenData->{'psp'} = {};
  foreach my $z ( keys %Z ) {
    $newScreenData->{'psp'}->{$z} = {};
    $newScreenData->{'psp'}->{$z}->{'input_hash'} = $opfData->{'completed'}->{$z}->{'input_hash'};
    $newScreenData->{'psp'}->{$z}->{'psp_hash'} = $opfData->{'completed'}->{$z}->{'psp_hash'};
  }
  

  #### RPA calculation

  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;

  writeExtraFiles( $newScreenData->{'structure'}, $newScreenData->{'screen'}, $newScreenData->{'general'} );
#                   $newScreenData->{'model'} );

  $newScreenData->{'combine'} = {} unless( exists $newScreenData->{'combine'} );
  $newScreenData->{'combine'}->{'complete'} = JSON::PP::true;
  unless( $newScreenData->{'model'}->{'complete'} && $newScreenData->{'screen'}->{'complete'} 
        && exists $screenData->{'combine'}->{'complete'} && $screenData->{'combine'}->{'complete'} )
  {
    $newScreenData->{'combine'}->{'complete'} = JSON::PP::false;
  }

  ##TODO:
#  unless( $newScreenData->{'combine'}->{'complete'} && $newScreenData->{'density'}->{'complete'} 
#        && $newScreenData->{'screen'}->{'complete'} )
#  {
#    print "Partial runs not enabled yet!\n";
#    $newScreenData->{'combine'}->{'complete'} = JSON::PP::false;
#    $newScreenData->{'density'}->{'complete'} = JSON::PP::false;
#    $newScreenData->{'screen'}->{'complete'} = JSON::PP::false;
#  }

#  $newScreenData->{'offset'} = {} unless ( exists $newScreenData->{'offset'} );
#  $newScreenData->{'offset'}->{'complete'} = JSON::PP::true;
#  unless( $newScreenData->{'density'}->{'complete'} && $newScreenData->{'screen'}->{'complete'} 
#        && exists $screenData->{'offset'}->{'complete'} && $screenData->{'offset'}->{'complete'} ) {
#    $newScreenData->{'offset'}->{'complete'} = JSON::PP::false;
#  }

  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;

  unless( $newScreenData->{'density'}->{'complete'} )
  {
    print "Running average\n";
    my $t0 = [gettimeofday];
    runRhoOfG( $newScreenData->{'model'} );
    my $errorCode = runDensityAverage( $newScreenData );
    if( $errorCode )
    { 
      die "Failed in runDensityAverage with code: $errorCode\n";
    }

    cleanAverage( $newScreenData->{'general'} );
    $newScreenData->{'density'}->{'complete'} = JSON::PP::true;
    $newScreenData->{'density'}->{'time'} = tv_interval( $t0 );
  } 
  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;

  unless( $newScreenData->{'model'}->{'complete'} )
  {
    my $t0 = [gettimeofday];
#    print "MODEL: " . $newScreenData->{'model'}->{'complete'} . "\n";
    if( lc($newScreenData->{'model'}->{'flavor'}) eq 'sll') {
      runVhommod($newScreenData->{'model'}->{'SLL'}) 
    } else {
      die "Unrecognized model flavor!\n";
    }

    $newScreenData->{'model'}->{'complete'} = JSON::PP::true;
    $newScreenData->{'model'}->{'time'} = tv_interval( $t0 );
  } 
  
  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;


  unless( $newScreenData->{'screen'}->{'complete'} )
  {
    my $t0 = [gettimeofday];
    print "SCREEN\n";
    grabOPF();

    $newScreenData->{'screen'}->{'grid2'} = dclone $newScreenData->{'screen'}->{'grid'};
    buildMKRB( $newScreenData->{'screen'}, $newScreenData->{'general'} );
    grabAngFile( $newScreenData->{'screen'} );
    prepWvfn(  $newScreenData->{'screen'},  $newScreenData->{'general'} );

    runScreenDriver( $newScreenData->{'computer'} );

    $newScreenData->{'screen'}->{'complete'} = JSON::PP::true;


    cleanScreen( $newScreenData->{'general'}, $newScreenData->{'screen'} );
    $newScreenData->{'screen'}->{'time'} = tv_interval( $t0 );

    open OUT, ">", "screen.json" or die;
    print OUT $json->encode($newScreenData);
    close OUT;
  } else {
    $newScreenData->{'screen'}->{'grid2'} = dclone $screenData->{'screen'}->{'grid2'};
    open OUT, ">", "screen.json" or die;
    print OUT $json->encode($newScreenData);
    close OUT;
  }

  unless( $newScreenData->{'combine'}->{'complete'} )
  {
    my $t0 = [gettimeofday];
    print "COMBINE\n";

#    cleanScreen( $newScreenData->{'general'}, $newScreenData->{'screen'} );
    finishCorePotentials( $newScreenData->{'general'}, $newScreenData->{'screen'} );


    $newScreenData->{'combine'}->{'complete'} = JSON::PP::true;
    $newScreenData->{'combine'}->{'time'} = tv_interval( $t0 );
    open OUT, ">", "screen.json" or die;
    print OUT $json->encode($newScreenData);
    close OUT;

  }

  foreach my $sec (@timeSections) {
    printf "Time %s: %f\n", $sec, $newScreenData->{$sec}->{'time'};
    $newScreenData->{'time'} += $newScreenData->{$sec}->{'time'};
  }
  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;

  exit 0;
}

print "Non-json version no longer supported\n";
exit 1;


# Currently WIP
sub interp
{
  my ( $xInRef, $yInRef, $xOutRef, $yOutRef ) = @_;

  my $xstart = 0;
  my $xstop;
  my $xmid;
  for ( my $i=0; $i < scalar @{ $xOutRef }; $i++ )
  {
    $xstop = scalar @{ $xInRef };
    
    while( $xstop - $xstart > 2 )
    {
      $xmid = floor($xstart + $xstop ) / 2;
      if( ${ $xInRef }[$xmid] < ${ $xOutRef }[$i] )
      {
        $xstart = $xmid;
      }
      else #( ${ $xInRef }[$xmid] > ${ $xOutRef }[$i] )
      {
        $xstop = $xmid;
      }
    }
    print "${ $xInRef }[$xstart]\t${ $xOutRef }[$i]\t${ $xInRef }[$xstop]\n";
    my $run = ${ $xOutRef }[$i] - ${ $xInRef }[$xstart];
    ${ $yOutRef }[$i] = ${ $yInRef }[$xstart] 
                      + (${ $yInRef }[$xstop] - ${ $yInRef }[$xstart] )
                        / (${ $xInRef }[$xstop] - ${ $xInRef }[$xstart] ) * $run;
  }
}

sub buildSiteEdgeWyck
{
  my ($genRef, $structRef) = @_;

  $genRef->{'sitelist'} = [];
  $genRef->{'edgelist'} = [];
  $genRef->{'wyck'} = [];
  $genRef->{'fulllist'} = [];

  my %countByName;
  my @index;
  for( my $i=0; $i< scalar @{$structRef->{'typat'}}; $i++ )
  {
    my $s = sprintf " %s %20.16f %20.16f %20.16f", $structRef->{'elname'}[$structRef->{'typat'}[$i]-1], 
                  $structRef->{'xred'}[$i][0], $structRef->{'xred'}[$i][1], $structRef->{'xred'}[$i][2];
    push @{$genRef->{'wyck'}}, $s;
    if( exists $countByName{$structRef->{'elname'}[$structRef->{'typat'}[$i]-1]} )
    {
      $countByName{$structRef->{'elname'}[$structRef->{'typat'}[$i]-1]}++;
    } else
    {  
      $countByName{$structRef->{'elname'}[$structRef->{'typat'}[$i]-1]}=1;
    }
    push @index, $countByName{$structRef->{'elname'}[$structRef->{'typat'}[$i]-1]};
#    print "$structRef->{'typat'}[$i]  $structRef->{'elname'}[$structRef->{'typat'}[$i]-1]  $countByName{$structRef->{'elname'}[$structRef->{'typat'}[$i]-1]}\n";
#    print $index[$i] . "\n";
  }

  my %edges;
  my %sites;
  foreach( @{$genRef->{'edges'}} )
  {
    $_ =~ m/^(\d+)\s+(\d+)\s+(\d+)/;
    my ($a, $n, $l ) = ($1, $2, $3 );
    $sites{ $1 } = $structRef->{'znucl'}[$structRef->{'typat'}[$1-1]-1];
    my $s = sprintf "%3i %2i %2i", $structRef->{'znucl'}[$structRef->{'typat'}[$1-1]-1], $2, $3;
    $edges{ $s } = 1;

    $s = sprintf "%s %i %i %i %i", $structRef->{'elname'}[$structRef->{'typat'}[$a-1]-1], 
         $structRef->{'znucl'}[$structRef->{'typat'}[$a-1]-1], $index[$a-1], 
         $n, $l;
    push @{$genRef->{'fulllist'}}, $s;
        
  }
  foreach  my $key ( sort { $a <=> $b } keys %sites )
  {
#    push $genRef->{'sitelist'}, sprintf "%s    %s   %s %s %s",  $structRef->{'sites'}[$key-1],
#    push $genRef->{'sitelist'}, sprintf "%s    %s   %s %s %s", $index[$key-1], 
#            $sites{ $key }, $structRef->{'xred'}[$key-1][0],
#            $structRef->{'xred'}[$key-1][1], $structRef->{'xred'}[$key-1][2];
    push @{$genRef->{'sitelist'}}, sprintf "%s %i %i", $structRef->{'elname'}[$structRef->{'typat'}[$key-1]-1],
            $structRef->{'znucl'}[$structRef->{'typat'}[$key-1]-1], $index[$key-1];
#    print "$key\n";
  }

  my %zee;
  foreach  my $key ( sort { $a <=> $b } keys %edges )
  {
    push @{$genRef->{'edgelist'}}, $key;
    $key =~ m/^\s+(\d+)/;
    $zee{$1} = 1;
  }


  ######
  open OUT, ">", "sitelist" or die $!;
  print OUT scalar @{$genRef->{'sitelist'}};
  print OUT "\n";
  foreach (@{$genRef->{'sitelist'}}) {
    print OUT $_ . "\n";
  }
  close OUT;

  open OUT, ">", "edgelist" or die $!;
  foreach (@{$genRef->{'edgelist'}}) {
    print OUT $_ . "\n";
  }
  close OUT;

  open OUT, ">", "xyz.wyck" or die $!;
  print OUT scalar @{$genRef->{'wyck'}};
  print OUT "\n";
  foreach (@{$genRef->{'wyck'}}) {
    print OUT $_ . "\n";
  }
  close OUT;

  open OUT, ">", "zeelist" or die $!;
  print OUT ( scalar keys %zee ) . "\n";
  foreach my $zee (keys %zee )
  {
    print OUT "$zee\n";
  }
  close OUT;
}

sub runDensityAverage
{
  my ($hashRef) = @_;

#  print "!!!!! Move file handling to correct order!!!!!\n";
#  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
#  for( my $i = 0; $i < 3; $i++ )
#  {
#    printf  OUT "%s  %s  %s\n", $hashRef->{'structure'}->{'avecs'}[$i][0],
#                                $hashRef->{'structure'}->{'avecs'}[$i][1],
#                                $hashRef->{'structure'}->{'avecs'}[$i][2];
#
#  }
#  close OUT;
#
#  open OUT, ">", "bvecs" or die "Failed to open bvecs\n$!";
#  for( my $i = 0; $i < 3; $i++ )
#  {
#    printf  OUT "%s  %s  %s\n", $hashRef->{'structure'}->{'bvecs'}[$i][0],
#                                $hashRef->{'structure'}->{'bvecs'}[$i][1],
#                                $hashRef->{'structure'}->{'bvecs'}[$i][2];
#
#  }
#  close OUT;

  open OUT, ">", "screen.mode" or die $!;
  print OUT $hashRef->{'general'}->{'mode'} . "\n";
  close OUT;


  if( $hashRef->{'general'}->{'mode'} eq 'core' ) {
#    my %sites;
#    foreach( @{$hashRef->{'density'}->{'edges'}} )
#    {
#      $_ =~ m/^(\d+)/;
#      $sites{ $1 } = 1;
#    }
#    open OUT, ">", "sitelist.new" or die $!;
#    my $n = keys %sites;
#    print OUT "$n\n";
#    foreach  my $key ( sort { $a <=> $b } keys %sites )
#    {
#      printf OUT "%s    %s %s %s\n",  $hashRef->{'structure'}->{'sites'}[$key-1], 
#            $hashRef->{'structure'}->{'xred'}[$key-1][0], 
#            $hashRef->{'structure'}->{'xred'}[$key-1][1], $hashRef->{'structure'}->{'xred'}[$key-1][2];
#    }
#    close OUT;
    
  } elsif ( $hashRef->{'general'}->{'mode'} eq 'grid' ) {
    open OUT, ">", "xmesh.ipt" or die $!;
    printf OUT "%i  %i  %i\n", $hashRef->{'general'}->{'xmesh'}[0], 
            $hashRef->{'general'}->{'xmesh'}[1], $hashRef->{'general'}->{'xmesh'}[2];
    close OUT;
  }
  else {
    return 101;
  }

#  copy( catfile( updir(), "DFT", "rhoofg" ), "rhoofg" );

  print "$hashRef->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1\n";
  system("$hashRef->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1" );
  if ($? == -1) {
      print "failed to execute: $!\n";
      die;
  }
  elsif ($? & 127) {
      printf "mpi_avg died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
      die;
  }
  else {
    my $errorCode = $? >> 8;
    if( $errorCode != 0 ) {
      die "CALCULATION FAILED\n  mpi_avg exited with value $errorCode\n";
    }
    else {
      printf "mpi_avg exited successfully with value %d\n", $errorCode;
    }
  }
  

  return 0;
}

sub buildMKRB
{
  my ($hashRef, $genRef) = @_;

  my $n = 0;
  $n = scalar @{$hashRef->{'grid2'}->{'ang'}} if( scalar @{$hashRef->{'grid2'}->{'ang'}} > $n );
  $n = scalar @{$hashRef->{'grid2'}->{'deltar'}} if( scalar @{$hashRef->{'grid2'}->{'deltar'}} > $n );
  $n = scalar @{$hashRef->{'grid2'}->{'rmode'}} if( scalar @{$hashRef->{'grid2'}->{'rmode'}} > $n );
  $n = scalar @{$hashRef->{'grid2'}->{'shells'}} if( scalar @{$hashRef->{'grid2'}->{'shells'}} > $n );

  #TODO: Move this into screen_driver to support multiple Z's in a single run
  if( $hashRef->{'grid2'}->{'shells'}[0] <= 0 )
  {
    if( $hashRef->{'screen'}->{'mode'} eq 'grid2' )
    {
      print "Warning. Negative value in screen.grid.shells doesn't make sense with valence grid calculation\n";
      $hashRef->{'grid2'}->{'shells'}[0] = 2 if( $hashRef->{'grid2'}->{'shells'}[1] > 2 );
      $hashRef->{'grid2'}->{'shells'}[0] = $hashRef->{'grid2'}->{'shells'}[1]/2
        if( $hashRef->{'grid2'}->{'shells'}[1] <= 2 );
      printf "   Used default of %.2f\n", $hashRef->{'grid2'}->{'shells'}[0];
    }
    else
    {
      $genRef->{'edgelist'}[0] =~ m/^\s+(\d+)/ or die;
      my $zeeName = sprintf "radfilez%03i", $1;
#      print catfile( "zpawinfo", $zeeName ) . "\n";
      open IN, "<", catfile( "zpawinfo", $zeeName ) or die $!;
      <IN> =~ m/^\s+(\S+)/ or die;
      $hashRef->{'grid2'}->{'shells'}[0] = $1*1;
    }
  }

  unless( $hashRef->{'grid2'}->{'rmode'}[0] =~ m/legendre/ || $hashRef->{'grid2'}->{'rmode'}[0] =~ m/uniform/ )
  {
    $hashRef->{'grid2'}->{'rmode'}[0] = 'legendre';
  }

  $hashRef->{'grid2'}->{'ang'}[0] = 5 if( $hashRef->{'grid2'}->{'ang'} < 5 );
  $hashRef->{'grid2'}->{'deltar'}[0] = 0.2 if( $hashRef->{'grid2'}->{'deltar'} < 0 );

  for( my $i = 1; $i < $n; $i++ )
  {
    $hashRef->{'grid2'}->{'ang'}[$i] = $hashRef->{'grid2'}->{'ang'}[$i-1] 
        if( scalar @{$hashRef->{'grid2'}->{'ang'}} <= $i );
    $hashRef->{'grid2'}->{'ang'}[$i] = 7
        if( $hashRef->{'grid2'}->{'ang'}[$i] < 0 || ! $hashRef->{'grid2'}->{'ang'}[$i] =~ m/\d/ );
    $hashRef->{'grid2'}->{'deltar'}[$i] = $hashRef->{'grid2'}->{'deltar'}[$i-1]
        if( scalar @{$hashRef->{'grid2'}->{'deltar'}} <= $i || $hashRef->{'grid2'}->{'deltar'}[$i] <= 0 );
    $hashRef->{'grid2'}->{'rmode'}[$i] = 'uniform'
        unless( $hashRef->{'grid2'}->{'rmode'}[$i] =~ m/uniform/ || $hashRef->{'grid2'}->{'rmode'}[$i] =~ m/legendre/);

    if( scalar @{$hashRef->{'grid2'}->{'shells'}} <= $i || $hashRef->{'grid2'}->{'shells'}[$i] < 0
        || $hashRef->{'grid2'}->{'shells'}[$i] 
              >= ( $hashRef->{'grid2'}->{'rmax'} - $hashRef->{'grid2'}->{'deltar'}[$i]/2) ) 
    {
      print "Ran out of shells, might be making fewer than expected!\n";
      $hashRef->{'grid2'}->{'shells'}[$i] = $hashRef->{'grid2'}->{'rmax'};
      $n = $i+1;
      last;
    }
  }
  while( scalar @{$hashRef->{'grid2'}->{'shells'}} > $n ) { pop @{$hashRef->{'grid2'}->{'shells'}} }
  while( scalar @{$hashRef->{'grid2'}->{'ang'}} > $n ) { pop @{$hashRef->{'grid2'}->{'ang'}} }
  while( scalar @{$hashRef->{'grid2'}->{'deltar'}} > $n ) { pop @{$hashRef->{'grid2'}->{'deltar'}} }
  while( scalar @{$hashRef->{'grid2'}->{'ang'}} > $n ) { pop @{$hashRef->{'grid2'}->{'ang'}} }

  open OUT, ">", "mkrb_control" or die "Failed to open mkrb_control for writing\n$!";
  printf OUT "%.16g  %i\n", $hashRef->{'grid2'}->{'rmax'}, $n;
  for( my $i = 0; $i < $n; $i++ )
  {
    printf OUT "%s %.16g %.16g %i specpnt\n", $hashRef->{'grid2'}->{'rmode'}[$i], 
              $hashRef->{'grid2'}->{'shells'}[$i], $hashRef->{'grid2'}->{'deltar'}[$i], 
              $hashRef->{'grid2'}->{'ang'}[$i];
#    print OUT "$hashRef->{'grid2'}->{'rmode'}[$i] $hashRef->{'grid2'}->{'shells'}[$i] "
#            . "$hashRef->{'grid2'}->{'deltar'}[$i] $hashRef->{'grid2'}->{'ang'}[$i] specpnt\n";
  }
  close OUT;

  return 0;
}

sub grabAngFile
{
  my ($hashRef) = @_;

  my %ang;
  foreach( @{$hashRef->{'grid2'}->{'ang'}} )
  {
    $ang{ $_ } = 1;
  }
  foreach my $key (keys %ang)
  {
    copy( catfile( $ENV{'OCEAN_BIN'}, "specpnt.$key"), "specpnt.$key" ) or die "Failed to get specpnt.$key\n";
  }

  return 0;
}

sub prepWvfn
{
  my ($screenHash, $genHash ) = @_;

  my $dirname = sprintf "k%i_%i_%iq%f_%f_%f", $screenHash->{'kmesh'}[0], $screenHash->{'kmesh'}[1],
                $screenHash->{'kmesh'}[2], $screenHash->{'kshift'}[0],
                $screenHash->{'kshift'}[1], $screenHash->{'kshift'}[2];

  $dirname = catdir( updir(), "DFT", $dirname );

  if( $genHash->{'program'} eq 'qe' )
  {
    unlink "Out" if( -l "Out" );
    symlink( catdir( $dirname, "Out" ), "Out" );
    if( $screenHash->{'isGamma'} )
    {
      open OUT, ">", "gamma" or die $!;
      print OUT "T\n";
      close OUT;
    }
    copy( catfile( $dirname, "wvfn.ipt" ), "wvfn.ipt" );
    copy( catfile( $dirname, "QE_EIGS.txt" ), "QE_EIGS.txt" );

    open OUT, ">", "prefix" or die $!;
    print OUT "system\n";
    close OUT;
  }
  else
  {
    unlink "RUN0001_WFK" if( -l "RUN0001_WFK" || -e "RUN0001_WFK" );
    symlink( catfile( $dirname, "NSCFx_WFK"), "RUN0001_WFK" ); 
    open OUT, ">", "wvfn.ipt";
    print OUT "abinit\n";
    close OUT;
  }
}

sub grabOPF
{
  symlink( catdir( updir(), "OPF", "zpawinfo" ), "zpawinfo" )
}

sub writeExtraFiles
{
  my ($structureRef, $screenRef, $genRef) = @_;

  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "%.16g  %.16g  %.16g\n", $structureRef->{'avecs'}[$i][0],
                                $structureRef->{'avecs'}[$i][1],
                                $structureRef->{'avecs'}[$i][2];

  }
  close OUT;

  open OUT, ">", "bvecs" or die "Failed to open bvecs\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "%.16g  %.16g  %.16g\n", $structureRef->{'bvecs'}[$i][0],
                                $structureRef->{'bvecs'}[$i][1],
                                $structureRef->{'bvecs'}[$i][2];

  }
  close OUT;
  
  open OUT, ">", "xmesh.ipt" or die "Failed to open xmesh.ipt\n$!";
  printf OUT "%i  %i  %i\n", $genRef->{'xmesh'}[0],
                                $genRef->{'xmesh'}[1],
                                $genRef->{'xmesh'}[2];
  close OUT;

  open OUT, ">", "brange.ipt" or die $!;
  printf OUT "%i  %i\n%i  %i\n", $screenRef->{'brange'}[0], $screenRef->{'brange'}[1],
                                 $screenRef->{'brange'}[2], $screenRef->{'brange'}[3];
  close OUT;

  open OUT, ">", "bands.ipt" or die $!;
  printf OUT "1  %i\n", $screenRef->{'nbands'};
  close OUT;

  open OUT, ">", "k0.ipt" or die $!;
  printf OUT "%f %f %f\n", $screenRef->{'kshift'}[0], $screenRef->{'kshift'}[1], $screenRef->{'kshift'}[2];
  close OUT;

  open OUT, ">", "kmesh.ipt" or die $!;
  printf OUT "%i %i %i\n", $screenRef->{'kmesh'}[0], $screenRef->{'kmesh'}[1], $screenRef->{'kmesh'}[2];
  close OUT;

  open OUT, ">", "nspin" or die $!;
  printf OUT "%i\n", $screenRef->{'nspin'};
  close OUT;
  
  open OUT, ">", "screen.augment" or die $!;
  if( $screenRef->{'augment'}) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

  open OUT, ">", "screen.mode" or die $!;
  print OUT $genRef->{'mode'} . "\n";
  close OUT;

  open OUT, ">", "screen.inversionstyle" or die $!;
  print OUT ($screenRef->{'inversionstyle'}) ."\n";
  close OUT;

  open OUT, ">", "screen.convertstyle" or die $!;
  print OUT ($screenRef->{'convertstyle'}) ."\n";
  close OUT;

  open OUT, ">", 'screen.allaug' or die $!;
  if( $screenRef->{'all_augment'}) {
    print OUT ".true.\n";
  } else {
    print OUT ".false.\n";
  }
  close OUT;

  open OUT, ">", 'screen.lmax' or die $!;
  printf OUT "%i\n", $screenRef->{'grid'}->{'lmax'};
  close OUT;

  open OUT, ">", 'shells' or die $!;
  printf OUT "%i\n", scalar @{$screenRef->{'shells'}};
  foreach ( @{$screenRef->{'shells'}} )
  {
    printf OUT "%f\n", $_;
  }
  close OUT;

  open OUT, ">", "efermiinrydberg.ipt" or die $!;
  print OUT ( $screenRef->{'fermi'} * 2 ) . "\n";
  close OUT;

  open OUT, ">", "screen.final.rmax" or die $!;
  print OUT ($screenRef->{'final'}->{'rmax'}) . "\n";
  close OUT;

  open OUT, ">", "screen.final.dr" or die $!;
  print OUT ($screenRef->{'final'}->{'dr'}) . "\n";
  close OUT;

  open OUT, ">", "screen.vext" or die $!;
  print OUT $screenRef->{'vext'} . "\n";
  close OUT;

  open OUT, ">", "epsilon" or die $!;
  print OUT ($structureRef->{'epsilon'}) . "\n";
  close OUT;

}

sub runRhoOfG
{
  my ($modelRef) = @_;

  my $rhofile = catfile( updir(), "DFT", "val.rhoofr" );
  if( (not -e $rhofile) || $modelRef->{'SLL'}->{'semicore_density'} ) {
    $rhofile = catfile( updir(), "DFT", "rhoofr" );
  }

  copy( $rhofile , "rhoofr" ) or die $!;
  copy( catfile( updir(), "DFT", "nfft" ), "nfft" ) or die $!;
  
  system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0  or die "Failed to run rhoofg.x\n";
  system("wc -l rhoG2 > rhoofg") == 0 or die "$!\n";
  system("sort -n -k 6 rhoG2 >> rhoofg") == 0 or die "$!\n";
  unlink( "rhoG2" );
  
}

sub runScreenDriver
{
  my $ref = $_[0];

  print "$ref->{'para_prefix'} $ENV{'OCEAN_BIN'}/screen_driver.x > screen_driver.log 2>&1\n";
  system("$ref->{'para_prefix'} $ENV{'OCEAN_BIN'}/screen_driver.x > screen_driver.log 2>&1\n");
  if ($? == -1) {
      print "failed to execute: $!\n";
      die;
  }
  elsif ($? & 127) {
      printf "screen_driver died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
      die;
  }
  else {
    my $errorCode = $? >> 8;
    if( $errorCode != 0 ) {
      die "CALCULATION FAILED\n  screen_driver exited with value $errorCode\n";
    }
    else {
      printf "screen_driver exited successfully with value %d\n", $errorCode;
    }
  }

}

sub runVhommod
{
  my $ref = $_[0];
  open OUT, ">", 'screen.model.dq' or die;
  printf OUT "%g\n", $ref->{'dq'};
  close OUT;


  open OUT, ">", 'screen.model.qmax' or die;
  printf OUT "%g\n", $ref->{'qmax'};
  close OUT;

  if( exists $ref->{'nav'} && $ref->{'nav'} >= 0 ) {
    open OUT, ">", 'fake_nav.ipt' or die;
    printf OUT "%g\n", $ref->{'nav'};
    close OUT;
  } elsif( -e 'fake_nav.ipt' ) {
    unlink 'fake_nav.ipt';
  }

  print "$ENV{'OCEAN_BIN'}/vhommod.x\n";
  system("$ENV{'OCEAN_BIN'}/vhommod.x" );
  if ($? == -1) {
      print "failed to execute: $!\n";
      die;
  }
  elsif ($? & 127) {
      printf "vhommod died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
      die;
  }
  else {
    my $errorCode = $? >> 8;
    if( $errorCode != 0 ) {
      die "CALCULATION FAILED\n  vhommod exited with value $errorCode\n";
    }
    else {
      printf "vhommod exited successfully with value %d\n", $errorCode;
    }
  }
}

sub buildAndCompareSitedir {
  my ($genHash, $oldGenHash) = @_;
  
  my @sitelist;

  if( $genHash->{'mode'} eq 'core' ) {
    foreach my $s ( @{$genHash->{'sitelist'}} ) {
      $s =~ m/(\S+)\s+\d+\s+(\d+)/ or die;
      push @sitelist, sprintf("z%2s%04i", $1, $2 );
    }
  } elsif( $genHash->{'mode'} eq 'grid' ) {
    my $nx = $genHash->{'xmesh'}[0]*$genHash->{'xmesh'}[1]*$genHash->{'xmesh'}[2];
    for( my $i = 1; $i <= $nx; $i++ ) {
      push @sitelist, sprintf( "x%06i", $i );
    }
  } else {
    die "Bad value of mode in general: $genHash->{'mode'}\n";
  }

  for( my $i = 0; $i < scalar @sitelist; $i++ ) {
    my $d = $sitelist[$i];
    unless( -d $d ) {
      mkdir $d or die "Failed to make directory $d\n$!";
    }
  }

  my $tmpHash->{'sitedir'} = [@sitelist];

  copyAndCompare( $genHash, $tmpHash, $oldGenHash, $genHash, ['sitedir'] );
}

sub cleanAverage {
  my ($genHash) = @_;

#  for( my $i = 0; $i < scalar @{$genHash->{'sitedir'}}; $i++ ) {
#    my $f = substr $sitelist[$i], 1;
  foreach my $site (@{$genHash->{'sitedir'}}) {
    my $f = substr $site, 1;
    move( "avg".$f, catfile( $site, "avg" ) ) or die "$!";
  }
}

#TODO hoist avg out of here and put it in cleanAverage routine
# TODO modify vhommod to look for chi in subfolder
#TODO make clean screen_driver routine to store chi0, grid, and chi
#TODO change sitename to just trim first character off of sitelist
sub cleanScreen
{
  my ($genHash, $screenHash) = @_;

  my @siteFiles = ("chi0", "grid", "chi", "avg" );
  my @potFiles = ( "vind", "vind0", "nind", "nind0" );

  my @sitelist;
  my %corelist;
  my @sitename;

  # Build sitelist for core or val
  if( $genHash->{'mode'} eq 'core' ) {
    foreach my $s ( @{$genHash->{'sitelist'}} ) {
      $s =~ m/(\S+)\s+\d+\s+(\d+)/ or die;
      push @sitelist, sprintf("z%2s%04i", $1, $2 );
      push @sitename, sprintf( "%2s%04i", $1, $2 );
    }
    foreach my $s (@{$genHash->{'fulllist'}}) {
      $s =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ or die;
      my $site = sprintf( "z%2s%04i", $1, $3 );
      my $edge = sprintf( "n%02il%02i", $4, $5) ;
      my $entry = sprintf "%3i %2i %2i", $2, $4, $5;
      $corelist{$site} = [] unless( exists $corelist{ $site} ) ;
      push @{$corelist{ $site}}, [$edge,$entry];
#      push @sitelist, catdir( sprintf( "z%2s%04i", $1, $2 ), sprintf( "n%02il%02i", $3, $4) );
    }
  }
  elsif( $genHash->{'mode'} eq 'grid' ) {
    my $nx = $genHash->{'xmesh'}[0]*$genHash->{'xmesh'}[1]*$genHash->{'xmesh'}[2];
    for( my $i = 1; $i <= $nx; $i++ ) {
      push @sitelist, sprintf( "x%06i", $i );
      push @sitename, sprintf( "%06i", $i );
    }
  }
  else {
    die "Bad value of mode in general: $genHash->{'mode'}\n";
  } 


  for( my $i = 0; $i < scalar @sitelist; $i++ ) {
    my $d = $sitelist[$i];
    unless( -d $d ) {
      mkdir $d or die "Failed to make directory $d\n$!";
    }

    foreach my $f (@siteFiles) {
       move( $f.$sitename[$i], catfile( $d, $f ) );
    }

    #TODO: decouple vind, etc from creation of chi & chi0
    #
    if( $genHash->{'mode'} eq 'core' ) {
      foreach my $e (@{$corelist{$d}} ) {
        my @e = @$e;
        unless( -d catdir( $d, $e[0] ) ) {
          mkdir catdir( $d, $e[0] );
        }
        foreach my $r (@{$screenHash->{'shells'}}) {
          my $radNam = sprintf( "zR%.2f", $r );
          my $radDir = catdir( $d, $e[0], $radNam );
          print sprintf("reopt%s.R%.2f", $sitename[$i], $r ) . "\n";
          move( sprintf("reopt%s.R%.2f", $sitename[$i], $r ), catfile( $d, sprintf("reopt.R%.2f",$r) ) );
          unless( -d $radDir ) {
            mkdir $radDir or die $!;
          }
          foreach my $f (@potFiles) {
#            copy( sprintf( "%s_%s.%.zR2f", $sitename[$i], $e, $f ), catfile( 
            print $sitelist[$i].'_'.$e[0].'.'.$radNam.$f . "\n";
            move( $sitelist[$i].'_'.$e[0].'.'.$radNam.$f, catfile( $radDir, $f ) );
          }
        }        
      }

#      finishCorePotentials( $genHash, $screenHash, \%corelist );
    } else {
      foreach my $r (@{$screenHash->{'shells'}}) {
        my $radNam = sprintf( "zR%.2f", $r );
        my $radDir = catdir( $d, $radNam );
        unless( -d $radDir ) {
          mkdir $radDir or die $!;
        }
        foreach my $f (@potFiles) {
          print $sitelist[$i].'.'.$radNam.$f . "\n";
          move( $sitelist[$i].'.'.$radNam.$f, catfile( $radDir, $f ) );
        }
      }  
    }
  }

  
}

sub finishCorePotentials
{
  my ($genHash, $screenHash ) = @_;

  #TODO check on valence paths
  return unless( $genHash->{'mode'} eq 'core' );
  my %completeList;
  foreach my $s (@{$genHash->{'fulllist'}}) {
    $s =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ or die;
    my $site = sprintf( "z%2s%04i", $1, $3 ); 
    my $edge = sprintf( "n%02il%02i", $4, $5) ;
    my $entry = sprintf "%3i %2i %2i", $2, $4, $5;
    $completeList{$site} = [] unless( exists $completeList{ $site} ) ;
    push @{$completeList{ $site}}, [$edge,$entry];
  }


  my $final_nr = sprintf("%.0f", $screenHash->{'final'}->{'rmax'} / $screenHash->{'final'}->{'dr'} );

#  my %completeList = %{$completeListRef};
  print Dumper( %completeList );
  # First load up core-only screened potential files
  #  These are the core-hole potential including the effective screening
  #  of the other core-level electrons. 
  my %vc_bare;
  my %vpseud1;
  my %vvallel;
  # This allows us to loop over these three potentials and read them into their 
  #  own hash locations using the same code
  my %potTypes = ( 'vc_bare' => \%vc_bare, 'vpseud1' => \%vpseud1, 'vvallel' => \%vvallel );
  if( $screenHash->{ "augment" } ) { %potTypes = ( 'vc_bare' => \%vc_bare ); }
  foreach my $edgeEntry (@{$genHash->{'edgelist'}})
  {
    print "$edgeEntry\n";
    my @edgeentry = split ' ', $edgeEntry;
    my $edgename2 = sprintf("z%03in%02il%02i",$edgeentry[0], $edgeentry[1], $edgeentry[2]);

    foreach my $potType (keys %potTypes )
    {
      my $potfile;
      my @pot;
      my @rad;
      $potfile =  catfile( "zpawinfo", "${potType}${edgename2}" );
      open IN, "<", $potfile or die "Failed to open $potfile for reading\n$!";
      while( my $line = <IN> )
      {
        $line =~ m/(\d\.\d+[Ee][+-]?\d+)\s+(-?\d\.\d+[Ee][+-]?\d+)/ or die "Failed to parse $potfile\n$line";
        push @rad, $1;
        push @pot, $2;
      }
      ${$potTypes{ $potType }}{ "$edgeEntry" } = [ \@rad, \@pot ];
      print ${$potTypes{ $potType }}{ "$edgeEntry" }[0][0] . "\t" .
            ${$potTypes{ $potType }}{ "$edgeEntry" }[1][0] . "\n";
    }
  }

# This framework can walk through every site and then every edge w/i that site
  foreach my $currentSite (keys %completeList)
  {
    print "$currentSite\n";

    # The modeled shell is for a given site and radius
    #  (doesn't depend on edge)
    my @reoptArray;
    foreach my $rad (@{$screenHash->{'shells'}})
    {
      my $reoptName = catfile( $currentSite, sprintf("reopt.R%03.2f",$rad) );
      open IN, "<", $reoptName or die "Failed to open $reoptName\n$!";
      # 1 3
      my @reoptRad; my @reoptPot;
      while( my $line = <IN> )
      {
        $line =~ m/(\d+\.\d+[Ee][+-]\d+)\s+(-?\d+\.\d+[Ee][+-]\d+)\s+(-?\d+\.\d+[Ee][+-]\d+)/
              or die "Failed to parse $reoptName\n\t\t$line";

        push @reoptRad, $1;
        push @reoptPot, $3;
      }
      close IN;
      push @reoptArray, [ \@reoptRad, \@reoptPot ];
    }

    # Now walk through each edge and each radius
    #  1. Read in RPA-screened response
    #  2. Add (interpolated) model
    #  3. Optionally add all-ell - pseduo atomic calc
    #  4. write out
#    for( my $i = 0; $i < scalar( @{ $completeList{ $currentSite }} ); $i++ )

    foreach my $currentEdge (  @{ $completeList{ $currentSite }} )
    {
      my @currentEdge = @{ $currentEdge };
#      my @currentEdge = split ' ', $currentEdge;

#      print "Comp list cur site\t" . ${ $completeList{ $currentSite } }[$i] . "\n";
      print "\t" . $currentEdge[0] . "\t" . $currentEdge[1]  . "\n";
      for( my $r = 0; $r < scalar @{$screenHash->{'shells'}}; $r++ )
      {
        my @vindRad; my @vindPot;

        my $radName = sprintf("zR%03.2f",$screenHash->{'shells'}[$r]) ;
#        print "\t\t$radName\t$screenHash->{'shells'}[$r]\t$r\n";
        printf "\t\t%s\t%.16g\t%i\n", $radName, $screenHash->{'shells'}[$r], $r;

        my $vindName = catfile( $currentSite, $currentEdge[0], $radName, "vind" );
#        print "vind: $vindName\n";
        open IN, "<", $vindName or die "Failed to open $vindName\n$!";
        while( my $line = <IN> )
        {
          $line =~ m/(\d+\.\d+[Ee][+-]\d+)\s+(-?\d+\.\d+[Ee][+-]\d+)/ or die "Failed to parse $vindName\n";
          push @vindRad, $1;
          push @vindPot, $2;
        }
        close IN;

        my $rundir = catfile( $currentSite, $currentEdge[0], $radName );
        print "$rundir\n";
        chdir $rundir or die;

        # If all-electron augmentation then can only make augmented=true version of screened potential
        #  but shouldn't use the faked atomic calculation version
        my $recon_start = 0;
        my $recon_stop = 1;
        if( $screenHash->{ "augment" } )
        {
          $recon_start = -1;
          $recon_stop = -1;
          print "screen_driver used augmentation\n";
        }
        else
        {
          print "screen_driver did not use augmentation\n";
        }

        for( my $reconstruct = $recon_start; $reconstruct <= $recon_stop; $reconstruct++ )
        {
          open OUT, ">", "ipt1" or die "Failed to open ipt\n$!";

          my $len = scalar @{ $reoptArray[$r][0] };
          print OUT "1 2\n$len\n";
          for( my $i = 0; $i < $len; $i++ )
          {
            print OUT "$reoptArray[$r][0][$i]  $reoptArray[$r][1][$i]\n";
          }
          $len = scalar @vindPot;
          print OUT "1 2\n$len\n";
          for( my $i = 0; $i < $len; $i++ )
          {
            print OUT "$vindRad[$i]  $vindPot[$i]\n";
          }

          print "$currentEdge[1]\n";
          $len = scalar @{ $vc_bare{ "$currentEdge[1]" }[0] };
          print $len . "\n";
          print OUT "1 2\n$len\n";
          for( my $i = 0; $i < $len; $i++ )
          {
            print OUT "$vc_bare{ $currentEdge[1] }[0][$i]  $vc_bare{ $currentEdge[1] }[1][$i]\n";
#            my $inv = -1 / $vc_bare{ $currentEdge[1] }[0][$i];
#            print OUT "$vc_bare{ $currentEdge[1] }[0][$i]  $inv\n";
          }

          # True reconstruction of wavefunctions
          if( $reconstruct == -1 )
          {
            print OUT ".false.\n$screenHash->{'final'}->{'dr'} $final_nr\n";
            close OUT;
            system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt") == 0 or die;
          }
          # No reconstruction
          elsif( $reconstruct == 0 )
          {
            print OUT ".false.\n$screenHash->{'final'}->{'dr'} $final_nr\n";
            close OUT;
            system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt_false") == 0 or die;
            move( "rpot", "rpot_false" ) or die "rpot\n$!";
            move( "rpothires", "rpothires_false" ) or die "rpothires\n$!";
            next;
          }
          # Fake reconstruction using atomic all-electron/pseudo difference
          else
          {
            print OUT ".true.\n";
            $len = scalar @{ $vpseud1{ "$currentEdge[1]" }[0] };
            print $len . "\n";
            print OUT "$len\n";
            for( my $i = 0; $i < $len; $i++ )
            {
              print OUT "$vpseud1{ $currentEdge[1] }[0][$i]  $vpseud1{ $currentEdge[1] }[1][$i]\n";
            }
            $len = scalar @{ $vvallel{ "$currentEdge[1]" }[0] };
            print $len . "\n";
            print OUT "$len\n";
            for( my $i = 0; $i < $len; $i++ )
            {
              print OUT "$vvallel{ $currentEdge[1] }[0][$i]  $vvallel{ $currentEdge[1] }[1][$i]\n";
            }
            print OUT "$screenHash->{'final'}->{'dr'} $final_nr\n";
            close OUT;
            system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt") == 0 or die;
          }
#          open IN, "<", "rpot" or die "Failed to open rpot\n$!";
#          while (<IN>) {
#            if( $_ =~ m/^\s*-\d/ ) {
#              print "WARNING bad screening! Repulsive core-hole potentials\n"
#                  . "  BSE stage may fail if this potential is used\n";
#              last;
#            }
#          }

        }

        # back out 3 levels
        $rundir = catfile( updir(), updir(),updir() );
        chdir $rundir;
      }
    }
  }

}

sub runCoreOffset
{
  my ($screenRef, $hashRef) = @_;

  unlink( "core_shift.txt") if( -e "core_shift.txt" );
  return unless( $screenRef->{'core_offset'}->{'enable'} );

#  if( $hashRef->{'general'}->{'program'} ne 'qe' ) {
#    print "WARNING!!!! Possibly units are wrong for Vxc for ABINIT!!!\n\n";
#  }
  my $vxc_factor = 1;
  if( $hashRef->{'general'}->{'program'} eq 'abi' ) {
    $vxc_factor = 2;
  }

  unless( -d "vxc_test" ) {
    mkdir "vxc_test" or die $!;
  }

  my @files = ( "avecsinbohr.ipt", "bvecs", "sitelist", "xyz.wyck" );
  foreach my $f (@files) {
    copy( $f, catfile( "vxc_test", "$f" ) ) or die $!;
  }
  copy( catfile( updir(), "DFT", "potofr" ), catfile( "vxc_test", "rhoofr" ) ) or die $!;
  copy( catfile( updir(), "DFT", "nfft.pot" ), catfile( "vxc_test", "nfft" ) ) or die $!;
  
  chdir "vxc_test" or die $!;

  system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0  or die "Failed to run rhoofg.x\n";
  system("wc -l rhoG2 > rhoofg") == 0 or die "$!\n";
  system("sort -n -k 6 rhoG2 >> rhoofg") == 0 or die "$!\n";
  

  open OUT, ">", "avg.ipt" or die $!;
  print OUT "500 0.01\n";
  close OUT;

  print "$hashRef->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1\n";
  system("$hashRef->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1" );
  if ($? == -1) {
      print "failed to execute: $!\n";
      die;
  }
  elsif ($? & 127) {
      printf "mpi_avg died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
      die;
  }
  else {
    my $errorCode = $? >> 8;
    if( $errorCode != 0 ) {
      die "CALCULATION FAILED\n  mpi_avg exited with value $errorCode\n";
    }
    else {
      printf "mpi_avg exited successfully with value %d\n", $errorCode;
    }
  }

  my @hfin;
  open OUT, ">", "hfinlist" or die $!;
  foreach my $line (@{$hashRef->{'general'}->{'fulllist'}}) {
    print $line ."\tHFIN\n";
    my @temp = split ' ', $line;
    printf OUT "%s %i %i %i %s %i\n", $temp[0], $temp[1], $temp[3], $temp[4], $temp[0], $temp[2];
    push @hfin, [ $temp[3], $temp[4], $temp[0],  $temp[2] ];
  }
  close OUT;

  system("$ENV{'OCEAN_BIN'}/projectVxc.pl") == 0 or die "Failed to run projectVxc.pl\n$!";

  my @newPot;
  my $Vsum = 0;

  open IN, "pot.txt";
  while( my $line = <IN> )
  {
    $line =~ m/^\s*(\S+)/ or die "Failed to parse pot.txt\n$line";
    my $tmpPot = $1 * $vxc_factor;
    push @newPot, $tmpPot;
    $Vsum += $tmpPot;
  }
  close IN;
  chdir updir();

  my @newWsum;
  my @newWshift;

  copy( catfile( "vxc_test", "hfinlist" ), "hfinlist" );
  open OUT, ">", "screen.shells" or die;
  for( my $j = 0; $j < scalar @{$screenRef->{'shells'}}; $j++ ) {
    printf OUT "%.2f\n", $screenRef->{'shells'}[$j];
  }
  close OUT;

  system("$ENV{'OCEAN_BIN'}/projectW.pl") == 0 or die "Failed to run projectW.pl\n$!";
  open IN, "W.txt" or die "Failed to open W.txt\n$!";

  my $Ry2eV = 13.605698066;
  for( my $i = 0; $i < scalar @hfin; $i++ )
  {
    for( my $j = 0; $j < scalar @{$screenRef->{'shells'}}; $j++ )
    {
      my $line = <IN> or die "W.txt was not long enough!";
      $line =~ m/^\s*(\S+)/ or die "Failed to parse W.txt\n$line";
      $newWshift[$i][$j] = $1;
      $newWsum[$j] += $1;
    }
  }
  close IN;

  open OUT, ">", "core_shift.log" or die "Failed to open core_shift.log\n$!";

  # Loop over radii and then hfin
  my $offset;
  my @shiftArray;
  for( my $i = 0; $i < scalar @{$screenRef->{'shells'}}; $i++ )
  {
    my $rad_dir = sprintf("zR%03.2f", $screenRef->{'shells'}[$i] );

    printf OUT "\nRadius = %03.2f Bohr\n", $screenRef->{'shells'}[$i];

    # If we are averaging, new shift by radius
    if( $screenRef->{'core_offset'}->{'average'} )
    {
      $offset = -( $Vsum + $newWsum[$i] ) * $Ry2eV / ( scalar @hfin );
  #    print "$rad_dir\t$offset\n";
      print OUT "  core_offset was set to true. Now set to $offset  \n";
    } else
    {
      $offset = $screenRef->{'core_offset'}->{'energy'};
    }

    print OUT "Site index    New potential   new1/2 Screening   core_offset       total offset\n";
    print OUT "                  (eV)             (eV)              (eV)              (eV)\n";
  # print  "   iiiiiii  -xxxxx.yyyyyyyyy  -xxxxx.yyyyyyyyy  -xxxx.yyyyyyyyy  -xxxx.yyyyyyyyy  -xxxx.yyyyyyyyy  -xxxx.yyyyyyyyy\n";

    # Loop over each atom in hfin
    for( my $j = 0; $j < scalar @hfin; $j++ )
    {
      my $nn = $hfin[$j][0];
      my $ll = $hfin[$j][1];
      my $el = $hfin[$j][2];
      my $el_rank = $hfin[$j][3];

      # Wshift is actually in Ha (convert to Ryd and multiply by 1/2 and nothing happens)
  #    my $shift = ( $Vshift[$j] + $Wshift[$j][$i] ) * $Ry2eV;
      my $shift;
      $shift = ( $newPot[$j] + $newWshift[$j][$i] ) * $Ry2eV;

      $shift += $offset;
      $shift *= -1;
      printf OUT "   %7i   %16.9f  %15.9f  %15.9f  %16.7f\n", $el_rank, $newPot[$j]*$Ry2eV, $newWshift[$j][$i]*$Ry2eV, $offset, $shift;
      $shiftArray[$j][$i] = $shift;

      my $string = sprintf("z%s%04d/n%02dl%02d",$el, $el_rank,$nn,$ll);
      open TMP, ">$string/$rad_dir/cls" or die "Failed to open $string/$rad_dir/cls\n$!";
      print TMP $shift . "\n";
      close TMP;
    }

    print OUT "\n";

  }

  close OUT;

  # Create summaries
  for( my $j = 0; $j < scalar @hfin; $j++ )
  {
    my $nn = $hfin[$j][0];
    my $ll = $hfin[$j][1];
    my $el = $hfin[$j][2];
    my $el_rank = $hfin[$j][3];

    my $dir1 = sprintf "z%2s%04i", $el, $el_rank;
    my $dir2 = sprintf "n%02il%02i", $nn, $ll;
    my $file = catfile( $dir1, $dir2, "cls_rad.txt" );

    print $file . "\n";
    open OUT, ">", $file or die "$!\nFailed to open $file\n";
    print OUT "# rad (Bohr)  CLS (ev)  V_{ind} (Ha)\n";

    for( my $i = 0; $i < scalar @{$screenRef->{'shells'}}; $i++ )
    {
      my $rad_dir = sprintf("zR%03.2f", $screenRef->{'shells'}[$i] );
      open IN, catfile( $dir1, $dir2, $rad_dir, "vind" ) or die "$!\nFailed to open vind: " . catfile( $dir1, $dir2, $rad_dir, "vind" );
      <IN> =~ m/^\s*(\S+)\s+(\S+)/ or die;
      my $vind = $2;

      printf OUT "%03.2f   %16.8f  %16.8f\n", $screenRef->{'shells'}[$i], $shiftArray[$j][$i], $vind;
    }
    close OUT;
  }


#  
#  open OUT, ">", "core_offset" or die $!;
#  if( $screenRef->{'core_offset'}->{'average'} ) {
#    print OUT "true\n";
#  } else {
#    printf OUT "%g\n", $screenRef->{'core_offset'}->{'energy'};
#  }
#  close OUT;
#
#  open OUT, ">", "screen.shells" or die $!;
#  foreach (@{$screenRef->{'shells'}}) {
#    printf OUt "%.2f\n", $_;
#  }
#  close OUT;
#
#  open OUT, ">", "para_prefix" or die $!;
#  print OUT $computerRef->{'para_prefix'} . "\n" ;
#  close OUT;

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
    unless( $complete->{'complete'} )
    {
      print "DIFF:   $t\n" ;
      print Dumper( $newRef->{$t} );
      print Dumper( $oldRef->{$t} );
    }
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
  elsif( ref( $newRef ) eq 'HASH' )
  {
    foreach my $key (keys %$newRef )
    {
      unless( exists $oldRef->{$key} )
      {
        $complete->{'complete'} = JSON::PP::false;
        return;
      }
      recursiveCompare( $newRef->{$key}, $oldRef->{$key}, $complete );
      return unless( $complete->{'complete'} );
    }
  }
  else
  {
#    print "#!  $newRef  $oldRef\n";
    if( looks_like_number( $newRef ) )
    {
      unless( $newRef == $oldRef ) {
        $complete->{'complete'} = JSON::PP::false;
        print $newRef;
        print " ";
        print $oldRef;
        print "  number\n";
      }
#      $complete->{'complete'} = JSON::PP::false unless( $newRef == $oldRef );
    }
    else
    {
#      $complete->{'complete'} = JSON::PP::false unless( $newRef eq $oldRef );
      unless( $newRef eq $oldRef )
      {
        $complete->{'complete'} = JSON::PP::false;
        print $newRef;
        print " ";
        print $oldRef;
        print "  number\n";
      }
    }
  }
}
