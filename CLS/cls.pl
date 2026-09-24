#!/usr/bin/perl
# Copyright (C) 2025 OCEAN collaboration
#   
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#   
# John Vinson, Jan 2025
    
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
  $0 =~ m/(.*)\/cls\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################

my $override = 0;
if( scalar @ARGV > 0 ) {
  $override = $ARGV[0];
  # Any non-zero value will trigger EXX calculation, but EXX will still only
  #  be added to the total if enabled.
  # A value of 2 will cause mpi_avg.x to be run w/o MPI, allowing interactive
}

my @spdf = ( 's', 'p', 'd', 'f' );

my $json = JSON::PP->new;
$json->canonical([1]);
$json->pretty([1]);

###
my $dataFile = catfile( updir(), "Common", "postDefaultsOceanDatafile" );
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

####
# Compatibility with older versions
unless( exists $commonOceanData->{"cls"} ) {
  print "FALLBACK!\n";
  $commonOceanData->{"cls"} = {};
  $commonOceanData->{"cls"}->{'enable'} = $commonOceanData->{"screen"}->{"core_offset"}->{"enable"};
  $commonOceanData->{"cls"}->{'average'} = $commonOceanData->{"screen"}->{"core_offset"}->{"average"};
  $commonOceanData->{"cls"}->{'energy'} = [$commonOceanData->{"screen"}->{"core_offset"}->{"energy"}];
  $override = 1 if( $override == 0 );
  $commonOceanData->{"cls"}->{'do_exx'} = JSON::PP::true;
} 
####
# Early exit if not core-level calculation
my $earlyExit = 0;
$earlyExit = 1 unless( $commonOceanData->{'cls'}->{'enable'} );
$earlyExit = 1 if( $commonOceanData->{'calc'}->{'mode'} eq 'val' );
if( $earlyExit != 0 && $override == 0 )
{
  print "EXIT EARLY\n";
  exit 0;
}


###

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

$dataFile = catfile( updir(), "OPF", "opf.json" );
my $opfData;
if( open( my $in, "<", $dataFile ))
{
  local $/ = undef;
  $opfData = $json->decode(<$in>);
  close($in);
} else {
  die "Failed to open config file $dataFile\n$!";
}

$dataFile = catfile( updir(), "SCREEN", "screen.json" );
my $screenData;
if( open( my $in, "<", $dataFile ))
{ 
  local $/ = undef;
  $screenData = $json->decode(<$in>);
  close($in);
} else { 
  die "Failed to open config file $dataFile\n$!";
} 

# Done loading JSON of previous sections

$dataFile = "cls.json";
my $clsData = {};
if( open( my $in, "<", $dataFile ))
{ 
  local $/ = undef;
  $clsData = $json->decode(<$in>);
  close($in);
} 


runVOffset(  $commonOceanData, $dftData, $clsData, $override);
runWOffset(  $commonOceanData, $screenData, $clsData);
runEXX( $commonOceanData, $dftData, $clsData, $override);

runVWSum( $commonOceanData, $clsData);


exit 0;
################
################

sub runVWSum
{
  my ( $cod, $cls ) = @_;

  my %ZNL;
  foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
    my ($i, $n, $l) = split ' ', $site;
    my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
    my $znl = sprintf "%2s%1i%01i", $z, $n, $l;
    $ZNL{ $znl } = 0;
  }

  
  # In the future, should support setting energies for several edges in the 
  # OCEAN input file, but this is not yet done
  unless( $cod->{'cls'}->{'average'} ) {
    foreach my $znl (keys %ZNL ) {
      $ZNL{ $znl } = $cod->{'cls'}->{'energy'}[0];
    }
  }

  my @shells;
  foreach my $r (@{$cod->{'screen'}->{'shells'}}) {
    my $rad = sprintf "%03.2f", $r;
    push @shells, $rad;
  }

  open OUT, ">", "core_shift.log" or die "Failed to open core_shift.log\n$!";

  #my $Ry2eV = 13.605698066;
  my $Ry2eV = 13.605693122990; #2022 CODATA
  my @ibe = indxByElement( $cod );
  foreach my $rad (@shells) {
    printf OUT "Radius = %.2f Bohr\n", $rad;
    print OUT "Site index    New potential   new1/2 Screening       EXX            core_offset       total offset\n";
    print OUT "                  (eV)             (eV)              (eV)              (eV)               (eV)\n";
    my %avg;
    my %count;
    foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
      my ($i, $n, $l) = split ' ', $site;
      my ($el, $j ) = split ' ', $ibe[$i-1];
      my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
      my $znl = sprintf "%2s%1i%01i", $z, $n, $l;
      my $nl = sprintf "%1i%1s", $n, $spdf[$l];
      
      my $V = $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'};
      my $EXX = 0;
      if( $cod->{'cls'}->{'do_exx'} ) {
        $EXX = $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'};
      }
      my $W = $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'}->{$rad};

      # V is typically negative (though pseudopotentials could be strange), however, for the 
      #  all-electron system, the total potential is going to be attractive for the core-level 
      #  orbitals. The more negative, the more bound the core level, the more energy it will
      #  take to excite it, hence the negative sign out front. 
      # W is positive, the screening reduces the energy needed to excite, so it should have 
      #  the opposite sign of V
      # X is defined to be negative
      $cls->{'total'}->{$el}->{$j}->{$nl}->{$rad} = -( ($V + $W)*$Ry2eV + $EXX );
      
      if( $cod->{'cls'}->{'average'} ) {
        $count{ $znl } += 1;
        $avg{ $znl } += $cls->{'total'}->{$el}->{$j}->{$nl}->{$rad} ;
      }
    }


    foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
      my ($i, $n, $l) = split ' ', $site;
      my ($el, $j ) = split ' ', $ibe[$i-1];
      my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
      my $znl = sprintf "%2s%1i%01i", $z, $n, $l;
      my $nl = sprintf "%1i%1s", $n, $spdf[$l];
      if( $cod->{'cls'}->{'average'} ) {
        $ZNL{ $znl } = ($avg{ $znl })/($count{ $znl });
      }
      my $EXX = 0;
      if( $cod->{'cls'}->{'do_exx'} ) {
        $EXX = $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'};
      }
      $cls->{'total'}->{$el}->{$j}->{$nl}->{$rad} -= $ZNL{ $znl };
      printf OUT  "   %7i   %16.9f  %15.9f  %15.9f  %15.9f %16.7f\n", $j, 
                  $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'}*$Ry2eV, 
                  $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'}->{$rad}*$Ry2eV, $EXX,
                  $ZNL{ $znl }, $cls->{'total'}->{$el}->{$j}->{$nl}->{$rad};
                  
    }
  }
  close OUT;

  $cls->{'total'}->{'units'} = 'eV';

  my $jsonFile = "cls.json";
  open OUT, ">", $jsonFile or die "$!";
  print OUT $json->encode( $cls );
  close OUT;
}

sub runEXX
{
  my ( $cod, $dft, $cls, $or ) = @_;

  unless ( $or ) {
    return unless( $cod->{'cls'}->{'do_exx'} );
  }


  print "EXX\n";

#  if( $cod->{'structure'}->{'metal'} ) {
#  }

  $cls->{'EXX'} = {} unless exists $cls->{'EXX'};
  $cls->{'EXX'}->{'edge'} = {} unless exists $cls->{'EXX'}->{'edge'};
  $cls->{'EXX'}->{'units'} = "eV";

  my @ibe = indxByElement( $cod );
  my @runXtot;
  my %uniqueZ;
  my %uniqueZNL;
  foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
    my ($i, $n, $l) = split ' ', $site;
    my ($el, $j ) = split ' ', $ibe[$i-1];
    my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
    my $nl = sprintf "%1i%1s", $n, $spdf[$l];
    $cls->{'EXX'}->{'edge'}->{$el} = {} unless( exists( $cls->{'EXX'}->{'edge'}->{$el} ) );
    $cls->{'EXX'}->{'edge'}->{$el}->{$j} = {} unless( exists( $cls->{'EXX'}->{'edge'}->{$el}->{$j} ) );
    $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl} = {} unless( exists( $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl} ) );
    my $s = sprintf "%2s %3i %1i %1i %04i", $el, $z, $n, $l, $j;
    unless( $dft->{'bse'}->{'hash'} eq $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl}->{'BSE hash'} ) {
      push @runXtot, $s;
      $uniqueZ{ $z } = 1;
      my $znl = sprintf "%02i %1i %1i", $z, $n, $l;
      $uniqueZNL{ $znl } = 1;
    }
  }

  if( scalar @runXtot > 0 ) {
    unless( -d "exx" ) {
      mkdir "exx" or die "$!";
    }
    chdir "exx";

    writeAvecs( $cod->{'structure'} );
    writeKmesh( $cod->{'bse'} );
    copy( catfile( $ENV{"OCEAN_BIN"}, 'Pquadrature' ), 'Pquadrature') or die "$!";
    copy( catfile( $ENV{"OCEAN_BIN"}, 'sphpts' ), 'sphpts') or die "$!";
    
    foreach my $z (keys %uniqueZ ) {
      my $prjFile = sprintf( "prjfilez%03i", $z );
      copy catfile( updir(), updir(), 'OPF', 'zpawinfo', $prjFile ), $prjFile or die "$!";
    }
    my @gk;
    my $dirname = catdir( updir(), updir(), "OPF", "zpawinfo" );
    opendir (my $dir, $dirname ) or die "$!";
    while ( my $file = readdir( $dir ) ) {
      next unless -f catfile( $dirname, $file );
      if( $file =~ m/^gk/ ) {
#        print catfile( $dirname, $file ) . "\n";
#        $gk{ catfile( $dirname, $file ) } = 1;
        push @gk, $file;
      } elsif( $file =~ m/^fk/ ) {
        push @gk, $file;
      }
    }
    closedir $dir;

    open OUT, ">", "edgelist" or die "Failed to open edgelist\n$!";
    foreach my $znl (keys %uniqueZNL) {
      my @znl = split ' ', $znl;
      printf OUT "%i %i %i\n", $znl[0], $znl[1], $znl[2];
      my $add10 = sprintf "z%03in%02il%02i", $znl[0], $znl[1], $znl[2];
      foreach my $file (@gk) {
        if( $file =~ m/$add10/ ) {
          copy catfile( $dirname, $file), $file;
        }
      }      
    }
    close OUT;
    
    open OUT, ">", "exx.inp" or die "Failed to open sitelist\n$!";
    printf OUT "%i\n", scalar @runXtot;
    foreach my $s (@runXtot) {
      my @s = split ' ', $s;
      my $par = sprintf "parcksv.%2s%04i", $s[0], $s[4];
      copy catfile( updir(), updir(), 'PREP', 'BSE', $par), $par or die "$!";
      print OUT "$s\n";
    }
    close OUT;
    if( $cod->{'structure'}->{'metal'} ) {
      copy catfile( updir(), updir(), 'PREP', 'BSE', "efermiinrydberg.ipt" ), "efermiinrydberg.ipt" or die "$!";
      copy catfile( updir(), updir(), 'PREP', 'BSE', "enkfile" ), "enkfile" or die "$!";
      copy catfile( updir(), updir(), 'PREP', 'BSE', "brange.ipt" ), "brange.ipt" or die "$!";
      copy catfile( updir(), updir(), 'PREP', 'BSE', "nspin" ), "nspin" or die "$!";
#      die "Metals not yet supported for exact exchange\n";
    }

    system("$ENV{'OCEAN_BIN'}/corex.x > corex.log") == 0 or die "Failed to run corex.x\n$!";
    open IN, "<", "corex.log" or die "Failed to open corex.log\n$!";
    while (my $line = <IN> ) {
      if( $line =~ m/AVG\s+(\w+\s+\d+\s+\d+\s+\d+\s+\S+)\s+(-?\d+\.\d+)/ ) {
        my $exx = $2;
        my ($el, $z, $n, $l, $j)  = split ' ', $1;
        my $nl = sprintf "%1i%1s", $n, $spdf[$l];
        $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl}->{'BSE hash'} = $dft->{'bse'}->{'hash'};
        $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'} = $exx;
        print "##### $line\n";
      }
    }
    close IN;
#    <IN>;
#    foreach my $s (@runXtot) {
#      my ($el, $z, $n, $l, $j)  = split ' ', $s;
#      my $nl = sprintf "%1i%1s", $n, $spdf[$l];
#      <IN> =~ m/^(\w+\s+\d+\s+\d+\s+\d+)\s+(-?\d+\.\d+)/ or die "$_";
#      $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl}->{'BSE hash'} = $dft->{'bse'}->{'hash'};
#      $cls->{'EXX'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'} = $2;
#      print $s . "  " . $1 . "\n";
#    } 
#    close IN;
    

    chdir updir();
    my $jsonFile = "cls.json" ;
    open OUT, ">", $jsonFile or die "$!";
    print OUT $json->encode( $cls );
    close OUT; 
  }
}
  
# 
sub runVOffset
{
  my ( $cod, $dft, $cls, $or ) = @_;

  my $t0 = [gettimeofday];
  print "V offset\n";

  $cls->{'V'} = {} unless exists $cls->{'V'};
  $cls->{'V'}->{'site'} = {} unless exists $cls->{'V'}->{'site'};
  $cls->{'V'}->{'edge'} = {} unless exists $cls->{'V'}->{'edge'};
  $cls->{'V'}->{'units'} = "Ryd.";


  my @ibe = indxByElement( $cod );
#  foreach ( indxByElement( $cod ) ) {
#    print "$_ \n";
#  }

  # Do we need to calc total potential?
  # 1. Do the previous calculations exist for these atomic sites?
  # 2. Did the SCF change?

  # List of sites (ignores edge)
  my %uniqueAtomicSite;
  foreach my $edge (@{$cod->{'calc'}->{'edges'}}) {
#    print "$edge \n";
    my $s = (split ' ', $edge )[0];
    $uniqueAtomicSite{ $s } = 1;  # First element of the string
#    $uniqueAtomicSite{ @{ split ' ', $edge }[0] } = 1;  # First element of the string
  }

  my @runVtot;
  foreach my $site (keys %uniqueAtomicSite) { 
    my ( $el, $i ) = split ' ', $ibe[$site-1];
    $cls->{'V'}->{'site'}->{$el} = {} unless exists( $cls->{'V'}->{'site'}->{$el} );
    if( exists $cls->{'V'}->{'site'}->{$el}->{$i} ) {
      unless( $dft->{'scf'}->{'hash'} eq $cls->{'V'}->{'site'}->{$el}->{$i}->{'SCF hash'} ) {
        push @runVtot, $site;
#        print "$site 0 \n";
      }
    } else {
      push @runVtot, $site;
#        print "$site 1 \n";
    }
  }

  if( scalar @runVtot > 0 ) {
    unless( -d "pot" ) {
      mkdir "pot" or die "$!";
    }

    copy( catfile( updir(), "DFT", "potofr" ), catfile( "pot", "rhoofr" ) ) or die $!;
    copy( catfile( updir(), "DFT", "nfft.pot" ), catfile( "pot", "nfft" ) ) or die $!;

    chdir "pot" or die "$!";

    writeAvecs( $cod->{'structure'} );
    writeBvecs( $cod->{'structure'} );
    writeSitelistNew(  $cod, \@runVtot );
    
#    writeSitelist( $screen->{'general'} );
#    writeXYZ( $screen->{'general'} );

    system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0  or die "Failed to run rhoofg.x\n$!";
    system("wc -l rhoG2 > rhoofg") == 0 or die "$!\n";
    system("sort -n -k 6 rhoG2 >> rhoofg") == 0 or die "$!\n";

    open OUT, ">", "avg.ipt" or die $!;
    print OUT "500 0.01\n";
    close OUT;

    if( $or == 2 ) {
      print "$ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1\n";
      system("$ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1" );
    } else {
      print "$cod->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1\n";
      system("$cod->{'computer'}->{'para_prefix'} $ENV{'OCEAN_BIN'}/mpi_avg.x > mpi_avg.log 2>&1" );
    }
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

    foreach my $site (@runVtot) {
      my ( $el, $i ) = split ' ', $ibe[$site-1];
      $cls->{'V'}->{'site'}->{$el}->{$i}->{'SCF hash'} = $dft->{'scf'}->{'hash'};
    }

    my $jsonFile = catfile( updir(), "cls.json" );
    open OUT, ">", $jsonFile or die "$!";
    print OUT $json->encode( $cls );
    close OUT;

    chdir updir();
  }


  # Now that the real-space projection is done, make a complete list
  #  (In the future when multiple edges (e.g. K + L23 ) are supported)
  @runVtot = ();
  foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
    my ($i, $n, $l) = split ' ', $site;
    my ($el, $j ) = split ' ', $ibe[$i-1];
    my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
    my $nl = sprintf "%1i%1s", $n, $spdf[$l];
    $cls->{'V'}->{'edge'}->{$el} = {} unless( exists( $cls->{'V'}->{'edge'}->{$el} ) );
    $cls->{'V'}->{'edge'}->{$el}->{$j} = {} unless( exists( $cls->{'V'}->{'edge'}->{$el}->{$j} ) );
    $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl} = {} unless( exists( $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl} ) );
    my $s = sprintf "%2s %3i %1i %1i %04i", $el, $z, $n, $l, $j;
#    if( exists $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}  ) {
      unless( $dft->{'scf'}->{'hash'} eq $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'SCF hash'} ) {
        push @runVtot, $s;
      }
#    } else {
#      push @runVtot, $s;
#    }
  }

  foreach (@runVtot) {
    print $_ . "\n";
  }


  if( scalar @runVtot > 0 ) {

    my $pot_factor = 1;
    if( $cod->{'general'}->{'program'} eq 'abi' ) {
      $pot_factor = 2;
    }

    my @Vtot = projectVtot( \@runVtot );

    for( my $i = 0; $i < scalar @Vtot; $i ++ ) {
      my @hfin = split ' ', $runVtot[$i];
      my $el = $hfin[0];
      my $j = $hfin[4];
      my $n = $hfin[2];
      my $l = $hfin[3];
      my $nl = sprintf "%1i%1s", $n, $spdf[$l];
      $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'SCF hash'} = $dft->{'scf'}->{'hash'};
      $cls->{'V'}->{'edge'}->{$el}->{$j}->{$nl}->{'pot'} = $Vtot[$i] * $pot_factor;

    }

    my $jsonFile = "cls.json" ;
    open OUT, ">", $jsonFile or die "$!";
    print OUT $json->encode( $cls );
    close OUT; 

  }
}

sub runWOffset
{
  my ( $cod, $screen, $cls ) = @_;

  my $t0 = [gettimeofday];
  print "W offset\n";

  $cls->{'W'} = {} unless exists $cls->{'W'};
  $cls->{'W'}->{'units'} = "Ryd.";

  my @ibe = indxByElement( $cod );

  my @runWtot;
  foreach my $site (@{$cod->{'calc'}->{'edges'}}) {
    my ($i, $n, $l) = split ' ', $site;
    my ($el, $j ) = split ' ', $ibe[$i-1];
    my $z = $cod->{'structure'}->{'znucl'}[$cod->{'structure'}->{'typat'}[$i-1]-1];
    my $nl = sprintf "%1i%1s", $n, $spdf[$l];

    $cls->{'W'}->{$el} = {} unless( exists $cls->{'W'}->{$el});
    $cls->{'W'}->{$el}->{$j} = {} unless( exists $cls->{'W'}->{$el}->{$j});
    $cls->{'W'}->{$el}->{$j}->{$nl} = {} unless( exists $cls->{'W'}->{$el}->{$j}->{$nl});
    my $s = sprintf "%2s %3i %1i %1i %04i", $el, $z, $n, $l, $j;
#    my $s = sprintf "%2s%04i", $el, $j;
    my $clean = 0;

    if( $screen->{'screen'}->{'hash'} eq $cls->{'W'}->{$el}->{$j}->{$nl}->{'hash'} ) {
      if( exists $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'} ) {
        foreach my $r (@{$screen->{'screen'}->{'shells'}}) {
          my $rad = sprintf "03.2f", $r;
          unless( exists $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'}->{$rad} ) {
            $clean = 1;
            last;
          }
        }
      } else {
        $clean = 1;
      }
    } else {
      $clean = 1;
    }
    if( $clean ) {
      $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'} = {};
      push @runWtot, $s;
    }
  }

  return if( scalar @runWtot == 0 );
  my @rad;
  foreach my $r (@{$screen->{'screen'}->{'shells'}}) {
    push @rad,sprintf( "%03.2f", $r);
  }

  foreach my $s (@runWtot) {
    my @edge = split ' ', $s;
    my $el = $edge[0];
    my $j = $edge[4];
    my $nl = sprintf "%1i%1s", $edge[2], $spdf[$edge[3]];
    my @W = projectW( $s, \@rad );
    for( my $i = 0; $i < scalar @rad; $i++) {
      $cls->{'W'}->{$el}->{$j}->{$nl}->{'pot'}->{$rad[$i]} = $W[$i]
    }
    $cls->{'W'}->{$el}->{$j}->{$nl}->{'hash'} = $screen->{'screen'}->{'hash'};
  }
  
  my $jsonFile = "cls.json" ;
  open OUT, ">", $jsonFile or die "$!";
  print OUT $json->encode( $cls );
  close OUT; 

}


sub writeAvecs
{
  my $structureRef = $_[0];

  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ ) {
    printf  OUT "%.16g  %.16g  %.16g\n", $structureRef->{'avecs'}[$i][0],
                                $structureRef->{'avecs'}[$i][1],
                                $structureRef->{'avecs'}[$i][2];

  }
  close OUT;
}

sub writeKmesh
{
  my $ref = $_[0];

  open OUT, ">", "kmesh.ipt" or die "Failed to open kmesh.ipt\n$!";
  printf OUT "%i %i %i\n", $ref->{'kmesh'}[0], $ref->{'kmesh'}[1], $ref->{'kmesh'}[2];
  close OUT;
}

sub writeBvecs
{
  my $structureRef = $_[0];
  open OUT, ">", "bvecs" or die "Failed to open bvecs\n$!";
  for( my $i = 0; $i < 3; $i++ ) {
    printf  OUT "%.16g  %.16g  %.16g\n", $structureRef->{'bvecs'}[$i][0],
                                $structureRef->{'bvecs'}[$i][1],
                                $structureRef->{'bvecs'}[$i][2];

  }
  close OUT;
}

sub writeSitelistNew #(  $cod->{'structure'}, \@runVsite );
{
  my ( $genRef, @sitelist ) = ( $_[0], @{$_[1]} );

  my @siteCounter;
  open OUT, ">", "sitelist.new" or die "$!";
  print OUT scalar @sitelist . "\n";
  my @output;
  my $natom = scalar @{$genRef->{'structure'}->{'typat'}};
  for ( my $i = 0; $i < $natom; $i ++ ) {
    my $t = $genRef->{'structure'}->{'typat'}[$i];
    my $z = $genRef->{'structure'}->{'znucl'}[$t-1];
    $siteCounter[$z] ++ ;
    push @output, sprintf "%2s %8i %.16g %.16g %.16g\n", $genRef->{'structure'}->{'elname'}[$t-1], 
                $siteCounter[$z], $genRef->{'structure'}->{'xred'}[$i][0], 
                $genRef->{'structure'}->{'xred'}[$i][1], $genRef->{'structure'}->{'xred'}[$i][2];
  }

  foreach my $s (@sitelist ) {
    print OUT $output[$s-1];
  }
  close OUT;
}

sub indxByElement
{
  my ( $genRef ) = @_;

  my @indxByEl;
  my @siteCounter;
  my $natom = scalar @{$genRef->{'structure'}->{'typat'}};
  for ( my $i = 0; $i < $natom; $i ++ ) {
    # Need to de-reference back to z to account for how different spins/U/etc are done in QE
    my $t = $genRef->{'structure'}->{'typat'}[$i];
    my $z = $genRef->{'structure'}->{'znucl'}[$t-1];
    $siteCounter[$z] ++ ;
    push @indxByEl, sprintf "%2s %04i", $genRef->{'structure'}->{'elname'}[$t-1],
                $siteCounter[$z];
  }

  return @indxByEl;
}

sub writeXYZ 
{
  my $genRef = $_;
  open OUT, ">", "xyz.wyck" or die $!;
  print OUT scalar @{$genRef->{'wyck'}};
  print OUT "\n";
  foreach (@{$genRef->{'wyck'}}) {
    print OUT $_ . "\n";
  }
  close OUT;
}

sub writeSitelist
{
  my $genRef = $_;
  open OUT, ">", "sitelist" or die $!;
  print OUT scalar @{$genRef->{'sitelist'}};
  print OUT "\n";
  foreach (@{$genRef->{'sitelist'}}) {
    print OUT $_ . "\n";
  }
  close OUT;
}

sub projectVtot
{
  my @edges = @{$_[0]};

  my @V;
  foreach my $edge ( @edges ) {
#    print "$edge \n";

    my @hfin = split ' ', $edge;
    my $cf = sprintf "coreorbz%03in%02il%02i", $hfin[1], $hfin[2], $hfin[3];
    my $coreFile = catfile( updir(), "OPF", "zpawinfo", $cf );
    open IN, "<", $coreFile or die "Failed to open $coreFile\n$!";
    <IN>;
    my @rad;
    my @wvfn;
    while ( my $line = <IN> )
    {
      $line =~ m/(\S+)\s+(\S+)/;
      push @rad, $1;
      push @wvfn, $2;
    }
    close IN;

    my $rat = $rad[-1]/$rad[0];
    my $dl = log( $rat )/ ($#rad);
    my $xrat = exp( $dl );

    my $xr1 = sqrt( $xrat ) - sqrt( 1.0/$xrat);
    my $rmin = $rad[0] / $xrat;


    my $sum = 0;
    my $sum2 = 0;
    for( my $i = 0; $i <= $#rad; $i++ )
    {
      my $temp_rad = $rmin * $xrat**($i+1);
      $sum += $rad[$i] * $rad[$i] * $rad[$i] * $xr1 * $wvfn[$i]**2;
      $sum2 += $rad[$i] * $xr1 * $wvfn[$i]**2;
    }
#    print "$sum2\n";

    my $filename = catfile( "pot", sprintf("avg%2s%04i", $hfin[0], $hfin[4]) );
    open IN, $filename or die "Failed to open $filename\n$!";
    my @prad;
    my @pot;
    while( my $line = <IN> )
    {
      $line =~ m/\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+/;
      push @prad, $1;
      push @pot, $2;
    }
    close IN;

    $sum = 0;
    my $j = 0;
    for( my $i = 0; $i <= $#rad; $i++ )
    {
      while( $rad[$i] > $prad[$j] )
      {
        die if( $j > $#prad );
        $j++;
      }
      my $interp = @pot[$j];
      if( $j > 0 && $j < $#prad )
      {
        my $run = $prad[$j+1]-$prad[$j];
        my $slope = ($pot[$j+1]-$pot[$j]) / $run;
        $interp += $slope * ( $rad[$i] - $prad[$j] );
      }
      $sum += $rad[$i] * $xr1 * $wvfn[$i]**2 * $interp;
    }

    push @V, $sum;
  }

  return @V;
}

sub projectW
{
  my $edge = $_[0];
  my @shells = @{$_[1]};

  my @Wpot;
  
  my @hfin = split ' ', $edge;
  my $el = $hfin[0];
  my $z = $hfin[1];
  my $n = $hfin[2];
  my $l = $hfin[3];
  my $indx = $hfin[4];
  
  my $cf = sprintf "coreorbz%03in%02il%02i", $z, $n, $l;
  my $coreFile = catfile( updir(), "OPF", "zpawinfo", $cf );
  open IN, "<", $coreFile or die "Failed to open $coreFile\n$!";
  <IN>;
  my @rad;
  my @wvfn;
  while ( my $line = <IN> )
  {
    $line =~ m/(\S+)\s+(\S+)/;
    push @rad, $1;
    push @wvfn, $2;
  }
  close IN;

  my $rat = $rad[-1]/$rad[0];
  my $dl = log( $rat )/ ($#rad);
  my $xrat = exp( $dl );
  my $xr1 = sqrt( $xrat ) - sqrt( 1.0/$xrat);
  my $rmin = $rad[0] / $xrat;

  my $zel = sprintf "z%s%04d", $el, $indx;
  my $znl = sprintf "n%02il%02i", $n, $l;


  foreach my $radius (@shells) {
    my $rad_dir = sprintf("zR%03.2f", $radius );
    my $filename = catfile( updir(), 'SCREEN', $zel, $znl, $rad_dir, 'ropt' );
    open IN, $filename or die "Failed to open $filename\n$!";


    my @prad;
    my @pot;
    while( my $line = <IN> )
    {
      $line =~ m/(\S+)\s+\S+\s+\S+\s+(\S+)/;
      push @prad, $1;
      push @pot, $2;
    }
    close IN;

    my $sum = 0;
    my $j = 0;
    for( my $i = 0; $i <= $#rad; $i++ )
    {
      while( $rad[$i] > $prad[$j] )
      {
        die if( $j > $#prad );
        $j++;
      }
      my $interp = @pot[$j];
      if( $j > 1 && $j < $#prad - 1 )
      {
        $interp = $pot[$j-2] * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j+1])
                            / (($prad[$j-2]-$prad[$j-1])*($prad[$j-2]-$prad[$j])*($prad[$j-2]-$prad[$j+1]) )
                + $pot[$j-1] * ($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j+1])*($rad[$i]-$prad[$j-2])
                            / (($prad[$j-1]-$prad[$j])*($prad[$j-1]-$prad[$j+1])*($prad[$j-1]-$prad[$j-2]) )
                + $pot[$j]   * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j+1])*($rad[$i]-$prad[$j-2])
                            / (($prad[$j]-$prad[$j-1])*($prad[$j]-$prad[$j+1])*($prad[$j]-$prad[$j-2]) )
                + $pot[$j+1] * ($rad[$i]-$prad[$j-1])*($rad[$i]-$prad[$j])*($rad[$i]-$prad[$j-2])
                            / (($prad[$j+1]-$prad[$j-1])*($prad[$j+1]-$prad[$j])*($prad[$j+1]-$prad[$j-2]) );
      }
      elsif( $j == $#prad - 1 || $j == 1 || $j == 0 )
      {
        my $run = $prad[$j+1]-$prad[$j];
        my $slope = ($pot[$j+1]-$pot[$j]) / $run;
        $interp += $slope * ( $rad[$i] - $prad[$j] );
      }
      $sum += $rad[$i] * $xr1 * $wvfn[$i]**2 * $interp;
    }
    push @Wpot, $sum;
  }

  return @Wpot;
}
  

