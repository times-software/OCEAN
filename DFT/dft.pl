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
use Digest::MD5 qw(md5_hex);

use FindBin;
use lib $FindBin::Bin;
require 'QEdriver.pl';
require 'ABIdriver.pl';

use Time::HiRes qw( gettimeofday tv_interval );

print localtime() .  "\n";
my ( $startSeconds, $startMicroseconds) = gettimeofday;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/dft\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_ESPRESSO_PW"} ) {$ENV{"OCEAN_ESPRESSO_PW"} = $ENV{"OCEAN_BIN"} . "/pw.x"; }
if (! $ENV{"OCEAN_ESPRESSO_PP"} ) {$ENV{"OCEAN_ESPRESSO_PP"} = $ENV{"OCEAN_BIN"} . "/pp.x"; }
if (! $ENV{"OCEAN_ESPRESSO_PH"} ) {$ENV{"OCEAN_ESPRESSO_PH"} = $ENV{"OCEAN_BIN"} . "/ph.x"; }
if (! $ENV{"OCEAN_ABINIT"} ) {$ENV{"OCEAN_ABINIT"} = $ENV{"OCEAN_BIN"} . "/abinit"; }
if (! $ENV{"OCEAN_CUT3D"} ) {$ENV{"OCEAN_CUT3D"} = $ENV{"OCEAN_BIN"} . "/cut3d"; }

#my $driver = catdir( $ENV{"OCEAN_BIN"}, "QEdriver.pl" );
#require "$driver";
my @timeSections = ( 'scf', 'density', 'potential', 'bse', 'screen' );

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
my $fake->{ 'complete' } = JSON::PP::false;
foreach my $sec (@timeSections) {
  $newDftData->{$sec}->{'time'} = $dftData->{$sec}->{'time'} if( exists $dftData->{$sec}->{'time'} );
}

# First we check, using the SCF flag to store result
# 1) Was previous run?
# 2) Structure matches
# 3) PSP matches
# 
# Any failures skip future tests and SCF to not done (which in turn, cascades to all futher runs)
$newDftData->{'scf'}->{'complete'} = JSON::PP::true;

# (was previous run done)
$newDftData->{'scf'}->{'complete'} = JSON::PP::false 
    unless( exists $dftData->{'scf'} && exists $dftData->{'scf'}->{'complete'} && $dftData->{'scf'}->{'complete'});

# (build the structure and check )
$newDftData->{'scf'}->{'complete'} = JSON::PP::false unless( exists $dftData->{'structure'} );

$newDftData->{'structure'} = {};

my @structureList = ( "typat", "xred", "znucl", "avecs", "zsymb", "valence_electrons", "bvecs", 
                      "metal", "charge", "magnetization", "semicore_electrons" );
copyAndCompare( $newDftData->{'structure'}, $commonOceanData->{'structure'}, $dftData->{'structure'}, 
                $newDftData->{'scf'}, \@structureList );

# Additional items for later stages -- unlikely these changed w/o changing manditory ones, but 
@structureList = ( "elname" );
copyAndCompare( $newDftData->{'structure'}, $commonOceanData->{'structure'}, $dftData->{'structure'},
                $fake, \@structureList );

$newDftData->{'psp'} = {};
my @pspList = ( "pphash" );
copyAndCompare( $newDftData->{'psp'}, $commonOceanData->{'psp'}, $dftData->{'psp'},
                $newDftData->{'scf'}, \@pspList );

@pspList = ( "pp_list", "ppdir" );
copyAndCompare( $newDftData->{'psp'}, $commonOceanData->{'psp'}, $dftData->{'psp'},
                $fake, \@pspList );

# Now do the general DFT parts

# Only check the first list against previous runs
my @generalList = ( "degauss", "ecut", "fband", "functional", "noncolin", "nspin", "occopt", 
                    "program", "smag", "spinorb", "verbatim", "isolated" );
my @generalSecondaryList = ( "calc_force", "calc_stress", "diagonalization", "mixing", "mixing_mode",
                             "nstep", "redirect", "startingwfc", "startingpot", "tmp_dir", "abpad" );
$newDftData->{'general'} = {};
copyAndCompare( $newDftData->{'general'}, $commonOceanData->{'dft'}, $dftData->{'general'},
                $newDftData->{'scf'}, \@generalList );
copyAndCompare( $newDftData->{'general'}, $commonOceanData->{'dft'}, $dftData->{'general'},
                $fake, \@generalSecondaryList );

# EXX if functional is specified 
if( $newDftData->{'general'}->{'functional'} ne 'default' )
{
  $newDftData->{'general'}->{'exx'} = {};
  copyAndCompare( $newDftData->{'general'}->{'exx'}, $commonOceanData->{'dft'}->{'exx'}, 
                  $dftData->{'general'}->{'exx'},
                  $newDftData->{'scf'}, [ 'qmesh' ] );
}

# LDA+U
$newDftData->{'general'}->{'ldau'} = {};
copyAndCompare( $newDftData->{'general'}->{'ldau'}, $commonOceanData->{'dft'}->{'ldau'}, 
                $dftData->{'general'}->{'ldau'},
                $newDftData->{'scf'}, [ 'enable' ] );
if( $newDftData->{'general'}->{'ldau'}->{'enable'} )
{
  my @ldauList = ( "Hubbard_J", "Hubbard_J0", "Hubbard_U", "Hubbard_V", "U_projection_type", 
                   "lda_plus_u_kind", "scf_only" );
  copyAndCompare( $newDftData->{'general'}->{'ldau'}, $commonOceanData->{'dft'}->{'ldau'},
                $dftData->{'general'}->{'ldau'}, $newDftData->{'scf'}, \@ldauList );
}


# and finally density run information
my @scfList = ( "auto", "kmesh", "kshift", "toldfe" );
my @scfSecondaryList = ( "poolsize" );

copyAndCompare( $newDftData->{'scf'}, $commonOceanData->{'dft'}->{'den'}, $dftData->{'scf'}, 
                $newDftData->{'scf'}, \@scfList );
copyAndCompare( $newDftData->{'scf'}, $commonOceanData->{'dft'}->{'den'}, $dftData->{'scf'},
                $fake, \@scfSecondaryList );

checkSetGamma( $newDftData->{'scf'} );
# Computer information
$newDftData->{'computer'} = {};
my @computerList = ( "cpu_factors", "cpu_square_factors", "ncpus", "para_prefix", "ser_prefix" );
copyAndCompare( $newDftData->{'computer'}, $commonOceanData->{'computer'}, $dftData->{'computer'},
                $fake, \@computerList );

# At this point SCF is sorted out
# all subsequent stages can be complete if SCF isn't being re-run
$newDftData->{'density'}->{'complete'} = JSON::PP::false;
$newDftData->{'potential'}->{'complete'} = JSON::PP::false;
$newDftData->{'epsilon'}->{'complete'} = JSON::PP::false;
#$newDftData->{'screen'}->{'complete'} = JSON::PP::false;
#$newDftData->{'bse'}->{'complete'} = JSON::PP::false;


if( $newDftData->{'scf'}->{'complete'} ) {
  print "Re-using previous SCF run\n";
  
  $newDftData->{'density'}->{'complete'} = $dftData->{'density'}->{'complete'} 
      if( exists $dftData->{'density'}->{'complete'} );
  $newDftData->{'potential'}->{'complete'} = $dftData->{'potential'}->{'complete'}
      if( exists $dftData->{'potential'}->{'complete'} );
  $newDftData->{'epsilon'}->{'complete'} = $dftData->{'epsilon'}->{'complete'}
      if( exists $dftData->{'epsilon'}->{'complete'} );
#  $newDftData->{'screen'}->{'complete'} = $dftData->{'screen'}->{'complete'}
#      if( exists $dftData->{'screen'}->{'complete'} );
#  $newDftData->{'bse'}->{'complete'} = $dftData->{'bse'}->{'complete'}
#      if( exists $dftData->{'bse'}->{'complete'} );

  # If SCF already run, copy additional info from previous time
  my @scfCopyList = ( "npool", "ncpus", "fermi", "etot", "time", "version", "nelec", "lowest", "highest", "hash" );
  copyAndCompare( $newDftData->{'scf'}, $dftData->{'scf'}, $dftData->{'scf'}, $fake, \@scfCopyList );

#  my @bseCopyList = ( "")
  
#  copyAndCompare( $newDftData->{'bse'}, $dftData->{'bse'}, $dftData->{'bse'}, $fake, [ "completed" ] );

  #Copy record of all completed NSCF runs
  copyAndCompare( $newDftData, $dftData, $dftData, $fake, [ 'znscf' ] );

} else {
  print "Need SCF run\n";
  $newDftData->{'znscf'} = {};
}


### Determining epsilon w/ DFPT
my @epsList = ( "metal_max", "metal_min", "method", "min_gap", "thresh" );
copyAndCompare( $newDftData->{'epsilon'}, $commonOceanData->{'dft'}->{'epsilon'}, $dftData->{'epsilon'},
                $newDftData->{'epsilon'}, \@epsList );
if( $newDftData->{'epsilon'}->{'method'} eq "input" ) {
  $newDftData->{'epsilon'}->{'complete'} = JSON::PP::true; # if( $newDftData->{'epsilon'}->{'method'} eq "input" );
  $newDftData->{'structure'}->{'epsilon'} = $commonOceanData->{'structure'}->{'epsilon'};
} elsif( $newDftData->{'epsilon'}->{'complete'} ) {
  $newDftData->{'epsilon'}->{'epsilon'} = $dftData->{'epsilon'}->{'epsilon'};
}

my @nscf_InitialList = ();

# Step 0 -- what info needs to be available for PREP?
$newDftData->{'bse'} = {} unless exists $newDftData->{'bse'};
copyAndCompare( $newDftData->{'bse'}, $commonOceanData->{'calc'}, $dftData->{'bse'},
                $fake, [ 'photon_q', 'nonzero_q' ] );
copyAndCompare( $newDftData->{'bse'}, $commonOceanData->{'dft'}->{'bse'}, $dftData->{'bse'},
                $fake, [ 'split' ] );
copyAndCompare( $newDftData->{'bse'}, $commonOceanData->{'bse'}, $dftData->{'bse'},
                $fake, [ 'con_start', 'val_stop' ] );

#$newDftData->{'bse'}->{'nonzero_q'} = JSON::PP::true;
#$newDftData->{'bse'}->{'split'} = JSON::PP::true;
#$newDftData->{'bse'}->{'photon_q'} = [ 0.01, 0.01, 0.01 ];
$newDftData->{'bse'}->{'directories'} = [];

# Step 1 -- do we have split=false && photon_q ?
my $niter = 1;
my $nosplit = 0;
if( $newDftData->{'bse'}->{'nonzero_q'} ) { 
  if( $newDftData->{'bse'}->{'split'} ) {
    $niter = 2;
  } else {
    $nosplit = 1;
  }
}

#if( $nosplit ) {die "No split not implemented yet\n";}



# Step 1 -- load up general info
for( my $i = 0; $i < $niter; $i ++ ) {

  push @nscf_InitialList, {};
  my @bseList = ( "toldfe", "poolsize", "diagonalization" );
  copyNoCompare( $nscf_InitialList[$i], $commonOceanData->{'dft'}->{'bse'}, \@bseList );
  @bseList = ( "kmesh", "kshift", "con_start", "val_stop" );
  copyNoCompare( $nscf_InitialList[$i], $commonOceanData->{'bse'}, \@bseList );

  # First loop is conduction bands, second is valence only on -q grid
  if( $i == 0 ) {
    @bseList = ( "nbands" );
    copyNoCompare( $nscf_InitialList[$i], $commonOceanData->{'bse'}, \@bseList );
  } else {  
    if( $nscf_InitialList[$i]->{'val_stop'} >= 1 ) {
#      @bseList = ( "nbands" );
      $nscf_InitialList[$i]->{'nbands'} = $nscf_InitialList[$i]->{'val_stop'};
    } else {
#      @bseList = ( "fband" );
      copyNoCompare( $nscf_InitialList[$i], $commonOceanData->{'dft'}, ["fband"] );
    }
    printf "%.6f  %.6f  %.6f", $nscf_InitialList[$i]->{'kshift'}[0], 
                              $nscf_InitialList[$i]->{'kshift'}[1], $nscf_InitialList[$i]->{'kshift'}[2];
    shiftKpointsByPhoton( $nscf_InitialList[$i], $newDftData->{'bse'}->{'photon_q'} );
    printf "  %.6f  %.6f  %.6f\n", $nscf_InitialList[$i]->{'kshift'}[0], 
                              $nscf_InitialList[$i]->{'kshift'}[1], $nscf_InitialList[$i]->{'kshift'}[2];
  }
  my $dirname;
  if( $nosplit ) {
    $nscf_InitialList[$i]->{'nonzero_q'} = JSON::PP::true;
    $nscf_InitialList[$i]->{'photon_q'} = $newDftData->{'bse'}->{'photon_q'};
    $dirname = sprintf "ks%i_%i_%iq%.6f_%.6f_%.6f", $nscf_InitialList[$i]->{'kmesh'}[0],
                    $nscf_InitialList[$i]->{'kmesh'}[1], $nscf_InitialList[$i]->{'kmesh'}[2],
                    $nscf_InitialList[$i]->{'kshift'}[0], $nscf_InitialList[$i]->{'kshift'}[1],
                    $nscf_InitialList[$i]->{'kshift'}[2];
  } else {
    $nscf_InitialList[$i]->{'nonzero_q'} = JSON::PP::false;
    $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $nscf_InitialList[$i]->{'kmesh'}[0],
                    $nscf_InitialList[$i]->{'kmesh'}[1], $nscf_InitialList[$i]->{'kmesh'}[2],
                    $nscf_InitialList[$i]->{'kshift'}[0], $nscf_InitialList[$i]->{'kshift'}[1],
                    $nscf_InitialList[$i]->{'kshift'}[2];
  }
  push @{$newDftData->{'bse'}->{'directories'}}, $dirname;
  checkSetGamma( $nscf_InitialList[$i] );
}
$newDftData->{'bse'}->{'brange'} = [ 0, 0, 0, $commonOceanData->{'bse'}->{'nbands'} ];

# repeat but with screen info, $niter is location of last item in InitalList
unless( $commonOceanData->{'calc'}->{'mode'} eq 'val' && ! ( $commonOceanData->{'screen'}->{'mode'} eq 'grid' ) )
{
  push @nscf_InitialList, {};
  my $i = $niter;
  my @bseList = ( "toldfe", "poolsize", "diagonalization" );
  copyNoCompare( $nscf_InitialList[$i], $commonOceanData->{'dft'}->{'screen'}, \@bseList );
  @bseList = ( "kmesh", "kshift", "nbands" );
  copyNoCompare( $nscf_InitialList[$i], $commonOceanData->{'screen'}, \@bseList );
  my $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $nscf_InitialList[$i]->{'kmesh'}[0],
                    $nscf_InitialList[$i]->{'kmesh'}[1], $nscf_InitialList[$i]->{'kmesh'}[2],
                    $nscf_InitialList[$i]->{'kshift'}[0], $nscf_InitialList[$i]->{'kshift'}[1],
                    $nscf_InitialList[$i]->{'kshift'}[2];
  $newDftData->{'screen'}->{'directories'} = [ $dirname ];
  $newDftData->{'screen'}->{'enable'} = JSON::PP::true;
  $newDftData->{'screen'}->{'brange'} = [ 0, 0, 0, $commonOceanData->{'screen'}->{'nbands'} ];
} else {
  $newDftData->{'screen'}->{'enable'} = JSON::PP::false;
}

my $nscf_TodoList = {};
foreach my $hashRef (@nscf_InitialList) {
  my $dirname = sprintf "%i_%i_%iq%.6f_%.6f_%.6f", $hashRef->{'kmesh'}[0],
                    $hashRef->{'kmesh'}[1], $hashRef->{'kmesh'}[2],
                    $hashRef->{'kshift'}[0], $hashRef->{'kshift'}[1],
                    $hashRef->{'kshift'}[2];
  if( $hashRef->{'nonzero_q'} ) {
    $dirname = 'ks' . $dirname;
  } else {
    $dirname = 'k' . $dirname;
  }

  my $addThisCalculation = 1;
  if( defined $newDftData->{'znscf'}->{$dirname} ) {
    print "$dirname exists\n";
    $addThisCalculation = 0;
    if( $hashRef->{'toldfe'} < $newDftData->{'znscf'}->{ $dirname }->{'toldfe'} ) {
      $addThisCalculation = 1;
    }
    if( defined( $hashRef->{'nbands'} ) ) {
      if( defined( $newDftData->{'znscf'}->{ $dirname }->{'nbands' } ) ) {
        if( $hashRef->{'nbands'} > $newDftData->{'znscf'}->{ $dirname }->{'nbands' } ) {
          $addThisCalculation = 1;
        }
      } else {  # new calc has nbands while old didn't
        $addThisCalculation = 1;
      }
    } 
    else {  # else, using fband
      unless( defined( $newDftData->{'znscf'}->{ $dirname }->{'fband'} ) &&
              $newDftData->{'znscf'}->{ $dirname }->{'fband'} >= $hashRef->{'fband'} ) {
      $addThisCalculation = 1;
      }
    }
#    if( $newDftData->{'znscf'}->{ $dirname }->{'toldfe'} <= $hashRef->{'toldfe'} &&
#        $newDftData->{'znscf'}->{ $dirname }->{'nbands' } >= $hashRef->{'nbands'} ) {
#      $addThisCalculation = 0;
#    }
  } else {
    print "$dirname is new\n";
  }
  if( $addThisCalculation ) {
    print "Will run $dirname\n";
    if( defined( $nscf_TodoList->{ $dirname } ) ) {
      print "Condensing two runs\n";
      $hashRef->{'toldfe'} = $nscf_TodoList->{ $dirname }->{'toldfe'} 
          if( $nscf_TodoList->{ $dirname }->{'toldfe'} < $hashRef->{'toldfe'} );
      $hashRef->{'nbands'} = $nscf_TodoList->{ $dirname }->{'nbands'} 
          if( $nscf_TodoList->{ $dirname }->{'nbands'} > $hashRef->{'nbands'} );
    }
    $nscf_TodoList->{ $dirname } = dclone( $hashRef );
  }
  else {   # check energy parsing
    if( ( $newDftData->{'bse'}->{'con_start'} != $newDftData->{'znscf'}->{ $dirname }->{'con_start'} 
        && ( defined( $newDftData->{'bse'}->{'con_start'}) || 
             defined( $newDftData->{'znscf'}->{ $dirname }->{'con_start'} ) ) ) 
   || ($newDftData->{'bse'}->{'val_stop'} != $newDftData->{'znscf'}->{ $dirname }->{'val_stop'} 
        && ( defined( $newDftData->{'bse'}->{'val_stop'}) || 
             defined($newDftData->{'znscf'}->{ $dirname }->{'val_stop'} ) ) ) )  {
      my $errorCode;
      if( $newDftData->{'general'}->{'program'} eq "qe" ) {
        my $temp = $newDftData->{'bse'}->{'split'};
        $newDftData->{'bse'}->{'split'} = JSON::PP::false;
        $newDftData->{'znscf'}->{ $dirname }->{'con_start'} = $hashRef->{'con_start'};
        $newDftData->{'znscf'}->{ $dirname }->{'val_stop'} = $hashRef->{'val_stop'};
        $errorCode = QErunParseEnergies( $newDftData, $newDftData->{'znscf'}->{ $dirname }, 0 );
        $newDftData->{'bse'}->{'split'} = $temp;
#        $newDftData->{'znscf'}->{ $dirname }->{'con_start'} = $newDftData->{'bse'}->{'con_start'};
#        $newDftData->{'znscf'}->{ $dirname }->{'val_stop'} = $newDftData->{'bse'}->{'val_stop'};
      }
      exit $errorCode if( $errorCode );
    }
  }
}

my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);
open OUT, ">", "derp.json" or die;
print OUT $json->encode($nscf_TodoList);
close OUT;


my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);
open OUT, ">", "dft2.json" or die;
print OUT $json->encode($newDftData);
close OUT;

#exit 0;

### BSE


#if( $newDftData->{'bse'}->{'complete'} )
#{   
#  my @bseCopyList = ( "npool", "ncpus", "fermi", "etot", "time", "version", "nelec", "lowest", "highest", "hash", "brange" );
#  copyAndCompare( $newDftData->{'bse'}, $dftData->{'bse'}, $dftData->{'bse'}, $fake, \@bseCopyList );
#}


### SCREEN

#if( $newDftData->{'screen'}->{'complete'} )
#{
#  my @screenCopyList = ( "npool", "ncpus", "fermi", "etot", "time", "version", "nelec", "lowest", "highest", "hash", "brange" );
#  copyAndCompare( $newDftData->{'screen'}, $dftData->{'screen'}, $dftData->{'screen'}, $fake, \@screenCopyList );
#}
#
#$newDftData->{'screen'}->{'enable'} = JSON::PP::true;
#if( $commonOceanData->{'calc'}->{'mode'} eq 'val' )
#{
#  $newDftData->{'screen'}->{'enable'} = JSON::PP::false unless( $commonOceanData->{'screen'}->{'mode'} eq 'grid' );
#}



print "Done parsing input for DFT stage\n";

# Save current outlook
my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);
open OUT, ">", "dft.json" or die;
print OUT $json->encode($newDftData);
close OUT;

### Need to write abstraction to support multiple DFT codes
unless( $newDftData->{'general'}->{'program'} eq "qe" ||
        $newDftData->{'general'}->{'program'} eq "abi" ) {
  print "Only QE and ABINIT supported at the moment!\t\t" . $newDftData->{'general'}->{'program'} . "\n";
  exit 1;
}

#call density stage
unless( $newDftData->{'scf'}->{'complete'} )
{
  my $errorCode = 0;
  my $t0 = [gettimeofday];
  if( $newDftData->{'general'}->{'program'} eq "qe" ) {
    $errorCode = QErunDensity( $newDftData );
    print "$errorCode\n";
  } elsif ( $newDftData->{'general'}->{'program'} eq "abi" ) {
    $errorCode = ABIrunDensity( $newDftData );
    print "$errorCode\n";
  } else {
    $errorCode = 1;
  }
  exit $errorCode if( $errorCode != 0 ) ;

  $newDftData->{'scf'}->{'complete'} = JSON::PP::true;

  my $s = $json->encode($newDftData->{'psp'}->{'pphash'});
  foreach ( 'general', 'scf', 'structure' )
  {
    $s .= $json->encode($newDftData->{$_});
  }
#  print "$s\n\n\n";
  $newDftData->{'scf'}->{'hash'} = md5_hex( $s );

  $newDftData->{'scf'}->{'time'} = tv_interval( $t0 );
  
  #TODO: clean old nscf runs here
  $newDftData->{'znscf'} = {};

  open OUT, ">", "dft.json" or die;
  print OUT $json->encode($newDftData);
  close OUT;
  print "SCF stage complete, total energy: $newDftData->{'scf'}->{'etot'}\n";

} else {
  $newDftData->{'scf'}->{'time'} = $dftData->{'scf'}->{'time'};
}

# Re-format density
unless( $newDftData->{'density'}->{'complete'} ) {
  my $t0 = [gettimeofday];
  print "Exporting density from SCF\n";
  my $errorCode;
  if( $newDftData->{'general'}->{'program'} eq "qe" ) {
    my ($emin, $emax) = QEconvertEnergies( $newDftData->{'structure'}->{'semicore_electrons'});
    $errorCode = QEparseDensityPotential( $newDftData, "density", $emin, $emax );
    if( $errorCode == 0 && abs($emin-$emax) > 0.1 ) {
      $errorCode = QEparseDensityPotential( $newDftData, "sc-density", $emin, $emax );
    }
  } elsif ( $newDftData->{'general'}->{'program'} eq "abi" ) {
     $errorCode = ABIparseDensityPotential( $newDftData, "density" );
  }
  exit $errorCode if( $errorCode != 0 );


  open OUT, ">", "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "%s  %s  %s\n", $commonOceanData->{'structure'}->{'avecs'}[$i][0], 
                                $commonOceanData->{'structure'}->{'avecs'}[$i][1], 
                                $commonOceanData->{'structure'}->{'avecs'}[$i][2];
  }
  close OUT;

  open OUT, ">", "bvecs" or die "Failed to open bvecs\n$!";
  for( my $i = 0; $i < 3; $i++ )
  {
    printf  OUT "%s  %s  %s\n", $commonOceanData->{'structure'}->{'bvecs'}[$i][0], 
                                $commonOceanData->{'structure'}->{'bvecs'}[$i][1], 
                                $commonOceanData->{'structure'}->{'bvecs'}[$i][2];
  }
  close OUT;

  system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0  or die "Failed to run rhoofg.x\n";
  system("wc -l rhoG2 > rhoofg") == 0 or die "$!\n";
  system("sort -n -k 6 rhoG2 >> rhoofg") == 0 or die "$!\n";

  unlink( "avecsinbohr.ipt" );
  unlink( "bvecs" );
  unlink( "rhoG2" );

  $newDftData->{'density'}->{'complete'} = JSON::PP::true;
  $newDftData->{'density'}->{'time'} = tv_interval( $t0 );
  open OUT, ">", "dft.json" or die;
  print OUT $json->encode($newDftData);
  close OUT;
  print "Density export complete\n";
} else {
  $newDftData->{'density'}->{'time'} = $dftData->{'density'}->{'time'};
}
  


# Re-format potential
unless( $newDftData->{'potential'}->{'complete'} ) {
  my $t0 = [gettimeofday];
  print "Exporting potential from SCF\n";
#  my $errorCode = QEparseDensityPotential( $newDftData, "potential" );
  my $errorCode;
  if( $newDftData->{'general'}->{'program'} eq "qe" ) {
     $errorCode = QEparseDensityPotential( $newDftData, "potential", 0, 0 );
  } elsif ( $newDftData->{'general'}->{'program'} eq "abi" ) {
     $errorCode = ABIparseDensityPotential( $newDftData, "potential" );
  }
  exit $errorCode if( $errorCode != 0 );


  $newDftData->{'potential'}->{'complete'} = JSON::PP::true;
  $newDftData->{'potential'}->{'time'} = tv_interval( $t0 );
  open OUT, ">", "dft.json" or die;
  print OUT $json->encode($newDftData);
  close OUT;
  print "Potential export complete\n";
} else {
  $newDftData->{'potential'}->{'time'} = $dftData->{'potential'}->{'time'};
}

unless( $newDftData->{'epsilon'}->{'complete'} ) {
  my $t0 = [gettimeofday];

  my $errorCode;
  if( $newDftData->{'general'}->{'program'} eq "qe" ) {
    $errorCode  = QErunDFPT(  $newDftData );
  } else {
    die "DFPT not enabled for ABINIT yet\n";
  }
  die "Failed to run DFPT\n" if( $errorCode != 0 );

  $newDftData->{'epsilon'}->{'complete'} = JSON::PP::true;
  $newDftData->{'epsilon'}->{'time'} = tv_interval( $t0 );
  open OUT, ">", "dft.json" or die;
  print OUT $json->encode($newDftData);
  close OUT;
  print "Epsilon calculation complete\n";
} else {
  $newDftData->{'epsilon'}->{'time'} = $dftData->{'epsilon'}->{'time'};
  unless ( $newDftData->{'epsilon'}->{'method'} eq "input" ) {
    $newDftData->{'structure'}->{'epsilon'} = $newDftData->{'epsilon'}->{'epsilon'};
  }
}


# Time for NSCF runs
my %sorted_nscf_TodoList;
foreach my $dirname (keys %{$nscf_TodoList}) {
  my $nb = $nscf_TodoList->{$dirname}->{'nbands'};
  my $nk = $nscf_TodoList->{$dirname}->{'kmesh'}[0]
         * $nscf_TodoList->{$dirname}->{'kmesh'}[1]
         * $nscf_TodoList->{$dirname}->{'kmesh'}[2];
  my $val = $nb*$nb*$nb*$nk;
  while( exists $sorted_nscf_TodoList{ $val } ) {
    $val *= (1 + rand()/50.0);
  }
  $sorted_nscf_TodoList{ $val } = $dirname;
#  print "$dirname:  $val\n";
}
foreach my $val (sort {$b <=> $a} keys %sorted_nscf_TodoList ) {
  my $dirname = $sorted_nscf_TodoList{ $val };
  print "$dirname:  $val\n";
}


#foreach my $dirname (keys %{$nscf_TodoList}) {
foreach my $val (sort {$b <=> $a} keys %sorted_nscf_TodoList ) {
  my $dirname = $sorted_nscf_TodoList{ $val };
  print "Running NSCF run for: " . $dirname . "\n";
  my $t0 = [gettimeofday];
  my $errorCode;
  if( $newDftData->{'general'}->{'program'} eq "qe" ) {
     $errorCode = QErunNSCF($newDftData, $nscf_TodoList->{$dirname}, 0 );
  } elsif ( $newDftData->{'general'}->{'program'} eq "abi" ) {
     $errorCode = ABIrunNSCF($newDftData, $nscf_TodoList->{$dirname}, 0 );
  }
  
  exit $errorCode if( $errorCode );
  my $s = $json->encode($newDftData->{'psp'}->{'pphash'});
  foreach ( 'general', 'structure', 'scf' )
  {
    $s .= $json->encode($newDftData->{$_});
  }
  $s .= $json->encode($nscf_TodoList->{$dirname});
  $nscf_TodoList->{$dirname}->{'hash'} = md5_hex( $s );
  $nscf_TodoList->{$dirname}->{'time'} = tv_interval( $t0 );


  $newDftData->{'znscf'}->{ $dirname } = dclone( $nscf_TodoList->{$dirname} );
  open OUT, ">", "dft.json" or die;
  print OUT $json->encode($newDftData);
  close OUT;
}

if( 0 ) {
# Time for SCREENING states
if( $newDftData->{'screen'}->{'enable'} ) {
  # Search for completed runs
  unless( $newDftData->{'screen'}->{'complete'} ) {
    my $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $newDftData->{'screen'}->{'kmesh'}[0],
                    $newDftData->{'screen'}->{'kmesh'}[1], $newDftData->{'screen'}->{'kmesh'}[2],
                    $newDftData->{'screen'}->{'kshift'}[0], $newDftData->{'screen'}->{'kshift'}[1],
                    $newDftData->{'screen'}->{'kshift'}[2];
    if( exists $newDftData->{'znscf'}->{ $dirname } ) {
      if( $newDftData->{'znscf'}->{ $dirname }->{'toldfe'} <= $newDftData->{'screen'}->{'toldfe'} &&
          $newDftData->{'znscf'}->{ $dirname }->{'nbands' } >= $newDftData->{'screen'}->{'nbands'} )
      {
        print "Found previous DFT NSCF run for the screening\n";
        $newDftData->{'screen'} = dclone( $newDftData->{'znscf'}->{ $dirname } );
      }
    }
  }
  unless( $newDftData->{'screen'}->{'complete'} ) {
    my $t0 = [gettimeofday];
    print "Running DFT for screening states\n";

    my $errorCode;
    if( $newDftData->{'general'}->{'program'} eq "qe" ) {
       $errorCode = QErunNSCF($newDftData, $newDftData->{'screen'}, 0 );
    } elsif ( $newDftData->{'general'}->{'program'} eq "abi" ) {
       $errorCode = ABIrunNSCF($newDftData, $newDftData->{'screen'}, 0 );
    }
    
    exit $errorCode if( $errorCode );

    $newDftData->{'screen'}->{'complete'} = JSON::PP::true;
    my $s = $json->encode($newDftData->{'psp'}->{'pphash'});
    foreach ( 'general', 'screen', 'structure', 'scf' )
    {
      $s .= $json->encode($newDftData->{$_});
    }
#    print "$s\n\n\n";
    $newDftData->{'screen'}->{'hash'} = md5_hex( $s );
    $newDftData->{'screen'}->{'time'} = tv_interval( $t0 );

    my $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $newDftData->{'screen'}->{'kmesh'}[0],
                    $newDftData->{'screen'}->{'kmesh'}[1], $newDftData->{'screen'}->{'kmesh'}[2],
                    $newDftData->{'screen'}->{'kshift'}[0], $newDftData->{'screen'}->{'kshift'}[1],
                    $newDftData->{'screen'}->{'kshift'}[2];
    $newDftData->{'znscf'}->{ $dirname } = dclone( $newDftData->{'screen'} );

    open OUT, ">", "dft.json" or die;
    print OUT $json->encode($newDftData);
    close OUT;
    print "DFT for screening states complete\n";
  } #else {
  #  $newDftData->{'screen'}->{'time'} = $dftData->{'screen'}->{'time'};
  #}
}


# Time for BSE final states
unless( $newDftData->{'bse'}->{'complete'} ) {

  my $t0 = [gettimeofday];
  print "Running DFT for BSE basis states\n";

  
  my $errorCode;
  if( $newDftData->{'general'}->{'program'} eq "qe" ) {
    $errorCode = QErunNSCF($newDftData, $newDftData->{'bse'}, 0 );
  } elsif ( $newDftData->{'general'}->{'program'} eq "abi" ) {
    $errorCode = ABIrunNSCF($newDftData, $newDftData->{'bse'}, 0 );
  }

  
  exit $errorCode if( $errorCode );

  $newDftData->{'bse'}->{'complete'} = JSON::PP::true;

  my $s = $json->encode($newDftData->{'psp'}->{'pphash'});
  foreach ( 'general', 'bse', 'structure', 'scf' )
  {
    $s .= $json->encode($newDftData->{$_});
  }
#  print "$s\n\n\n";
  $newDftData->{'bse'}->{'hash'} = md5_hex( $s );
  $newDftData->{'bse'}->{'time'} = tv_interval( $t0 );

  my $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $newDftData->{'bse'}->{'kmesh'}[0],
                  $newDftData->{'bse'}->{'kmesh'}[1], $newDftData->{'bse'}->{'kmesh'}[2],
                  $newDftData->{'bse'}->{'kshift'}[0], $newDftData->{'bse'}->{'kshift'}[1],
                  $newDftData->{'bse'}->{'kshift'}[2];
  $newDftData->{'znscf'}->{ $dirname } = dclone( $newDftData->{'bse'} );

  open OUT, ">", "dft.json" or die;
  print OUT $json->encode($newDftData);
  close OUT;
  print "DFT for BSE final states complete\n";
} else {
  $newDftData->{'bse'}->{'time'} = $dftData->{'bse'}->{'time'};
  if( $newDftData->{'bse'}->{'con_start'} != $dftData->{'bse'}->{'con_start'} 
        && ( defined( $newDftData->{'bse'}->{'con_start'}) || defined($dftData->{'bse'}->{'con_start'} ) ) ) {
#  if( ( $newDftData->{'bse'}->{'con_start'} != $dftData->{'bse'}->{'con_start'} 
#        && ( defined( $newDftData->{'bse'}->{'con_start'}) || defined($dftData->{'bse'}->{'con_start'} ) ) ) 
#   || ($newDftData->{'bse'}->{'val_stop'} != $dftData->{'bse'}->{'val_stop'}
#        && ( defined( $newDftData->{'bse'}->{'val_stop'}) || defined($dftData->{'bse'}->{'val_stop'} ) ) ) ) {
    my $errorCode;
    if( $newDftData->{'general'}->{'program'} eq "qe" ) {
      $errorCode = QErunParseEnergies( $newDftData, $newDftData->{'bse'}, 0 );
    }
    exit $errorCode if( $errorCode );
    open OUT, ">", "dft.json" or die;
    print OUT $json->encode($newDftData);
    close OUT;
  }
}
}

if( -e "redo_energies" ) {
  print "Re-calculating energies\n";
  foreach (@{$newDftData->{'bse'}->{'directories'}}) {
    QErunParseEnergies( $newDftData, $newDftData->{'znscf'}->{$_}, 0 );
  }
  QErunParseEnergies( $newDftData, $newDftData->{'znscf'}->{$newDftData->{'screen'}->{'directories'}[0]}, 0 );
#  QErunParseEnergies( $newDftData, $newDftData->{'bse'}, 0 );
#  QErunParseEnergies( $newDftData, $newDftData->{'screen'}, 0 );
}

# touch up Fermi if insulator
if( $newDftData->{'general'}->{'occopt'} == 1 && $newDftData->{'general'}->{'program'} eq "qe" )
{
  my $low = $newDftData->{'scf'}->{'lowest'};
  my $high = $newDftData->{'scf'}->{'highest'};
  foreach (@{$newDftData->{'bse'}->{'directories'}} ) {
    $low = $newDftData->{'znscf'}->{$_}->{'lowest'} if ( $newDftData->{'znscf'}->{$_}->{'lowest'} < $low );
    $high = $newDftData->{'znscf'}->{$_}->{'highest'} if ( $newDftData->{'znscf'}->{$_}->{'highest'} > $high );
  }
#  $low = $newDftData->{'bse'}->{'lowest'} if ( $newDftData->{'bse'}->{'lowest'} < $low );

#  $high = $newDftData->{'bse'}->{'highest'} if ( $newDftData->{'bse'}->{'highest'} > $high );

#  if( $newDftData->{'screen'}->{'enable'} )
##  {
#    $low = $newDftData->{'screen'}->{'lowest'} if ( $newDftData->{'screen'}->{'lowest'} < $low );
#    $high = $newDftData->{'screen'}->{'highest'} if ( $newDftData->{'screen'}->{'highest'} > $high );
#  }

  $newDftData->{'scf'}->{'fermi'} = ($low + $high)/2;
  open OUT, ">", "dft.json" or die;
  print OUT $json->encode($newDftData);
  close OUT;
  print "DFT for BSE final states complete\n";
}

# Fix up brange
if( $newDftData->{'screen'}->{'enable'} ) {
  my $d = $newDftData->{'screen'}->{'directories'}[0];
  for( my $i = 0; $i < 3; $i++ ) {
    $newDftData->{'screen'}->{'brange'}[$i] = $newDftData->{'znscf'}->{ $d }->{'brange'}[$i];
  }
}
$newDftData->{'bse'}->{'brange'}[0] = 
      $newDftData->{'znscf'}->{$newDftData->{'bse'}->{'directories'}[-1]}->{'brange'}[0];
$newDftData->{'bse'}->{'brange'}[1] = 
      $newDftData->{'znscf'}->{$newDftData->{'bse'}->{'directories'}[-1]}->{'brange'}[1];
$newDftData->{'bse'}->{'brange'}[2] = 
      $newDftData->{'znscf'}->{$newDftData->{'bse'}->{'directories'}[0]}->{'brange'}[2];

if( $commonOceanData->{'bse'}->{'mimic_exciting_bands'} ) {
  $newDftData->{'bse'}->{'brange'}[3] = $newDftData->{'bse'}->{'brange'}[2] 
                                      + $commonOceanData->{'bse'}->{'mimic_exciting_bands'} - 1;
}

if( exists $newDftData->{'znscf'}->{$newDftData->{'bse'}->{'directories'}[-1]}->{'max_occ_band'} ) {
  $newDftData->{'bse'}->{'max_occ_band'} = $newDftData->{'znscf'}->{$newDftData->{'bse'}->{'directories'}[-1]}->{'max_occ_band'};
} else {
  $newDftData->{'bse'}->{'max_occ_band'} = $newDftData->{'bse'}->{'brange'}[1];
}
if( exists $newDftData->{'znscf'}->{$newDftData->{'bse'}->{'directories'}[0]}->{'min_unocc_band'} ) {
  $newDftData->{'bse'}->{'min_unocc_band'} = $newDftData->{'znscf'}->{$newDftData->{'bse'}->{'directories'}[0]}->{'min_unocc_band'};
} else {
  $newDftData->{'bse'}->{'min_unocc_band'} = $newDftData->{'bse'}->{'brange'}[2];
}



# Grab times and hashes
$newDftData->{'screen'}->{'time'} = $newDftData->{'znscf'}->{ $newDftData->{'screen'}->{'directories'}[0] }->{'time'};
$newDftData->{'screen'}->{'hash'} = $newDftData->{'znscf'}->{ $newDftData->{'screen'}->{'directories'}[0] }->{'hash'};
$newDftData->{'bse'}->{'time'} = 0;
$newDftData->{'bse'}->{'hash'} = '';
foreach (@{$newDftData->{'bse'}->{'directories'}} ) {
  $newDftData->{'bse'}->{'time'} += $newDftData->{'znscf'}->{ $_ }->{'time'};
  $newDftData->{'bse'}->{'hash'} .= $newDftData->{'znscf'}->{ $_ }->{'hash'};
  if( $newDftData->{'screen'}->{'directories'}[0] eq $_ ) {
    $newDftData->{'screen'}->{'time'} = 0;
  } 
}

my ( $endSeconds, $endMicroseconds) = gettimeofday;
print localtime() .  "\n";

my $elapsedSeconds = $endSeconds - $startSeconds;
my $elapsedMicroseconds = $endMicroseconds - $startMicroseconds;
if( $elapsedMicroseconds < 0 ) {
  $elapsedMicroseconds += 1000000;
  $elapsedSeconds -= 1;
}
$newDftData->{'time_script'} = sprintf "%i.%06i", $elapsedSeconds, $elapsedMicroseconds;
$newDftData->{'time'} = 0;
foreach my $sec ( 'scf', 'density', 'potential', 'bse', 'screen', 'epsilon' ) {
  printf "Time %s: %f\n", $sec, $newDftData->{$sec}->{'time'};
  $newDftData->{'time'} += $newDftData->{$sec}->{'time'};
}
open OUT, ">", "dft.json" or die;
print OUT $json->encode($newDftData);
close OUT;
print "DFT section is complete\n\n";

exit 0;

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
    if( ref( $commonRef->{ $t } ) eq '' )
    {
#      print "$commonRef->{ $t } ---\n";
      $newRef->{ $t } = $commonRef->{ $t };
    }
    else
    {
      $newRef->{ $t } = dclone $commonRef->{ $t };
#      recursiveTouch( $newRef->{ $t }  );
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

sub copyNoCompare
{
  my $newRef = $_[0];
  my $commonRef = $_[1];
  my @tags = @{$_[2]};

  foreach my $t (@tags)
  {
    if( ref( $commonRef->{ $t } ) eq '' )
    {
      $newRef->{ $t } = $commonRef->{ $t };
    }
    else
    {
      $newRef->{ $t } = dclone $commonRef->{ $t };
    }
  }
}


sub recursiveTouch
{
  my $newRef = $_[0];
  if( ref( $newRef ) eq 'ARRAY' )
  {
    for( my $i = 0; $i < scalar @{ $newRef }; $i++ )
    {
      recursiveTouch( @{$newRef}[$i] );
    }
  }
  else
  {
    $newRef*=1 if( looks_like_number( $newRef ) );
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
    unless( ref( $oldRef ) eq 'HASH' ) {
      $complete->{'complete'} = JSON::PP::false;
      return;
    }
    if( scalar keys %{$newRef} != scalar keys %{$oldRef} ) {
      $complete->{'complete'} = JSON::PP::false;
      return;
    }
    foreach my $key (keys %{$newRef} ) {
      unless( exists $oldRef->{$key} ) {
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
      $complete->{'complete'} = JSON::PP::false unless( $newRef == $oldRef );
    }
    else
    {
      $complete->{'complete'} = JSON::PP::false unless( $newRef eq $oldRef );
    }
  }
}

sub checkSetGamma
{
  my $hashRef = $_[0];

  my $minQ = 0.0000001;

  if( $hashRef->{'kmesh'}[0] == 1 && $hashRef->{'kmesh'}[1] == 1 && $hashRef->{'kmesh'}[1] == 1 
      && abs($hashRef->{'kshift'}[0]) < $minQ && abs($hashRef->{'kshift'}[1]) < $minQ 
      && abs($hashRef->{'kshift'}[0]) < $minQ )
  {
    $hashRef->{'isGamma'} = JSON::PP::true;
  } 
  else
  {
    $hashRef->{'isGamma'} = JSON::PP::false;
  }

}

sub shiftKpointsByPhoton 
{
  my $hashRef = $_[0];
  my $photon_q = $_[1];
  
  for( my $i = 0; $i < 3; $i ++ ) {
    my $kactual = $hashRef->{'kshift'}[$i]/$hashRef->{'kmesh'}[$i];
    $kactual -= @{$photon_q}[$i];
    while( $kactual < 0 ) { $kactual += 1; }
    while( $kactual >= 1 ) { $kactual -= 1; }
    $kactual *= $hashRef->{'kmesh'}[$i];
    $hashRef->{'kshift'}[$i] = $kactual;
  }
}
