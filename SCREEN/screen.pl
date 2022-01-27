#!/usr/bin/perl
# Copyright (C) 2010, 2013 - 2021 OCEAN collaboration
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
  $newScreenData->{'density'} = {};
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

  $newScreenData->{'model'} = {};
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
  
  print "MODEL: " . $newScreenData->{'model'}->{'complete'} . "\n";

  
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

  my @screenList = ( "all_augment", "augment", "convertstyle", "grid", "inversionstyle", "kmesh", 
                     "kshift", "mode", "nbands", "shells", "final", "core_offset" );
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
  

  #### RPA calculation

  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;

  writeExtraFiles( $newScreenData->{'structure'}, $newScreenData->{'screen'}, $newScreenData->{'general'} );

  $newScreenData->{'combine'} = {};
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

  $newScreenData->{'offset'} = {};
  $newScreenData->{'offset'}->{'complete'} = JSON::PP::true;
  unless( $newScreenData->{'density'}->{'complete'} && $newScreenData->{'screen'}->{'complete'} 
        && exists $screenData->{'offset'}->{'complete'} && $screenData->{'offset'}->{'complete'} ) {
    $newScreenData->{'offset'}->{'complete'} = JSON::PP::false;
  }

  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;

  unless( $newScreenData->{'density'}->{'complete'} )
  {
    my $errorCode = runDensityAverage( $newScreenData );
    if( $errorCode )
    { 
      die "Failed in runDensityAverage with code: $errorCode\n";
    }

    cleanAverage( $newScreenData->{'general'} );
    $newScreenData->{'density'}->{'complete'} = JSON::PP::true;
  }
  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;

  unless( $newScreenData->{'model'}->{'complete'} )
  {
#    print "MODEL: " . $newScreenData->{'model'}->{'complete'} . "\n";
    if( $newScreenData->{'model'}->{'flavor'} eq 'SLL') {
      runVhommod($newScreenData->{'model'}->{'SLL'}) 
    } else {
      die "Unrecognized model flavor!\n";
    }

    $newScreenData->{'model'}->{'complete'} = JSON::PP::true;
  }
  
  open OUT, ">", "screen.json" or die;
  print OUT $json->encode($newScreenData);
  close OUT;


  unless( $newScreenData->{'screen'}->{'complete'} )
  {
    print "SCREEN\n";
    grabOPF();

    $newScreenData->{'screen'}->{'grid2'} = dclone $newScreenData->{'screen'}->{'grid'};
    buildMKRB( $newScreenData->{'screen'}, $newScreenData->{'general'} );
    grabAngFile( $newScreenData->{'screen'} );
    prepWvfn(  $newScreenData->{'screen'},  $newScreenData->{'general'} );

    runScreenDriver( $newScreenData->{'computer'} );

    $newScreenData->{'screen'}->{'complete'} = JSON::PP::true;

    cleanScreen( $newScreenData->{'general'}, $newScreenData->{'screen'} );

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
    print "COMBINE\n";

#    cleanScreen( $newScreenData->{'general'}, $newScreenData->{'screen'} );
    finishCorePotentials( $newScreenData->{'general'}, $newScreenData->{'screen'} );


    $newScreenData->{'combine'}->{'complete'} = JSON::PP::true;
    open OUT, ">", "screen.json" or die;
    print OUT $json->encode($newScreenData);
    close OUT;

  }

  unless( $newScreenData->{'offset'}->{'complete'} )
  {
    print "OFFSET\n";

    runCoreOffset( $newScreenData->{'screen'}, $newScreenData);

    $newScreenData->{'offset'}->{'complete'} = JSON::PP::true;

    open OUT, ">", "screen.json" or die;
    print OUT $json->encode($newScreenData);
    close OUT;
  }

  exit 0;
}



my @CommonFiles = ("znucl", "opf.hfkgrid", "opf.fill", "opf.opts", "pplist", "screen.shells", 
                   "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution",  
                   "screen.k0", "scfac", "core_offset", "dft", "avecsinbohr.ipt", "coord", 
                   "nspin", "prefix", "work_dir" );

my @CommonFiles2 = ("para_prefix", "calc" );

my @ScreenFiles = ("screen.grid.scheme", "screen.grid.rmode", "screen.grid.ninter", 
                   "screen.grid.shells", "screen.grid.xyz", "screen.grid.rmax", "screen.grid.ang",
                   "screen.lmax", "screen.grid.nb", "screen.grid.nr", "screen.final.rmax", 
                   "screen.final.dr", "screen.legacy", "screen.model.dq", "screen.model.qmax", 
                   "screen.augment", "screen.wvfn", "screen.convertstyle", "screen.inversionstyle", 
                   "screen.mode", "screen.grid.deltar" );

my @DenDipFiles = ("rhoofg", "bvecs", "efermiinrydberg.ipt", "xmesh.ipt", "epsilon");
my @DenDipFiles2 = ( "masterwfile", "listwfile", "enkfile", "kmesh.ipt", "brange.ipt" );

my @ExtraFiles = ("specpnt", "Pquadrature", "hqp", "lqp", "EvenQuadHalf.txt" );

my @DFTFiles = ( "potofr" );

my $runSCREEN = 1;
if (-e "../PREP/SCREEN/old" && -e "done" ) {
  $runSCREEN = 0;
  foreach (@CommonFiles) {
#    if (`diff -q $_ ../Common/$_`) {
    if( compare("$_", "../Common/$_" ) != 0 )   # should get diff or non-exist
    {
      print "$_ differs\n";
      $runSCREEN = 1;
      last;
    }
  }
  if( $runSCREEN == 0 )
  {
    foreach (@ScreenFiles)
    {
#      if (`diff -q $_ ../Common/$_`) {
    if( compare("$_", "../Common/$_" ) != 0 )   # should get diff or non-exist
    {
        print "$_ differs\n";
        $runSCREEN = 1;
        last;
      }
    }
  }
}
else
{
  if( -e "done" )
  {
    print "PREP was updated. Will re-run\n";
  }
  else
  {
    print "Screening wasn't previously completed\n";
  }
}

# Secondary check
if( $runSCREEN == 0 )
{
  open IN, "screen.mode" or die "Failed to open screen.mode\n$!";
  if( <IN> =~ m/grid/i )
  {
    if( $runSCREEN == 0 ) 
    {
      if( compare("xmesh.ipt", "../Common/xmesh.ipt" ) != 0 )
      {
        print "Using screen.mode = grid and xmesh.ipt has changed\n";
        $runSCREEN = 1;
      }
    }
  }
  close IN;
}

# OPF check
if( $runSCREEN == 0 )
{
  open IN, "screen.augment" or die "Failed to open screen.augment\n";
  if( <IN> =~ m/t/i )
  {
    unless( -e "../OPF/old" )
    {
      print "OPFs were updated recently. Will re-run\n";
      $runSCREEN = 1;
    }
  }
}

if( $runSCREEN == 0 )
{
  if( compare("epsilon", "../PREP/epsilon" ) != 0 )
  {
    print "epsilon has change";
    $runSCREEN = 1;
  }
}

if ($runSCREEN == 0 ) {
  print "Nothing new needed for SCREEN stage\n";
  open OUT, ">", "old" or die;
  print OUT "1\n";
  close OUT;  
  exit 0;
}

unlink "done";
unlink "old";


foreach (@CommonFiles) {
  copy( "../Common/$_", $_ ) or die "Failed to get $_ from Common/\n$!";
}

copy( "screen.k0", "k0.ipt") or die "$!";

my %screen_data_files = {};
foreach my $filename (@ScreenFiles)
{
  copy( "../Common/$filename", $filename ) or die "Failed to get $filename from Common\n$!";
  open IN, $filename or die "filename: $!";
  my $string;
    while( my $a = <IN> )
    {
      $string .= $a;
    }
    chomp $string;

    close IN;
    $filename =~ m/screen\.(\w+\.?\w+)/ or die;
    my $store_name = $1;
    $screen_data_files{ "$store_name" } = $string;
}

# Attempt to parse and copy in angular grid file
#if( $screen_data_files{ 'grid.ang' } =~ m/(\w+)\s+(\d+)/ )
#{
#  my $angularGridFile = $1 . '.' . $2;

my @angList    = split ' ', $screen_data_files{'grid.ang'};
foreach (@angList)
{
  my $angularGridFile = 'specpnt.' . $_;
  if( -e  "$ENV{'OCEAN_BIN'}/$angularGridFile" )
  {
    copy( "$ENV{'OCEAN_BIN'}/$angularGridFile", "$angularGridFile" );
  }
  else
  {
    print "Couldn't find requested angular grid file: $angularGridFile\n";
  }
}
  

# Need calc
foreach( @CommonFiles2 )
{
  copy( "../Common/$_", "$_" ) or die "Failed to copy $_ from Common\n$!";
}

my $valenceGrid;
open CALC, "calc" or die "Failed to open calc\n";
if( <CALC> =~ m/VAL/i )
{
  if( $screen_data_files{ "mode" } =~ m/grid/i )
  {
    print "Running valence grid\n";
    $valenceGrid = 1;
  }
  else
  {
    print "No (core-level) SCREEN calc for valence run\n";
    close CALC;
    exit 0;
  }
}
else
{
  # if we aren't running valence set valenceGrid to zero no matter what
  $valenceGrid = 0;
}
close CALC;

foreach (@DFTFiles) {
  copy( "../DFT/$_", $_ ) or die "Failed to get $_ from DFT/\n$!";
}
foreach (@DenDipFiles) {
  copy( "../PREP/$_", $_ ) or die "Failed to get $_ from PREP/\n$!";
}
foreach (@DenDipFiles2) {
  copy( "../PREP/SCREEN/$_", $_ ) or die "Failed to get $_ from PREP/SCREEN/\n$!";
}

foreach (@ExtraFiles) {
  copy( "$ENV{'OCEAN_BIN'}/$_", $_ ) or die "Failed to get $_ from $ENV{'OCEAN_BIN'}/\n$!";
}

my $dft = `cat dft`;
if( $dft =~ m/qe/i )
{
  `ln -s ../DFT/Out SCF`;
  `ln -s ../DFT/SCREEN/Out .`;
  copy( "../DFT/SCREEN/gamma", "gamma" ) if -e "../DFT/SCREEN/gamma";
}
elsif( $dft =~ m/abi/i )
{
  `ln -s ../DFT/SCREEN/RUN0001_WFK .`;
}

# Detect qe version ( or, in the future, abinit )
if( $screen_data_files{ 'wvfn' } =~ m/qe/ || $screen_data_files{ 'wvfn' } =~ m/new/ )
{
  open IN, "prefix" or die "Failed to open prefix\n$!";
  my $prefix = <IN>;
  chomp( $prefix );
  if( -e "Out/$prefix.save/data-file.xml" )
  {
    print "Detected QE54-style DFT run\n";
    $screen_data_files{ 'wvfn' } = "qe54";
  }
  elsif( -e "Out/$prefix.save/data-file-schema.xml" )
  {
    print "Detected QE62-style DFT run\n";
    $screen_data_files{ 'wvfn' } = "qe62";
    copy "../DFT/SCREEN/eig62.txt", "eig62.txt";
    copy "../DFT/SCREEN/enkfile", "enkfile_raw";
  }
  else
  {
    print "WARNING! Failed to detect style of QE output\nWill attempt to continue\n";
  }
}
open WVFN, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
print WVFN $screen_data_files{ 'wvfn' } . "\n";
close WVFN;

open LISTW, "listwfile" or die "Failed to open listwfile\n";
while (<LISTW> ) {
  m/(\d+)\s+(\S+)/ or die "Malformed listwfile\n";
  `ln -sf ../PREP/SCREEN/$2 $2`;
}

unless( -e 'clipbands' ) {
print "Creating clipbands\n";
open BRANGE, "brange.ipt" or die "Failed to open brange.ipt\n";
open CLIPS, ">clipbands" or die "Failed to open clipbands for writing\n";
<BRANGE> =~ m/(\d+)\s+\d+/ or die "Malformed brange\n";
print CLIPS "$1   ";
<BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Malformed brange\n";
print CLIPS "$2\n";
close BRANGE;
close CLIPS;
}

my $para_prefix = "";
if( open PARA_PREFIX, "para_prefix" )
{
  $para_prefix = <PARA_PREFIX>;
  chomp($para_prefix);
  close( PARA_PREFIX);
}

###################################

# Setup
###################################
print "Running SCREEN Setup\n";
system("$ENV{'OCEAN_BIN'}/pawsetup.x") == 0 or die "Failed to run pawsetup.x\n";

`ln -sf ../OPF/zpawinfo zpawinfo`;
###################################



# shells
##################################
open SHELLS, "screen.shells" or die "Failed to open screen.shells\n";
my $numshells = 0;
my $allshells = '';
while (<SHELLS>) {
  chomp;
  $allshells .= $_ ." ";
}
close SHELLS;
my @rads = split ' ', $allshells;
$numshells = $#rads + 1;
open SHELLS, ">shells" or die "Failed to open shells for writing\n";
print SHELLS "$numshells\n";
for( my $i = 0; $i < scalar @rads; $i++ )
{
  print SHELLS "$rads[$i]\n";
}
#print SHELLS "$allshells\n";
close SHELLS;

######################################

print "Starting xipp section\n";

if( -e "$ENV{'OCEAN_BIN'}/mpi_avg.x" )
{
  print "Running mpi_avg.x\n";
  system("$para_prefix $ENV{'OCEAN_BIN'}/mpi_avg.x") == 0 or die "$!\nFailed to run mpi_avg.x\n";
}
else
{
  print "Running avg.x\n";
  system("$ENV{'OCEAN_BIN'}/avg.x") == 0 or die "$!\nFailed to run avg.x\n";
}


# Pre-process some screen params here
$screen_data_files{'final.rmax'} =~ m/(\d+\.?\d?)/ or die "Failed to parse screen.final.rmax\n";
$screen_data_files{'final.rmax'} = $1;
$screen_data_files{'final.dr'} =~ m/(\d+\.?\d*)/ or die "Failed to parse screen.final.dr\n";
$screen_data_files{'final.dr'} = $1;
my $final_nr = sprintf("%.i", $screen_data_files{'final.rmax'} / $screen_data_files{'final.dr'} );
print "Final grid: " . $screen_data_files{'final.dr'} . " " . $final_nr . "\n";

if( $valenceGrid == 1 ) # valence grid, must use screen_driver.x
{

  ###MOVE up to universal at some point ###
  my %requiredSpecpnt;

  # Start putting together defaults
  #By default 2 shells, one for OPF and one for the rest
  my $ninter = 2;
  my @angList    = split ' ', $screen_data_files{'grid.ang'};
  $ninter = scalar @angList if( scalar @angList > $ninter );

  my @rmodeList  = split ' ', $screen_data_files{'grid.rmode'};
  $ninter = scalar @rmodeList if( scalar @rmodeList > $ninter );

  my @deltarList = split ' ', $screen_data_files{'grid.deltar'};
  $ninter = scalar @rmodeList if( scalar @deltarList > $ninter );

  my @shellsList = split ' ', $screen_data_files{'grid.shells'};
  $ninter = scalar @shellsList if( scalar @shellsList > $ninter );

  # 
  if( $shellsList[0] <= 0 )
  { 
    open IN, "hfinlist" or die;
    <IN> =~ m/^\s*\S+\s+(\d+)/ or die;
    my $zee = $1;
    close IN;
    my $zeeName = sprintf("z%03i",$zee);
    open IN, "zpawinfo/radfile$zeeName" or die;
    <IN> =~ m/^\s+(\S+)/ or die;
    $shellsList[0] = $1;
  }

  unless( $rmodeList[0] =~ m/legendre/ || $rmodeList[0] =~ m/uniform/ )
  {
    $rmodeList[0] = 'legendre';
  }

  $angList[0] = 5 if( $angList[0] < 0 );
  $deltarList[0] = 0.2 if( $deltarList[0] < 0 );

  for( my $i = 1; $i < $ninter; $i++ )
  {
    $angList[$i] = $angList[$i-1] if( scalar @angList <= $i );
    $angList[$i] = 7 if( $angList[$i] < 0 || ! $angList[$i] =~ m/\d/ ); #|| scalar @angList <= $i );
    print "$deltarList[$i] !!  ";
    $deltarList[$i] = $deltarList[$i-1] if( $deltarList[$i] < 0 || scalar @deltarList <= $i );
    print "$deltarList[$i]\n";
    $rmodeList[$i] = 'uniform' unless( $rmodeList[$i] =~ m/legendre/ || $rmodeList[$i] =~ m/uniform/ );

    if( $shellsList[$i] < 0 || $shellsList[$i] <= $shellsList[$i-1]
                            || $shellsList[$i] >= ( $screen_data_files{'grid.rmax'} - $deltarList[$i]/2 ) )
    {
      print "Ran out of shells, might be making fewer than expected!\n";
      $shellsList[$i] = $screen_data_files{'grid.rmax'};
      $ninter = $i+1;
      last;
    }
  }
  print "$ninter\n\n";

  for( my $i = 0; $i < $ninter; $i++ )
  {
    $requiredSpecpnt{"$angList[$i]"} = 1;
  }
  foreach my $key (sort(keys %requiredSpecpnt))
  {
    print "specpnt.$key\n";
    copy( "$ENV{'OCEAN_BIN'}/specpnt.$key", "specpnt.$key" ) or die "Failed to get specpnt.$key\n";
  }

  open MKRB, ">", "mkrb_control" or die "Failed to open mkrb_control for writing\n$!";
#    print MKRB "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $screen_data_files{'grid.ninter'}\n" 
#             . "$screen_data_files{'grid.scheme'} $screen_data_files{'grid.rmode'}\n";
#    close MKRB;
  print MKRB "$screen_data_files{'grid.rmax'} $ninter\n";
  for( my $i = 0; $i < $ninter; $i++ )
  {
    print MKRB "$rmodeList[$i] $shellsList[$i] $deltarList[$i] $angList[$i] specpnt\n";
  }
  close MKRB;



#  copy "specpnt", "specpnt.5";
#  open MKRB, ">", "mkrb_control" or die "Failed to open mkrb_control for writing\n$!";
##  print MKRB "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $screen_data_files{'grid.ninter'}\n"
##           . "$screen_data_files{'grid.scheme'} $screen_data_files{'grid.rmode'}\n";
#  print MKRB "$screen_data_files{'grid.rmax'} 1\n" 
#           . "$screen_data_files{'grid.rmode'}  $screen_data_files{'grid.rmax'} 0.1 5 specpnt\n";
#  close MKRB;

  # Make directory structure
  open XMESH, "xmesh.ipt" or die "Failed to open xmesh.ipt\n$!";
  <XMESH> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse xmesh.ipt\n$_";
  my $nx = $1*$2*$3;
#  my @xmesh = ( $1, $2, $3 );
  close XMESH;

  for( my $i = 1; $i <= $nx; $i++ )
  {
    my $sitename = sprintf("x%06i", $i );
    unless( -d $sitename )
    {
      mkdir $sitename or die "Failed to make directory $sitename\n$!";
    }
    foreach my $rad (@rads)
    {
      my $radName = sprintf("zR%03.2f",$rad);
      mkdir catfile( $sitename, $radName );
    }
  }
  
  my $allaug = 0;
  if( -e "screen.allaug" ) {
    open IN, "screen.allaug" or die;
    if( <IN> =~ m/t/i ) {
      $allaug = 1;
    }
    close IN;
  }
  open ZEELIST, ">", "zeelist" or die "Failed to open zeelist for writing\n$!";
  if( $allaug == 0 ) {
    print ZEELIST "0\n";
  }
  else
  { print ZEELIST "1\n9"; }
  close ZEELIST;

  system("$para_prefix $ENV{'OCEAN_BIN'}/screen_driver.x > screen_driver.log") == 0 or die "$!\nFailed to run screen_driver.x\n";
  print "screen_driver.x done\n";

  system("$ENV{'OCEAN_BIN'}/vhommod.x") == 0 or die "$!\nFailed to run vhommod.x\n";

  # Begin file cleaning
  my @siteFiles = ("chi0", "grid", "chi", "avg" );
  my @potFiles = ( "vind", "vind0", "nind", "nind0" );
  my @sitelist;
  for( my $i = 1; $i <= $nx; $i++ )
  {
    my $siteName = sprintf("%06i", $i );
    foreach my $fileName (@siteFiles)
    {
      my $newFileName = catfile( "x$siteName", $fileName );
      move( "${fileName}${siteName}", $newFileName ) or die "${fileName}${siteName}\n$newFileName\n$!";
    }

    foreach my $rad (@rads)
    {
      my $radName = sprintf("R%03.2f",$rad);
      my $newFileName = catfile( "x$siteName", "reopt.$radName" );
      move( "reopt${siteName}.$radName", $newFileName );
    }
  
    
    foreach my $rad (@rads)
    {
      my $radName = sprintf("zR%03.2f",$rad);
      foreach my $fileName (@potFiles)
      {
        my $fullFileName = "x" . $siteName . "." . $radName . $fileName;
        my $target = catfile( "x$siteName", $radName, $fileName );
        print "$fullFileName\t\t$target\n";
        move( $fullFileName, $target ) or die "Failed to move $fullFileName to $target\n$!";
      }
    }
  }
  #Files have been cleaned and put in directories
  

  # Now build total potential
  my @avec;
  open AVEC, "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i ++ )
  {
    <AVEC> =~ m/(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/ or die "Couldn't parse avecsinbohr.ipt$_\n";
    $avec[0][$i] = $1;
    $avec[1][$i] = $2;
    $avec[2][$i] = $3;
  }
  close AVEC;

  my $volume = $avec[0][0] * ($avec[1][1] * $avec[2][2] - $avec[2][1] * $avec[1][2] )
             - $avec[1][0] * ($avec[0][1] * $avec[2][2] - $avec[2][1] * $avec[0][2] )
             + $avec[2][0] * ($avec[0][1] * $avec[1][2] - $avec[1][1] * $avec[0][2] );

  # 3/(4 pi ) Volume/nx
  my $valRadius = ( 0.238732414637843007 * $volume / $nx ) ** (1/3);

  my @barePotential;
  my @bareRadius;
  my $smallDelta = $screen_data_files{'final.dr'}/10;  
  my $saveRadius;
  for( $saveRadius = 0.0; $saveRadius <= $valRadius; $saveRadius += $smallDelta )
  {
    push @bareRadius, $saveRadius;
    my $pot = -( 0.5 / $valRadius ) * ( 3.0 - $saveRadius**2/ ( $valRadius**2 ) );
    push @barePotential, $pot;
  }
    
  for( my $radius = $saveRadius; $radius <= $screen_data_files{'final.rmax'}; $radius += $smallDelta )
  {
    push @bareRadius, $radius;
    my $pot = -1.0/$radius;
    push @barePotential, $pot;
  }

  open POT, ">", "valencePotential" or die "$!";
  for( my $i = 0; $i < scalar @barePotential; $i ++ )
  {
    print POT $bareRadius[$i] . "   " . $barePotential[$i] . "\n";
  }
  close POT;
  
  for( my $j = 1; $j <= $nx; $j++ )
  {
    my $currentSite = sprintf("x%06i", $j );

    my @reoptArray;
    foreach my $rad (@rads)
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
    #     no aug or fake-aug 
    #  3. write out

    for( my $r = 0; $r < scalar @rads; $r++ )
    {
      my @vindRad; my @vindPot;

      my $radName = sprintf("zR%03.2f",$rads[$r]) ;
      print "\t\t$radName\t$rads[$r]\t$r\n";

      my $vindName = catfile( $currentSite, $radName, "vind" );
      open IN, "<", $vindName or die "Failed to open $vindName\n$!";
      while( my $line = <IN> )
      {
        $line =~ m/(\d+\.\d+[Ee][+-]\d+)\s+(-?\d+\.\d+[Ee][+-]\d+)/ or die "Failed to parse $vindName\n";
        push @vindRad, $1;
        push @vindPot, $2;
      }
      close IN;

      my $rundir = catfile( $currentSite, $radName );
      chdir $rundir or die;

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

      $len = scalar @barePotential;
      print $len . "\n";
      print OUT "1 2\n$len\n";
      for( my $i = 0; $i < $len; $i++ )
      {
        print OUT "$bareRadius[$i]  $barePotential[$i]\n";
      }

      print OUT ".false.\n$screen_data_files{'final.dr'} $final_nr\n";
      close OUT;
      system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt") == 0 or die;

      # back out 2 levels
      $rundir = catfile( updir(), updir() );
      chdir $rundir;    

    }
  }

}
else
{
  if( $screen_data_files{ 'legacy' } != 0 )  #legacy-style run for core
  {
    open HFINLIST, "hfinlist" or die "Failed to open hfinlist\n";
    print "Running $ENV{'OCEAN_BIN'}/vhommod.x\n";
    system( "$ENV{'OCEAN_BIN'}/vhommod.x" ) == 0 or die "Failed to run vhommod.x\n$!";

    my $rad;
    my $edgename; 
    my $hfinline; my $ppfilename; my $znucl; my $nnum; my $lnum; my $elname; my $elnum;
    while ($hfinline = <HFINLIST>) {
    #  print $hfinline;
      ($hfinline =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)/) or die "Malformed hfinlist\t$1 $2 $3 $4 $5 $6\n";
      $ppfilename = $1;
      $znucl = $2;
      $nnum = $3;
      $lnum = $4;
      $elname = $5;
      $elnum = $6;

      $edgename = sprintf("z%2s%04i/n%02il%02i", $elname, $elnum, $nnum, $lnum);
      print "$edgename\n";
      `mkdir -p $edgename` == 0 or die "Failed to make dir $edgename\n";

      my $avden =  sprintf("avg%2s%04i",$elname,$elnum);
      copy( $avden,  "avden" ) or die "Failed to copy density $avden\n$!";


      my $edgename2 = sprintf("z%03in%02il%02i",$znucl, $nnum, $lnum);
    #  `mkdir -p $edgename` == 0 or die "Failed to make dir $edgename\n";


      foreach $rad (@rads) {
        my $fullrad = sprintf("%03.2f",$rad);
    #      `echo "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $elname $elnum" | $ENV{'OCEAN_BIN'}/mkrbfile.x`;
          `echo "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $screen_data_files{'grid.ninter'}" > mkrb_control`;
          `echo "$screen_data_files{'grid.scheme'} $screen_data_files{'grid.rmode'}" >> mkrb_control`;
          `echo 1 >> mkrb_control`;
          `echo "$elname $elnum" >> mkrb_control`;
          system( " $ENV{'OCEAN_BIN'}/mkrbfile_mult.x" ) == 0 or die "Failed to run mkrbfile_mult.x\n$!";
          `mkdir -p ${edgename}/zRXT${fullrad}`;
          `mkdir -p ${edgename}/zRXF${fullrad}`;
          `mkdir -p ${edgename}/zRXS${fullrad}`;
          chdir "$edgename";
          `ln -s -f zRXT${fullrad} zR${fullrad}`;
          chdir "../../";
          copy( "zpawinfo/vc_bare${edgename2}R${fullrad}", "tmp" ) 
              or die "Failed to copy zpawinfo/vc_bare${edgename2}R${fullrad}\n$!";
          `wc tmp > vpert`;
          `cat tmp >> vpert`;
          system("$ENV{'OCEAN_BIN'}/builder.x") == 0 or die;
          copy( "ximat", "${edgename}/zR${fullrad}/ximat" ) or die "Failed to copy ximat to ${edgename}/\n$!";
          copy( "ximat_small", "${edgename}/zR${fullrad}/ximat_small" ) or die "Failed to copy ximat to ${edgename}/\n$!";

          `echo "$screen_data_files{'grid.nb'}" > ipt`;
          `$ENV{'OCEAN_BIN'}/xipps.x < ipt`;
          move( "ninduced", "nin" ) or die "Failed to move ninduced.\n$!";
  #        `echo $fullrad > ipt`;
  #        `cat ibase epsilon >> ipt`;
  #        `time $ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
          my $reoptName = sprintf( "reopt%2s%04i.R%s", $elname, $elnum, $fullrad );
          unless( -e $reoptName )
          {
            die "Failed to find file: $reoptName\n";
          }
  #        move( "reopt", "rom" ) or die "Failed to move reopt.\n$!";
          `echo 1 3 > ipt`;
          `wc $reoptName >> ipt`;
          `cat $reoptName >> ipt`;
  #        `wc rom >> ipt`;
  #        `cat rom >> ipt`;
          `echo 1 4 >> ipt`;
          `wc nin >> ipt`;
          `cat nin >> ipt`;
          `echo 1 2 >> ipt`;
          `wc zpawinfo/vc_bare${edgename2} >> ipt`;
          `cat zpawinfo/vc_bare${edgename2} >> ipt`;
      
          # Full ximat, but using false/older style in rscombine
          copy( "ipt", "ipt1" ) or die;
          `echo .false. >> ipt1`;
          `echo 0.1 100 >> ipt1`;
          `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXF${fullrad}/ropt`;
          move( "rpot", "$edgename/zRXF$fullrad/" ) or die "Failed to move rpot\n$!";
          move( "rpothires", "$edgename/zRXF$fullrad/" ) or die "Failed to move rpothires\n$!";

          # Full ximat and most up-to-date rscombine settings
          copy( "ipt", "ipt1" ) or die;
          `echo .true. >> ipt1`;
          `wc zpawinfo/vpseud1${edgename2} >> ipt1`;
          `cat zpawinfo/vpseud1${edgename2} >> ipt1`;
          `wc zpawinfo/vvallel${edgename2} >> ipt1`;
          `cat zpawinfo/vvallel${edgename2} >> ipt1`;
          `echo "$screen_data_files{'final.dr'} $final_nr" >> ipt1`;
          `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXT${fullrad}/ropt`;
          move( "rpot", "$edgename/zRXT$fullrad/" ) or die "Failed to move rpot\n$!";
          move( "rpothires", "$edgename/zRXT$fullrad/" ) or die "Failed to move rpothires\n$!";
          copy( $reoptName, "$edgename/zRXT$fullrad/rom" ) or die "Failed to copy rom\n$!";
  #        move( "rom", "$edgename/zRXT$fullrad/" ) or die "Failed to move rom\n$!";
          move( "nin", "$edgename/zRXT$fullrad/" ) or die "Failed to move nin\n$!";

          #####################
          # Now use ximat_small
          move( "ximat_small", "ximat" ) or die "$!";
          `echo "$screen_data_files{'grid.nb'}" > ipt`;
          `$ENV{'OCEAN_BIN'}/xipps.x < ipt`;
          move( "ninduced", "nin" ) or die "Failed to move ninduced.\n$!";
  #        `echo $fullrad > ipt`;
  #        `cat ibase epsilon >> ipt`;
  #        `time $ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
  #        move( "reopt", "rom" ) or die "Failed to move reopt.\n$!";
          `echo 1 3 > ipt`;
          `wc $reoptName >> ipt`;
          `cat $reoptName >> ipt`;
  #        `wc rom >> ipt`;
  #        `cat rom >> ipt`;
          `echo 1 4 >> ipt`;
          `wc nin >> ipt`;
          `cat nin >> ipt`;
          `echo 1 2 >> ipt`;
          `wc zpawinfo/vc_bare${edgename2} >> ipt`;
          `cat zpawinfo/vc_bare${edgename2} >> ipt`;

          copy( "ipt", "ipt1" ) or die;
          `echo .true. >> ipt1`;
          `wc zpawinfo/vpseud1${edgename2} >> ipt1`;
          `cat zpawinfo/vpseud1${edgename2} >> ipt1`;
          `wc zpawinfo/vvallel${edgename2} >> ipt1`;
          `cat zpawinfo/vvallel${edgename2} >> ipt1`;
          `echo "$screen_data_files{'final.dr'} $final_nr" >> ipt1`;
          `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXS${fullrad}/ropt`;
          move( "rpot", "$edgename/zRXS$fullrad/" ) or die "Failed to move rpot\n$!";
          move( "rpothires", "$edgename/zRXS$fullrad/" ) or die "Failed to move rpothires\n$!";
  #        move( "rom", "$edgename/zRXS$fullrad/" ) or die "Failed to move rom\n$!";
          copy( $reoptName, "$edgename/zRXT$fullrad/rom" ) or die "Failed to copy rom\n$!";
          move( "nin", "$edgename/zRXS$fullrad/" ) or die "Failed to move nin\n$!";
          # finished with ximat small
          #######################

      }
    }
    close HFINLIST;
  }
  else  #new run for core
  {
#    copy "specpnt", "specpnt.5";
    my %requiredSpecpnt;

    # Start putting together defaults
    #By default 2 shells, one for OPF and one for the rest
    my $ninter = 2;
    my @angList    = split ' ', $screen_data_files{'grid.ang'};
    $ninter = scalar @angList if( scalar @angList > $ninter );

    my @rmodeList  = split ' ', $screen_data_files{'grid.rmode'};
    $ninter = scalar @rmodeList if( scalar @rmodeList > $ninter );

    my @deltarList = split ' ', $screen_data_files{'grid.deltar'};
    $ninter = scalar @rmodeList if( scalar @deltarList > $ninter );

    my @shellsList = split ' ', $screen_data_files{'grid.shells'};
    $ninter = scalar @shellsList if( scalar @shellsList > $ninter );

    # 
    if( $shellsList[0] <= 0 )
    {
      open IN, "hfinlist" or die;
      <IN> =~ m/^\s*\S+\s+(\d+)/ or die;
      my $zee = $1;
      close IN;
      my $zeeName = sprintf("z%03i",$zee);
      open IN, "zpawinfo/radfile$zeeName" or die;
      <IN> =~ m/^\s+(\S+)/ or die;
      $shellsList[0] = $1;
    }

    unless( $rmodeList[0] =~ m/legendre/ || $rmodeList[0] =~ m/uniform/ )
    {
      $rmodeList[0] = 'legendre';
    }

    $angList[0] = 5 if( $angList[0] < 0 );
    $deltarList[0] = 0.2 if( $deltarList[0] < 0 );

    for( my $i = 1; $i < $ninter; $i++ )
    {
      $angList[$i] = $angList[$i-1] if( scalar @angList <= $i );
      $angList[$i] = 7 if( $angList[$i] < 0 || ! $angList[$i] =~ m/\d/ ); #|| scalar @angList <= $i );
      print "$deltarList[$i] !!  ";
      $deltarList[$i] = $deltarList[$i-1] if( $deltarList[$i] < 0 || scalar @deltarList <= $i );
      print "$deltarList[$i]\n";
      $rmodeList[$i] = 'uniform' unless( $rmodeList[$i] =~ m/legendre/ || $rmodeList[$i] =~ m/uniform/ );

      if( $shellsList[$i] < 0 || $shellsList[$i] <= $shellsList[$i-1] 
                              || $shellsList[$i] >= ( $screen_data_files{'grid.rmax'} - $deltarList[$i]/2 ) )
      {
        print "Ran out of shells, might be making fewer than expected!\n";
        $shellsList[$i] = $screen_data_files{'grid.rmax'};
        $ninter = $i+1;
        last;
      }
    }
    print "$ninter\n\n";

    for( my $i = 0; $i < $ninter; $i++ )
    {
      $requiredSpecpnt{"$angList[$i]"} = 1;
    }
    foreach my $key (sort(keys %requiredSpecpnt)) 
    {
      print "specpnt.$key\n";
      copy( "$ENV{'OCEAN_BIN'}/specpnt.$key", "specpnt.$key" ) or die "Failed to get specpnt.$key\n";
    }

    open MKRB, ">", "mkrb_control" or die "Failed to open mkrb_control for writing\n$!";
#    print MKRB "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $screen_data_files{'grid.ninter'}\n" 
#             . "$screen_data_files{'grid.scheme'} $screen_data_files{'grid.rmode'}\n";
#    close MKRB;
    print MKRB "$screen_data_files{'grid.rmax'} $ninter\n";
    for( my $i = 0; $i < $ninter; $i++ )
    {
      print MKRB "$rmodeList[$i] $shellsList[$i] $deltarList[$i] $angList[$i] specpnt\n";
    }
    close MKRB;
    


    # Identify unique ZNL combos
    my %completeList;
    my %edgelist;
    my %zeelist;
    my @hfinlist;
    open HFINLIST, "hfinlist" or die "Failed to open hfinlist\n";

    while (my $hfinline = <HFINLIST>) 
    {
      ($hfinline =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)/) or die "Malformed hfinlist\t$1 $2 $3 $4 $5 $6\n";
      push @hfinlist, [( $1, $2, $3, $4, $5, $6 ) ];
  #    $ppfilename = $1;
  #    $znucl = $2;
  #    $nnum = $3;
  #    $lnum = $4;
  #    $elname = $5;
  #    $elnum = $6;

      my $siteName = sprintf("z%2s%04i", $5, $6 );
      my $edgeName = sprintf("n%02il%02i", $3, $4 );
  #    my $edgename = sprintf("z%2s%04i_n%02il%02i", $5, $6, $3, $4);
      my $edgeEntry = sprintf "%3i %2i %2i", $2, $3, $4;

      $edgelist{ "$edgeEntry" } = '1';
      $zeelist{ "$2" } = '1';

      push @{ $completeList{ $siteName } }, [( $edgeName, $edgeEntry) ];
      
      unless( -d $siteName )
      {
        mkdir $siteName or die "Failed to make directory $siteName\n$!";
      }
      my $dirName = catfile( $siteName, $edgeName );
      unless( -d $dirName )
      {
        mkdir $dirName or die "Failed to make directory $dirName\n$!";
      }

      foreach my $rad (@rads) 
      {
        my $radName = sprintf("zR%03.2f",$rad);
        mkdir catfile( $siteName, $edgeName, $radName );
      }

    }
    close HFINLIST;

    open EDGELIST, ">", "edgelist" or die "Failed to open edgelist for writing\n$!";
    foreach my $edgeEntry (keys %edgelist )
    {
      print EDGELIST "$edgeEntry\n";
    }
    close EDGELIST;

    open ZEELIST, ">", "zeelist" or die "Failed to open zeelist for writing\n$!";
    print ZEELIST ( scalar keys %zeelist ) . "\n";
    foreach my $Zee (keys %zeelist )
    {
      print ZEELIST "$Zee\n";
    }
    close ZEELIST;

    print "$para_prefix $ENV{'OCEAN_BIN'}/screen_driver.x > screen_driver.log 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_BIN'}/screen_driver.x > screen_driver.log 2>&1") == 0 or die "$!\nFailed to run screen_driver.x\n";
    print "screen_driver.x done\n";

    system("$ENV{'OCEAN_BIN'}/vhommod.x") == 0 or die "$!\nFailed to run vhommod.x\n";

    # Begin file cleaning
    open IN, "<", "sitelist" or die "Failed to open sitelist\n$!";
    <IN>;
    my @siteFiles = ("chi0", "grid", "chi", "avg" );
    my @sitelist;
    while( my $siteline = <IN> )
    {
      $siteline =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/ or die "Malformed sitelist\t$1 $2 $3 $4\n";
      push @sitelist, [( $1, $2 )];
      my $siteName = sprintf("%2s%04i", $1, $3 );
      foreach my $fileName (@siteFiles)
      {
        my $newFileName = catfile( "z$siteName", $fileName );
        move( "${fileName}${siteName}", $newFileName ) or die "${fileName}${siteName}\n$newFileName\n$!";
      }
      
      foreach my $rad (@rads)
      {
        my $radName = sprintf("R%03.2f",$rad);
        my $newFileName = catfile( "z$siteName", "reopt.$radName" );
        move( "reopt${siteName}.$radName", $newFileName );
      }
    }
    close IN;

    my @potFiles = ( "vind", "vind0", "nind", "nind0" );
    for( my $i = 0; $i <= $#hfinlist; $i++ )
    {
      my $siteName = sprintf("z%2s%04i", $hfinlist[$i][4], $hfinlist[$i][5] );
      my $edgeName = sprintf("n%02il%02i", $hfinlist[$i][2], $hfinlist[$i][3] );
      foreach my $rad (@rads)
      {
        my $radName = sprintf("zR%03.2f",$rad);
        foreach my $fileName (@potFiles)
        {
          my $fullFileName = $siteName . "_" . $edgeName . "." . $radName . $fileName;
          my $target = catfile( $siteName, $edgeName, $radName, $fileName );
          print "$fullFileName\t\t$target\n";
          move( $fullFileName, $target ) or die "Failed to move $fullFileName to $target\n$!";
        }
      }
    }
    #Files have been cleaned and put in directories

    # Now start preparing to bring together total potential

    # First load up core-only screened potential files
    #  These are the core-hole potential including the effective screening
    #  of the other core-level electrons. 
    my %vc_bare;
    my %vpseud1;
    my %vvallel;
    # This allows us to loop over these three potentials and read them into their 
    #  own hash locations using the same code
    my %potTypes = ( 'vc_bare' => \%vc_bare, 'vpseud1' => \%vpseud1, 'vvallel' => \%vvallel );
    foreach my $edgeEntry (keys %edgelist )
    {
      print "$edgeEntry\n";
      my @edgeentry = split ' ', $edgeEntry;
      my $edgename2 = sprintf("z%03in%02il%02i",$edgeentry[0], $edgeentry[1], $edgeentry[2]);
      
      foreach my $potType (keys %potTypes )
      {
        my $potfile;
        my @pot;
        my @rad;
  #    $potfile = catfile( "zpawinfo", "vc_bare${edgename2}" );
        $potfile =  catfile( "zpawinfo", "${potType}${edgename2}" );
        open IN, "<", $potfile or die "Failed to open $potfile for reading\n$!";
        while( my $line = <IN> ) 
        {
          $line =~ m/(\d\.\d+[Ee][+-]?\d+)\s+(-?\d\.\d+[Ee][+-]?\d+)/ or die "Failed to parse $potfile\n$line";
          push @rad, $1;
          push @pot, $2;
        }
  #    $vc_bare{ "$edgeEntry" } = [ \@rad, \@pot ];
        ${$potTypes{ $potType }}{ "$edgeEntry" } = [ \@rad, \@pot ];
  #    print $vc_bare{ "$edgeEntry" }[0][0] . "\t" . $vc_bare{ "$edgeEntry" }[1][0] . "\n";
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
      foreach my $rad (@rads)
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
  #      print "Comp list cur site\t" . ${ $completeList{ $currentSite } }[$i] . "\n";
        print "\t" . $currentEdge[0] . "\t" . $currentEdge[1]  . "\n";
        for( my $r = 0; $r < scalar @rads; $r++ )
        {
          my @vindRad; my @vindPot;

          my $radName = sprintf("zR%03.2f",$rads[$r]) ;
          print "\t\t$radName\t$rads[$r]\t$r\n";

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
          chdir $rundir or die;

          # If all-electron augmentation then can only make augmented=true version of screened potential
          #  but shouldn't use the faked atomic calculation version
          my $recon_start = 0; 
          my $recon_stop = 1;
          if( $screen_data_files{ "augment" } =~ m/t/i )
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
              print OUT ".false.\n$screen_data_files{'final.dr'} $final_nr\n";
              close OUT;
              system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt") == 0 or die;
            }
            # No reconstruction
            elsif( $reconstruct == 0 )
            {
              print OUT ".false.\n$screen_data_files{'final.dr'} $final_nr\n";
              close OUT;
              system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt_false") == 0 or die;
              move( "rpot", "rpot_false" ) or die "rpot\n$!";
              move( "rpothires", "rpothires_false" ) or die "rpothires\n$!";
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
              print OUT "$screen_data_files{'final.dr'} $final_nr\n";
              close OUT;
              system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt") == 0 or die;
            }

          }
          
          # back out 3 levels
          $rundir = catfile( updir(), updir(),updir() );
          chdir $rundir;
        }
      }

  #    print $reoptArray[0][0][0] . "\t" . $reoptArray[0][1][0] . "\n";
  #    print $reoptArray[0][0][-1] . "\t" . $reoptArray[0][1][-1] . "\n";
    }

    # then for each site and core level
  }
  # core offsets
  #my $dft = `cat dft`;
  #if( $dft =~ m/qe/i )
  #{
  #  `ln -s ../DFT/Out SCF`;
  #}
  my $core_offset = `cat core_offset`;
  chomp $core_offset;
  if( $core_offset =~ m/false/i )
  {
    print "No core shift\n";
    `rm core_shift.txt` if( -e "core_shift.txt" );
  } else
  {
    `perl $ENV{'OCEAN_BIN'}/core_shift.pl > core_shift.log`;
  }
}

`touch done`;

exit 0;


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
    print "$structRef->{'typat'}[$i]  $structRef->{'elname'}[$structRef->{'typat'}[$i]-1]  $countByName{$structRef->{'elname'}[$structRef->{'typat'}[$i]-1]}\n";
    print $index[$i] . "\n";
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
    push $genRef->{'fulllist'}, $s;
        
  }
  foreach  my $key ( sort { $a <=> $b } keys %sites )
  {
#    push $genRef->{'sitelist'}, sprintf "%s    %s   %s %s %s",  $structRef->{'sites'}[$key-1],
#    push $genRef->{'sitelist'}, sprintf "%s    %s   %s %s %s", $index[$key-1], 
#            $sites{ $key }, $structRef->{'xred'}[$key-1][0],
#            $structRef->{'xred'}[$key-1][1], $structRef->{'xred'}[$key-1][2];
    push $genRef->{'sitelist'}, sprintf "%s %i %i", $structRef->{'elname'}[$structRef->{'typat'}[$key-1]-1],
            $structRef->{'znucl'}[$structRef->{'typat'}[$key]-1], $index[$key-1];
    print "$key\n";
  }

  my %zee;
  foreach  my $key ( sort { $a <=> $b } keys %edges )
  {
    push $genRef->{'edgelist'}, $key;
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

  copy( catfile( updir(), "DFT", "rhoofg" ), "rhoofg" );

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
        unless( $hashRef->{'grid2'}->{'rmode'} =~ m/uniform/ || $hashRef->{'grid2'}->{'rmode'} =~ m/legendre/);

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
  printf OUT "%i\n", $screenRef->{'lmax'};
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

  open OUT, ">", "epsilon" or die $!;
  print OUT ($structureRef->{'epsilon'}) . "\n";
  close OUT;

  # 'screen.quadorder' 'screen.chi0integrand'  'screen.appx'
  
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
            print OUT "$$screenHash->{'final'}->{'dr'} $final_nr\n";
            close OUT;
            system( "$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ropt") == 0 or die;
          }

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

  # Loop over radii and then hfin
  my $offset;
  for( my $i = 0; $i < scalar @{$screenRef->{'shells'}}; $i++ )
  {
    my $rad_dir = sprintf("zR%03.2f", $screenRef->{'shells'}[$i] );

    printf "\nRadius = %03.2f Bohr\n", $screenRef->{'shells'}[$i];

    # If we are averaging, new shift by radius
    if( $screenRef->{'core_offset'}->{'average'} )
    {
      $offset = -( $Vsum + $newWsum[$i] ) * $Ry2eV / ( scalar @hfin );
  #    print "$rad_dir\t$offset\n";
      print "  core_offset was set to true. Now set to $offset  \n";
    } else
    {
      $offset = $screenRef->{'core_offset'}->{'energy'};
    }

    print  "Site index    New potential   new1/2 Screening   core_offset       total offset\n";
    print  "                  (eV)             (eV)              (eV)              (eV)\n";
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
      printf "   %7i   %16.9f  %15.9f  %15.9f  %16.7f\n", $el_rank, $newPot[$j]*$Ry2eV, $newWshift[$j][$i]*$Ry2eV, $offset, $shift;

      my $string = sprintf("z%s%04d/n%02dl%02d",$el, $el_rank,$nn,$ll);
      open OUT, ">$string/$rad_dir/cls" or die "Failed to open $string/$rad_dir/cls\n$!";
      print OUT $shift . "\n";
      close OUT;
    }

    print "\n";

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
    foreach my $key (keys $newRef )
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

