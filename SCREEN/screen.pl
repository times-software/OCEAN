#!/usr/bin/perl
# Copyright (C) 2010, 2013 - 2019 OCEAN collaboration
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


my @CommonFiles = ("znucl", "opf.hfkgrid", "opf.fill", "opf.opts", "pplist", "screen.shells", 
                   "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution", "epsilon", 
                   "screen.k0", "scfac", "core_offset", "dft", "avecsinbohr.ipt", 
                   "nspin", "prefix" );

my @CommonFiles2 = ("para_prefix", "calc" );

my @ScreenFiles = ("screen.grid.scheme", "screen.grid.rmode", "screen.grid.ninter", 
                   "screen.grid.shells", "screen.grid.xyz", "screen.grid.rmax", "screen.grid.ang",
                   "screen.lmax", "screen.grid.nb", "screen.grid.nr", "screen.final.rmax", 
                   "screen.final.dr", "screen.legacy", "screen.model.dq", "screen.model.qmax", 
                   "screen.augment", "screen.wvfn", "screen.convertstyle", "screen.inversionstyle", 
                   "screen.mode" );

my @DenDipFiles = ("rhoofg", "bvecs", "efermiinrydberg.ipt", "xmesh.ipt");
my @DenDipFiles2 = ( "masterwfile", "listwfile", "enkfile", "kmesh.ipt", "brange.ipt" );

my @ExtraFiles = ("specpnt", "Pquadrature", "hqp", "lqp", "gauss16", "EvenQuadHalf.txt" );

my @DFTFiles = ( "potofr" );

my $runSCREEN = 1;
if (-e "../PREP/PAW/old" && -e "done" ) {
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
if( $screen_data_files{ 'grid.ang' } =~ m/(\w+)\s+(\d+)/ )
{
  my $angularGridFile = $1 . '.' . $2;
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
  copy( "../PREP/PAW/$_", $_ ) or die "Failed to get $_ from PREP/PAW/\n$!";
}

foreach (@ExtraFiles) {
  copy( "$ENV{'OCEAN_BIN'}/$_", $_ ) or die "Failed to get $_ from $ENV{'OCEAN_BIN'}/\n$!";
}

my $dft = `cat dft`;
if( $dft =~ m/qe/i )
{
  `ln -s ../DFT/Out SCF`;
  `ln -s ../DFT/SCREEN/Out .`;
}

open WVFN, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
print WVFN $screen_data_files{ 'wvfn' } . "\n";
close WVFN;

open LISTW, "listwfile" or die "Failed to open listwfile\n";
while (<LISTW> ) {
  m/(\d+)\s+(\S+)/ or die "Malformed listwfile\n";
  `ln -sf ../PREP/PAW/$2 $2`;
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
print "Running PAW Setup\n";
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
  copy "specpnt", "specpnt.5";
  open MKRB, ">", "mkrb_control" or die "Failed to open mkrb_control for writing\n$!";
  print MKRB "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $screen_data_files{'grid.ninter'}\n"
           . "$screen_data_files{'grid.scheme'} $screen_data_files{'grid.rmode'}\n";
  close MKRB;

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
  

  open ZEELIST, ">", "zeelist" or die "Failed to open zeelist for writing\n$!";
  print ZEELIST "0\n";
  close ZEELIST;

  system("$para_prefix $ENV{'OCEAN_BIN'}/screen_driver.x") == 0 or die "$!\nFailed to run screen_driver.x\n";
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
          `time $ENV{'OCEAN_BIN'}/xipps.x < ipt`;
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
          `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXF${fullrad}/ropt`;
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
          `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXT${fullrad}/ropt`;
          move( "rpot", "$edgename/zRXT$fullrad/" ) or die "Failed to move rpot\n$!";
          move( "rpothires", "$edgename/zRXT$fullrad/" ) or die "Failed to move rpothires\n$!";
          copy( $reoptName, "$edgename/zRXT$fullrad/rom" ) or die "Failed to copy rom\n$!";
  #        move( "rom", "$edgename/zRXT$fullrad/" ) or die "Failed to move rom\n$!";
          move( "nin", "$edgename/zRXT$fullrad/" ) or die "Failed to move nin\n$!";

          #####################
          # Now use ximat_small
          move( "ximat_small", "ximat" ) or die "$!";
          `echo "$screen_data_files{'grid.nb'}" > ipt`;
          `time $ENV{'OCEAN_BIN'}/xipps.x < ipt`;
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
          `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXS${fullrad}/ropt`;
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
    copy "specpnt", "specpnt.5";
    open MKRB, ">", "mkrb_control" or die "Failed to open mkrb_control for writing\n$!";
    print MKRB "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $screen_data_files{'grid.ninter'}\n" 
             . "$screen_data_files{'grid.scheme'} $screen_data_files{'grid.rmode'}\n";
    close MKRB;
  #  `echo "$screen_data_files{'grid.rmax'} $screen_data_files{'grid.nr'} $screen_data_files{'grid.ninter'}" > mkrb_control`;
  #  `echo "$screen_data_files{'grid.scheme'} $screen_data_files{'grid.rmode'}" >> mkrb_control`;


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

    print "$para_prefix $ENV{'OCEAN_BIN'}/screen_driver.x\n";
    system("$para_prefix $ENV{'OCEAN_BIN'}/screen_driver.x") == 0 or die "$!\nFailed to run screen_driver.x\n";
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
    `time perl $ENV{'OCEAN_BIN'}/core_shift.pl > core_shift.log`;
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
