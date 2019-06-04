#!/usr/bin/perl
# Copyright (C) 2014, 2016, 2017, 2019 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Path qw( rmtree );
use Cwd 'abs_path';
use File::Compare;
use File::Copy;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/dendip\.pl/;
#  $ENV{"OCEAN_BIN"} = $1;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################
my $RunDenDip = 0;

#open STATUS, "../DFT/abinitstage.stat" or die;
#<STATUS> =~ m/(\d)/;
#if ($1 == 1) {
#  print "Wavefunctions unchanged. Nothing to be done for DenDip\n";
#  exit 0;
#}

my $split_dft = 0;
my $stat = 0;
$stat = 1 if (-e "done");
`rm -f done`;
my $oldden = 0;
$oldden = 1 if (-e "../DFT/old");


my @AbFiles = ( "rhoofr", "density.out", "nkpt", "screen.nkpt", "qinunitsofbvectors.ipt", "efermiinrydberg.ipt");
my @CommonFiles = ( "avecsinbohr.ipt", "nspin", "xmesh.ipt", "dft", "nspin", "dft.split", 
                    "screen.mode", "bse.wvfn", "k0.ipt", "calc", "screen.wvfn" );
my @NewMethodFiles = ( "ntype", "typat", "natoms", "znucl", "taulist", "edges", "core_offset", "metal", "cksshift",
                       "cksstretch", "nedges", "edges", "pplist", "opf.opts", "opf.fill" );
my @ExtraFiles = ("specpnt.5", "Pquadrature", "sphpts" );
my @OPFFiles = ("hfinlist", "xyz.wyck");

foreach (@AbFiles) {
  system("cp ../DFT/$_ .") == 0 or die "Failed to copy $_\n";
}
foreach (@CommonFiles)
{
  system("cp ../Common/$_ .") == 0 or die "Failed to copy $_\n";
}
foreach (@NewMethodFiles) {
  system("cp ../Common/$_ .") == 0 or die "Failed to copy $_\n";
}
foreach (@OPFFiles)
{
  system("cp ../OPF/$_ .") == 0 or die "Failed to copy $_\n";
}




print "$stat  $oldden\n";
unless ($stat && $oldden) {
  -e "../DFT/SCx_DEN" or die "SCx_DEN not found\n";
  `ln -sf ../DFT/SCx_DEN SCx_DEN`;

  open IN, "dft.split" or die "Failed to open dft.split\n$!";
  if( <IN> =~ m/t/i )
  {
    $split_dft = 1;
    print "DFT has been split into to runs\n";
  }

  #`ln -sf ../DFT/RUN????_WFK .`;

  `tail -n 1 rhoofr > nfft`;

  system("$ENV{'OCEAN_BIN'}/nelectron.x") == 0 
    or die "Failed to run nelectron.x\n";

  system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0
    or die "Failed to run bvecs.pl\n";

  system("$ENV{'OCEAN_BIN'}/gvecs2.pl") == 0
    or die "Failed to run gvecs2.pl\n";
}

open IN, "calc" or die "Failed to open calc\n";
<IN> =~ m/(\w+)/ or die "Failed to parse calc\n";
my $calc = $1;
close IN;

open IN, "screen.mode" or die "Failed to open screen.mode\n";
<IN> =~ m/(\w+)/ or die "Failed to parse screen.mode\n";
my $screen_mode = $1;
close IN;
my $run_screen = 1;
if( $calc =~ m/val/i )
{
  $run_screen = 0 unless( $screen_mode =~ m/grid/i );
}

my $rundir;
my @nkpt;
## process screen wf files ##
if( $run_screen == 1 )
{

  open NKPT, "screen.nkpt" or die "Failed to open screen.nkpt";
  <NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse screen.nkpt\n";
  @nkpt = ($1, $2, $3);
  close NKPT;
  $rundir = "../DFT/SCREEN";

  unless( -e "PAW/done" && -e "${rundir}/old" ) {
  `rm -r PAW` if (-e "PAW");
  mkdir "PAW"; 
  chdir "PAW";

  open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
  print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
  close NKPT;


  foreach ("kmesh.ipt", "brange.ipt", "qinunitsofbvectors.ipt" ) {
    system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
  }
  #`cp ../qinunitsofbvectors.ipt .`;
  `cp ../bvecs .`;
  `cp ../dft .`;
  `cp ../nspin .`;
  `cp ../${rundir}/umklapp .`;

  my $runfile;
  $runfile = sprintf("../${rundir}/RUN%04u_WFK", 1 );
  system("ln -s $runfile .") == 0  or die "Failed to link $runfile\n";

  system("$ENV{'OCEAN_BIN'}/wfconvert.x") == 0 
    or die "Failed to run wfconvert.x\n";

  `touch done`;
  chdir "../";
  }
  else {
    `touch PAW/old`;
    print  "Nothing needed for SCREEN wfns\n";
  }
  print "Done with SCREEN files\n";
}
else
{
  print "Nothing needed for SCREEN wvfn\n";
}

## process bse wf files ##

open NKPT, "nkpt" or die "Failed to open nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse nkpt\n";
@nkpt = ($1, $2, $3);
close NKPT;

$rundir = sprintf("../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);
my @BSECommonFiles = ( "qinunitsofbvectors.ipt", "bvecs", "dft", "nelectron", "avecsinbohr.ipt",
                       "nspin", "dft.split", "prefix", "natoms", "typat", "ntype","znucl", "taulist",
                       "edges", "k0.ipt", "core_offset", "metal", "cksshift", "cksstretch" );
my @rundirFiles = ( "kmesh.ipt", "brange.ipt", "umklapp" );
my @BSEBonusFiles = ("xmesh.ipt", "calc", "hfinlist", "xyz.wyck" );

my $runBSE = 1;
if( -e "BSE/done" && -e "${rundir}/old" )
{
  $runBSE = 0;
  foreach( @BSECommonFiles )
  {
    if( compare( "$_", "BSE/$_" ) != 0 )
    {
      print "Differences found in $_\nWill re-run BSE prep\n";
      last;
    }
  }
  if( $runBSE == 0 )
  {
    foreach( @rundirFiles )
    {
      if( compare( "$rundir/$_", "BSE/$_" ) != 0 )
      {
        print "Differences found in $_\nWill re-run BSE prep\n";
        last;
      }
    }
  }
  if( $runBSE == 0 )
  {   
    $runBSE = 2 if( compare( "xmesh.ipt", "BSE/xmesh.ipt" ) != 0 );
    if( compare( "hfinlist", "BSE/hfinlist" ) != 0 || compare( "calc", "BSE/calc" ) != 0 )
    {
      $runBSE = 3;
    }
  } 
} 

  
if( $runBSE != 0 )
{
  unlink "BSE/done";
  if( $runBSE == 1 )
  {
    rmtree( "BSE" ) if (-e "BSE");
    mkdir "BSE";
  }
  chdir "BSE";

  open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
  print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
  close NKPT;

  foreach( @rundirFiles )
  { 
    copy( "../${rundir}/$_", $_ );
  }
  foreach( @BSECommonFiles )
  {
    copy ( "../$_", $_ );
  }
  foreach( @BSEBonusFiles )
  {
    copy ( "../$_", $_ );
  }

  foreach( @ExtraFiles )
  {
    copy("$ENV{'OCEAN_BIN'}/$_", $_ ) or die;
  }

  if( $runBSE == 1 )
  {

    my $runfile;
    my $Nfiles = 1;

    $Nfiles = 2 if( $split_dft == 1);

    for (my $i = 1; $i <= $Nfiles; $i++) {
      $runfile = sprintf("../${rundir}/RUN%04u_WFK", $i );
      system("ln -s $runfile .") == 0 or die "Failed to link $runfile\n";
    }

    system("$ENV{'OCEAN_BIN'}/wfconvert.x") == 0 
      or die "Failed to run wfconvert.x\n";

    system("$ENV{'OCEAN_BIN'}/ofermi.pl") == 0
      or die "Failed to run ofermi.pl\n";
    `cp eshift.ipt ../`;
    system("cp ../efermiinrydberg.ipt ./") == 0 
      or die "Failed to copy efermiinrydberg.ipt\n";
  }
  if( $runBSE == 1 || $runBSE == 3 )
  {
    copy( "kmesh.ipt", "kgrid" ) or die "$!";
    copy( "k0.ipt", "scaledkzero.ipt" ) or die "$!";
    copy( "qinunitsofbvectors.ipt", "cksdq" ) or die "$!";
    my $is_xas;
    open TMPFILE, "calc" or die "Failed to open calc\n";
    my $mode = <TMPFILE>;
    close TMPFILE;
    chomp($mode);
    if( lc($mode) =~ m/xes/ )
    {
      print "Calculating XES\n";
      $is_xas = 0;
    }
    elsif( lc($mode) =~ m/xas/ )
    {
      print "Calculating XAS\n";
      $is_xas = 1;
    }
    else
    {
      print "Unrecognized mode. Calculating XAS\n";
      $is_xas = 1;
    }

    open TMPFILE, ">cks.normal" or die "Failed to open cks.normal for writing\n$!";
    if( $is_xas == 1 )
    {
      print TMPFILE ".true.\n";
    }
    else
    {
      print TMPFILE ".false.\n";
    }
    close TMPFILE;

#      open CKS, ">cks.in" or die "Failed to open cks.in\n";
    my $znl_string = 0;
    my $ncks = 0;
    my $cks_string;

    my %unique_z;
    open EDGE, "hfinlist" or die "Failed to open hfinlist\n$!";
    while (<EDGE>)
    {
      $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
      my $ppname = $1;
      my $znum = $2;
      my $nnum = $3;
      my $lnum = $4;
      my $elname = $5;
      my $elnum = $6;

      my $cks;
      if( $is_xas == 1  ) {
        $cks = "cksc.${elname}";
      }
      else {
        $cks = "cksv.${elname}"
      }

      # For each unique Z we need to grab some files from OPF
      unless( exists $unique_z{ "$znum" } )
      {
        my $zstring = sprintf("z%03in%02il%02i", $znum, $nnum, $lnum);
        print $zstring ."\n";
        
        `ln -sf ../../OPF/zpawinfo/gk*${zstring} .`;
        `ln -sf ../../OPF/zpawinfo/fk*${zstring} .`;
        `ln -sf ../../OPF/zpawinfo/melfile${zstring} .`;
        `ln -sf ../../OPF/zpawinfo/coreorb${zstring} .`;
        
        
        my $z3 = sprintf("z%03i", $znum);
        `ln -sf ../../OPF/zpawinfo/phrc?${z3} .`;
        `ln -sf ../../OPF/zpawinfo/prjfile${z3} .`;
        `ln -sf ../../OPF/zpawinfo/ft?${z3} .`;
        `ln -sf ../../OPF/zpawinfo/ae?${z3} .`;
        `ln -sf ../../OPF/zpawinfo/ps?${z3} .`;
        `ln -sf ../../OPF/zpawinfo/corezeta${z3} .`;
        `ln -sf ../../OPF/zpawinfo/radfile${z3} .`;
        
      } 
      
      print "CKS NAME = $cks\n";
      my $temp_znl = sprintf "%i  %i  %i", $znum, $nnum, $lnum;
      if( $znl_string == 0 ) 
      { 
        $znl_string = $temp_znl;
        open ZNL, ">ZNL" or die;
        print ZNL "$znl_string\n";
        close ZNL;
      } 
      
      if( $znl_string eq $temp_znl )
      { 
        $ncks++;
        $cks_string .= "$elname  $elnum  $cks\n";
      } 
      else
      { 
        print "New ZNL!\nRunning $ncks through cks\n";
        unless ( $ncks == 0 )
        { 
          open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
          print CKSIN "$ncks\n$cks_string";
          close CKSIN;
          print "cks\n";
          system("$ENV{'OCEAN_BIN'}/cks.x < cks.in > cks.log") == 0 or die;
        } 
        $znl_string = $temp_znl;
        open ZNL, ">ZNL" or die;
        print ZNL "$znl_string\n";
        close ZNL;
        $ncks = 1;
        $cks_string = "$elname  $elnum  $cks\n";
      } 
      
    } 
    close EDGE;

    unless ( $ncks == 0 )
    {
      print "Final cks: $ncks\n";
      open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
      print CKSIN "$ncks\n$cks_string";
      close CKSIN;
      print "cks\n";
      system("$ENV{'OCEAN_BIN'}/cks.x < cks.in > cks.log") == 0 or die;
    }

  }

  # 3 is cks-only option, 1 & 2 will make new u2
  unless( $runBSE == 3 )
  {
    print "Running setup\n";
    system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0
      or die "Failed to run setup\n";

    print "conugtoux\n";
    system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
    print "orthog\n";
    system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;
  }

  `touch done`;
  chdir "../";
}
else {
  `touch BSE/old`;
  print "Nothing needed for bse wfns\n";
}

unless ($stat && $oldden) {
system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0
  or die "Failed to run rhoofg.x\n";
`wc -l rhoG2 > rhoofg`;
`sort -n -k 6 rhoG2 >> rhoofg`;
}

`touch done`;
exit 0;

