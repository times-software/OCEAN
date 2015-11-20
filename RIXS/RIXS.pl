#!/usr/bin/perl

# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Copy;

if (! $ENV{"AI2NBSE_BIN"} ) {
  $0 =~ m/(.*)\/RIXS\.pl/;
  $ENV{"AI2NBSE_BIN"} = $1;
  print "AI2NBSE_BIN not set. Setting it to $1\n";
}
if (! $ENV{"AI2NBSE_WORKDIR"}){ $ENV{"AI2NBSE_WORKDIR"} = `pwd` . "../" ; }

my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );


my @CommonFiles = ("epsilon", "xmesh.ipt", "nedges", "k0.ipt", "nbuse.ipt", 
  "cnbse.rad", "cnbse.ways", "metal", "cksshift", "cksstretch", "cksdq", 
  "cnbse.niter", "cnbse.spect_range", "cnbse.broaden", "cnbse.mode", "nphoton", "dft", 
  "para_prefix", "cnbse.strength", "serbse", "core_offset", "avecsinbohr.ipt" );

my @DFTFiles = ("nelectron");

my @DenDipFiles = ("kmesh.ipt", "masterwfile", "listwfile", "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "tmels", "nelectron", "eshift.ipt", "ldagap" );

my @WFNFiles = ("kmesh.ipt",  "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", "nbuse.ipt", "wvfcninfo", "wvfvainfo", "nbuse_xes.ipt", "obf_control", "ibeg.h");

my @ExtraFiles = ("Pquadrature", "sphpts" );

my @PawFiles = ("hfinlist" , "xyz.wyck");

# At the moment this is minimum viable edits
#  since you have to run OCEAN anyway this is just a hack job 
#  of the cnbse_mpi.pl script
my @rixsFiles = ("gwgap", "gwcstr", "gwvstr", "bflag", "lflag", "niter", "backf", "aldaf", "qpflg", "bwflg", "bande", "abs_gap", "spect.h");

my @prepFiles = ("gvectors", "rhoofr", "nfft" );

my @cnbseFiles = ("wvfcninfo", "wvfvainfo", "ZNL" );

foreach (@CommonFiles) {
  copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}

open DFT, "dft" or die "Failed to open dft\n";
<DFT> =~ m/(\w+)/ or die;
my $dft_type = $1;
close DTF;
my $obf;
if( $dft_type =~ m/obf/i )
{
  $obf = 1;
}
else
{
  $obf = 0;
}


if( $obf == 1 ) 
{
  foreach (@DFTFiles) {
    copy( "../DFT/$_", $_) or die "Failed to get DFT/$_\n$!";
  }
  foreach (@WFNFiles) {
    copy( "../zWFN/$_", $_ ) or die "Failed to get zWFN/$_\n$!";
  }
}
else
{
  foreach (@DenDipFiles) {
    copy( "../PREP/BSE/$_", $_ ) or die "Failed to get PREP/BSE/$_\n$!" ;
  }
  foreach (@prepFiles) {
    copy( "../PREP/$_", $_ ) or die "Failed to get PREP/BSE/$_\n$!" ;
  }
}

foreach (@ExtraFiles) {
  copy( "$ENV{'AI2NBSE_BIN'}/$_", $_ ) or die "Failed to get ../$_\n$!";
}

foreach (@PawFiles) {
  copy( "../SCREEN/$_", $_ ) or die "Failed to get ../SCREEN/$_\n$!";
}

foreach (@rixsFiles) {
  copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}

foreach (@cnbseFiles) {
  copy( "../CNBSE/$_", $_ ) or die "Failed to get CNBSE/$_\n$!";
}

copy( "$ENV{'AI2NBSE_BIN'}/stepper.h", "stepper.h");
copy( "$ENV{'AI2NBSE_BIN'}/rixs.h", "rixs.h");

##### Trigger serial bse fallback here
# later we should remove this and fold these two perl scripts together
my $run_serial = 0;
if( -e "serbse" )
{
  open IN, "serbse" or die "$!";
  <IN> =~ m/(\d)/;
  if( $1 == 1 )
  {
    print "Serial BSE requested!!\n";
#    exit system("$ENV{'AI2NBSE_BIN'}/OBF_cnbse.pl");
    $run_serial = 1;
  }
  close IN;
}
unless( -e "$ENV{'AI2NBSE_BIN'}/ocean.x" )
{
  print "Parallel BSE executable not present in $ENV{'AI2NBSE_BIN'}\nAttempting serial run ...\n";
#  exit system("$ENV{'AI2NBSE_BIN'}/OBF_cnbse.pl");
  $run_serial = 1;
}

# Grab the needed photon files, copy them into the CNBSE directory,
# and store their names into the array @photon_files
my $nphoton = -1;
open NPHOTON, "nphoton" or die "Failed to open nphoton\n$!";
while (<NPHOTON>)
{
  if( $_ =~ m/(-?\d+)/ )
  {
    $nphoton = $1;
    last;
  }
}
close NPHOTON;

my @photon_files;
if( $nphoton > 0 )
{
  for( my $i = 1; $i <= $nphoton; $i++ )
  {
    if( -e "../photon${i}" )
    {
      copy( "../photon${i}", "photon${i}") or die "Failed to copy photon$i\n$!";
      push @photon_files, "photon${i}";
    }
    elsif( -e "../jtv${i}" )
    {
       copy( "../jtv${i}", "jtv{$i}") or die "Failed to copy jtv$1\n$!";
      push @photon_files, "jtv${i}";
    }
    else
    {
      print "Could not find photon file # $i\n";
    }
  }
}
else
{
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
  my @sorted_photon_files = sort{ my ($anum,$bnum); $a =~ m/\w+(\d+)/; $anum = $1; $b =~ m/\w+(\d+)/; $bnum = $1; $anum <=> $bnum } @photon_files;
  @photon_files = @sorted_photon_files;

  foreach( @photon_files )
  {
    print "        $_\n";
    copy( "../$_", $_ ) or die "$!";
  }
}



##### misc other setup
#`echo gmanner > format65`;
#`cp nspin nspn.ipt`;
copy( "kmesh.ipt", "kgrid" ) or die "$!";
copy( "k0.ipt", "scaledkzero.ipt" ) or die "$!";
copy( "qinunitsofbvectors.ipt", "cksdq" ) or die "$!";

my $para_prefix = "";
if( open PARA_PREFIX, "para_prefix" )
{
  $para_prefix = <PARA_PREFIX>;
  chomp($para_prefix);
  close( PARA_PREFIX);
} else
{
  print "Failed to open para_prefix. Error: $!\nRunning serially\n";
}


# Set up mode
  my $is_xas;
  open TMPFILE, "cnbse.niter" or die "Failed to open cnbse.niter\n$!";
  <TMPFILE> =~ m/(\d+)/ or die "Failed to parse cnbse.niter";
  my $num_haydock_iterations = $1;
  close TMPFILE;
  
  open TMPFILE, "cnbse.strength" or die "Failed to open cnbse.strength\n$!";
  <TMPFILE> =~ m/([0-9]*\.?[0-9]+)/ or die "Failed to parse cnbse.strength\n";
  my $interaction_strength = $1; 
  close TMPFILE;
    
  print "Unrecognized mode. Calculating RIXS\n";

# write cks.normal file
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

#write mode file
  open TMPFILE, ">mode" or die "Failed to open mode for writing\n$!";
  print TMPFILE "$interaction_strength    $num_haydock_iterations\n";
  close TMPFILE;

###############
# If we are using QE/ABI w/o OBFs we need to set nbuse
my $run_text = '';
if( $obf == 1 )
{
  close RUNTYPE;
  if ($is_xas == 1 )
  {
    $run_text = 'XAS';
  } 
  else
  {
    move( "nbuse.ipt", "nbuse_xas.ipt" ) or die "$!";
    copy( "nbuse_xes.ipt", "nbuse.ipt" ) or die "$!";
#    `mv nbuse.ipt nbuse_xas.ipt`;
#    `cp nbuse_xes.ipt nbuse.ipt`;
    $run_text = 'XES';
    print "XES!\n";
  }
}
else  ### Abi/QE w/o obf
{ 
  open NBUSE, "nbuse.ipt" or die "Failed to open nbuse.ipt\n";
  <NBUSE> =~ m/(\d+)/ or die "Failed to parse nbuse.ipt\n";
  my $nbuse = $1;
  close NBUSE;
  my @brange;
  if ($nbuse == 0) {
    open BRANGE, "brange.ipt" or die "Failed to open brange.ipt\n";
    <BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Failed to parse brange.ipt\n";
    $brange[0] = $1;
    $brange[1] = $2;
    <BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Failed to parse brange.ipt\n";
    $brange[2] = $1;
    $brange[3] = $2;
    close BRANGE;

    if( $is_xas == 1 )
    {
      $run_text = 'XAS';
      $nbuse = $brange[3] - $brange[2] + 1;
    }
    else
    {
      print "XES!\n";
      $run_text = 'XES';
      $nbuse = $brange[1] - $brange[0] + 1;
    }
    open NBUSE, ">nbuse.ipt" or die "Failed to open nbuse.ipt\n";
    print NBUSE "$nbuse\n";
    close NBUSE;
  }
}


system("$ENV{'AI2NBSE_BIN'}/getnval.x") == 0 or die "Failed to get nval\n";


#####################
if( $obf == 1 )
{
  if( -e "../zWFN/u2par.dat" )
  {
    `ln -s ../zWFN/u2par.dat`;
    open OUT, ">bloch_type" or die;
    print OUT "new\n";
    close OUT;
  }
  else
  {
    `ln -s ../zWFN/u2.dat`;
  }
  die "OBFs not supported yet\n";
}
else  # We are using abi/qe path w/o obfs
{
  # grab .Psi
  `touch .Psi`;
  system("rm .Psi*");
  open LISTW, "listwfile" or die "Failed to open listwfile\n";
  while (<LISTW>) 
  {
    $_ =~ m/(\d+)\s+(\S+)/ or die "Failed to parse listwfile\n";
    system("ln -sf ../PREP/BSE/$2 .") == 0 or die "Failed to link $2\n";
  }  

  print "Running setup\n";
  system("$ENV{'AI2NBSE_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";

  if (-e "../PREP/BSE/u2.dat")
  {
    `ln -s ../PREP/BSE/u2.dat`;
  }
  else
  {
    print "conugtoux\n";
    system("$ENV{'AI2NBSE_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
    print "orthog\n";
    system("$ENV{'AI2NBSE_BIN'}/orthog.x > orthog.log") == 0 or die;
  }
}

#write decut
open OUT, ">decut" or die;
print OUT "10 2\n";
close OUT;

# populate gwipt
open IN, "gwcstr" or die;
my $gwcstr = <IN>;
chomp($gwcstr);
close IN;

open IN, "gwvstr" or die;
my $gwvstr = <IN>;
chomp $gwvstr;
close IN;

open IN, "abs_gap" or die;
my $absgap = <IN>;
chomp $absgap;
close IN;

open IN, "gwgap" or die;
my $gwgap = <IN>;
chomp $gwgap;
close IN;

open IN, "ldagap" or die;
my $ldagap = <IN>;
chomp $ldagap;
close IN;

unless( $absgap)
{
  $gwgap = $gwgap + $ldagap;
}
open OUT, ">gwipt" or die;
print OUT "$gwgap $ldagap $gwcstr $gwvstr\n";
close OUT;

# build niter.h
my $bflag = `cat bflag`;
my $lflag = `cat lflag`;
my $niter = `cat niter`;
my $backf = `cat backf`;
my $aldaf = `cat aldaf`;
my $qpflg = `cat qpflg`;
my $bwflg = `cat bwflg`;
my $bande = `cat bande`;
open OUT, ">niter.h" or die;
printf OUT "  %5i          niter\n", $niter;
printf OUT "  %5i          bflag\n", $bflag;
printf OUT "  %5i          lflag\n", $lflag;
printf OUT "  %5i          backf\n", $backf;
printf OUT "  %5i          aldaf\n", $aldaf;
printf OUT "  %5i          qpflg\n", $qpflg; 
printf OUT "  %5i          bwflg\n", $bwflg; 
printf OUT "  %5i          bande\n", $bande; 

`$ENV{'AI2NBSE_BIN'}/bridgegw.x > bridgegw.log`;
`$ENV{'AI2NBSE_BIN'}/condens2.x             > condens.log`;
`$ENV{'AI2NBSE_BIN'}/ll.x        < epsilon  > ll.log`;
`$ENV{'AI2NBSE_BIN'}/whom0.x                > whom0.log`;
`$ENV{'AI2NBSE_BIN'}/bridgelad.x            > bridgelad.log`;


# 
open IN, "spect.h" or die;
my $spect = <IN>;
chomp $spect;
close IN;

open IN, "cnbse.broaden" or die "$!";
my $broaden = <IN>;
chomp $broaden;
close IN;

open SPECT, ">spect.in" or die "Failed to open spect.in";
print SPECT "$broaden $spect -1\n";
close SPECT;

#
open EDGE, "hfinlist" or die "$!";
while (<EDGE>)
{
  $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
  my $ppname = $1;
  my $znum = $2;
  my $nnum = $3;
  my $lnum = $4;
  my $elname = $5;
  my $elnum = $6;

  my $echamp = sprintf("../CNBSE/echamp_$elname.%04i_${nnum}$alphal{$lnum}_01.0001", $elnum );
  copy( $echamp, "echamp.$elnum" ) or die "$echamp\n$!";

  my $mels = sprintf("../CNBSE/mels.z%03un%02ul%02up%02u", $znum, $nnum, $lnum, 2 );
  copy( $mels, "mels" ) or die "$mels\n$!";

  my $cks = sprintf("../CNBSE/cksv.${elname}%04u", $elnum );
  copy( $cks, "cbinf" ) or die "$cks\n$!";

  system("echo cbinf | $ENV{'AI2NBSE_BIN'}/dotter.x") == 0 or die "Failed to run dotter\n";

  for( my $i = 0; $i <= $lnum; $i++ )
  {
    my $beff = sprintf("beff%02i", $i+1 );
    copy( $beff, "$beff.$elnum" ) or die "$beff\n$!";
  }
}



`$ENV{'AI2NBSE_BIN'}/smnanl.x    < rixs.h > smnanl.log`;

`$ENV{'AI2NBSE_BIN'}/spect.x < spect.in`;
