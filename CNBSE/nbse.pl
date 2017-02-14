#!/usr/bin/perl

# Copyright (C) 2015 - 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
# Adapted from cnbse_mpi.pl

use strict;
use File::Copy;
use POSIX;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/nbse\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }

my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

my @CommonFiles = ("epsilon", "xmesh.ipt", "nedges", "k0.ipt", "nbuse.ipt", 
  "cnbse.rad", "cnbse.ways", "metal", "cksshift", "cksstretch", "cksdq", 
  "cnbse.niter", "cnbse.spect_range", "cnbse.broaden", "cnbse.mode", "nphoton", "dft", 
  "para_prefix", "cnbse.strength", "serbse", "core_offset", "avecsinbohr.ipt", 
  "cnbse.solver", "cnbse.gmres.elist", "cnbse.gmres.erange", "cnbse.gmres.nloop", 
  "cnbse.gmres.gprc", "cnbse.gmres.ffff", "cnbse.write_rhs", "spin_orbit", "nspin", 
  "niter", "backf", "aldaf", "bwflg", "bande", "bflag", "lflag", "decut", "spect.h" );

my @DFTFiles = ("nelectron", "rhoofr");

my @DenDipFiles = ("kmesh.ipt", "masterwfile", "listwfile", "efermiinrydberg.ipt", 
                   "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "tmels", "nelectron", 
                   "eshift.ipt" );

my @WFNFiles = ("kmesh.ipt",  "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", 
                "wvfcninfo", "wvfvainfo", "obf_control", "ibeg.h", "q.out", "tmels.info", 
                "val_energies.dat", "con_energies.dat" );

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
  open OUT, ">tmel_selector" or die;
  print OUT "1\n";
  close OUT;
  open OUT, ">enk_selector" or die;
  print OUT "1\n";
  close OUT;
  open OUT, ">bloch_selector" or die;
  print OUT "1\n";
  close OUT;
}
else
{
  foreach (@DFTFiles) {
    copy( "../PREP/$_", $_) or die "Failed to get DFT/$_\n$!";
  }

  foreach (@DenDipFiles) {
    copy( "../PREP/BSE/$_", $_ ) or die "Failed to get PREP/BSE/$_\n$!" ;
  }
  open OUT, ">tmel_selector" or die;
  print OUT "0\n";
  close OUT;
  open OUT, ">enk_selector" or die;
  print OUT "0\n";
  close OUT;
  open OUT, ">bloch_selector" or die;
  print OUT "0\n";
  close OUT;

}



##### Rixs requires that gmres be used #####
open IN, "cnbse.solver" or die "Failed to open cnbse.solver!\n$!";
my $line = <IN>;
close IN;
my $solver;
if( lc($line) =~ m/hay/ )
{
  $solver = 'hay';
}
elsif( lc($line) =~ m/gmres/ )
{
  $solver = 'gmres';
}
else
{
  print "Trouble parsing cnbse.solver!!\n*** Will default to  Haydock recursion ***\n";
  $solver = 'hay';
}
## Now if gmres we need to parse the inputs for that
my $gmres_footer = "";
my $gmres_count = 0;
if( $solver eq 'gmres' )
{
  open IN, "cnbse.gmres.nloop" or die "Failed to open cnbse.gmres.nloop\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  my $gmres_header = $line;

  open IN, "cnbse.broaden" or die "Failed to open cnbse.broaden\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  $line /= 27.2114;
  $gmres_header .= " " . $line;

  open IN, "cnbse.gmres.gprc" or die "Failed to open cnbse.gmres.gprc\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  $gmres_header .= " " . $line;

  open IN, "cnbse.gmres.ffff" or die "Failed to open cnbse.gmres.ffff\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  $gmres_header .= " " . $line;

  $gmres_header .= "  0.0\n";

  my $have_elist = 0;
  my $have_erange = 0;
  open IN, "cnbse.gmres.elist" or die "Failed to open cnbse.gmres.elist\n$!";
  $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    $gmres_footer = $gmres_header . "list\n";
    my $temp .= $line;
    my $i = 1;
    while( $line = <IN> )
    {
      $temp .= $line;
      $i++;
    }
    $gmres_count = $i;
    $gmres_footer .= "$i\n";
    $gmres_footer .= "$temp";
    $have_elist = 1;
    close IN;
  }

  open IN, "cnbse.gmres.erange" or die "Failed to open cnbse.gmres.erange\n$!";
  $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    if( $have_erange == 1 )
    {
      print "Both erange and elist were specified for GMRES. We are using erange\n";
    }
    else
    {
      $gmres_footer = $gmres_header . "loop\n";
      $line =~ m/([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/ or die "Failed to parse erange setting\n$line";
      my $estart = $1; my $estop = $3; my $estep = $5;
      $gmres_count = floor( ($estop - $estart + $estep*.9 ) / $estep );
      $gmres_footer .= $line;
      $have_erange = 1;
    }
    close IN;
  }

#  if( $have_erange + $have_elist == 2 )
#  {
#    print "Both erange and elist were specified for GMRES. We are using erange\n";
#  }
  if( $have_erange + $have_elist == 0 )
  {
    print "Neither elist nor erange were specified for GMRES!\n";
    print "*** WARNING ***\n";
    $have_elist = 1;
    $gmres_footer = $gmres_header . "list\n1\n0\n";
    $gmres_count = 1;
  }
}

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
    print "Parallel required\nCannot continue\n";
    exit 1;
  }
  close IN;
}
unless( -e "$ENV{'OCEAN_BIN'}/ocean.x" )
{
  print "Parallel BSE executable not present in $ENV{'OCEAN_BIN'}\n";
  print "Parallel required\nCannot continue\n";
  exit 1;
}



##### misc other setup
#`echo gmanner > format65`;
copy( "kmesh.ipt", "kgrid" ) or die "$!";
copy( "k0.ipt", "scaledkzero.ipt" ) or die "$!";
copy( "qinunitsofbvectors.ipt", "cksdq" ) or die "$!";

##### qin must be non-zero
open QIN, "qinunitsofbvectors.ipt" or die "Failed to open qinunitsofbvectors.ipt\n$!";
<QIN> =~ m/(\d*\.?\d+)\s+(\d*\.?\d+)\s+(\d*\.?\d+)/ or die "Failed to parse qinunitsofbvectors.ipt\n";
close QIN;
if( abs( $1 ) + abs( $2 ) + abs( $3 ) < 0.0000001 ) 
{
  print "Non-zero q required for VAL\n";
  exit 1;
}

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


# Set up mode -- need to run as xas which means ignore cnbse.mode
  my $is_xas = 1;
  open TMPFILE, "cnbse.niter" or die "Failed to open cnbse.niter\n$!";
  <TMPFILE> =~ m/(\d+)/ or die "Failed to parse cnbse.niter";
  my $num_haydock_iterations = $1;
  close TMPFILE;
  
  open TMPFILE, "cnbse.strength" or die "Failed to open cnbse.strength\n$!";
  <TMPFILE> =~ m/([0-9]*\.?[0-9]+)/ or die "Failed to parse cnbse.strength\n";
  my $interaction_strength = $1; 
  close TMPFILE;
    
# write cks.normal file
  open TMPFILE, ">cks.normal" or die "Failed to open cks.normal for writing\n$!";
  print TMPFILE ".true.\n";
  close TMPFILE;

#write mode file
  open TMPFILE, ">mode" or die "Failed to open mode for writing\n$!";
  print TMPFILE "$interaction_strength    $num_haydock_iterations\n";
  close TMPFILE;

###############
# If we are using QE/ABI w/o OBFs we need to set nbuse
my @brange;
my $run_text = '';
open NBUSE, "nbuse.ipt" or die "Failed to open nbuse.ipt\n";
<NBUSE> =~ m/(\d+)/ or die "Failed to parse nbuse.ipt\n";
my $nbuse = $1;
close NBUSE;
if( $obf == 1 )
{
  close RUNTYPE;
  if ($is_xas == 1 )
  {
    $run_text = 'XAS';
    if( $nbuse == 0 )
    {
      copy( "../zWFN/nbuse.ipt", "nbuse.ipt" ) or die "$!";
    }
    print "XAS!\n";
  } 
  else
  {
    if( $nbuse == 0 )
    {
      copy( "../zWFN/nbuse_xes.ipt", "nbuse.ipt" ) or die "$!";
    }
    $run_text = 'XES';
    print "XES!\n";
  }
}
else  ### Abi/QE w/o obf
{ 
#  my @brange;
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
      $nbuse = $brange[3] - $brange[1];
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
  else
  {
    if( $is_xas == 1 )
    {
      $run_text = 'XAS';
    }
    else
    {
      print "XES!\n";
      $run_text = 'XES';
    }
  }
}


system("$ENV{'OCEAN_BIN'}/getnval.x") == 0 or die "Failed to get nval\n";


#####################
if( $obf == 1 )
{
  if( -e "../zWFN/u2par.dat" )
  {
    `ln -s ../zWFN/u2par.dat`;
    `ln -s ../zWFN/ptmels.dat`;
    `cp ../zWFN/tmels.info .`;
    open OUT, ">bloch_type" or die;
    print OUT "new\n";
    close OUT;
  }
  else
  {
    `ln -s ../zWFN/u2.dat`;
  }
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
  system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";

  if (-e "../PREP/BSE/u2.dat")
  {
    `ln -s ../PREP/BSE/u2.dat`;
  }
  else
  {
    print "conugtoux\n";
    system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
    print "orthog\n";
    system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;
  }
}

my $gamma0 = `cat cnbse.broaden`;
chomp($gamma0);


open RUNLIST, ">runlist";
print RUNLIST "1\n0  0  0  __  __  0  0  VAL\n";
close RUNLIST;

# Set up for running the valence pathway for RIXS

print "\nSetting up valence section\n";


`tail -n 1 rhoofr > nfft`;

open INFILE, ">bse.in" or die "Failed to open bse.in\n";
print INFILE "0 0 0\n";
print INFILE "0 0 0\n";
my $spectrange = `cat spect.h`;
chomp($spectrange);
my $num_haydock_iterations = `cat niter`;
chomp($num_haydock_iterations);


print INFILE "hay\n";
print INFILE "$num_haydock_iterations  $spectrange  $gamma0  0.000\n";
close INFILE;

`cat bse.in`;

`echo 0 0 0 > ZNL`;

print "Running valence\n";
print "time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log";
system("time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log") == 0 or die "Failed to finish\n";


exit 0;

