#!/usr/bin/perl
# Copyright (C) 2014 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use Getopt::Long;
use Cwd 'abs_path';
use Cwd 'getcwd';
use File::Copy;
use strict;


my @OceanFolders = ("Common", "DFT", "zWFN", "OPF", "SCREEN", "CNBSE", "PREP", "NBSE" );

print "Welcome to OCEAN\n";

##########################################
#???
my $ColWidth = sprintf('%i',`tput cols` );
if ($ColWidth < 0 ) { $ColWidth=10; }
if ($ColWidth > 100 ) { $ColWidth=100; }
my $Separator =  "";
for (my $i = 0; $i < $ColWidth; $i++ ) {
  $Separator  .= "#";
}

if ($ColWidth > 60 ) {
  print $Separator . "\n";
  print "#                 OOO   CC  EEEE  AA  N  N                 #\n";
  print "#                O   O C  C E    A  A NN N                 #\n";
  print "#      **        O   O C    E    A  A NN N         **      #\n";
  print "#     **         O   O C    EE   AAAA N NN          **     #\n";
  print "#    ***    **   O   O C    E    A  A N NN    **    ***    #\n";
  print "#   *****  ***   O   O C  C E    A  A N  N    ***  *****   #\n";
  print "# **************  OOO   CC  EEEE A  A N  N  ************** #\n";
  print $Separator . "\n";
}

if (! $ENV{"OCEAN_BIN"} ) {
  abs_path($0) =~ m/(.*)\/ocean\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
}
my $OCEAN_BIN = $ENV{"OCEAN_BIN"};

$ENV{"OCEAN_WORKDIR"} = getcwd();
if ( open VERSION, "$ENV{'OCEAN_BIN'}/Version" )
{
  my $version = <VERSION>;
  chomp( $version );
  $ENV{"OCEAN_VERSION"} = $version;
  if( $version = <VERSION> )
  {
    chomp( $version );
    $ENV{"OCEAN_VERSION_HASH"} = $version;
  }
}
else
{
  print "Version file not found\n";
  $ENV{"OCEAN_VERSION"} = "NaN";
}
close VERSION;

##########################################
print "Version $ENV{'OCEAN_VERSION'}\n";
if( length $ENV{"OCEAN_VERSION_HASH"} > 1 )
{
  print "  commit hash: $ENV{'OCEAN_VERSION_HASH'}\n";
}


# Process Options
##########################################
my $OPT_Clean = '1';
my $OPT_CleanSaved = '';
my $OPT_SaveMin = '';

GetOptions ('clean' => \$OPT_Clean, 'cleansaved' => \$OPT_CleanSaved, 'min' => \$OPT_SaveMin);
##########################################

# Clean
##########################################
if ($OPT_CleanSaved) {
  print "Clear all saved folders? (y/n)\n";
  if (getc eq "y") {
    `rm -r Save.????`;
  }
  else {
    exit 1;
  }
  $OPT_Clean = 1;
}
##########################################


# Get input file
##########################################
my $InputFile = $ARGV[0];
unless ( -e $InputFile ) { die "$InputFile not found\n"; }
##########################################

# Test level of differ, move to save folder (depending on flags)
##########################################
unless ($OPT_Clean) {
  if (-e "Common/") {
    my $tmpline = `ls -ld Save.*`;
    my @SaveFolder = split(/'\n'/, $tmpline);
    $SaveFolder[-1] =~ m/Save\.(\d+)/;
    my $NewSave = sprintf("Save.%04i",$1+1);
    print "$NewSave\n";
    foreach (@OceanFolders) {
      `cp -r $_ $NewSave`;
    }
  }
}

if (-e "Common") {
  if (-e "Common/$InputFile" && -e "Common/Inp_Clean") {
    `mv Common/Inp_Clean Common/Inp_old`;
    `cp $InputFile Common/`;
  }
  else { print "Run differs\n"; }
} 
else {
  print "No previous run found\n";
  foreach (@OceanFolders) {
    `rm -r $_`;
  }
}

foreach (@OceanFolders) {
  `mkdir -p $_`;
}

##########################################
print "Setup complete, parsing ...\n";

# Common stage
##########################################

chdir "Common";
#`cp ../$InputFile .`;
copy("../$InputFile","$InputFile");

system("$ENV{'OCEAN_BIN'}/parse.pl --all $InputFile $ENV{'OCEAN_BIN'}/oparse.h") == 0 
  or die "Failed to parse the input file\n$!";

open DFT_TYPE, "dft";
<DFT_TYPE> =~ m/(\w+)/ or die;
my $dft_type = $1;
close DTF_TYPE;
my $script_pre;
if( $dft_type =~ m/abinit/i )
{
  $script_pre = 'abi';
} 
elsif( $dft_type =~ m/qe/i )
{
  $script_pre = 'qe';
} 
elsif( $dft_type =~ m/obf/i )
{
  $script_pre = 'OBF';
}
else
{
  print "WARNING!!! DFT helper program unspecified. Using ABINIT\n";
  $script_pre = 'abi';
}

system("$ENV{'OCEAN_BIN'}/defaults.pl") == 0 or die "Failed to run defaults.pl\n$!";

system("$ENV{'OCEAN_BIN'}/structure.pl") == 0 or die "Failed to run structure.pl\n$!";

system("$ENV{'OCEAN_BIN'}/edges.pl") == 0 or die "Failed to run edges.pl\n$!";

### CALC ###

# The type of calculation will determine a bit about what scripts/stages get run
# In the future this may be more complicated
open CALC, "calc" or die "Failed to open file calc: $!\n";
<CALC> =~ m/(\w+)/ or die "Failed to read file calc\n";
my $calc = lc($1);
close CALC;
my $run_opf;
my $run_screen;
if( $calc =~ m/val/i )
{
  $run_opf = 0;
  $run_screen = 0;
}
else
{
  $run_opf = 1;
  $run_screen = 1;
}

my $run_rixs = 0;
if( $calc =~ m/rixs/i )
{
  $run_rixs = 1;
}

### \CALC ###


### EXCITON ###

open EXC, "plotexc" or die "Failed to open file plotexc: $!\n";
<EXC> =~ m/(\w+)/ or die "Failed to read file plotexc\n";
my $excp = lc($1);
close EXC;
my $run_val_exc = 0;
my $run_core_exc = 0;

if( $excp =~ m/true/i )
{
    if( $calc =~ m/val/i )
    {
        $run_val_exc = 1;
    }
    elsif( $calc =~ m/xas/i )
    {
        $run_core_exc = 1;
    }
    elsif( $calc =~ m/xes/i )
    {
        $run_core_exc = 1;
    }
    elsif( $calc =~ m/rixs/i )
    {
        $run_core_exc = 1;
        $run_val_exc = 1;
    }
}

### \EXCITON ###


### CHECK STATUS ###

my $run_paw = 0;
if ( open STATUS, "../OPF/status" )
{
   $run_paw = <STATUS>;
   chomp( $run_paw );
   print "DFT STATUS : $run_paw\n";
   close STATUS;
}


my $run_dft = 0;
if ( open STATUS, "../DFT/status" )
{
   $run_dft = <STATUS>;
   chomp( $run_dft );
   print "DFT STATUS : $run_dft\n";
   close STATUS;
}


my $run_wfn = 0;
if ( open STATUS, "../PREP/status" )
{
   $run_wfn = <STATUS>;
   chomp( $run_wfn );
   print "PREP STATUS : $run_wfn\n";
   close STATUS;
}

if ( open STATUS, "../zWFN/status" )
{
   $run_wfn = <STATUS>;
   chomp( $run_wfn );
   print "zWFN STATUS : $run_wfn\n";
   close STATUS;
}


my $run_scr = 0;
if ( open STATUS, "../SCREEN/status" )
{
   $run_scr = <STATUS>;
   chomp( $run_scr );
   print "SCREEN STATUS : $run_scr\n";
   close STATUS;
}


my $run_vbse = 0;
if ( open STATUS, "../NBSE/status" )
{
   $run_vbse = <STATUS>;
   chomp( $run_vbse );
   print "NBSE STATUS : $run_vbse\n";
   close STATUS;
}

my $run_cbse = 0;
if ( open STATUS, "../CNBSE/status" )
{
   $run_cbse = <STATUS>;
   chomp( $run_cbse );
   print "CNBSE STATUS : $run_cbse\n";
   close STATUS;
}


# if RIXS folder exists see if it ran successfully
# i.e. skip if we are plotting excitons after running RIXS
my $run_rxs = 0;
if ( -d "../RIXS/") {
   if ( open STATUS, "../RIXS/status" )
   {
      $run_rxs = <STATUS>;
      chomp( $run_rxs );
      print "RIXS STATUS : $run_rxs\n";
      close STATUS;
   }
}

### \CHECK STATUS ###



chdir "../";
print "Done with parsing\n";
##########################################
#
# OPF stage
##########################################
print "$Separator\n";
if( $run_paw == 0 )
{
  if( $run_opf )
  {
    print "Entering OPF stage\n";
    chdir "OPF";
    print_status(0);
    system("$OCEAN_BIN/opf.pl") == 0 or die "OPF stage failed\n";
    print_status(1);
    chdir "../";
  }
}
else
{
    print "No need to run PAW stage\n";
}
##########################################
#
# DFT stage
##########################################
if( $run_dft == 0 )
{
  if( $script_pre eq 'OBF' )
  {
    print "$Separator\n";
    print "Entering DFT stage\n";
    chdir "DFT";
    print_status(0);
    system("$OCEAN_BIN/dft.pl") == 0 or die "DFT Stage Failed\n$!";
    print_status(1);
    chdir "../";
  }
  elsif( $script_pre eq 'qe' )
  {
    print "$Separator\n";
    print "Entering QESPRESSO stage\n";
    chdir "DFT";
    print_status(0);
    system("$OCEAN_BIN/dft.pl") == 0 or die "Qespresso Stage Failed\n";
    print_status(1);
    chdir "../";
  }
  else
  {
    print "$Separator\n";
    print "Entering ABINIT stage\n";
    chdir "DFT";
    print_status(0);
    system("$OCEAN_BIN/AbinitDriver.pl") == 0 or die "Abinit Stage Failed\n";
    print_status(1);
    chdir "../";
  }
}
else
{
    print "No need to run DFT stage\n";
}
##########################################
#
# zWFN stage
##########################################
if( $run_wfn == 0 )
{
  if( $script_pre eq 'OBF' )
  {
    print "$Separator\n";
    print "Entering zWFN stage\n";
    chdir "zWFN" or die "$!\n";
    print_status(0);
    system("$OCEAN_BIN/${script_pre}_wfn.pl") == 0 or die "zWFN Stage Failed\n$!";
    print_status(1);
    chdir "../";
  }
  else
  {
    print "$Separator\n";
    print "Entering PREP stage\n";
    chdir "PREP" or die "$!\n";
    print_status(0);
    if( $script_pre eq 'qe' )
    {
        system("$OCEAN_BIN/qe_dendip.pl") == 0 or die "PREP Stage Failed\n$!";
    }
    else
    {
        system("$OCEAN_BIN/dendip.pl") == 0 or die "PREP Stage Failed\n$!";
    }
    print_status(1);
    chdir "../";
  }
}
else
{
      print "No need to run WFN stage\n";
}
##########################################
#
# SCREENing stage
##########################################
print "$Separator\n";
if( $run_screen)
{
  if( $run_scr == 0 )
  {
    print "Entering SCREENing stage\n";
    chdir "SCREEN";
    print_status(0);
    if( $script_pre eq 'OBF' )
    {
      system("$OCEAN_BIN/${script_pre}_screen_multi.pl") == 0 or die "SCREEN stage failed\n$!";
    }
    else
    {
      system("$OCEAN_BIN/screen.pl") == 0 or die "SCREEN stage failed\n$!";
    }
    print_status(1);
    chdir "../";
  }
  else
  {
    print "No need to run core SCREEN stage\n";
  }
}
##########################################
#
# CNBSE stage
##########################################
print "$Separator\n";
if( $calc =~ m/val/ )
{
  if( $run_vbse == 0 )
  {
    print "Entering NBSE stage\n";
    chdir "NBSE";
    print_status(0);
    system("$OCEAN_BIN/nbse.pl") == 0 or die "NBSE stage failed\n$!";
    print_status(1);
    chdir "../";
  }
  else
  {
    print "No need to run valence BSE stage\n";
  }
}
else
{
  if( $run_cbse == 0 )
  {
    print "Entering CNBSE stage\n";
    chdir "CNBSE";
    print_status(0);
    system("$OCEAN_BIN/cnbse_mpi.pl") == 0 or die "CNBSE stage failed\n$!";
    print_status(1);
    chdir "../";
  }
  else
  {
    print "No need to run core BSE stage\n";
  }
}
##########################################
#
# RIXS stage
##########################################
print "$Separator\n";
if( $run_rixs == 1 )
{
  if ( $run_rxs == 0 )
  {
    print "Entering RIXS stage\n";
    mkdir "RIXS" unless -d "RIXS";
    chdir "RIXS";
    print_status(0);
    `echo gmres > ../Common/cnbse.solver`;
    `echo xas > ../Common/calc`;
    system("$OCEAN_BIN/rixs.pl > rixs.log") == 0 or die "RIXS stage failed\n$!";
    print_status(1);
    chdir "../";
  }
  else
  {
    print "No need to run RIXS stage\n";
  }
}
##########################################
#
# EXCITON stage
##########################################
print "$Separator\n";
if( $run_core_exc == 1 )
{
  print "Entering core EXCITON stage\n";
  `mkdir -p EXCITON`;
  chdir "EXCITON";
#  print_status(0);
  system("$OCEAN_BIN/core-exciton.pl > exciton.log") == 0 or die "EXCITON stage failed\n$!";
#  print_status(1);
  chdir "../";
}
if( $run_val_exc == 1 )
{
  print "Entering valence EXCITON stage\n";
  `mkdir -p EXCITON`;
  chdir "EXCITON";
#  print_status(0);
  if( $calc =~ m/val/i )
  {
    system("$OCEAN_BIN/optical-exc.pl > exciton.log") == 0 or die "EXCITON stage failed\n$!";
  }
  elsif( $calc =~ m/rixs/i )
  {
    system("$OCEAN_BIN/rixs-val-exc.pl > exciton.log") == 0 or die "EXCITON stage failed\n$!";
  }
#  print_status(1);
  chdir "../";
}
##########################################
print "$Separator\n";
print "Ocean is done\n";
print "$Separator\n";



exit 0;


sub print_status {
   if (-e "status") {
   unlink('status') or die "Failed to delete status!\n";
   }
   open(my $runfh, '>', "status") or die "Could not open file 'status' $!";
   print $runfh "$_[0]\n";
   close $runfh;
}

