#!/usr/bin/perl

use Getopt::Long;
use Switch;
use strict;


my @OceanFolders = ("Common", "DFT", "zWFN", "PAW", "SCREEN", "CNBSE", "ABINIT", "QESPRESSO", "PREP");

print "Welcome to OCEAN\n";

##########################################
my $ColWidth = `tput cols`;
if ($ColWidth < 0 ) { $ColWidth=10; }
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
  $0 =~ m/(.*)\/ocean\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
my $OCEAN_BIN = $ENV{"OCEAN_BIN"};

$ENV{"OCEAN_WORKDIR"} = `pwd`;
$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`;
##########################################
#print ', version "$ENV{OCEAN_VERSION}"\n';


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
`cp ../$InputFile .`;
system("$ENV{'OCEAN_BIN'}/parse --all $InputFile $ENV{'OCEAN_BIN'}/oparse.h") == 0 
  or die "Failed to parse the input file\n$!";
#if ($? != 0 ) {
#  print "Error parsing $InputFile: $?\n";
#}
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


### remove print statments in this section after testing
##
print "Opening Stages\n";

open(STAGES, "<stages");
my $tmp;
my @stages;
while(<STAGES>)
{
  chomp;
  $tmp = $_;
  @stages = split(/ /,$tmp);
}

print "@stages\n";
print scalar @stages;
print "\n";

my $runstage = "00000\n";
print $runstage;

foreach my $stage (@stages) {
  print "$stage\n";
  switch ($stage) {
    case "paw" {substr($runstage,0,1) = "1";}
    case "dft" {substr($runstage,1,1) = "1";}
    case "prep" {substr($runstage,2,1) = "1";}
    case "screen" {substr($runstage,3,1) = "1";}
    case "bse" {substr($runstage,4,1) = "1\n";}
    case "all" {$runstage = "11111\n";}
  }
  print "$runstage";
}

print "Closing Stages\n";
##
###


chdir "../";
print "Done with parsing\n";
##########################################
#
# PAW stage
##########################################
if (substr($runstage,0,1) eq "1") {
  print "$Separator\n";
  print "Entering PAW stage\n";
  chdir "PAW";
  system("$OCEAN_BIN/paw.pl >& paw.log") == 0 or die "PAW stage failed\n";
  chdir "../";
}
##########################################
#
# DFT stage
##########################################
if (substr($runstage,1,1) eq "1") {

  if( $script_pre eq 'OBF' ) 
  {
	print "$Separator\n";
	print "Entering DFT stage\n";
	chdir "DFT";
	system("$OCEAN_BIN/${script_pre}_dft.pl >& dft.log") == 0 or die "DFT Stage Failed\n$!";
	chdir "../";
  } 
  elsif( $script_pre eq 'qe' )
  {
	print "$Separator\n";
        print "Entering QESPRESSO stage\n";
        chdir "QESPRESSO";
        system("$OCEAN_BIN/QespressoDriver.pl") == 0 or die "Qespresso Stage Failed\n";
        chdir "../";
  }
  else
  {
	print "$Separator\n";
	print "Entering ABINIT stage\n";
	chdir "ABINIT";
	system("$OCEAN_BIN/AbinitDriver.pl") == 0 or die "Abinit Stage Failed\n";
        chdir "../";
  }

}
##########################################
#
# zWFN stage
##########################################
if (substr($runstage,2,1) eq "1") {

  if( $script_pre eq 'OBF' ) 
  {
	print "$Separator\n";
	print "Entering zWFN stage\n";
	chdir "zWFN" or die "$!\n";
	system("$OCEAN_BIN/${script_pre}_wfn.pl >& wfn.log") == 0 or die "zWFN Stage Failed\n$!";
  }
  else
  {
	print "$Separator\n";
	print "Entering PREP stage\n";
	chdir "PREP" or die "$!\n";
	if( $script_pre eq 'qe' )
	{
  		system("$OCEAN_BIN/qe_dendip.pl >& prep.log ") == 0 or die "PREP Stage Failed\n$!";
	}
	else
	{
		system("$OCEAN_BIN/dendip.pl >& prep.log ") == 0 or die "PREP Stage Failed\n$!";
	}
  }

}
##########################################
#
# SCREENing stage
##########################################
if (substr($runstage,3,1) eq "1") {

  print "$Separator\n";
  print "Entering SCREENing stage\n";
  chdir "../SCREEN";

  if( $script_pre eq 'abi' )
  {
	system("$OCEAN_BIN/screen.pl >& screen.log") == 0 or die "SCREEN stage failed\n$!";
  }
  else
  {
        system("$OCEAN_BIN/screen.pl >& screen.log") == 0 or die "SCREEN stage failed\n$!";
#	system("$OCEAN_BIN/${script_pre}_screen.pl >& screen.log") == 0 or die "SCREEN stage failed\n$!";
  }

}
##########################################
#
# CNBSE stage
##########################################
if (substr($runstage,4,1) eq "1") {

  print "$Separator\n";
  print "Entering CNBSE stage\n";
  chdir "../CNBSE";
  if( $script_pre eq 'OBF' )
  {
        system("$OCEAN_BIN/cnbse.pl >& cnbse.log") == 0 or die "CNBSE stage failed\n$!";
#	system("$OCEAN_BIN/${script_pre}_cnbse.pl >& cnbse.log") == 0 or die "CNBSE stage failed\n$!";
  }
  else
  {
	system("$OCEAN_BIN/cnbse.pl >& cnbse.log") == 0 or die "CNBSE stage failed\n$!";
  }

}
##########################################
print "$Separator\n";
print "Ocean is done\n";
print "$Separator\n";



exit 0;
