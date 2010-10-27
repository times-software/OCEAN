#!/usr/bin/perl

use Getopt::Long;
use strict;


my @OceanFolders = ("Common", "ABINIT", "PREP", "PAW", "SCREEN", "CNBSE");

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
  $0 =~ m/(.*)\/par_ocean\.pl/;
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
system("$ENV{'OCEAN_BIN'}/parse --all $InputFile $ENV{'OCEAN_BIN'}/oparse.h") == 0 or die;
#if ($? != 0 ) {
#  print "Error parsing $InputFile: $?\n";
#}

chdir "../";
print "Done with parsing\n";
##########################################
#
# Abinit stage
##########################################
print "$Separator\n";
print "Entering ABINIT stage\n";
chdir "ABINIT";
system("$OCEAN_BIN/NewAbinitDriver.pl") == 0 or die "Abinit Stage Failed\n";
chdir "../";
exit;
##########################################
#
# DenDip stage
##########################################
print "$Separator\n";
print "Entering DenDip stage\n";
chdir "DenDip" or die "$!\n";
system("$OCEAN_BIN/dendip.pl") == 0 or die "DenDip Stage Failed\n";
##########################################
#
# PAW stage
##########################################
print "$Separator\n";
print "Entering PAW stage\n";
chdir "../PAW";
system("$OCEAN_BIN/paw2.pl") == 0 or die "PAW stage failed\n";
##########################################
#
# CNBSE stage
##########################################
print "$Separator\n";
print "Entering CNBSE stage\n";
chdir "../CNBSE";
system("$OCEAN_BIN/cnbse.pl") == 0 or die "CNBSE stage failed\n";

##########################################
print "$Separator\n";
print "Ocean is done\n";
print "$Separator\n";



exit 0;
