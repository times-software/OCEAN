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
  $0 =~ m/(.*)\/ocean_finish\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
my $OCEAN_BIN = $ENV{"OCEAN_BIN"};

$ENV{"OCEAN_WORKDIR"} = `pwd`;
$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`;
##########################################
#print ', version "$ENV{OCEAN_VERSION}"\n';


##########################################
#
# Prep stage
##########################################
print "$Separator\n";
print "Entering PREP stage\n";
chdir "PREP" or die "$!\n";
system("$OCEAN_BIN/par_dendip.pl") == 0 or die "PREP Stage Failed\n";
##########################################
#
# PAW stage
##########################################
print "$Separator\n";
print "Entering PAW stage\n";
chdir "../PAW";
system("$OCEAN_BIN/paw.pl") == 0 or die "PAW stage failed\n";
##########################################
#
# SCRENing stage
##########################################
print "$Separator\n";
print "Entering SCREENing stage\n";
chdir "../SCREEN";
system("$OCEAN_BIN/screen.pl") == 0 or die "SCREEN stage failed\n";
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
