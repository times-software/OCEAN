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

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/core_shift\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################

# Figure out our plan
my $offset;
my $control
if( -e "core_offset" )
{
  open IN, "core_offset" or die "Failed to open core_offset\n$!";
  $control = <IN>;
	chomp($control);
  if( $control =~ m/false/ )
	{
		print "How did I get here?\n";
    exit 0;
  }
  elsif( $control =~ m/true/ )
  {
    print "Offset given as true.\nWill return an average offsetof 0\n";
	}
  else
  {
    $offset = $control;
  }
}
else
{
  print "Can't find core_offset control file\n";
  exit 0;
}


# Load up all the radii we need
open RAD, "paw.shells" or die "Failed to open paw.shells\n";
my $line;
while( <RAD> )
{
  chomp;
  $line .= $_ . " ";
}
close RAD;
my @rads = split( /\s+/, $line );

#my $rad_dir = sprintf("zR%03.2f",$rad);

# Para prefix
open IN, "para_prefix" or die "Failed to open para_prefix\n$!";
my $para_prefix = <IN>;
chomp($para_prefix);
close IN;

# DFT flavor. Right now ABINIT or QE (OBF gets treated as QE)
open IN, "dft" or die "Failed to open dft\n$!";
my $line = <IN>;
close IN;
my $dft;
if( $line =~ m/abi/i )
{
  $dft = "abi";
}
elsif( $line =~ m/qe/i )
{
  $dft = "qe";
}
elsif( $line = ~m/obf/i )
{
  $dft = "qe"; # OBF if treated as QE for CLS!
}
else
{
  die "Failed to parse dft flavor\n";
}


open IN, "natoms" or die "Failed to open natoms\n$!";
my $natoms = <IN>;
chomp $natoms;
close IN;


# Bring in all of hfinlist
my @hfin;
open HFIN, "hfinlist" or die "Failed to open hfinlist";
while ( my $line = <HFIN>) 
{
# 07-n.lda.fhi                                         7   1   0 N_   1
  $line =~ m/\S+\s+\d+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die "Failed to parse hfin.\n$line";
  push @hfin, ( $1, $2, $3, $4 );
}
close HFIN;

# After this section DFT_pot will contain the full DFT potential evaluated for site
my @DFT_pot;

# For QE need to find their ridiculous coordinate formating
if( $dft eq 'qe' )
{
  my @coords;
  copy( "../DFT/scf.out", "scf.out" );

#`grep -A $natom "site" scf.out | tail -n $natom | awk '{print \$2, \$7, \$8, \$9}'   > xyz.alat`;
  open SCF, "scf.out" or die "Failed to open scf.out\n$!";

# skip down until we find 'site'
  while (<SCF>)
  {
    last if ($_ =~ m/site/ );
  }
# read in the next natom lines
  for( my $i=0; $i < $natom; $i++ )
  {
    $line = <SCF>;
    $line =~ m/\d+\s+(\w+)\s+tau\(\s*\d+\)\s+=\s+\(\s+(\S+)\s+(\S+)\s+(\S+)/ 
                or die "Failed to parse scf.out\n$line\nAtom index : " . $i+1 . "\n";
    print "$1\t$2\t$3\t$4\n";
    $coords[$i] = "$1\t$2\t$3\t$4\n";
  }
  close SCF;

  print "Pre-comp\n";
  open OUT, ">pot_prep.in";
  print OUT "&inputpp\n"
     .  "  prefix = 'system'\n"
     .  "  outdir = './Out'\n"
     .  "  filplot = 'system.pot'\n"
     .  "  plot_num = 1\n"
     .  "/\n";
  close OUT;
  system("$para_prefix $ENV{'OCEAN_BIN'}/pp.x < pot_prep.in > pot_prep.out");
  print "Pre-comp complete.\n";


  for( my $i = 0; $i < scalar @hfin; $i++ )
  {
  
    my $nn = $hfin[$i][0];
    my $ll = $hfin[$i][1];
    my $el = $hfin[$i][2];
    my $el_rank = $hfin[$i][3];

#  my $taustring = `grep $el xyz.wyck | head -n $el_rank | tail -n 1`;
    my $small_el = $el;
    $small_el =~ s/_//;

    my $taustring;
    my $count = 0;
    foreach $taustring ( @coords )
    {
      # For each element = small_el iterate count
      $count++ if( $taustring =~ m/$small_el/ );
      last if ( $count == $el_rank );
    }
#  my $taustring = `grep $small_el xyz.alat |  head -n $el_rank | tail -n 1`;
    print "$el_rank, $small_el, $taustring\n";
     $taustring =~ m/\S+\s+(\S+)\s+(\S+)\s+(\S+)/;
    my $x = $1;
    my $y = $2;
    my $z = $3;

    open OUT, ">pot.in" ;
    print OUT   "&inputpp\n/\n&plot\n"
     .  "  nfile = 1\n"
     .  "  filepp(1) = 'system.pot', weight(1) = 1\n"
     .  "  iflag = 1\n"
     .  "  output_format = 0\n"
     .  "  fileout = 'system.pot.$el_rank'\n"
     .  "  e1(1) = 0, e1(2) = 0, e1(3) = 1\n"
     .  "  x0(1) = $x"
     .  "  x0(2) = $y"
     .  "  x0(3) = $z"
     .  "  nx = 2\n"
     .  "/\n";
    close OUT;
   
    `cp pot.in pot.in.$el_rank`;
    system("$para_prefix $ENV{'OCEAN_BIN'}/pp.x < pot.in > pot.out.$el_rank");

# Vshift here is in Rydberg
  my $Vshift = `head -n 1 system.pot.$el_rank | awk '{print \$2}'`;


  my $string = sprintf("z%s%02d_n%02dl%02d",$el, $el_rank,$nn,$ll);
  print "$string\n";
# W shift is in Ha., but we want to multiple by 1/2 anyway, so the units work out
  my $Wshift = `head -n 1 $string/$rad_dir/ropt | awk '{print \$4}'`;

  my $shift = $Vshift + $Wshift;
  $shift *= 13.605;
#  if( $offset eq "no" )
#  {
#	$offset = -$shift;
#  }
  $shift += $offset;
  $shift *= -1;
  print "$el_rank\t$Vshift\t$Wshift\t$shift\n";

  print CORESHIFT "$shift\n";

}
close HFIN;
close CORESHIFT;

#`cp core_shift.txt ../`;
