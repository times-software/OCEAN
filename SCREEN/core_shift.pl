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

my $Ry2eV = 13.605698066;

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
my $control;
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
    print "Offset given as true.\nWill return an average offset of 0\n";
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
#  die "Core level shift support not written for abinit yet!\n";
}
elsif( $line =~ m/qe/i )
{
  $dft = "qe";
}
elsif( $line = ~m/obf/i )
{
  $dft = "obf"; 
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
  push @hfin, [( $1, $2, $3, $4 )];
}
close HFIN;

# After this section DFT_pot will contain the full DFT potential evaluated for site
my @Vshift;
my @Wshift;
my $Vsum = 0;
my @Wsum;
for( my $i = 0; $i < scalar @rads; $i++ )
{
  $Wsum[$i] = 0;
}

# For QE need to find their ridiculous coordinate formating
if( $dft eq 'qe' || $dft eq 'obf' )
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
  for( my $i=0; $i < $natoms; $i++ )
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
  if( $dft eq 'qe' )
  {
    system("$para_prefix $ENV{'OCEAN_BIN'}/pp.x < pot_prep.in > pot_prep.out") == 0 
        or die "Failed to run pp.x for the pre-computation\n";
  }
  else
  {
    system("$para_prefix $ENV{'OCEAN_BIN'}/obf_pp.x < pot_prep.in > pot_prep.out") == 0 
        or die "Failed to run pp.x for the pre-computation\n";
  }
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
    for( my $j = 0; $j < scalar @coords; $j++ )
    {
      $taustring = $coords[$j];
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
    print "$x\t$y\t$z\n";

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
    if( $dft eq 'qe' )
    {
      system("$para_prefix $ENV{'OCEAN_BIN'}/pp.x < pot.in > pot.out.$el_rank") == 0
          or die "Failed to run pp.x for $el_rank\n$!";
    }
    else
    {
      system("$para_prefix $ENV{'OCEAN_BIN'}/obf_pp.x < pot.in > pot.out.$el_rank") == 0
          or die "Failed to run pp.x for $el_rank\n$!";
    }

  # Vshift here is in Rydberg
    $Vshift[$i] = `head -n 1 system.pot.$el_rank | awk '{print \$2}'`;
    chomp( $Vshift[$i] );
    $Vsum += $Vshift[$i];


    my $string = sprintf("z%s%04d_n%02dl%02d",$el, $el_rank,$nn,$ll);
    print "$string\n";
  # W shift is in Ha., but we want to multiple by 1/2 anyway, so the units work out

    for( my $j = 0; $j < scalar @rads; $j++ )
    {
      my $rad_dir = sprintf("zR%03.2f", $rads[$j] );
      my $temp = `head -n 1 $string/$rad_dir/ropt | awk '{print \$4}'`;
      chomp( $temp );
      $Wshift[$i][$j] = $temp;
      $Wsum[$j] += $temp;
    }
  }
}
else
{
  copy( "../DFT/SCx_POT", "SCx_POT" );
  my @coords;

  open IN, "xyz.wyck" or die "Failed to open xyz.wyck.\n$!";
  <IN>;
  while (<IN>)
  {
    push @coords, $_;
  }
  close IN;

  
  open OUT, ">pot.in" or die "Failed top open pot.in for writing.\n$!";

  print OUT "SCx_POT\n1\n";

  for( my $i = 0; $i < scalar @hfin; $i++ )
  {

    my $nn = $hfin[$i][0];
    my $ll = $hfin[$i][1];
    my $el = $hfin[$i][2];
    my $el_rank = $hfin[$i][3];

    my $taustring;
    my $count = 0;
    for( my $j = 0; $j < scalar @coords; $j++ )
    {
      $taustring = $coords[$j];
      # For each element = small_el iterate count
      $count++ if( $taustring =~ m/$el/ );
      last if ( $count == $el_rank );
    }
    print "$el_rank, $el, $taustring\n";
     $taustring =~ m/\S+\s+(\S+)\s+(\S+)\s+(\S+)/;
    my $x = $1; 
    my $y = $2;
    my $z = $3;
    print "$x\t$y\t$z\n";
    
    print OUT "1\n"
            . "2\n"
            . "$x $y $z\n";
    if( ($i+1) < scalar @hfin )
    {
      print OUT "1\n";
    }
    else
    {
      print OUT "0\n";
    }
  }
  close OUT;

  system( "$ENV{'OCEAN_BIN'}/cut3d < pot.in > pot.out") == 0 or die "Failed to run cut3d\n$!";

  open IN, "pot.out" or die "Failed to open pot.out\n$!";
  my $line;
  while( $line = <IN> )
  {
    if( $line =~ /value=\s+(\S+)/ )
    {
      # Vshift is in Ha not Ry for ABINIT
      push @Vshift, 2*$1;
      $Vsum += 2*$1;
    }
  }
  close IN;
  

  for( my $i = 0; $i < scalar @hfin; $i++ )
  {

    my $nn = $hfin[$i][0];
    my $ll = $hfin[$i][1];
    my $el = $hfin[$i][2];
    my $el_rank = $hfin[$i][3];

    my $string = sprintf("z%s%04d_n%02dl%02d",$el, $el_rank,$nn,$ll);
    print "$string\n";

    for( my $j = 0; $j < scalar @rads; $j++ )
    {
      my $rad_dir = sprintf("zR%03.2f", $rads[$j] );
      my $temp = `head -n 1 $string/$rad_dir/ropt | awk '{print \$4}'`;
      chomp( $temp );
      $Wshift[$i][$j] = $temp;
      $Wsum[$j] += $temp;
    }
  }
  
}


# Loop over radii and then hfin
for( my $i = 0; $i < scalar @rads; $i++ )
{
  my $rad_dir = sprintf("zR%03.2f", $rads[$i] );

  # If we are averaging, new shift by radius
  if( $control =~ m/true/ )
  {
    $offset = -( $Vsum + $Wsum[$i] ) * $Ry2eV / ( scalar @hfin );
    print "$rad_dir\t$offset\n";
  }

  # Loop over each atom in hfin
  for( my $j = 0; $j < scalar @hfin; $j++ )
  {
    my $nn = $hfin[$j][0];
    my $ll = $hfin[$j][1];
    my $el = $hfin[$j][2];
    my $el_rank = $hfin[$j][3];

    my $shift = ( $Vshift[$j] + $Wshift[$j][$i] ) * $Ry2eV;

    $shift += $offset;
    $shift *= -1;
    print "$el_rank\t$Vshift[$j]\t$Wshift[$j][$i]\t$shift\n";

    my $string = sprintf("z%s%04d_n%02dl%02d",$el, $el_rank,$nn,$ll);
    open OUT, ">$string/$rad_dir/cls" or die "Failed to open $string/$rad_dir/cls\n$!";
    print OUT $shift . "\n";
    close OUT;
  }

}

exit 0;

