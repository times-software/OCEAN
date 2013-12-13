#!/usr/bin/perl

use strict;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/core_shift\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################



#my $offset = "no"; #121.611826250456;
my $para_prefix = `cat para_prefix`;
chomp($para_prefix);
my $offset = 121.611826250456;
open HFIN, "hfinlist" or die "Failed to open hfinlist";


if( -e "core_offset" )
{
	$offset = `cat core_offset`;
	chomp($offset);
        if( $offset =~ m/false/ )
	{
		print "How did I get here?\n";
		die;
	}
}

my $yderp = 1;   #0.946557;
my $zderp = 1;   #0.860222;

open CORESHIFT, ">core_shift.txt";

`cp ../DFT/scf.out .`;

my $natom = `cat natoms`;
chomp($natom);
#`grep -A $natom "site" scf.out | tail -n $natom | awk '{print \$2, \$7, \$8, \$9}'   > xyz.alat`;
open SCF, "scf.out" or die "$!\n";
while (<SCF>)
{
  last if ($_ =~ m/site/ );
}
open ALAT, ">xyz.alat" or die "$_\n";
for( my $i=0; $i < $natom; $i++ )
{
  my $line = <SCF>;
  $line =~ m/\d+\s+(\w+)\s+tau\(\s*\d+\)\s+=\s+\(\s+(\S+)\s+(\S+)\s+(\S+)/;
  print "$1\t$2\t$3\t$4\n";
  print  ALAT "$1\t$2\t$3\t$4\n";
}
close ALAT;
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
system("mpirun -n 8 $ENV{'OCEAN_BIN'}/pp.x < pot_prep.in >& pot_prep.out");

while ( my $line = <HFIN>) 
{
# 07-n.lda.fhi                                         7   1   0 N_   1
  $line =~ m/\S+\s+\d+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die "Failed to parse hfin.\n$line";
  my $nn = $1;
  my $ll = $2;
  my $el = $3;
  my $el_rank = $4;

#  my $taustring = `grep $el xyz.wyck | head -n $el_rank | tail -n 1`;
  my $small_el = $el;
  $small_el =~ s/_//;
  my $taustring = `grep $small_el xyz.alat |  head -n $el_rank | tail -n 1`;
  print "$el_rank, $small_el, $taustring\n";
   $taustring =~ m/\S+\s+(\S+)\s+(\S+)\s+(\S+)/;
  my $x = $1;
  my $y = $2;
  my $z = $3;

  $y *=$yderp;
  $z *= $zderp;

  open OUT, ">pot.in" ;
  print OUT   "&inputpp\n/\n&plot\n"
#"&inputpp\n"
#   .  "  prefix = 'system'\n"
#   .  "  outdir = './Out'\n"
#   .  "  filplot = 'system.pot'\n"
#   .  "  plot_num = 1\n"
#   .  "/\n"
#   .  "&plot\n"
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
 
#  system("$para_prefix ~/bin/gnu/OCEAN_atlas/par_pp.x < pot.in >& pot.out.$el_rank"); 
#  system("~/bin/gnu/OCEAN_atlas/pp.x < pot.in >& pot.out.$el_rank");
  `cp pot.in pot.in.$el_rank`;
  system("mpirun -n 8 $ENV{'OCEAN_BIN'}/pp.x < pot.in >& pot.out.$el_rank");

# Vshift here is in Rydberg
  my $Vshift = `head -n 1 system.pot.$el_rank | awk '{print \$2}'`;


  my $string = sprintf("z%s%02d_n%02dl%02d",$el, $el_rank,$nn,$ll);
  print "$string\n";
# W shift is in Ha., but we want to multiple by 1/2 anyway, so the units work out
  my $Wshift = `head -n 1 $string/zR4.00/ropt | awk '{print \$4}'`;

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
