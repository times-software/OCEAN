use strict;

#my $offset = "no"; #121.611826250456;
my $offset = 121.611826250456;
open HFIN, "hfinlist" or die "Failed to open hfinlist";


if( -e "offset_overide" )
{
	$offset = `cat offset_overide`;
	chomp($offset);
}

my $yderp = 1;   #0.946557;
my $zderp = 1;   #0.860222;

open CORESHIFT, ">core_shift.txt";

while ( my $line = <HFIN>) 
{
# 07-n.lda.fhi                                         7   1   0 N_   1
  $line =~ m/\S+\s+\d+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die "Failed to parse hfin.\n$line";
  my $nn = $1;
  my $ll = $2;
  my $el = $3;
  my $el_rank = $4;

  my $taustring = `grep $el xyz.wyck | head -n $el_rank | tail -n 1`;
  print "$el_rank, $taustring";
   $taustring =~ m/\S+\s+(\S+)\s+(\S+)\s+(\S+)/;
  my $x = $1;
  my $y = $2;
  my $z = $3;

  $y *=$yderp;
  $z *= $zderp;

  open OUT, ">pot.in" ;
  print OUT "&inputpp\n"
   .  "  prefix = 'system'\n"
   .  "  outdir = './Out'\n"
   .  "  filplot = 'system.pot'\n"
   .  "  plot_num = 1\n"
   .  "/\n"
   .  "&plot\n"
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
 
  system("~/bin/gnu/OCEAN/pp.x < pot.in >& pot.out.$el_rank"); 

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
  print "$el_rank\t$Vshift\t$Wshift\t$shift\n";

  print CORESHIFT "$shift\n";

}
close HFIN;
close CORESHIFT;

#`cp core_shift.txt ../`;
