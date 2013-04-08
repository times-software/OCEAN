#!/usr/bin/perl

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/OBF_cnbse\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }

my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

my @CommonFiles = ("epsilon", "xmesh.ipt", "nedges", "k0.ipt", #"nbuse.ipt", 
  "cnbse.rad", "cnbse.ways", "metal", "cksshift", "cksstretch", "cksdq", "cks.normal",
  "cnbse.niter", "cnbse.spect_range", "cnbse.broaden", "cnbse.mode");

my @AbinitFiles = ("avecsinbohr.ipt");

my @DFTFiles = ("nelectron");

my @DenDipFiles = ("kmesh.ipt", "masterwfile", "listwfile", "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "tmels", "nelectron", "eshift.ipt" );

my @WFNFiles = ("kmesh.ipt",  "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", "avecsinbohr.ipt", "nbuse.ipt", "wvfcninfo", "wvfvainfo", "nbuse_xes.ipt");

my @ExtraFiles = ("Pquadrature", "sphpts" );

my @jtv = ("jtv1");
my @ways = ( 1, 2, 3 );

my @PawFiles = ("hfinlist" , "xyz.wyck");
#`cp ../xyz.wyck .`;
foreach (@CommonFiles) {
  `cp ../Common/$_ .` == 0 or die "Failed to get Common/$_\n";
}
#foreach (@AbinitFiles) {
#  `cp ../SCREEN/$_ .` == 0 or die "Failed to get ABINIT/$_\n";
#}
#foreach (@DenDipFiles) {
#  `cp ../PREP/BSE/$_ .` == 0 or die "Failed to get PREP/BSE/$_\n" ;
#}
foreach (@DFTFiles) {
  `cp ../DFT/$_ . ` == 0 or die "Failed to get DFT/$_\n";
}
foreach (@WFNFiles) {
  `cp ../zWFN/$_ .`== 0 or die "Failed to get zWFN/$_\n";
}
foreach (@ExtraFiles) {
  `cp $ENV{'OCEAN_BIN'}/$_ .` == 0 or die "Failed to get ../$_\n";
}
foreach (@ways) {
  `cp ../jtv$_ .` == 0 or die "Failed to get ../$_\n";
}

foreach (@PawFiles) {
  `cp ../SCREEN/$_ .` == 0 or die "Failed to get ../SCREEN/$_\n";
}

##### misc other setup
#`echo gmanner > format65`;
`cp kmesh.ipt kgrid`;
`cp k0.ipt scaledkzero.ipt`;
`mv cnbse.mode mode`;
`cp qinunitsofbvectors.ipt cksdq`;




system("$ENV{'OCEAN_BIN'}/getnval.x") == 0 or die "Failed to get nval\n";


open RUNTYPE, "cks.normal" or die;
my $runtype = <RUNTYPE>;
close RUNTYPE;
if ($runtype =~ m/true/ ) 
{
  $runtype = 1;
} elsif ($runtype =~ m/false/ )
{
  $runtype = 0;
  `mv nbuse.ipt nbuse_xas.ipt`;
  `cp nbuse_xes.ipt nbuse.ipt`;
  print "XES!\n";
} else {
  die "Failed to parse cks.normal\n";
}


#####################
#print "Running setup\n";
#system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";


`ln -s ../zWFN/u2.dat`;

my $pawrad = `cat cnbse.rad`;
chomp($pawrad);
$pawrad = sprintf("%.2f", $pawrad);

open EDGE, "hfinlist";
while (<EDGE>) {
  $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
  my $ppname = $1;
  my $znum = $2;
  my $nnum = $3;
  my $lnum = $4;
  my $elname = $5;
  my $elnum = $6;

  my $cks;
  if( $runtype ) {
    $cks = sprintf("cksc.${elname}%04u", $elnum );
  } 
  else {
    $cks = sprintf("cksv.${elname}%04u", $elnum );
  }

  print "CKS NAME = $cks\n";
  `ln -sf ../zWFN/$cks cbinf0001`;
#  `ln -sf ../zWFN/test_cks.bin cbinf0001`;
  `ln -sf cbinf0001 ufmi` == 0 or die;
#  my $zstring = sprintf("z%03u", $znum);
#  `ln -sf ../PAW/zpawinfo/*${zstring}* .` == 0 or die;
  my $zstring = sprintf("z%03i", $znum);
  print $zstring ."\n";
  `ln -sf ../PAW/zpawinfo/*${zstring}* .`;
#  `ln -sf ../PAW/zpawinfo/${zstring} .` == 0 or die;
  my $templine = `ls ../PAW/zpawinfo/*$zstring`;
  chomp($templine);
  my @pslist = split(/\s+/, $templine);
  foreach (@pslist) 
  {
    $_ =~ m/ae(\S+)/;
    `ln -sf ../PAW/zpawinfo/ae$1 .`;
    `ln -sf ae$1 ps$1`;
  }
  
  open ZNL, ">ZNL" or die;
#  `echo "$znum  $nnum  $lnum$" > ZNL`;
  print ZNL "$znum  $nnum  $lnum\n";
  close ZNL;
  $zstring = sprintf("z%03un%02ul%02u", $znum, $nnum, $lnum);
  $zstring = sprintf("z%2s%02i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
  system("cp ../SCREEN/${zstring}/zR${pawrad}/rpot ./rpotfull") == 0 
    or die "Failed to grab rpot\n../SCREEN/${zstring}/zR${pawrad}/rpot ./rpotfull\n";
#
#
  foreach my $way (@ways) {
    system("cp jtv${way} spectfile") ;#== 0 or die;
    system("$ENV{'OCEAN_BIN'}/meljtv.x");
    `mv mels jtvmels${way}`;
  }  


#  open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
#  print CKSIN "1\n$elname  $elnum  cbinf\n";
#  close CKSIN;


  foreach my $way ( @ways ) {

  `cp jtvmels${way} mels`;

  print "dotter\n";
  system("echo 'cbinf0001' | $ENV{'OCEAN_BIN'}/dotter.x") == 0 or die;


  open INFILE, ">bse.in" or die "Failed to open bse.in\n";
  my $filename = sprintf("deflinz%03un%02ul%02u", $znum, $nnum, $lnum);

  open TMPFILE, $filename or die "Failed to open $filename\n";
  my $line = <TMPFILE>;
  close TMPFILE;
  
  print INFILE $line;
  my $lookup = sprintf("%1u%1s", $nnum, $alphal{$lnum}) or die;
  my $filename = sprintf("corezetaz%03u", $znum);
  print "$lookup\t$filename\n";
  my $line = `grep $lookup $filename`;
  print INFILE $line;
  
  print INFILE "hay\n";
  open TMPFILE, "cnbse.niter" or die "Failed to open niter\n";
  <TMPFILE> =~ m/(\d+)/ or die "Failed to parse niter\n";
  my $niter = $1;
  close TMPFILE;
  my $spectrange = `cat cnbse.spect_range`;
  chomp($spectrange);
  my $gamma0 = `cat cnbse.broaden`;
  chomp($gamma0);
  
  print INFILE "$niter  $spectrange  $gamma0  0.000\n";
  close INFILE;

  if( -e "../SCREEN/core_shift.txt" )
  {
    `head -n $elnum ../core_shift.txt | tail -n 1 > core_offset `;
  } else
  {
     `rm -f core_offset`;
  }
#  system("../swbsys3.job") == 0 or die;
  system("$ENV{'OCEAN_BIN'}/cainmultip.x < bse.in >& cm.log") == 0 or die "Failed to finish\n"; 
#  my $absspct = sprintf("absspct_%2s.%u_%2s", $elname
  `mkdir -p ${zstring}_${way}/`;
  `cp {absspct,lanceigs,mulfile} ${zstring}_${way}/`;
  `cp absspct "absspct_${elname}.${elnum}_${lookup}_${way}"`;
}
}

exit 0;

