#!/usr/bin/perl

# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/OBF_cnbse_mpi\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }

my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

my @CommonFiles = ("epsilon", "xmesh.ipt", "nedges", "k0.ipt", #"nbuse.ipt", 
  "cnbse.rad", "cnbse.ways", "metal", "cksshift", "cksstretch", "cksdq", "cks.normal",
  "cnbse.niter", "cnbse.spect_range", "cnbse.broaden", "cnbse.mode", "para_prefix", "cnbse.strength" );

my @AbinitFiles = ("avecsinbohr.ipt");

my @DFTFiles = ("nelectron");

my @DenDipFiles = ("kmesh.ipt", "masterwfile", "listwfile", "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "tmels", "nelectron", "eshift.ipt" );

my @WFNFiles = ("kmesh.ipt",  "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", "avecsinbohr.ipt", "nbuse.ipt", "wvfcninfo", "wvfvainfo", "nbuse_xes.ipt", "obf_control", "ibeg.h");

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
`cp qinunitsofbvectors.ipt cksdq`;

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


# Set up mode
  my $is_xas;
  open TMPFILE, "cnbse.niter" or die "Failed to open cnbse.niter\n$!";
  <TMPFILE> =~ m/(\d+)/ or die "Failed to parse cnbse.niter";
  my $num_haydock_iterations = $1;
  close TMPFILE;
  
  open TMPFILE, "cnbse.strength" or die "Failed to open cnbse.strength\n$!";
  <TMPFILE> =~ m/([0-9]*\.?[0-9]+)/ or die "Failed to parse cnbse.strength\n";
  my $interaction_strength = $1; 
  close TMPFILE;
    
  open TMPFILE, "cnbse.mode" or die "Failed to open cnbse.mode\n";
  my $mode = <TMPFILE>;
  close TMPFILE;
  chomp($mode);
  if( lc($mode) eq 'xes' )
  {
    print "Calculating XES\n";
    $interaction_strength = 0.0;
    $is_xas = ".false.";
  } 
  elsif( lc($mode) eq 'xas' )
  {
    print "Calculating XAS\n";
    $is_xas = ".true.";
  }
  else
  {
    print "Unrecognized mode. Calculating XAS\n";
    $is_xas = ".true.";
  }

# write cks.normal file
  open TMPFILE, ">cks.normal" or die "Failed to open cks.normal for writing\n$!";
  print TMPFILE "$is_xas\n";
  close TMPFILE;

#write mode file
  open TMPFILE, ">mode" or die "Failed to open mode for writing\n$!";
  print TMPFILE "$interaction_strength    $num_haydock_iterations\n";
  close TMPFILE;





system("$ENV{'OCEAN_BIN'}/getnval.x") == 0 or die "Failed to get nval\n";


open RUNTYPE, "cks.normal" or die;
my $runtype = <RUNTYPE>;
my $run_text = '';
close RUNTYPE;
if ($runtype =~ m/true/ ) 
{
  $runtype = 1;
  $run_text = 'XAS';
} elsif ($runtype =~ m/false/ )
{
  $runtype = 0;
  `mv nbuse.ipt nbuse_xas.ipt`;
  `cp nbuse_xes.ipt nbuse.ipt`;
  $run_text = 'XES';
  print "XES!\n";
} else {
  die "Failed to parse cks.normal\n";
}


#####################
#print "Running setup\n";
#system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";

if( -e "../zWFN/u2par.dat" )
{
  `ln -s ../zWFN/u2par.dat`;
  open OUT, ">bloch_type" or die;
  print OUT "new\n";
  close OUT;
}
else
{
  `ln -s ../zWFN/u2.dat`;
}
my $pawrad = `cat cnbse.rad`;
chomp($pawrad);
$pawrad = sprintf("%.2f", $pawrad);


my %unique_z;
my %unique_znl;

open RUNLIST, ">runlist";
my $hfinlength = `wc -l hfinlist`;
chomp($hfinlength);
$hfinlength *= ($#ways + 1 );
print "$hfinlength\n";
print RUNLIST "$hfinlength\n";



open EDGE, "hfinlist";
while (<EDGE>) {
  $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
  my $ppname = $1;
  my $znum = $2;
  my $nnum = $3;
  my $lnum = $4;
  my $elname = $5;
  my $elnum = $6;

  foreach my $way (@ways) {
    print RUNLIST "$znum  $nnum  $lnum  $elname  ${nnum}$alphal{$lnum}  $elnum  $way  $run_text\n";
  }


  my $cks;
  if( $runtype ) {
    $cks = sprintf("cksc.${elname}%04u", $elnum );
  } 
  else {
    $cks = sprintf("cksv.${elname}%04u", $elnum );
  }

  print "CKS NAME = $cks\n";
  `cp ../zWFN/$cks .`;

#  my $add10_zstring = sprintf("z%03un%02ul%02u", $znum, $nnum, $lnum);
  my $zstring = sprintf("z%2s%02i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
  system("cp ../SCREEN/${zstring}/zR${pawrad}/rpot ./rpot.${zstring}") == 0 
    or die "Failed to grab rpot\n../SCREEN/${zstring}/zR${pawrad}/rpot ./rpot.${zstring}\n";
#

  $unique_z{ "$znum" } = 1;  
  $unique_znl{ "$znum $nnum $lnum" } = 1;
}
close EDGE;
close RUNLIST;

while ( my ($key, $value ) = each(%unique_z) )
{
  my $zstring = sprintf("z%03i", $key);
  print $zstring ."\n";
  `ln -sf ../PAW/zpawinfo/*${zstring}* .`;
  my $templine = `ls ../PAW/zpawinfo/*$zstring`;
  chomp($templine);
  my @pslist = split(/\s+/, $templine);
  foreach (@pslist) 
  {
    $_ =~ m/ae(\S+)/;
    `ln -sf ../PAW/zpawinfo/ae$1 .`;
    `ln -sf ae$1 ps$1`;
  }
}

while ( my ($key, $value ) = each(%unique_znl) )
{ 
  $key =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die;
  my $znum = $1;
  my $nnum = $2;
  my $lnum = $3;

  open ZNL, ">ZNL" or die;
  print ZNL "$znum  $nnum  $lnum\n";
  close ZNL;

#
  foreach my $way (@ways) {
    system("cp jtv${way} spectfile") ;#== 0 or die;
    system("$ENV{'OCEAN_BIN'}/meljtv.x");
#    `mv mels jtvmels${way}`;
    my $mel_targ = sprintf("mels.z%03un%02ul%02up%02u", $znum, $nnum, $lnum, $way );
    `mv mels $mel_targ`;
  }  


}

#  foreach my $way ( @ways ) {
#
#  `cp jtvmels${way} mels`;
#
#  print "dotter\n";
#  system("echo 'cbinf0001' | $ENV{'OCEAN_BIN'}/dotter.x") == 0 or die;
#  open EDGE, "<hfinlist";
#  my $tmp_line = <EDGE>;
#  $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
#  my $ppname = $1;
#  my $znum = $2;
#  my $nnum = $3;
#  my $lnum = $4;
#  my $elname = $5;
#  my $elnum = $6;
  my $ZNL = `cat ZNL`;
  $ZNL =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die;
  my $znum = $1;
  my $nnum = $2;
  my $lnum = $3;


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
    `cp ../SCREEN//core_shift.txt .`;
    `head -n 1 ../SCREEN/core_shift.txt  > core_offset`;
  } else
  {
     `rm -f core_offset`;
    `rm -f core_shift.txt`;
  }

  $ENV{"OMP_NUM_THREADS"}=1;

#  print "time mpirun -n 64 $ENV{'OCEAN_BIN'}/ocean.x >& cm.log";
#  system("time mpirun -n 64 $ENV{'OCEAN_BIN'}/ocean.x >& cm.log") == 0 or die "Failed to finish\n"; 


  print "time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > cm.log";
  system("time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > cm.log") == 0 or die "Failed to finish\n"; 

exit 0;

