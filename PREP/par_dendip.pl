#!/usr/bin/perl

use strict;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/par_dendip\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################
my $RunDenDip = 0;

my $stat = 0;
$stat = 1 if (-e "done");
`rm -f done`;
my $oldden = 0;
$oldden = 1 if (-e "../ABINIT/old");


my @AbFiles = ( "nkpt", "paw.nkpt",
                "avecsinbohr.ipt", "qinunitsofbvectors.ipt", "nbands", "paw.nbands", "ndtset", "gs.out");

my @CommonFiles = ( "occopt", "xmesh.ipt", "natoms", "fband", "metal" )

foreach (@AbFiles) {
  system("cp ../ABINIT/$_ .") == 0 or die "Failed to copy $_\n$!\n";
}

foreach (@CommonFile) {
  system("cp ../Common/$_ .") == 0 or die "Failed to copy $_\n$!\n";
}


my $occopt = `cat occopt`;
chomp($occopt);
my $metal; = `cat metal`;
if( $occopt == 1 ){
  if( $metal eq ".true." ) 
  {
    print "Metal set to true, occopt set to insulator. Changing metal to false\n";
    `echo .false. > metal`;
  }
}
else
{
  if( $metal != 1 ) 
  {
    print "Metal set to false, occopt set to metallic. Changin mertal to true\n";
    `echo .true. > metal`;
  }
}


open LOG, "gs.out";
my $vb;
while (<LOG>) {
  if ($_ =~ m/mband\s+(\d+)/) {
    $vb = $1;
    last;
  }
}
close LOG;
my $natoms = `cat natoms`;
my $fband = `cat fband`;
$pawnbands = `cat paw.nbands`;
my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
$cb = 1 if ($cb < 1);



#print "$stat  $oldden\n";
unless ($stat && $oldden) {
#-e "../ABINIT/RUN0001_DS1_DEN" or die "SCx_DEN not found\n";
#`ln -sf ../ABINIT/RUN0001_DS1_DEN SCx_DEN`;
`ln -sf ../ABINIT/SCx_DEN .`;

system("echo 'SCx_DEN\n1\n6\nrhoofr\n0' | $ENV{'OCEAN_BIN'}/cut3d > cut3d.log") == 0
      or die "Failed to run cut3d\n";

`tail -n 1 rhoofr > nfft`;

system("$ENV{'OCEAN_BIN'}/nelectron.x") == 0 
  or die "Failed to run nelectron.x\n";

system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0
  or die "Failed to run bvecs.pl\n";

system("$ENV{'OCEAN_BIN'}/gvecs2.pl") == 0
  or die "Failed to run gvecs2.pl\n";
}

my $rundir;
## process paw wf files ##

open NKPT, "paw.nkpt" or die "Failed to open paw.nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse paw.nkpt\n";
my @nkpt = ($1, $2, $3);
close NKPT;
$rundir = sprintf("../ABINIT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

unless( -e "PAW/done" ) { #&& -e "${rundir}/old" ) {
`rm -r PAW` if (-e "PAW");
mkdir "PAW"; 
chdir "PAW";

open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
close NKPT;


`cp ../qinunitsofbvectors.ipt .`;
`cp ../bvecs .`;
`cp ../paw.nbands .`;
`cp ../paw.nkpt .`;
`cp paw.nkpt kmesh.ipt`;


my $bandmax = `cat paw.nbands`;
chomp($bandmax);
my $nelectron = `cat ../nelectron`;

open BRANGE, ">brange.ipt";

if( $occopt == 1 ) {
  print BRANGE "1  " . $nelectron/2 . "\n";
  print BRANGE $nelectron/2+1 . "    $bandmax\n";
}
else {
  print BRANGE "1  $vb\n$cb  $bandmax\n";
}
  
close BRANGE;

#my $Nfiles = `cat Nfiles`;
`echo 1 > Nfiles`;
#my $runfile = "../../ABINIT/RUN0001_DS2_WFK";
my $runfile = "../../ABINIT/RUN0002_WFK";
#for (my $i = 1; $i <= $Nfiles; $i++) {
#  $runfile = sprintf("../${rundir}/RUN%04u_WFK", $i );
  system("ln -s $runfile RUN0001_WFK") == 0  or die "Failed to link $runfile\n";
#}

system("$ENV{'OCEAN_BIN'}/wfconvert.x") == 0 
  or die "Failed to run wfconvert.x\n";

`touch done`;
chdir "../";
}
else {
  `touch PAW/old`;
  print  "Nothing needed for PAW wfns\n";
}

## process bse wf files ##

open NKPT, "nkpt" or die "Failed to open nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse nkpt\n";
@nkpt = ($1, $2, $3);
close NKPT;

$rundir = sprintf("../ABINIT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

unless( -e "BSE/done" ) { #&& -e "${rundir}/old" ) {
`rm -r BSE` if (-e "BSE");
mkdir "BSE";
chdir "BSE";

open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
close NKPT;

#foreach ("Nfiles", "kmesh.ipt", "brange.ipt" ) {
#  system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
#}
`cp ../qinunitsofbvectors.ipt .`;
`cp ../bvecs .`;
`cp ../nelectron .`;
#my $Nfiles = `cat Nfiles`;
`cp ../nbands .`;
`echo 1 > Nfiles`;
`cp ../nkpt kmesh.ipt`;

#open BRANGE, "brange.ipt";
#my @brange;
#<BRANGE> =~ m/(\d+)\s+(\d+)/;
#$brange[0] = $1;
#$brange[1] = $2;
#<BRANGE> =~ m/(\d+)\s+(\d+)/;
#$brange[2] = $1;
#$brange[3] = $2;
#close BRANGE;
my $bandmax = `cat nbands`;
chomp($bandmax);
my $nelectron = `cat ../nelectron`;
open BRANGE, ">brange.ipt";
print BRANGE "1  " . $nelectron/2 . "\n";
print BRANGE $nelectron/2+1 . "    $bandmax\n";
close BRANGE;

my $ndtset = `cat ../ndtset`;
chomp($ndtset);
#my $runfile = "../../ABINIT/RUN0001_DS${ndtset}_WFK";
my $runfile = "../../ABINIT/RUN0001_WFK";
#for (my $i = 1; $i <= $Nfiles; $i++) {
#  $runfile = sprintf("../${rundir}/RUN%04u_WFK", $i );
  system("ln -s $runfile RUN0001_WFK") == 0 or die "Failed to link $runfile\n";
#}

system("$ENV{'OCEAN_BIN'}/wfconvert.x") == 0 
  or die "Failed to run wfconvert.x\n";

system("$ENV{'OCEAN_BIN'}/ofermi.pl") == 0
  or die "Failed to run ofermi.pl\n";
`cp eshift.ipt ../`;
system("cp efermiinrydberg.ipt ../") == 0 
  or die "Failed to copy efermiinrydberg.ipt\n";

`touch done`;
chdir "../";
}
else {
  `touch BSE/old`;
  print "Nothing needed for bse wfns\n";
}

unless ($stat && $oldden) {
system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0
  or die "Failed to run rhoofg.x\n";
`wc -l rhoG2 > rhoofg`;
`sort -n -k 6 rhoG2 >> rhoofg`;
}

chdir "BSE/";
`cp ../avecsinbohr.ipt .`;
`cp ../../Common/xmesh.ipt .`;
print "Running setup\n";
system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0
  or die "Failed to run setup\n";

  print "conugtoux\n";
  system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
  print "orthog\n";
  system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;




`touch done`;
exit 0;

