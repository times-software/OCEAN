#!/usr/bin/perl

use strict;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/dendip\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################
my $RunDenDip = 0;

#open STATUS, "../DFT/abinitstage.stat" or die;
#<STATUS> =~ m/(\d)/;
#if ($1 == 1) {
#  print "Wavefunctions unchanged. Nothing to be done for DenDip\n";
#  exit 0;
#}

my $stat = 0;
$stat = 1 if (-e "done");
`rm -f done`;
my $oldden = 0;
$oldden = 1 if (-e "../DFT/old");


my @AbFiles = ( "rhoofr", "density.out", "nkpt", "paw.nkpt",
                "avecsinbohr.ipt", "qinunitsofbvectors.ipt");

foreach (@AbFiles) {
  system("cp ../DFT/$_ .") == 0 or die "Failed to copy $_\n";
}
print "$stat  $oldden\n";
unless ($stat && $oldden) {
-e "../DFT/SCx_DEN" or die "SCx_DEN not found\n";
`ln -sf ../DFT/SCx_DEN SCx_DEN`;

#`ln -sf ../DFT/RUN????_WFK .`;

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
$rundir = sprintf("../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

unless( -e "PAW/done" && -e "${rundir}/old" ) {
`rm -r PAW` if (-e "PAW");
mkdir "PAW"; 
chdir "PAW";

open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
close NKPT;


foreach ("Nfiles", "kmesh.ipt", "brange.ipt", "qinunitsofbvectors.ipt" ) {
  system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
}
#`cp ../qinunitsofbvectors.ipt .`;
`cp ../bvecs .`;
`cp ../${rundir}/umklapp .`;
open BRANGE, "brange.ipt";
my @brange;
<BRANGE> =~ m/(\d+)\s+(\d+)/;
$brange[0] = $1;
$brange[1] = $2;
<BRANGE> =~ m/(\d+)\s+(\d+)/;
$brange[2] = $1;
$brange[3] = $2;
close BRANGE;
my $nelectron = `cat ../nelectron`;
open BRANGE, ">brange.ipt";
print BRANGE "1  " . $nelectron/2 . "\n";
print BRANGE $nelectron/2+1 . "    $brange[3]\n";
close BRANGE;

my $Nfiles = `cat Nfiles`;
my $runfile;
for (my $i = 1; $i <= $Nfiles; $i++) {
  $runfile = sprintf("../${rundir}/RUN%04u_WFK", $i );
  system("ln -s $runfile .") == 0  or die "Failed to link $runfile\n";
}

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

$rundir = sprintf("../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

unless( -e "BSE/done" && -e "${rundir}/old" ) {
`rm -r BSE` if (-e "BSE");
mkdir "BSE";
chdir "BSE";

open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
close NKPT;

foreach ("Nfiles", "kmesh.ipt", "brange.ipt") {
  system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
}
`cp ../qinunitsofbvectors.ipt .`;
`cp ../bvecs .`;
`cp ../nelectron .`;
`cp ../${rundir}/umklapp .`;
my $Nfiles = `cat Nfiles`;

open BRANGE, "brange.ipt";
my @brange;
<BRANGE> =~ m/(\d+)\s+(\d+)/;
$brange[0] = $1;
$brange[1] = $2;
<BRANGE> =~ m/(\d+)\s+(\d+)/;
$brange[2] = $1;
$brange[3] = $2;
close BRANGE;
my $nelectron = `cat ../nelectron`;
open BRANGE, ">brange.ipt";
print BRANGE "1  " . $nelectron/2 . "\n";
print BRANGE $nelectron/2+1 . "    $brange[3]\n";
close BRANGE;


my $runfile;
for (my $i = 1; $i <= $Nfiles; $i++) {
  $runfile = sprintf("../${rundir}/RUN%04u_WFK", $i );
  system("ln -s $runfile .") == 0 or die "Failed to link $runfile\n";
}

system("$ENV{'OCEAN_BIN'}/wfconvert.x") == 0 
  or die "Failed to run wfconvert.x\n";

system("$ENV{'OCEAN_BIN'}/ofermi.pl") == 0
  or die "Failed to run ofermi.pl\n";
`cp eshift.ipt ../`;
system("cp efermiinrydberg.ipt ../") == 0 
  or die "Failed to copy efermiinrydberg.ipt\n";


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

`touch done`;
exit 0;

