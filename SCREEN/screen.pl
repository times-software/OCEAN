#!/usr/bin/perl

use strict;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/screen\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################


my @CommonFiles = ("znucl", "paw.hfkgrid", "paw.fill", "paw.opts", "pplist", "paw.shells", "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution", "epsilon", "k0.ipt", "ibase", "scfac" );

my @AbinitFiles = ("avecsinbohr.ipt");

my @DenDipFiles = ("rhoofg", "bvecs", "efermiinrydberg.ipt");
my @DenDipFiles2 = ( "masterwfile", "listwfile", "enkfile", "kmesh.ipt", "brange.ipt" );

my @ExtraFiles = ("specpnt", "Pquadrature" );

my $runPAW = 1;
if (-e "../DenDip/PAW/old" && -e "done" ) {
  $runPAW = 0;
  foreach (@CommonFiles) {
    if (`diff -q $_ ../Common/$_`) {
      $runPAW = 1;
      last;
    }
  }
}

if ($runPAW == 0 ) {
  print "Nothing new needed for PAW stage\n";
  exit 0;
}
`rm -f done`;


foreach (@CommonFiles) {
  `cp ../Common/$_ .` == 0 or die "Failed to get $_ from Common/\n";
}
foreach (@AbinitFiles) {
  `cp ../ABINIT/$_ .` == 0 or die "Failed to get $_ from ABINIT/\n";
}
foreach (@DenDipFiles) {
  `cp ../DenDip/$_ .` == 0 or die "Failed to get $_ from DenDip/\n";
}
foreach (@DenDipFiles2) {
  `cp ../DenDip/PAW/$_ .` == 0 or die "Failed to get $_ from DenDip/PAW/\n";
}

foreach (@ExtraFiles) {
  `cp $ENV{'OCEAN_BIN'}/$_ .` == 0 or die;
}

open LISTW, "listwfile" or die "Failed to open listwfile\n";
while (<LISTW> ) {
  m/(\d+)\s+(\S+)/ or die "Malformed listwfile\n";
  `ln -sf ../DenDip/PAW/$2 $2`;
}

#open PPLIST, "pplist" or die "Failed to open pplist\n";
#while (<PPLIST>) {
#  chomp;
#  `cp "../$_" .`;
#}

unless( -e 'clipbands' ) {
print "Creating clipbands\n";
open BRANGE, "brange.ipt" or die "Failed to open brange.ipt\n";
open CLIPS, ">clipbands" or die "Failed to open clipbands for writing\n";
<BRANGE> =~ m/(\d+)\s+\d+/ or die "Malformed brange\n";
print CLIPS "$1   ";
<BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Malformed brange\n";
print CLIPS "$2\n";
close BRANGE;
close CLIPS;
}

###################################

# Setup
###################################
print "Running PAW Setup\n";
system("$ENV{'OCEAN_BIN'}/pawsetup.x") == 0 or die "Failed to run pawsetup.x\n";

#`mkdir -p zdiag/ zpawinfo/`;
`ln -sf ../PAW/zpawinfo zpawinfo`;
###################################



# shells
##################################
open SHELLS, "paw.shells" or die "Failed to open paw.shells\n";
my $numshells = 0;
my $allshells = '';
while (<SHELLS>) {
  chomp;
  $allshells .= $_ ." ";
}
close SHELLS;
my @rads = split(/ /, $allshells);
$numshells = $#rads + 1;
open SHELLS, ">shells" or die "Failed to open shells for writing\n";
print SHELLS "$numshells\n";
print SHELLS "$allshells\n";
close SHELLS;

######################################

print "Starting xipp section\n";

system("$ENV{'OCEAN_BIN'}/avg.x") == 0 or die "Failed to run avg.x\n";

open HFINLIST, "hfinlist" or die "Failed to open hfinlist\n";
#seek (HFINLIST, 0, 0);
#open SHELLS, "paw.shells" or die "Failed to open shells\n";
#my $line = <SHELLS>;
#my @rads = split(/ /, <SHELLS>);
#close SHELLS;

my $rad;
my $edgename; 
my $hfinline; my $ppfilename; my $znucl; my $nnum; my $lnum; my $elname; my $elnum;
while ($hfinline = <HFINLIST>) {
#  print $hfinline;
  ($hfinline =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)/) or die "Malformed hfinlist\t$1 $2 $3 $4 $5 $6\n";
  $ppfilename = $1;
  $znucl = $2;
  $nnum = $3;
  $lnum = $4;
  $elname = $5;
  $elnum = $6;

  $edgename = sprintf("z%2s%02i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
  print "$edgename\n";
  `mkdir -p $edgename` == 0 or die "Failed to make dir $edgename\n";

  my $avden =  sprintf("avg%2s%02i",$elname,$elnum);
  system("cp $avden avden") == 0 or die "Failed to copy density $avden\n";


  my $edgename2 = sprintf("z%03in%02il%02i",$znucl, $nnum, $lnum);
#  `mkdir -p $edgename` == 0 or die "Failed to make dir $edgename\n";

  foreach $rad (@rads) {
    my $fullrad = sprintf("%03.2f",$rad);
      `echo "8 25 $elname $elnum" | $ENV{'OCEAN_BIN'}/mkrbfile.x`;
      `mkdir -p ${edgename}/zRXT${fullrad}`;
      `mkdir -p ${edgename}/zRXF${fullrad}`;
      chdir "$edgename";
      `ln -s -f zRXT${fullrad} zR${fullrad}`;
      chdir "../";
      `cp zpawinfo/vcxxxxx${edgename2}R${fullrad} ./tmp`;
      `wc tmp > vpert`;
      `cat tmp >> vpert`;
#      system("builder.x < builder.in") == 0 or die;
      system("$ENV{'OCEAN_BIN'}/builder.x") == 0 or die;
      `echo 24 > ipt`;
      `time $ENV{'OCEAN_BIN'}/xipps.x < ipt`;
      `mv ninduced nin`;
      `echo $fullrad > ipt`;
      `cat ibase epsilon >> ipt`;
      `time $ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
      `mv reopt rom`;
      `echo 1 3 > ipt`;
      `wc rom >> ipt`;
      `cat rom >> ipt`;
      `echo 1 4 >> ipt`;
      `wc nin >> ipt`;
      `cat nin >> ipt`;
      `echo 1 2 >> ipt`;
      `wc zpawinfo/vcxxxxx${edgename2} >> ipt`;
      `cat zpawinfo/vcxxxxx${edgename2} >> ipt`;
  
      `cp ipt ipt1`;
      `echo .false. >> ipt1`;
      `echo 0.1 100 >> ipt1`;
      `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXF${fullrad}/ropt`;
      `mv {rpot,rpothires} ${edgename}/zRXF${fullrad}/`;

      `cp ipt ipt1`;
      `echo .true. >> ipt1`;
      `wc zpawinfo/vvpseud${edgename2} >> ipt1`;
      `cat zpawinfo/vvpseud${edgename2} >> ipt1`;
      `wc zpawinfo/vvallel${edgename2} >> ipt1`;
      `cat zpawinfo/vvallel${edgename2} >> ipt1`;
      `echo 0.1 100 >> ipt1`;
      `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXT${fullrad}/ropt`;
      `mv {rpot,rpothires,rom,nin} ${edgename}/zRXT${fullrad}/`;
  }
}
close HFINLIST;

`touch done`;

exit 0;
