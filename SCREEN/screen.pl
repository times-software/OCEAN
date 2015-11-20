#!/usr/bin/perl
# Copyright (C) 2014 OCEAN collaboration
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
  $0 =~ m/(.*)\/screen\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################


my @CommonFiles = ("znucl", "paw.hfkgrid", "paw.fill", "paw.opts", "pplist", "paw.shells", "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution", "epsilon", "k0.ipt", "ibase", "scfac", "core_offset", "dft", "avecsinbohr.ipt", "para_prefix" );

my @DenDipFiles = ("rhoofg", "bvecs", "efermiinrydberg.ipt");
my @DenDipFiles2 = ( "masterwfile", "listwfile", "enkfile", "kmesh.ipt", "brange.ipt" );

my @ExtraFiles = ("specpnt", "Pquadrature" );

my $runPAW = 1;
if (-e "../PREP/PAW/old" && -e "done" ) {
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
  copy( "../Common/$_", $_ ) or die "Failed to get $_ from Common/\n$!";
}
#foreach (@DFTFiles) {
#  copy( "../DFT/$_", $_ ) or die "Failed to get $_ from DFT/\n$!";
#}
foreach (@DenDipFiles) {
  copy( "../PREP/$_", $_ ) or die "Failed to get $_ from PREP/\n$!";
}
foreach (@DenDipFiles2) {
  copy( "../PREP/PAW/$_", $_ ) or die "Failed to get $_ from PREP/PAW/\n$!";
}

foreach (@ExtraFiles) {
  copy( "$ENV{'OCEAN_BIN'}/$_", $_ ) or die "Failed to get $_ from $ENV{'OCEAN_BIN'}/\n$!";
}

open LISTW, "listwfile" or die "Failed to open listwfile\n";
while (<LISTW> ) {
  m/(\d+)\s+(\S+)/ or die "Malformed listwfile\n";
  `ln -sf ../PREP/PAW/$2 $2`;
}

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
  copy( $avden,  "avden" ) or die "Failed to copy density $avden\n$!";


  my $edgename2 = sprintf("z%03in%02il%02i",$znucl, $nnum, $lnum);
#  `mkdir -p $edgename` == 0 or die "Failed to make dir $edgename\n";

  foreach $rad (@rads) {
    my $fullrad = sprintf("%03.2f",$rad);
      `echo "8 25 $elname $elnum" | $ENV{'OCEAN_BIN'}/mkrbfile.x`;
      `mkdir -p ${edgename}/zRXT${fullrad}`;
      `mkdir -p ${edgename}/zRXF${fullrad}`;
      `mkdir -p ${edgename}/zRXS${fullrad}`;
      chdir "$edgename";
      `ln -s -f zRXT${fullrad} zR${fullrad}`;
      chdir "../";
      copy( "zpawinfo/vcxxxxx${edgename2}R${fullrad}", "tmp" ) 
          or die "Failed to copy zpawinfo/vcxxxxx${edgename2}R${fullrad}\n$!";
      `wc tmp > vpert`;
      `cat tmp >> vpert`;
      system("$ENV{'OCEAN_BIN'}/builder.x") == 0 or die;
      copy( "ximat", "${edgename}/zR${fullrad}/ximat" ) or die "Failed to copy ximat to ${edgename}/\n$!";
      copy( "ximat_small", "${edgename}/zR${fullrad}/ximat_small" ) or die "Failed to copy ximat to ${edgename}/\n$!";

      `echo 24 > ipt`;
      `time $ENV{'OCEAN_BIN'}/xipps.x < ipt`;
      move( "ninduced", "nin" ) or die "Failed to move ninduced.\n$!";
      `echo $fullrad > ipt`;
      `cat ibase epsilon >> ipt`;
      `time $ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
      move( "reopt", "rom" ) or die "Failed to move reopt.\n$!";
      `echo 1 3 > ipt`;
      `wc rom >> ipt`;
      `cat rom >> ipt`;
      `echo 1 4 >> ipt`;
      `wc nin >> ipt`;
      `cat nin >> ipt`;
      `echo 1 2 >> ipt`;
      `wc zpawinfo/vcxxxxx${edgename2} >> ipt`;
      `cat zpawinfo/vcxxxxx${edgename2} >> ipt`;
  
			# Full ximat, but using false/older style in rscombine
      copy( "ipt", "ipt1" ) or die;
      `echo .false. >> ipt1`;
      `echo 0.1 100 >> ipt1`;
      `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXF${fullrad}/ropt`;
      move( "rpot", "$edgename/zRXF$fullrad/" ) or die "Failed to move rpot\n$!";
      move( "rpothires", "$edgename/zRXF$fullrad/" ) or die "Failed to move rpothires\n$!";

			# Full ximat and most up-to-date rscombine settings
      copy( "ipt", "ipt1" ) or die;
      `echo .true. >> ipt1`;
      `wc zpawinfo/vvpseud${edgename2} >> ipt1`;
      `cat zpawinfo/vvpseud${edgename2} >> ipt1`;
      `wc zpawinfo/vvallel${edgename2} >> ipt1`;
      `cat zpawinfo/vvallel${edgename2} >> ipt1`;
      `echo 0.1 100 >> ipt1`;
      `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXT${fullrad}/ropt`;
      move( "rpot", "$edgename/zRXT$fullrad/" ) or die "Failed to move rpot\n$!";
      move( "rpothires", "$edgename/zRXT$fullrad/" ) or die "Failed to move rpothires\n$!";
      move( "rom", "$edgename/zRXT$fullrad/" ) or die "Failed to move rom\n$!";
      move( "nin", "$edgename/zRXT$fullrad/" ) or die "Failed to move nin\n$!";

      #####################
      # Now use ximat_small
			move( "ximat_small", "ximat" ) or die "$!";
      `echo 24 > ipt`;
      `time $ENV{'OCEAN_BIN'}/xipps.x < ipt`;
      move( "ninduced", "nin" ) or die "Failed to move ninduced.\n$!";
      `echo $fullrad > ipt`;
      `cat ibase epsilon >> ipt`;
      `time $ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
      move( "reopt", "rom" ) or die "Failed to move reopt.\n$!";
      `echo 1 3 > ipt`;
      `wc rom >> ipt`;
      `cat rom >> ipt`;
      `echo 1 4 >> ipt`;
      `wc nin >> ipt`;
      `cat nin >> ipt`;
      `echo 1 2 >> ipt`;
      `wc zpawinfo/vcxxxxx${edgename2} >> ipt`;
      `cat zpawinfo/vcxxxxx${edgename2} >> ipt`;

      copy( "ipt", "ipt1" ) or die;
      `echo .true. >> ipt1`;
      `wc zpawinfo/vvpseud${edgename2} >> ipt1`;
      `cat zpawinfo/vvpseud${edgename2} >> ipt1`;
      `wc zpawinfo/vvallel${edgename2} >> ipt1`;
      `cat zpawinfo/vvallel${edgename2} >> ipt1`;
      `echo 0.1 100 >> ipt1`;
      `time $ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXS${fullrad}/ropt`;
      move( "rpot", "$edgename/zRXS$fullrad/" ) or die "Failed to move rpot\n$!";
      move( "rpothires", "$edgename/zRXS$fullrad/" ) or die "Failed to move rpothires\n$!";
      move( "rom", "$edgename/zRXS$fullrad/" ) or die "Failed to move rom\n$!";
      move( "nin", "$edgename/zRXS$fullrad/" ) or die "Failed to move nin\n$!";
      # finished with ximat small
      #######################

  }
}
close HFINLIST;

# core offsets for QE
my $dft = `cat dft`;
if( $dft =~ m/qe/i )
{
	`ln -s ../DFT/Out .`;
	my $core_offset = `cat core_offset`;
	chomp $core_offset;
	if( $core_offset =~ m/false/i )
	{
	        print "No core shift\n";
		`rm core_shift.txt` if( -e "core_shift.txt" );
	} else
	{
	        `time perl $ENV{'OCEAN_BIN'}/core_shift.pl > core_shift.log`;
	}
}

`touch done`;

exit 0;
