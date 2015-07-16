#!/usr/bin/perl

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/AbinitDriver\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }
if (! $ENV{"OCEAN_ABINIT"} ) {$ENV{"OCEAN_ABINIT"} = $ENV{"OCEAN_BIN"} . "/abinit"; }
if (! $ENV{"OCEAN_CUT3D"} ) {$ENV{"OCEAN_CU3D"} = $ENV{"OCEAN_BIN"} . "/cut3d"; }

####################################
# Executables to be run
# kgen.x
# PP
# ABINIT
# avecs?

my $RunKGen = 0;
my $RunPP = 0;
my $RunABINIT = 0;
my $pawRUN = 0;
my $bseRUN = 0;

my @GeneralFiles = ("core", "para_prefix" );

my @KgenFiles = ("nkpt", "k0.ipt", "qinunitsofbvectors.ipt", "paw.nkpt");
my @BandFiles = ("nbands", "paw.nbands");
my @AbinitFiles = ( "rscale", "rprim", "ntype", "natoms", "typat",
    "verbatim", "coord", "taulist", "ecut", "etol", "nrun", "wftol", 
    "fband", "occopt", "ngkpt", "abpad");
my @PPFiles = ("pplist", "znucl");
my @OtherFiles = ("epsilon");


foreach (@PPFiles) {
  if ( -e $_ ) {
    if ( `diff -q $_ ../Common/$_` ) {
      $RunPP = 1;
      print "$_ differs\n";
      last;
    }
  }
  else {
    $RunPP = 1;
    print "$_ not found\n";
    last;
  }
}
unless ($RunPP) {
  $RunPP = 1;
  if (open STATUS, "pp.stat" ) {
    if (<STATUS> == 1) { $RunPP = 0; }
  }
  close STATUS;
}
if ($RunPP != 0) {
  `rm -f pp.stat`;
} 

if ( $RunPP ) {
  $RunABINIT = 1;
}
else {
  foreach (@AbinitFiles) {
    if ( -e $_ ) {
      if ( `diff -q $_ ../Common/$_` ) {
        $RunABINIT = 1;
        last;
      }
    }
    else {
      $RunABINIT = 1;
      last;
    }
  }
}
unless ($RunABINIT) {
 $RunABINIT = 1;
  if (open STATUS, "abinit.stat" ) {
    if (<STATUS> == 1) { $RunABINIT = 0; }
  }
  close STATUS;
}
if ($RunABINIT) {
  print "Differences found for density run. Clearing all old data\n";
#  die;
  my @dirlisting = <*>;
  foreach my $file (@dirlisting) {
    chomp($file);
    `rm -r $file`;
  }
  $RunPP = 1;
}
else {
  `touch old`;
}



#unless ($RunABINIT || $RunPP || $RunKGen ) {
#  print "Nothing needed in ABINIT Stage\n";
#  open GOUT, ">abinitstage.stat" or die;
#  print GOUT "1";
#  close GOUT;
#  exit 0;
#}

open GOUT, ">abinitstage.stat" or die;
print GOUT "0";
close GOUT;

foreach (@GeneralFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@KgenFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@BandFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

open PARA, "para_prefix" or die "$!";
my $para_prefix = <PARA>;
chomp($para_prefix);
close PARA;

#if  ($RunKGen) {
#  foreach (@KgenFiles) {
#    system("cp ../Common/$_ .") == 0 or die;
#  } 
#`echo "1" > kgen.stat`;
#}

#if ($RunPP) {
  foreach (@PPFiles) {
    system("cp ../Common/$_ .") == 0 or die;
  } 
#`echo "1" > pp.stat`;
#}

#if ($RunABINIT) {
  foreach (@AbinitFiles, @OtherFiles) {
    system("cp ../Common/$_ .") == 0 or die;
  } 
#}

#############################################
my $ecut = `cat ecut`;
chomp($ecut);
open OUT, ">ecutRy" or die;
print OUT "$ecut Ry\n";
close OUT;

# test paw.nkpt, paw.nbands
open NKPT, "paw.nkpt" or die "Failed to open paw.nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. paw.nkpt\n";
my @pawnkpt = ($1, $2, $3);
close NKPT;
open NKPT, "nkpt" or die "Failed to open nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. nkpt\n";
my @nkpt = ($1, $2, $3);
close NKPT;
my $pawnbands;
my $nbands;
open NBANDS, "paw.nbands" or die "Failed to open paw.nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse paw.nbands\n";
$pawnbands = $1;
close NBANDS;
open NBANDS, "nbands" or die "Failed to open nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
$nbands = $1;
close NBANDS;

if ( $nkpt[0] + $nkpt[1] + $nkpt[2] == 0 ) {
  `cp nkpt paw.nkpt`;
  @pawnkpt = @nkpt;
  if ( $pawnbands == 0 ) {
    `cp nbands paw.nbands`;
    $pawnbands = $nbands;
  }
  elsif ( $nbands > $pawnbands ) {
    die "paw.nbands must be larger than nbands\b";
  }
}


# test the directory for the PAW run first
my $pawDIR = sprintf("%03u%03u%03u", $pawnkpt[0], $pawnkpt[1], $pawnkpt[2] );
if ( -d $pawDIR ) {
  chdir $pawDIR;
  if (-e "abinit.stat") {
    if ( `diff -q nbands ../paw.nbands`) {
      open NBANDS, "nbands" or die "Failed to open `pwd`/nbands\n";
      <NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
      my $tmpnbands = $1;
      close NBANDS;
      $pawRUN = 1 if ( $tmpnbands < $pawnbands);
    }
    $pawRUN = 1 if ( `diff -q k0.ipt ../k0.ipt` );
  }
  else {
    $pawRUN = 1;
  }
  chdir "../"
}
else {
  $pawRUN = 1;
}

if ($pawRUN == 1) {
  print "Need to run for PAW\n";
#  die;
  `rm -rf $pawDIR`;
  mkdir $pawDIR;
}
else {
  `touch $pawDIR/old`;
}

# test the directory for the NBSE run
my $bseDIR = sprintf("%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2] );
if ( -d $bseDIR) {
  chdir $bseDIR;
  if (-e "abinit.stat") {
    if ( `diff -q nbands ../nbands`) {
      open NBANDS, "nbands" or die "Failed to open `pwd`/nbands\n";
      <NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
      my $tmpnbands = $1;
      close NBANDS;
      $bseRUN = 1 if ( $tmpnbands < $nbands);
    }
    $bseRUN = 1 if ( `diff -q k0.ipt ../k0.ipt` );
  }
  else {
    $bseRUN = 1;
  }
  chdir "../"
}
else {
  $bseRUN = 1;
}

if ( $pawRUN == 1 && $nkpt[0] == $pawnkpt[0] && $nkpt[1] == $pawnkpt[1] && $nkpt[2] == $pawnkpt[2] ) {
  $bseRUN = 0;
}

if ($bseRUN == 1 ) {
  print "Need run for the BSE\n";
  #die;
  `rm -rf $bseDIR`;
  mkdir $bseDIR;
}
else {
  `touch $bseDIR/old`;
}

#if ($RunKGen) {
#  print "Running kgen.x\n";
#  `cp nkpt kmesh.ipt`;
#  system("$ENV{'OCEAN_BIN'}/kgen.x") == 0 or die "KGEN.X Failed\n";
#  `echo "1" > kgen.stat`;
#}

if ($RunPP) {
  system("$ENV{'OCEAN_BIN'}/pp.pl znucl pplist finalpplist") == 0
    or die "Failed to run pp.pl\n";
  `echo "1" > pp.stat`;
}

#############################################
# Run abinit
# Determin what type of run is being done, par or seq

my $AbinitType = "seq";

if ($RunABINIT) {
  `echo symmorphi 0 > abfile`;
  `echo chksymbreak 0 >> abfile`;
  `echo 'acell ' >> abfile`;
  `cat rscale >> abfile`;
  `echo rprim >> abfile`;
  `cat rprim >> abfile`;
  `echo 'ntypat ' >> abfile`;
  `cat ntype >> abfile`;
  `echo 'znucl ' >> abfile`;
  `cat znucl >> abfile`;
  `echo 'natom ' >> abfile`;
  `cat natoms >> abfile`;
  `echo 'typat ' >> abfile`;
  `cat typat >> abfile`;
  `cat coord >> abfile`;
  `cat taulist >> abfile`;
  `echo 'ecut ' >> abfile`;
  `cat ecutRy >> abfile`;
  `echo 'nstep ' >> abfile`;
  `cat nrun >> abfile`;
  `echo 'diemac ' >> abfile`;
  `cat epsilon >> abfile`;
  `cat verbatim >> abfile`;
  `echo 'occopt ' >> abfile`;
  `cat occopt >> abfile`;
  `echo 'npfft 1' >> abfile`;

#if ($AbinitType eq "par" ) {
#  die "not an option\n";
#}
#########################
# seq Abinit run
#else {

### Clean ####
  `rm -f density.out`;
  `rm -f SCx_DEN`;
  `rm -f SCx_EIG`;
  `rm -f SCx_WFK`;
#  `rm -f Run.????.out`;
### Done cleaning ###


  open FILES, ">denout.files";
  print FILES "inai.denout\n"
            . "density.out\n"
            . "SC\n"
            . "SCx\n"
            . "Scxx\n";
  close FILES;
  `cat finalpplist >> denout.files`;
  
  `cat abfile > inai.denout`;
  `echo 'fband ' >> inai.denout`;
  `cat fband >> inai.denout`;
  `echo prtden 1 >> inai.denout`;
  `echo kptopt 1 >> inai.denout`;
  `echo 'ngkpt ' >> inai.denout`;
  `cat ngkpt >> inai.denout`;
  `echo 'toldfe ' >> inai.denout`;
  `cat etol >> inai.denout`;
#  `echo prtdos 3 >> inai.denout`;
#  `echo prtdosm 1 >> inai.denout`;


  print "Self-Consisten Density Run\n";
  system("$para_prefix $ENV{'OCEAN_ABINIT'} < denout.files > density.log 2> density.err") == 0
    or die "Failed to run initial density stage\n$para_prefix $ENV{'OCEAN_ABINIT'}\n";
  `echo 1 > den.stat`;

  `ln -s SCx_DEN SCx_DS0_DEN`;

  open CUTIN, ">cut3d.in" or die "Failed to open cut3d.in for writing.\n$!\n";
  print CUTIN "SCx_DEN\n1\n6\nrhoofr\n0\n";
  close CUTIN;

  system("$ENV{'OCEAN_BIN'}/cut3d < cut3d.in > cut3d.log 2> cut3d.err") == 0
      or die "Failed to run cut3d\n";

  
  `echo "1" > abinit.stat`;

}

open LOG, "density.log";
my $vb;
while (<LOG>) {
  if ($_ =~ m/nband\s+(\d+)/) {
    $vb = $1;
    last;
  } 
}   
close LOG;


if ( $pawRUN ) {
  print "PAW run\n";
  chdir $pawDIR;   
  `cp ../abfile .`;
 # copy all files over
  foreach ( @GeneralFiles, @AbinitFiles, @PPFiles, @OtherFiles) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  foreach ( "paw.nkpt", "paw.nbands", "k0.ipt", "qinunitsofbvectors.ipt", "finalpplist" ) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  `cp paw.nkpt nkpt`;
  `cp paw.nbands nbands`;
 # run KGEN
  print "Running kgen2.x\n";
  `cp paw.nkpt kmesh.ipt`;
  `echo 0.0 0.0 0.0 > qinunitsofbvectors.ipt`;
  system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";
  `echo "1" > kgen.stat`;

 
  open ABPAD, "abpad" or die;
  my $abpad = <ABPAD>;
  close ABPAD;
  $pawnbands += $abpad; 
  
  `echo "nband $pawnbands" >> abfile`;
  `echo "nbdbuf $abpad" >> abfile`;
  `echo 'iscf -2' >> abfile`;
  `echo 'tolwfr ' >> abfile`;
  `cat wftol >> abfile`;
  `echo getden -1 >> abfile`;
  `echo kptopt 0 >> abfile`;
  `echo "istwfk *1" >> abfile`;
  
  open NRUNS, "Nfiles" or die;
  my $nfiles = <NRUNS>;
  close NRUNS;
  for (my $runcount = 1; $runcount <= $nfiles; $runcount++ ) 
  {
    my $abfilename = sprintf("inabinit.%04i", $runcount );
    `cp abfile "$abfilename"`;
    my $kptfile =  sprintf("kpts.%04i", $runcount);
    `cat $kptfile >> $abfilename`;

    my $deninFiles = sprintf("denin.files.%04i", $runcount );
    open FILES, ">$deninFiles";
    my $Runout = sprintf('Run.%04i.out', $runcount );
    my $RUN = sprintf('RUN%04i', $runcount);
    print FILES "$abfilename\n"
              . "$Runout\n"
              . "../SCx\n"
              . "$RUN\n"
              . "SCxx" . $runcount . "\n";
    close FILES;
    `cat finalpplist >> $deninFiles`;

    my $denin = sprintf("denin.files.%04i", $runcount);
    print "$para_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log";
    system("$para_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log 2> ABINIT.$runcount.err") == 0 or
      die "$!\n";
    print "\n";
  }



  my $natoms = `cat natoms`;
  my $fband = `cat fband`;
  $pawnbands = `cat paw.nbands`;
  my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
  $cb = 1 if ($cb < 1);
  open BRANGE, ">brange.ipt" or die;
  print BRANGE "1  $vb\n"
             . "$cb $pawnbands\n";
  close BRANGE;

  open STATUS, ">abinit.stat" or die;
  print STATUS "1";
  close STATUS;

####
  chdir "../";
}

if ( $bseRUN ) {
  chdir $bseDIR;
  print "BSE run\n";
  `cp ../abfile .`;
 # copy all files over
  foreach ( @GeneralFiles, @AbinitFiles, @PPFiles, @OtherFiles) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  foreach ( "nkpt", "nbands", "k0.ipt", "qinunitsofbvectors.ipt", "finalpplist" ) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
 # run KGEN
  print "Running kgen2.x\n";
  `cp nkpt kmesh.ipt`;
  system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";
  `echo "1" > kgen.stat`;

 
  open ABPAD, "abpad" or die;
  my $abpad = <ABPAD>;
  close ABPAD;
  $nbands += $abpad; 
  
  `echo "nband $nbands" >> abfile`;
  `echo "nbdbuf $abpad" >> abfile`;
  `echo 'iscf -2' >> abfile`;
  `echo 'tolwfr ' >> abfile`;
  `cat wftol >> abfile`;
  `echo getden -1 >> abfile`;
  `echo kptopt 0 >> abfile`;
  `echo "istwfk *1" >> abfile`;

  open NRUNS, "Nfiles" or die;
  my $nfiles = <NRUNS>;
  close NRUNS;
  for (my $runcount = 1; $runcount <= $nfiles; $runcount++ )
  {
    my $abfilename = sprintf("inabinit.%04i", $runcount );
    `cp abfile "$abfilename"`;
    my $kptfile =  sprintf("kpts.%04i", $runcount);
    `cat $kptfile >> $abfilename`;
  
    my $deninFiles = sprintf("denin.files.%04i", $runcount );
    open FILES, ">$deninFiles";
    my $Runout = sprintf('Run.%04i.out', $runcount );
    my $RUN = sprintf('RUN%04i', $runcount);
    print FILES "$abfilename\n"
              . "$Runout\n"
              . "../SCx\n"
              . "$RUN\n"
              . "SCxx" . $runcount . "\n";
    close FILES;
    `cat finalpplist >> $deninFiles`;

    my $denin = sprintf("denin.files.%04i", $runcount);
    print "$para_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log";
    system("$para_prefix $ENV{'OCEAN_ABINIT'} < $denin > ABINIT.$runcount.log 2> ABINIT.$runcount.err") == 0 or
      die "$!\n";
    print "\n";
  }



  my $natoms = `cat natoms`;
  my $fband = `cat fband`;
  my $bands = `cat nbands`;
  my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
  $cb = 1 if ($cb < 1);
  open BRANGE, ">brange.ipt" or die;
  print BRANGE "1  $vb\n"
             . "$cb $bands\n";
  close BRANGE;


  
  open STATUS, ">abinit.stat" or die;
  print STATUS "1";
  close STATUS;
####
  chdir "../";
}

print `pwd`;
open RSCALE, "rscale" or die;
open RPRIM, "rprim" or die;
<RSCALE> =~  m/(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)/ or die;
my @rscale = ($1, $3, $5);
print "$1\t$3\t$5\n";
close RSCALE;

open AVECS, ">avecsinbohr.ipt" or die;
for (my $i = 0; $i < 3; $i++ ) {
  <RPRIM> =~  m/([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)/ or die "$_";
  print AVECS $1*$rscale[0] . "  " . $3*$rscale[1] .  "  " . $5*$rscale[2] . "\n";
  print "$1\t$3\t$5\n";
}
close RPRIM;
close AVECS;


print "Abinit stage complete\n";


exit 0;

