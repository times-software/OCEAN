#!/usr/bin/perl

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/QespressoDriver\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (! $ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }


my $RunKGen = 0;
my $RunPP = 0;
my $RunESPRESSO = 0;
my $pawRUN = 0;
my $bseRUN = 0;

my @GeneralFiles = ("core" );

my @KgenFiles = ("k0.ipt", "qinunitsofbvectors.ipt", "paw.nkpt", "nkpt");
#my @KgenFiles = ("k0.ipt", "qinunitsofbvectors.ipt", "scf.nkpt", "paw.nkpt", "bse.nkpt");
my @BandFiles = ("paw.nbands", "nbands");
#my @BandFiles = ("paw.nbands", "bse.nbands");
my @EspressoFiles = ( "coord", "degauss", "ecut", "etol", "fband", "ibrav",
    "isolated", "mixing", "natoms", "ngkpt", "noncolin", "nrun", "nscf.kshift", "ntype",
    "occopt", "occtype", "prefix", "ppdir", "rprim", "rscale", "scf.kshift", "smearing",
    "spinorb", "taulist", "typat", "verbatim", "work_dir", "wftol");
#my @EspressoFiles = ( "coord", "degauss", "ecut", "etol", "fband", "ibrav", 
#    "isolated", "mixing", "natoms", "ngkpt", "noncolin", "nrun", "nscf.kshift", "ntype", 
#    "occopt", "occtype", "prefix", "ppdir", "rprim", "rscale", "scf.kshift", "smearing", 
#    "spinorb", "taulist", "typat", "verbatim", "work_dir", "wftol");
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
  $RunESPRESSO = 1;
}
else {
  system("touch ../Common/scf.kshift");
  system("echo '0  0  0' > ../Common/scf.kshift");
  system("cp ../Common/scf.kshift ../Common/nscf.kshift");
  foreach (@EspressoFiles) {
    if ( -e $_ ) {
      if ( `diff -q $_ ../Common/$_` ) {
        $RunESPRESSO = 1;
        last;
      }
    }
    else {
      $RunESPRESSO = 1;
      last;
    }
  }
}
unless ($RunESPRESSO) {
 $RunESPRESSO = 1;
  if (open STATUS, "espresso.stat") {
    if (<STATUS> == 1) { $RunESPRESSO = 0; }
  }
  close STATUS;
}
if ($RunESPRESSO) {
  print "Differences found for density run. Clearing all old data\n";
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


open GOUT, ">espressostage.stat" or die;
print GOUT "0";
close GOUT;

foreach (@GeneralFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@KgenFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
system("cp paw.nkpt scf.nkpt");
system("mv nkpt bse.nkpt");
foreach (@BandFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
system("mv nbands bse.nbands");
foreach (@PPFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
system("echo '0  0  0' > ../Common/scf.kshift");
system("cp ../Common/scf.kshift ../Common/nscf.kshift");
foreach (@EspressoFiles, @OtherFiles) {
#  system("cp ../Common/$_ .") == 0 or die;
  system("cp ../Common/$_ .");
} 

#############################################

# test paw.nkpt, paw.nbands
open NKPT, "paw.nkpt" or die "Failed to open paw.nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. paw.nkpt\n";
my @pawnkpt = ($1, $2, $3);
close NKPT;
open NKPT, "bse.nkpt" or die "Failed to open bse.nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. bse.nkpt\n";
my @bsenkpt = ($1, $2, $3);
close NKPT;
my $pawnbands;
my $bsenbands;
open NBANDS, "paw.nbands" or die "Failed to open paw.nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse paw.nbands\n";
$pawnbands = $1;
close NBANDS;
open NBANDS, "bse.nbands" or die "Failed to open nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse bse.nbands\n";
$bsenbands = $1;
close NBANDS;

my $pawnkpts = ( $pawnkpt[0] * $pawnkpt[1] * $pawnkpt[2] );
`echo pawnkpts`;
`echo $pawnkpts`;

if ( $bsenkpt[0] + $bsenkpt[1] + $bsenkpt[2] == 0 ) {
  `cp bse.nkpt paw.nkpt`;
  @pawnkpt = @bsenkpt;
  if ( $pawnbands == 0 ) {
    `cp bse.nbands paw.nbands`;
    $pawnbands = $bsenbands;
  }
  elsif ( $bsenbands > $pawnbands ) {
    die "paw.nbands must be larger than bse.nbands\b";
  }
}


# test the directory for the PAW run first
my $pawDIR = sprintf("%03u%03u%03u", $pawnkpt[0], $pawnkpt[1], $pawnkpt[2] );
if ( -d $pawDIR ) {
  chdir $pawDIR;
  if (-e "espresso.stat") {
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
  `rm -rf $pawDIR`;
  mkdir $pawDIR;
}
else {
  `touch $pawDIR/old`;
}

# test the directory for the NBSE run
my $bseDIR = sprintf("%03u%03u%03u", $bsenkpt[0], $bsenkpt[1], $bsenkpt[2] );
if ( -d $bseDIR) {
  chdir $bseDIR;
  if (-e "espresso.stat") {
    if ( `diff -q nbands ../bse.nbands`) {
      open NBANDS, "nbands" or die "Failed to open `pwd`/nbands\n";
      <NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
      my $tmpnbands = $1;
      close NBANDS;
      $bseRUN = 1 if ( $tmpnbands < $bsenbands);
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

my $bsenkpts = ( $bsenkpt[0] * $bsenkpt[1] * $bsenkpt[2] );
`echo bsenkpts`;
`echo $bsenkpts`;

if ( $pawRUN == 1 && $bsenkpt[0] == $pawnkpt[0] && $bsenkpt[1] == $pawnkpt[1] && $bsenkpt[2] == $pawnkpt[2] ) {
  $bseRUN = 0;
}

if ($bseRUN == 1 ) {
  print "Need run for the BSE\n";
  `rm -rf $bseDIR`;
  mkdir $bseDIR;
}
else {
  `touch $bseDIR/old`;
}


if ($RunPP) {
  system("$ENV{'OCEAN_BIN'}/pp.pl znucl pplist finalpplist") == 0
    or die "Failed to run pp.pl\n";
  `echo "1" > pp.stat`;
}

#############################################
# Run espresso
# Determine what type of run is being done, par or seq
# > make this unecessary

# remove this switch
my $EspressoType = "seq"; 



# make additional files for QE input card

print "making the coordinates\n";
system("$ENV{'OCEAN_BIN'}/makecoords.x") == 0
    or die "Failed to make coordinates\n";

print "making acell\n";
system("$ENV{'OCEAN_BIN'}/makeacell.x") == 0
    or die "Failed to make acell\n";

print "making atompp\n";
system("$ENV{'OCEAN_BIN'}/makeatompp.x") == 0
    or die "Failed to make acell\n";

my $kshift = " 0 0 0\n";

open(KSHIFT, '>kshift') or die "couldn't open kshift\n";
print KSHIFT $kshift;
close(KSHIFT);




if ($RunESPRESSO) {


 # collect information about the unit cell
 # 
 my $ibrav = 0;
 open(IBRAV, 'ibrav') or die "couldn't open ibrav";
 $ibrav = <IBRAV>;
 close(IBRAV);

 my $line = "";
 my $celldm1 = 0;
 my $celldm2 = 0;
 my $celldm3 = 0;
 open(RSCALE, 'rscale') or die "couldn't open rscale";
 foreach $line (<RSCALE>) {
   ($celldm1, $celldm2, $celldm3) = split(' ' ,$line);
 }
 close(RSCALE);

# `echo "We are in "`;
# `pwd`;
 system("python $ENV{'OCEAN_BIN'}/celldm.py") == 0
     or die "Celldm failed\n";

 ### write SCF input card for initial density
 #   need to move this to a separate script for tidiness
 #
 `echo "&control" > qefile`;
 `echo "   calculation = 'scf'" >> qefile`;
 `echo -n "   prefix = '" >> qefile`;
 `head -c -1 prefix >> qefile`;
 `echo "'" >> qefile`;
 `echo -n "   pseudo_dir = " >> qefile`;
 `cat ppdir >> qefile`;
 `echo -n "   outdir = " >> qefile`;
 `cat work_dir >> qefile`;
 `echo "   tstress = .true." >> qefile`;
 `echo "   tprnfor = .true." >> qefile`;
 `echo "   wf_collect = .true." >> qefile`;
 `echo "/" >> qefile`;

 `echo "&system" >> qefile`;
 `echo -n "   ibrav = " >> qefile`;    
 `cat ibrav >> qefile`;                
 if ($ibrav != 0) {
  `cat cellinfo >> qefile`;
#   `echo -n "   celldm(1) = " >> qefile`;
#   `echo ${celldm1} >> qefile`;
#   `echo -n "   celldm(2) = " >> qefile`;
#   `echo ${celldm2} >> qefile`;
#   `echo -n "   celldm(3) = " >> qefile`;
#   `echo ${celldm3} >> qefile`;
#   `echo -n "   celldm(4) = " >> qefile`;
#   `echo ${celldm4} >> qefile`;
#   `echo -n "   celldm(5) = " >> qefile`;
#   `echo ${celldm5} >> qefile`;
#   `echo -n "   celldm(6) = " >> qefile`;
#   `echo ${celldm6} >> qefile`;
 }
# `echo -n "   celldm(2) = " >> qefile`;
# `echo ${celldm2} >> qefile`;
# `echo -n "   celldm(3) = " >> qefile`;
# `echo ${celldm3} >> qefile`;
 `echo -n "   nat = " >> qefile`;
 `cat natoms >> qefile`;
 `echo -n "   ntyp = " >> qefile`;
 `cat ntype >> qefile`;
# `echo -n "   nbnd = " >> qefile`;
# `cat paw.nbands >> qefile`;
 `echo -n "   noncolin = " >> qefile`;
 `cat noncolin >> qefile`;
 `echo -n "   lspinorb = " >> qefile`;
 `cat spinorb >> qefile`;
 `echo -n "   ecutwfc = " >> qefile`;
 `cat ecut >> qefile`;
 #`echo "   ecutrho = " >> qefile`;
 #`cat ecutrho >> qefile`;
 `echo -n "   occupations = " >> qefile`;
 `cat occtype >> qefile`;
 `echo -n "   smearing = " >> qefile`;
 `cat smearing >> qefile`;
 `echo -n "   degauss = " >> qefile`;
 `cat degauss >> qefile`;
 `echo "   nosym = .true." >> qefile`;
 `echo "   noinv = .true." >> qefile`;
# `echo -n "   assume_isolated = " >> qefile`;
# `cat isolated >> qefile`;
 `echo "/" >> qefile`;

 `echo "&electrons" >> qefile`;
 `echo -n "   conv_thr    = " >> qefile`;
 `cat etol >> qefile`;
 `echo -n "   mixing_beta = " >> qefile`;
 `cat mixing >> qefile`;
 `echo -n "   electron_maxstep = " >> qefile`;
 `cat nrun >> qefile`;
 `echo "/" >> qefile`;

 `echo "&ions" >> qefile`;
 `echo "/" >> qefile`;

 `echo "ATOMIC_SPECIES" >> qefile`;
 `cat atompp >> qefile`;

 if ($ibrav == 0) {
    `echo "CELL_PARAMETERS cubic" >> qefile`;
   `cat acell >> qefile`;
 }

 `echo "ATOMIC_POSITIONS crystal" >> qefile`;
 `cat coords >> qefile`;

 `echo "K_POINTS automatic" >> qefile`;
 `tr '\n' ' ' < scf.nkpt >> qefile`;
 `cat scf.kshift >> qefile`;

 # mv qefile to appropriate location
 `cp qefile scf.in`;



 ## SCF Post-Processing initialize and set defaults
 
 my $outdir = "'./Out/'\n";
 
 my $filplot = "'system.rho'\n";
 my $plotnum = "0\n";
 
 my $nfile = "1\n";
 my $filepp = "'system.rho', ";
 my $weight = "1.0\n";
 my $iflag = "3\n";
 my $outformat = "6\n";
 my $fileout = "'system.rho.dat'\n";
 
 
 ### write Post-Processing input card for density
 
 `echo "&inputpp" > ppfile`;
 `echo -n "   prefix = " >> ppfile`;
 `cat prefix >> ppfile`;
 `echo -n "   outdir = " >> ppfile`;
 `cat work_dir >> ppfile`;
 `echo "   filplot = 'system.rho'" >> ppfile`;
 `echo "   plot_num = 0" >> ppfile`;
 `echo "/" >> ppfile`;
 
 `echo "&plot" >> ppfile`;
 `echo "   nfile = 1" >> ppfile`;
 `echo "   filepp(1) = 'system.rho', weight(1) = 1" >> ppfile`;
 `echo "   iflag = 3" >> ppfile`;
 `echo "   output_format = 6" >> ppfile`;
 `echo "   fileout = 'system.rho.dat'" >> ppfile`;
 `echo "/" >> ppfile`;
 
 `mv ppfile pp.in`;
 

### get the number of cores we are going to run on
#
 my $nc=0;
 open CORES, "core";
 $nc = int(<CORES>);
 `echo $nc >> testlog`;
 close CORES;
          
# `echo NodeFile`;
# `echo $PBS_NODEFILE`;

 ### the SCF run for initial density
 ##
 print "Density SCF Run\n";
  system("mpirun -np $nc $ENV{'OCEAN_BIN'}/pw.x < scf.in > scf.out 2>&1") == 0
      or die "Failed to run scf stage for Density\n";
 `echo 1 > scf.stat`;

 print "Density PP Run\n";
  system("mpirun -np $nc $ENV{'OCEAN_BIN'}/pp.x < pp.in > pp.out 2>&1") == 0
     or die "Failed to run post-process for density stage\n";
 `echo 1 > den.stat`;

 ## convert the density file to proper format
 print "Density conversion\n";
 system("$ENV{'OCEAN_BIN'}/converter.py system.rho.dat") == 0
     or die "Failed to convert density\n";

 ## find the top of the valence bance
 system("$ENV{'OCEAN_BIN'}/qeband.x") == 0
     or die "Failed to count bands\n";

 ## obtain the number of occupied bands
 system("$ENV{'OCEAN_BIN'}/qebocc.x") == 0
     or die "Failed to count bands\n";
      

 ## get the number of occupied bands
 `grep 1.000 bands.out | wc -l > vb`;
 my $vb = 0;
 open BANDS, "vb" or die;
 $vb = <BANDS>;
 close BANDS;

 ## make brange file
 my $natoms = `cat natoms`;
 my $fband = `cat fband`;
 $pawnbands = `cat paw.nbands`;
 my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
 $cb = 1 if ($cb < 1);
 open BRANGE, ">brange.ipt" or die;
 print BRANGE "1  $vb"
            . "$cb $pawnbands";
 close BRANGE;

 open STATUS, ">espresso.stat" or die;
 print STATUS "1";
 close STATUS;

 `cp scf.out density.out`;

 ####
 # chdir "../";

} # end SCF for density
      





### Do NSCF for PAW run

if ( $pawRUN ) {
  print "PAW run\n";
  chdir $pawDIR;   
# `cp ../paw.scf.in scf.in`;
# `cp ../paw.pp.in pp.in`;
 # copy all files over
  `cp -r ../Out/ .`; # try just linking this
  `cp ../acell .`;
  `cp ../atompp .`;
  `cp ../coords .`;
  foreach ( @GeneralFiles, @EspressoFiles, @PPFiles, @OtherFiles) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  foreach ( "paw.nkpt", "nscf.kshift", "paw.nbands", "k0.ipt", "qinunitsofbvectors.ipt", "finalpplist", "cellinfo" ) {
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


 my $ibrav = 0;
 open(IBRAV, 'ibrav') or die "couldn't open ibrav";
 $ibrav = <IBRAV>;
 close(IBRAV);


 my $line = "";
 my $celldm1 = 0;
 my $celldm2 = 0;
 my $celldm3 = 0;
 open(RSCALE, 'rscale') or die "couldn't open rscale";
 foreach $line (<RSCALE>) {
   ($celldm1, $celldm2, $celldm3) = split(' ' ,$line);
 }
 close(RSCALE);


### PAW NSCF write input card
#   push this into a separate script for cleanliness
#    
 `echo "&control" > qefile`;
 `echo "   calculation = 'nscf'" >> qefile`;
 `echo -n "   prefix = '" >> qefile`;
 `head -c -1 prefix >> qefile`;
 `echo "'" >> qefile`;
 `echo -n "   pseudo_dir = " >> qefile`;
 `cat ppdir >> qefile`;
 `echo -n "   outdir = " >> qefile`;
 `cat work_dir >> qefile`;
 `echo "   tstress = .true." >> qefile`;
 `echo "   tprnfor = .true." >> qefile`;
 `echo "   wf_collect = .true." >> qefile`;
 `echo "/" >> qefile`;

 `echo "&system" >> qefile`;
 `echo -n "   ibrav = " >> qefile`;
 `cat ibrav >> qefile`;
 if ($ibrav != 0) {
   `cat cellinfo >> qefile`;
#   `echo -n "   celldm(1) = " >> qefile`;
#   `echo ${celldm1} >> qefile`;
#   `echo -n "   celldm(2) = " >> qefile`;
#   `echo ${celldm2} >> qefile`;
#   `echo -n "   celldm(3) = " >> qefile`;
#   `echo ${celldm3} >> qefile`;
#   `echo -n "   celldm(4) = " >> qefile`;
#   `echo ${celldm4} >> qefile`;
#   `echo -n "   celldm(5) = " >> qefile`;
#   `echo ${celldm5} >> qefile`;
#   `echo -n "   celldm(6) = " >> qefile`;
#   `echo ${celldm6} >> qefile`;
 }
# `echo -n "   celldm(2) = " >> qefile`;
# `echo ${celldm2} >> qefile`;
# `echo -n "   celldm(3) = " >> qefile`;
# `echo ${celldm3} >> qefile`;
 `echo -n "   nat = " >> qefile`;
 `cat natoms >> qefile`;
 `echo -n "   ntyp = " >> qefile`;
 `cat ntype >> qefile`;
 `echo -n "   nbnd = " >> qefile`;
 `cat paw.nbands >> qefile`;
 `echo -n "   noncolin = " >> qefile`;
 `cat noncolin >> qefile`;
 `echo -n "   lspinorb = " >> qefile`;
 `cat spinorb >> qefile`;
 `echo -n "   ecutwfc = " >> qefile`;
 `cat ecut >> qefile`;
 #`echo "   ecutrho = " >> qefile`;
 #`cat ecutrho >> qefile`;
 `echo -n "   occupations = " >> qefile`;
 `cat occtype >> qefile`;
 `echo -n "   smearing = " >> qefile`;
 `cat smearing >> qefile`;
 `echo -n "   degauss = " >> qefile`;
 `cat degauss >> qefile`;
 `echo "   nosym = .true." >> qefile`;
 `echo "   noinv = .true." >> qefile`;
# `echo -n "   assume_isolated = " >> qefile`;
# `cat isolated >> qefile`;
 `echo "/" >> qefile`;

 `echo "&electrons" >> qefile`;
 `echo -n "   conv_thr    = " >> qefile`;
 `cat etol >> qefile`;
 `echo -n "   mixing_beta = " >> qefile`;
 `cat mixing >> qefile`;
 `echo -n "   electron_maxstep = " >> qefile`;
 `cat nrun >> qefile`;
 `echo "/" >> qefile`;

 `echo "&ions" >> qefile`;
 `echo "/" >> qefile`;

 `echo "ATOMIC_SPECIES" >> qefile`;
 `cat atompp >> qefile`;

 if ($ibrav == 0) {
    `echo "CELL_PARAMETERS cubic" >> qefile`;
    `cat acell >> qefile`;
 }

 `echo "ATOMIC_POSITIONS crystal" >> qefile`;
 `cat coords >> qefile`;

 #`echo "K_POINTS automatic" >> qefile`;
 #`tr '\n' ' ' < paw.nkpt >> qefile`;
 #`cat kshift >> qefile`;

 # paste k-point list onto end of scf.in file
  `echo "K_POINTS crystal " >> qefile`;
  `echo  $pawnkpts >> qefile`;
  `tail -$pawnkpts kpts4qe.0001 >> qefile`;

 # mv qefile to appropriate location
  `mv qefile nscf.in`;
   

### get the number of cores we are going to run on
##
  my $nc=0;
  open CORES, "core";
  $nc = int(<CORES>);
  `echo $nc >> testlog`;
  close CORES;


### the PAW NSCF run
#
  print "PAW NSCF Run\n";
  system("mpirun -np $nc $ENV{'OCEAN_BIN'}/pw.x < nscf.in > nscf.out 2>&1") == 0
      or die "Failed to run nscf for PAW\n";

  `echo 1 > scf.stat`;


 ## find the top of the valence bance
 system("$ENV{'OCEAN_BIN'}/qeband.x") == 0
     or die "Failed to count bands\n";

 ## obtain the number of occupied bands
 system("$ENV{'OCEAN_BIN'}/qebocc.x") == 0
     or die "Failed to count bands\n";


 ## get the number of occupied bands
 `grep 1.000 bands.out | wc -l > vb`;
 my $vb = 0;
 open BANDS, "vb" or die;
 $vb = <BANDS>;
 close BANDS;
                   

  my $natoms = `cat natoms`;
  my $fband = `cat fband`;
  $pawnbands = `cat paw.nbands`;
  my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
  $cb = 1 if ($cb < 1);
  open BRANGE, ">brange.ipt" or die;
  print BRANGE "1  $vb"
             . "$cb $pawnbands";
  close BRANGE;

  open STATUS, ">espresso.stat" or die;
  print STATUS "1";
  close STATUS;

####
  chdir "../";

} # end NSCF for PAW run




### Do NSCF for BSE run

if ( $bseRUN ) {
  chdir $bseDIR;
  print "BSE run\n";
 #  `cp ../bse.scf.in scf.in`;
 #  `cp ../bse.pp.in pp.in`;
 # copy all files over
  `cp -r ../Out/ .`;
  `cp ../acell .`;
  `cp ../atompp .`;
  `cp ../coords .`;
  foreach ( @GeneralFiles, @EspressoFiles, @PPFiles, @OtherFiles) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
  foreach ( "bse.nkpt", "nscf.kshift", "bse.nbands", "k0.ipt", "qinunitsofbvectors.ipt", "finalpplist", "cellinfo" ) {
    system("cp ../$_ .") == 0 or die "Failed to copy $_\n";
  }
 # run KGEN
  `cp bse.nkpt nkpt`;
  `cp bse.nbands nbands`;
  print "Running kgen2.x\n";
  `cp bse.nkpt kmesh.ipt`;
  system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";
  `echo "1" > kgen.stat`;


 my $ibrav = 0;
 open(IBRAV, 'ibrav') or die "couldn't open ibrav";
 $ibrav = <IBRAV>;
 close(IBRAV);


 my $line = "";
 my $celldm1 = 0;
 my $celldm2 = 0;
 my $celldm3 = 0;
 open(RSCALE, 'rscale') or die "couldn't open rscale";
 foreach $line (<RSCALE>) {
   ($celldm1, $celldm2, $celldm3) = split(' ' ,$line);
 }
 close(RSCALE);


 `echo "&control" > qefile`;
 `echo "   calculation = 'nscf'" >> qefile`;
 `echo -n "   prefix = '" >> qefile`;
 `head -c -1 prefix >> qefile`;
 `echo "'" >> qefile`;
 `echo -n "   pseudo_dir = " >> qefile`;
 `cat ppdir >> qefile`;
 `echo -n "   outdir = " >> qefile`;
 `cat work_dir >> qefile`;
 `echo "   tstress = .true." >> qefile`;
 `echo "   tprnfor = .true." >> qefile`;
 `echo "   wf_collect = .true." >> qefile`;
 `echo "/" >> qefile`;

 `echo "&system" >> qefile`;
 `echo -n "   ibrav = " >> qefile`;
 `cat ibrav >> qefile`;
 if ($ibrav != 0) {
   `cat cellinfo >> qefile`;
#   `echo -n "   celldm(1) = " >> qefile`;
#   `echo ${celldm1} >> qefile`;
 }
# `echo -n "   celldm(2) = " >> qefile`;
# `echo ${celldm2} >> qefile`;
# `echo -n "   celldm(3) = " >> qefile`;
# `echo ${celldm3} >> qefile`;
 `echo -n "   nat = " >> qefile`;
 `cat natoms >> qefile`;
 `echo -n "   ntyp = " >> qefile`;
 `cat ntype >> qefile`;
 `echo -n "   nbnd = " >> qefile`;
 `cat nbands >> qefile`;
 `echo -n "   noncolin = " >> qefile`;
 `cat noncolin >> qefile`;
 `echo -n "   lspinorb = " >> qefile`;
 `cat spinorb >> qefile`;
 `echo -n "   ecutwfc = " >> qefile`;
 `cat ecut >> qefile`;
 #`echo "   ecutrho = " >> qefile`;
 #`cat ecutrho >> qefile`;
 `echo -n "   occupations = " >> qefile`;
 `cat occtype >> qefile`;
 `echo -n "   smearing = " >> qefile`;
 `cat smearing >> qefile`;
 `echo -n "   degauss = " >> qefile`;
 `cat degauss >> qefile`;
 `echo "   nosym = .true." >> qefile`;
 `echo "   noinv = .true." >> qefile`;
# `echo -n "   assume_isolated = " >> qefile`;
# `cat isolated >> qefile`;
 `echo "/" >> qefile`;
 
 `echo "&electrons" >> qefile`;
 `echo -n "   conv_thr    = " >> qefile`;
 `cat etol >> qefile`;
 `echo -n "   mixing_beta = " >> qefile`;
 `cat mixing >> qefile`;
 `echo -n "   electron_maxstep = " >> qefile`;
 `cat nrun >> qefile`;
 `echo "/" >> qefile`;
 
 `echo "&ions" >> qefile`;
 `echo "/" >> qefile`;

 `echo "ATOMIC_SPECIES" >> qefile`;
 `cat atompp >> qefile`;

 if ($ibrav == 0) {
    `echo "CELL_PARAMETERS cubic" >> qefile`;
    `cat acell >> qefile`;
 }

 `echo "ATOMIC_POSITIONS crystal" >> qefile`;
 `cat coords >> qefile`;

 #`echo "K_POINTS automatic" >> qefile`;
 # paste k-point list onto end of scf.in file
 `echo "K_POINTS crystal " >> qefile`;
 `echo  $bsenkpts >> qefile`;
 `tail -$bsenkpts kpts4qe.0001 >> qefile`;
 #`tr '\n' ' ' < nkpt >> qefile`;
 #`cat kshift >> qefile`;
 
 # mv qefile to appropriate location
 `mv qefile nscf.in`;
 


### get the number of cores we are going to run on
##
  my $nc=0;
  open CORES, "core";
  $nc = int(<CORES>);
  `echo $nc >> testlog`;
  close CORES;
##


### running NSCF for the BSE stage
  print "BSE NSCF Run\n";

  system("mpirun -np $nc $ENV{'OCEAN_BIN'}/pw.x < nscf.in > nscf.out 2>&1") == 0
      or die "Failed to run nscf for BSE\n";
  `echo 1 > nscf.stat`;

 ## obtain the number of occupied bands
 system("$ENV{'OCEAN_BIN'}/qebocc.x") == 0
     or die "Failed to count bands\n";

  `grep 1.000 bands.out | wc -l > vb`;
  my $vb = 0;
  open BANDS, "vb" or die;
  $vb = <BANDS>;
  close BANDS;

  my $natoms = `cat natoms`;
  my $fband = `cat fband`;
  my $bands = `cat nbands`;
  my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
  $cb = 1 if ($cb < 1);
  open BRANGE, ">brange.ipt" or die;
  print BRANGE "1  $vb"
             . "$cb $bands";
  close BRANGE;

  open STATUS, ">espresso.stat" or die;
  print STATUS "1";
  close STATUS;
####
  chdir "../";

} # end NSCF for BSE run




#  system("$ENV{'OCEAN_BIN'}/avec.x") == 0 or die "Failed to run avec.x\n";
#  
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


print "Espresso stage complete\n";

exit 0;

