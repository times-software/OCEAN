#!/usr/bin/perl

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/OBF_dft\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (! $ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }
if (! $ENV{"OCEAN_ESPRESSO_PW"} ) {$ENV{"OCEAN_ESPRESSO_PW"} = $ENV{"OCEAN_BIN"} . "/pw.x"; }
if (! $ENV{"OCEAN_ESPRESSO_PP"} ) {$ENV{"OCEAN_ESPRESSO_PP"} = $ENV{"OCEAN_BIN"} . "/pp.x"; }

####################################
# Executables to be run
# kgen.x
# PP
# ABINIT
# avecs?

my $RunKGen = 0;
my $RunPP = 0;
my $RunESPRESSO = 0;
my $nscfRUN = 1;

my @GeneralFiles = ("para_prefix" );#, "getden");

my @KgenFiles = ("nkpt", "k0.ipt", "qinunitsofbvectors.ipt", "paw.nkpt");
my @BandFiles = ("nbands", "paw.nbands");
my @EspressoFiles = ( "coord", "degauss", "ecut", "etol", "fband", "ibrav", 
    "isolated", "mixing", "natoms", "ngkpt", "noncolin", "nrun", "ntype", 
    "occopt", "occtype", "prefix", "ppdir", "rprim", "rscale", "smearing", 
    "spinorb", "taulist", "typat", "verbatim", "work_dir", "wftol", 
    "den.kshift", "obkpt.ipt", "trace_tol");
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
#    `rm -r $file`;
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
foreach (@BandFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

foreach (@PPFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

foreach (@EspressoFiles, @OtherFiles) {
  system("cp ../Common/$_ .") == 0 or die;
} 

#############################################


if ($RunPP) {
  system("$ENV{'OCEAN_BIN'}/pp.pl znucl pplist finalpplist") == 0
    or die "Failed to run pp.pl\n";
  `echo "1" > pp.stat`;
}

#############################################
# Run espresso
# Determin what type of run is being done, par or seq

my $EspressoType = "seq"; 

### load up the para_prefix
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
#`echo 1 > core`;



# make additional files for QE input card

# Coords are wrong, currently
print "making the coordinates";
system("$ENV{'OCEAN_BIN'}/makecoords.x") == 0
    or die "Failed to make coordinates\n";

print "making acell";
system("$ENV{'OCEAN_BIN'}/makeacell.x") == 0
    or die "Failed to make acell\n";

print "making atompp";
system("$ENV{'OCEAN_BIN'}/makeatompp.x") == 0
    or die "Failed to make acell\n";


if ($RunESPRESSO) {

 open(NGKP, 'ngkpt') or die "couldn't open ngkpt\n$!";
 my $ngkpt = <NGKP>;
 chomp($ngkpt);
 close(NGKPT);

 open(KSHIFT, 'den.kshift') or die "couldn't open den.kshift\n$!";
 my $kshift = <KSHIFT>;
 chomp($kshift);
 close(KSHIFT);


 my $ibrav = 0;
 open(IBRAV, 'ibrav') or die "couldn't open ibrav\n$!";
 $ibrav = <IBRAV>;
 close(IBRAV);

 my $line = "";
 my $celldm1 = 0;
 my $celldm2 = 0;
 my $celldm3 = 0;
 open(RSCALE, 'rscale') or die "couldn't open rscale\n$!";
 foreach $line (<RSCALE>) {
   ($celldm1, $celldm2, $celldm3) = split(' ' ,$line);
 }
 close(RSCALE);

 ### write SCF input card for initial density

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
   `echo -n "   celldm(1) = " >> qefile`;
   `echo ${celldm1} >> qefile`;
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

#test 2013-4-5
# `echo "ATOMIC_POSITIONS angstrom" >> qefile`;
`echo "ATOMIC_POSITIONS crystal" >> qefile`;
 `cat coords >> qefile`;


 `echo "K_POINTS automatic" >> qefile`;
 #`tr '\n' ' ' < paw.nkpt >> qefile`;
 #`cat kshift >> qefile`;
# `echo "2 2 1 1 1 0" >> qefile`;
  `echo "$ngkpt $kshift" >> qefile`;
 #KG#

 # mv qefile to appropriate location
 `mv qefile scf.in`;



 ## SCF PP initialize and set defaults
 
 #KG# my $prefix = "'atom'\n";
 my $outdir = "'./Out'\n";
 
 my $filplot = "'system.rho'\n";
 my $plotnum = "0\n";
 
 my $nfile = "1\n";
 my $filepp = "'system.rho', ";
 my $weight = "1.0\n";
 my $iflag = "3\n";
 my $outformat = "6\n";
 my $fileout = "'system.rho.dat'\n";
 
 
 ### write PP input card for density
 
 `echo "&inputpp" > ppfile`;
 `echo -n "   prefix = '" >> ppfile`;
 `head -c -1 prefix >> ppfile`;
 `echo "'" >> ppfile`;
# `echo -n "   prefix = " >> ppfile`;
# `cat prefix >> ppfile`;
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
 

          
# `echo NodeFile`;
# `echo $PBS_NODEFILE`;

 ### the SCF run for initial density
 ##
 print "Density SCF Run\n";
# system("mpiexec /global/home/users/kgilmore/Code/QEspresso-4.2.1/bin/pw.x < scf.in >& scf.out") == 0
#     or die "Failed to run scf stage for Density\n";
#  system("mpirun -np $nc /global/home/users/kgilmore/Code/QEspresso-4.2.1/bin/pw.x < scf.in > scf.out 2>&1") == 0
#      or die "Failed to run scf stage for Density\n";
#  system("mpirun -np $nc $ENV{'OCEAN_ESPRESSO_PW'}  < scf.in > scf.out 2>&1") == 0
  system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  < scf.in > scf.out 2>&1") == 0
      or die "Failed to run scf stage for Density\n";
 `echo 1 > scf.stat`;

 print "Density PP Run\n";
  system("$para_prefix $ENV{'OCEAN_ESPRESSO_PP'} < pp.in > pp.out 2>&1") == 0
     or die "Failed to run density stage for PAW\n";
 `echo 1 > den.stat`;

 ## convert the density file to proper format
 print "Density conversion\n";
 system("$ENV{'OCEAN_BIN'}/converter.py system.rho.dat") == 0
     or die "Failed to convert density\n$!\n";

# ## find the top of the valence bance
# system("$ENV{'OCEAN_BIN'}/qeband.x") == 0
#     or die "Failed to count bands\n$!\n";
#
# ## obtain the number of occupied bands
# system("$ENV{'OCEAN_BIN'}/qebocc.x") == 0
#     or die "Failed to count bands\n";
#      
#
# ## get the number of occupied bands
# `grep 1.000 bands.out | wc -l > vb`;
# my $vb = 0;
# open BANDS, "vb" or die;
# $vb = <BANDS>;
# close BANDS;
#
# ## make brange file
# my $natoms = `cat natoms`;
# my $fband = `cat fband`;
# my $nbands = `cat nbands`;
# my $cb = sprintf("%.0f", $vb - 2*$natoms*$fband);
# $cb = 1 if ($cb < 1);
# open BRANGE, ">brange.ipt" or die;
# print BRANGE "1  $vb"
#            . "$cb $nbands";
# close BRANGE;

 open STATUS, ">espresso.stat" or die;
 print STATUS "1";
 close STATUS;

 `cp scf.out density.out`;


  my $fermi = 'no';
  my $nelectron = 'no';

  open SCF, "scf.out" or die "$!";
  while( my $line = <SCF> )
  {
    if( $line  =~  m/the Fermi energy is\s+([+-]?\d?\.?\d+)/ )
    {
      $fermi = $1;
      print "Fermi level found at $fermi eV\n";
      $fermi = $fermi/13.60569252;
    }
    if( $line =~ m/number of electrons\s+=\s+(\d+)/ )
    {
      $nelectron = $1;
    }
  }
  close SCF;

  die "Fermi level not found in scf.out\n" if( $fermi eq 'no' ) ;
  die "Number of electrons not found in scf.out\n" if( $nelectron eq 'no' );

  open FERMI, ">efermiinrydberg.ipt" or die "Failed to open efermiinrydberg\n$!";
  print FERMI "$fermi\n";
  close FERMI;

  open NELECTRON, ">nelectron" or die "Failed to open nelectron\n$!";
  print NELECTRON "$nelectron\n";
  close NELECTRON;


} # end SCF for density
      



### Do NSCF run

if ( $nscfRUN ) {
  print "NSCF run\n";

  my $ibrav = 0;
  open(IBRAV, 'ibrav') or die "couldn't open ibrav\n$!";
  $ibrav = <IBRAV>;
  close(IBRAV);


  my $line = "";
  my $celldm1 = 0;
  my $celldm2 = 0;
  my $celldm3 = 0;
  open(RSCALE, 'rscale') or die "couldn't open rscale\n$!";
  foreach $line (<RSCALE>) {
    ($celldm1, $celldm2, $celldm3) = split(' ' ,$line);
  }
  close(RSCALE);

  open(OBFKPT, 'obkpt.ipt') or die "couldn't open obkpt.ipt\n$!";
  my $obfkpt = <OBFKPT>;
  chomp($obfkpt);
  close(OBFKPT);


 ### PAW NSCF write input card
  
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
    `echo -n "   celldm(1) = " >> qefile`;
    `echo ${celldm1} >> qefile`;
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

  `echo "ATOMIC_POSITIONS angstrom" >> qefile`;
  `cat coords >> qefile`;

  #`echo "K_POINTS automatic" >> qefile`;
  #`tr '\n' ' ' < paw.nkpt >> qefile`;
  #`cat kshift >> qefile`;

 # paste k-point list onto end of scf.in file
  `echo "K_POINTS automatic" >> qefile`;
  `echo "$obfkpt 0 0 0 " >> qefile`;

 # mv qefile to appropriate location
  `mv qefile nscf.in`;
   


### the PAW NSCF run
#
  system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} < nscf.in > nscf.out 2>&1") == 0
     or die "Failed to run nscf stage for PAW\n";

  `echo 1 > nscf.stat`;
} # end NSCF for PAW run


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


print "Create Basis\n";
`echo "&input" > basis.in`;
`echo -n "   prefix = '" >> basis.in`;
`head -c -1 prefix >> basis.in`;
`echo "'" >> basis.in`;
`echo -n "   outdir = " >> basis.in`;
`cat work_dir >> basis.in`;
`echo -n "   trace_tol = " >> basis.in`;
`cat trace_tol >> basis.in`;
`echo "/" >> basis.in`;

system("time $para_prefix $ENV{'OCEAN_BIN'}/shirley_basis.x  < basis.in > basis.out 2>&1") == 0
      or die "Failed to run shirley_basis.x\n$!";

print "Create Shirley Hamiltonian\n";
`echo "&input" > ham.in`;
`echo "   prefix = 'system_opt'" >> ham.in`;
`echo -n "   outdir = " >> ham.in`;
`cat work_dir >> ham.in`;
`echo "   updatepp = .false." >> ham.in`;
`echo "   ncpp = .true." >> ham.in`;
`echo "/" >> ham.in`;
`echo " K_POINTS" >> ham.in`;
`echo "4 4 4 0 0 0" >> ham.in`;


system("time $para_prefix $ENV{'OCEAN_BIN'}/shirley_ham_o.x  < ham.in > ham.out 2>&1") == 0
      or die "Failed to run shirley_ham_o.x\n$!";

print "Espresso stage complete\n";

exit 0;


