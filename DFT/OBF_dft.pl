#!/usr/bin/perl

use strict;
use File::Copy;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/OBF_dft\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (! $ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }
if (! $ENV{"OCEAN_ESPRESSO_PW"} ) {$ENV{"OCEAN_ESPRESSO_PW"} = $ENV{"OCEAN_BIN"} . "/obf_pw.x"; }
if (! $ENV{"OCEAN_ESPRESSO_PP"} ) {$ENV{"OCEAN_ESPRESSO_PP"} = $ENV{"OCEAN_BIN"} . "/obf_pp.x"; }

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
    "occopt", "occtype", "prefix", "ppdir", "rprim", "rscale", 
    "spinorb", "taulist", "typat", "verbatim", "work_dir", "wftol", 
    "den.kshift", "obkpt.ipt", "trace_tol", "ham_kpoints", "obf.nbands","tot_charge", "nspin", "smag", "ldau", "zsymb");
my @PPFiles = ("pplist", "znucl");
my @OtherFiles = ("epsilon", "pool_control");


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

my $coord_type = `cat coord`;
chomp($coord_type);



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

  my @qe_data_files = ('prefix', 'ppdir', 'work_dir', 'ibrav', 'natoms', 'ntype', 'noncolin',
                       'spinorb', 'ecut', 'occtype', 'degauss', 'etol', 'mixing', 'nrun', 'trace_tol', 'tot_charge', 'nspin', 'smag', 'ldau' );
  my %qe_data_files = {};
  foreach my $file_name (@qe_data_files)
  {
    open IN, $file_name or die "$file_name:  $!";
    my $string = <IN>;
    chomp $string;
    # Trim ', " and also leading or trailing spaces
    $string =~ s/\'//g;
    $string =~ s/\"//g;
    $string =~ s/^\s+//g;
    $string =~ s/\s+$//g;
    close IN;
    $qe_data_files{ "$file_name" } = $string;
  }

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

  my $string;
  open QE, ">scf.in" or die "Failed to open scf.in.\n$!";
  print QE "&control\n" 
        .  "  calculation = 'scf'\n"
        .  "  prefix = \'$qe_data_files{'prefix'}\'\n"
        .  "  pseudo_dir = \'$qe_data_files{'ppdir'}\'\n"
        .  "  outdir = \'$qe_data_files{'work_dir'}\'\n"
        .  "  tstress = .true.\n"
        .  "  tprnfor = .true.\n"
        .  "  wf_collect = .true.\n"
#        .  "  disk_io = 'low'\n"
        .  "/\n";
  print QE "&system\n"
        .  "  ibrav = $qe_data_files{'ibrav'}\n"
        .  "  nat = $qe_data_files{'natoms'}\n"
        .  "  ntyp = $qe_data_files{'ntype'}\n"
        .  "  noncolin = $qe_data_files{'noncolin'}\n"
        .  "  lspinorb = $qe_data_files{'spinorb'}\n"
        .  "  ecutwfc = $qe_data_files{'ecut'}\n"
        .  "  occupations = '$qe_data_files{'occtype'}'\n"
        .  "  degauss = $qe_data_files{'degauss'}\n"
        .  "  nspin  = $qe_data_files{'nspin'}\n"
        .  "  tot_charge  = $qe_data_files{'tot_charge'}\n"
        .  "  nosym = .true.\n"
        .  "  noinv = .true.\n";
  if( $qe_data_files{'smag'}  ne "" )
  {
    print QE "$qe_data_files{'smag'}\n";
  }
  if( $qe_data_files{'ldau'}  ne "" )
  {
    print QE "$qe_data_files{'ldau'}\n";
  }
  if( $qe_data_files{'ibrav'} != 0 ) 
  {
    print QE "  celldim(1) = ${celldm1}\n";
  }
  print QE "/\n"
        .  "&electrons\n"
        .  "  conv_thr = $qe_data_files{'etol'}\n"
        .  "  mixing_beta = $qe_data_files{'mixing'}\n"
        .  "  electron_maxstep = $qe_data_files{'nrun'}\n"
        .  "/\n"
        .  "&ions\n"
        .  "/\n";

  open IN, "atompp" or die "$!";
  my $atompp;
  while (<IN>) { $atompp .= $_; }
  close IN;
  chomp $atompp;
  print QE "ATOMIC_SPECIES\n" . $atompp . "\n";

  if ($ibrav == 0) {
    open IN, "acell" or die "$!";
    my $acell;
    while (<IN>) { $acell .= $_; }
    close IN;
    chomp $acell;
    print QE "CELL_PARAMETERS cubic\n" . $acell . "\n";
  }

  if( $coord_type =~ m/angst/ )
  {
    print QE "ATOMIC_POSITIONS angstrom\n";
  }
  else
  {
    print QE "ATOMIC_POSITIONS crystal\n";
  }
  open IN, "coords" or die;
  while (<IN>) { print QE $_;}
  close IN;

  print QE "K_POINTS automatic\n$ngkpt $kshift\n";

  close QE;

#  copy( "qefile", "scf.in");



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
  open PP, ">pp.in";
  print PP "&inputpp\n"
          . "  prefix = \'$qe_data_files{'prefix'}\'\n" 
          . "  outdir = \'$qe_data_files{'work_dir'}\'\n"
          . "  filplot= 'system.rho'\n"
          . "  plot_num = 0\n"
          . "/\n"
          . "&plot\n"
          . "  plot_center_atom = -1\n"
          . "  nfile = 1\n"
          . "  filepp(1) = 'system.rho', weight(1) = 1\n"
          . "  iflag = 3\n"
          . "  output_format = 6\n"
          . "  fileout = 'system.rho.dat'\n"
          . "  /\n";
  close PP;

# `echo "&inputpp" > ppfile`;
# `echo -n "   prefix = '" >> ppfile`;
# `head -c -1 prefix >> ppfile`;
# `echo "'" >> ppfile`;
# `echo -n "   prefix = " >> ppfile`;
# `cat prefix >> ppfile`;
# `echo -n "   outdir = " >> ppfile`;
# `cat work_dir >> ppfile`;
# `echo "   filplot = 'system.rho'" >> ppfile`;
# `echo "   plot_num = 0" >> ppfile`;
# `echo "/" >> ppfile`;
 
# `echo "&plot" >> ppfile`;
# `echo "   nfile = 1" >> ppfile`;
# `echo "   filepp(1) = 'system.rho', weight(1) = 1" >> ppfile`;
# `echo "   iflag = 3" >> ppfile`;
# `echo "   output_format = 6" >> ppfile`;
# `echo "   fileout = 'system.rho.dat'" >> ppfile`;
# `echo "/" >> ppfile`;
 
# `mv ppfile pp.in`;
 

          
# `echo NodeFile`;
# `echo $PBS_NODEFILE`;

  my $npool = 1;
  open INPUT, "pool_control" or die;
  while (<INPUT>)
  {
    if( $_ =~ m/^scf\s+(\d+)/ )
    {
      $npool = $1;
      last;
    }
  }
  close INPUT;


 ### the SCF run for initial density
 ##
 print "Density SCF Run\n";
  print "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < scf.in > scf.out 2>&1\n";
  system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < scf.in > scf.out 2>&1") == 0
      or die "Failed to run scf stage for Density\n";
 `echo 1 > scf.stat`;
 print "SCF complete\n";
# sleep( 30 );

 print "Density PP Run\n";
  print "$para_prefix $ENV{'OCEAN_ESPRESSO_PP'}  -npool $npool < pp.in > pp.out 2>&1\n";
  system("$para_prefix $ENV{'OCEAN_ESPRESSO_PP'} -npool $npool < pp.in > pp.out 2>&1") == 0
#  system("~/espresso-5.0.1/bin/pp.x < pp.in >& pp.out") == 0 
     or die "Failed to run density stage for PAW\n";
 `echo 1 > den.stat`;
# sleep( 30 );

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
    if( $line  =~  m/the Fermi energy is\s+([+-]?\d+\.?\d+)/ )
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

  my $nbands = `cat obf.nbands`;
  chomp($nbands);
  if( $nbands < 0 ) 
  {
    $nbands = `cat nbands`;
    chomp( $nbands );
  }
 ### PAW NSCF write input card
  
  open QE, ">nscf.in" or die "Failed to open nscf.in\n$!";
  print QE "&control\n"
        .  "  calculation = 'nscf'\n"
        .  "  prefix = \'$qe_data_files{'prefix'}\'\n"
        .  "  pseudo_dir = \'$qe_data_files{'ppdir'}\'\n"
        .  "  outdir = \'$qe_data_files{'work_dir'}\'\n"
        .  "  tstress = .true.\n"
        .  "  tprnfor = .true.\n"
        .  "  wf_collect = .true.\n"
#        .  "  disk_io = 'low'\n"
        .  "/\n";
  print QE "&system\n"
        .  "  ibrav = $qe_data_files{'ibrav'}\n"
        .  "  nat = $qe_data_files{'natoms'}\n"
        .  "  ntyp = $qe_data_files{'ntype'}\n"
        .  "  noncolin = $qe_data_files{'noncolin'}\n"
        .  "  lspinorb = $qe_data_files{'spinorb'}\n"
        .  "  ecutwfc = $qe_data_files{'ecut'}\n"
        .  "  occupations = '$qe_data_files{'occtype'}'\n"
        .  "  degauss = $qe_data_files{'degauss'}\n"
        .  "  nspin  = $qe_data_files{'nspin'}\n"
        .  "  tot_charge  = $qe_data_files{'tot_charge'}\n"
        .  "  nosym = .true.\n"
        .  "  noinv = .true.\n"
        .  "  nbnd = $nbands\n";
  if( $qe_data_files{'smag'}  ne "" )
  {
    print QE "$qe_data_files{'smag'}\n";
  }
  if( $qe_data_files{'ldau'}  ne "" )
  {
    print QE "$qe_data_files{'ldau'}\n";
  }
  if( $qe_data_files{'ibrav'} != 0 )
  {
    print QE "  celldim(1) = ${celldm1}\n";
  }
  print QE "/\n"
        .  "&electrons\n"
        .  "  conv_thr = $qe_data_files{'etol'}\n"
        .  "  mixing_beta = $qe_data_files{'mixing'}\n"
        .  "  electron_maxstep = $qe_data_files{'nrun'}\n"
        .  "/\n"
        .  "&ions\n"
        .  "/\n";

  open IN, "atompp" or die "$!";
  my $atompp;
  while (<IN>) { $atompp .= $_; }
  close IN;
  chomp $atompp;
  print QE "ATOMIC_SPECIES\n" . $atompp . "\n";

  if ($ibrav == 0) {
    open IN, "acell" or die "$!";
    my $acell;
    while (<IN>) { $acell .= $_; }
    close IN;
    chomp $acell;
    print QE "CELL_PARAMETERS cubic\n" . $acell . "\n";
  }

  if( $coord_type =~ m/angst/ )
  {
    print QE "ATOMIC_POSITIONS angstrom\n";
  }
  else
  {
    print QE "ATOMIC_POSITIONS crystal\n";
  }
  open IN, "coords" or die;
  while (<IN>) { print QE $_;}
  close IN;

  print QE  "K_POINTS automatic\n"
          . "$obfkpt 0 0 0\n";

#  `echo "&control" > qefile`;
#  `echo "   calculation = 'nscf'" >> qefile`;
#  `echo -n "   prefix = '" >> qefile`;
#  `head -c -1 prefix >> qefile`;
#  `echo "'" >> qefile`;
#  `echo -n "   pseudo_dir = " >> qefile`;
#  `cat ppdir >> qefile`;
#  `echo -n "   outdir = " >> qefile`;
#  `cat work_dir >> qefile`;
#  `echo "   tstress = .true." >> qefile`;
#  `echo "   tprnfor = .true." >> qefile`;
#  `echo "   wf_collect = .true." >> qefile`;
#  `echo "/" >> qefile`;
#
#  `echo "&system" >> qefile`;
#  `echo -n "   ibrav = " >> qefile`;
#  `cat ibrav >> qefile`;
#  if ($ibrav != 0) {
#    `echo -n "   celldm(1) = " >> qefile`;
#    `echo ${celldm1} >> qefile`;
#  }
# # `echo -n "   celldm(2) = " >> qefile`;
# # `echo ${celldm2} >> qefile`;
# # `echo -n "   celldm(3) = " >> qefile`;
# # `echo ${celldm3} >> qefile`;
#  `echo -n "   nat = " >> qefile`;
#  `cat natoms >> qefile`;
#  `echo -n "   ntyp = " >> qefile`;
#  `cat ntype >> qefile`;
#  `echo -n "   nbnd = " >> qefile`;
#  my $nbands = `cat obf.nbands`;
#  chomp($nbands);
#  if( $nbands < 0 ) 
#  {
#  	`cat nbands >> qefile`;
#  }
#  else
#  {
#	`cat obf.nbands >> qefile`;
#  }
#  `echo -n "   noncolin = " >> qefile`;
#  `cat noncolin >> qefile`;
#  `echo -n "   lspinorb = " >> qefile`;
#  `cat spinorb >> qefile`;
#  `echo -n "   ecutwfc = " >> qefile`;
#  `cat ecut >> qefile`;
#  #`echo "   ecutrho = " >> qefile`;
#  #`cat ecutrho >> qefile`;
#  `echo -n "   occupations = " >> qefile`;
#  `cat occtype >> qefile`;
#  `echo -n "   smearing = " >> qefile`;
#  `cat smearing >> qefile`;
#  `echo -n "   degauss = " >> qefile`;
#  `cat degauss >> qefile`;
#  `echo "   nosym = .true." >> qefile`;
#  `echo "   noinv = .true." >> qefile`;
# # `echo -n "   assume_isolated = " >> qefile`;
# # `cat isolated >> qefile`;
#  `echo "/" >> qefile`;

#  `echo "&electrons" >> qefile`;
#  `echo -n "   conv_thr    = " >> qefile`;
#  `cat etol >> qefile`;
#  `echo -n "   mixing_beta = " >> qefile`;
#  `cat mixing >> qefile`;
#  `echo -n "   electron_maxstep = " >> qefile`;
#  `cat nrun >> qefile`;
#  `echo "/" >> qefile`;
#
#  `echo "&ions" >> qefile`;
#  `echo "/" >> qefile`;

#  `echo "ATOMIC_SPECIES" >> qefile`;
#  `cat atompp >> qefile`;
#
#  if ($ibrav == 0) {
#     `echo "CELL_PARAMETERS cubic" >> qefile`;
#     `cat acell >> qefile`;
#  }

# `echo "ATOMIC_POSITIONS angstrom" >> qefile`;
#  `echo "ATOMIC_POSITIONS crystal" >> qefile`;
#if( $coord_type =~ m/angst/ )
#{
#        `echo "ATOMIC_POSITIONS angstrom" >> qefile`;
#}
#else
#{
#        `echo "ATOMIC_POSITIONS crystal" >> qefile`;
#}
#
#  `cat coords >> qefile`;

  #`echo "K_POINTS automatic" >> qefile`;
  #`tr '\n' ' ' < paw.nkpt >> qefile`;
  #`cat kshift >> qefile`;

 # paste k-point list onto end of scf.in file
#  `echo "K_POINTS automatic" >> qefile`;
#  `echo "$obfkpt 0 0 0 " >> qefile`;

 # mv qefile to appropriate location
#  `mv qefile nscf.in`;
   
  my $npool = 1;
  open INPUT, "pool_control" or die;
  while (<INPUT>)
  {
    if( $_ =~ m/obf\s+(\d+)/ )
    {
      $npool = $1;
      last;
    }
  }
  close INPUT;


### the PAW NSCF run
#
  print "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} -npool $npool < nscf.in > nscf.out 2>&1\n";
  system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} -npool $npool < nscf.in > nscf.out 2>&1") == 0
     or die "Failed to run nscf stage for PAW\n";
  print "NSCF complete\n";
# sleep( 30 );

  `echo 1 > nscf.stat`;
} # end NSCF for PAW run


#  
#print `pwd`;
#open RSCALE, "rscale" or die;
#open RPRIM, "rprim" or die;
#<RSCALE> =~  m/(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)/ or die;
#my @rscale = ($1, $3, $5);
#print "$1\t$3\t$5\n";
#close RSCALE;

#open AVECS, ">avecsinbohr.ipt" or die;
#for (my $i = 0; $i < 3; $i++ ) {
#  <RPRIM> =~  m/([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)/ or die "$_";
#  print AVECS $1*$rscale[0] . "  " . $3*$rscale[1] .  "  " . $5*$rscale[2] . "\n";
#  print "$1\t$3\t$5\n";
#}
#close RPRIM;
#close AVECS;


print "Create Basis\n";
open BASIS, ">basis.in" or die "Failed top open basis.in\n$!";
print BASIS "&input\n"
        .  "  prefix = \'$qe_data_files{'prefix'}\'\n"
        .  "  outdir = \'$qe_data_files{'work_dir'}\'\n"
        .  "  trace_tol = $qe_data_files{'trace_tol'}\n"
        .  "/\n";
close BASIS;

#`echo "&input" > basis.in`;
#`echo -n "   prefix = '" >> basis.in`;
#`head -c -1 prefix >> basis.in`;
#`echo "'" >> basis.in`;
#`echo -n "   outdir = " >> basis.in`;
#`cat work_dir >> basis.in`;
#`echo -n "   trace_tol = " >> basis.in`;
#`cat trace_tol >> basis.in`;
#`echo "/" >> basis.in`;

system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_basis.x  < basis.in > basis.out 2>&1") == 0
      or die "Failed to run shirley_basis.x\n$!";

my $ham_kpoints = `cat ham_kpoints`;
chomp $ham_kpoints;

print "Create Shirley Hamiltonian\n";
open HAM, ">ham.in" or die "Failed to open ham.in\n$!";
print HAM "&input\n"
        . "  prefix = 'system_opt'\n"
        . "  outdir = \'$qe_data_files{'work_dir'}\'\n"
        . "  updatepp = .false.\n"
        . "  ncpp = .true.\n"
        . "  nspin_ham = $qe_data_files{'nspin'}\n"
        . "/\n"
        . "K_POINTS\n"
        . "$ham_kpoints 0 0 0\n";
close HAM;

#`echo "&input" > ham.in`;
#`echo "   prefix = 'system_opt'" >> ham.in`;
#`echo -n "   outdir = " >> ham.in`;
#`cat work_dir >> ham.in`;
#`echo "   updatepp = .false." >> ham.in`;
#`echo "   ncpp = .true." >> ham.in`;
#`echo "/" >> ham.in`;
#`echo " K_POINTS" >> ham.in`;
#`echo "$ham_kpoints 0 0 0" >> ham.in`;


system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_ham_o.x  < ham.in > ham.out 2>&1") == 0
      or die "Failed to run shirley_ham_o.x\n$!";

print "Espresso stage complete\n";

exit 0;


