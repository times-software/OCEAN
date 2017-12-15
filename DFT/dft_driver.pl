#!/usr/bin/perl
# Copyright (C) 2015 - 2017 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
# At present this runs JDFTx, but the plan is to move all supported DFT codes to 
# be controlled by this script


use strict;
use File::Copy;
use Cwd 'abs_path';

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/dft_driver\.pl/;
  $ENV{"OCEAN_BIN"} = abs_path($1);
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}


if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (! $ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }
if (! $ENV{"OCEAN_ESPRESSO_PW"} ) {$ENV{"OCEAN_ESPRESSO_PW"} = $ENV{"OCEAN_BIN"} . "/pw.x"; }
if (! $ENV{"OCEAN_ESPRESSO_PP"} ) {$ENV{"OCEAN_ESPRESSO_PP"} = $ENV{"OCEAN_BIN"} . "/pp.x"; }
if (! $ENV{"OCEAN_ESPRESSO_OBF_PW"} ) 
    {$ENV{"OCEAN_ESPRESSO_OBF_PW"} = $ENV{"OCEAN_BIN"} . "/obf_pw.x"; }
if (! $ENV{"OCEAN_ESPRESSO_OBF_PP"} ) 
    {$ENV{"OCEAN_ESPRESSO_OBF_PP"} = $ENV{"OCEAN_BIN"} . "/obf_pp.x"; }
if (! $ENV{"OCEAN_JDFTX"} ) {$ENV{"OCEAN_JDFTX"} = $ENV{"OCEAN_BIN"} . "/jdftx"; }

####################################

my $RunKGen = 0;
my $RunPP = 0;
my $runDFT = 0;
my $nscfRUN = 1;

my @GeneralFiles = ("para_prefix", "calc");

my @KgenFiles = ("nkpt", "k0.ipt", "qinunitsofbvectors.ipt", "screen.nkpt");
my @BandFiles = ("nbands", "screen.nbands");
my @EspressoFiles = ( "coord", "degauss", "ecut", "etol", "fband", "ibrav", 
    "isolated", "mixing", "natoms", "ngkpt", "noncolin", "nrun", "ntype", 
    "occopt", "prefix", "ppdir", "rprim", "rscale", "metal",
    "spinorb", "taulist", "typat", "verbatim", "work_dir", "tmp_dir", "wftol", 
    "den.kshift", "obkpt.ipt", "trace_tol", "ham_kpoints", "obf.nbands","tot_charge", 
    "nspin", "smag", "ldau", "qe_scissor", "zsymb", "dft.calc_stress", "dft.calc_force", "dft.split", "dft",
    "dft.startingwfc", "dft.diagonalization" );
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
  $runDFT = 1;
}
else {
  foreach (@EspressoFiles) {
    if ( -e $_ ) {
      if ( `diff -q $_ ../Common/$_` ) {
        $runDFT = 1;
        last;
      }
    }
    else {
      $runDFT = 1;
      last;
    }
  }
}
unless ($runDFT) {
 $runDFT = 1;
  if (open STATUS, "espresso.stat") {
    if (<STATUS> == 1) { $runDFT = 0; }
  }
  close STATUS;
}
if ($runDFT) {
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


open GOUT, ">dft.stat" or die;
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

open DFT, "dft" or die "Failed to open dft\n";
<DFT> =~ m/(\w+)/ or die "Failed to parse dft\n";
my $input_dft_type = $1;
close DTF;
my $dft_type = 0;
my $obf;
if( $input_dft_type =~ m/obf/i )
{
  $dft_type = 'obf';
  $obf = 1;
  print "Running DFT calculation with OBF extension\n"
}
elsif( $input_dft_type =~ m/qe/i )
{
  $dft_type = 'qe';
  print "Running DFT calculation using QuantumESPRESSO\n";  
}
elsif( $input_dft_type =~ m/abi/i )
{
  $dft_type = 'abi';
  print "Running DFT calculation using ABINIT\n";
  exit 1;
}
elsif( $input_dft_type =~ m/jdf/i )
{
  $dft_type = 'jdf';
  print "Running DFT calculation using JDFTx\n";
}
else
{
  print "Unrecognized DFT option. Will not continue\n";
  exit 2;
}

open IN, "calc" or die "Failed to open calc\n";
<IN> =~m/(\w+)/ or die "Failed to parse calc\n";
my $calc = $1;
close IN;
my $run_screen = 1;
if( $calc =~ m/val/i )
{
  $run_screen = 0;
}


if ($RunPP) {
  system("$ENV{'OCEAN_BIN'}/pp.pl znucl pplist finalpplist") == 0
    or die "Failed to run pp.pl\n";
  `echo "1" > pp.stat`;
}


#############################################

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



my @dft_data_files = ('prefix', 'ppdir', 'work_dir', 'tmp_dir', 'ibrav', 'natoms', 'ntype', 'noncolin',
                     'spinorb', 'ecut', 'degauss', 'etol', 'mixing', 'nrun', 'occopt',
                     'trace_tol', 'tot_charge', 'nspin', 'ngkpt', 'k0.ipt', 'metal',
                     'den.kshift', 'obkpt.ipt', 'obf.nbands', 'nkpt', 'nbands', 'screen.nbands',
                     'screen.nkpt', 'dft.calc_stress', 'dft.calc_force', 'dft.startingwfc', 
                     'dft.diagonalization' );



my %dft_data_files = {};
foreach my $file_name (@dft_data_files)
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
    $dft_data_files{ "$file_name" } = $string;
}
my $line = "";
my $celldm1 = 0;
my $celldm2 = 0;
my $celldm3 = 0;
open(RSCALE, 'rscale') or die "couldn't open rscale\n$!";
foreach $line (<RSCALE>) {
 ($celldm1, $celldm2, $celldm3) = split(' ' ,$line);
}
close(RSCALE);
$dft_data_files{ "celldm1" } = $celldm1;
$dft_data_files{ "celldm2" } = $celldm2;
$dft_data_files{ "celldm3" } = $celldm3;


# Switch ppdir to absolute path
$dft_data_files{ "ppdir" } = abs_path( $dft_data_files{ "ppdir" } ) . "/";


#QE optional files
my @qe_opt_files = ('acell', 'coords', 'atompp', 'smag', 'ldau', 'qe_scissor' );
foreach my $file_name (@qe_opt_files)
{
    open IN, $file_name or die "$file_name:  $!";
    my $string;
    while( my $a = <IN> ) 
    {
      $string .= $a;
    }
    chomp $string;
    
    close IN;
    $dft_data_files{ "$file_name" } = $string;
}
##################

# Map QE/Abinit occupation options
if( $dft_data_files{ "occopt" } < 1 || $dft_data_files{ "occopt" } > 7 )
{
  print "WARNING! Occopt set to a non-sensical value. Changing to 3";
  $dft_data_files{ "occopt" } = 3;
}
# Don't support abinit 2
$dft_data_files{ "occopt" } = 1 if( $dft_data_files{ "occopt" } == 2 );
# Override occopt if metal was specified
if( $dft_data_files{ 'metal' } =~ m/true/i )
{
  if( $dft_data_files{ "occopt" } == 1 )
  {
    print "WARNING! Mismatch between occopt and metal flags.\n  Setting occopt to 3\n";
    $dft_data_files{ "occopt" } = 3;
  }
}

$dft_data_files{'occtype'} = 'smearing';
# At the moment we are leaving QE as smearing even if occopt = 1
#   therefore we want to clamp down the smearing a bunch
if( $dft_data_files{ "occopt" } == 1 )
{
  $dft_data_files{'occtype'} = 'fixed';
  $dft_data_files{ 'degauss' } = 0.002;
}

# Array of QE names for smearing by occopt
my @QE_smear;
$QE_smear[1] = "'gaussian'";     # Still need to fix to be insulator
$QE_smear[3] = "'fermi-dirac'";  # ABINIT = fermi-dirac
$QE_smear[4] = "'marzari-vanderbilt'";  # ABINIT = Marzari cold smearing a = -0.5634
$QE_smear[5] = "'marzari-vanderbilt'";  # ABINIT = Marzari a = -0.8165
$QE_smear[6] = "'methfessel-paxton'";  # ABINIT = Methfessel and Paxton PRB 40, 3616 (1989)
$QE_smear[7] = "'gaussian'";     # ABINIT = Gaussian



if ($runDFT) {


 ### write SCF input card for initial density

  open my $QE, ">scf.in" or die "Failed to open scf.in.\n$!";

  # Set the flags that change for each input/dft run
  $dft_data_files{'calctype'} = 'scf';
  $dft_data_files{'print kpts'} = "K_POINTS automatic\n$dft_data_files{'ngkpt'} $dft_data_files{'den.kshift'}\n";
  $dft_data_files{'print nbands'} = -1;


  &print_qe( $QE, %dft_data_files );

  close $QE;


 ## SCF PP initialize and set defaults
 
 ### write PP input card for density
  open PP, ">pp.in";
  print PP "&inputpp\n"
          . "  prefix = \'$dft_data_files{'prefix'}\'\n" 
          . "  outdir = \'$dft_data_files{'work_dir'}\'\n"
          . "  filplot= 'system.rho'\n"
          . "  plot_num = 0\n"
          . "/\n";
  close PP;


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
  if( $obf == 1 )
  {
    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'}  -npool $npool < scf.in > scf.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'}  -npool $npool < scf.in > scf.out 2>&1") == 0
        or die "Failed to run scf stage for Density\n";
  }
  else
  {
    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < scf.in > scf.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < scf.in > scf.out 2>&1") == 0
        or die "Failed to run scf stage for Density\n";
  }
  open OUT, ">scf.stat" or die "Failed to open scf.stat\n$!";
  print OUT "1\n";
  close OUT;
  print "SCF complete\n";

  print "Density PP Run\n";
  if( $obf == 1 )
  {  
    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'}  -npool $npool < pp.in > pp.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'} -npool $npool < pp.in > pp.out 2>&1") == 0
       or die "Failed to run density stage for SCREENING\n";
  }
  else
  {
    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PP'}  -npool $npool < pp.in > pp.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_PP'} -npool $npool < pp.in > pp.out 2>&1") == 0
       or die "Failed to run density stage for SCREENING\n";
  }
  open OUT, ">den.stat" or die "Failed to open scf.stat\n$!";
  print OUT "1\n";
  close OUT;

  ## convert the density file to proper format
  print "Density conversion\n";
  system("$ENV{'OCEAN_BIN'}/qe2rhoofr.pl" ) == 0 
    or die "Failed to convert density\n$!\n";



  open STATUS, ">espresso.stat" or die;
  print STATUS "1";
  close STATUS;


  # Find Fermi level and number of electrons
  my $fermi = 'no';
  my $nelectron = 'no';
  my $units;

  # First attempt to grab from outfile
  my $data_file = $dft_data_files{'work_dir'} . "/" . $dft_data_files{'prefix'} . ".save/data-file.xml";
  print "Looking for $data_file \n";
  if( open SCF, $data_file )
  {
    while( my $scf_line = <SCF> )
    {
      if( $scf_line =~ m/\<UNITS_FOR_ENERGIES UNITS=\"(\w+)/ )
      {
        $units = $1;
      }
      if( $scf_line =~ m/\<FERMI_ENERGY/ )
      {
        $scf_line = <SCF>;
        $scf_line =~ m/([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ or die "$scf_line";
        $fermi = $1;
      }
      if( $scf_line =~m/\<NUMBER_OF_ELECTRONS/ )
      {
        $scf_line = <SCF>;
        $scf_line =~ m/(\d+\.\d+[Ee]?[-+]?(\d+)?)/ or die "$scf_line";
        $nelectron = $1;
      }
    }
    close SCF;
  }
  
  if( $fermi =~ m/no/ || $nelectron =~ m/no/ )
  {
    open SCF, "scf.out" or die "$!";
    while( my $line = <SCF> )
    {
      if( $line  =~  m/the Fermi energy is\s+([+-]?\d+\.?\d+)/ )
      {
        $fermi = $1;
        print "Fermi level found at $fermi eV\n";
        $fermi = $fermi/13.60569253;
      }
      if( $line =~ m/number of electrons\s+=\s+(\d+)/ )
      {
        $nelectron = $1;
      }
    }
    close SCF;
  }
  else
  {
    if( $units =~ m/hartree/i )
    {
      $fermi *= 2;
    }
    elsif( $units =~ m/eV/i )
    {
      $fermi /= 13.60569253;
    }
    my $eVfermi = $fermi * 13.60569253;
    print "Fermi level found at $eVfermi eV\n";
  }

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

#JTV
  my $line = "";

  #IF( OBF ) then "Single run, main directory"

  #ELSE( is QE ) then "2 runs for screening and BSE wavefunctions

  my $nbands = $dft_data_files{'obf.nbands'};
  $nbands = $dft_data_files{'nbands'} if ( $nbands < 1 );
  
  if( $obf == 1 )
  {
    open QE, ">nscf.in" or die "Failed to open nscf.in\n$!";
    print QE "&control\n"
          .  "  calculation = 'nscf'\n"
          .  "  prefix = \'$dft_data_files{'prefix'}\'\n"
          .  "  pseudo_dir = \'$dft_data_files{'ppdir'}\'\n"
          .  "  outdir = \'$dft_data_files{'work_dir'}\'\n"
          .  "  wfcdir = \'$dft_data_files{'tmp_dir'}\'\n"
#          .  "  tstress = $dft_data_files{'dft.calc_stress'}\n"
#          .  "  tprnfor = $dft_data_files{'dft.calc_force'}\n"
          .  "  wf_collect = .true.\n"
  #        .  "  disk_io = 'low'\n"
          .  "/\n";
    print QE "&system\n"
          .  "  ibrav = $dft_data_files{'ibrav'}\n"
          .  "  nat = $dft_data_files{'natoms'}\n"
          .  "  ntyp = $dft_data_files{'ntype'}\n"
          .  "  noncolin = $dft_data_files{'noncolin'}\n"
          .  "  lspinorb = $dft_data_files{'spinorb'}\n"
          .  "  ecutwfc = $dft_data_files{'ecut'}\n"
          .  "  occupations = '$dft_data_files{'occtype'}'\n"
          .  "  degauss = $dft_data_files{'degauss'}\n"
          .  "  nspin  = $dft_data_files{'nspin'}\n"
          .  "  tot_charge  = $dft_data_files{'tot_charge'}\n"
          .  "  nosym = .true.\n"
          .  "  noinv = .true.\n"
          .  "  nbnd = $nbands\n";
    if( $dft_data_files{'smag'}  ne "" )
    {
      print QE "$dft_data_files{'smag'}\n";
    }
    if( $dft_data_files{'ldau'}  ne "" )
    {
      print QE "$dft_data_files{'ldau'}\n";
    }
    if( $dft_data_files{'qe_scissor'}  ne "" )
    {
      print QE "$dft_data_files{'qe_scissor'}\n";
    }
    if( $dft_data_files{'ibrav'} != 0 )
    {
      print QE "  celldim(1) = ${celldm1}\n";
    }
    print QE "/\n"
          .  "&electrons\n"
          .  "  conv_thr = $dft_data_files{'etol'}\n"
          .  "  mixing_beta = $dft_data_files{'mixing'}\n"
          .  "  electron_maxstep = $dft_data_files{'nrun'}\n"
          .  "  startingwfc = \'$dft_data_files{'dft.startingwfc'}\'\n"
          .  "  diagonalization = \'$dft_data_files{'dft.diagonalization'}\'\n"
          .  "/\n"
          .  "&ions\n"
          .  "/\n";

    open IN, "atompp" or die "$!";
    my $atompp;
    while (<IN>) { $atompp .= $_; }
    close IN;
    chomp $atompp;
    print QE "ATOMIC_SPECIES\n" . $atompp . "\n";

    if ($dft_data_files{'ibrav'} == 0) {
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
    elsif( $coord_type =~ m/bohr/ )
    {
      print QE "ATOMIC_POSITIONS bohr\n";
    }
    else
    {
      print QE "ATOMIC_POSITIONS crystal\n";
    }
    open IN, "coords" or die;
    while (<IN>) { print QE $_;}
    close IN;

    print QE  "K_POINTS automatic\n"
            . "$dft_data_files{'obkpt.ipt'} 0 0 0\n";
    close QE;

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


    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} -npool $npool < nscf.in > nscf.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} -npool $npool < nscf.in > nscf.out 2>&1") == 0
       or die "Failed to run nscf stage for OBFs\n";
    print "NSCF complete\n";

    open OUT, ">nscf.stat" or die "Failed to open nscf.stat\n$!";
    print OUT "1\n";
    close OUT;

    print "Create Basis\n";
    open BASIS, ">basis.in" or die "Failed top open basis.in\n$!";
    print BASIS "&input\n"
             .  "  prefix = \'$dft_data_files{'prefix'}\'\n"
             .  "  outdir = \'$dft_data_files{'work_dir'}\'\n";
    unless( $dft_data_files{'tmp_dir'} =~ m/undefined/ )
    {
#      print "$dft_data_files{'tmp_dir'}\n";
      print BASIS "  wfcdir = \'$dft_data_files{'tmp_dir'}\'\n";
    }
    print BASIS "  trace_tol = $dft_data_files{'trace_tol'}\n"
             .  "/\n";
    close BASIS;

    system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_basis.x  < basis.in > basis.out 2>&1") == 0
          or die "Failed to run shirley_basis.x\n$!";

    my $ham_kpoints = `cat ham_kpoints`;
    chomp $ham_kpoints;

    print "Create Shirley Hamiltonian\n";
    open HAM, ">ham.in" or die "Failed to open ham.in\n$!";
    print HAM "&input\n"
            . "  prefix = 'system_opt'\n"
            . "  outdir = \'$dft_data_files{'work_dir'}\'\n";
    unless( $dft_data_files{'tmp_dir'} =~ m/undefined/ )
    {
#      print "$dft_data_files{'tmp_dir'}\n";
      print HAM "  wfcdir = \'$dft_data_files{'tmp_dir'}\'\n";
    }
    print HAM "  updatepp = .false.\n"
            . "  ncpp = .true.\n"
            . "  nspin_ham = $dft_data_files{'nspin'}\n"
            . "/\n"
            . "K_POINTS\n"
            . "$ham_kpoints 0 0 0\n";
    close HAM;

    system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_ham_o.x  < ham.in > ham.out 2>&1") == 0
          or die "Failed to run shirley_ham_o.x\n$!";

  }
  else
  {

    my $bseDIR = sprintf("%03u%03u%03u", split( /\s+/,$dft_data_files{'nkpt'})); 
    mkdir $bseDIR unless ( -d $bseDIR );
    chdir $bseDIR;

    # kpts
    copy "../nkpt", "nkpt";
    copy "../qinunitsofbvectors.ipt", "qinunitsofbvectors.ipt";
    copy "../k0.ipt", "k0.ipt";
    copy "../dft.split", "dft.split";

    my $split_dft;
    if( open IN, "dft.split" )
    {
      $split_dft = 0;
      if( <IN> =~ m/t/i )
      {
        $split_dft = 1;
      }
      close IN;
      # Only makes sense if we have q
      if( $split_dft )
      {
        open IN, "qinunitsofbvectors.ipt" or die "Failed to open qinunitofbvectors\n$!";
        <IN> =~ m/([+-]?\d+\.?\d*([eE][+-]?\d+)?)\s+([+-]?\d+\.?\d*([eE][+-]?\d+)?)\s+([+-]?\d+\.?\d*([eE][+-]?\d+)?)/ 
                    or die "Failed to parse qinunitsofbvectors.ipt\n";
        my $fake_qmag = abs($1) + abs($3) + abs($5);
        close IN;
        $split_dft = 0 if( $fake_qmag < 0.000000001 );
      }
    }
    else
    {
      $split_dft = 0;
    } 

    print "DFT runs will be split\n" if( $split_dft );

    $dft_data_files{'prefix_shift'} = $dft_data_files{'prefix'} . "_shift";

    mkdir "Out" unless ( -d "Out" );

    # This will loop back and do everything for prefix_shift if we have split
    my $repeat = 0;
    $repeat = 1 if( $split_dft );
    my $prefix = 'prefix';
    for( my $i = 0; $i <= $repeat; $i++ )
    {
      mkdir "Out/$dft_data_files{$prefix}.save" unless ( -d "Out/$dft_data_files{$prefix}.save" );

      copy "../Out/$dft_data_files{'prefix'}.save/charge-density.dat", "Out/$dft_data_files{$prefix}.save/charge-density.dat";
      copy "../Out/$dft_data_files{'prefix'}.save/data-file.xml", "Out/$dft_data_files{$prefix}.save/data-file.xml";


      if( $dft_data_files{'nspin'} == 2 )
      {
        copy "../Out/$dft_data_files{'prefix'}.save/spin-polarization.dat", 
             "Out/$dft_data_files{$prefix}.save/spin-polarization.dat";
      }

      if( $dft_data_files{'ldau'}  ne "" )
      {
        # Starting w/ QE-6.0 this is the DFT+U info from the SCF
        if( -e "../Out/$dft_data_files{'prefix'}.save/occup.txt" )
        {
          copy "../Out/$dft_data_files{'prefix'}.save/occup.txt", "Out/$dft_data_files{$prefix}.save/occup.txt";
        }
        # QE 4.3-5.x
        elsif( -e "../Out/$dft_data_files{'prefix'}.occup" )
        {
          copy "../Out/$dft_data_files{'prefix'}.occup", "Out/$dft_data_files{$prefix}.occup";
        }
      }
      $prefix = "prefix_shift";
    }


#    copy "../acell", "acell";
#    copy "../atompp", "atompp";
#    copy "../coords", "coords";
#    open OUT, ">core" or die;
#    print OUT "1\n";
#    close OUT;
    system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";

    open my $QE, ">nscf.in" or die "Failed to open nscf.in\n$!";

    # Set the flags that change for each input/dft run
    $dft_data_files{'calctype'} = 'nscf';
    my $kpt_text = "K_POINTS crystal\n";
    open IN, "nkpts" or die;
    my $nkpts = <IN>;
    close IN;
    $kpt_text .= $nkpts;
    open IN, "kpts4qe.0001" or die;
    while(<IN>)
    {
      $kpt_text .= $_;
    }
    close IN;

    $dft_data_files{'print kpts'} = $kpt_text;
    unless( $split_dft ) 
    {
      $dft_data_files{'print nbands'} = $dft_data_files{'nbands'};
    }

    &print_qe( $QE, %dft_data_files );

    close $QE;

    my $npool = 1;
    open INPUT, "../pool_control" or die;
    while (<INPUT>)
    {
      if( $_ =~ m/^nscf\s+(\d+)/ )
      {
        $npool = $1;
        last;
      }
    }
    close INPUT;

    print "BSE NSCF Run\n";
    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < nscf.in > nscf.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < nscf.in > nscf.out 2>&1") == 0
        or die "Failed to run nscf stage for BSE wavefunctions\n";


    if( $split_dft )
    {
      open my $QE, ">nscf_shift.in" or die "Failed to open nscf_shift.in\n$!";

      $prefix = $dft_data_files{'prefix'};
      $dft_data_files{'prefix'} = $dft_data_files{'prefix_shift'};

      $kpt_text = "K_POINTS crystal\n";
      open IN, "nkpts" or die;
      my $nkpts = <IN>;
      close IN;
      $kpt_text .= $nkpts;
      open IN, "kpts4qe.0002" or die;
      while(<IN>)
      {
        $kpt_text .= $_;
      }
      close IN;
      $dft_data_files{'print kpts'} = $kpt_text;
      $dft_data_files{'print nbands'} = $dft_data_files{'nbands'};

      &print_qe( $QE, %dft_data_files );

      close $QE;

      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < nscf_shift.in > nscf_shift.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < nscf_shift.in > nscf_shift.out 2>&1") == 0
          or die "Failed to run nscf stage for shifted BSE wavefunctions\n";

      $dft_data_files{'prefix'} = $prefix;

    }



    open OUT, ">nscf.stat" or die "Failed to open nscf.stat\n$!";
    print OUT "1\n";
    close OUT;
    print "BSE NSCF complete\n";

    ## find the top of the valence bance
    system("$ENV{'OCEAN_BIN'}/qeband.pl") == 0
       or die "Failed to count bands\n$!\n";

    open IN, "brange.stub" or die;
    open OUT, ">brange.ipt" or die;
    while(<IN>)
    {
      print OUT $_;
    }
    print OUT "$dft_data_files{'nbands'}\n";
    close IN;
    close OUT;

    copy "nkpt", "kmesh.ipt";

    chdir "../";
  }
}

if( $obf == 0 && $run_screen == 1 )
{
  
#  my $bseDIR = sprintf("%03u%03u%03u", split( /\s+/,$dft_data_files{'screen.nkpt'}));
  my $bseDIR = "SCREEN";
  mkdir $bseDIR unless ( -d $bseDIR );
  chdir $bseDIR;

  mkdir "Out" unless ( -d "Out" );
  mkdir "Out/$dft_data_files{'prefix'}.save" unless ( -d "Out/$dft_data_files{'prefix'}.save" );

  copy "../Out/$dft_data_files{'prefix'}.save/charge-density.dat", "Out/$dft_data_files{'prefix'}.save/charge-density.dat";
  copy "../Out/$dft_data_files{'prefix'}.save/data-file.xml", "Out/$dft_data_files{'prefix'}.save/data-file.xml";

  if( $dft_data_files{'nspin'} == 2 )
  {
    copy "../Out/$dft_data_files{'prefix'}.save/spin-polarization.dat", 
         "Out/$dft_data_files{'prefix'}.save/spin-polarization.dat";
  }

  if( $dft_data_files{'ldau'}  ne "" )
  {
    # Starting w/ QE-6.0 this is the DFT+U info from the SCF
    if( -e "../Out/$dft_data_files{'prefix'}.save/occup.txt" )
    {
      copy "../Out/$dft_data_files{'prefix'}.save/occup.txt", "Out/$dft_data_files{'prefix'}.save/occup.txt";
    }
    # QE 4.3-5.x
    elsif( -e "../Out/$dft_data_files{'prefix'}.occup" )
    {
      copy "../Out/$dft_data_files{'prefix'}.occup", "Out/$dft_data_files{'prefix'}.occup";
    }
  }


  # kpts
  copy "../screen.nkpt", "nkpt";
  copy "../k0.ipt", "k0.ipt";

#  copy "../acell", "acell";
#  copy "../atompp", "atompp";
#  copy "../coords", "coords";

  # QINB is 0 for screening
  open OUT, ">qinunitsofbvectors.ipt" or die;
  print OUT "0.0 0.0 0.0\n";
  close OUT;
  open OUT, ">core" or die;
  print OUT "1\n";
  close OUT;
  system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";

  open my $QE, ">nscf.in" or die "Failed to open nscf.in\n$!";

  $dft_data_files{'calctype'} = 'nscf';
  my $kpt_text = "K_POINTS crystal\n";
  open IN, "nkpts" or die;
  my $nkpts = <IN>;
  close IN;
  $kpt_text .= $nkpts;
  open IN, "kpts4qe.0001" or die;
  while(<IN>)
  {
    $kpt_text .= $_;
  }
  close IN;
  
  $dft_data_files{'print kpts'} = $kpt_text;
  $dft_data_files{'print nbands'} = $dft_data_files{'screen.nbands'};

  &print_qe( $QE, %dft_data_files );

  close QE;

  my $npool = 1;
  open INPUT, "../pool_control" or die;
  while (<INPUT>)
  {
    if( $_ =~ m/^screen\s+(\d+)/ )
    {
      $npool = $1;
      last;
    }
  }
  close INPUT;

  print "Screening NSCF Run\n";
  print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < nscf.in > nscf.out 2>&1\n";
  system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'}  -npool $npool < nscf.in > nscf.out 2>&1") == 0
      or die "Failed to run nscf stage for SCREENing wavefunctions\n";
  open OUT, ">nscf.stat" or die "Failed to open nscf.stat\n$!";
  print OUT "1\n";
  close OUT;
  print "Screening NSCF complete\n";

    ## find the top of the valence bance
    system("$ENV{'OCEAN_BIN'}/qeband.pl") == 0
       or die "Failed to count bands\n$!\n";

  open IN, "brange.stub" or die;
  open OUT, ">brange.ipt" or die;
  while(<IN>)
  {
    print OUT $_;
  }
  print OUT "$dft_data_files{'screen.nbands'}\n";
  close IN;
  close OUT;

  copy "nkpt", "kmesh.ipt";

  chdir "../";
}



print "Espresso stage complete\n";
exit 0;


# Subroutine to print a QE-style input file
#   Pass in a file handle and a hashtable
sub print_qe 
{
  my ($fh, %inputs ) = @_;

  print $fh "&control\n"
        .  "  calculation = \'$inputs{'calctype'}\'\n"
        .  "  prefix = \'$inputs{'prefix'}\'\n"
        .  "  pseudo_dir = \'$inputs{'ppdir'}\'\n"
        .  "  outdir = \'$inputs{'work_dir'}\'\n"
        .  "  wfcdir = \'$inputs{'tmp_dir'}\'\n"
        .  "  tstress = $inputs{'dft.calc_stress'}\n"
        .  "  tprnfor = $inputs{'dft.calc_force'}\n"
        .  "  wf_collect = .true.\n"
        .  "/\n";
  print $fh "&system\n"
        .  "  ibrav = $inputs{'ibrav'}\n"
        .  "  nat = $inputs{'natoms'}\n"
        .  "  ntyp = $inputs{'ntype'}\n"
        .  "  noncolin = $inputs{'noncolin'}\n"
        .  "  lspinorb = $inputs{'spinorb'}\n"
        .  "  ecutwfc = $inputs{'ecut'}\n"
        .  "  occupations = '$inputs{'occtype'}'\n"
        .  "  smearing = $QE_smear[$inputs{'occopt'}]\n"
        .  "  degauss = $inputs{'degauss'}\n"
        .  "  nspin  = $inputs{'nspin'}\n"
        .  "  tot_charge  = $inputs{'tot_charge'}\n"
        .  "  nosym = .true.\n"
        .  "  noinv = .true.\n";
  if( $inputs{'print nbands'} > 0 ) # for scf no nbnd is set. 
                                    # Therefore -1 is passed in and nothing is written to the input file
  {
    print $fh "  nbnd = $inputs{'print nbands'}\n";
  }
  if( $inputs{'smag'}  ne "" )
  {
    print $fh "$inputs{'smag'}\n";
  }
  if( $inputs{'ldau'}  ne "" )
  {
    print $fh "$inputs{'ldau'}\n";
  }
  if( $inputs{'qe_scissor'}  ne "" )
  {
    print $fh "$inputs{'qe_scissor'}\n";
  }
  if( $inputs{'ibrav'} != 0 )
  {
    print $fh "  celldim(1) = $inputs{'celldm1'}\n";
  }
  print $fh "/\n"
        .  "&electrons\n"
        .  "  conv_thr = $inputs{'etol'}\n"
        .  "  mixing_beta = $inputs{'mixing'}\n"
        .  "  electron_maxstep = $inputs{'nrun'}\n"
        .  "  startingwfc = \'$dft_data_files{'dft.startingwfc'}\'\n"
        .  "  diagonalization = \'$dft_data_files{'dft.diagonalization'}\'\n"
        .  "/\n"
        .  "&ions\n"
        .  "/\n";

#  open IN, "atompp" or die "$!";
#  my $atompp;
#  while (<IN>) { $atompp .= $_; }
#  close IN;
#  chomp $atompp;
#  print $fh "ATOMIC_SPECIES\n" . $atompp . "\n";
  print $fh "ATOMIC_SPECIES\n" . $inputs{'atompp'} . "\n";

  if ($inputs{'ibrav'} == 0) {
#    open IN, "acell" or die "$!";
#    my $acell;
#    while (<IN>) { $acell .= $_; }
#    close IN;
#    chomp $acell;
#    print $fh "CELL_PARAMETERS cubic\n" . $acell . "\n";
    print $fh "CELL_PARAMETERS cubic\n" . $inputs{'acell'} . "\n";
  }

  if( $coord_type =~ m/angst/ )
  {
    print $fh "ATOMIC_POSITIONS angstrom\n";
  }
  elsif( $coord_type =~ m/bohr/ )
  {
    print $fh "ATOMIC_POSITIONS bohr\n";
  }
  else
  {
    print $fh "ATOMIC_POSITIONS crystal\n";
  }
#  open IN, "coords" or die;
#  while (<IN>) { print $fh $_;}
#  close IN;
  print $fh $inputs{'coords'} . "\n";

  print $fh $inputs{'print kpts'};

}

  
