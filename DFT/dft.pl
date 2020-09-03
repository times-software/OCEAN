#!/usr/bin/perl
# Copyright (C) 2015 - 2019 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#


use strict;
use File::Copy;
use Cwd 'abs_path';
use File::Compare;
use File::Spec::Functions;
use POSIX;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/dft\.pl/;
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

####################################

my $RunKGen = 0;
my $RunPP = 0;
my $RunESPRESSO = 0;
my $nscfRUN = 0;
my $run_screen = 0;

my @GeneralFiles = ("para_prefix", "calc");

my @KgenFiles = ("nkpt", "k0.ipt", "qinunitsofbvectors.ipt", "screen.nkpt", "screen.k0", "dft.split");
my @BandFiles = ("nbands", "screen.nbands");
my @EspressoFiles = ( "coord", "degauss", "ecut", "etol", "fband", "ibrav", 
    "isolated", "mixing", "natoms", "ngkpt", "noncolin", "nrun", "ntype", 
    "occopt", "prefix", "ppdir", "rprim", "rscale", "metal",
    "spinorb", "taulist", "typat", "verbatim", "work_dir", "tmp_dir", "wftol", 
    "den.kshift", "obkpt.ipt", "trace_tol", "ham_kpoints", "obf.nbands","tot_charge", 
    "nspin", "smag", "ldau", "qe_scissor", "zsymb", "dft.calc_stress", "dft.calc_force", "dft",
    "dft.startingwfc", "dft.diagonalization", "dft.qe_redirect", "dft.ndiag", "dft.functional", "dft.exx.qmesh", 
    "ngkpt.auto" );
my @PPFiles = ("pplist", "znucl");
my @OtherFiles = ("epsilon", "pool_control", "screen.mode");

my @SCFBonus = ("charge-density.kin.dat");
my @exx = ("hse");

unless( -e "scf.stat" )
{
  $RunPP = 1;
}


foreach (@PPFiles) {
  if ( -e $_ ) {
    if( compare( "$_", "../Common/$_") != 0 )
    {
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

if ( $RunPP ) {
  $RunESPRESSO = 1;
}
else {
  foreach (@EspressoFiles) {
    if ( -e $_ ) {
      if( compare( "$_", "../Common/$_") != 0 )
      {
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

if ($RunESPRESSO) {
  print "Differences found for density run. Clearing all old data\n";
  my @dirlisting = <*>;
  foreach my $file (@dirlisting) {
    chomp($file);
#    `rm -r $file`;
  }
  $RunPP = 1;
  $nscfRUN = 1;
  $run_screen = 1;
  unlink "scf.stat";
}
else {
  `touch old`;
}

unless( $nscfRUN == 1)
{
  foreach( "nkpt", "k0.ipt", "nbands" )
  {
    if( compare( "$_", "../Common/$_") != 0 )
    {
      $nscfRUN = 1;
      print "Difference found in $_\n";
      last;
    }
  }
  unless( $nscfRUN == 1 )
  {
    if( compare( "qinunitsofbvectors.ipt", "../Common/qinunitsofbvectors.ipt" ) != 0 )
    {
      $nscfRUN = 2;
      print "Difference found in qinunitsofbvectors.ipt\n";
    }
  }
}
unless( $run_screen == 1)
{
  foreach( "screen.nkpt", "screen.k0", "screen.nbands" )
  {
    if( compare( "$_", "../Common/$_") != 0 )
    {
      $run_screen = 1;
      print "Difference found in $_\n";
      last;
    }
  }
}

if( $nscfRUN == 0 )
{
  $nscfRUN = 1 unless( -e "bse.stat" );
}
if( $run_screen == 0 )
{
  $run_screen = 1 unless( -e "screen.stat" );
}

foreach (@GeneralFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@KgenFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@BandFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

if( $RunPP == 1 )
{
  foreach (@PPFiles) {
    system("cp ../Common/$_ .") == 0 or die;
  }
}

open IN, "calc" or die "Failed to open calc\n";
<IN> =~m/(\w+)/ or die "Failed to parse calc\n";
my $calc = $1;
close IN;

my $old_screen_mode;
if( -e "screen.mode" )
{
  open IN, "screen.mode" or die "Failed to open screen.mode\n$!";
  <IN> =~m/(\w+)/ or die "Failed to parse screen.mode\n";
  $old_screen_mode = $1;
  close IN;
}
else
{
  $old_screen_mode = '';
}


foreach (@EspressoFiles, @OtherFiles) {
  system("cp ../Common/$_ .") == 0 or die;
} 

open IN, "screen.mode" or die "Failed to open screen.mode";
<IN> =~m/(\w+)/ or die "Failed to parse screen.mode\n";
my $screen_mode = $1;
close IN;

open IN, "screen.mode" or die "Failed to open screen.mode\n";
<IN> =~m/(\w+)/ or die "Failed to parse screen.mode\n";
my $screen_mode = $1;
close IN;
if( $calc =~ m/val/i )
{
  $run_screen = 0 unless( $screen_mode =~ m/grid/i );
}
if( $run_screen == 0 && $screen_mode =~ m/grid/i )
{
  unless( $old_screen_mode =~ m/grid/i )
  {
    print "Need screening for valence: $old_screen_mode\n";
    $run_screen = 1;
  }
}

if( $nscfRUN != 0 )
{
  unlink "bse.stat";
}
if( $run_screen == 1 )
{
   unlink "screen.stat";
}
#############################################

open DFT, "dft" or die "Failed to open dft\n";
<DFT> =~ m/(\w+)/ or die "Failed to parse dft\n";
my $dft_type = $1;
close DTF;
my $obf;
if( $dft_type =~ m/obf/i )
{
  $obf = 1;
  print "Running DFT calculation with OBF extension\n"
}
else
{
  $obf = 0;
  print "Running DFT calculation using QE\n";
}


# Input to QE can be done via redirect (legacy) or -inp (more stable)
open IN, "dft.qe_redirect" or die "Failed to open dft.qe_redirect\n$!";
my $qe_redirect = <IN>;
close IN;
chomp( $qe_redirect );
if( $qe_redirect =~ m/f/i )
{
  $qe_redirect = 0;
}
elsif( $qe_redirect =~ m/t/i )
{
  $qe_redirect = 1;
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

if( -e "../Common/atompp" ) {
  copy "../Common/atompp", "atompp";
}
else {
  print "making atompp";
  move( "pplist", "pplist.hold");
  open IN, "pplist.hold" or die;
  open OUT, ">", "pplist" or die;
  while( my $line = <IN> )
  {
    chomp $line;
    $line =~ s/.upf$//i;
    print OUT $line . "\n";
  }
  close OUT;
  close IN;
  system("$ENV{'OCEAN_BIN'}/makeatompp.x") == 0
     or die "Failed to make acell\n";
  move( "pplist.hold", "pplist");
}



my @qe_data_files = ('prefix', 'ppdir', 'work_dir', 'tmp_dir', 'ibrav', 'natoms', 'ntype', 'noncolin',
                     'spinorb', 'ecut', 'degauss', 'etol', 'mixing', 'nrun', 'occopt',
                     'trace_tol', 'tot_charge', 'nspin', 'ngkpt', 'k0.ipt', 'metal',
                     'den.kshift', 'obkpt.ipt', 'obf.nbands', 'nkpt', 'nbands', 'screen.nbands',
                     'screen.nkpt', 'dft.calc_stress', 'dft.calc_force', 'dft.startingwfc', 
                     'dft.diagonalization', 'dft.ndiag', 'dft.functional' );



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
my $line = "";
my $celldm1 = 0;
my $celldm2 = 0;
my $celldm3 = 0;
open(RSCALE, 'rscale') or die "couldn't open rscale\n$!";
foreach $line (<RSCALE>) {
 ($celldm1, $celldm2, $celldm3) = split(' ' ,$line);
}
close(RSCALE);
$qe_data_files{ "celldm1" } = $celldm1;
$qe_data_files{ "celldm2" } = $celldm2;
$qe_data_files{ "celldm3" } = $celldm3;

#Set startingpot
$qe_data_files{ "dft.startingpot" } = 'atomic';

# Switch ppdir to absolute path
$qe_data_files{ "ppdir" } = abs_path( $qe_data_files{ "ppdir" } ) . "/";


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
    $qe_data_files{ "$file_name" } = $string;
}

# Load up qmesh for EXX
open EXX, "dft.exx.qmesh" or die "Failed to open dft.exx.qmesh\n$!";
<EXX> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse dft.exx.qmesh\n";
$qe_data_files{'nqx1'} = $1;
$qe_data_files{'nqx2'} = $2;
$qe_data_files{'nqx3'} = $3;
close EXX;
##################

# Map QE/Abinit occupation options
if( $qe_data_files{ "occopt" } < 1 || $qe_data_files{ "occopt" } > 7 )
{
  print "WARNING! Occopt set to a non-sensical value. Changing to 3";
  $qe_data_files{ "occopt" } = 3;
}
# Don't support abinit 2
$qe_data_files{ "occopt" } = 1 if( $qe_data_files{ "occopt" } == 2 );
# Override occopt if metal was specified
if( $qe_data_files{ 'metal' } =~ m/true/i )
{
  if( $qe_data_files{ "occopt" } == 1 )
  {
    print "WARNING! Mismatch between occopt and metal flags.\n  Setting occopt to 3\n";
    $qe_data_files{ "occopt" } = 3;
  }
}

$qe_data_files{'occtype'} = 'smearing';
# At the moment we are leaving QE as smearing even if occopt = 1
#   therefore we want to clamp down the smearing a bunch
if( $qe_data_files{ "occopt" } == 1 )
{
  $qe_data_files{'occtype'} = 'fixed';
  $qe_data_files{ 'degauss' } = 0.002;
}

# Array of QE names for smearing by occopt
my @QE_smear;
$QE_smear[1] = "'gaussian'";     # Still need to fix to be insulator
$QE_smear[3] = "'fermi-dirac'";  # ABINIT = fermi-dirac
$QE_smear[4] = "'marzari-vanderbilt'";  # ABINIT = Marzari cold smearing a = -0.5634
$QE_smear[5] = "'marzari-vanderbilt'";  # ABINIT = Marzari a = -0.8165
$QE_smear[6] = "'methfessel-paxton'";  # ABINIT = Methfessel and Paxton PRB 40, 3616 (1989)
$QE_smear[7] = "'gaussian'";     # ABINIT = Gaussian



if ($RunESPRESSO) {


  unlink "scf.stat";


 ## SCF PP initialize and set defaults
 
 ### write PP input card for density
  open PP, ">pp.in";
  print PP "&inputpp\n"
          . "  prefix = \'$qe_data_files{'prefix'}\'\n" 
          . "  outdir = \'$qe_data_files{'work_dir'}\'\n"
          . "  filplot= 'system.rho'\n"
          . "  plot_num = 0\n"
          . "/\n";
  close PP;

 ### write PP input card for total potential
  open PP, ">pp2.in";
  print PP "&inputpp\n"
          . "  prefix = \'$qe_data_files{'prefix'}\'\n"
          . "  outdir = \'$qe_data_files{'work_dir'}\'\n"
          . "  filplot= 'system.pot'\n"
          . "  plot_num = 1\n"
          . "/\n";
  close PP;



  my $npool = 1;
  my $ncpus = 1;
  open INPUT, "pool_control" or die;
  while (<INPUT>)
  {
    if( $_ =~ m/^scf\s+(\d+)/ )
    {
      $npool = $1;
    }
    elsif( $_ =~ m/^total\s+(\d+)/ )
    { 
      $ncpus = $1;
    }

  }
  close INPUT;

  print "TEST\n";
  print $qe_data_files{'dft.ndiag'} . "\n";
  if( $qe_data_files{'dft.ndiag'} =~ m/(-?\d+)/ )
  {
    if( $1 > 0 )
    {
      $qe_data_files{'dft.ndiag'} = $1;
    }
    else
    {
      $qe_data_files{'dft.ndiag'} = $ncpus;
    }
  }
  else
  {
    $qe_data_files{'dft.ndiag'} = 4;
  }

  my $scfConv = 0.0000073502388828 * $qe_data_files{'natoms'};
  if( $qe_data_files{'etol'} =~ m/(\d+\.?\d?([edED][+-]?\d+)?)/ )
  {
    my $conv_thr = $1;
    $conv_thr =~ s/d/e/i;
    print "$conv_thr\n";
    $scfConv = 20.0*$conv_thr if ( $scfConv < 20.0*$conv_thr );
  }
  print "$scfConv\n";
  my $oldSCFEnergy = 0;
  my $SCFEnergy;

  my $scfcountmax = 0;
  if( open INPUT, "ngkpt.auto" )
  {
    $scfcountmax = 6 if( <INPUT> =~ m/T/i );
    close INPUT;
  }

  for( my $scfcount = 0; $scfcount < $scfcountmax; $scfcount++ )
  {
   ### write SCF input card for initial density

    open my $QE, ">scf.in" or die "Failed to open scf.in.\n$!";

    # Set the flags that change for each input/dft run
    $qe_data_files{'calctype'} = 'scf';
    $qe_data_files{'print kpts'} = "K_POINTS automatic\n$qe_data_files{'ngkpt'} $qe_data_files{'den.kshift'}\n";
    $qe_data_files{'print nbands'} = -1;
    if( $obf == 1 ) 
    {
      $qe_data_files{'nosym'} = '.true.';
      $qe_data_files{'noinv'} = '.true.';
    }
    else
    {
      $qe_data_files{'nosym'} = '.false.';
      $qe_data_files{'noinv'} = '.false.';
    }

    # Check for Gamma-only, and over-write 'print kpts'
    $qe_data_files{'ngkpt'} =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "$qe_data_files{'ngkpt'}";
    if( $1 * $2 * $3 == 1 )
    {
      $qe_data_files{'den.kshift'} =~ m/(\S+)\s+(\S+)\s+(\S+)/ or die "$qe_data_files{'den.kshift'}";
      unless ( abs($1) > 0.000001 || abs($2) > 0.000001 || abs($3) > 0.000001 )
      {
        $qe_data_files{'print kpts'} = "K_POINTS gamma\n";
      }
      else { print "KSHIFT: $1  $2  $3\n"; }
    }
    else
    { print "KPOINTS: $1  $2  $3\n"; }

    &print_qe( $QE, %qe_data_files );

    close $QE;



    my $scf_prefix = $para_prefix;
    if( $obf != 1 ) 
    {
      print "Testing parallel QE execution\n";
      my $ser_prefix = $para_prefix;
      $ser_prefix =~ s/\d+/1/;
      open TMP, '>', "$qe_data_files{'prefix'}.EXIT" or die "Failed to open file $qe_data_files{'prefix'}.EXIT\n$!";
      close TMP;
      if( $qe_redirect )
      {
        print  "$ser_prefix $ENV{'OCEAN_ESPRESSO_PW'} < scf.in > test.out 2>&1\n";
        system("$ser_prefix $ENV{'OCEAN_ESPRESSO_PW'} < scf.in >test.out 2>&1");
      }
      else
      {
        print  "$ser_prefix $ENV{'OCEAN_ESPRESSO_PW'} -inp scf.in > test.out 2>&1\n";
        system("$ser_prefix $ENV{'OCEAN_ESPRESSO_PW'} -inp scf.in > test.out 2>&1");
      }

      if( open TMP, "test.out" )
      {
        my $actualKpts = -1;
        my $numKS;
        while (<TMP>)
        {
          if( $_ =~ m/number of Kohn-Sham states=\s+(\d+)/ )
          {
            $numKS = $1;
          }
          if( $_ =~ m/number of k points=\s+(\d+)/ )
          {
            $actualKpts = $1;
            last;
          }
        }
        close TMP;
        if( $actualKpts == -1 )
        {
          print "Had trouble parsing test.out\nDidn't find number of k points\n";
        }
        else
        {
          if( $actualKpts > $ncpus )
          {
            $npool = $ncpus;
          }
          else
  #        if( $npool > $actualKpts )
          {
            for( my $i = 1; $i <= $actualKpts; $i++ )
            {
              $npool = $i if(  $ncpus % $i == 0 );
            }
          }
          print "SCF has $actualKpts k-points\nWill use $npool pools\n";
        }
        if( defined( $numKS ) )
        {
          my $maxProcs = $numKS * $npool;
          print "   $ncpus  $maxProcs\n";
          if( $maxProcs < $ncpus )
          {
            $scf_prefix =~ s/\d+/$maxProcs/;
          }
        }
      }
      else
      {
        print "Had trouble parsing test.out\n. Will attempt to continue.\n";
      }
    }

   ### the SCF run for initial density
   ##
    print "Density SCF Run\n";
    my $qeCommandLine = "-ndiag $qe_data_files{'dft.ndiag'} -npool $npool";
    if( $obf == 1 )
    {
      if( $qe_redirect ) 
      {
        print  "$scf_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} $qeCommandLine < scf.in > scf.out 2>&1\n";
        system("$scf_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} $qeCommandLine < scf.in > scf.out 2>&1") == 0
            or die "Failed to run scf stage for Density\n";
      }
      else
      {
        print  "$scf_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} $qeCommandLine -inp scf.in > scf.out 2>&1\n";
        system("$scf_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} $qeCommandLine -inp scf.in > scf.out 2>&1") == 0
            or die "Failed to run scf stage for Density\n";
      }
    }
    else
    {
      if( $qe_redirect )
      {    
        print  "$scf_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < scf.in > scf.out 2>&1\n";
        system("$scf_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < scf.in > scf.out 2>&1") == 0
            or die "Failed to run scf stage for Density\n";
      } 
      else
      {
        print  "$scf_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp scf.in > scf.out 2>&1\n";
        system("$scf_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp scf.in > scf.out 2>&1") == 0
            or die "Failed to run scf stage for Density\n";
      }
    }

#    my $SCFEnergy = `grep ! scf.out`;
    `grep ! scf.out` =~ m/(-?\d+\.\d+)\s+Ry/ or die "Failed to parse scf.out\n";
    $SCFEnergy = $1;
    copy( "scf.out", "scf.out.$scfcount" );
    copy( "scf.in", "scf.in.$scfcount" );
    if( $scfcount > 1 && abs( $SCFEnergy - $oldSCFEnergy ) < $scfConv )
    {
      print abs( $SCFEnergy - $oldSCFEnergy ) . "   $scfConv\n";
      last;
    }
    my @ngkpt = split ' ', $qe_data_files{'ngkpt'};
    my @acell = split ' ', $qe_data_files{'acell'};
    my $kgden;
    my $testden;
    print "$qe_data_files{'ngkpt'}\n$ngkpt[0]  $ngkpt[1]  $ngkpt[2]\n";
    # Length of bvector is 1/avector (ignoring 2pi)
    # denisty is Ng / length(b) = Ng * length(a)
    $testden = $ngkpt[0] * sqrt( $acell[0]**2 + $acell[1]**2 + $acell[2]**2 );
    $kgden = $testden;
    $testden = $ngkpt[1] * sqrt( $acell[3]**2 + $acell[4]**2 + $acell[5]**2 );
    $kgden = $testden if( $testden > $kgden );
    $testden = $ngkpt[2] * sqrt( $acell[6]**2 + $acell[7]**2 + $acell[8]**2 );
    $kgden = $testden if( $testden > $kgden );
    $kgden += 0.1;
    $ngkpt[0] = ceil( $kgden/sqrt( $acell[0]**2 + $acell[1]**2 + $acell[2]**2 ) );
    $ngkpt[1] = ceil( $kgden/sqrt( $acell[3]**2 + $acell[4]**2 + $acell[5]**2 ) );
    $ngkpt[2] = ceil( $kgden/sqrt( $acell[6]**2 + $acell[7]**2 + $acell[8]**2 ) );
    $qe_data_files{'ngkpt'} = "$ngkpt[0] $ngkpt[1] $ngkpt[2]";
#    copy( "scf.out", "scf.out.$scfcount" );
#    copy( "scf.in", "scf.in.$scfcount" );
    $qe_data_files{'dft.startingpot'} = 'file';
    print "Re-running SCF: " . abs( $SCFEnergy - $oldSCFEnergy ) . "   $scfConv\n";
    $oldSCFEnergy = $SCFEnergy;
  }
  open OUT, ">scf.stat" or die "Failed to open scf.stat\n$!";
  print OUT "1\n";
  close OUT;
  print "SCF complete\n";

  print "Density PP Run\n";
  if( $obf == 1 )
  {  
    if( $qe_redirect )
    {
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'}  -npool $npool < pp.in > pp.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'} -npool $npool < pp.in > pp.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'}  -npool $npool < pp2.in > pp2.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'} -npool $npool < pp2.in > pp2.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
    } else
    {
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'}  -npool $npool -inp pp.in > pp.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'} -npool $npool -inp pp.in > pp.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'}  -npool $npool -inp pp2.in > pp2.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PP'} -npool $npool -inp pp2.in > pp2.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
    }
  }
  else
  {
    if( $qe_redirect )
    {  
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PP'}  -npool $npool < pp.in > pp.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_PP'} -npool $npool < pp.in > pp.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PP'}  -npool $npool < pp2.in > pp2.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_PP'} -npool $npool < pp2.in > pp2.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
    } else
    {
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PP'}  -npool $npool -inp pp.in > pp.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_PP'} -npool $npool -inp pp.in > pp.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PP'}  -npool $npool -inp pp2.in > pp2.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_PP'} -npool $npool -inp pp2.in > pp2.out 2>&1") == 0
         or die "Failed to run density stage for SCREENING\n";
    }
  }
  open OUT, ">den.stat" or die "Failed to open scf.stat\n$!";
  print OUT "1\n";
  close OUT;

  ## convert the density file to proper format
  print "Density conversion\n";
  system("$ENV{'OCEAN_BIN'}/qe2rhoofr.pl system.rho rhoofr" ) == 0 
    or die "Failed to convert density\n$!\n";

  print "Potential conversion\n";
  system("$ENV{'OCEAN_BIN'}/qe2rhoofr.pl system.pot potofr" ) == 0
    or die "Failed to convert potential\n$!\n";


  open STATUS, ">scf.stat" or die;
  print STATUS "1";
  close STATUS;


  # Find Fermi level and number of electrons
  my $fermi = 'no';
  my $nelectron = 'no';
  my $units;

  # First attempt to grab from outfile (works for 5.4 >= QE <= 6.2 (OLD_XML) )
  my $qe54_file = catfile( $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save", "data-file.xml" );
  my $qe62_file = catfile( $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save", "data-file-schema.xml" );
#  my $data_file = $qe_data_files{'work_dir'} . "/" . $qe_data_files{'prefix'} . ".save/data-file.xml";
  if( -e $qe54_file )
  {
    print "Looking for $qe54_file \n";
    open SCF, $qe54_file or die "Failed to open $qe54_file\n$!";
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
    if( $units =~ m/hartree/i )
    {
      $fermi *= 2;
    }
    elsif( $units =~ m/eV/i )
    {
      $fermi /= 13.60569253;
    }

    open OUT, '>', 'dftVersion' or die "Failed to open dftVersion for writing\n$!";
    print OUT "qe54\n";
    close OUT;
  }
  if( -e $qe62_file )  # Starting in QE6.5 it looks like both xml files are written 
  {
    print "$qe62_file\n";
    open SCF, $qe62_file or die "Failed to open $qe62_file\n$!";

    #Assume Hartree!
    my $highest;
    my $lowest = 'cow';
    while( my $scf_line = <SCF> )
    { 
      if( $scf_line =~ m/\<highestOccupiedLevel\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $highest = $1; 
      }
      elsif( $scf_line =~ m/\<lowestUnoccupiedLevel\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $lowest = $1;
      }
      elsif( $scf_line =~ m/\<fermi_energy\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $fermi = $1;
      }
      # We just average the two for spin=2 
      elsif( $scf_line =~ m/\<two_fermi_energies\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)\s+([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $fermi = ($1+$3)/2;
      }
      elsif( $scf_line =~ m/\<nelec\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $nelectron = $1;
      }
    }
    close SCF;
    if( $fermi eq 'no' )
    {
      if( $lowest eq 'cow' )
      { # Assumed Hartree
        $fermi = $highest * 2
      }
      else
      {
        $fermi = $highest + $lowest;
      }
    }
    else
    {
      # Move from Ha to Ry
      $fermi *= 2;
    }
    open OUT, '>', 'dftVersion' or die "Failed to open dftVersion for writing\n$!";
    print OUT "qe62\n";
    close OUT;
  }
  else # last shot
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
      elsif( $line  =~  m/Fermi energies are\s+([+-]?\d+\.?\d+)\s+([+-]?\d+\.?\d+)/ )
      {
        $fermi = ($1+$2)/2;
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

  my $eVfermi = $fermi * 13.60569253;
  print "Fermi level found at $eVfermi eV\n";

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

  my $nbands = $qe_data_files{'obf.nbands'};
  $nbands = $qe_data_files{'nbands'} if ( $nbands < 1 );
  
  if( $obf == 1 )
  {
    open QE, ">nscf.in" or die "Failed to open nscf.in\n$!";
    print QE "&control\n"
          .  "  calculation = 'nscf'\n"
          .  "  prefix = \'$qe_data_files{'prefix'}\'\n"
          .  "  pseudo_dir = \'$qe_data_files{'ppdir'}\'\n"
          .  "  outdir = \'$qe_data_files{'work_dir'}\'\n"
          .  "  wfcdir = \'$qe_data_files{'tmp_dir'}\'\n"
#          .  "  tstress = $qe_data_files{'dft.calc_stress'}\n"
#          .  "  tprnfor = $qe_data_files{'dft.calc_force'}\n"
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
    if( $qe_data_files{'qe_scissor'}  ne "" )
    {
      print QE "$qe_data_files{'qe_scissor'}\n";
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
          .  "  startingwfc = \'$qe_data_files{'dft.startingwfc'}\'\n"
          .  "  diagonalization = \'$qe_data_files{'dft.diagonalization'}\'\n"
          .  "/\n"
          .  "&ions\n"
          .  "/\n";

    open IN, "atompp" or die "$!";
    my $atompp;
    while (<IN>) { $atompp .= $_; }
    close IN;
    chomp $atompp;
    print QE "ATOMIC_SPECIES\n" . $atompp . "\n";

    if ($qe_data_files{'ibrav'} == 0) {
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
    elsif( $coord_type =~ m/bohr/ || $coord_type =~ m/xcart/ )
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
            . "$qe_data_files{'obkpt.ipt'} 0 0 0\n";
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

    if( $qe_redirect )
    {  
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} -npool $npool < nscf.in > nscf.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} -npool $npool < nscf.in > nscf.out 2>&1") == 0
         or die "Failed to run nscf stage for OBFs\n";
    } else
    {
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} -npool $npool -inp nscf.in > nscf.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_OBF_PW'} -npool $npool -inp nscf.in > nscf.out 2>&1") == 0
         or die "Failed to run nscf stage for OBFs\n";
    }
    print "NSCF complete\n";

    open OUT, ">nscf.stat" or die "Failed to open nscf.stat\n$!";
    print OUT "1\n";
    close OUT;

    print "Create Basis\n";
    open BASIS, ">basis.in" or die "Failed top open basis.in\n$!";
    print BASIS "&input\n"
             .  "  prefix = \'$qe_data_files{'prefix'}\'\n"
             .  "  outdir = \'$qe_data_files{'work_dir'}\'\n";
    unless( $qe_data_files{'tmp_dir'} =~ m/undefined/ )
    {
#      print "$qe_data_files{'tmp_dir'}\n";
      print BASIS "  wfcdir = \'$qe_data_files{'tmp_dir'}\'\n";
    }
    print BASIS "  trace_tol = $qe_data_files{'trace_tol'}\n"
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
            . "  outdir = \'$qe_data_files{'work_dir'}\'\n";
    unless( $qe_data_files{'tmp_dir'} =~ m/undefined/ )
    {
#      print "$qe_data_files{'tmp_dir'}\n";
      print HAM "  wfcdir = \'$qe_data_files{'tmp_dir'}\'\n";
    }
    print HAM "  updatepp = .false.\n"
            . "  ncpp = .true.\n"
            . "  nspin_ham = $qe_data_files{'nspin'}\n"
            . "/\n"
            . "K_POINTS\n"
            . "$ham_kpoints 0 0 0\n";
    close HAM;

    system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_ham_o.x  < ham.in > ham.out 2>&1") == 0
          or die "Failed to run shirley_ham_o.x\n$!";

  }
  else
  {

    my $bseDIR = sprintf("%03u%03u%03u", split( /\s+/,$qe_data_files{'nkpt'})); 
    mkdir $bseDIR unless ( -d $bseDIR );
    chdir $bseDIR;

    unlink "old" if( -e "old" );

    # kpts
    copy "../nkpt", "nkpt";
    copy "../qinunitsofbvectors.ipt", "qinunitsofbvectors.ipt";
    copy "../k0.ipt", "k0.ipt";
    copy "../dft.split", "dft.split";

    my $split_dft = 0;
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

    $qe_data_files{'prefix_shift'} = $qe_data_files{'prefix'} . "_shift";

    my $qeVersion;
#    mkdir "Out" unless ( -d "Out" );
 
    mkdir $qe_data_files{'work_dir'} unless( -d $qe_data_files{'work_dir'} );

    # This will loop back and do everything for prefix_shift if we have split
    my $repeat = 0;
    $repeat = 1 if( $split_dft );
    my $prefix = 'prefix';
    my $startRepeat = 0;

    $startRepeat = 1 if( $nscfRUN == 2 );

    for( my $i = $startRepeat; $i <= $repeat; $i++ )
    {
      my $savedir = catdir( $qe_data_files{'work_dir'}, $qe_data_files{$prefix} . ".save" ); 
#      mkdir "Out/$qe_data_files{$prefix}.save" unless ( -d "Out/$qe_data_files{$prefix}.save" );
      mkdir $savedir unless( -d $savedir );

      my $chargeDensity = catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save", 
                                   "charge-density.dat" );

      die "Couldn't find SCF charge density: $chargeDensity" unless( -e $chargeDensity );
      copy $chargeDensity, catfile( $savedir, "charge-density.dat");
#      copy "../Out/$qe_data_files{'prefix'}.save/charge-density.dat", "Out/$qe_data_files{$prefix}.save/charge-density.dat";

      if( -e catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save", "data-file.xml" ) )
      {
        copy "../Out/$qe_data_files{'prefix'}.save/data-file.xml", "Out/$qe_data_files{$prefix}.save/data-file.xml";
        $qeVersion = 54;
      }
      elsif( -e catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save", 
                         "data-file-schema.xml" ) )
      {
        copy catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save", 
                         "data-file-schema.xml" ), $savedir;
        $qeVersion = 62;
      }
      else
      {
        die "Couldn't find data-file or data-file-schema\n";
      }

      if( $qe_data_files{'nspin'} == 2 )
      {
        copy "../Out/$qe_data_files{'prefix'}.save/spin-polarization.dat", 
             "Out/$qe_data_files{$prefix}.save/spin-polarization.dat";
      }

      if( $qe_data_files{'ldau'}  ne "" )
      {
        # Starting w/ QE-6.0 this is the DFT+U info from the SCF
        if( -e "../Out/$qe_data_files{'prefix'}.save/occup.txt" )
        {
          copy "../Out/$qe_data_files{'prefix'}.save/occup.txt", "Out/$qe_data_files{$prefix}.save/occup.txt";
        }
        # QE 4.3-5.x
        elsif( -e "../Out/$qe_data_files{'prefix'}.occup" )
        {
          copy "../Out/$qe_data_files{'prefix'}.occup", "Out/$qe_data_files{$prefix}.occup";
        }
      }
      foreach my $bonusFile ( @SCFBonus )
      {
        my $tempFile = "../Out/$qe_data_files{'prefix'}.save/$bonusFile";
        copy $tempFile, "Out/$qe_data_files{$prefix}.save/$bonusFile" if( -e $tempFile );
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
    $qe_data_files{'calctype'} = 'nscf';
    $qe_data_files{'dft.startingpot'} = 'file';
#    # if have exact exchange flip back to scf
#    print "$qe_data_files{'dft.functional'}\n";

    # some of this needs to be moved up
    foreach( @exx )
    {
#      print "$_\n";
      if( $qe_data_files{'dft.functional'} =~ m/$_/i )
      {
        $qe_data_files{'calctype'} = 'scf';
        $qe_data_files{'nscfEXX'} = 1;
        last;
      }
    }
    $qe_data_files{'nosym'} = '.true.';
    $qe_data_files{'noinv'} = '.true.';
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

    $qe_data_files{'print kpts'} = $kpt_text;
    # QE behaves cnoverges incorrectly if only give occupied states
    if( $split_dft && $qe_data_files{ "occopt" } == 1 ) 
    {
      open NEL, "../nelectron" or die "Failed top open ../nelectron for reading\n$!";
      my $nelectron = <NEL>;
      close NEL;
      my $tempBand = $nelectron / 2 + 1;
      $tempBand++ if( $tempBand%2 == 1 );
      $qe_data_files{'print nbands'} = $tempBand;
    }
    elsif( $split_dft == 0)
#    unless( $split_dft ) 
    {
      $qe_data_files{'print nbands'} = $qe_data_files{'nbands'};
    }

    &print_qe( $QE, %qe_data_files );

    close $QE;

    my $npool = 1;
    my $ncpus = 1;
    open INPUT, "../pool_control" or die;
    while (<INPUT>)
    {
      if( $_ =~ m/^nscf\s+(\d+)/ )
      {
        $npool = $1;
      }
      elsif( $_ =~ m/^total\s+(\d+)/ )
      {
        $ncpus = $1;
      }

    }
    close INPUT;

    print "TEST\n";
    print $qe_data_files{'dft.ndiag'} . "\n";
    if( $qe_data_files{'dft.ndiag'} =~ m/(-?\d+)/ )
    {
      if( $1 > 0 )
      {
        $qe_data_files{'dft.ndiag'} = $1;
      }
      else
      {
        $qe_data_files{'dft.ndiag'} = $ncpus;
      }
    }
    else
    {
      $qe_data_files{'dft.ndiag'} = $ncpus;
    }


    my $qeCommandLine = "-ndiag $qe_data_files{'dft.ndiag'} -npool $npool";

    print "BSE NSCF Run\n";
    if( $qe_redirect )
    {  
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < nscf.in > nscf.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < nscf.in > nscf.out 2>&1") == 0
          or die "Failed to run nscf stage for BSE wavefunctions\n";
    } else
    {
      print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp nscf.in > nscf.out 2>&1\n";
      system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp nscf.in > nscf.out 2>&1") == 0
          or die "Failed to run nscf stage for BSE wavefunctions\n";
    }


    if( $split_dft && $nscfRUN == 2 )
    {
      print "Unoccupied states re-used from previous calculation\n";
    }
    elsif( $split_dft )
    {
      open my $QE, ">nscf_shift.in" or die "Failed to open nscf_shift.in\n$!";

      $prefix = $qe_data_files{'prefix'};
      $qe_data_files{'prefix'} = $qe_data_files{'prefix_shift'};

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
      $qe_data_files{'print kpts'} = $kpt_text;
      $qe_data_files{'print nbands'} = $qe_data_files{'nbands'};

      &print_qe( $QE, %qe_data_files );

      close $QE;

      if( $qe_redirect )
      {  
        print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < nscf_shift.in > nscf_shift.out 2>&1\n";
        system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < nscf_shift.in > nscf_shift.out 2>&1") == 0
            or die "Failed to run nscf stage for shifted BSE wavefunctions\n";
      } else
      {
        print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp nscf_shift.in > nscf_shift.out 2>&1\n";
        system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp nscf_shift.in > nscf_shift.out 2>&1") == 0
            or die "Failed to run nscf stage for shifted BSE wavefunctions\n";
      }

      $qe_data_files{'prefix'} = $prefix;

    }



    open OUT, ">nscf.stat" or die "Failed to open nscf.stat\n$!";
    print OUT "1\n";
    close OUT;
    print "BSE NSCF complete\n";

    ## find the top of the valence bance
    if( $qeVersion == 54 )
    {
      system("$ENV{'OCEAN_BIN'}/qeband.pl") == 0
         or die "Failed to count bands\n$!\n";
    }
    elsif( $qeVersion == 62 )
    {
      system("$ENV{'OCEAN_BIN'}/qe62band.pl") == 0
         or die "Failed to count bands\n$!\n";
    }
    else
    { die "qeVersion wasn't set\n"; }

    open IN, "brange.stub" or die;
    open OUT, ">brange.ipt" or die;
    while(<IN>)
    {
      print OUT $_;
    }
    print OUT "$qe_data_files{'nbands'}\n";
    close IN;
    close OUT;

    copy "nkpt", "kmesh.ipt";

    chdir "../";

    open OUT, ">", "bse.stat" or die;
    print OUT "1\n";
    close OUT;
  }
}
else
{
  my $bseDIR = sprintf("%03u%03u%03u", split( /\s+/,$qe_data_files{'nkpt'}));
  die "Problem with $bseDIR\n" unless( chdir $bseDIR );
  open OUT, ">", "old";
  print OUT "old\n";
  close OUT;
  chdir "../";
}

if( $obf == 0 && $run_screen == 1 )
{
  
#  my $bseDIR = sprintf("%03u%03u%03u", split( /\s+/,$qe_data_files{'screen.nkpt'}));
  my $bseDIR = "SCREEN";
  mkdir $bseDIR unless ( -d $bseDIR );
  chdir $bseDIR;

  unlink "old";

  mkdir "Out" unless ( -d "Out" );
  mkdir "Out/$qe_data_files{'prefix'}.save" unless ( -d "Out/$qe_data_files{'prefix'}.save" );

#  copy "../Out/$qe_data_files{'prefix'}.save/charge-density.dat", "Out/$qe_data_files{'prefix'}.save/charge-density.dat";
#  copy "../Out/$qe_data_files{'prefix'}.save/data-file.xml", "Out/$qe_data_files{'prefix'}.save/data-file.xml";


  my $savedir = catdir( $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save" );
#      mkdir "Out/$qe_data_files{$prefix}.save" unless ( -d "Out/$qe_data_files{$prefix}.save" );
  mkdir $savedir unless( -d $savedir );

  my $qeVersion;

  my $chargeDensity = catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save",
                               "charge-density.dat" );

  die "Couldn't find SCF charge density: $chargeDensity" unless( -e $chargeDensity );
  copy $chargeDensity, catfile( $savedir, "charge-density.dat");
#      copy "../Out/$qe_data_files{'prefix'}.save/charge-density.dat", "Out/$qe_data_files{$prefix}.save/charge-density.dat";

  if( -e catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save", "data-file.xml" ) )
  {
    copy "../Out/$qe_data_files{'prefix'}.save/data-file.xml", "Out/$qe_data_files{'prefix'}.save/data-file.xml";
    $qeVersion = 54;
  }
  elsif( -e catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save",
                     "data-file-schema.xml" ) )
  {
    copy catfile( updir(), $qe_data_files{'work_dir'}, $qe_data_files{'prefix'} . ".save",
                     "data-file-schema.xml" ), $savedir;
    $qeVersion = 62;
  }
  else
  {
    die "Couldn't find data-file or data-file-schema\n";
  }


  if( $qe_data_files{'nspin'} == 2 )
  {
    copy "../Out/$qe_data_files{'prefix'}.save/spin-polarization.dat", 
         "Out/$qe_data_files{'prefix'}.save/spin-polarization.dat";
  }

  if( $qe_data_files{'ldau'}  ne "" )
  {
    # Starting w/ QE-6.0 this is the DFT+U info from the SCF
    if( -e "../Out/$qe_data_files{'prefix'}.save/occup.txt" )
    {
      copy "../Out/$qe_data_files{'prefix'}.save/occup.txt", "Out/$qe_data_files{'prefix'}.save/occup.txt";
    }
    # QE 4.3-5.x
    elsif( -e "../Out/$qe_data_files{'prefix'}.occup" )
    {
      copy "../Out/$qe_data_files{'prefix'}.occup", "Out/$qe_data_files{'prefix'}.occup";
    }
  }

  foreach my $bonusFile ( @SCFBonus )
  {
    my $tempFile = "../Out/$qe_data_files{'prefix'}.save/$bonusFile";
    copy $tempFile, "Out/$qe_data_files{'prefix'}.save/$bonusFile" if( -e $tempFile );
  }


  # kpts
  copy "../screen.nkpt", "nkpt";
  copy "../screen.k0", "k0.ipt";


  # QINB is 0 for screening
  open OUT, ">qinunitsofbvectors.ipt" or die;
  print OUT "0.0 0.0 0.0\n";
  close OUT;
  open OUT, ">core" or die;
  print OUT "1\n";
  close OUT;
  system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 or die "KGEN.X Failed\n";

  open my $QE, ">nscf.in" or die "Failed to open nscf.in\n$!";


  $qe_data_files{'calctype'} = 'nscf';
  $qe_data_files{'dft.startingpot'} = 'file';
  # if have exact exchange flip back to scf
  foreach( @exx )
  {
    if( $qe_data_files{'dft.functional'} =~ m/$_/i )
    {
      $qe_data_files{'calctype'} = 'scf';
    }
  }

  $qe_data_files{'nosym'} = '.true.';
  $qe_data_files{'noinv'} = '.true.';

  my $gamma = 1;

  open IN, "nkpt" or die "Failed to open nkpt (screening)\n$!";
  <IN> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse screen.nkpt.\n";
  $gamma *= 0 if ( $1 * $2 * $3 > 1 ); 
  close IN;

  open IN, "k0.ipt"  or die "Failed to open k0.ipt (screening)\n$!";
  <IN> =~ m/(\S+)\s+(\S+)\s+(\S+)/ or die "Failed to parse k0.ipt\n";
  $gamma *= 0 if( abs($1) > 0.000001 || abs($2) > 0.000001 || abs($3) > 0.000001 );
  close IN;

  my $kpt_text;
  if( $gamma == 1 ) 
  {
    $kpt_text = "K_POINTS gamma\n";
  }
  else
  {
    $kpt_text = "K_POINTS crystal\n";
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
  }
  
  $qe_data_files{'print kpts'} = $kpt_text;
  $qe_data_files{'print nbands'} = $qe_data_files{'screen.nbands'};

  &print_qe( $QE, %qe_data_files );

  close QE;

  my $npool = 1;
  my $ncpus = 1;
  open INPUT, "../pool_control" or die;
  while (<INPUT>)
  {
    if( $_ =~ m/^screen\s+(\d+)/ )
    {
      $npool = $1;
    }
    elsif( $_ =~ m/^total\s+(\d+)/ )
    {
      $ncpus = $1;
    }

  }
  close INPUT;

  print "TEST\n";
  print $qe_data_files{'dft.ndiag'} . "\n";
  if( $qe_data_files{'dft.ndiag'} =~ m/(-?\d+)/ )
  {
    if( $1 > 0 )
    {
      $qe_data_files{'dft.ndiag'} = $1;
    }
    else
    {
      $qe_data_files{'dft.ndiag'} = $ncpus;
    }
  }
  else
  {
    $qe_data_files{'dft.ndiag'} = $ncpus;
  }


  my $qeCommandLine = "-ndiag $qe_data_files{'dft.ndiag'} -npool $npool";

  print "Screening NSCF Run\n";
  if( $qe_redirect )
  {  
    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < nscf.in > nscf.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine < nscf.in > nscf.out 2>&1") == 0
        or die "Failed to run nscf stage for SCREENing wavefunctions\n";
  } else
  {
    print  "$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp nscf.in > nscf.out 2>&1\n";
    system("$para_prefix $ENV{'OCEAN_ESPRESSO_PW'} $qeCommandLine -inp nscf.in > nscf.out 2>&1") == 0
        or die "Failed to run nscf stage for SCREENing wavefunctions\n";
  }
  open OUT, ">nscf.stat" or die "Failed to open nscf.stat\n$!";
  print OUT "1\n";
  close OUT;
  print "Screening NSCF complete\n";

  ## find the top of the valence bands
  if( $qeVersion == 54 )
  {
    system("$ENV{'OCEAN_BIN'}/qeband.pl") == 0
       or die "Failed to count bands\n$!\n";
  }
  elsif( $qeVersion == 62 )
  {
    system("$ENV{'OCEAN_BIN'}/qe62band.pl") == 0
       or die "Failed to count bands\n$!\n";
  }
  else
  { die "qeVersion wasn't set\n"; }

  # Figure out Gamma point usage within QE
  my $dataFileName;
  my $workdir = $qe_data_files{"work_dir"};
  $workdir =~ s/\'//g;
  $workdir =~ s/^\.//;
  $workdir =~ s/^\///;
  my $prefix = $qe_data_files{"prefix"};
  if(  $qeVersion == 54 )
  {
    $dataFileName = $workdir . '/' . $prefix . ".save/data-file.xml";
  }
  else
  {
    $dataFileName = $workdir . '/' . $prefix . ".save/data-file-schema.xml";
  }
  
  open IN, $dataFileName or die "Failed to open $dataFileName\n$!\n";
  open GAMMA, ">", "gamma" or die "$!";

  while (my $line = <IN>)
  {

    if ( $line =~ m/GAMMA_ONLY/i )
    {
      chomp($line);
      $line .= <IN>;
      chomp($line);
      $line .= <IN>;
      chomp($line);

      if( $line =~ m/>([ft])(alse)?(rue)?\s*</i )
      {
        my $gamma = $1;
        print $gamma . "\n";
        print GAMMA $gamma . "\n";
      }
      else
      {
        print "Nope!\n$line\n";
        print GAMMA "F" . "\n";
      }

    last;
    }
  }
  close IN;
  close GAMMA;

  open IN, "brange.stub" or die;
  open OUT, ">brange.ipt" or die;
  while(<IN>)
  {
    print OUT $_;
  }
  print OUT "$qe_data_files{'screen.nbands'}\n";
  close IN;
  close OUT;

  copy "nkpt", "kmesh.ipt";

  chdir "../";

  open OUT, ">", "screen.stat" or die;
  print OUT "1\n";
  close OUT;
}
else
{
  my $bseDIR = "SCREEN";
  if( -d $bseDIR )
  {
    die "Problem with $bseDIR\n" unless( chdir $bseDIR );
    open OUT, ">", "old";
    print OUT "old\n";
    close OUT;
    chdir "../";
  }
}

# For occopt = 1, the SCF run only gives the highest occupied
# With a sparse k-point grid this doesn't give a good position for the Fermi
# So, if we have re-run any segment then figure out a better Fermi level
# We do this by taking the highest occupied from SCREEN and BSE and the lowest unoccupied
# and then setting the Fermi to be the midpoint
if( $qe_data_files{ "occopt" } == 1 && ( $RunESPRESSO + $nscfRUN + $run_screen > 0 ) )
{
  print "Fixing Fermi level for occopt=1, insulating system\n";

  my $valenceE;
  my $conductionE;

  if( -e "SCREEN/nscf.out" )
  {
    open IN, "SCREEN/nscf.out" or die "Failed to open SCREEN/nscf.out\n$!";
    while( my $line = <IN> )
    {
      if( $line =~ m/highest occupied, lowest unoccupied level\s\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/ )
      {
        $valenceE = $1;
        $conductionE = $2;
        print "$valenceE\t$conductionE\n";
        last;
      }
    }
    close IN;
  }
  my $bseDIR = sprintf("%03u%03u%03u", split( /\s+/,$qe_data_files{'nkpt'}));
  if( -e "$bseDIR/nscf_shift.out" )
  {
    open IN, "$bseDIR/nscf_shift.out" or die "Failed to open $bseDIR/nscf_shift.out\n$!";
    while( my $line = <IN> )
    {
      if( $line =~ m/highest occupied, lowest unoccupied level\s\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/ )
      {
        print $line;
        my $v = $1;
        my $c = $2;
        unless( defined( $valenceE ) )
        {
          $valenceE = $v;
          $conductionE = $c;
        }
        else
        {
          $valenceE = $v if( $v > $valenceE );
          $conductionE = $c if( $c < $conductionE );
        }
        last;
        print "$valenceE\t$conductionE\n";
      }
    }
    close IN;
  }
  if( -e "$bseDIR/nscf.out" )
  {
    open IN, "$bseDIR/nscf.out" or die "Failed to open $bseDIR/nscf.out\n$!";
    while( my $line = <IN> )
    {
      if( $line =~ m/highest occupied, lowest unoccupied level\s\S+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/ )
      {
        print $line;
        my $v = $1;
        my $c = $2;
        unless( defined( $valenceE ) )
        {
          $valenceE = $v;
          $conductionE = $c;
        }
        else
        { 
          $valenceE = $v if( $v > $valenceE );
          $conductionE = $c if( $c < $conductionE );
        }
        last;
        print "$valenceE\t$conductionE\n";
      }
    }
    close IN;
  }
  else
  {
    print "Couldn't fine $bseDIR/nscf.out\n";
  }

  unless( defined( $valenceE ) )
  {
    print "Failed to correct Fermi level!\nLikely DFT runs didn't finish correctly!!";
  }
  else
  {
    if( $valenceE > $conductionE )
    {
      print "WARNING!!!! Highest occupied is greater than lowest unoccupied!\n";
      print "Likely you specified metal = .false. for a metallic system or your structure is incorrect\n";
      print "OCEAN will continue, but the results are probably bad!\nWARNING!!!!!\n";
    }
    else
    {
      print "$valenceE\t$conductionE\n";
      my $eVfermi = ( $valenceE + $conductionE ) / 2;
      my $fermi = $eVfermi / 13.60569253;
      print "Fermi level found at $eVfermi eV\n";

      open FERMI, ">efermiinrydberg.ipt" or die "Failed to open efermiinrydberg\n$!";
      print FERMI "$fermi\n";
      close FERMI;

    }
  }
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
        .  "  disk_io = 'low'\n"
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
        .  "  nosym = $inputs{'nosym'}\n"
        .  "  noinv = $inputs{'noinv'}\n";
#        .  "  nosym = .true.\n"
#        .  "  noinv = .true.\n";
  unless( $inputs{'dft.functional'} =~ m/none/ )
  {
    print $fh "  input_dft = \'$inputs{'dft.functional'}\'\n";
    foreach( @exx )
    {
      if( $inputs{'dft.functional'} =~ m/$_/i )
      {
        print $fh "  nqx1 = $inputs{'nqx1'}, nqx2 = $inputs{'nqx2'}, nqx3 = $inputs{'nqx3'}\n";
        last;
      }
    }
  }
    


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
        .  "  startingwfc = \'$inputs{'dft.startingwfc'}\'\n"
        .  "  startingpot = \'$inputs{'dft.startingpot'}\'\n"
        .  "  diagonalization = \'$inputs{'dft.diagonalization'}\'\n";
  if( $inputs{'print nbands'} > 100 && $inputs{'calctype'} =~ m/nscf/i )
  {
    print $fh "  diago_david_ndim = 2\n";
  }
  if( $inputs{'nscfEXX'} == 1 )
  {
    # Since (at the moment) we are loading the SCF density
    #  don't converge the density for the first iteration w/o EXX
    print $fh "  adaptive_thr = .true., conv_thr_init = 1\n";
  }
  print $fh "/\n"
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
  elsif( $coord_type =~ m/bohr/ || $coord_type =~ m/cart/ )
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

  
