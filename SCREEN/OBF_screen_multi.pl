#!/usr/bin/perl
# Copyright (C) 2015 OCEAN collaboration
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
  $0 =~ m/(.*)\/OBF_screen_multi\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################

#FAKE INPUTS FOR NOW

#my $para_prefix = "mpirun -n 16 ";
#my $trace_tolerance = "5.0d-15";

# Spline for Hamiltonian
#my $ham_kpoints = "4 4 4";

my $band_start = 1;
# A value of minus one should cause the max number to be used
# (ie. band max = number of OBFs )
my $band_stop  = -1; #800;

#my $screen_nkpt = "2 2 2";



# Step 1: Create support files
my @CommonFiles = ("znucl", "paw.hfkgrid", "paw.fill", "paw.opts", "pplist", "paw.shells", 
                   "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution", 
                   "epsilon", "k0.ipt", "ibase", "scfac", "para_prefix", "tmp_dir", 
                   "paw.nbands", "core_offset", "paw.nkpt", "pool_control", "ham_kpoints", 
                   "cnbse.rad", "dft", "avecsinbohr.ipt", 
                   "screen.grid.rmax", "screen.grid.nr", "screen.grid.ang", "screen.grid.lmax", 
                   "screen.grid.nb", "screen.final.rmax", "screen.final.dr", "screen.model.dq", "screen.model.qmax" );

my @DFTFiles = ("rhoofr", "nscf.out", "system.rho.dat", "efermiinrydberg.ipt");
my @ExtraFiles = ("Pquadrature" );

foreach(@ExtraFiles)
{
  `cp $ENV{'OCEAN_BIN'}/$_ .` == 0 or die;
}

foreach (@DFTFiles)
{
  `cp ../DFT/$_ .` == 0 or die "Failed to get $_ from DFT/\n";
}

foreach (@CommonFiles) {
  `cp ../Common/$_ .` == 0 or die "Failed to get $_ from Common/\n";
}

my $tmpdir = "undefined";
if( open TMPDIR, "tmp_dir" )
{
  $tmpdir = <TMPDIR>;
  chomp( $tmpdir );
  close TMPDIR;
}

open IN, "epsilon" or die "Failed to open epsilon\n$!";
my $epsilon = <IN>;
chomp $epsilon;
close IN;

#open IN, "ibase" or die "Failed to open ibase\n$!";
#my $ibase;
#while( <IN> ) { $ibase .= $_; }

my $pool_size = 1;
open INPUT, "pool_control" or die;
while (<INPUT>)
{
  if( $_ =~ m/interpolate paw\s+(\d+)/ )
  {
    $pool_size = $1;
    last;
  }
}
close INPUT;

my $para_prefix = "";
if( open PARA_PREFIX, "para_prefix" )
{
  $para_prefix = <PARA_PREFIX>;
  chomp($para_prefix);
  close( PARA_PREFIX);
} else
{
  print "Failed to open para_prefix. Error: $!\nRunning serially\n";
  $pool_size = 1;
}

#if( -e "../DFT/ham_kpoints" )
#{
#	`cp ../DFT/ham_kpoints .`;
#	$ham_kpoints = `cat ham_kpoints`;
#	chomp($ham_kpoints);
#}
open IN, "ham_kpoints" or die;
<IN> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse ham_kpoints\n$_";
my $ham_kpoints = "$1  $2  $3";



my $screen_nkpt = `cat paw.nkpt`;
chomp($screen_nkpt);


`ln -s ../DFT/Out .`;

#my $fermi = 0;
#open SCF, "nscf.out" or die "$!";
#while( my $line = <SCF> )
#{
#  if( $line  =~  m/the Fermi energy is\s+([+-]?\d\S+)/ )
#  {
#    $fermi = $1;
#    print "Fermi level found at $fermi eV\n";
#    last;
#  }
#}
#$fermi = $fermi/13.60569252;
#`echo "$fermi" > efermiinrydberg.ipt`;
open EFERMI, "efermiinrydberg.ipt" or die "$!\nFailed to open eferiinrydberg.ipt\n";
my $fermi = <EFERMI>;
chomp( $fermi);
close EFERMI;
print "Fermi level found at " . $fermi*13.60569252 . " eV\n";


system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0
    or die "Failed to run bvecs.pl\n";

`tail -n 1 rhoofr > nfft` if( -e "rhoofr" );
system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0
  or die "Failed to run rhoofg.x\n";
`wc -l rhoG2 > rhoofg`;
`sort -n -k 6 rhoG2 >> rhoofg`;


print "Running PAW Setup\n";
system("$ENV{'OCEAN_BIN'}/pawsetup.x") == 0 or die "$!\nFailed to run pawsetup.x\n";

if( -e "$ENV{'OCEAN_BIN'}/mpi_avg.x" )
{
  print "Running mpi_avg.x\n";
  system("$para_prefix $ENV{'OCEAN_BIN'}/mpi_avg.x") == 0 or die "$!\nFailed to run mpi_avg.x\n";
}
else
{
  print "Running avg.x\n";
  system("$ENV{'OCEAN_BIN'}/avg.x") == 0 or die "$!\nFailed to run avg.x\n";
}


`ln -sf ../OPF/zpawinfo zpawinfo`;
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

open NBAND, "paw.nbands" or die "$!";
if( <NBAND> =~ m/(\d+)/ )
{
  $band_stop = $1;
} else
{
  print "Failed to parse paw.nbands, running with $band_stop\n";
}
close( NBAND );

# Step 4: Loop over every core site, building W(r)

# Prep input file
open BOFR, ">bofr.in" or die "$!\nFailed to open bofr.in for writing\n";
print BOFR "&input\n" .
          "  prefix = 'system_opt'\n" .
          "  outdir = './Out'\n" .
          "  wfcdir = '$tmpdir'\n" .
          "  updatepp = .false.\n" .
          "  ncpp = .true.\n" .
          "  calculation = 'ocean_bofr_multi'\n" .
          "/\n" .
          " K_POINTS\n" .
          "$ham_kpoints 0 0 0\n";
close BOFR;

open BUILDER, ">builder.in" or die "$!\nFailed to open builder.in for writing\n";
print BUILDER "&input\n" .
          "  prefix = 'system_opt'\n" .
          "  outdir = './Out'\n" .
          "  band_subset = $band_start  $band_stop\n" .
          "/\n" .
          " K_POINTS\n" .
          " automatic\n $screen_nkpt 1 1 1 \n";
close BUILDER;

open HFINLIST, "hfinlist" or die "Failed to open hfinlist\n";

my $rad;
my $edgename;
my $hfinline; my $ppfilename; my $znucl; my $nnum; my $lnum; my $elname; my $elnum;

my $nedges = `cat nedges`;
chomp($nedges);


##### Grab inputs that control screening grid
open IN, "screen.grid.rmax" or die "Failed to open screen.grid.rmax\n$!";
<IN> =~ m/([+-]?\d+\.?\d*)/ or die "Failed to parse screen.grid.rmax\n";
my $grid_rmax = $1;
close IN;
$grid_rmax = 8 if( $grid_rmax < 0 );

open IN, "screen.grid.nr" or die "Failed to open screen.grid.nr\n$!";
<IN> =~ m/([+-]?\d+)/ or die "Failed to parse screen.grid.nr\n";
my $grid_nr = $1;
close IN;
$grid_nr = 24 if( $grid_nr < 0 );

open IN, "screen.grid.ang" or die "Failed to open screen.grid.ang\n$!";
<IN> =~ m/(\S+)\s+(\d+)/ or die "Failed to parse screen.grid.ang\n";
my $grid_ang_type = $1;
my $grid_ang_max = $2;
close IN;
if( lc($grid_ang_type) ne "lebdev" )
{
  print "Un-supported grid type requested: $grid_ang_type.\nOCEAN will default to Lebdev grid\n";
  $grid_ang_type = "lebdev"
}
my $specpnt;
if( lc($grid_ang_type) eq "lebdev" )
{
  $specpnt = sprintf("specpnt.%i",$grid_ang_max);
}
else
{
  die "Unsupported angular grid requested\n";
}

unless( -e "$ENV{'OCEAN_BIN'}/$specpnt" )
{
  die "Requested angular grid not found in OCEAN installation: $ENV{'OCEAN_BIN'}/$specpnt\n";
}
print "$specpnt\n";
copy "$ENV{'OCEAN_BIN'}/$specpnt", "specpnt" or die "$!";

open IN, "screen.grid.lmax" or die "Failed to open screen.grid.lmax\n$!";
<IN> =~ m/([+-]?\d+)/ or die "Failed to parse screen.grid.lmax\n";
my $grid_lmax = $1;
close IN;
if( $grid_lmax < 0 )
{
  print "WARNING: Negative LMAX requested for screening calculation. Using 0 instead\n"
}
if( $grid_lmax > $grid_ang_max/2 )
{
  $grid_lmax = sprintf("%i",($grid_ang_max/2) );
  print "WARNING: Requested LMAX is too large for grid. Using $grid_lmax instead\n";
}

# For now!
if( $grid_lmax != 0 ) 
{
  print "WARNING: Only LMAX = 0 is programmed for the screening\n";
  $grid_lmax = 0;
}
# \For now

open IN, "screen.grid.nb" or die "Failed to open screen.grid.nb\n$!";
<IN> =~ m/([+-]?\d+)/ or die "Failed to parse screen.grid.nb\n";
my $grid_nb = $1;
close IN;
if( $grid_nb < 0 )
{
  $grid_nb = ($grid_nr - 1);
  print "WARNING: Negative number of basis functions requested. Using $grid_nb.\n";
}
if( $grid_nb > ($grid_nr - 1) ) 
{
  $grid_nb = ($grid_nr - 1);
  print "WARNING: For screening calc the number of radial grid points \n"
      . "  must exceed the number of radial basis functions.\n"
      . "  Reducing this to $grid_nb\n";
}

open IN, "screen.final.rmax" or die "Failed to open screen.final.rmax\n$!";
<IN> =~ m/(\d+\.?\d?)/ or die "Failed to parse screen.final.rmax\n";
my $final_rmax = $1;
close IN;

open IN, "screen.final.dr" or die "Failed to open screen.final.dr\n$!";
<IN> =~ m/(\d+\.?\d*)/ or die "Failed to parse screen.final.dr\n";
my $final_dr = $1;
close IN;
print "$final_dr  $final_rmax\n";
my $final_nr = sprintf("%.i",$final_rmax/$final_dr);
print "$final_dr  $final_nr\n";

open IN, "screen.model.dq" or die "Failed to open screen.model.dq\n$!";
<IN> =~ m/(\d+\.?\d*)/ or die "Failed to parse screen.model.dq\n";
my $model_dq = $1;
close IN;

open IN, "screen.model.qmax" or die "Failed to open screen.model.qmax\n$!";
<IN> =~ m/(\d+\.?\d*)/ or die "Failed to parse screen.model.qmax\n";
my $model_qmax = $1;
close IN;

my $ibase = "$final_dr  $final_nr\n$model_dq  $model_qmax\n";
print $ibase;

open IBASE, ">ibase" or die "Failed to open to open ibase for writing\n";
print IBASE $ibase;
close IBASE;

#### Have everything needed for screening grid




open MKRB_CONTROL, ">mkrb_control" or die "Failed to open mkrb_control\n$!\n";
#print MKRB_CONTROL "8 25\n";
print MKRB_CONTROL "$grid_rmax $grid_nr\n";
print MKRB_CONTROL "$nedges\n";


my $temp_edgename;

while ($hfinline = <HFINLIST>) {
  chomp $hfinline;
  print $hfinline . "\n";
  ($hfinline =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)/) 
                  or die "Malformed hfinlist\t$1 $2 $3 $4 $5 $6\n";
  $ppfilename = $1;
  $znucl = $2;
  $nnum = $3;
  $lnum = $4;
  $elname = $5;
  $elnum = $6;

  print MKRB_CONTROL "$elname $elnum\n";
  $temp_edgename = sprintf("z%03in%02il%02i",$znucl, $nnum, $lnum);

}
close MKRB_CONTROL;

open BOFR_CONTROL, ">bofr_tau_control" or die "Failed to open bofr_tau_control\n$!\n";
print BOFR_CONTROL "$nedges 0\n";
close BOFR_CONTROL;

seek(HFINLIST,0,0);

  # Step 4.1: Build the radial grid for this element
print "mkrbfile.x\n";
system("$ENV{'OCEAN_BIN'}/mkrbfile_mult.x") == 0 
  or die "Failed to run mkrbfile_mult.x\n$!\n";

  # Step 4.2: Project the basis functions onto this radial grid
print "bofr\n";
system("$para_prefix $ENV{'OCEAN_BIN'}/shirley_ham_o.x < bofr.in > bofr.out") == 0 
        or die "$!\nFailed to run shirley_ham from bofr.in\n";
#JTV only works for all the same core hole potential right now
#JTV and all the same radius too
foreach $rad (@rads)
{
#  my $temp_rad = sprintf("%03.2f",$rads[0]);
  my $temp_rad = sprintf("%03.2f",$rad);
  open VCx, "zpawinfo/vc_bare${temp_edgename}R${temp_rad}" or die "Failed to open vc_bare${temp_edgename}R${temp_rad}\n$!";
  open VPERT, ">vpert.$rad" or die;
  while (<VCx>){}
  my $vpert_length = $.;
  seek(VCx,0,0);
  print VPERT "$vpert_length\n";
  while (<VCx>)
  {
    print VPERT $_;
  }
  close VCx;
  close VPERT;


}
unlink "vpert" if (-e "vpert");
symlink "vpert.$rads[0]", "vpert";

#`cp zpawinfo/vcxxxxx${temp_edgename}R${temp_rad} ./tmp`;
#`wc tmp > vpert`;
#`cat tmp >> vpert`;

print "ocean_builder.x\n";
system("$para_prefix $ENV{'OCEAN_BIN'}/ocean_builder.x  $pool_size  < builder.in 1> builder.out 2> builder.err") == 0
        or die "$!\nFailed to run ocean_builder.x\n";

my $itau = 0;

while ($hfinline = <HFINLIST>) {

  $itau++;
  chomp $hfinline;
  print $hfinline . "\n";
  $hfinline =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\d+)/ 
        or die "Malformed hfinlist\t$1 $2 $3 $4 $5 $6\n";
  $ppfilename = $1;
  $znucl = $2;
  $nnum = $3;
  $lnum = $4;
  $elname = $5;
  $elnum = $6;


  $edgename = sprintf("z%2s%04i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
  print "$edgename\n";
  unless( -d $edgename )
  {
    mkdir $edgename or die "Failed to make dir $edgename\n";
  }

  my $avden =  sprintf("avg%2s%04i",$elname,$elnum);
  copy( $avden, "avden" ) or die "Failed to copy density $avden\n$!";

  my $edgename2 = sprintf("z%03in%02il%02i",$znucl, $nnum, $lnum);


# Step 5: For each core site, loop over radius
#           This radius is for the neutralizing charge
  foreach $rad (@rads) {
    my $fullrad = sprintf("%03.2f",$rad);

    chdir "$edgename" or die "$!";
    foreach my $t_dir ( "zRXT${fullrad}", "zRXF${fullrad}", "zRXS${fullrad}" )
    {
      unless( -d $t_dir )
      {
        mkdir "$t_dir" or die "$!";
      }
    }
    unlink "zR${fullrad}" if ( -e "zR${fullrad}" );
    symlink "zRXT${fullrad}",  "zR${fullrad}" ;
    chdir "../";
    

    my $ximat_name = sprintf("ximat%04i",$itau);
    unlink "ximat" if ( -e "ximat" );
    symlink $ximat_name, "ximat";
#    `ln -sf $ximat_name ximat`;
    
    open IPT, ">ipt" or die "Failed to open ipt\n$!";
    # Number of basis functions for chi
    print IPT "$grid_nb\n";
#    print IPT "24\n";
    close IPT;


    unlink "vpert" if( -e "vpert");
    symlink "vpert.$rad", "vpert";
    `$ENV{'OCEAN_BIN'}/xipps.x < ipt`;
    move( "ninduced", "nin" );
#    `mv ninduced nin`;

    open IPT, ">ipt" or die "Failed to open ipt\n$!";
    print IPT "$fullrad\n";
#    print IPT $ibase;
#    print IPT "$epsilon\n";
    close IPT;
    `$ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
    move( "reopt", "rom" );
#    `mv reopt rom`;

#    open IPT, ">ipt" or die "Failed to open ipt\n$!";
#    print IPT "1 3\n";
    
    `echo 1 3 > ipt`;
    `wc rom >> ipt`;
    `cat rom >> ipt`;
    `echo 1 4 >> ipt`;
    `wc nin >> ipt`;
    `cat nin >> ipt`;
    `echo 1 2 >> ipt`;
    `wc zpawinfo/vc_bare${edgename2} >> ipt`;
    `cat zpawinfo/vc_bare${edgename2} >> ipt`;

    `cp ipt ipt1`;
    `echo .false. >> ipt1`;
#
#    `echo 0.1 100 >> ipt1`;
    `echo "$final_dr $final_nr" >> ipt1`;
    `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXF${fullrad}/ropt`;
    move( "rpot", "${edgename}/zRXF${fullrad}/" );
    move( "rpothires", "${edgename}/zRXF${fullrad}/" );
#    `mv {rpot,rpothires} ${edgename}/zRXF${fullrad}/`;

    `cp ipt ipt1`;
    `echo .true. >> ipt1`;
    `wc zpawinfo/vpseud1${edgename2} >> ipt1`;
    `cat zpawinfo/vpseud1${edgename2} >> ipt1`;
    `wc zpawinfo/vvallel${edgename2} >> ipt1`;
    `cat zpawinfo/vvallel${edgename2} >> ipt1`;
#    `echo 0.1 100 >> ipt1`;
    `echo "$final_dr $final_nr" >> ipt1`;
    `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXT${fullrad}/ropt`;

    foreach( "rpot","rpothires","rom","nin")
      { move( $_ , "${edgename}/zRXT${fullrad}/"); }
    
#    `mv {rpot,rpothires,rom,nin} ${edgename}/zRXT${fullrad}/`;



#    `mv ximat ximat_full`;
#    `cp ximat_small ximat`;
     my $ximat_name = sprintf("ximat_small%04i",$itau);
    `ln -sf $ximat_name ximat`;

    open IPT, ">ipt" or die "Failed to open ipt\n$!";
    # Number of basis functions for chi
    print IPT "$grid_nb\n";
    close IPT;

    `$ENV{'OCEAN_BIN'}/xipps.x < ipt`;
    move( "ninduced", "nin" );
#    `mv ninduced nin`;
    open IPT, ">ipt" or die "Failed to open ipt\n$!";
    print IPT "$fullrad\n";
#    print IPT $ibase;
#    print IPT "$epsilon\n";
    close IPT;

#    `echo $fullrad > ipt`;
#    `echo $ibase >> ipt`;
#    `cat epsilon >> ipt`;
##    `cat ibase epsilon >> ipt`;
    `$ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
    move( "reopt", "rom" );
#    `mv reopt rom`;
    `echo 1 3 > ipt`;
    `wc rom >> ipt`;
    `cat rom >> ipt`;
    `echo 1 4 >> ipt`;
    `wc nin >> ipt`;
    `cat nin >> ipt`;
    `echo 1 2 >> ipt`;
    `wc zpawinfo/vc_bare${edgename2} >> ipt`;
    `cat zpawinfo/vc_bare${edgename2} >> ipt`;
    `cp ipt ipt1`;
    `echo .true. >> ipt1`;
    `wc zpawinfo/vpseud1${edgename2} >> ipt1`;
    `cat zpawinfo/vpseud1${edgename2} >> ipt1`;
    `wc zpawinfo/vvallel${edgename2} >> ipt1`;
    `cat zpawinfo/vvallel${edgename2} >> ipt1`;
#    `echo 0.1 100 >> ipt1`;
    `echo "$final_dr $final_nr" >> ipt1`;
    `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXS${fullrad}/ropt`;
#    `mv {rpot,rpothires,rom,nin} ${edgename}/zRXS${fullrad}/`;
    foreach( "rpot","rpothires","rom","nin")
      { move( $_ , "${edgename}/zRXS${fullrad}/"); }

  }
}
close HFINLIST;

`touch done`;

my $core_offset = `cat core_offset`;
chomp $core_offset;
if( $core_offset =~ m/false/i )
{
	print "No core shift\n";
	unlink "core_shift.txt" if( -e "core_shift.txt" );
} else
{
	`perl $ENV{'OCEAN_BIN'}/core_shift.pl > core_shift.log 2> core_shift.err`;
}

exit 0;
