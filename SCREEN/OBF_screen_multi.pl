#!/usr/bin/perl

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
my $trace_tolerance = "5.0d-15";

# Spline for Hamiltonian
my $ham_kpoints = "4 4 4";

my $band_start = 1;
my $band_stop  = 800;

my $screen_nkpt = "2 2 2";



# Step 1: Create support files
my @CommonFiles = ("znucl", "paw.hfkgrid", "paw.fill", "paw.opts", "pplist", "paw.shells", "ntype",
                   "natoms", "typat", "taulist", "nedges", "edges", "caution", "epsilon", "k0.ipt", 
                   "ibase", "scfac", "rscale", "rprim", "para_prefix", "paw.nbands", "core_offset",
                   "paw.nkpt", "pool_control");
my @ExtraFiles = ("specpnt", "Pquadrature" );
my @DFTFiles = ("rhoofr", "nscf.out", "system.rho.dat");

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

open IN, "epsilon" or die "Failed to open epsilon\n$!";
my $epsilon = <IN>;
chomp $epsilon;
close IN;

open IN, "ibase" or die "Failed to open ibase\n$!";
my $ibase;
while( <IN> ) { $ibase .= $_; }

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

if( -e "../DFT/ham_kpoints" )
{
	`cp ../DFT/ham_kpoints .`;
	$ham_kpoints = `cat ham_kpoints`;
	chomp($ham_kpoints);
}

$screen_nkpt = `cat paw.nkpt`;
chomp($screen_nkpt);


`ln -s ../DFT/Out .`;

my $fermi = 0;
open SCF, "nscf.out" or die "$!";
while( my $line = <SCF> )
{
  if( $line  =~  m/the Fermi energy is\s+([+-]?\d\S+)/ )
  {
    $fermi = $1;
    print "Fermi level found at $fermi eV\n";
    last;
  }
}
$fermi = $fermi/13.60569252;
`echo "$fermi" > efermiinrydberg.ipt`;

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

system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0
    or die "Failed to run bvecs.pl\n";

`tail -n 1 rhoofr > nfft` if( -e "rhoofr" );
system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0
  or die "Failed to run rhoofg.x\n";
`wc -l rhoG2 > rhoofg`;
`sort -n -k 6 rhoG2 >> rhoofg`;


print "Running PAW Setup\n";
system("$ENV{'OCEAN_BIN'}/pawsetup.x") == 0 or die "$!\nFailed to run pawsetup.x\n";

print "Running avg.x\n";
system("$ENV{'OCEAN_BIN'}/avg.x") == 0 or die "$!\nFailed to run avg.x\n";


`ln -sf ../PAW/zpawinfo zpawinfo`;
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
open MKRB_CONTROL, ">mkrb_control" or die "Failed to open mkrb_control\n$!\n";
print MKRB_CONTROL "8 25\n";
print MKRB_CONTROL "$nedges\n";


my $temp_edgename;
my $temp_rad = sprintf("%03.2f",$rads[0]);

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

open VCx, "zpawinfo/vcxxxxx${temp_edgename}R${temp_rad}" or die "Failed to open vcxxxxxx\n$!";
open VPERT, ">vpert" or die;
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
#`cp zpawinfo/vcxxxxx${temp_edgename}R${temp_rad} ./tmp`;
#`wc tmp > vpert`;
#`cat tmp >> vpert`;

print "ocean_builder.x\n";
system("$para_prefix $ENV{'OCEAN_BIN'}/ocean_builder_mult.x  $pool_size  < builder.in 1> builder.out 2> builder.err") == 0
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


  $edgename = sprintf("z%2s%02i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
  print "$edgename\n";
  unless( -d $edgename )
  {
    mkdir $edgename or die "Failed to make dir $edgename\n";
  }

  my $avden =  sprintf("avg%2s%02i",$elname,$elnum);
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
    
#    `mkdir -p ${edgename}/zRXT${fullrad}`;
#    `mkdir -p ${edgename}/zRXF${fullrad}`;
#    `mkdir -p ${edgename}/zRXS${fullrad}`;
#    chdir "$edgename";
#    `ln -s -f zRXT${fullrad} zR${fullrad}`;
#    chdir "../";

#    $screen_nkpt =~ m/(\d+)\s+(\d+)\s+(\d+)/;
#    my $np_builder = $1*$2*$3;
##      system("builder.x < builder.in") == 0 or die;
##      system("$ENV{'OCEAN_BIN'}/builder.x") == 0 or die;
#    system("$para_prefix $ENV{'OCEAN_BIN'}/ocean_builder.x  $pool_size < builder.in >& builder.out ") == 0
#        or die "$!\nFailed to run ocean_builder.x\n";


    my $ximat_name = sprintf("ximat%04i",$itau);
    unlink "ximat" if ( -e "ximat" );
    symlink $ximat_name, "ximat";
#    `ln -sf $ximat_name ximat`;
    
    open IPT, ">ipt" or die "Failed to open ipt\n$!";
    # Number of basis functions for chi
#    `echo 24 > ipt`;
    print IPT "24\n";
    close IPT;
    `$ENV{'OCEAN_BIN'}/xipps.x < ipt`;
    move( "ninduced", "nin" );
#    `mv ninduced nin`;

    open IPT, ">ipt" or die "Failed to open ipt\n$!";
    print IPT "$fullrad\n";
#    `echo $fullrad > ipt`;
#    `cat ibase epsilon >> ipt`;
    print IPT $ibase;
    print IPT "$epsilon\n";
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
    `wc zpawinfo/vcxxxxx${edgename2} >> ipt`;
    `cat zpawinfo/vcxxxxx${edgename2} >> ipt`;

    `cp ipt ipt1`;
    `echo .false. >> ipt1`;
    `echo 0.1 100 >> ipt1`;
    `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXF${fullrad}/ropt`;
    move( "rpot", "${edgename}/zRXF${fullrad}/" );
    move( "rpothires", "${edgename}/zRXF${fullrad}/" );
#    `mv {rpot,rpothires} ${edgename}/zRXF${fullrad}/`;

    `cp ipt ipt1`;
    `echo .true. >> ipt1`;
    `wc zpawinfo/vvpseud${edgename2} >> ipt1`;
    `cat zpawinfo/vvpseud${edgename2} >> ipt1`;
    `wc zpawinfo/vvallel${edgename2} >> ipt1`;
    `cat zpawinfo/vvallel${edgename2} >> ipt1`;
    `echo 0.1 100 >> ipt1`;
    `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXT${fullrad}/ropt`;

    foreach( "rpot","rpothires","rom","nin")
      { move( $_ , "${edgename}/zRXT${fullrad}/"); }
    
#    `mv {rpot,rpothires,rom,nin} ${edgename}/zRXT${fullrad}/`;



#    `mv ximat ximat_full`;
#    `cp ximat_small ximat`;
     my $ximat_name = sprintf("ximat_small%04i",$itau);
    `ln -sf $ximat_name ximat`;
    `echo 24 > ipt`;
    `$ENV{'OCEAN_BIN'}/xipps.x < ipt`;
    `mv ninduced nin`;
    `echo $fullrad > ipt`;
    `cat ibase epsilon >> ipt`;
    `$ENV{'OCEAN_BIN'}/vhommod.x < ipt`;
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
    `echo .true. >> ipt1`;
    `wc zpawinfo/vvpseud${edgename2} >> ipt1`;
    `cat zpawinfo/vvpseud${edgename2} >> ipt1`;
    `wc zpawinfo/vvallel${edgename2} >> ipt1`;
    `cat zpawinfo/vvallel${edgename2} >> ipt1`;
    `echo 0.1 100 >> ipt1`;
    `$ENV{'OCEAN_BIN'}/rscombine.x < ipt1 > ./${edgename}/zRXS${fullrad}/ropt`;
#    `mv {rpot,rpothires,rom,nin} ${edgename}/zRXS${fullrad}/`;
    foreach( "rpot","rpothires","rom","nin")
      { move( $_ , "${edgename}/zRXT${fullrad}/"); }

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
