#!/usr/bin/perl

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/cnbse\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }

my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

my @CommonFiles = ("epsilon", "xmesh.ipt", "nedges", "k0.ipt", "nbuse.ipt", 
  "cnbse.rad", "cnbse.ways", "metal", "cksshift", "cksstretch", "cksdq", "cks.normal",
  "cnbse.niter", "cnbse.spect_range", "cnbse.broaden", "cnbse.mode", "nphoton");

my @AbinitFiles = ("avecsinbohr.ipt");

my @DenDipFiles = ("kmesh.ipt", "masterwfile", "listwfile", "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "tmels", "nelectron", "eshift.ipt" );

my @ExtraFiles = ("Pquadrature", "sphpts" );

my @PawFiles = ("hfinlist", "xyz.wyck");

foreach (@CommonFiles) {
  `cp ../Common/$_ .` == 0 or die "Failed to get Common/$_\n";
}
foreach (@AbinitFiles) {
  `cp ../DFT/$_ .` == 0 or die "Failed to get DFT/$_\n";
}
foreach (@DenDipFiles) {
  `cp ../PREP/BSE/$_ .` == 0 or die "Failed to get PREP/BSE/$_\n" ;
}
foreach (@ExtraFiles) {
  `cp $ENV{'OCEAN_BIN'}/$_ .` == 0 or die "Failed to get ../$_\n";
}

foreach (@PawFiles) {
  `cp ../SCREEN/$_ .` == 0 or die "Failed to get ../SCREEN/$_\n";
}


# Grab the needed photon files, copy them into the CNBSE directory,
# and store their names into the array @photon_files
my $nphoton = -1;
open NPHOTON, "nphoton" or die "Failed to open nphoton\n$!";
while (<NPHOTON>)
{
  if( $_ =~ m/(-?\d+)/ )
  {
    $nphoton = $1;
    last;
  }
}
close NPHOTON;

my @photon_files;
if( $nphoton > 0 ) 
{
  for( my $i = 1; $i <= $nphoton; $i++ )
  {
    if( -e "../photon${i}" )
    {
      `cp ../photon${i} .`;
      push @photon_files, "photon${i}";
    }
    elsif( -e "../jtv${i}" )
    {
      `cp ../jtv${i} .`;
      push @photon_files, "jtv${i}";
    }
    else
    {
      print "Could not find photon file # $i\n";
    }
  }
}
else
{
  print "Looking for available photon files\n";
  opendir DIR, "../" or die $!;
  while( my $file = readdir( DIR ) )
  {
    if( $file =~ m/^photon\d+$/ )
    {
      push @photon_files, $file;
    }
  }
  closedir DIR;
  
  if( $#photon_files == -1 )  # no photon files, fall back to jtv for now
  {
    opendir DIR, "../" or die $!;
    while( my $file = readdir( DIR ) )
    {
      if( $file =~ m/^jtv\d+$/ )
      {
        push @photon_files, $file;
      }
    }
    closedir DIR;
  }
}

if( $#photon_files == -1 )
{
  print "Could not find any photon files!\nWill have to quit :(\n";
  exit 1;
}
else
{
  $nphoton = $#photon_files+1;
  if( $nphoton > 1 ) 
  {
    print "\nRunning with $nphoton photon files\n";
  }
  else
  {
    print "\nRunning with $nphoton photon file\n";
  }
  foreach( @photon_files )
  {
    print "$_\n";
    `cp ../$_ .`;
  }
}


##### misc other setup
#`echo gmanner > format65`;
`cp kmesh.ipt kgrid`;
`cp k0.ipt scaledkzero.ipt`;
`mv cnbse.mode mode`;
`cp qinunitsofbvectors.ipt cksdq`;

system("rm .Psi*");
open LISTW, "listwfile" or die "Failed to open listwfile\n";
while (<LISTW>) {
  $_ =~ m/(\d+)\s+(\S+)/ or die "Failed to parse listwfile\n";
  system("ln -sf ../PREP/BSE/$2 .") == 0 or die "Failed to link $2\n";
}

system("$ENV{'OCEAN_BIN'}/getnval.x") == 0 or die "Failed to get nval\n";

open NBUSE, "nbuse.ipt" or die "Failed to open nbuse.ipt\n";
<NBUSE> =~ m/(\d+)/ or die "Failed to parse nbuse.ipt\n";
my $nbuse = $1;
close NBUSE;
my @brange;
if ($nbuse == 0) {
  open COND, "cks.normal" or die "Failed to open cks.normal\n";
  open BRANGE, "brange.ipt" or die "Failed to open brange.ipt\n";
  <BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Failed to parse brange.ipt\n";
  $brange[0] = $1;
  $brange[1] = $2;
  <BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Failed to parse brange.ipt\n";
  $brange[2] = $1;
  $brange[3] = $2;
  close BRANGE;

  if( <COND> =~ m/\.true\./ ) {
    print "cks in normal\n";
    $nbuse = $brange[3] - $brange[2] + 1;
  }
  else{
    print "val\n";
    $nbuse = $brange[1] - $brange[0] + 1;
  }
  open NBUSE, ">nbuse.ipt" or die "Failed to open nbuse.ipt\n";
  print NBUSE "$nbuse\n";
  close NBUSE;
}


#####################
print "Running setup\n";
system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";
`touch cbinf0001`;
`ln -sf cbinf0001 ufmi` == 0 or die;

if (-e "../PREP/BSE/u2.dat") 
{
  `ln -s ../PREP/BSE/u2.dat`;
}
else
{
  print "conugtoux\n";
  system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
  print "orthog\n";
  system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;
}

my $pawrad = `cat cnbse.rad`;
chomp($pawrad);
$pawrad = sprintf("%.2f", $pawrad);

open EDGE, "hfinlist";
while (<EDGE>) {
  $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
  my $ppname = $1;
  my $znum = $2;
  my $nnum = $3;
  my $lnum = $4;
  my $elname = $5;
  my $elnum = $6;

  my $zstring = sprintf("z%03i", $znum);
  print $zstring ."\n";
  `ln -sf ../PAW/zpawinfo/*${zstring}* .`;
#  `ln -sf ../PAW/zpawinfo/${zstring} .` == 0 or die;
  my $templine = `ls ../PAW/zpawinfo/*$zstring`;
  chomp($templine);
  my @pslist = split(/\s+/, $templine);
  foreach (@pslist) 
  {
    $_ =~ m/ae(\S+)/;
    `ln -sf ../PAW/zpawinfo/ae$1 .`;
    `ln -sf ae$1 ps$1`;
  }
  
  open ZNL, ">ZNL" or die;
#  `echo "$znum  $nnum  $lnum$" > ZNL`;
  print ZNL "$znum  $nnum  $lnum\n";
  close ZNL;
  $zstring = sprintf("z%03un%02ul%02u", $znum, $nnum, $lnum);
  $zstring = sprintf("z%2s%02i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
  system("cp ../SCREEN/${zstring}/zR${pawrad}/rpot ./rpotfull") == 0 
    or die "Failed to grab rpot\n../SCREEN/${zstring}/zR${pawrad}/rpot ./rpotfull\n";
#
#
  foreach my $way (@photon_files) {
    system("cp ${way} spectfile") ;#== 0 or die;
    system("$ENV{'OCEAN_BIN'}/meljtv.x");
    `mv mels mels.${way}`;
  }  

#  print "conugtoux\n";
#  system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
#  print "orthog\n";
#  system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;

  open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
  print CKSIN "1\n$elname  $elnum  cbinf\n";
  close CKSIN;

  print "cks\n";
  system("$ENV{'OCEAN_BIN'}/cks.x < cks.in > cks.log") == 0 or die;

  foreach my $file (@photon_files) {

    `cp mels.${file} mels`;

    $file =~ m/(\d+)/ or die;
    my $way = $1;

    print "dotter\n";
    system("echo 'cbinf0001' | $ENV{'OCEAN_BIN'}/dotter.x") == 0 or die;


    open INFILE, ">bse.in" or die "Failed to open bse.in\n";
    my $filename = sprintf("deflinz%03un%02ul%02u", $znum, $nnum, $lnum);

    open TMPFILE, $filename or die "Failed to open $filename\n";
    my $line = <TMPFILE>;
    close TMPFILE;
    
    print INFILE $line;
    my $lookup = sprintf("%1u%1s", $nnum, $alphal{$lnum}) or die;
    my $filename = sprintf("corezetaz%03u", $znum);
    print "$lookup\t$filename\n";
    my $line = `grep $lookup $filename`;
    print INFILE $line;
    
    print INFILE "hay\n";
    open TMPFILE, "cnbse.niter" or die "Failed to open niter\n";
    <TMPFILE> =~ m/(\d+)/ or die "Failed to parse niter\n";
    my $niter = $1;
    close TMPFILE;
    my $spectrange = `cat cnbse.spect_range`;
    chomp($spectrange);
    my $gamma0 = `cat cnbse.broaden`;
    chomp($gamma0);
    
    print INFILE "$niter  $spectrange  $gamma0  0.000\n";
    close INFILE;

    if( -e "../SCREEN/core_shift.txt" )
    {
      `head -n $elnum ../SCREEN/core_shift.txt | tail -n 1 > core_offset `;
    } else
    {
       `rm -f core_offset`;
    }

    system("$ENV{'OCEAN_BIN'}/cainmultip.x < bse.in > cm.log") == 0 or die "Failed to finish\n"; 

    my $store_string = sprintf("%2s.%04i_%2s_%02i", $elname, $elnum, $lookup, $way); 
    
    `mkdir -p ${store_string}/`;
    `cp absspct lanceigs mulfile nval.h ${store_string}/`;
    `cp absspct "absspct_${store_string}"`;
  }
}

exit 0;

