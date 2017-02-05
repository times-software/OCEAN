#!/usr/bin/perl

# Copyright (C) 2015, 2016 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Copy;
use POSIX;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/rixs\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }

my %alphal = ( "0" => "s", "1" => "p", "2" => "d", "3" => "f" );

my @CommonFiles = ("epsilon", "xmesh.ipt", "nedges", "k0.ipt", "nbuse.ipt", 
  "cnbse.rad", "cnbse.ways", "metal", "cksshift", "cksstretch", "cksdq", 
  "cnbse.niter", "cnbse.spect_range", "cnbse.broaden", "cnbse.mode", "nphoton", "dft", 
  "para_prefix", "cnbse.strength", "serbse", "core_offset", "avecsinbohr.ipt", 
  "cnbse.solver", "cnbse.gmres.elist", "cnbse.gmres.erange", "cnbse.gmres.nloop", 
  "cnbse.gmres.gprc", "cnbse.gmres.ffff", "cnbse.write_rhs", "spin_orbit", "nspin", 
  "niter", "backf", "aldaf", "bwflg", "bande", "bflag", "lflag", "decut", "spect.h" );

my @DFTFiles = ("nelectron", "rhoofr");

my @DenDipFiles = ("kmesh.ipt", "masterwfile", "listwfile", "efermiinrydberg.ipt", 
                   "qinunitsofbvectors.ipt", "brange.ipt", "enkfile", "tmels", "nelectron", 
                   "eshift.ipt" );

my @WFNFiles = ("kmesh.ipt",  "efermiinrydberg.ipt", "qinunitsofbvectors.ipt", "brange.ipt", 
                "wvfcninfo", "wvfvainfo", "obf_control", "ibeg.h", "q.out", "tmels.info", 
                "val_energies.dat", "con_energies.dat" );

my @ExtraFiles = ("Pquadrature", "sphpts" );

my @PawFiles = ("hfinlist" , "xyz.wyck");

foreach (@CommonFiles) {
  copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}

open DFT, "dft" or die "Failed to open dft\n";
<DFT> =~ m/(\w+)/ or die;
my $dft_type = $1;
close DTF;
my $obf;
if( $dft_type =~ m/obf/i )
{
  $obf = 1;
}
else
{
  $obf = 0;
}


if( $obf == 1 ) 
{
  foreach (@DFTFiles) {
    copy( "../DFT/$_", $_) or die "Failed to get DFT/$_\n$!";
  }
  foreach (@WFNFiles) {
    copy( "../zWFN/$_", $_ ) or die "Failed to get zWFN/$_\n$!";
  }
  open OUT, ">tmel_selector" or die;
  print OUT "1\n";
  close OUT;
  open OUT, ">enk_selector" or die;
  print OUT "1\n";
  close OUT;
  open OUT, ">bloch_selector" or die;
  print OUT "1\n";
  close OUT;
}
else
{
  foreach (@DFTFiles) {
    copy( "../PREP/$_", $_) or die "Failed to get DFT/$_\n$!";
  }

  foreach (@DenDipFiles) {
    copy( "../PREP/BSE/$_", $_ ) or die "Failed to get PREP/BSE/$_\n$!" ;
  }
  open OUT, ">tmel_selector" or die;
  print OUT "0\n";
  close OUT;
  open OUT, ">enk_selector" or die;
  print OUT "0\n";
  close OUT;
  open OUT, ">bloch_selector" or die;
  print OUT "0\n";
  close OUT;

}

foreach (@ExtraFiles) {
  copy( "$ENV{'OCEAN_BIN'}/$_", $_ ) or die "Failed to get ../$_\n$!";
}

foreach (@PawFiles) {
  copy( "../SCREEN/$_", $_ ) or die "Failed to get ../SCREEN/$_\n$!";
}

##### Rixs requires that gmres be used #####
open IN, "cnbse.solver" or die "Failed to open cnbse.solver!\n$!";
my $line = <IN>;
close IN;
my $solver;
if( lc($line) =~ m/hay/ )
{
  print "Must select GMRES for use with rixs\n";
  print "Cannot continue\n";
  exit 1;
}
elsif( lc($line) =~ m/gmres/ )
{
  $solver = 'gmres';
}
else
{
  print "Trouble parsing cnbse.solver!!\n*** Will default to  Haydock recursion ***\n";
  print "Must select GMRES for use with rixs\n";
  print "Cannot continue\n";
  exit 1;
}
## Now if gmres we need to parse the inputs for that
my $gmres_footer = "";
my $gmres_count = 0;
if( $solver eq 'gmres' )
{
  open IN, "cnbse.gmres.nloop" or die "Failed to open cnbse.gmres.nloop\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  my $gmres_header = $line;

  open IN, "cnbse.broaden" or die "Failed to open cnbse.broaden\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  $line /= 27.2114;
  $gmres_header .= " " . $line;

  open IN, "cnbse.gmres.gprc" or die "Failed to open cnbse.gmres.gprc\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  $gmres_header .= " " . $line;

  open IN, "cnbse.gmres.ffff" or die "Failed to open cnbse.gmres.ffff\n$!";
  $line = <IN>;
  close IN;
  chomp $line;
  $gmres_header .= " " . $line;

  $gmres_header .= "  0.0\n";

  my $have_elist = 0;
  my $have_erange = 0;
  open IN, "cnbse.gmres.elist" or die "Failed to open cnbse.gmres.elist\n$!";
  $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    $gmres_footer = $gmres_header . "list\n";
    my $temp .= $line;
    my $i = 1;
    while( $line = <IN> )
    {
      $temp .= $line;
      $i++;
    }
    $gmres_count = $i;
    $gmres_footer .= "$i\n";
    $gmres_footer .= "$temp";
    $have_elist = 1;
    close IN;
  }

  open IN, "cnbse.gmres.erange" or die "Failed to open cnbse.gmres.erange\n$!";
  $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    if( $have_erange == 1 )
    {
      print "Both erange and elist were specified for GMRES. We are using erange\n";
    }
    else
    {
      $gmres_footer = $gmres_header . "loop\n";
      $line =~ m/([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/ or die "Failed to parse erange setting\n$line";
      my $estart = $1; my $estop = $3; my $estep = $5;
      $gmres_count = floor( ($estop - $estart + $estep*.9 ) / $estep );
      $gmres_footer .= $line;
      $have_erange = 1;
    }
    close IN;
  }

#  if( $have_erange + $have_elist == 2 )
#  {
#    print "Both erange and elist were specified for GMRES. We are using erange\n";
#  }
  if( $have_erange + $have_elist == 0 )
  {
    print "Neither elist nor erange were specified for GMRES!\n";
    print "*** WARNING ***\n";
    $have_elist = 1;
    $gmres_footer = $gmres_header . "list\n1\n0\n";
    $gmres_count = 1;
  }
}

##### Trigger serial bse fallback here
# later we should remove this and fold these two perl scripts together
my $run_serial = 0;
if( -e "serbse" )
{
  open IN, "serbse" or die "$!";
  <IN> =~ m/(\d)/;
  if( $1 == 1 )
  {
    print "Serial BSE requested!!\n";
    print "Parallel required\nCannot continue\n";
    exit 1;
  }
  close IN;
}
unless( -e "$ENV{'OCEAN_BIN'}/ocean.x" )
{
  print "Parallel BSE executable not present in $ENV{'OCEAN_BIN'}\n";
  print "Parallel required\nCannot continue\n";
  exit 1;
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
      copy( "../photon${i}", "photon${i}") or die "Failed to copy photon$i\n$!";
      push @photon_files, "photon${i}";
    }
    elsif( -e "../jtv${i}" )
    {
       copy( "../jtv${i}", "jtv{$i}") or die "Failed to copy jtv$1\n$!";
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
  print "Looking for available photon files:\n";
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
  print "!!!  Could not find any photon files  !!!\n  Will have to quit :(\n";
  exit 1;
}
else
{
  $nphoton = $#photon_files+1;
  if( $nphoton > 1 )
  {
    print "    Running with $nphoton photon files\n";
  }
  else
  {
    print "    Running with $nphoton photon file\n";
  }
  # Want to sort these back to numerical order (to avoid confusing/annoying users)
  my @sorted_photon_files = sort{ my ($anum,$bnum); $a =~ m/\w+(\d+)/; $anum = $1; $b =~ m/\w+(\d+)/; $bnum = $1; $anum <=> $bnum } @photon_files;
  @photon_files = @sorted_photon_files;

  foreach( @photon_files )
  {
    print "        $_\n";
    copy( "../$_", $_ ) or die "$!";
  }
}



##### misc other setup
#`echo gmanner > format65`;
copy( "kmesh.ipt", "kgrid" ) or die "$!";
copy( "k0.ipt", "scaledkzero.ipt" ) or die "$!";
copy( "qinunitsofbvectors.ipt", "cksdq" ) or die "$!";

##### qin must be non-zero
open QIN, "qinunitsofbvectors.ipt" or die "Failed to open qinunitsofbvectors.ipt\n$!";
<QIN> =~ m/(\d*\.?\d+)\s+(\d*\.?\d+)\s+(\d*\.?\d+)/ or die "Failed to parse qinunitsofbvectors.ipt\n";
close QIN;
if( abs( $1 ) + abs( $2 ) + abs( $3 ) < 0.0000001 ) 
{
  print "Non-zero q required for RIXS\n";
  exit 1;
}

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


# Set up mode -- need to run as xas which means ignore cnbse.mode
  my $is_xas = 1;
  open TMPFILE, "cnbse.niter" or die "Failed to open cnbse.niter\n$!";
  <TMPFILE> =~ m/(\d+)/ or die "Failed to parse cnbse.niter";
  my $num_haydock_iterations = $1;
  close TMPFILE;
  
  open TMPFILE, "cnbse.strength" or die "Failed to open cnbse.strength\n$!";
  <TMPFILE> =~ m/([0-9]*\.?[0-9]+)/ or die "Failed to parse cnbse.strength\n";
  my $interaction_strength = $1; 
  close TMPFILE;
    
# write cks.normal file
  open TMPFILE, ">cks.normal" or die "Failed to open cks.normal for writing\n$!";
  print TMPFILE ".true.\n";
  close TMPFILE;

#write mode file
  open TMPFILE, ">mode" or die "Failed to open mode for writing\n$!";
  print TMPFILE "$interaction_strength    $num_haydock_iterations\n";
  close TMPFILE;

###############
# If we are using QE/ABI w/o OBFs we need to set nbuse
my @brange;
my $run_text = '';
open NBUSE, "nbuse.ipt" or die "Failed to open nbuse.ipt\n";
<NBUSE> =~ m/(\d+)/ or die "Failed to parse nbuse.ipt\n";
my $nbuse = $1;
close NBUSE;
if( $obf == 1 )
{
  close RUNTYPE;
  if ($is_xas == 1 )
  {
    $run_text = 'XAS';
    if( $nbuse == 0 )
    {
      copy( "../zWFN/nbuse.ipt", "nbuse.ipt" ) or die "$!";
    }
    print "XAS!\n";
  } 
  else
  {
    if( $nbuse == 0 )
    {
      copy( "../zWFN/nbuse_xes.ipt", "nbuse.ipt" ) or die "$!";
    }
    $run_text = 'XES';
    print "XES!\n";
  }
}
else  ### Abi/QE w/o obf
{ 
#  my @brange;
  if ($nbuse == 0) {
    open BRANGE, "brange.ipt" or die "Failed to open brange.ipt\n";
    <BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Failed to parse brange.ipt\n";
    $brange[0] = $1;
    $brange[1] = $2;
    <BRANGE> =~ m/(\d+)\s+(\d+)/ or die "Failed to parse brange.ipt\n";
    $brange[2] = $1;
    $brange[3] = $2;
    close BRANGE;

    if( $is_xas == 1 )
    {
      $run_text = 'XAS';
      $nbuse = $brange[3] - $brange[1];
    }
    else
    {
      print "XES!\n";
      $run_text = 'XES';
      $nbuse = $brange[1] - $brange[0] + 1;
    }
    open NBUSE, ">nbuse.ipt" or die "Failed to open nbuse.ipt\n";
    print NBUSE "$nbuse\n";
    close NBUSE;
  }
  else
  {
    if( $is_xas == 1 )
    {
      $run_text = 'XAS';
    }
    else
    {
      print "XES!\n";
      $run_text = 'XES';
    }
  }
}


system("$ENV{'OCEAN_BIN'}/getnval.x") == 0 or die "Failed to get nval\n";


#####################
if( $obf == 1 )
{
  if( -e "../zWFN/u2par.dat" )
  {
    `ln -s ../zWFN/u2par.dat`;
    open OUT, ">bloch_type" or die;
    print OUT "new\n";
    close OUT;
  }
  else
  {
    `ln -s ../zWFN/u2.dat`;
  }
}
else  # We are using abi/qe path w/o obfs
{
  # grab .Psi
  `touch .Psi`;
  system("rm .Psi*");
  open LISTW, "listwfile" or die "Failed to open listwfile\n";
  while (<LISTW>) 
  {
    $_ =~ m/(\d+)\s+(\S+)/ or die "Failed to parse listwfile\n";
    system("ln -sf ../PREP/BSE/$2 .") == 0 or die "Failed to link $2\n";
  }  

  print "Running setup\n";
  system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0 or die "Setup failed\n";

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
}


my $core_offset;
open IN, "core_offset" or die "Failed to open core_offset\n$!";
$core_offset = <IN>;
chomp($core_offset);
close IN;


my $pawrad = `cat cnbse.rad`;
chomp($pawrad);
$pawrad = sprintf("%.2f", $pawrad);


my %unique_z;
my %unique_znl;

open RUNLIST, ">runlist";
my $hfinlength = `wc -l hfinlist`;
chomp($hfinlength);
$hfinlength *= ($#photon_files + 1 );
print "$hfinlength\n";
print RUNLIST "$hfinlength\n";


my $cls_average = 0;
my $cls_count = 0;
open EDGE, "hfinlist";
while (<EDGE>) {
  $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
  my $ppname = $1;
  my $znum = $2;
  my $nnum = $3;
  my $lnum = $4;
  my $elname = $5;
  my $elnum = $6;

#  for( my $i = 1; $i <= $#photon_files+1; $i++ ) 
  foreach my $way (@photon_files) 
  {
    $way =~ m/(\d+)$/ or die "Malformed photon file name:\t$way\n";
    my $i = $1;
    print RUNLIST "$znum  $nnum  $lnum  $elname  ${nnum}$alphal{$lnum}  $elnum  $i  $run_text\n";
  }



  # For each unique Z we need to grab some files from OPF
  unless( exists $unique_z{ "$znum" } )
  {
    my $zstring = sprintf("z%03i", $znum);
    print $zstring ."\n";
    `ln -sf ../OPF/zpawinfo/*${zstring}* .`;
    `ln -sf ../OPF/zpawinfo/phrc? .`;
    my $templine = `ls ../OPF/zpawinfo/*$zstring`;
    chomp($templine);
    my @pslist = split(/\s+/, $templine);
  }

  for( my $cks_iter = 1; $cks_iter <= 2; $cks_iter++ )
  {
    my $cks;
    if( $cks_iter == 2  ) {
      $cks = sprintf("cksc.${elname}%04u", $elnum );
    } 
    else {
      $cks = sprintf("cksv.${elname}%04u", $elnum );
    }
    print "CKS NAME = $cks\n";
    if( $obf == 1 )
    {
      copy( "../zWFN/$cks", $cks ) or die "Failed to grab $cks\n$!";
    }
    else # qe/abi w/o obf need to calculate cainkset
    {
      if( $cks_iter == 2  ) {
      open CKSNORM, ">cks.normal" or die "Failed to open cks.normal\n$!";
      print CKSNORM ".true.\n";
      close CKSNORM;
      $nbuse = $brange[3] - $brange[1];

      }
      else {
        open CKSNORM, ">cks.normal" or die "Failed to open cks.normal\n$!";
        print CKSNORM ".false.\n";
        close CKSNORM;
        $nbuse = $brange[1] - $brange[0] + 1;
      }

      open NBUSE, ">nbuse.ipt" or die "Failed to open nbuse.ipt\n";
      print NBUSE "$nbuse\n";
      close NBUSE;


      open ZNL, ">ZNL" or die;
      print ZNL "$znum  $nnum  $lnum\n";
      close ZNL;

      open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
      print CKSIN "1\n$elname  $elnum  cbinf\n";
      close CKSIN;


      print "cks\n";
      system("$ENV{'OCEAN_BIN'}/cks.x < cks.in > cks.log") == 0 or die;
      move( "cbinf0001", $cks ) or die "Failed to move cbinf0001 to $cks\n$!";
    }
  }

  my $zstring = sprintf("z%2s%04i_n%02il%02i", $elname, $elnum, $nnum, $lnum);
  copy( "../SCREEN/${zstring}/zR${pawrad}/rpot", "rpot.${zstring}" )
    or die "Failed to grab rpot\n../SCREEN/${zstring}/zR${pawrad}/rpot ./rpot.${zstring}\n";

  # If we don't want CLS then make sure the file is not here
  if( $core_offset =~ m/false/i )
  {
    if( -e "cls.${zstring}" ) 
    {
      unlink "cls.${zstring}" or die "Failed to remove cls.${zstring}\n$!";
    }
  }
  else
  {
    copy( "../SCREEN/${zstring}/zR${pawrad}/cls", "cls.${zstring}" ) 
      or warn "WARNING!\nCore-level shift support requested, but could not find ../SCREEN/${zstring}/zR${pawrad}/cls\n"
            . "No CLS will be done for this site!\n";
    $cls_count++;
    open IN, "cls.${zstring}" or die "Failed to open cls.${zstring}\n$!";
    my $cls = <IN>;
    chomp $cls;
    $cls_average += $cls;
  }


  $unique_z{ "$znum" } = 1;  
  $unique_znl{ "$znum $nnum $lnum" } = 1;
}
close EDGE;
close RUNLIST;
if( $cls_count > 0 )
{
  $cls_average /= $cls_count;
}
else
{
  $cls_average = 0;
}
open OUT, ">cls_average" or die "$!";
print OUT "$cls_average\n";
close OUT;


while ( my ($key, $value ) = each(%unique_znl) )
{ 
  $key =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die;
  my $znum = $1;
  my $nnum = $2;
  my $lnum = $3;

  open ZNL, ">ZNL" or die;
  print ZNL "$znum  $nnum  $lnum\n";
  close ZNL;

#
  foreach my $way (@photon_files) {
    copy( $way, "spectfile" ) or die "Failed to copy $way\n$!";
    system("$ENV{'OCEAN_BIN'}/meljtv.x");
    $way =~ m/(\d+)$/ or die "Misformed photon file name\n";
    my $i = $1;
    my $mel_targ = sprintf("mels.z%03un%02ul%02up%02u", $znum, $nnum, $lnum, $i );
    move( "mels", $mel_targ ) or die "Failed to move mels to $mel_targ\n$!";
  }  


}

my $ZNL = `cat ZNL`;
$ZNL =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die;
my $znum = $1;
my $nnum = $2;
my $lnum = $3;


open INFILE, ">bse.in" or die "Failed to open bse.in\n";
my $filename = sprintf("prjfilez%03u", $znum );

open TMPFILE, $filename or die "Failed to open $filename\n";
$line = <TMPFILE>;
close TMPFILE;
$line =~ m/(\d+)\s+(\d+)\s+\d+\s+\S+/ or die "Failed to match first line of $filename\n";

print INFILE "$lnum\t$1\t$2\n";

# spin orbit splitting
open IN, "spin_orbit" or die "Failed to open spin_orbit\n$!";
my $spin_orbit = <IN>;
chomp $spin_orbit;

if( $spin_orbit < 0 )
{
	my $lookup = sprintf("%1u%1s", $nnum, $alphal{$lnum}) or die;
	my $filename = sprintf("corezetaz%03u", $znum);
	print "$lookup\t$filename\n";
	$line = `grep $lookup $filename`;
}
else
{
	print "Overiding spin-orbit splitting! Using: $spin_orbit (eV)\n";
  $line = "$spin_orbit\n";
}
print INFILE $line;
$line =~ m/^\s*(\d+(\.\d+)?)/ or die;
$spin_orbit = $1;
open OUT, ">so.ipt" or die "$!";
print OUT "$spin_orbit\n";
close OUT;

close TMPFILE;
my $spectrange = `cat cnbse.spect_range`;
chomp($spectrange);
# check if spectrange needs correcting
unless( $spectrange =~ m/\d+\s+[-]?\d+(\.\d+)?\s+[-]?\d+(\.\d+)?/ )
{
  system("$ENV{'OCEAN_BIN'}/spect_range.pl") == 0 or die "Failed to run spect_range.pl";
  $spectrange = `cat cnbse.spect_range`;
  chomp($spectrange);
}
my $gamma0 = `cat cnbse.broaden`;
chomp($gamma0);

if( $solver eq 'gmres' )
{
  print INFILE "inv\n";
  print INFILE $gmres_footer . "\n";
}
else
{
  print INFILE "hay\n";
  print INFILE "$num_haydock_iterations  $spectrange  $gamma0  0.000\n";
}
close INFILE;

# open echamp.inp
open INFILE, ">echamp.inp" or die "$!";
print INFILE ".true.\n";
close INFILE;

# Run ocean.x once for all the selected energies to get the excitons
$ENV{"OMP_NUM_THREADS"}=1;

print "Running XAS through GMRES\n";
print "time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > core.log";
system("time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > core.log") == 0 or die "Failed to finish\n"; 


# Set up for running the valence pathway for RIXS

print "\nSetting up valence section\n";
move("runlist", "runlist.xas");

open IN, "runlist.xas" or die "Failed to open runlist.xas\n$!";
<IN>;
<IN> =~ m/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\w\w)\s+(\w\w)/ or die;
my $znum = $1;
my $nnum = $2;
my $lnum = $3;
my $elname = $4;
my $corelevel = $5;
close IN;

open RUNLIST, ">runlist" or die "Failed to open runlist for writing\n$!";
my $photon_combo = $nphoton * ($nphoton-1);
my $tot_gmres_count = $gmres_count * $photon_combo;
print RUNLIST "$tot_gmres_count\n";
for( my $e = 1; $e <= $gmres_count; $e++ )
{
  for( my $i = 1; $i <= $nphoton; $i++ )
  {
    for( my $j = 1; $j <= $nphoton; $j++ )
    {
      next if( $i == $j );
      print RUNLIST "$znum  $nnum  $lnum  $elname  $corelevel  0  $i  RXS  $e  $j\n";
    }
  }
}
close RUNLIST;

`tail -n 1 rhoofr > nfft`;

open INFILE, ">bse.in" or die "Failed to open bse.in\n";
print INFILE "0 0 0\n";
print INFILE "0 0 0\n";
my $spectrange = `cat spect.h`;
chomp($spectrange);
my $num_haydock_iterations = `cat niter`;
chomp($num_haydock_iterations);


print INFILE "hay\n";
print INFILE "$num_haydock_iterations  $spectrange  $gamma0  0.000\n";
close INFILE;

`cat bse.in`;
sleep 10;

print "Running valence for RIXS\n";
print "time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log";
system("time $para_prefix $ENV{'OCEAN_BIN'}/ocean.x > val.log") == 0 or die "Failed to finish\n";


exit 0;

