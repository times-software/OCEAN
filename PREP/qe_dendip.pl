#!/usr/bin/perl
# Copyright (C) 2014, 2016 - 2019 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Path qw( rmtree );
use Cwd 'abs_path';
use File::Compare;
use File::Copy;

###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/qe_dendip\.pl/;
#  $ENV{"OCEAN_BIN"} = $1;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################
my $RunDenDip = 0;

my $stat = 0;
$stat = 1 if (-e "done");
`rm -f done`;
my $oldden = 0;
$oldden = 1 if (-e "../DFT/old");


my @QEFiles     = ( "rhoofr", "efermiinrydberg.ipt" );
my @CommonFiles = ( "screen.nkpt", "nkpt", "qinunitsofbvectors.ipt", "avecsinbohr.ipt", "dft", 
                    "nspin", "xmesh.ipt", "dft.split", "prefix", "calc", "screen.wvfn", "screen.legacy", 
                    "screen.mode", "bse.wvfn", "k0.ipt" );
my @NewMethodFiles = ( "ntype", "typat", "natoms", "znucl", "taulist", "edges", "core_offset", "metal", "cksshift", 
                       "cksstretch", "nedges", "edges", "pplist", "opf.opts", "opf.fill" );
my @ExtraFiles = ("specpnt.5", "Pquadrature", "sphpts" );

foreach (@QEFiles) {
  system("cp ../DFT/$_ .") == 0 or die "Failed to copy $_\n";
}
foreach (@CommonFiles) {
  system("cp ../Common/$_ .") == 0 or die "Failed to copy $_\n";
}
foreach (@NewMethodFiles) {
  system("cp ../Common/$_ .") == 0 or die "Failed to copy $_\n";
}

#JTV
( copy "../OPF/hfinlist", "hfinlist" ) == 1 or die "Failed to copy hfinlist\n$!";
( copy "../OPF/xyz.wyck", "xyz.wyck" ) == 1 or die "Failed to copy hfinlist\n$!";


system("mv nkpt bse.nkpt") == 0 or die "Failed to rename nkpt $_\n";
print "$stat  $oldden\n";
unless ($stat && $oldden) {


  `tail -n 1 rhoofr > nfft`;

  system("$ENV{'OCEAN_BIN'}/nelectron.x") == 0 
    or die "Failed to run nelectron.x\n";

  system("$ENV{'OCEAN_BIN'}/bvecs.pl") == 0
    or die "Failed to run bvecs.pl\n";

  system("$ENV{'OCEAN_BIN'}/gvecs2.pl") == 0
    or die "Failed to run gvecs2.pl\n";
}


open IN, "calc" or die "Failed to open calc\n";
<IN> =~ m/(\w+)/ or die "Failed to parse calc\n";
my $calc = $1;
close IN;

open IN, "screen.mode" or die "Failed to open screen.mode\n";
<IN> =~ m/(\w+)/ or die "Failed to parse screen.mode\n";
my $screen_mode = $1;
close IN;
my $run_screen = 1;
if( $calc =~ m/val/i )
{
  $run_screen = 0 unless( $screen_mode =~ m/grid/i );
}

open IN, "screen.wvfn" or die "Failed to open screen.wvfn\n$!";
my $screenWvfn = <IN>;
close IN;

open IN, "screen.legacy" or die "Failed to open screen.legacy\n$!";
my $screenLegacy = <IN>;
chomp $screenLegacy;
close IN;

open IN, "bse.wvfn" or die "Failed to open bse.wvfn\n$!";
my $bseWvfn = <IN>;
close IN;

my $rundir;
my @nkpt;

## process screen wf files ##
if( $run_screen == 1 )
{
  open NKPT, "screen.nkpt" or die "Failed to open screen.nkpt";
  <NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse screen.nkpt\n";
  @nkpt = ($1, $2, $3);
  close NKPT;
  $rundir = "../DFT/SCREEN";

  unless( -e "PAW/done" && -e "${rundir}/old" ) 
  {
    `rm -r PAW` if (-e "PAW");
    mkdir "PAW"; 
    chdir "PAW";

    open NKPTS, ">nkpts" or die "Failed to open nkpts for writing\n";
    print NKPTS $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
    close NKPTS;


    foreach ("kmesh.ipt", "brange.ipt", "qinunitsofbvectors.ipt" ) {
      system("cp ../${rundir}/$_ .") == 0 or die "Failed to copy $_\n";
    }
    #`cp ../qinunitsofbvectors.ipt .`;
    `cp ../bvecs .`;
    `cp ../dft .`;
    `cp ../nspin .`;
    `cp ../${rundir}/umklapp .`;
    `cp ../prefix .`;

    my $prefix;
    open PREFIX, "prefix";
    $prefix = <PREFIX>;
    close (PREFIX);
    chomp( $prefix );


    #`cp -r ../${rundir}/Out .`;
    if( -l "Out" )  # Out is an existing link
    {
      unlink "Out" or die "Problem cleaning old 'Out' link\n$!";
    }
    elsif(  -d "Out" ) #or Out is existing directory
    {
      rmtree( "Out" );
    }
    elsif( -e "Out" ) #or Out is some other file
    {
      unlink "Out";
    }
    print "../$rundir/Out\n";
    symlink ("../$rundir/Out", "Out") == 1 or die "Failed to link Out\n$!";

    # New methods for skipping wfconvert
    if( $screenWvfn =~ m/qe(\d+)/ && $screenLegacy == 0 )
    {
      unlink "PAW/old";
      my $qeVersion = $1;
      print "Don't convert DFT. Using new method for screening wavefunctions: $qeVersion\n";
      `touch listwfile masterwfile enkfile`;
      if( $qeVersion == 62 )
      {
  #      print "../$rundir/enkfile\n";
        copy "../$rundir/enkfile", "enkfile";
      }
  #    else{ print $qeVersion . "\n"; }
    }
    else  # old method, run wfconvert
    {
      print "$screenWvfn $screenLegacy\n";
      print "$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml\n";
      system("$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml") == 0 
        or die "Failed to run qe_data_file.pl\n$!";

      system("$ENV{'OCEAN_BIN'}/wfconvert.x") == 0 
        or die "Failed to run wfconvert.x\n$!";
    }


    `touch done`;
    chdir "../";
  }
  else {
    `touch PAW/old`;
    print  "Nothing needed for SCREEN wfns\n";
  }


  print "Done with SCREEN files\n";
}
else
{
  print "Nothing needed for SCREEN wvfn\n";
}


## process bse wf files ##
my $split_dft = 0;
open IN, "dft.split" or die "$!\n";
if( <IN> =~ m/t/i )
{
  $split_dft = 1;
}
close IN;



open NKPT, "bse.nkpt" or die "Failed to open nkpt";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse nkpt\n";
@nkpt = ($1, $2, $3);
close NKPT;

$rundir = sprintf("../DFT/%03u%03u%03u", $nkpt[0], $nkpt[1], $nkpt[2]);

my @BSECommonFiles = ( "qinunitsofbvectors.ipt", "bvecs", "dft", "nelectron", "avecsinbohr.ipt", 
                       "nspin", "dft.split", "prefix", "natoms", "typat", "ntype","znucl", "taulist", 
                       "edges", "k0.ipt", "core_offset", "metal", "cksshift", "cksstretch", "bse.wvfn" );
my @rundirFiles = ( "kmesh.ipt", "brange.ipt", "umklapp" );
my @BSEBonusFiles = ("xmesh.ipt", "calc", "hfinlist", "xyz.wyck" );
#Checks for BSE prep
my $runBSE = 1;
if( -e "BSE/done" && -e "${rundir}/old" )
{
  $runBSE = 0;
  foreach (@BSECommonFiles)
  {
    if( compare( "$_", "BSE/$_") != 0 )
    {
      $runBSE = 1;
      print "Difference found in $_\n";
      last;
    }
  }
  if( $runBSE == 0 )
  {
    foreach( @rundirFiles )
    {
      if( compare( "${rundir}/$_", "BSE/$_") != 0 )
      {
        $runBSE = 1;
        print "Difference found in $_\n";
        last;
      }
    }
  }
  # Check xmesh changes only
  if( $runBSE == 0 )
  {
    $runBSE = 2 if( compare( "xmesh.ipt", "BSE/xmesh.ipt" ) != 0 );
    if( compare( "hfinlist", "BSE/hfinlist" ) != 0 || compare( "calc", "BSE/calc" ) != 0 )
    {
      $runBSE = 3;
    }
  }
}

#$runBSE = 0 if( $bseWvfn =~ m/qe62/ );

if( $runBSE != 0 )
{
  unlink "BSE/done";
  if( $runBSE == 1 )
  {
    rmtree( "BSE" ) if (-e "BSE");
    mkdir "BSE";
  }
  chdir "BSE";

  open NKPT, ">nkpts" or die "Failed to open nkpts for writing\n";
  print NKPT $nkpt[0]*$nkpt[1]*$nkpt[2] . "\n";
  close NKPT;

  foreach( @rundirFiles )
  {
    copy( "../${rundir}/$_", $_ );
  }
  foreach( @BSECommonFiles )
  {
    copy ( "../$_", $_ );
  }
  foreach( @BSEBonusFiles )
  {
    copy ( "../$_", $_ );
  }

  my $prefix;
  open PREFIX, "prefix" or die "$!\n";
  $prefix = <PREFIX>;
  chomp($prefix);
  close (PREFIX);



  if( -l "Out" )  # Out is an existing link
  {
    unlink "Out" or die "Problem cleaning old 'Out' link\n$!";
  } 
  elsif(  -d "Out" ) #or Out is existing directory
  {
    rmtree( "Out" );
  } 
  elsif( -e "Out" ) #or Out is some other file
  {
    unlink "Out";
  } 
  print "../$rundir/Out\n";
  symlink ("../$rundir/Out", "Out") == 1 or die "Failed to link Out\n$!";

  if( $bseWvfn =~ m/qe/ || $bseWvfn =~ m/new/ )
  {
    if( -e "Out/$prefix.save/data-file.xml" )
    {
      print "Detected QE54-style DFT run\n";
      open TMP, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
      print TMP "qe54\n";
      close TMP;
    }
    elsif( -e "Out/$prefix.save/data-file-schema.xml" )
    {
      print "Detected QE62-style DFT run\n";
      open TMP, ">", "wvfn.ipt" or die "Failed to open wvfn.ipt for writing\n$!";
      print TMP "qe62\n";
      close TMP;
    }
    else
    { 
      print "WARNING! Failed to detect style of QE output\nWill attempt to continue\n";
      copy "../bse.wvfn", "wvfn.ipt";
    }
    system("cp ../efermiinrydberg.ipt ./") == 0 
      or die "Failed to copy efermiinrydberg.ipt\n";
    unless( -e "zpawinfo" )
    {
      symlink ("../../OPF/zpawinfo", "zpawinfo" ) == 1 or die "Failed to link zpawinfo\n$!";
    }

    foreach( @ExtraFiles )
    {
      copy("$ENV{'OCEAN_BIN'}/$_", $_ ) or die;
    }

#    system( "/users/jtv1/cluster/Software/OCEAN/PREP/src/ocean_prep.x" );
    print "$ENV{'OCEAN_BIN'}/ocean_prep.x > ocean_prep.log 2>&1\n";
    system("$ENV{'OCEAN_BIN'}/ocean_prep.x > ocean_prep.log 2>&1" ) == 0
          or die "Failed to run ocean_prep.x\n$!";
  }
  else
  {

    foreach( @ExtraFiles )
    {
      copy("$ENV{'OCEAN_BIN'}/$_", $_ ) or die;
    }
    if( $runBSE == 1 )
    {
      if( $split_dft ) 
      {
        print "$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml Out/${prefix}_shift.save/data-file.xml\n";
        system("$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml Out/${prefix}_shift.save/data-file.xml") == 0
          or die "Failed to run qe_data_file.pl\n$!";
      }
      else
      {
        print "$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml\n";
        system("$ENV{'OCEAN_BIN'}/qe_data_file.pl Out/$prefix.save/data-file.xml") == 0
          or die "Failed to run qe_data_file.pl\n$!";
      }

      system("$ENV{'OCEAN_BIN'}/wfconvert.x system") == 0 
        or die "Failed to run wfconvert.x\n";

      system("$ENV{'OCEAN_BIN'}/ofermi.pl") == 0
        or die "Failed to run ofermi.pl\n";

      `cp eshift.ipt ../`;
      system("cp ../efermiinrydberg.ipt ./") == 0 
        or die "Failed to copy efermiinrydberg.ipt\n";
    }
    
    open TMPFILE, "calc" or die "Failed to open calc\n";
    my $mode = <TMPFILE>;
    close TMPFILE;
    chomp($mode);
    my $runCKS = 1;
    $runCKS = 0 if( lc($mode) =~ m/val/ );

    if( ( $runBSE == 1 || $runBSE == 3 ) && ( $runCKS == 1 ) )
    {
      copy( "kmesh.ipt", "kgrid" ) or die "$!";
      copy( "k0.ipt", "scaledkzero.ipt" ) or die "$!";
      copy( "qinunitsofbvectors.ipt", "cksdq" ) or die "$!";
      my $is_xas;
      open TMPFILE, "calc" or die "Failed to open calc\n";
      my $mode = <TMPFILE>;
      close TMPFILE;
      chomp($mode);
      if( lc($mode) =~ m/xes/ )
      {
        print "Calculating XES\n";
        $is_xas = 0;
      }
      elsif( lc($mode) =~ m/xas/ )
      {
        print "Calculating XAS\n";
        $is_xas = 1;
      }
      else
      {
        print "Unrecognized mode. Calculating XAS\n";
        $is_xas = 1;
      }

      open TMPFILE, ">cks.normal" or die "Failed to open cks.normal for writing\n$!";
      if( $is_xas == 1 )
      {
        print TMPFILE ".true.\n";
      }
      else
      {
        print TMPFILE ".false.\n";
      }
      close TMPFILE;

#      open CKS, ">cks.in" or die "Failed to open cks.in\n";
      my $znl_string = 0;
      my $ncks = 0;
      my $cks_string;

      my %unique_z;
      open EDGE, "hfinlist" or die "Failed to open hfinlist\n$!";
      while (<EDGE>) 
      {
        $_ =~ m/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/ or die;
        my $ppname = $1;
        my $znum = $2;
        my $nnum = $3;
        my $lnum = $4;
        my $elname = $5;
        my $elnum = $6;
      
        my $cks;
        if( $is_xas == 1  ) {
          $cks = "cksc.${elname}";
        }
        else {
          $cks = "cksv.${elname}"
        }

        # For each unique Z we need to grab some files from OPF
        unless( exists $unique_z{ "$znum" } )
        {
          my $zstring = sprintf("z%03in%02il%02i", $znum, $nnum, $lnum);
          print $zstring ."\n";

          `ln -sf ../../OPF/zpawinfo/gk*${zstring} .`;
          `ln -sf ../../OPF/zpawinfo/fk*${zstring} .`;
          `ln -sf ../../OPF/zpawinfo/melfile${zstring} .`;
          `ln -sf ../../OPF/zpawinfo/coreorb${zstring} .`;


          my $z3 = sprintf("z%03i", $znum);
          `ln -sf ../../OPF/zpawinfo/phrc?${z3} .`;
          `ln -sf ../../OPF/zpawinfo/prjfile${z3} .`;
          `ln -sf ../../OPF/zpawinfo/ft?${z3} .`;
          `ln -sf ../../OPF/zpawinfo/ae?${z3} .`;
          `ln -sf ../../OPF/zpawinfo/ps?${z3} .`;
          `ln -sf ../../OPF/zpawinfo/corezeta${z3} .`;
          `ln -sf ../../OPF/zpawinfo/radfile${z3} .`;

        }

        print "CKS NAME = $cks\n";
        my $temp_znl = sprintf "%i  %i  %i", $znum, $nnum, $lnum;
        if( $znl_string == 0 )
        { 
          $znl_string = $temp_znl;
          open ZNL, ">ZNL" or die;
          print ZNL "$znl_string\n";
          close ZNL;
        }
        
        if( $znl_string eq $temp_znl )
        { 
          $ncks++;
          $cks_string .= "$elname  $elnum  $cks\n";
        }
        else
        { 
          print "New ZNL!\nRunning $ncks through cks\n";
          unless ( $ncks == 0 )
          { 
            open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
            print CKSIN "$ncks\n$cks_string";
            close CKSIN;
            print "cks\n";
            system("$ENV{'OCEAN_BIN'}/cks.x < cks.in > cks.log") == 0 or die;
          }
          $znl_string = $temp_znl;
          open ZNL, ">ZNL" or die;
          print ZNL "$znl_string\n";
          close ZNL;
          $ncks = 1;
          $cks_string = "$elname  $elnum  $cks\n";
        }

      }
      close EDGE;
      
      unless ( $ncks == 0 )
      {
        print "Final cks: $ncks\n";
        open CKSIN, ">cks.in" or die "Failed to open cks.in\n";
        print CKSIN "$ncks\n$cks_string";
        close CKSIN;
        print "cks\n";
        system("$ENV{'OCEAN_BIN'}/cks.x < cks.in > cks.log") == 0 or die;
      }

    }

    # 3 is cks-only option, 1 & 2 will make new u2
    unless( $runBSE == 3 )
    {
      print "Running setup\n";
      system("$ENV{'OCEAN_BIN'}/setup2.x > setup.log") == 0
        or die "Failed to run setup\n";

      print "conugtoux\n";
      system("$ENV{'OCEAN_BIN'}/conugtoux.x > conugtoux.log");# == 0 or die;
      print "orthog\n";
      system("$ENV{'OCEAN_BIN'}/orthog.x > orthog.log") == 0 or die;
    }

  }

  `touch done`;
  chdir "../";
}


else {
  `touch BSE/old`;
  print "Nothing needed for bse wfns\n";
}

unless ($stat && $oldden) {
system("$ENV{'OCEAN_BIN'}/rhoofg.x") == 0
  or die "Failed to run rhoofg.x\n";
`wc -l rhoG2 > rhoofg`;
`sort -n -k 6 rhoG2 >> rhoofg`;
}

`touch done`;
exit 0;

