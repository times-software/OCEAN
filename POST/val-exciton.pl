#!/usr/bin/perl
#
## Copyright (C) 2015 OCEAN collaboration
##
## This file is part of the OCEAN project and distributed under the terms 
## of the University of Illinois/NCSA Open Source License. See the file 
## `License' in the root directory of the present distribution.
##
##
#
use strict;
use File::Copy;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/val-exciton\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }


my @CommonFiles = ("calc", "natoms", "ntype", "typat", "taulist", "znucl" );

foreach (@CommonFiles) {
   copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}


### run valence BSE with full diagonalization
## if RIXS was run we need to take things from there
open CALC, "calc" or die "Failed to open file calc: $!\n";
<CALC> =~ m/(\w+)/ or die "Failed to read file calc\n";
my $calc = lc($1);
close CALC;

if( $calc =~ m/val/i )
{
   open SOLV, "../NBSE/cnbse.solver" or die "Failed to open file cnbse.solver: $!\n";
   <SOLV> =~ m/(\w+)/ or die "Failed to read file cnbse.solver\n";
   my $solver = lc($1);
   close SOLV;
   # if gmres, copy everything from NBSE, otherwise, rerun as below
   if( $solver =~ m/gmres/i )
   {
      print "Found sovler is GMRES\n";
      chdir "../NBSE";
      map { copy($_, "../EXCITON/$_") or die "Failed to copy echamp" } glob("ehamp_*");
      chdir "../EXCITON";
      my @BseFiles = ("avecsinbohr.ipt", "bloch_selector", "brange.ipt", "kmesh.ipt", "k0.ipt", "nspin", "nbuse.ipt", "opcons", "qinunitsofbvectors.ipt", "runlist", "u2.dat", "xmesh.ipt", "ZNL" );
      foreach (@BseFiles) {
         copy( "../NBSE/$_", $_ ) or die "Failed to get NBSE/$_\n$!";
      }
   }
   else
   {
      `echo .true. > echamp.inp`;
      `echo gmres > ../Common/cnbse.solver`;
      system("$ENV{'OCEAN_BIN'}/nbse.pl") == 0 or die "Failed to run valence BSE\n$!";
   }
}
elsif( $calc =~ m/rixs/i )
{
    print "Still working on RIXS final state exciton plotting\n";
    # let's try rerunning full RIXS calculation, but using gmres at each step
    system("$ENV{'OCEAN_BIN'}/rixs_exc.pl > rixs.log") == 0 or die "Failed to run RIXS for final exciton states\n$!";
    print "Made it through RIXS rerun\n";

    # now we test whether the rest could work
    copy("rxsspct_Si.2p_01.00001.02","opcons");

#    # try the brute force approach
#    opendir my $rixsdir, "../RIXS" or die "Cannot open RIXS directory$!\n";
#    my @files = readdir $rixsdir;
#    closedir $rixsdir;
#    foreach my $file (@files) {
#       copy "../RIXS/$file", $file;
#    }
#    `echo .true. > echamp.inp`;
#    `echo gmres > ../Common/cnbse.solver`;
#    system("$ENV{'OCEAN_BIN'}/rixs_nbse_exc.pl") == 0 or die "Failed to run valence BSE for RIXS exciton states\n$!";
}


# make xyz.wyck
system("$ENV{'OCEAN_BIN'}/makewyck.x") == 0 or die "Failed to make xyz.wyck\n";


# make ehcoor.ipt
# does this matter for valence excitons?
`echo 0.1 0.1 0.1 > ehcoor.ipt`;


### want to loop over all ehamp files
# get the number of excitons we generated
`wc -l < opcons > nexc`;
open NEX, "nexc" or die "Failed to open nexc\n$!";
my $line = <NEX>;
close NEX;
chomp $line;
my $nexc = $line;


my $zi;
# Select the hole part of the exciton
`echo 0 > ehflag.ipt`;

for( my $i = 1; $i <= $nexc; $i++ )
{

   $zi = sprintf("%04d",$i);
   # make exciton_plot.ipt
   open EXCPLOT, ">exciton_plot.ipt";
   print EXCPLOT "ehamp_$zi\n";
   print EXCPLOT "hole_$zi.cube\n";
   print EXCPLOT "3 3 3\n";
   print EXCPLOT "-1 -1 -1\n";
   close EXCPLOT;

   # run cube generator
   system("$ENV{'OCEAN_BIN'}/val_exciton_plot.x") == 0 or die "Failed to run exciton plotter\n";

   # rename input file to store
   `mv exciton_plot.ipt exc_$zi.ipt`;
}


# Select the electron part of the exciton
`echo 3 > ehflag.ipt`;

for( my $i = 1; $i <= $nexc; $i++ )
{

   $zi = sprintf("%04d",$i);
   # make exciton_plot.ipt
   open EXCPLOT, ">exciton_plot.ipt";
   print EXCPLOT "ehamp_$zi\n";
   print EXCPLOT "electron_$zi.cube\n";
   print EXCPLOT "3 3 3\n";
   print EXCPLOT "-1 -1 -1\n";
   close EXCPLOT;

   # run cube generator
   system("$ENV{'OCEAN_BIN'}/val_exciton_plot.x") == 0 or die "Failed to run exciton plotter\n";

   # rename input file to store
   `mv exciton_plot.ipt exc_$zi.ipt`;
}



chdir "../";

exit 0;
