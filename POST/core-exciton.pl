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
  $0 =~ m/(.*)\/core-exciton\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }


my @CommonFiles = ("calc", "natoms", "ntype", "typat", "taulist", "znucl" );

foreach (@CommonFiles) {
   copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}


### run core BSE with full diagonalization
## if RIXS was run we can just take everything from there
open CALC, "calc" or die "Failed to open file calc: $!\n";
<CALC> =~ m/(\w+)/ or die "Failed to read file calc\n";
my $calc = lc($1);
close CALC;
if( $calc =~ m/rixs/i )
{
   chdir "../RIXS";
   map { copy($_, "../EXCITON/$_") or die "Failed to copy echamp" } glob("echamp_*"); 
   map { copy($_, "../EXCITON/$_") or die "Failed to copy absspct" } glob("abss*");
   chdir "../EXCITON";
   my @RixsFiles = ("avecsinbohr.ipt", "bloch_selector", "brange.ipt", "kmesh.ipt", "k0.ipt", "nspin", "nbuse.ipt", "qinunitsofbvectors.ipt", "runlist.xas", "u2.dat", "xmesh.ipt", "xyz.wyck", "ZNL" );
   foreach (@RixsFiles) {
      copy( "../RIXS/$_", $_ ) or die "Failed to get RIXS/$_\n$!";
   }
   move("runlist.xas","runlist");
}
elsif( $calc =~ m/xas/i || $calc =~ m/xes/i )
{
   print "Found calc is XAS or XES\n";
# check cnbse.solver in CNBSE
   open SOLV, "../CNBSE/cnbse.solver" or die "Failed to open file cnbse.solver: $!\n";
   <SOLV> =~ m/(\w+)/ or die "Failed to read file cnbse.solver\n";
   my $solver = lc($1);
   close SOLV;
# if gmres, copy everything from CNBSE, otherwise, rerun as below
   if( $solver =~ m/gmres/i )
   {
      print "Found sovler is GMRES\n";
      chdir "../CNBSE";
      map { copy($_, "../EXCITON/$_") or die "Failed to copy echamp" } glob("echamp_*");
      if ( $calc =~ m/xas/i )
      {
          map { copy($_, "../EXCITON/$_") or die "Failed to copy absspct" } glob("abss*");
      }
      else
      {
          map { copy($_, "../EXCITON/$_") or die "Failed to copy xesspct" } glob("xess*");
      }
      chdir "../EXCITON";
      my @BseFiles = ("avecsinbohr.ipt", "bloch_selector", "brange.ipt", "kmesh.ipt", "k0.ipt", "nspin", "nbuse.ipt", "qinunitsofbvectors.ipt", "runlist", "u2.dat", "xmesh.ipt", "xyz.wyck", "ZNL" );
      foreach (@BseFiles) {
         copy( "../CNBSE/$_", $_ ) or die "Failed to get CNBSE/$_\n$!";
      }
    }
    else
    {
        `echo .true. > echamp.inp`;
        `echo gmres > ../Common/cnbse.solver`;
        system("$ENV{'OCEAN_BIN'}/cnbse_mpi.pl") == 0 or die "Failed to run core BSE\n$!";
    }
}

### want to loop over all echamp files
### need to collect some information on absspct list
#
# parse runlist
my @elm;
my @site;
my @core;
my @ph;
open RL, "runlist" or die "Failed to open runlist\n$!";
my $nspec = <RL>;
chomp $nspec;
# print "Number of absspct files : $nspec\n";
# need to build arrays of length nspec for elm name, core state, site num, photon num
my $line;
my @entry;
my $ic=1;
while ( $ic <= $nspec ) {
    $line = <RL>;
    chomp $line;
    @entry = split /\s+/, $line;
    push @elm, $entry[3];
    push @core, $entry[4];
    push @site, $entry[5];
    push @ph, $entry[6];
    $ic += 1;
}
close RL;

# declare variables for labelling each exciton
my $label;
my $zsite;
my $zph;
my $zj;


# get the number of exciton energies
$zsite = sprintf("%04d_",$site[0]);
$zph = sprintf("%02d",$ph[0]);
if ( $calc =~ m/xas/i || $calc =~ m/rixs/i )
{
    `wc -l < absspct_$elm[0].$zsite$core[0]_$zph > nexc`;
}
else
{
    `wc -l < xesspct_$elm[0].$zsite$core[0]_$zph > nexc`;
}

open NEX, "nexc" or die "Failed to open nexc\n$!";
my $line = <NEX>;
close NEX;
chomp $line;
my $nexc = $line;

### end of collecting information on absspct list


# Select the hole part of the exciton
print "Starting on the hole part of the exciton\n";
`echo 0 > ehflag.ipt`;


# declare variables for obtaining exciton origin
my @tau;


## presently hole part is not supported for core excitons
#  runs electron part instead
my $irun=0;
if ( $irun )
{
for( my $j = 1; $j <= $nexc; $j++ )
{
   for( my $i = 1; $i <= $nspec; $i++ )
   {

      # creation $label for each exciton file
      $zsite = sprintf("%04d_",$site[$i-1]);
      $zph = sprintf("%02d",$ph[$i-1]);
      $zj = sprintf("%04d",$j);
      $label = "$elm[$i-1].$zsite$core[$i-1]_$zph.$zj";

      # get the location of the absorbing atom, write to ehcoor.ipt
      @tau = get_tau($elm[$i-1],$i);
      unlink "ehcoor.ipt";
      open EHCOOR, ">ehcoor.ipt";
      print EHCOOR join("  ", @tau) . "\n";
      close EHCOOR;

      # create the exciton plotting input file
      open EXCPLOT, ">exciton_plot.ipt";
      print EXCPLOT "echamp_$label\n";
      print EXCPLOT "hole_$label.cube\n";
      print EXCPLOT "3 3 3\n";      # could be made an input in the future
      print EXCPLOT "-1 -1 -1\n";   # could be made an input in the future
      close EXCPLOT;

      # run cube generator
#      print "Starting exciton_plot\n";
      system("$ENV{'OCEAN_BIN'}/exciton_plot.x") == 0 or die "Failed to run exciton plotter\n";
#      print "Back from exciton_plot\n";

      # rename input file
#      `mv exciton_plot.ipt exc_$label.ipt`;

   }
}
}


# Select the electron part of the exciton
print "Now for the electron part of the exciton\n";
`echo 3 > ehflag.ipt`;

for( my $j = 1; $j <= $nexc; $j++ )
{
   for( my $i = 1; $i <= $nspec; $i++ )
   {

      # creation $label for each exciton file
      $zsite = sprintf("%04d_",$site[$i-1]);
      $zph = sprintf("%02d",$ph[$i-1]);
      $zj = sprintf("%04d",$j);
      $label = "$elm[$i-1].$zsite$core[$i-1]_$zph.$zj";

      # get the location of the absorbing atom, write to ehcoor.ipt
      @tau = get_tau($elm[$i-1],$i);
      unlink "ehcoor.ipt";
      open EHCOOR, ">ehcoor.ipt";
      print EHCOOR join("  ", @tau) . "\n";
      close EHCOOR;

      # create the exciton plotting file
      open EXCPLOT, ">exciton_plot.ipt";
      print EXCPLOT "echamp_$label\n";
      print EXCPLOT "electron_$label.cube\n";
      print EXCPLOT "3 3 3\n";      # could be made an input in the future
      print EXCPLOT "-1 -1 -1\n";   # could be made an input in the future
      close EXCPLOT;

      # run cube generator
      system("$ENV{'OCEAN_BIN'}/exciton_plot.x") == 0 or die "Failed to run exciton plotter\n";

      # rename input file
      `mv exciton_plot.ipt exc_$label.ipt`;

   }
}

## move everything into a new folder Core/ in case we also have a valence exciton (i.e. for RIXS)
## might want to do this only in the case of RIXS
`mkdir -p Core/`;
opendir my $excdir, "./" or die "Cannot open directory$!\n";
my @files = readdir $excdir;
closedir $excdir;
foreach my $file (@files) {
    move $file, "Core/$file";
}


chdir "../";

exit 0;


sub get_tau {
   my ($elmx,$ix) = @_;
   my @tau;

   my @taux;
   my @tauy;
   my @tauz;

   # parse xyz.wyck for elements and coordinates
   open XYZ, "xyz.wyck" or die "Failed to open xyz.wyck$!\n";
   my $ntau = <XYZ>;
   chomp $ntau;
   my $line;
   my $i = 1;
   while ( $i <= $ntau ) {
      $line = <XYZ>;
      chomp $line;
      @entry = split /\s+/, $line;
      push @elm, $entry[1];
      push @taux, $entry[2];
      push @tauy, $entry[3];
      push @tauz, $entry[4];
      $i += 1;
   }
   close XYZ;

#   print "Looking for element $elmx, instance $ix\n";

   # extract the coordinates of the ix^th entry for element elmx
   my $ic = 0;
   $i = 1;
   while ( $i <= $ntau ) {
      if ( $elm[$i-1] eq $elmx ) {
         $ic += 1;
         if ( $ic == $ix){
             push @tau, $taux[$i-1];
             push @tau, $tauy[$i-1];
             push @tau, $tauz[$i-1];
         } 
      }
      $i += 1;
   }

   return @tau;

}

