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
  $0 =~ m/(.*)\/optical-exc\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }


my @CommonFiles = ("calc", "natoms", "ntype", "typat", "taulist", "znucl", "vnbse.gmres.erange", "vnbse.gmres.elist" );

foreach (@CommonFiles) {
   copy( "../Common/$_", $_ ) or die "Failed to get Common/$_\n$!";
}


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
   #   map { copy($_, "../EXCITON/$_") or die "Failed to copy absspct" } glob("abss*");
      chdir "../EXCITON";
      my @BseFiles = ("avecsinbohr.ipt", "bloch_selector", "brange.ipt", "kmesh.ipt", "k0.ipt", "nspin", "nbuse.ipt", "opcons", "qinunitsofbvectors.ipt", "runlist", "u2.dat", "xmesh.ipt", "ZNL" );
      foreach (@BseFiles) {
         copy( "../NBSE/$_", $_ ) or die "Failed to get CNBSE/$_\n$!";
      }
   }
   else
   {
      `echo .true. > echamp.inp`;
      `echo gmres > ../Common/cnbse.solver`;
      system("$ENV{'OCEAN_BIN'}/nbse.pl") == 0 or die "Failed to run valence BSE\n$!";
   }


# make xyz.wyck
system("$ENV{'OCEAN_BIN'}/makewyck.x") == 0 or die "Failed to make xyz.wyck\n";


# make ehcoor.ipt
# does this matter for valence excitons?
`echo 0.1 0.1 0.1 > ehcoor.ipt`;


   my $zi;
   my $nve = get_nve();


# Select the hole part of the exciton
print "Working on the hole part\n";
`echo 0 > ehflag.ipt`;

#if( $calc =~ m/val/i )
#{
#   my $zi;
#   my $nve = get_nve();
   for( my $i = 1; $i <= $nve; $i++ )
   {
      $zi = sprintf("%04d",$i);
      # make exciton_plot.ipt
      open EXCPLOT, ">exciton_plot.ipt";
      print EXCPLOT "ehamp_$zi \n";
      print EXCPLOT "hole_$zi.cube\n";
      print EXCPLOT "3 3 3\n";
      print EXCPLOT "-1 -1 -1\n";
      close EXCPLOT;

      # run cube generator
      system("$ENV{'OCEAN_BIN'}/val_exciton_plot.x") == 0 or die "Failed to run exciton plotter\n";

      # rename input file to store
      `mv exciton_plot.ipt exc_$zi.ipt`;
   }
#}



# Select the electron part of the exciton
print "Working on the electron part\n";
`echo 3 > ehflag.ipt`;

#if( $calc =~ m/val/i )
#{
#   my $zi;
#   my $nve = get_nve();
    for( my $i = 1; $i <= $nve; $i++ )
    {
        $zi = sprintf("%04d",$i);
        # make exciton_plot.ipt
        open EXCPLOT, ">exciton_plot.ipt";
        print EXCPLOT "ehamp_$zi \n";
        print EXCPLOT "electron_$zi.cube\n";
        print EXCPLOT "3 3 3\n";
        print EXCPLOT "-1 -1 -1\n";
        close EXCPLOT;

        # run cube generator
        system("$ENV{'OCEAN_BIN'}/val_exciton_plot.x") == 0 or die "Failed to run exciton plotter\n";

        # rename input file to store
        `mv exciton_plot.ipt exc_$zi.ipt`;

    }
#}


`mkdir -p Valence/`;
opendir my $excdir, "./" or die "Cannot open directory$!\n";
my @files = readdir $excdir;
closedir $excdir;
foreach my $file (@files) {
    move $file, "Valence/$file";
}
# may want to introduce a check here (if Valence/Core exists)
move "./Valence/Core/", "./Core";

chdir "../";

exit 0;


sub get_nce {

  my $nce = 0;
  my $have_elist = 0;
  my $have_erange = 0;

  open IN, "vnbse.gmres.elist" or die "Failed to open vnbse.gmres.elist\n$!";
  my $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    $have_elist = 1;
    $nce = 1;
    while( $line = <IN> )
    {
      $nce++;
    }
    close IN;
  }

  open IN, "vnbse.gmres.erange" or die "Failed to open vnbse.gmres.erange\n$!";
  $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    $have_erange = 1;
    my @erange = split ' ', $line;
    my $emin = @erange[0];
    my $emax = @erange[1];
    my $estep = @erange[2];

    for( my $c = $emin; $c <= $emax; $c += $estep ) {
            $nce++;
    }
#    print "We have $nve energies\n";
    close IN;
  }

  if( $have_erange + $have_elist == 2 )
  {
    print "Both erange and elist were specified for GMRES. We are using erange\n";
  }
  elsif( $have_erange + $have_elist == 0 )
  {
    print "Neither elist nor erange were specified for GMRES!\nFalling back to Haydock\n";
  }

  return $nce;
}


sub get_nve {

  my $nve = 0;
  my $have_elist = 0;
  my $have_erange = 0;

  open IN, "vnbse.gmres.elist" or die "Failed to open vnbse.gmres.elist\n$!";
  my $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    $have_elist = 1;
    $nve = 1;
    while( $line = <IN> )
    {
      $nve++;
    }
    close IN;
  }

  open IN, "vnbse.gmres.erange" or die "Failed to open vnbse.gmres.erange\n$!";
  $line = <IN>;
  if( $line =~ m/false/ )
  {
    close IN;
  }
  else
  {
    $have_erange = 1;
    my @erange = split ' ', $line;
    my $emin = @erange[0];
    my $emax = @erange[1];
    my $estep = @erange[2];

    for( my $c = $emin; $c <= $emax; $c += $estep ) {
            $nve++;
    }
#    print "We have $nve energies\n";
    close IN;
  }

  if( $have_erange + $have_elist == 2 )
  {
    print "Both erange and elist were specified for GMRES. We are using erange\n";
  }
  elsif( $have_erange + $have_elist == 0 )
  {
    print "Neither elist nor erange were specified for GMRES!\nFalling back to Haydock\n";
  }

  return $nve;
}
