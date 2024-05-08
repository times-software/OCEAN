#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use POSIX;
use File::Spec;
use File::Copy;
use Math::Trig;

#TODO: should loop here over split if we have it
sub QErunNSCF
{
  my ( $hashRef, $specificHashRef, $shift ) = @_;

#  my $dirname = "k" . $specificHashRef->{'kmesh'}[0] . '_' . $specificHashRef->{'kmesh'}[1] . '_' 
#              . $specificHashRef->{'kmesh'}[2] . "q" . $specificHashRef->{'kshift'}[0] . '_' 
#              . $specificHashRef->{'kshift'}[1] . '_'  . $specificHashRef->{'kshift'}[2];

  my $dirname;
  my $dualKpt = 0;
  if( $specificHashRef->{'nonzero_q'} && not $specificHashRef->{'split'}) {
    $dirname = sprintf "ks%i_%i_%iq%.6f_%.6f_%.6f", $specificHashRef->{'kmesh'}[0],
                    $specificHashRef->{'kmesh'}[1], $specificHashRef->{'kmesh'}[2], 
                    $specificHashRef->{'kshift'}[0], $specificHashRef->{'kshift'}[1],
                    $specificHashRef->{'kshift'}[2];
    $dualKpt = 1;
  } else {
    $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $specificHashRef->{'kmesh'}[0], 
                    $specificHashRef->{'kmesh'}[1], $specificHashRef->{'kmesh'}[2], 
                    $specificHashRef->{'kshift'}[0], $specificHashRef->{'kshift'}[1],
                    $specificHashRef->{'kshift'}[2];
  }

  print "$dirname\n";

  mkdir $dirname unless( -d $dirname );
  chdir $dirname;

  QEsetupNSCF( $shift );


  # open input file
  my $fileName = "nscf.in";
  $fileName = "nscf_shift.in" if( $shift );

  open my $fh, ">", $fileName or die "Failed to open $fileName.\n$!";

  my $calcFlag = 0;

  my $kptString;
#  if( $specificHashRef->{'isGamma'} ) 
  if( $specificHashRef->{'kmesh'}[0] == 1 && $specificHashRef->{'kmesh'}[1] == 1 
   && $specificHashRef->{'kmesh'}[2] == 1 && $specificHashRef->{'kshift'}[0] == 0
   && $specificHashRef->{'kshift'}[1] == 0 && $specificHashRef->{'kshift'}[2] == 0 )  
  {
    $kptString = "K_POINTS gamma\n";
  }
  else 
  {
    my $nk = $specificHashRef->{'kmesh'}[0] * $specificHashRef->{'kmesh'}[1] 
           * $specificHashRef->{'kmesh'}[2];
    $nk *=2  if( $dualKpt );
    $kptString = "K_POINTS crystal\n  $nk\n" . kptGen( $specificHashRef, $dualKpt );  
  }

  # Modify QEprintInput to take two different hashes and a string 
  QEprintInput( $fh, $hashRef, $specificHashRef, $calcFlag, $kptString ) ;
  
  close $fh;

  my $ncpus = 1;
  my $npool = 1;
  my $nbd = 1;
  my $nk = $specificHashRef->{'kmesh'}[0] * $specificHashRef->{'kmesh'}[1] * $specificHashRef->{'kmesh'}[2];
  $ncpus = $hashRef->{'computer'}->{'ncpus'};
  $ncpus = $nk if( $nk < $ncpus );
  $npool = $ncpus;

  ($ncpus, $npool, $nbd) = QErunTest( $hashRef, "nscf.in", $ncpus, $npool, $nbd );

  QEsetupNSCF( $shift );

  my $prevNpool = $npool;
  ($ncpus, $npool, $nbd) = QErunTest( $hashRef, "nscf.in", $ncpus, $npool, $nbd);

  QEsetupNSCF( $shift );

  if( $prevNpool * 2 < $npool ) {
    ($ncpus, $npool, $nbd) = QErunTest( $hashRef, "nscf.in", $ncpus, $npool, $nbd);

    QEsetupNSCF( $shift );
  }



  # full run

  print "$ncpus  $npool\n";
  $specificHashRef->{'ncpus'} = $ncpus*1;
  $specificHashRef->{'npool'} = $npool*1;

  my $prefix = $hashRef->{'computer'}->{'para_prefix'};
  $prefix =~ s/$hashRef->{'computer'}->{'ncpus'}/$ncpus/ if( $hashRef->{'computer'}->{'ncpus'} != $ncpus );

  my $cmdLine = " -npool $npool";
  if( $nbd > 1 ) { 
    $cmdLine .= " -pd true -ntg $nbd";
  }

  QErunPW( $hashRef->{'general'}->{'redirect'}, $prefix, $cmdLine, "nscf.in", "nscf.out" );


  my $errorCode = QEparseOut( "nscf.out", $specificHashRef );
  return $errorCode if( $errorCode );

  $errorCode = QEparseEnergies( $hashRef, $specificHashRef );

  chdir updir();
  return $errorCode;
}

sub QErunParseEnergies
{ 
  my ( $hashRef, $specificHashRef, $shift ) = @_;

  my $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $specificHashRef->{'kmesh'}[0],
                    $specificHashRef->{'kmesh'}[1], $specificHashRef->{'kmesh'}[2], 
                    $specificHashRef->{'kshift'}[0], $specificHashRef->{'kshift'}[1],
                    $specificHashRef->{'kshift'}[2];
  
  chdir $dirname;

  my $errorCode = QEparseEnergies( $hashRef, $specificHashRef );
  
  chdir updir();
  return $errorCode;
}

sub QEsetupNSCF
{
  my ($shift ) = @_;

  my $outDir = "Out";
  $outDir = "Out_shift" if ( $shift );
    
  mkdir $outDir unless( -d $outDir );
  $outDir = catdir( $outDir, "system.save" );
  mkdir $outDir unless( -d $outDir );
  
  my $targDir = catdir( updir(), "Out", "system.save" );
  print "$targDir   $outDir\n";
  print `pwd`;
  my @fileList = ( "charge-density.dat", "spin-polarization.dat", "occup.txt", 
                   "charge-density.kin.dat", "ekin-density.dat", 
                   "data-file.xml", "data-file-schema.xml" );

  foreach my $f (@fileList) {
    my $file = catfile( $targDir, $f );
    if( -e $file ) {
      copy $file, catfile( $outDir, $f );
    }
  }
}

sub QErunDensity
{
  my $hashRef = $_[0];

  # Make SCF k-point string
  my $kptString;
  if( $hashRef->{'scf'}->{'kmesh'}[0] * $hashRef->{'scf'}->{'kmesh'}[1] * $hashRef->{'scf'}->{'kmesh'}[2] == 1 
      && $hashRef->{'scf'}->{'kshift'}[0] == 0 && $hashRef->{'scf'}->{'kshift'}[1] == 0 
      && $hashRef->{'scf'}->{'kshift'}[2] == 0 ) {
   $kptString = "K_POINTS gamma\n";
  } else { 
    $kptString = "K_POINTS automatic\n" 
               . sprintf "%i %i %i %i %i %i\n", $hashRef->{'scf'}->{'kmesh'}[0], $hashRef->{'scf'}->{'kmesh'}[1],
                                                $hashRef->{'scf'}->{'kmesh'}[2], $hashRef->{'scf'}->{'kshift'}[0],
                                                $hashRef->{'scf'}->{'kshift'}[1], $hashRef->{'scf'}->{'kshift'}[2];
  }
  

  # open input file
  open my $fh, ">scf.in" or die "Failed to open scf.in.\n$!";

  my $calcFlag = 1;

  # Modify QEprintInput to take two different hashes and a string 
  QEprintInput( $fh, $hashRef, $hashRef->{'scf'}, $calcFlag, $kptString ) ;

  close $fh;

#  # test run for k-points
#  print "Testing parallel QE execution\n";
#  open TMP, '>', "system.EXIT" or die "Failed to open file system.EXIT\n$!";
#  close TMP;
#
#  QErunPW( $hashRef->{'general'}->{'redirect'}, $hashRef->{'computer'}->{'ser_prefix'}, "" , "scf.in", "test.out" );
#  
#  # parse test results
#  my $npool = 1;
#  my $ncpus = $hashRef->{'computer'}->{'ncpus'};
#  if( open TMP, "test.out" )
#  {
#    my $actualKpts = -1;
#    my $numKS = -1;
#    my $mem = -1;
#    while (<TMP>)
#    {
#      if( $_ =~ m/number of Kohn-Sham states=\s+(\d+)/ )
#      {
#        $numKS = $1;
#      }
#      elsif( $_ =~ m/number of k points=\s+(\d+)/ )
#      {
#        $actualKpts = $1;
#      }
#      elsif( $_ =~ m/Estimated max dynamical RAM per process\s+>\s+(\d+\.\d+)/ )
#      {
#        $mem = $1;
#      }
#      last if( $actualKpts > 0 && $numKS > 0 && $mem > 0 );
#    }
#    close TMP;
#
#    if( $actualKpts == -1 )
#    {
#      print "Had trouble parsing test.out\nDidn't find number of k points\n";
#    }
#    else
#    {
#      ($ncpus, $npool) = QEPoolControl( $actualKpts, $numKS, $mem, $hashRef->{'computer'} );
#    }
#  }
#  else
#  {
#    print "Had trouble parsing test.out\n. Will attempt to continue.\n";
#  }

  my $ncpus = 1;
  my $npool = 1;
  my $nbd = 1;
  my $nk = $hashRef->{'kmesh'}[0] * $hashRef->{'kmesh'}[1] * $hashRef->{'kmesh'}[2];
  $ncpus = $hashRef->{'computer'}->{'ncpus'};
  $ncpus = $nk if( $nk < $ncpus );
  $npool = $ncpus;

  ($ncpus, $npool, $nbd) = QErunTest( $hashRef, "scf.in", 1, 1, 1 );

  my $prevNpool = $npool;
  ($ncpus, $npool, $nbd) = QErunTest( $hashRef, "scf.in", $ncpus, $npool, $nbd);
  
  if( $prevNpool * 2 < $npool ) {
    ($ncpus, $npool, $nbd) = QErunTest( $hashRef, "scf.in", $ncpus, $npool, $nbd);
  }

  # full run

  print "$ncpus  $npool\n";
  $hashRef->{'scf'}->{'ncpus'} = $ncpus*1;
  $hashRef->{'scf'}->{'npool'} = $npool*1;

  my $prefix = $hashRef->{'computer'}->{'para_prefix'};
  $prefix =~ s/$hashRef->{'computer'}->{'ncpus'}/$ncpus/ if( $hashRef->{'computer'}->{'ncpus'} != $ncpus );

  my $cmdLine = " -npool $npool ";
  if( $nbd > 1 ) {
    $cmdLine .= " -pd true -ntg $nbd";
  }

  QErunPW( $hashRef->{'general'}->{'redirect'}, $prefix, $cmdLine, "scf.in", "scf.out" );


  my $errorCode = QEparseOut( "scf.out", $hashRef->{'scf'} );

  return $errorCode;
}

sub QErunTest
{
  my ( $hashRef, $fileName, $previousNcpus, $previousNpool, $previousNbd ) = @_;
    # test run for k-points
  print "Testing parallel QE execution\n";
  open TMP, '>', "system.EXIT" or die "Failed to open file system.EXIT\n$!";
  close TMP;

  if( $previousNcpus == 1 ) {
    QErunPW( $hashRef->{'general'}->{'redirect'}, $hashRef->{'computer'}->{'ser_prefix'}, "" , $fileName, "test.out" );
  } else {
    my $prefix = $hashRef->{'computer'}->{'para_prefix'};
    $prefix =~ s/$hashRef->{'computer'}->{'ncpus'}/$previousNcpus/ 
        if( $hashRef->{'computer'}->{'ncpus'} != $previousNcpus );
  
    my $cmdLine = " -npool $previousNpool ";
    if( $previousNbd > 1 ) {
      $cmdLine .= " -pd true -ntg $previousNbd";
    }
    QErunPW( $hashRef->{'general'}->{'redirect'}, $prefix, $cmdLine, $fileName, "test.out" );
  }

  # parse test results
  my $npool = 1;
  my $ncpus = $hashRef->{'computer'}->{'ncpus'};
  my $nbd = 1;
  if( open TMP, "test.out" )
  { 
    my $actualKpts = -1;
    my $numKS = -1;
    my $mem = -1;
    my $nspin = 1; 
    my $FFT = -1;
    my $OMP = -1;
    $nspin = $hashRef->{'general'}->{'nspin'} if( exists $hashRef->{'general'}->{'nspin'} );
    while (<TMP>)
    { 
      if( $_ =~ m/number of Kohn-Sham states=\s+(\d+)/ )
      { 
        $numKS = $1;
      }
      elsif( $_ =~ m/number of k points=\s+(\d+)/ )
      { 
        $actualKpts = $1 * $nspin;
      }
      elsif( $_ =~ m/Estimated max dynamical RAM per process\s+>\s+(\d+\.\d+)\s(\w+)/ )
      { 
        $mem = $1;
        $mem *= 1000 if (lc($2) eq 'gb' );
      }
      elsif( $_ =~ m/FFT dimensions:\s\(\s*(\d+),\s*(\d+),\s*(\d+)/ ) {
        $FFT = $3;
      }
      elsif( $_ =~ m/Threads\/MPI process:\s+(\d+)/ ) {
        $OMP = $1;
      }
      last if( $actualKpts > 0 && $numKS > 0 && $mem > 0 && $FFT > 0 && $OMP > 0 );
    }
    close TMP;
    
    if( $actualKpts < 0 )
    { 
      print "Had trouble parsing test.out\nDidn't find number of k points\n";
    }
    else
    { 
      ($ncpus, $npool, $nbd) = QEPoolControl( $actualKpts, $numKS, $mem, $FFT, $OMP, $hashRef->{'computer'} );
    }
  }
  else
  {
    print "Had trouble parsing test.out\n. Will attempt to continue.\n";
  }

  return ( $ncpus, $npool, $nbd);
}

# parse various info from the 
sub QEparseOut
{
  my( $outFile, $hashRef ) = @_;

  my $errorCode = 1;
  

  my $fermi;
  my $highest;
  my $lowest;
  my $xmlFile = catfile( "Out", "system.save", "data-file-schema.xml" );
  if( -e $xmlFile ) {
    open SCF, "<", $xmlFile or die "Failed to open $xmlFile\n";

    print "QE62!!\n";
    $hashRef->{'version'} = '6.2';

    open OUT, ">", "wvfn.ipt" or die $!;
    print OUT "qe62\n";
    close OUT;

    while( my $scf_line = <SCF> )
    {
      if( $scf_line =~ m/\<highestOccupiedLevel\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $highest = $1;
        $hashRef->{'highest'} = $1*1;
      }
      elsif( $scf_line =~ m/\<lowestUnoccupiedLevel\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $lowest = $1;
        $hashRef->{'lowest'} = $1*1;
      }
      elsif( $scf_line =~ m/\<fermi_energy\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $fermi = $1*1;
      }
      # We just average the two for spin=2 
      elsif( $scf_line =~ m/\<two_fermi_energies\>([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)\s+([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
      {
        $fermi = ($1+$3)/2;
      }
      elsif( $scf_line =~ m/\<two_fermi_energies\>\s*$/ ) {
        $scf_line = <SCF>;
        if( $scf_line =~ m/([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)\s+([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ )
        {
          $fermi = ($1+$3)/2;
        }
      }
#      elsif( $scf_line =~ m/\<two_fermi_energies\>\s*
      elsif( $scf_line =~ m/NAME=\"PWSCF\" VERSION=\"([\d\.\w]+)\"/ )
      {
        $hashRef->{'versionString'} = $1;
        $hashRef->{'versionString'} =~ m/(^\d+\.?\d*)/;
        $hashRef->{'version'} = $1;
      }
      elsif( $scf_line =~ m/wall>(\d+\.\d+[eE]?\+?\d+)/ )
      {
        $hashRef->{'time'} = $1*1;
      }
      elsif( $scf_line =~ m/etot>(-?\d+\.\d+[eE]?\+?\d?)/ )
      {
        $hashRef->{'etot'} = $1*1;
      }
      elsif( $scf_line =~ m/convergence_achieved>(\w+)/ )
      {
        $errorCode = 0 if( $1 =~ m/t/i );
      }
      elsif( $scf_line =~ m/<calculation>nscf/ )
      {
        $errorCode = 0;
      }
      elsif( $scf_line =~ m/<nelec>(\d\.\d+([eE]\+?\d+)?)/ )
      {
        $hashRef->{'nelec'} = $1*1;
      }
      elsif( $scf_line =~ m/<n_scf_steps>(\d+)</ ) 
      {
        $errorCode = 0 if( $1 > 0 );
      }
    }
    close SCF;
    if( defined $highest && defined $lowest ) {
      $fermi = ($highest + $lowest)/2; # Move from Ha to Ryd
      $hashRef->{'highest'} = $highest;
      $hashRef->{'lowest'} = $lowest;
      $hashRef->{'fermi'} = $fermi;
      printf "%g  %g  %g\n", $highest*13.60569253*2, $lowest*13.60569253*2, $fermi*13.60569253*2;
      printf "%g  %g  %g\n", $highest, $lowest, $fermi;
    } else {
      $fermi *= 1; # Move from Ha to Ryd
      $hashRef->{'fermi'} = $fermi;
    }

    print "Convergence not achieved\n" if( $errorCode == 1 );
    return $errorCode; 
  } elsif( -e catfile( "Out", "system.save", "data-file.xml" ) ) {

    print "QE54!!\n";
    $hashRef->{'version'} = '5.4';
    
    open OUT, ">", "wvfn.ipt" or die $!;
    print OUT "qe54\n";
    close OUT;
    my $units;

    open SCF, "<", catfile( "Out", "system.save", "data-file.xml" ) or die "Failed to open $xmlFile\n";

    while( my $scf_line = <SCF> )
    {
      if( $scf_line =~ m/\<UNITS_FOR_ENERGIES UNITS=\"(\w+)/ )
      {
        $units = $1;
      }
      elsif( $scf_line =~ m/\<FERMI_ENERGY/ )
      {
        $scf_line = <SCF>;
        $scf_line =~ m/([-+]?\d+\.\d+[Ee]?[-+]?(\d+)?)/ or die "$scf_line";
        $fermi = $1;
      }
      elsif(  $scf_line =~ m/NAME=\"PWSCF\" VERSION=\"([\d\.\w]+)\"/ )
      {
        $hashRef->{'versionString'} = $1;
        $hashRef->{'versionString'} =~ m/(^\d+\.?\d*)/;
#        $hashRef->{'versionString'} =~ m/([\d\.]+)/;
        $hashRef->{'version'} = $1;
#        print ">>>>>>> $hashRef->{'version'} <<<<<<<<< \n";
      }

    }
    close SCF;
    if( $hashRef->{'version'} < 1 ) {
      die "Failed to match version : $hashRef->{'version'}\n";
    }
    if( defined( $fermi ) ) {
      if( $units =~ m/hartree/i )
      {
        $fermi *= 1;
      }
      elsif( $units =~ m/eV/i )
      {
        $fermi /= 13.60569253/2;
      }

      $hashRef->{'fermi'} = $fermi;
      $errorCode = 0;
    }

    open OUTFILE, "<", $outFile or die "Failed to open $outFile\n";
    while( my $line=<OUTFILE> )
    {
      if( $line =~ m/number of electrons\s+=\s+(\d+\.\d+)/ )
      {
        $hashRef->{'nelec'} = $1*1;
      }
    }
    close OUTFILE;

  } else {
  
    open OUTFILE, "<", $outFile or die "Failed to open $outFile\n";
    while( my $line=<OUTFILE> )
    {
      if( $line =~ m/!\s+total energy\s+=\s+(-?\d+\.\d+)/ )
      {
        $hashRef->{'etot'} = $1*1;
#        $totEnergy = $1;
  #      print ">>>>> $totEnergy\n";
      }
      elsif( $line =~ m/convergence has been achieved/ )
      {
        $errorCode = 0;
      }
      elsif( $line =~ m/End of band/ )
      {
        $errorCode = 0;
      }
      elsif( $line  =~  m/the Fermi energy is\s+([+-]?\d+\.?\d+)/ )
      {
        $fermi = $1;
        print "Fermi level found at $fermi eV\n";
        $fermi = $fermi/13.60569253/2;
      }
      elsif( $line  =~  m/Fermi energies are\s+([+-]?\d+\.?\d+)\s+([+-]?\d+\.?\d+)/ )
      {
        $fermi = ($1+$2)/2;
        print "Fermi level found at $fermi eV\n";
        $fermi = $fermi/13.60569253/2;
      }
      elsif( $line =~ m/highest occupied, lowest unoccupied level \(ev\):\s+([+-]?\d+\.?\d+)\s+([+-]?\d+\.?\d+)/ )
      {
        $fermi = ($1+$2)/2;
        $hashRef->{'highest'} = $1/13.60569253/2;
        $hashRef->{'lowest'} = $2/13.60569253/2;
        print "Fermi level found at $fermi eV\n";
        $fermi = $fermi/13.60569253/2;
      }
      elsif( $line =~ m/number of electrons\s+=\s+(\d+\.\d+)/ )
      {
        $hashRef->{'nelec'} = $1*1;
      }
      elsif(  $line =~ m/PWSCF\s+v\.([\d\.\w]+)/ )
      {
        $hashRef->{'versionString'} = $1;
#        $hashRef->{'versionString'} =~ m/([\d\.]+)/;
        $hashRef->{'versionString'} =~ m/(^\d+\.?\d*)/;
        $hashRef->{'version'} = $1;
#        print ">>>>>>> $hashRef->{'version'} <<<<<<<<< \n";
      } 

    }
    $hashRef->{'fermi'} = $fermi if( defined( $fermi ) );
    close OUTFILE;
  }
  return 2 unless( defined( $fermi ) );
  print "FERMI: $hashRef->{'fermi'}\n";

  if( exists( $hashRef->{'lowest'} ) ){
    if( $hashRef->{'lowest'} <= $hashRef->{'highest'} ) {
      print "ERROR!!\n  Insulator selected, but there is no gap.\n"
                    ."  This run will not continue\n";
      $errorCode = 10;
      return $errorCode;
    }
  }
  return $errorCode;
}

sub QErunDFPT
{
  my ($hashRef) = @_;

  if( $hashRef->{'structure'}->{'metal'} ) {
    $hashRef->{'structure'}->{'epsilon'} = $hashRef->{'epsilon'}->{'metal_max'};
    return 0;
  }

  my $nnode = 1;
  my $npool = 1;
  my $ncpus = 1;
  if( open IN, "<", "scf.out" ) {
    while(<IN>) {
      $ncpus = $1 if( $_ =~ m/Number of MPI processes:\s+(\d+)/ );
      $nnode = $1 if( $_ =~ m/(\d+)\s+nodes/ );
      if( $_ =~ m/npool\s+=\s+(\d+)/ )
      {
        $npool = $1;
        last;
      }
    }
    close IN;
  }

  open PH, ">", "ph.in" or die "$!";
  print PH "title\n&inputph\n"
          . "  prefix = \'system\'\n"
          . "  outdir = \'Out\'\n"
          .  "  epsil = .true.\n"
          .  "  start_irr = 1\n"
          .  "  last_irr = 0\n"
          .  "  trans = .false\n"
          .  "  reduce_io = .true.\n"
          .  "/\n0 0 0\n";
  close PH;

  my $n = $nnode;
  $n = $npool if( $npool > $nnode );
  $ncpus = $n if( $n > $ncpus );
  while( $ncpus > 1 && ( $ncpus % $n != 0 ) ) {
    $ncpus --;
  }
  my $prefix = $hashRef->{'computer'}->{'para_prefix'};
  $prefix =~ s/$hashRef->{'computer'}->{'ncpus'}/$ncpus/ if( $hashRef->{'computer'}->{'ncpus'} != $ncpus );
  
  print  "$prefix $ENV{'OCEAN_ESPRESSO_PH'} -npool $n  -inp ph.in > ph.out 2>&1\n";
  system("$prefix $ENV{'OCEAN_ESPRESSO_PH'} -npool $n  -inp ph.in > ph.out 2>&1\n") == 0
    or die "Failed to run ph.x\n";
  open IN, "ph.out" or die;

  my @epsilon;
  $epsilon[0] = -1;
  while (<IN>)
  {
    if( $_ =~ m/Dielectric constant in cartesian axis/ )
    {
      <IN>;
      <IN> =~ m/(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/;
      $epsilon[0] = $1;
      <IN> =~ m/(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/;
      $epsilon[1] = $2;
      <IN> =~ m/(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)/;
      $epsilon[2] = $3;
      last;
    }
  }
  close IN;
  $hashRef->{'structure'}->{'epsilon'} = ( $epsilon[0] + $epsilon[1] + $epsilon[2] ) / 3;
  $hashRef->{'epsilon'}->{'epsilon'} = $hashRef->{'structure'}->{'epsilon'};

  return 9 if( $epsilon[0] < 0 );

  printf "DFPT done, epsilon = %f\n", $hashRef->{'structure'}->{'epsilon'};
  return 0;
}


sub QEparseDensityPotential
{
  my ($hashRef, $type, $emin, $emax ) = @_;

  printf "DENPOT VERSION %s\n", $hashRef->{'scf'}->{'version'};

  my $flag;
  my $filplot;
  my $infile;
  my $convert;
  my $bonusInputs ="";
  if( $type eq 'density' ) {
    $flag = 0;
    $filplot = 'system.rho';
    $infile = 'pp.in';
    $convert = "system.rho rhoofr";
  } elsif( $type eq 'sc-density' ) {
    # emin and emax are set to be equal to each other in case of errors
    $flag = 0;
    $filplot = 'system.val.rho';
    $infile = 'pp3.in';
    $convert = "system.val.rho val.rhoofr";
    if( abs( $emin - $emax ) > 0.1 && $hashRef->{'scf'}->{'version'} > 6.95) {
      $flag = 23;
      $bonusInputs = sprintf "  emin = %.15f\n  emax = %.15f\n", $emin, $emax;
    }
  } elsif( $type eq 'potential' ) {
    $flag = 11;
    $filplot = 'system.pot';
    $infile = 'pp2.in';
    $convert = "system.pot potofr";
  } else {
    return 1;
  }

  my $npool = 1;
  my $ncpus = $hashRef->{'computer'}->{'ncpus'};
  $ncpus = $hashRef->{'scf'}->{'ncpus'} if( exists $hashRef->{'scf'}->{'ncpus'} && $hashRef->{'scf'}->{'ncpus'} >= 1 );
  $npool = $hashRef->{'scf'}->{'npool'} if( exists $hashRef->{'scf'}->{'npool'} && $hashRef->{'scf'}->{'npool'} >= 1 );
  if( $ncpus % $npool == 0 ) {
    $ncpus = $ncpus / $npool;
    $npool = 1;
    $ncpus = 1;
  }
  if( $type eq 'potential' ) {
    $ncpus = 1;
    $npool = 1;
  }

  ### write PP input card for density
  open PP, ">", "$infile" or die "$!";
  print PP "&inputpp\n"
          . "  prefix = \'system\'\n"
          . "  outdir = \'Out\'\n"
          . "  filplot= \'$filplot\'\n"
          . $bonusInputs
          . "  plot_num = $flag\n"
          . "/\n";
  close PP;

  my $prefix = $hashRef->{'computer'}->{'para_prefix'};
  $prefix =~ s/$hashRef->{'computer'}->{'ncpus'}/$ncpus/ if( $hashRef->{'computer'}->{'ncpus'} != $ncpus );

  my $cmdLine = " -npool $npool ";

  my $outfile = $infile;
  $outfile =~ s/\.in/.out/;
  
#  my $didSwap = 0;
#  if($type eq 'potential' ) {
#    $didSwap = QEmggaFix();
#  }
  QErunPP( $hashRef->{'general'}->{'redirect'}, $prefix, $cmdLine, "$infile", "$outfile" );
#  if( $didSwap ) {
#    QEmggaUnFix();
#  }

  # Check for NaN with metaGGA
  if( $type eq 'potential' ) {
    my $error = QEfixPP( $hashRef->{'general'}->{'redirect'}, $prefix, $cmdLine, "$infile", "$outfile", "system.pot" );
    print "POT: $error\n";
    return $error if( $error != 0 );
  }

  system("$ENV{'OCEAN_BIN'}/qe2rhoofr.pl $convert" ) == 0 
    or die "Failed to convert $type\n$!\n";

  return 0;
}

sub QEfixPP
{
  my ( $redirect, $prefix, $cmdLine, $infile, $outfile, $potfile ) = @_;

  print "QE fix 1 $potfile\n";
  open IN, "<", $potfile or die "Failed to open $potfile\n$!";
  my $NAN = 0;
  while( my $line = <IN> ) {
    if( $line =~ m/nan/i ) {
      $NAN = 1;
      last;
    }
  }
  close IN;
  return 0 if( $NAN == 0 );
  print "QE fix 2\n";

  my $qe62File = catfile( "Out", "system.save", "data-file-schema.xml" );
  if( -e $qe62File ) {
    my @file;
    open IN, "<", $qe62File or die "$!";
    while (<IN>) {
      push @file, $_;
    }
    close IN;
    my $funct;
    open OUT, ">", $qe62File or die "$!";
    foreach my $line (@file) {
      if( $line =~ m/<functional>(.*)<\/functional>/ ) {
        $funct = $1;
        $line =~ s/$funct/PBE/;
      }
      print OUT $line;
    }
    close OUT;
    print "QE fix $funct\n";

    QErunPP( $redirect, $prefix, $cmdLine, $infile, $outfile );

    open IN, "<", $potfile or die "Failed to open $potfile\n$!";
    $NAN = 0;
    while( my $line = <IN> ) {
      if( $line =~ m/nan/i ) {
        $NAN = 1;
        last;
      }
    }
    close IN;
  
    open OUT, ">", $qe62File or die "$!";
    foreach my $line (@file) {
      if( $line =~ m/<functional>PBE<\/functional>/ ) {
        $line =~ s/PBE/$funct/;
      }
      print OUT $line;
    }
    close OUT;
  }

  return $NAN;
}

# Figure out how many pools to run with
sub QEPoolControl
{
  my ( $actualKpts, $numKS, $mem, $FFT, $OMP, $hashRef ) = @_;

  $OMP = 1 if( $OMP < 1 );
  my $maxMem = 2000*$OMP;
  my $minPool = 0.9;
  $minPool = $mem / $maxMem if( $mem > $maxMem );

  print "Memory estimate $mem\n";
  print "Min pool size $minPool\n";
  

#  print "$hashRef->{'ncpus'}\n";

  my %cpusAndPools;
  # j is number of processors
  for( my $j = $hashRef->{'ncpus'}; $j > 0; $j-- )
  {
    # i is number of pools
    my $maxNpool = $j;
    $maxNpool = $actualKpts if ( $actualKpts < $j );
    for( my $i = 1; $i <= $maxNpool; $i++ )
    {
      # gotta be an even factor
      next if ( $j % $i );
      # pool size (j/i) must be larger than minimum size
      next if ( $j / $i < $minPool );

      next if ( $numKS > 0 && $j / $i > $numKS );

      my $cpuPerPool = $j / $i;
      my $nbd = 1;
      if( $cpuPerPool > $FFT && $FFT > 0 ) {
        for( my $ii = 2; $ii <= $cpuPerPool/2; $ii++ ) {
          next if( $numKS % $ii || $cpuPerPool % $ii);
          $nbd = $ii;
          last if( $cpuPerPool/$nbd <= $FFT );
        }
        next if( $nbd == 1 );
      }
          
      my $kPerPool = ceil( $actualKpts / $i );
      my $cost = $kPerPool / $cpuPerPool;
      # Penalize multi-procs per pool since k-point parallelization is most efficient
#      print "$j  $i  $nbd $kPerPool  $cost";
      $cost /= ( 0.999**$cpuPerPool );
#      print "  $cost";
      # Penalize is not many bands per 
      $cost /= ( atan($numKS/$cpuPerPool) * 2.0/pi());
#      print "  $cost\n";
#      my $cost = $kPerPool / ( $cpuPerPool * ( 0.999**$cpuPerPool ) );
#      print "$j  $i  $nbd $kPerPool  $cost\n";
      $cpusAndPools{ $cost } = [ $j, $i, $nbd ];
    }
  }

  if( scalar keys %cpusAndPools < 1 ) {
    my $j = $hashRef->{'ncpus'};
    my $i = 1;
    my $cpuPerPool = $j;
    my $nbd = 1;
    if ( $numKS > 0 && $j / $i > $numKS ) { die "Memory won't fit" }
    if( $cpuPerPool > $FFT && $FFT > 0 ) {
      for( my $ii = 2; $ii <= $cpuPerPool/2; $ii++ ) {
        next if( $numKS % $ii || $cpuPerPool % $ii);
        $nbd = $ii;
        last if( $cpuPerPool/$nbd <= $FFT );
      }
      next if( $nbd == 1 );
    }

    my $kPerPool = ceil( $actualKpts / $i );
    my $cost = $kPerPool / ( $cpuPerPool * ( 0.999**$cpuPerPool ) );
#    print "$j  $i  $nbd $kPerPool  $cost\n";
    $cpusAndPools{ $cost } = [ $j, $i, $nbd ];
  }

  my $ncpus; my $npool; my $nbd;
  my @sortedCost = sort {$a <=> $b} keys %cpusAndPools;
  print " N procs     Pool     Band        Cost\n";
  for( my $i = 0; $i < 5; $i ++ ) {
    last if( $i >= scalar @sortedCost );
    printf "%8d %8d %8d       %g\n", $cpusAndPools{$sortedCost[$i]}[0], $cpusAndPools{$sortedCost[$i]}[1], $cpusAndPools{$sortedCost[$i]}[2], $sortedCost[$i]; 
  }
  foreach my $i (sort {$b <=> $a} keys %cpusAndPools )
  {
#    print "$i  $cpusAndPools{$i}[0]  $cpusAndPools{$i}[1]\n";
    $ncpus = $cpusAndPools{$i}[0];
    $npool = $cpusAndPools{$i}[1];
    $nbd = $cpusAndPools{$i}[2];
  }

  return ( $ncpus, $npool, $nbd );
}

# subroutine to handle running each QE step
# make one of these for pp and ph too!
sub QErunPW
{
  my( $redirect, $prefix, $cmdLine, $in, $out ) = @_;

  if( $redirect )
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine < $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine < $in > $out 2>&1");
  }
  else
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine -inp $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PW'} $cmdLine -inp $in > $out 2>&1");
  }
}


# subroutine to handle running each QE step
# make one of these for ph too!
sub QErunPP
{
  my( $redirect, $prefix, $cmdLine, $in, $out ) = @_;

  if( $redirect )
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PP'} $cmdLine < $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PP'} $cmdLine < $in > $out 2>&1");
  }
  else
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PP'} $cmdLine -inp $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PP'} $cmdLine -inp $in > $out 2>&1");
  }
}

# subroutine to handle running each QE step
sub QErunPH
{
  my( $redirect, $prefix, $cmdLine, $in, $out ) = @_;

  if( $redirect )
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PH'} $cmdLine < $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PH'} $cmdLine < $in > $out 2>&1");
  }
  else
  {
    print  "$prefix $ENV{'OCEAN_ESPRESSO_PH'} $cmdLine -inp $in > $out 2>&1\n";
    system("$prefix $ENV{'OCEAN_ESPRESSO_PH'} $cmdLine -inp $in > $out 2>&1");
  }
}


# Subroutine to print a QE-style input file
#   Pass in a file handle and a hashtable
sub QEprintInput
{
  my ($fh, $generalRef, $specificRef, $calcFlag, $kptString ) = @_;

  my $calc;
  my $tstress = '.false.';
  my $tprnfor = '.false.';
  my $nosyminv;
  my $startingPot;
  my $diagonalization;
  if( $calcFlag )
  {
    $calc = 'scf';
    $tstress = 'true' if( $generalRef->{'general'}->{'calc_stress'} );
    $tprnfor = 'true' if( $generalRef->{'general'}->{'calc_force'} );
    $nosyminv = 'false';
    $startingPot = 'atomic';
    if( $generalRef->{'general'}->{'startingpot'} =~ m/file/ )
    {
      if( -e catfile( "Out", "system.save", "charge-density.dat" ) || 
        -e catfile( "Out", "system.save", "charge-density.xml" ) ) {
        $startingPot = 'file';
        print "  WARNING: restarting from existing charge density!\n";
      }
    }
    if( $generalRef->{'general'}->{'diagonalization'} eq 'default' ) {
      $diagonalization = 'david';
    } else {
      $diagonalization = $generalRef->{'general'}->{'diagonalization'};
    }
  }
  else
  {
    $calc = 'nscf';
    $tstress = 'false';
    $tprnfor = 'false';
    $nosyminv = 'true';
    $startingPot = 'file';
    $diagonalization = $specificRef->{'diagonalization'};
    if( $generalRef->{'scf'}->{'version'} > 0 && $generalRef->{'scf'}->{'version'} < 6.5 ) {
      unless( $diagonalization eq 'cg' || $diagonalization eq 'david' ) {
        $diagonalization = 'david';
      }
    } else {
      printf "QE version %g\n", $generalRef->{'scf'}->{'version'};
    }
  }

  my $noncolin = '.false.';
  $noncolin = '.true.' if( $generalRef->{'general'}->{'noncolin'} );
  my $spinorb = '.false.';
  $spinorb = '.true.' if( $generalRef->{'general'}->{'spinorb'} );

  my $nbnd;
  if( exists $specificRef->{'nbands'} )
  {
    $nbnd = $specificRef->{'nbands'};
  } else {
    $nbnd = $generalRef->{'structure'}->{'valence_electrons'} / 2
          + $generalRef->{'general'}->{'fband'} * scalar @{$generalRef->{'structure'}->{'xred'}};
    $nbnd = ceil( $nbnd );
    $nbnd++ if( $nbnd % 2 );
  }
    

  # Array of QE names for smearing by occopt
  my @QE_smear;
  $QE_smear[1] = "'gaussian'";     # Still need to fix to be insulator
  $QE_smear[3] = "'fermi-dirac'";  # ABINIT = fermi-dirac
  $QE_smear[4] = "'marzari-vanderbilt'";  # ABINIT = Marzari cold smearing a = -0.5634
  $QE_smear[5] = "'marzari-vanderbilt'";  # ABINIT = Marzari a = -0.8165
  $QE_smear[6] = "'methfessel-paxton'";  # ABINIT = Methfessel and Paxton PRB 40, 3616 (1989)
  $QE_smear[7] = "'gaussian'";     # ABINIT = Gaussian

  my $occopt = "fixed";
  $occopt = "smearing" if( $generalRef->{'general'}->{'occopt'} != 1 );

  print $fh "&control\n"
        .  "  calculation = \'$calc\'\n"
        .  "  prefix = \'system\'\n"
        .  "  pseudo_dir = \'$generalRef->{'psp'}->{'ppdir'}\'\n"
        .  "  outdir = \'Out\'\n"
        .  "  wfcdir = \'$generalRef->{'general'}->{'tmp_dir'}\'\n"
        .  "  tstress = $tstress\n"
        .  "  tprnfor = $tprnfor\n"
        .  "  wf_collect = .true.\n"
        .  "  disk_io = 'low'\n";
  if( $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'control'} ne '' ) {
    my $temp = $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'control'};
    $temp =~ s/,/,\n/g;
    print $fh $temp . "\n";
#    print $fh $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'control'} . "\n";
  }
  print $fh "/\n"
        .  "&system\n"
        .  "  ibrav = 0\n"
        .  "  nat = " . scalar @{$generalRef->{'structure'}->{'xred'}} . "\n"
        .  "  ntyp = " . scalar @{$generalRef->{'structure'}->{'znucl'}} . "\n"
        .  "  noncolin = $noncolin\n" 
        .  "  lspinorb = $spinorb\n"
#        .  "  ecutwfc = $generalRef->{'general'}->{'ecut'}\n"
        . (sprintf "  ecutwfc = %i\n", $generalRef->{'general'}->{'ecut'})
        .  "  occupations = \'$occopt\'\n"
        .  "  smearing = $QE_smear[$generalRef->{'general'}->{'occopt'}]\n"
        . (sprintf "  degauss = %g\n  nspin  = %i\n  tot_charge = %g\n", 
              $generalRef->{'general'}->{'degauss'}, $generalRef->{'general'}->{'nspin'},
              $generalRef->{'structure'}->{'charge'})
        .  "  nosym = $nosyminv\n"
        .  "  noinv = $nosyminv\n";
  unless( $generalRef->{'general'}->{'functional'} =~ m/default/ )
  {
    print $fh "  input_dft = \'$generalRef->{'general'}->{'functional'}\'\n";
    print $fh "  nqx1 = $generalRef->{'general'}->{'exx'}->{'qmesh'}[0],"
             . " nqx2 = $generalRef->{'general'}->{'exx'}->{'qmesh'}[1],"
             . " nqx3 = $generalRef->{'general'}->{'exx'}->{'qmesh'}[2]\n";
  }
  if( exists $generalRef->{'general'}->{'isolated'} && 
      $generalRef->{'general'}->{'isolated'} ne 'none' ) {
    print $fh "  assume_isolated = \'$generalRef->{'general'}->{'isolated'}\'\n";
  }
  if( exists $generalRef->{'structure'}->{'magnetization'} &&
      $generalRef->{'structure'}->{'magnetization'} ne 'none' ) {
    printf $fh "  tot_magnetization = %g\n", $generalRef->{'structure'}->{'magnetization'};
  }

  print $fh "  nbnd = $nbnd\n";

  if( $generalRef->{'general'}->{'smag'}  ne "" )
  {
    print $fh "$generalRef->{'general'}->{'smag'}\n";
  }
  if( $generalRef->{'general'}->{'ldau'}->{'enable'} && 
      ( $calcFlag || not $generalRef->{'general'}->{'ldau'}->{'scf_only'} ) )
  {
    print $fh "lda_plus_u = true\n" 
            . "lda_plus_u_kind = $generalRef->{'general'}->{'ldau'}->{'lda_plus_u_kind'}\n"
            . "U_projection_type = '$generalRef->{'general'}->{'ldau'}->{'U_projection_type'}'\n";
    print $fh "$generalRef->{'general'}->{'ldau'}->{'Hubbard_U'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_U'} ne "" );
    print $fh "$generalRef->{'general'}->{'ldau'}->{'Hubbard_V'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_V'} ne "" );
    print $fh "$generalRef->{'general'}->{'ldau'}->{'Hubbard_J'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_J'} ne "" );
    print $fh "$generalRef->{'general'}->{'ldau'}->{'Hubbard_J0'}\n" 
        if( $generalRef->{'general'}->{'ldau'}->{'Hubbard_J0'} ne "" );
  }
#  if( $inputs{'qe_scissor'}  ne "" )
#  {
#    print $fh "$inputs{'qe_scissor'}\n";
#  }
#  if( $inputs{'ibrav'} != 0 )
#  {
#    print $fh "  celldim(1) = $inputs{'celldm1'}\n";
#  }

  if( $generalRef->{'general'}->{'mixing_mode'} =~ m/local-tf/i ) {
    $generalRef->{'general'}->{'mixing_mode'} = 'local-TF';
    print "local-TF\n";
  } elsif( $generalRef->{'general'}->{'mixing_mode'} =~ m/tf/i ) {
    $generalRef->{'general'}->{'mixing_mode'} = 'TF';
    print "TF\n";
  } else {
    $generalRef->{'general'}->{'mixing_mode'} = 'plain';
    print "local\n";
  }

  if( $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'system'} ne '' ) {
    my $temp = $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'system'};
    $temp =~ s/,/,\n/g;
    print $fh $temp . "\n";
#    print $fh $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'system'} . "\n";
  }

  print $fh "/\n"
        .  "&electrons\n"
        .  (sprintf "  conv_thr = %g\n", $specificRef->{'toldfe'})
        .  (sprintf "  mixing_beta = %g\n", $generalRef->{'general'}->{'mixing'}) 
        .  (sprintf "  mixing_mode = '%s'\n", $generalRef->{'general'}->{'mixing_mode'})
        .  (sprintf "  electron_maxstep = %i\n", $generalRef->{'general'}->{'nstep'})
        .  "  startingwfc = \'$generalRef->{'general'}->{'startingwfc'}\'\n"
        .  "  startingpot = \'$startingPot\'\n"
        .  "  diagonalization = \'$diagonalization\'\n";

  if( $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'electrons'} ne '' ) {
    my $temp = $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'electrons'};
    $temp =~ s/,/,\n/g;
    print $fh $temp . "\n";
#    print $fh $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'electrons'} . "\n";
  }

  print $fh "/\n";
#  if( $inputs{'nscfEXX'} == 1 )
#  {
#    # Since (at the moment) we are loading the SCF density
#    #  don't converge the density for the first iteration w/o EXX
#    print $fh "  adaptive_thr = .true., conv_thr_init = 1\n";
#  }

  print $fh "ATOMIC_SPECIES\n";# . $inputs{'atompp'} . "\n";
  for( my $a=0; $a < scalar @{$generalRef->{'structure'}->{'zsymb'}}; $a++ )
  {
    printf $fh "%4s  0.0  %s\n", $generalRef->{'structure'}->{'zsymb'}[$a], $generalRef->{'psp'}->{'pp_list'}[$a];
  }

#  if ($inputs{'ibrav'} == 0) {
    print $fh "CELL_PARAMETERS cubic\n";  # . $inputs{'acell'} . "\n";
    for( my $i = 0; $i<3; $i++ )
    {
      printf $fh " %.16f    %.16f   %.16f \n", $generalRef->{'structure'}->{'avecs'}[$i][0], 
          $generalRef->{'structure'}->{'avecs'}[$i][1], $generalRef->{'structure'}->{'avecs'}[$i][2];
    }
#  }

  print $fh "ATOMIC_POSITIONS crystal\n";

  for( my $a=0; $a < scalar @{$generalRef->{'structure'}->{'xred'}}; $a++ )
  {
    printf $fh "%4s  %.16f    %.16f    %.16f \n", 
          $generalRef->{'structure'}->{'zsymb'}[$generalRef->{'structure'}->{'typat'}[$a]-1], 
          $generalRef->{'structure'}->{'xred'}[$a][0],
          $generalRef->{'structure'}->{'xred'}[$a][1],
          $generalRef->{'structure'}->{'xred'}[$a][2];
  }

  if( exists $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'hubbard'}  &&
        $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'hubbard'}  ne '' ) 
  {
      my $s = $generalRef->{'general'}->{'verbatim'}->{'qe'}->{'hubbard'};
      $s =~ s/\s+V /\nV /ig;
      $s =~ s/\s+U /\nU /ig;
      $s =~ s/\s+J /\nJ /ig;
      $s =~ s/\s+J0 /\nJ0 /ig;
      $s =~ s/\s+B /\nB /ig;
      print $fh $s . "\n";
  }


  print $fh $kptString;

}

#Move this to standalone at some point for re-use with ABINIT
sub kptGen
{
  my ( $hashRef, $dualKpt ) = @_;

  my $kptText = "";
  my @q = ( 0, 0, 0 );
  if( $dualKpt ) {
    $q[0] = $hashRef->{'photon_q'}[0];
    $q[1] = $hashRef->{'photon_q'}[1];
    $q[2] = $hashRef->{'photon_q'}[2];
  }

 
  for( my $x = 0; $x < $hashRef->{'kmesh'}[0]; $x++ ) {
    my $xk = $hashRef->{'kshift'}[0]/$hashRef->{'kmesh'}[0] + $x/$hashRef->{'kmesh'}[0] - $q[0];
    while( $xk > 1 ) { $xk -= 1.0; }
    while( $xk < -1 ) { $xk += 1.0; }
    my $xk2 = $hashRef->{'kshift'}[0]/$hashRef->{'kmesh'}[0] + $x/$hashRef->{'kmesh'}[0];
    while( $xk2 > 1 ) { $xk2 -= 1.0; }
    while( $xk2 < -1 ) { $xk2 += 1.0; }
    for( my $y = 0; $y < $hashRef->{'kmesh'}[1]; $y++ ) {
      my $yk = $hashRef->{'kshift'}[1]/$hashRef->{'kmesh'}[1] + $y/$hashRef->{'kmesh'}[1] - $q[1];
      while( $yk > 1 ) { $yk -= 1.0; }
      while( $yk < -1 ) { $yk += 1.0; }
      my $yk2 = $hashRef->{'kshift'}[1]/$hashRef->{'kmesh'}[1] + $y/$hashRef->{'kmesh'}[1];
      while( $yk2 > 1 ) { $yk2 -= 1.0; }
      while( $yk2 < -1 ) { $yk2 += 1.0; }
      for( my $z = 0; $z < $hashRef->{'kmesh'}[2]; $z++ ) {
        my $zk = $hashRef->{'kshift'}[2]/$hashRef->{'kmesh'}[2] + $z/$hashRef->{'kmesh'}[2] - $q[2];
        while( $zk > 1 ) { $zk -= 1.0; }
        while( $zk < -1 ) { $zk += 1.0; }
        my $zk2 = $hashRef->{'kshift'}[2]/$hashRef->{'kmesh'}[2] + $z/$hashRef->{'kmesh'}[2];
        while( $zk2 > 1 ) { $zk2 -= 1.0; }
        while( $zk2 < -1 ) { $zk2 += 1.0; }

        $kptText .= sprintf "%19.15f  %19.15f  %19.15f  1\n", $xk, $yk, $zk;
        if( $dualKpt ) { 
          $kptText .= sprintf "%19.15f  %19.15f  %19.15f  1\n", $xk2, $yk2, $zk2;
        }
      }
    }
  }

  
  return $kptText;
}

sub QEparseEnergies
{
  my ($hashRef, $specificHashRef ) = @_;

  my $spin = $hashRef->{'general'}->{'nspin'};
  my $occopt = $hashRef->{'general'}->{'occopt'};
  my $insulator = 1;
  $insulator = 0 if( $spin > 1 || $occopt != 1 );
  my $split = $specificHashRef->{'split'};
#  $split = 1 if ( $split );
  my $nonzero_q = $specificHashRef->{'nonzero_q'};

  my $fermi = $hashRef->{'scf'}->{'fermi'};

  my @energies;
  
  my $qe62File = catfile( "Out", "system.save", "data-file-schema.xml" );
  my $qe54File = catfile( "Out", "system.save", "data-file.xml" );
  
  if( -e $qe62File )
  {
    my $splitFile = $qe62File;
    $splitFile = catfile( "Out_shift", "system.save", "data-file-schema.xml" ) if( $split );
    @energies = QEparseEnergies62( $qe62File, $splitFile, $spin );
  }
  elsif( -e $qe54File )
  {
    my $splitFile = $qe54File;
    $splitFile = catfile( "Out_shift", "system.save", "data-file.xml" ) if( $split );
    @energies = QEparseEnergies54( $qe54File, $splitFile, $spin );
  }
  else
  {
    die;
    return 1;
  }

  my $nks = $spin * $specificHashRef->{'kmesh'}[0] * $specificHashRef->{'kmesh'}[1]
                  * $specificHashRef->{'kmesh'}[2];
  $nks *=2 if( $nonzero_q );

  if( $nks != scalar @energies )
  {
    print "Wrong number of spins and k-points!  Expected: $nks . Found: " . scalar @energies . "\n";
    return 2;
  }

  my @b;
  $b[0] = 1;
  if( scalar @energies > 1 )
  {
    $b[3] = scalar @{$energies[1]};
  } else {
    $b[3] = scalar @{$energies[0]};
  }

  if( $insulator )
  {
    my $nelec = sprintf "%.0f", $hashRef->{'scf'}->{'nelec'};
    die "Fractional nuber of electrons, but fixed occupations\n" 
        if( abs( $nelec - $hashRef->{'scf'}->{'nelec'} ) > 0.001 );
    die "Odd electrons, but fixed occupations\n" if( $nelec % 2 );
    $b[1] = $nelec / 2 ;
    $b[2] = $b[1] + 1;
  }
  else
  {
    $b[1] = 1;
    $b[2] = $b[3];

    for( my $k = 0; $k < scalar @energies; $k++ )
    {
#      my $temp = 1; #$b[1] - 1;
      my $temp = $b[1] - 1;
      for( my $i = $temp; $i < scalar @{$energies[$k]}; $i++ )
      {
        $b[1] = $i+1 if( $energies[$k][$i] < $fermi );
        last if( $energies[$k][$i] > $fermi );
#        print "vvv $b[1]\n";
      }

      $temp = $b[2] - 1;
      for( my $i = $temp; $i >= 0; $i-- )
      {
        $b[2] = $i+1 if( $energies[$k][$i] > $fermi );
        last if( $energies[$k][$i] < $fermi );
#        print "ccc $b[2]\n";
      }
    }
  }

  # Store max band of occupied and min band of unoccupied
  #  before storing brange
  $specificHashRef->{'max_occ_band'} = $b[1];
  $specificHashRef->{'min_unocc_band'} = $b[2];

  if( exists( $specificHashRef->{'con_start'} ) ) {
    $b[2] = $specificHashRef->{'con_start'} if( $specificHashRef->{'con_start'} >= 1 );
  }

  if( exists( $specificHashRef->{'val_stop'} ) ) {
    $b[1] = $specificHashRef->{'val_stop'} if( $specificHashRef->{'val_stop'} >= 1 );
    $b[1] = $b[3] if( $b[1] > $b[3] );
  }

  $specificHashRef->{'brange'} = \@b ;
  
  open OUT, ">", "QE_EIGS.txt" or die;
  my $nk = $specificHashRef->{'kmesh'}[0] * $specificHashRef->{'kmesh'}[1] * $specificHashRef->{'kmesh'}[2];
  printf OUT "%i %i %i %i %i %i\n", $b[0], $b[1], $b[2], $b[3], $nk, $spin;

  my $delim = " ";
  my $n = 3;
  my $ik = 0;
  my $nq = 0;
  $nq ++ if( $nonzero_q );
  for( my $s = 0; $s < $spin; $s++ )
  {
    for( my $k = 0; $k < $nk; $k++ )
    {
      my @eslice = @{ $energies[$ik] }[ $b[0]-1 .. $b[1]-1 ];
      foreach my $x (@eslice) { $x = $x * 2; }
      while (my @x = splice @eslice, 0, $n)
      {
        print OUT join($delim, @x), "\n";
      }
      $ik++ if( $nonzero_q );
      my @eslice = @{ $energies[$ik] }[ $b[2]-1 .. $b[3]-1 ];
      foreach my $x (@eslice) { $x = $x * 2; } 
      while (my @x = splice @eslice, 0, $n)
      { 
        print OUT join($delim, @x), "\n";
      } 
      $ik++;
    }
  }
  close OUT;
  
  

  return 0;
}

sub QEparseEnergies62
{
  my ($f1, $f2, $spin) = @_;

  my $split = 0;
  open my $fh1, "<", $f1 or die "Failed to open $f1\n$!";

  my @energies1;
  my @energies1_spin;

  while( my $line = <$fh1> )
  {

    if( $line =~ m/<eigenvalues size=\"\d+\">/ )
    {
      $line =~ s/<eigenvalues size=\"\d+\">//;
      chomp $line;
      my $eigs = $line . ' ';
      until( $line =~ m/eigenvalues/ )
      { 
        $line = <$fh1>;
        chomp $line;
        $eigs .= $line . ' ';
      }
      $eigs =~ s/\s+<\/eigenvalues>//;

      my @eigs = split( ' ', $eigs );
      if( $spin == 2 )
      {
        my $half = scalar @eigs / 2;
        my @t1 = @eigs[ 0..$half-1 ];
        push @energies1, \@t1;
        my @t2 = @eigs[ $half..scalar @eigs-1 ];
        push @energies1_spin, \@t2;
#        print scalar @eigs . "  $half  " . scalar @t1 . "  " .scalar @t2 . "\n";
      }
      else
      {
        push @energies1, \@eigs;
      }
    }
  }
  close $fh1;

  unless( $f1 ne $f2 )
  {
#    push @energies1, @energies1_spin if( $spin == 2 );
    
    if( $spin == 2 ) {
      print scalar @energies1 . "\n";
      print scalar @energies1_spin . "\n";
      for( my $i = 0; $i < scalar @energies1_spin; $i++ ) {
        push @energies1, \@{ $energies1_spin[$i] };
      }
      print scalar @energies1 . "\n";
    }
    return @energies1;
  }


  
  open my $fh2, "<", $f2 or die "Failed to open $f2\n$!";
  my @energies2;
  my @energies2_spin;

  while( my $line = <$fh2> )
  {
    
    if( $line =~ m/<eigenvalues size=\"\d+\">/ )
    { 
      $line =~ s/<eigenvalues size=\"\d+\">//;
      chomp $line;
      my $eigs = $line . ' ';
      until( $line =~ m/eigenvalues/ )
      { 
        $line = <$fh2>;
        chomp $line;
        $eigs .= $line . ' ';
      }
      $eigs =~ s/\s+<\/eigenvalues>//;
      
      my @eigs = split( ' ', $eigs );
      if( $spin == 2 )
      {
        my $half = scalar @eigs / 2;
        my @t1 = @eigs[ 0..$half-1 ];
        push @energies2, \@t1;
        my @t2 = @eigs[ $half..scalar @eigs-1 ];
        push @energies2_spin, \@t2;
      }
      else
      {
        push @energies2, \@eigs;
      }
    }
  }
  close $fh2;

  if( scalar @energies1 != scalar @energies2 )
  { 
    die "K-point mismatch between split DFT runs\n";
  }
  my @energies;
  for( my $i = 0; $i < scalar @energies1; $i++ )
  {
    #TODO: Not sure this is necessary or even works
    push @energies, \@{ $energies1[$i] };
    push @energies, \@{ $energies2[$i] };
  }
  if( $spin == 2 )
  {
    for( my $i = 0; $i < scalar @energies1; $i++ )
    {
      #TODO: Not sure this is necessary or even works
      push @energies, \@{ $energies1_spin[$i] };
      push @energies, \@{ $energies2_spin[$i] };
    }
  }

  return @energies;
  
}

sub QEparseEnergies54
{
  my ($f1, $f2, $spin) = @_;

  my $split = 0;
  open my $fh1, "<", $f1 or die "Failed to open $f1\n$!";

#  print ">>>>>>>>>> $f1 <<<<<<<<<<\n";
  
  my @energies1;
  my @energies1_spin;

  while( my $line = <$fh1> ) 
  {
#    print $line;
#    if( $line =~ m/iotk_link=\"(\S+eigenval.xml)\"/ ) {
    if(  $line =~ m/\"(\.\/K\d+\/eigenval.xml)/ ) {
      my $k = $1;
      my $eigFile = $f1;
      $eigFile =~ s/data-file.xml/$k/;
      my $eigs = '';
      open IN, $eigFile or die "Failed to open $eigFile\n$!";
      while( $line = <IN> ) {
        if( $line =~ m/EIGENVALUES/ ) {
          $line = <IN>;
          until ( $line =~ m/\/EIGENVALUES/ ) {
            chomp $line;
            $eigs .= $line . ' ';
            $line = <IN>;
          }
          my @eigs = split( ' ', $eigs );
#          print "$eigFile  " . scalar @eigs . "\n";
          if( $spin == 2 )
          {
            my $half = scalar @eigs / 2;
            my @t1 = @eigs[ 0..$half-1 ];
            push @energies1, \@t1;
            my @t2 = @eigs[ $half..scalar @eigs-1 ];
            push @energies1_spin, \@t2;
          }
          else
          {
            push @energies1, \@eigs;
          }
        }
      }
      close IN;
    }
  }
  close $fh1;
  unless( $f1 ne $f2 )
  {
    if( $spin == 2 ) {
      print scalar @energies1 . "\n";
      print scalar @energies1_spin . "\n";
      for( my $i = 0; $i < scalar @energies1_spin; $i++ ) {
        push @energies1, \@{ $energies1_spin[$i] };
      }
      print scalar @energies1 . "\n";
    }
    return @energies1;
  }

  die "QE54 w/ shift not implemented. Upgrade your QE version to something from the past 5 years\n";
}


sub QEconvertEnergies
{
  my ($Nsemicore) = @_;
  my $emin = 0;
  my $emax = 0;

  print "########  QE convert\n";

  my $xmlFile = catfile( "Out", "system.save", "data-file-schema.xml" );
  if( -e $xmlFile ) {

    open IN, $xmlFile or die "Failed to open $xmlFile\n$!";
    while (my $line = <IN> ) {
      last if ( $line =~m/band_structure/ );
    }
    my $nk = 0;
    my $nspin = 0;
    my $nbands = 0;
    my $nband_up = 0;
    my $nband_dw = 0;
    my @energies;
    my @weights;
    my @kpts;
    my $fermi;
    while( my $line = <IN> ) { 
      if( $line =~ m/<lsda>(\w+)</ ) {
        if( $1 =~ m/true/ ) {
          $nspin = 2;
        } else {
          $nspin = 1;
        }
      } elsif( $line =~ m/<nbnd>(\d+)</ ) {
        $nbands = $1;
      } elsif( $line =~ m/<nbnd_up>(\d+)</ ) {
        $nband_up = $1;
      } elsif($line =~ m/<nbnd_dw>(\d+)</ ) {
        $nband_dw = $1;
      } elsif( $line =~ m/<nks>(\d+)</ ) {
        $nk = $1;
      } elsif( $line =~ m/<fermi_energy>(-?\d+\.\d+[eE]-?\d+)</ ) {
        $fermi = $1;
      } elsif( $line =~ m/<ks_energies>/ ) {
        <IN> =~ m/<k_point weight=\"(-?\d+\.\d+[eE]-?\d+)\">(-?\d+\.\d+[eE]-?\d+)\s+(-?\d+\.\d+[eE]-?\d+)\s+(-?\d+\.\d+[eE]-?\d+)</;
        push @weights, $1;
        push @kpts, [ $2, $3, $4];
        <IN>;
        <IN>;
        $line = <IN>;
        my $temp;
        until ($line =~m/<\/eigenvalues>/ ) {
          $temp .= $line;
          $line = <IN>;
        }
        my @temp = split ' ', $temp;
        push @energies, \@temp;
        until( $line =~m/<\/ks_energies>/ ) {
          $line = <IN>;
        }
      }

    }
    close IN;

    $nbands = $nband_up + $nband_dw if( $nbands == 0 );

    my @allEnergies;
    open OUT, ">QE.txt";
    printf OUT "%i %i %i\n", $nbands, $nk, $nspin;
#    for( my $spin = 0; $spin < $nspin; $spin++ ) {
      for( my $k = 0; $k < $nk; $k++ ) {
        printf OUT "%.15e %.15e %.15e %.15e\n", $kpts[$k][0], $kpts[$k][1],$kpts[$k][2], $weights[$k];
        for( my $b = 0; $b < $nbands; $b++ ) {
          printf OUT "%.15e\n", $energies[$k][$b];
          push @allEnergies, $energies[$k][$b];
        }
      }
#    } 
    close OUT;

    # keep emin and emax equal to each other to detect errors late
    if( defined $fermi ) {
      $emax = $fermi*27.211386245988;
      $emin = $emax;

      # here, if fermi isn't defined we won't bother detecting emin
      if( $Nsemicore > 0 ) {
        my @allSorted = sort { $a <=> $b } @allEnergies;
        my $n = $Nsemicore * $nk * $nspin / 2;
        printf "Semicore energy 'gap': %.4f  %.4f  %.4f  %.4f\n", $allSorted[$n-1]*27.211386245988, $allSorted[$n]*27.211386245988, $allSorted[$n-1], $allSorted[$n];
        $emin = ( $allSorted[$n-1] + $allSorted[$n] )*27.211386245988/2;
      }
    }
  } else {
    print "WARNING! QEconvertEnergiers called but $xmlFile is missing\n";
  }
  return( $emin, $emax);
}

sub QEmggaFix {
  my $sub = 'PBE';
  my $didSwap = 0;
  my $funct = '';

  my $xmlFile = catfile( 'Out', 'system.save', 'data-file-schema.xml' );
  my $xmlFile2 = catfile( 'Out', 'system.save', 'data-file-schema.xml_orig' );
  if( -e $xmlFile ) {
    copy $xmlFile, $xmlFile2;
    my $wholeFile;
    open IN, $xmlFile or die "$!";
    while( my $line = <IN> ) {
      if( $line =~ m/<functional>(.*)<\/functional>/ ) {
        $funct = $1;
        if( $funct =~ m/XC-(\d+)\w-(\d+)\w-(\d+)\w-(\d+)\w-(\d+)\w-(\d+)\w/ ) {
          if( $5 > 0 || $6 > 0 ) {
            $didSwap = 1;
            $line =~ s/$funct/$sub/;
          }
        } elsif( $funct =~ m/SCAN/i ) {
          $didSwap = 1;
          $line =~ s/$funct/$sub/;
        }
      }
      $wholeFile .= $line;
    }
    close IN;

    if( $didSwap ) {
      open OUT, ">", $xmlFile;
      print OUT $wholeFile;
      close OUT;

      print "Detected m-GGA functional $funct. Swapping with $sub to avoid QE bugs\n"; 
    }
  }

  return $didSwap;
}

sub QEmggaUnFix {
  my $xmlFile = catfile( 'Out', 'system.save', 'data-file-schema.xml' );
  my $xmlFile2 = catfile( 'Out', 'system.save', 'data-file-schema.xml_orig' );
  move $xmlFile2, $xmlFile;
}

1;
