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

sub ABIrunNSCF 
{
  my ( $hashRef, $specificHashRef, $shift ) = @_;

  my $dirname = sprintf "k%i_%i_%iq%.6f_%.6f_%.6f", $specificHashRef->{'kmesh'}[0],
                    $specificHashRef->{'kmesh'}[1], $specificHashRef->{'kmesh'}[2],
                    $specificHashRef->{'kshift'}[0], $specificHashRef->{'kshift'}[1],
                    $specificHashRef->{'kshift'}[2];

  print "$dirname\n";

  mkdir $dirname unless( -d $dirname );
  chdir $dirname;
  
  ABIsetupNSCF();

  my $nk = $specificHashRef->{'kmesh'}[0] * $specificHashRef->{'kmesh'}[1]
           * $specificHashRef->{'kmesh'}[2];
  $nk *=2  if( $specificHashRef->{'nonzero_q'} && not $specificHashRef->{'split'});
  my $kptString = sprintf "kptopt 0\nnkpt %i\nkpt\n", $nk;
  $kptString .= kptGenAbi( $specificHashRef, $shift );

  open my $input, ">", "nscf.in" or die "Failed to open nscf.in.\n$!";
  open my $files, ">", "nscf.files" or die "Failed to open nscf.files\n$!";

  my $calcFlag = 0;
  ABIprintFiles( $files, $hashRef, $calcFlag );
  ABIprintInput( $input, $hashRef, $specificHashRef, $calcFlag, "nscf", $kptString ) ;

  close $input;
  close $files;

  ABIrunABINIT( $hashRef->{'computer'}->{'para_prefix'}, "nscf.files", "nscf.log" );

  $specificHashRef->{'nelec'} = $hashRef->{'scf'}->{'nelec'};
  $specificHashRef->{'fermi'} = $hashRef->{'scf'}->{'fermi'};

  my $errorCode = ABIparseOut( "nscf.out", $specificHashRef );
  if( $errorCode != 0 ) {
    print "NSCF stage has error:  $errorCode\n";
    return $errorCode;
  }

  $errorCode = ABIbrange( "NSCFx_EIG", $specificHashRef, $hashRef->{'general'}->{'abpad'} );

  chdir updir();



  return $errorCode;
}

sub ABIsetupNSCF 
{
  my $prefix = 'SCFx_';
  my @fileList = ( 'DEN' );
  my @optfileList = ( 'KDEN' );

  foreach my $f (@fileList) {
    copy catfile( updir(), $prefix . $f ), $prefix . $f or die "Failed to copy ${prefix}${f}\n$!";
  }

  foreach my $f (@optfileList) {
    if( -e catfile( updir(), $prefix . $f ) ) {
      copy catfile( updir(), $prefix . $f ), $prefix . $f or die "Failed to copy ${prefix}${f}\n$!";
    }
  }
}

sub ABIrunDensity 
{
  my $hashRef = $_[0];

  open my $input, ">", "scf.in" or die "Failed to open scf.in.\n$!";
  open my $files, ">", "scf.files" or die "Failed to open scf.files\n$!";

  my $calcFlag = 1;

  my $kptString = sprintf "ngkpt %i %i %i\nnshiftk 1\nshiftk", $hashRef->{'scf'}->{'kmesh'}[0],
                        $hashRef->{'scf'}->{'kmesh'}[1], $hashRef->{'scf'}->{'kmesh'}[2];
  for( my $i = 0; $i < 3; $i++ ) {
    if( $hashRef->{'scf'}->{'kshift'}[$i]==1 ) {
      $kptString .= " 0.5";
    } else {
      $kptString .= " 0.0";
    }
  }
  $kptString .= "\n";

  # Modify QEprintInput to take two different hashes and a string 
  ABIprintFiles( $files, $hashRef, $calcFlag );
  ABIprintInput( $input, $hashRef, $hashRef->{'scf'}, $calcFlag, "scf", $kptString ) ;

  close $input;
  close $files;

  ABIrunABINIT( $hashRef->{'computer'}->{'para_prefix'}, "scf.files", "scf.log" );

  my $errorCode = ABIparseOut( "scf.out", $hashRef->{'scf'} );
  return $errorCode;
#  return 1;
}


#sub QErunTest {}

sub ABIbrange
{
  my ($file, $specificHashRef, $pad ) = @_;
  
  my $nkpt = $specificHashRef->{'kmesh'}[0]*$specificHashRef->{'kmesh'}[1]*$specificHashRef->{'kmesh'}[2];
  open IN, "<", $file or die "Failed to open $file\n$!";
  <IN> =~ m/nkpt=\s+(\d+)/ or die "Failed to parse $file\n";
  if( $1 != $nkpt ) {
    print "Mismatch between k-points in $file and expected $nkpt\n";
    return 11;
  }
  my @brange;
  $brange[0] = 1;
  $brange[1] = 1;
  
  <IN>;
  for( my $k = 0; $k < $nkpt; $k++ ) {
    my @eig;
    while( my $line = <IN> ) {
      last if ($line =~ m/kpt/ );
      @eig = (@eig, split ' ', $line );
    }
    if( $k==0 ) {
      $brange[2] = scalar @eig;
      $brange[3] = scalar @eig - $pad;
    }
    for( my $i = $brange[1]-1; $i < $brange[3]; $i++ ) {
      if( $eig[$i] < $specificHashRef->{'fermi'} ) {
        $brange[1] = $i+1;
      } else { 
        last; 
      }
    }
    for( my $i = $brange[2]-1; $i >= 0; $i-- ) {
      if( $eig[$i] > $specificHashRef->{'fermi'} ) {
        $brange[2] = $i+1;
      } else {
        last;
      }
    }
  }

  $specificHashRef->{'brange'} = \@brange;
  return 0;
}

sub ABIparseOut 
{
  my( $outFile, $hashRef ) = @_;

  my $errorCode = 1;

  my $fermi = 'no';
  my @occ;
  my $nband;

  open SCF, "<", $outFile or die "Failed to open $outFile\n$!";

  while( my $scf_line = <SCF> )
  { 
    if( $scf_line  =~  m/Fermi \(or HOMO\) energy \(hartree\) =\s+([+-]?\d+\.?\d+)/ )
    {
      $fermi = $1;
      $hashRef->{'fermi'} = $fermi*1;
    }
    elsif( $scf_line =~ m/Version\s+(\d+\.\d+\.\d+)\s+of ABINIT/ ) {
      $hashRef->{'version'} = $1;
    }
    elsif( $scf_line =~ m/wall=\s+(\d+\.\d+[eE]?[+-]?\d?)/ ) {
      $hashRef->{'time'} = $1*1;
    }
    elsif( $scf_line =~ m/etotal\s+(-?\d+\.\d+[eE]?[-+]?\d*)/ )
    { 
      $hashRef->{'etot'} = $1*1;
    }
#    elsif( $scf_line =~ m/is converged/ )
    elsif( $scf_line =~ m/Delivered\s+(\d+)\s+WARNING/ )
    {
      $errorCode = $1;
    }
    elsif( $scf_line =~ m/nband\s+(\d+)/ )
    {
      $nband = $1;
    }
    elsif( $scf_line =~ m/occ\s+\d/ && scalar @occ == 0 )
    {
      @occ = split ' ', $scf_line;
      shift @occ;
      $scf_line = <SCF>;
      while( $scf_line =~ m/^\s+\d/ ) {
        @occ = (@occ, split ' ', $scf_line );
        $scf_line = <SCF>;
      }
    }
#    elsif( $scf_line =~ m/nelect=\s+(\d+\.?\d*)/ ) 
#    {
#      $hashRef->{'nelec'} = $1*1;
#    }

  }

  close SCF;

  if( scalar @occ > 0 ) {
    $hashRef->{'nelec'} = 0;
    for( my $i = 0; $i<$nband; $i++ ) {
      $hashRef->{'nelec'} += $occ[$i];
    }
  }

  return $errorCode;
}

sub ABIparseDensityPotential
{
  my ( $hashRef, $type ) = @_;

  my $infile;
  my $datafile;
  my $outfile;
  my $logfile;
  my $nfft;
  if( $type eq 'density' ) {
    $infile = "cut3d.in";
    $datafile = "SCFx_DEN";
    $outfile = "rhoofr";
    $logfile = 'cut3d.log';
    $nfft = "nfft";
  } elsif( $type eq 'potential' ) {
    $infile = "cut3d2.in";
    $datafile = "SCFx_POT";
    $outfile = "potofr";
    $logfile = 'cut3d2.log';
    $nfft = "nfft.pot";
  } else {
    return 1;
  }
  
  
  my $AbiVersion = 9;
  if( $hashRef->{'scf'}->{'version'} =~ m/(\d+)\.(\d+)\.(\d+)/ ) {
    $AbiVersion = $1;
  }

  open CUTIN, ">", $infile or die "Failed to open $infile for writing.\n$!\n";
  if( $AbiVersion <= 7 )
  {
    if( $hashRef->{'general'}->{'nspin'} == 2 )
    {
      print CUTIN "$datafile\n1\n0\n6\n$outfile\n0\n";
    }
    else
    {
      print CUTIN "$datafile\n1\n6\n$outfile\n0\n";
    }
  }
  else
  {
    if(  $hashRef->{'general'}->{'nspin'}  == 2 )
    {
      print CUTIN "$datafile\n0\n6\n$outfile\n0\n";
    }
    else
    {
      print CUTIN "$datafile\n6\n$outfile\n0\n";
    }
  }
  close CUTIN;
  

  print "$hashRef->{'computer'}->{'ser_prefix'} $ENV{'OCEAN_CUT3D'} < $infile > $logfile 2>&1\n";
  unless (system("$hashRef->{'computer'}->{'ser_prefix'} $ENV{'OCEAN_CUT3D'} < $infile > $logfile 2>&1") == 0)
  {
    print "Failed to run cut3d\n";
    return 2;
  }

  unless( system( "tail -n 1 $outfile > $nfft" ) == 0 ) {
    print "Failed to make $nfft\n";
    return 3;
  }
   
      
  return 0;
}


#sub QEPoolControl {}

sub ABIrunABINIT
{
  my( $prefix, $in, $out ) = @_;

  print  "$prefix $ENV{'OCEAN_ABINIT'} < $in > $out 2>&1\n";
  system("$prefix $ENV{'OCEAN_ABINIT'} < $in > $out 2>&1");
  
}

#sub QErunPP {} 

sub ABIprintFiles
{
  my ($files, $generalRef, $calcFlag ) = @_;

  if( $calcFlag == 1 ) {
    print $files "scf.in\nscf.out\nSCF\nSCFx\nSCFxx\n"; 
  } else {
    print $files "nscf.in\nnscf.out\nSCFx\nNSCFx\nNSCFxx\n";
  }
  foreach my $p ( @{$generalRef->{'psp'}->{'pp_list'}}) {
    print $files catfile( $generalRef->{'psp'}->{'ppdir'}, $p ) . "\n";
  }
}

sub ABIprintInput 
{
  my ($input, $generalRef, $specificRef, $calcFlag, $baseString, $kptString ) = @_;


#  print $files $baseString . ".in\n" . $baseString . ".out\n" 
#          . uc($baseString) . "\n" . uc($baseString) . "x\n" 
#          . uc($baseString) . "xx\n";
#  foreach my $p ( @{$generalRef->{'psp'}->{'pp_list'}}) {
#    print $files catfile( $generalRef->{'psp'}->{'ppdir'}, $p ) . "\n";
#  }
  

  print $input "symmorphi 0\nautoparal 1\nchksymbreak 0\n"
             . "rprim\n";
  printf $input "% .16g  % .16g  % .16g\n% .16g  % .16g  % .16g\n% .16g  % .16g  % .16g\n", 
        $generalRef->{'structure'}->{'avecs'}[0][0], $generalRef->{'structure'}->{'avecs'}[0][1],
        $generalRef->{'structure'}->{'avecs'}[0][2], $generalRef->{'structure'}->{'avecs'}[1][0],
        $generalRef->{'structure'}->{'avecs'}[1][1], $generalRef->{'structure'}->{'avecs'}[1][2],
        $generalRef->{'structure'}->{'avecs'}[2][0], $generalRef->{'structure'}->{'avecs'}[2][1],
        $generalRef->{'structure'}->{'avecs'}[2][2];

  printf $input "ntypat %i\nznucl", scalar @{$generalRef->{'structure'}->{'znucl'}};
  foreach my $z (@{$generalRef->{'structure'}->{'znucl'}}) {
    printf $input " %i", $z;
  }
  printf $input "\nnatom %i\ntypat", scalar @{$generalRef->{'structure'}->{'typat'}};
  foreach my $z (@{$generalRef->{'structure'}->{'typat'}}) {
    printf $input " %i", $z;
  }
  print $input "\nxred\n";
  for( my $i = 0; $i < scalar @{$generalRef->{'structure'}->{'typat'}}; $i++ ) {
    printf $input " % .16g  % .16g  % .16g\n", $generalRef->{'structure'}->{'xred'}[$i][0], 
        $generalRef->{'structure'}->{'xred'}[$i][1], $generalRef->{'structure'}->{'xred'}[$i][2];
  }

  printf $input "ecut % .16g Ry\n", $generalRef->{'general'}->{'ecut'};
  printf $input "nstep %i\n", $generalRef->{'general'}->{'nstep'};
  if( $generalRef->{'epsilon'}->{'method'} eq "input" ) {
    printf $input "diemac % .16g\n", $generalRef->{'structure'}->{'epsilon'};
  } else {  # If no input dielectric then need to guess
    print $input "diemac 10\n";
  }

  printf $input "occopt %i\ntsmear % .16g Ry\n", 
      $generalRef->{'general'}->{'occopt'}, $generalRef->{'general'}->{'degauss'};

  print $input "npfft 1\n";
  printf $input "charge % .16g\n", $generalRef->{'structure'}->{'charge'};
  
  printf $input "nsppol %i\n", $generalRef->{'general'}->{'nspin'};
  if( $generalRef->{'general'}->{'nspin'} != 1 ) {
    printf $input "spinat\n%s", $generalRef->{'general'}->{'smag'};
  }

  if( $calcFlag == 1 ) {
    printf $input "fband % .16g\n", $generalRef->{'general'}->{'fband'};
    print $input "prtden 1\nprtpot 1\nkptopt 1\n";
#    printf $input "ngkpt %i %i %i\nnshiftk 1\nshiftk", $specificRef->{'kmesh'}[0], 
#          $specificRef->{'kmesh'}[1], $specificRef->{'kmesh'}[2];
#    for( my $i = 0; $i < 3; $i++ ) {
#      if( $specificRef->{'kshift'}[$i]==1 ) {
#        print $input " 0.5";  
#      } else {
#        print $input " 0.0";
#      }
#    }
#    print $input "\n";

    printf $input "toldfe  % .16g\n", $specificRef->{'toldfe'};
    if( $generalRef->{'calc_stress'} ) {
      print $input "optstress 1\n";
    } else {
      print $input "optstress 0\n";
    }

    if( $generalRef->{'calc_force'} ) {
      print $input "optforces 1\n";
    } else {
      print $input "optforces 0\n";
    }

  } else {
    # NOTE difference in definition between QE and ABINIT for tolerance leads to square here
    printf $input "nband %i\ntolwfr % .16g\n", $specificRef->{'nbands'}+$generalRef->{'general'}->{'abpad'}, 
                                               $specificRef->{'toldfe'};
    printf $input "iscf -2\ngetden -1\nistwfk *1\nnbdbuf %i\n", $generalRef->{'general'}->{'abpad'};
    printf $input "wfoptalg 114\n";
  }
  if( defined $generalRef->{'general'}->{'verbatim'}->{'abinit'} ) {
    my $s = $generalRef->{'general'}->{'verbatim'}->{'abinit'};
    $s =~ s/;/\n/g;
#    print $input $generalRef->{'general'}->{'verbatim'}->{'abinit'} . "\n";
    print $input $s . "\n";
  }

  print $input $kptString;
}

sub kptGenAbi
{
  my ( $hashRef, $split ) = @_;
      
  my $kptText = "";
  my @q = ( 0, 0, 0 ); 
  if( $hashRef->{'nonzero_q'} && not $split ) {
    $q[0] = $hashRef->{'photon_q'}[0];
    $q[1] = $hashRef->{'photon_q'}[1];
    $q[2] = $hashRef->{'photon_q'}[2];
  } 
  

  for( my $x = 0; $x < $hashRef->{'kmesh'}[0]; $x++ ) {
    my $xk = $hashRef->{'kshift'}[0]/$hashRef->{'kmesh'}[0] + $x/$hashRef->{'kmesh'}[0] - $q[0];
    while( $xk > 1 ) { $xk -= 1.0; }
    while( $xk < -1 ) { $xk += 1.0; }
    for( my $y = 0; $y < $hashRef->{'kmesh'}[1]; $y++ ) {
      my $yk = $hashRef->{'kshift'}[1]/$hashRef->{'kmesh'}[1] + $y/$hashRef->{'kmesh'}[1] - $q[1];
      while( $yk > 1 ) { $yk -= 1.0; }
      while( $yk < -1 ) { $yk += 1.0; }
      for( my $z = 0; $z < $hashRef->{'kmesh'}[2]; $z++ ) {
        my $zk = $hashRef->{'kshift'}[2]/$hashRef->{'kmesh'}[2] + $z/$hashRef->{'kmesh'}[2] - $q[2];
        while( $zk > 1 ) { $zk -= 1.0; }
        while( $zk < -1 ) { $zk += 1.0; }

        $kptText .= sprintf "%19.15f  %19.15f  %19.15f\n", $xk, $yk, $zk;
      }
    }
  }


  return $kptText;
}
#sub QEparseEnergies {}
#sub QEparseEnergies62  {}


1;
