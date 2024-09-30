#!/usr/bin/perl
# Copyright (C) 2019-2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
require JSON::PP;
JSON::PP->import;
use File::Copy;
use Cwd qw(cwd);
use Digest::MD5 qw(md5_hex);

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/psp_parser\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = cwd . "/../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }

my $dataFile = "oceanDatafile";
if( -e $dataFile )
{
  my $json = JSON::PP->new;
  my $oceanData;
  if( open( my $in, "<", $dataFile ))
  {
    local $/ = undef;
    $oceanData = $json->decode(<$in>);
    close($in);
  }
  else
  {
    die "Failed to open config file $dataFile\n$!";
  }
  my @typat = @{ $oceanData->{'structure'}->{'typat'} };
  my @countByTypat;
  foreach my $i (@typat)
  {
#    print $i . "\n";
    $countByTypat[$i]++;
  }
  for( my $i = 1; $i < scalar @countByTypat; $i++ )
  {
    printf "%i  : %i\n", $i, $countByTypat[$i];
  }
  my @psp = @{ $oceanData->{'psp'}->{'pp_list'} };
  my $suffix = '';
  my $dft = $oceanData->{'dft'}->{'program'};
  if( $dft =~ m/qe/ )
  {
    $suffix = '.UPF';
  }
  if( $oceanData->{'psp'}->{'source'} eq 'manual' )
  {
    my $ppdir = $ENV{"OCEAN_WORKDIR"} . '/' . $oceanData->{'psp'}->{'ppdir'};
    print "$ppdir \n";
    mkdir "psp";
    foreach my $p (@psp)
    {
      copy( $ppdir . '/' . $p . $suffix, "psp") if( -e $ppdir . '/' . $p . $suffix );
      copy( $ppdir . '/' . $p, "psp") if ( $suffix ne '' && -e $ppdir . '/' . $p, "psp" );
      my $q = $p;
      $q =~ s/\.psp8$//i;
      $q =~ s/\.upf$//i;
      $q .= '.in';
      copy( $ppdir . '/' . $q, "psp") if( -e  $ppdir . '/' . $q );
    }
  }
  my @electrons;
  my @semicore;
  my @ppHash;
  my $ppdir = './psp/';
  foreach my $ppfile (@psp)
  {
    $ppfile .= $suffix unless( $ppfile =~ m/$suffix$/i );
    my $pp;
  #  open $pp, "../$ppfile" or die "Failed to open ../$ppfile\n$!";
    open $pp, $ppdir . $ppfile or die "Failed to open $ppdir" . "$ppfile\n$!";

    my $nelectron = -1;
    my $semicore = -1;
    if( $suffix eq '' )
    {
      $nelectron = abinitParser( $pp );
    }
    elsif( $suffix eq '.UPF' )
    {
      ($nelectron, $semicore )= upfParser( $pp );
    }
    else
    {
      die "Un-supported psp\n";
    }
    die "Failed to parse $ppfile: $nelectron\n" if( $nelectron < 1 );

    $nelectron *= 1;
    print "$ppfile has $nelectron electrons\n";
    push @electrons, $nelectron;
    print "$ppfile has $semicore semicore electrons\n";
    push @semicore, $semicore;

    close $pp;
    if( open $pp, $ppdir . $ppfile  )
    {
      local $/ = undef;
      push @ppHash, md5_hex(<$pp>);
    }
  }

  my $totalElectrons = 0;
  my $totalSemicore = 0;
  foreach my $typat (@typat)
  {
    if( $typat > scalar @electrons )
    {  die "Found entry in typat that is higher than number of pseudopotentials\n";  }
    $totalElectrons += $electrons[$typat-1];
    $totalSemicore += $semicore[$typat-1] if( $semicore[$typat-1] >= 0 );
  }

  print "Total electrons in system: $totalElectrons\n";
  print "Total semi-core electrons in system: $totalSemicore\n";


  $oceanData->{'structure'}->{'valence_electrons'} = $totalElectrons; 
  $oceanData->{'structure'}->{'semicore_electrons'} = $totalSemicore; 
  # Adjust by charge (electrons are negative)
  $oceanData->{'structure'}->{'valence_electrons'} -= $oceanData->{'structure'}->{'charge'};

  $oceanData->{'psp'}->{'ppdir'} = cwd . '/psp/';

  $oceanData->{'psp'}->{'pphash'} = [@ppHash];
  
  my $enable = 1;
  $json->canonical([$enable]);
  $json->pretty([$enable]);
  open OUT, ">", $dataFile or die;
  print OUT $json->encode($oceanData);
  close OUT;

}
else
{
  open OUT, ">", "nelectrons" or die "Failed to open nelectrons for writing\n$!";
  print OUT "-1\n";
  close OUT;

  my @typat;
  open TYPAT, "typat" or die "Failed to open typat\n$!";
  while (<TYPAT>)
  {
    chomp;
    push @typat, split ' ';
  }
  close TYPAT;

  my @psp;
  open PSP, "pplist" or die "Failed top open pplist\n$!";
  while (<PSP>)
  {
    chomp;
    push @psp, split ' ';
  }
  close PSP;

  #This is to deal with the psp name difference between abinit and QE
  my $suffix = '';
  open DFT, "dft" or die "Failed to open file dft\n$!";
  my $dft = <DFT>;
  close DFT;
  if( $dft =~ m/qe/i )
  {
    $dft = 'qe';
    $suffix = '.UPF';
  }

  open PPDIR, "ppdir" or die "Failed to open file ppdir\n$!";
  my $ppdir = <PPDIR>;
  close PPDIR;
  chomp $ppdir;
  $ppdir =~ s/'//g;
  $ppdir .= '/';

  my @electrons;
  foreach my $ppfile (@psp)
  {
    $ppfile .= $suffix unless( $ppfile =~ m/$suffix$/ );
    my $pp;
  #  open $pp, "../$ppfile" or die "Failed to open ../$ppfile\n$!";
    open $pp, $ppdir . $ppfile or die "Failed to open $ppdir" . "$ppfile\n$!";

    my $nelectron = -1;
    if( $suffix eq '' )
    {
      $nelectron = abinitParser( $pp );
    }
    elsif( $suffix eq '.UPF' )
    {
      $nelectron = upfParser( $pp );
    }
    else
    {
      die "Un-supported psp\n";
    }
    die "Failed to parse $ppfile: $nelectron\n" if( $nelectron < 1 );

    print "$ppfile has $nelectron electrons\n";
    push @electrons, $nelectron;

    close $pp;
  }

  my $totalElectrons = 0;
  foreach my $typat (@typat)
  {
    if( $typat > scalar @electrons ) 
    {  die "Found entry in typat that is higher than number of pseudopotentials\n";  }
    $totalElectrons += $electrons[$typat-1];
  }

  print "Total electrons in system: $totalElectrons\n";


  open OUT, ">", "nelectrons" or die "Failed to open nelectrons for writing\n$!";
  print OUT "$totalElectrons\n";
  close OUT;
}
exit;

sub upfParser
{
  my $pp = shift;
  my $n = -1;
  my $element = "";
  while( my $line = <$pp> )
  {
    if( $line =~ m/z_valence/i || $line =~ m/z valence/i )
    {
      
      if( $line =~ m/(\d+.?\d*[Ee]?\+?\d*)/ )
      {
        $n = $1;
      }
      else
      {
        print "FAILED:   $line\n";
        $n = -2;
      }
      last;
    } 
    elsif( $line =~ m/element\s*="\s*(\w+)/i ) {
      $element = $1;
    }
#    last if( $n != -1 && $element ne "" );
  }

  my $semicore = -1;
  if( $element ne "" ) {
    my $z = ocean_el2z( $element );
    if( $z > 0 ) {
      my $core = 0;
      $core += 2 if( $z >= 3 );
      $core += 2 if ($z >= 11 );
      $core += 6 if( $z >= 12 ); # Na 2p is too shallow, overlaps with say O 2s
      $core += 2 if ($z >= 19 );
      $core += 6 if( $z >= 21 ); # K, Ca 3p is too shallow
      $core += 10 if( $z >= 33 );  # the 3d is tricky, but As should be ~40eV
      $core += 2 if ($z >= 37 );
      $core += 6 if( $z >= 39 );  # Rb/Sr 4p too shallow
      $core += 10 if( $z >= 52 );  # the 4d is tricky, but Te should be ~40eV
      $core += 2 if( $z >= 55 );
      $core += 6 if( $z >= 72  ); # 5p too shallow for most 4f elements
      $core += 14 if ( $z >= 74 ); # 4f hopefully semi core by W
      # don't know about 6d
      $core += 18 if ( $z >= 87 );
      $semicore = $n - ($z - $core );
    }
  }

  return ($n, $semicore);
}

sub abinitParser
{
  my $pp = shift;
  my $n = -1;
  <$pp>;
  my $line = <$pp>;
  $line =~ m/^\s*(\d+\.?\d*)\s+(\d+\.?\d*)/ or die "Failed to parse\n$line";
  $n = $2;
  return $n;
}


sub ocean_el2z 
{
  my ($el) = @_;

  my @chars = split("", $el);
  
  my $normEl;
  $chars[0] = uc( $chars[0] );
  if( $chars[1] =~ m/\w/ ) {
    $normEl = $chars[0] . lc( $chars[1] );
  } else {
    $normEl = $chars[0];
  }

  my %el = ( "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8,
             "F" => 9, "Ne" => 9, "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15, "S" => 16, 
             "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20, "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, 
             "Mn" => 25, "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30, "Ga" => 31, "Ge" => 32, 
             "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36, "Rb" => 37, "Sr" => 38, "Y" => 39, "Zr" => 40, 
             "Nb" => 41, "Mo" => 42, "Tc" => 43, "Ru" => 44, "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48,
             "In" => 49, "Sn" => 50, "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54, "Cs" => 55, "Ba" => 56, 
             "La" => 57, "Ce" => 58, "Pr" => 59, "Nd" => 60, "Pm" => 61, "Sm" => 62, "Eu" => 63, "Gd" => 64,
             "Tb" => 65, "Dy" => 66, "Ho" => 67, "Er" => 68, "Tm" => 69, "Yb" => 70, "Lu" => 71, "Hf" => 72,
             "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77, "Pt" => 78, "Au" => 79, "Hg" => 80,
             "Tl" => 81, "Pb" => 82, "Bi" => 83, "Po" => 84, "At" => 85, "Rn" => 86, "Fr" => 87, "Ra" => 88,
             "Ac" => 89, "Th" => 90, "Pa" => 91, "U" => 92, "Np" => 93, "Pu" => 94, "Am" => 95, "Cm" => 96,
             "Bk" => 97, "Cf" => 98, "Es" => 99, "Fm" => 100, "Md" => 101, "No" => 102, "Lr" => 103 ); 

  my $z = -1;
  $z = $el{$normEl}  if( exists $el{$normEl} );

  return $z;
}
