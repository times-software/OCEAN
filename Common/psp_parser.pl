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
  foreach my $i (@typat)
  {
    print $i . "\n";
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
      copy( $ppdir . '/' . $p . $suffix, "psp");
      copy( $ppdir . '/' . $p, "psp") if ( $suffix ne '' && -e $ppdir . '/' . $p, "psp" );
      my $q = $p;
      $q =~ s/\.psp8$//i;
      $q =~ s/\.upf$//i;
      $q .= '.in';
      copy( $ppdir . '/' . $q, "psp") if( -e  $ppdir . '/' . $q );
    }
  }
  my @electrons;
  my @ppHash;
  my $ppdir = './psp/';
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

    $nelectron *= 1;
    print "$ppfile has $nelectron electrons\n";
    push @electrons, $nelectron;

    close $pp;
    if( open $pp, $ppdir . $ppfile  )
    {
      local $/ = undef;
      push @ppHash, md5_hex(<$pp>);
    }
  }

  my $totalElectrons = 0;
  foreach my $typat (@typat)
  {
    if( $typat > scalar @electrons )
    {  die "Found entry in typat that is higher than number of pseudopotentials\n";  }
    $totalElectrons += $electrons[$typat-1];
  }

  print "Total electrons in system: $totalElectrons\n";


  $oceanData->{'structure'}->{'valence_electrons'} = $totalElectrons; 
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
  while( my $line = <$pp> )
  {
#    if( $line =~ m/z_valence=\"\s*(\d+.?\d*)\s*\"/i  || 
#        $line =~ m/z\s+valence=\"\s*(\d+.?\d*)\s*\"/i )
    if( $line =~ m/z_valence/i || $line =~ m/z valence/i )
    {
      
      if( $line =~ m/(\d+.?\d*[Ee]?\+?\d*)/ )
      {
#      print $line;
        $n = $1;
      }
      else
      {
        print "FAILED:   $line\n";
        $n = -2;
      }
      last;
    }
#    elsif( $line =~ m/z_valence/i )
#    {
#      print $line;
#      $n = -2;
#      last;
#    }
  }

  return $n;
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
