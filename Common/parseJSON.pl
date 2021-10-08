#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
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


my $input_filename = $ARGV[0];
my $config_filename = $ARGV[1];
my $type_filename = $ARGV[2];

unless(  -e $input_filename )
{
  die "Could not find file $input_filename!\n";
}

unless( -e $config_filename )
{
  die "Could not find config file $config_filename!\n";
}

my $json = JSON::PP->new;
my $config;
if( open( my $in, "<", $config_filename ))
{
  local $/ = undef;
  $config = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open config file $config_filename\n$!";
}

my $typeDef;
if( open( my $in, "<", $type_filename ) )
{ 
  local $/ = undef;
  $typeDef = $json->decode(<$in>);
  close($in);
}
else
{ 
  die "Failed to open $type_filename\n$!";
}




my %decoder = (
  'calc' => 'calc.mode',
  'dft' => 'dft.program',
  'nkpt' => 'bse.kmesh',
  'ngkpt' => 'dft.den.kmesh',
  'ngkpt.auto' => 'dft.den.auto',
  'photon_q' => 'calc.photon_q',
  'dft.split' => 'dft.bse.split',
  'dft.qe_redirect' => 'dft.redirect',
  'nbands' => 'bse.nbands',
  'dft_energy_range' => 'bse.dft_energy_range',
  'obf_energy_range' => 'nope.obf_energy_range',
  'obkpt' => 'nope.obkpt',
  'obf.nbands' => 'nope.obf_nbands',
  'trace_tol' => 'nope.trace_tol',
  'acc_level' => 'nope.acc_level',
  'k0' => 'bse.kshift',
  'fband' => 'dft.fband',
  'occopt' => 'dft.occopt',
  'mixing' => 'dft.mixing',
  'acell' => 'structure.rscale',
  'rprim' => 'structure.rprim',
  'ntypat' => 'nope.ntype',
  'typat' => 'structure.typat',
  'znucl' => 'structure.znucl',
  'zsymb' => 'structure.zsymb',
  'pp_list' => 'psp.pp_list',
  'pp_database' => 'psp.pp_database',
  'ecut.qualtiy' => 'psp.ecut_quality',
  'natom' => 'nope.natoms',
  'coord' => 'nope.structure.coord',
  'xred' => 'structure.xred',
  'ecut' => 'dft.ecut',
  'diemac' => 'structure.epsilon',
  'toldfe' => 'dft.toldfe',
  'tolwfr' => 'dft.tolwfr',
  'nstep' => 'dft.nstep',
  'dft.startingwfc' => 'dft.startingwfc',
  'dft.diagonalization' => 'dft.diagonalization',
  'dft.ndiag' => 'dft.ndiag',
  'dft.functional' => 'dft.functional',
  'dft.exx.qmesh' => 'dft.exx.qmesh',
  'dft.nscf.poolsize' => 'dft.bse.poolsize',
  'verbatim' => 'dft.verbatim',
  'para_prefix' => 'computer.para_prefix',
  'ser_prefix' => 'computer.ser_prefix',
  'abpad' => 'dft.abpad',
  'scfac' => 'bse.core.scfac',
  'screen.shells' => 'screen.shells',
  'opf.hfkgrid' => 'opf.shirley.hfkgrid',
  'opf.fill' => 'opf.shirley.fill',
  'opf.opts' => 'opf.shirley.opts',
  'opf.program' => 'opf.program',
  'screen.nkpt' => 'screen.kmesh',
  'screen.k0' => 'screen.kshift',
  'screen.nbands' => 'screen.nbands',
  'caution' => 'opf.shirley.caution',
  'nedges' => 'nope.nedges',
  'edges' => 'calc.edges',
  'cnbse.nbuse' => 'nope.nbuse',
  'cnbse.xmesh' => 'bse.xmesh',
  'cnbse.rad' => 'bse.core.screen_radius',
  'metal' => 'structure.metal',
  'cksshift' => 'nope.cksshift',
  'cksstretch' => 'nope.cksstretch',
  'cnbse.niter' => 'bse.core.haydock.niter',
  'haydock_convergence' => 'bse.core.haydock.converge.thresh',
  'cnbse.spect_range' => 'bse.core.plot.range',
  'cnbse.broaden' => 'bse.core.broaden',
  'cnbse.strength' => 'bse.core.strength',
  'cnbse.solver' => 'bse.core.solver',
  'cnbse.gmres.elist' => 'bse.core.gmres.elist',
  'cnbse.gmres.erange' => 'bse.core.gmres.erange',
  'cnbse.gmres.nloop' => 'bse.core.gmres.nloop',
  'cnbse.gmres.gprc' => 'bse.core.gmres.gprc',
  'cnbse.gmres.ffff' => 'bse.core.gmres.ffff',
  'cnbse.write_rhs' => 'bse.core.write_rhs',
  'cnbse.gw.control' => 'bse.core.gw.control',
  'bse.gw.cstr' => 'nope.bse_gw_cstr',
  'bse.gw.vstr' => 'nope.bse_gw_vstr',
  'bse.gw.gap' => 'nope.gwgap',
  'degauss' => 'dft.degauss',
  'ibrav' => 'nope.ibrav',
  'isolated' => 'nope.isolated',
  'noncolin' => 'dft.noncolin',
  'prefix' => 'nope.prefix',
  'ppdir' => 'psp.ppdir',
  'dft.calc_stress' => 'dft.calc_stress',
  'dft.calc_force' => 'dft.calc_force',
  'spinorb' => 'dft.spinorb',
  'work_dir' => 'nope.wordir',
  'tmp_dir' => 'dft.tmp_dir',
  'den.kshift' => 'dft.den.kshift',
  'core_offset' => 'screen.core_offset.enable',
  'ham_kpoints' => 'nope.ham_kpoints',
  'nbse.niter' => 'bse.val.haydock.niter',
  'nbse.backf' => 'bse.val.backf',
  'nbse.aldaf' => 'bse.val.aldaf',
  'nbse.qpflg' => 'bse.val.qpflg',
  'nbse.bwflg' => 'bse.val.bwflg',
  'nbse.bande' => 'bse.val.bande',
  'nbse.bflag' => 'bse.val.bflag',
  'nbse.lflag' => 'bse.val.lflag',
  'nbse.convergence' => 'nope.convergence',
  'nbse.decut' => 'bse.val.decut',
  'nbse.se_rs' => 'nope.se.rs',
  'nbse.se_metal' => 'nope.se.metal',
  'nbse.se_niter' => 'nope.se.niter',
  'nbse.spect_range' => 'bse.val.plot.range',
  'tot_charge' => 'dft.tot_charge',
  'nspin' => 'dft.nspin',
  'smag' => 'dft.smag',
  'ldau' => 'dft.ldau.Hubbard_U',
  'qe_scissor' => 'nope.qe_scissor',
  'nphoton' => 'nope.nphoton',
  'ser_bse' => 'nope.ser_bse',
  'spin_orbit' => 'bse.core.spin_orbit',
  'screen_energy_range' => 'screen.dft_energy_range',
  'screen.grid.scheme' => 'screen.grid.scheme',
  'screen.grid.rmode' => 'screen.grid.rmode',
  'screen.grid.ninter' => 'nope.screen_grid_ninter',
  'screen.grid.shells' => 'screen.grid.shells',
  'screen.grid.xyz' => 'nope.screen_grid_xyz',
  'screen.grid.rmax' => 'screen.grid.rmax',
  'screen.grid.nr' => 'nope.screen_grid_nr',
  'screen.grid.ang' => 'screen.grid.ang',
  'screen.grid.deltar' => 'screen.grid.deltar',
  'screen.lmax' => 'screen.grid.lmax',
  'screen.grid.nb' => 'nope.screen_grid_nb',
  'screen.final.rmax' => 'screen.final.rmax',
  'screen.final.dr' => 'screen.final.dr',
  'screen.model.dq' => 'screen.model.SLL.dq',
  'screen.model.qmax' => 'screen.model.SLL.qmax',
  'screen.legacy' => 'nope.screen_legacy',
  'screen.augment' => 'screen.augment',
  'screen.wvfn' => 'nope.screen_wvfn',
  'screen.convertstyle' => 'screen.convertstyle',
  'screen.inversionstyle' => 'screen.inversionstyle',
  'screen.mode' => 'screen.mode',
  'bse.wvfn' => 'nope.bse_wvfn',
  'hamnum' => 'nope.hamnum',
  'echamp' => 'bse.core.gmres.echamp',
  'bshift' => 'nope.bshift' );


open IN, "<", $input_filename or die "Failed to open input file $input_filename\n$!\n";

my $rawInputFile = '';
my $inputString = '';
while( my $line = <IN> )
{
  $rawInputFile .= $line;
  # if there are comment characters -- #, *, or ! --
  #   remove them and everything following
  $line =~ s/[#\*\\!].*/ /;
  # just pad with spaces, not that inefficient
  $line =~ s/\{/ \{ /;
  $line =~ s/\}/ \} /;
  $inputString .= $line;
}
close IN;


my @inputFile = split ' ', $inputString;

if( 0 ){
foreach my $i (@inputFile)
{
  print "$i\n";
}
}

my %inputHash;
my $i = 0;
my $tag = 1;
my $curly = 0;
my $key = '';
my $val = '';
my $errorBuffer = '';
while( $i < scalar @inputFile )
{
  if( $tag == 1 )
  {
    $errorBuffer .= $inputFile[$i] . "\n";
    die "Misplaced braces when expecting a tag\n>>>>\n$errorBuffer<<<<<\n" if( $inputFile[$i] =~ m/\{|\}/ );
#    print "$inputFile[$i] >>>> ";
    $tag = 0;
    $key = $inputFile[$i];
    $val = '';
    $errorBuffer = '';
    $errorBuffer .= $inputFile[$i-1] . "\n" if( $i > 0 );
    $errorBuffer .= $inputFile[$i] . "\n";
  }
  elsif ( $curly == 0 ) 
  {
    $errorBuffer .= $inputFile[$i] . "\n";
    if( $inputFile[$i] =~ m/\{/ )
    {
      $curly = 1;
    }
    elsif( $inputFile[$i] =~ m/\}/ )
    {
      die "Close brace when not expected\n>>>>\n$errorBuffer<<<<<\n";
    }
    else
    {
#      print "$inputFile[$i]\n";
      $tag = 1;
      $val .= $inputFile[$i] . " ";
    }
  } else
  {
    $errorBuffer .= $inputFile[$i] . "\n";
    die "Second open {\n>>>>\n$errorBuffer<<<<<\n" if( $inputFile[$i] =~ m/\{/ );
    if( $inputFile[$i] =~ m/\}/ )
    {
      $curly = 0;
      $tag = 1;
#      print "\n";
    }
    else
    {
      $val .= $inputFile[$i] . " ";
#      print "$inputFile[$i] ";
    }
  }
  $i++;
  # there will always be a trailing space
#  chop( $val );
  $inputHash{ $key } = $val;
}

my $haveLegacy = 0;
INPUT: foreach my $key ( keys %inputHash ) 
{
  my @newKey = split /\./, $key;
  my $ref = $config;
  for( my $i = 0; $i < scalar @newKey; $i++ )
  {
    if( exists $ref->{$newKey[$i]} )
    {
      $ref = $ref->{$newKey[$i]};
    }
    else
    {
      $haveLegacy = 1;
      print "Unrecognized input flag: $key\n  Attempting legacy conversion\n";
      last INPUT;
    }
  }
}

if( $haveLegacy == 1 )
{
  foreach my $key ( keys %inputHash )
  {
    die "Unrecognized input flag: $key\n No recovery possible!" unless( exists $decoder{$key} );
    my $newKey = $decoder{ $key };
    print "$key : $newKey  $inputHash{ $key }\n";
    $inputHash{ $newKey } = $inputHash{ $key };
    delete( $inputHash{ $key } ) unless( $key eq $newKey );
    $rawInputFile =~ s/$key/$newKey/;
    my @newKey = split /\./, $newKey;
    my $ref = $config;
    next if( $newKey[0] eq 'nope' );
    for( my $i = 0; $i < scalar @newKey; $i++ )
    {
      if( exists $ref->{$newKey[$i]} )
      {
        $ref = $ref->{$newKey[$i]};
      }
      else
      {
        die "Unrecognized input flag: $key\n  Legacy conversion failed\n";
      }
    }
  }

  my $newInputFile = $input_filename . ".mod3";
  open OUT, ">", $newInputFile or die "$!";
  print OUT $rawInputFile;
  close OUT;
}

print "Storing parsed data\n\n";
# If we made it here all the keys are valid
foreach my $key ( keys %inputHash )
{
  my $value = $inputHash{ $key };
  print "$key $value\n";
  my @newKey = split /\./, $key;
  next if( $newKey[0] eq 'nope' );


  my $type = $typeDef;
  for( my $i = 0; $i < scalar @newKey; $i++ )
  {
    $type = $type->{$newKey[$i]} ;
  } 

  my $hashref = $config;
  for( my $i = 0; $i < scalar @newKey - 1; $i++ )
  {
    $hashref = $hashref->{$newKey[$i]};
  }

  my $regex;
  $regex = '^(-?\d+)' if( $type =~ m/i/ );
  $regex = '^(-?\d*\.?\d+([eEdD][+-]?\d+)?)' if( $type =~ m/f/ );
  $regex = '^([\w\S\s]+)$' if ( $type =~ m/s|S/ );


  # if array
  if( $type =~ m/a/ )
  {
    my @rawArray = split ' ', $value;
    foreach my $i (@rawArray)
    {
      die "Failed to match: $i of type $type\n" unless( $i =~ m/$regex/ );
    }
    if( $type =~ m/[if]/ )
    {
      for( my $i = 0; $i < scalar @rawArray; $i++ )
      { $rawArray[$i] *= 1 }
    }
    elsif( $type =~ m/s/ )
    {
      for( my $i = 0; $i < scalar @rawArray; $i++ )
      { $rawArray[$i] = lc $rawArray[$i] }
    }

    # If we did legacy translation, patch up the incompatibilities
    if( $haveLegacy == 1 )
    {
      if( $key =~ m/bse.(core|val).plot.range/ )
      {
        my $points = shift @rawArray;
        $config->{'bse'}->{'core'}->{'plot'}->{'points'} = $points if( $key =~ m/core/ );
        $config->{'bse'}->{'val'}->{'plot'}->{'points'} = $points if( $key =~ m/val/ );
      }
      if( $key =~m/pp_list/ )
      {
        $config->{'psp'}->{'source'} = 'manual' if( scalar @rawArray > 0 );
      }
    }
    # End legacy fix
    $hashref->{$newKey[-1]} = [@rawArray];
  }
  else
  {
    unless( $value eq ' ' )
    {
      $value =~ s/^\s+//;
      $value =~ s/\s+$//;
    }
    # If we did legacy translation, patch up the incompatibilities
    if( $key =~ m/epsilon/ )
    {
      if( $value =~ m/dfpt/ )
      {
         $value = 0;
         $config->{'dft'}->{'epsilon'}->{'method'} = 'dfpt'
      }
      else
      {
         $config->{'dft'}->{'epsilon'}->{'method'} = 'input'
      }
    }
    elsif( $key =~ m/core_offset/ )
    {
      if( $value =~ m/\d/ )
      {
        $config->{'screen'}->{'core_offset'}->{'energy'} = $value;
        $value = 'true';
      }
    }
    # end fix  
    if( $type =~ m/b/ )
    {
      if( $value =~ m/t|true|1/i )
      {
        $value = $JSON::PP::true
      }
      elsif( $value =~ m/f|false|0/i )
      {
        $value = $JSON::PP::false
      }
      else
      {
        die "Failed to match: $value of type $type\n";
      }
    }
    else
    {
      if( $value =~ m/$regex/ )
      {
        $value = $1;
        $value *= 1 if( $type =~ m/[if]/ );
        $value = lc $value if( $type =~ m/s/ );
      }
      else
      {
        die "Failed to match: $value of type $type\n";
      }
    }
    $hashref->{$newKey[-1]} = $value;
  }
}

my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);
open OUT, ">", "parsedInputFile" or die;
print OUT $json->encode($config);
close OUT;


copy( "parsedInputFile", "oceanDatafile") ;


##### REMOVE IN FUTURE
open OUT, ">", "dft" or die;
print OUT $config->{'dft'}->{'program'} . "\n";
close OUT;

open OUT, ">", "calc" or die;
print OUT $config->{'calc'}->{'mode'} . "\n";
close OUT;
