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
use Scalar::Util qw( looks_like_number );


my $input_filename = $ARGV[0];
my $config_filename = $ARGV[1];
my $type_filename = $ARGV[2];
my $enum_filename = $ARGV[3];

unless(  -e $input_filename )
{
  die "Could not find file $input_filename!\n";
}

unless( -e $config_filename )
{
  die "Could not find config file $config_filename!\n";
}

unless( defined $enum_filename && -e $enum_filename )
{
  die "Could not find enum file " . ( defined $enum_filename ? $enum_filename : '<undef>' ) . "!\n";
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

my $enumDef = loadEnumDefinition( $enum_filename, $json );

validateSchemaTypeMatch( $config, $typeDef );
validateCaseInsensitivePaths( $config, '' );
validateTypeSpecs( $typeDef, '' );
validateEnumSpecs( $enumDef, $config, $typeDef, '' );




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
  'control' => 'nope.control',
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
  'ecut.quality' => 'psp.ecut_quality',
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
  'core_offset' => 'cls.enable',
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
my @inputOrder;
my %seenRawKey;
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
    if( exists $seenRawKey{ $key } )
    {
      die "Duplicate input flag: $key\n";
    }
    $seenRawKey{ $key } = 1;
    push @inputOrder, $key;
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
my $haveLegacyCorePlot = 0;
my $haveLegacyValPlot = 0;

INPUT: foreach my $key ( @inputOrder ) 
{
  unless( defined canonicalInputKey( $config, $key ) )
  { 
    $haveLegacy = 1;
    print "Unrecognized input flag: $key\n  Attempting legacy conversion\n";
    last INPUT;
  }
}

my %seenInputKey;
my %inputSource;
my %resolvedInputHash;

if( $haveLegacy == 1 )
{
  foreach my $key ( @inputOrder )
  {
    if( lc($key) eq 'ppdir' ) {
      if( $inputHash{ $key } =~ m/^\s*'(.+)'\s*$/ ) {
        print $inputHash{ $key };
        $inputHash{ $key } = $1;
        print "  " .$inputHash{ $key } . "\n";
      }
      if( $inputHash{ $key } =~ m/^\.\.\/$/ ) {
        print $inputHash{ $key };
        $inputHash{ $key } = './';
        print "  " .$inputHash{ $key } . "\n";
      }
    }
  }
  foreach my $key ( @inputOrder )
  {
    my $canonicalKey = canonicalInputKey( $config, $key );
    if( defined $canonicalKey )
    {
      if( exists $seenInputKey{ $canonicalKey } )
      {
        die "Duplicate input after legacy conversion: "
          . $seenInputKey{ $canonicalKey } . " and $key both set $canonicalKey\n";
      }
      $seenInputKey{ $canonicalKey } = $key;
      $inputSource{ $canonicalKey } = { raw => $key, legacy => 0 };
      $resolvedInputHash{ $canonicalKey } = $inputHash{ $key };
      print "Comment: Mixed new and legacy input:  $key\n";
      next;
    }

    my $lckey = lc($key);
#    $key = lc($key) unless( exists $decoder{$key} );
#    die "Unrecognized input flag: $key\n No recovery possible!" unless( exists $decoder{$lckey} );
    if( exists $decoder{$lckey} ) {
      my $newKey = $decoder{ $lckey };
      if( exists $seenInputKey{ $newKey } )
      {
        die "Duplicate input after legacy conversion: "
          . $seenInputKey{ $newKey } . " and $key both set $newKey\n";
      }
      $seenInputKey{ $newKey } = $key;
      $inputSource{ $newKey } = { raw => $key, legacy => ( $key ne $newKey ) };
      $resolvedInputHash{ $newKey } = $inputHash{ $key };
      print "$key : $newKey  $inputHash{ $key }\n";
      $rawInputFile =~ s/$key/$newKey/;
      my @newKey = split /\./, $newKey;
      my $ref = $config;
      if( $newKey[0] eq 'nope' )
      {
        print "COMMENT: Ignoring recognized legacy input flag: $key ($newKey)\n";
        next;
      }
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
      if( $lckey eq 'cnbse.spect_range' ) {
        $haveLegacyCorePlot = 1;
      } elsif ( $lckey eq 'nbse.spect_range' ) {
        $haveLegacyValPlot = 1;
      }
    } else {
      die "Unsupported input flag. Neither new nor legacy:  $key\n";
    }
  }

  my $newInputFile = $input_filename . ".mod3";
  open OUT, ">", $newInputFile or die "Failed to open $newInputFile \n$!";
  print OUT $rawInputFile;
  close OUT;
}
else
{
  foreach my $key ( @inputOrder )
  {
    my $canonicalKey = canonicalInputKey( $config, $key );
    if( exists $seenInputKey{ $canonicalKey } )
    {
      die "Duplicate input after case normalization: "
        . $seenInputKey{ $canonicalKey } . " and $key both set $canonicalKey\n";
    }
    $seenInputKey{ $canonicalKey } = $key;
    $inputSource{ $canonicalKey } = { raw => $key, legacy => 0 };
    $resolvedInputHash{ $canonicalKey } = $inputHash{ $key };
  }
}

%inputHash = %resolvedInputHash;

print "Storing parsed data\n\n";
# If we made it here all the keys are valid
my %suppliedInputKey;
my %suppliedInputSource;
foreach my $key ( keys %inputHash )
{
  my $value = $inputHash{ $key };
  print "$key $value\n";
  my @newKey = split /\./, $key;
  next if( $newKey[0] eq 'nope' );
  my $sourceRef = $inputSource{ $key };
  $suppliedInputKey{ $key } = 1;
  $suppliedInputSource{ $key } = $sourceRef;


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
  my ($baseType, $constraints, $unsupported) = parseTypeSpec( $type, $key );
  if( $unsupported )
  {
    my $keyContext = formatInputKeyContext( $key, $sourceRef );
    print "WARNING: Input flag $keyContext has no functionality in the current code.\n";
  }
  $regex = '^\s*(-?\d+)\s*$' if( $baseType =~ m/i/ );
  # Full-token floating point match:
  #   ^\s* and \s*$ allow only optional leading/trailing whitespace.
  #   -? allows an optional minus sign.
  #   (?:\d+(?:\.\d*)?|\.\d+) accepts either digits with an optional decimal
  #     point and optional following digits, or a leading decimal point followed
  #     by digits. This allows 1, 1., 1.0, and .1, but rejects bare ".".
  #   (?:[eEdD][+-]?\d+)? accepts an optional Fortran/C exponent with e, E, d,
  #     or D, an optional sign, and at least one exponent digit.
  $regex = '^\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?)\s*$' if( $baseType =~ m/f/ );
  $regex = '^([\w\S\s]+)$' if ( $baseType =~ m/s|S/ );


  # if array
  if( $baseType =~ m/a/ )
  {
    my @rawArray = split ' ', $value;
    for( my $j = 0; $j < scalar @rawArray; $j++ )
    {
      dieInputTypeError( $key, $sourceRef, $rawArray[$j], $type, undef, $j + 1 )
        unless( $rawArray[$j] =~ m/$regex/ );
    }
    if( $baseType =~ m/[if]/ )
    {
      for( my $i = 0; $i < scalar @rawArray; $i++ )
      { 
        $rawArray[$i] =~ s/[dD]/e/; # won't matter for int
        $rawArray[$i] *= 1 }
    }
    elsif( $baseType =~ m/s/ )
    {
      for( my $i = 0; $i < scalar @rawArray; $i++ )
      { $rawArray[$i] = lc $rawArray[$i] }
    }

    # If we did legacy translation, patch up the incompatibilities
    if( $key =~ m/bse.core.plot.range/ )
    {
      if( $haveLegacyCorePlot ) 
      {
        my $points = shift @rawArray;
        $config->{'bse'}->{'core'}->{'plot'}->{'points'} = $points if( $key =~ m/core/ );
      }
    } elsif( $key =~ m/bse.val.plot.range/ ) 
    {
      if( $haveLegacyValPlot ) {
        my $points = shift @rawArray;
        $config->{'bse'}->{'val'}->{'plot'}->{'points'} = $points if( $key =~ m/val/ );
      }
    } elsif( $haveLegacy == 1 && $key =~m/pp_list/ )
    {
      $config->{'psp'}->{'source'} = 'manual' if( scalar @rawArray > 0 );
    }
    # End legacy fix
    validateInputEnumValue( $enumDef, $key, $sourceRef, $rawArray[$_], $_ + 1 )
      for( 0 .. $#rawArray );
    validateArrayLength( $key, $sourceRef, $type, $constraints, \@rawArray );
    $hashref->{$newKey[-1]} = [@rawArray];
  }
  else
  {
#    if( $key =~ m/tol/ ) {
#      print "STOP" . $value . "\n" . $regex . "\n";
#      $value =~ s/d/e/;
#      $value =~ m/$regex/;
#      print $1 . "\n";
#      print $1*1 . "\n";
#      exit 1;
#    }
    unless( $value eq ' ' )
    {
      $value =~ s/^\s+//;
      $value =~ s/\s+$//;
    }
    # If we did legacy translation, patch up the incompatibilities
    if( $key =~ m/structure.epsilon/ )
    {
      if( $value =~ m/dfpt/ )
      {
         $value = 0;
         $config->{'dft'}->{'epsilon'}->{'method'} = 'dfpt';
      }
      else
      {
         $config->{'dft'}->{'epsilon'}->{'method'} = 'input';
      }
    }
    elsif( $key =~ m/core_offset/ )
    {
      if( $value =~ m/\d/ )
      {
        $config->{'cls'}->{'energy'} = $value;
        $value = 'true';
      } else {
        $config->{'cls'}->{'average'} = 'true';
      }
    }
    # end fix  
    if( $baseType =~ m/b/ )
    {
      my $boolValue = lc $value;
      if( $boolValue eq 't' || $boolValue eq 'true' || $boolValue eq '.t.'
          || $boolValue eq '.true.' || $boolValue eq '1' )
      {
        $value = $JSON::PP::true
      }
      elsif( $boolValue eq 'f' || $boolValue eq 'false' || $boolValue eq '.f.'
             || $boolValue eq '.false.' || $boolValue eq '0' )
      {
        $value = $JSON::PP::false
      }
      else
      {
        dieInputTypeError( $key, $sourceRef, $value, $type,
          "accepted forms: true/false, t/f, .true./.false., 1/0" );
      }
    }
    else
    {
      if( $value =~ m/$regex/ )
      {
        $value = $1;
        $value =~ s/[dD]/e/ if( $baseType =~ m/f/ );
        $value *= 1 if( $baseType =~ m/[if]/ );
        $value = lc $value if( $baseType =~ m/s/ );
        validateInputEnumValue( $enumDef, $key, $sourceRef, $value );
      }
      else
      {
        dieInputTypeError( $key, $sourceRef, $value, $type );
      }
    }
    $hashref->{$newKey[-1]} = $value;
  }
}

validateRelationalArrayLengths( $config, $typeDef, \%suppliedInputKey, \%suppliedInputSource );

my $enable = 1;
$json->canonical([$enable]);
$json->pretty([$enable]);
open OUT, ">", "parsedInputFile" or die "Failed to open parsedInputFile\n$!";
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



# Make sure case-insensitive input matching can choose between two schema keys.
sub validateCaseInsensitivePaths
{
  my ($ref, $prefix) = @_;
  return unless( ref( $ref ) eq 'HASH' );

  my %seenKey;
  foreach my $key ( keys %$ref )
  {
    my $lcKey = lc($key);
    if( exists $seenKey{ $lcKey } )
    {
      my $path = length $prefix ? $prefix : '<root>';
      die "Case-insensitive input schema collision under $path: "
        . "$seenKey{ $lcKey } and $key\n";
    }
    $seenKey{ $lcKey } = $key;
  }

  foreach my $key ( keys %$ref )
  {
    my $newPrefix = length $prefix ? "$prefix.$key" : $key;
    validateCaseInsensitivePaths( $ref->{$key}, $newPrefix );
  }
}


# Return the exact oparse.json leaf path for an input key, ignoring user case.
sub canonicalInputKey
{
  my ($config, $key) = @_;
  my @newKey = split /\./, $key;
  my @canonicalKey;
  my $ref = $config;

  for( my $i = 0; $i < scalar @newKey; $i++ )
  {
    return undef unless( ref( $ref ) eq 'HASH' );

    my $matchedKey;
    if( exists $ref->{$newKey[$i]} )
    {
      $matchedKey = $newKey[$i];
    }
    else
    {
      my $lcKey = lc($newKey[$i]);
      foreach my $possibleKey ( keys %$ref )
      {
        if( lc($possibleKey) eq $lcKey )
        {
          $matchedKey = $possibleKey;
          last;
        }
      }
    }

    return undef unless( defined $matchedKey );
    push @canonicalKey, $matchedKey;
    $ref = $ref->{$matchedKey};
  }

  return undef if( ref( $ref ) eq 'HASH' );
  return join '.', @canonicalKey;
}


# Format the effective input key, adding the original legacy key when useful.
sub formatInputKeyContext
{
  my ($effectiveKey, $sourceRef) = @_;
  if( defined $sourceRef && $sourceRef->{'legacy'} )
  {
    return "$effectiveKey (from legacy $sourceRef->{'raw'})";
  }
  return $effectiveKey;
}


# Quote a user input token for diagnostics without changing parser behavior.
sub quoteInputValue
{
  my ($value) = @_;
  $value = '' unless( defined $value );
  $value =~ s/\\/\\\\/g;
  $value =~ s/\n/\\n/g;
  $value =~ s/\t/\\t/g;
  $value =~ s/'/\\'/g;
  return "'$value'";
}


# Convert compact parser type strings into user-facing descriptions.
sub describeType
{
  my ($baseType, $asElement) = @_;
  my $type = 'value';

  if( $baseType =~ m/f/ )
  {
    $type = 'floating point';
  }
  elsif( $baseType =~ m/i/ )
  {
    $type = 'integer';
  }
  elsif( $baseType =~ m/b/ )
  {
    $type = 'boolean';
  }
  elsif( $baseType =~ m/S/ )
  {
    $type = 'case-sensitive string';
  }
  elsif( $baseType =~ m/s/ )
  {
    $type = 'case-insensitive string';
  }

  return $type if( $asElement || $baseType !~ m/a/ );
  return "array of $type values";
}


# Return "value" or "values" to make array length errors read naturally.
sub valueWord
{
  my ($count) = @_;
  return $count == 1 ? 'value' : 'values';
}


# Emit a consistent user-facing type conversion error.
sub dieInputTypeError
{
  my ($key, $sourceRef, $value, $type, $rule, $arrayIndex) = @_;
  my $keyContext = defined $arrayIndex ? "$key\[$arrayIndex\]" : $key;
  $keyContext = formatInputKeyContext( $keyContext, $sourceRef );
  my ($baseType) = split /:/, $type, 2;
  my $typeDescription = describeType( $baseType, defined $arrayIndex );
  my $message = "Invalid value for $keyContext: got " . quoteInputValue( $value )
              . ", expected type $typeDescription";
  $message .= " ($rule)" if( defined $rule && length $rule );
  die "$message\n";
}


# Emit a consistent parser-schema/type-file error.
sub dieTypeSpecError
{
  my ($path, $type, $message) = @_;
  $path = '<unknown>' unless( defined $path && length $path );
  die "Invalid type constraint for $path: $message in $type\n";
}


# Load the sparse enum definition. The enum file is required so a missing
# installation/configuration file cannot silently disable validation.
sub loadEnumDefinition
{
  my ($enumFilename, $json) = @_;

  my $enumDef;
  if( open( my $in, "<", $enumFilename ) )
  {
    local $/ = undef;
    $enumDef = $json->decode(<$in>);
    close($in);
  }
  else
  {
    die "Failed to open $enumFilename\n$!";
  }
  die "Invalid enum definition for <root>: expected JSON object\n"
    unless( ref( $enumDef ) eq 'HASH' );
  return $enumDef;
}


sub dieEnumSpecError
{
  my ($path, $message) = @_;
  $path = '<root>' unless( defined $path && length $path );
  die "Invalid enum definition for $path: $message\n";
}


# Fetch a dot-separated path while distinguishing missing paths from JSON null.
sub getPathValue
{
  my ($tree, $key) = @_;
  return (1, $tree) unless( defined $key && length $key );

  my @path = split /\./, $key;
  my $ref = $tree;
  foreach my $part ( @path )
  {
    return (0, undef) unless( ref( $ref ) eq 'HASH' && exists $ref->{$part} );
    $ref = $ref->{$part};
  }
  return (1, $ref);
}


# Validate the sparse enum tree against oparse.json and oparse.type.json.
sub validateEnumSpecs
{
  my ($enumRef, $config, $typeDef, $prefix) = @_;

  if( ref( $enumRef ) eq 'HASH' )
  {
    if( length $prefix )
    {
      my ($exists, $configValue) = getPathValue( $config, $prefix );
      dieEnumSpecError( $prefix, "path not found in oparse.json" ) unless( $exists );
      dieEnumSpecError( $prefix, "expected object path in oparse.json" )
        unless( ref( $configValue ) eq 'HASH' );
    }

    foreach my $key ( sort keys %$enumRef )
    {
      my $newPrefix = length $prefix ? "$prefix.$key" : $key;
      validateEnumSpecs( $enumRef->{$key}, $config, $typeDef, $newPrefix );
    }
    return;
  }

  dieEnumSpecError( $prefix, "expected allowed-value array or object" )
    unless( ref( $enumRef ) eq 'ARRAY' );
  dieEnumSpecError( $prefix, "allowed-value array cannot be empty" )
    unless( scalar @$enumRef );

  my ($configExists, $configValue) = getPathValue( $config, $prefix );
  dieEnumSpecError( $prefix, "path not found in oparse.json" ) unless( $configExists );

  my ($typeExists, $type) = getPathValue( $typeDef, $prefix );
  dieEnumSpecError( $prefix, "path not found in oparse.type.json" ) unless( $typeExists );
  dieEnumSpecError( $prefix, "enum entries are only valid on leaf type strings" )
    if( ref( $type ) );

  my ($baseType) = parseTypeSpec( $type, $prefix );
  dieEnumSpecError( $prefix, "enum entries are only valid for s or as types, got $baseType" )
    unless( $baseType eq 's' || $baseType eq 'as' );

  my %seenValue;
  foreach my $value ( @$enumRef )
  {
    dieEnumSpecError( $prefix, "allowed values must be strings" )
      if( ! defined $value || ref( $value ) || looks_like_number( $value ) );
    dieEnumSpecError( $prefix, "allowed value '$value' is not lowercase" )
      unless( lc($value) eq $value );
    dieEnumSpecError( $prefix, "duplicate allowed value '$value'" )
      if( exists $seenValue{$value} );
    $seenValue{$value} = 1;
  }
}


sub getEnumAllowedValues
{
  my ($enumDef, $key) = @_;
  my ($exists, $value) = getPathValue( $enumDef, $key );
  return undef unless( $exists );
  dieEnumSpecError( $key, "expected allowed-value array" )
    unless( ref( $value ) eq 'ARRAY' );
  return $value;
}


sub validateInputEnumValue
{
  my ($enumDef, $key, $sourceRef, $value, $arrayIndex) = @_;
  my $allowedRef = getEnumAllowedValues( $enumDef, $key );
  return unless( defined $allowedRef );

  foreach my $allowed ( @$allowedRef )
  {
    return if( defined $value && $value eq $allowed );
  }

  my $keyContext = defined $arrayIndex ? "$key\[$arrayIndex\]" : $key;
  $keyContext = formatInputKeyContext( $keyContext, $sourceRef );
  die "Invalid value for $keyContext: got " . quoteInputValue( $value )
    . ", allowed values: " . join( ', ', @$allowedRef ) . "\n";
}


# Split an oparse.type.json leaf into the original compact base type and any
# comma-separated constraint suffixes. For example:
#   af:len=3
# becomes base type "af" and constraints { len => 3 }.
sub parseTypeSpec
{
  my ($type, $path) = @_;
  my $unsupported = ( $type =~ s/u$// );
  my ($baseType, $constraintString) = split /:/, $type, 2;
  my %constraints;

  if( defined $constraintString && length $constraintString )
  {
    foreach my $constraint ( split /,/, $constraintString )
    {
      $constraint =~ m/^(\w+)=(.+)$/
        or dieTypeSpecError( $path, $type, "malformed constraint '$constraint'" );
      my $name = $1;
      my $value = $2;

      if( $name eq 'len' || $name eq 'signlen' )
      {
        dieTypeSpecError( $path, $type, "$name constraint '$constraint' requires an array type" )
          unless( $baseType =~ m/a/ );
        $value =~ m/^\d+$/
          or dieTypeSpecError( $path, $type, "invalid numeric constraint '$constraint'" );
        dieTypeSpecError( $path, $type, "invalid signlen constraint '$constraint'" )
          if( $name eq 'signlen' && $value == 0 );
        $constraints{$name} = $value * 1;
      }
      elsif( $name eq 'oneof' || $name eq 'lenmatch' )
      {
        dieTypeSpecError( $path, $type, "$name constraint '$constraint' requires an array type" )
          unless( $baseType =~ m/a/ );
        $value =~ m/^[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*$/
          or dieTypeSpecError( $path, $type, "invalid path constraint '$constraint'" );
        $constraints{$name} = $value;
      }
      elsif( $name eq 'lenmul' )
      {
        dieTypeSpecError( $path, $type, "$name constraint '$constraint' requires an array type" )
          unless( $baseType =~ m/a/ );
        $value =~ m/^[1-9]\d*\*[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*$/
          or dieTypeSpecError( $path, $type, "invalid lenmul constraint '$constraint'" );
        $constraints{$name} = $value;
      }
      else
      {
        dieTypeSpecError( $path, $type, "unknown constraint '$name'" );
      }
    }
  }

  return ($baseType, \%constraints, $unsupported);
}


# Walk the full type tree once at startup so malformed constraints in
# oparse.type.json fail even when the associated input is not supplied.
sub validateTypeSpecs
{
  my ($typeRef, $prefix) = @_;
  if( ref( $typeRef ) eq 'HASH' )
  {
    foreach my $key ( keys %$typeRef )
    {
      my $newPrefix = length $prefix ? "$prefix.$key" : $key;
      validateTypeSpecs( $typeRef->{$key}, $newPrefix );
    }
  }
  else
  {
    parseTypeSpec( $typeRef, $prefix );
  }
}


# Collect dot-separated leaf paths from a JSON tree so oparse.json and
# oparse.type.json can be checked for exact schema/type coverage.
sub collectLeafPaths
{
  my ($node, $prefix, $pathsRef) = @_;
  if( ref( $node ) eq 'HASH' )
  {
    foreach my $key ( keys %$node )
    {
      my $newPrefix = length $prefix ? "$prefix.$key" : $key;
      collectLeafPaths( $node->{$key}, $newPrefix, $pathsRef );
    }
  }
  else
  {
    $pathsRef->{$prefix} = 1;
  }
}


# Require oparse.json and oparse.type.json to have exactly the same leaves.
sub validateSchemaTypeMatch
{
  my ($config, $typeDef) = @_;
  my %configPaths;
  my %typePaths;
  my @errors;

  collectLeafPaths( $config, '', \%configPaths );
  collectLeafPaths( $typeDef, '', \%typePaths );

  foreach my $path ( sort keys %configPaths )
  {
    push @errors, "Missing type definition for $path" unless( exists $typePaths{$path} );
  }

  foreach my $path ( sort keys %typePaths )
  {
    push @errors, "Extra type definition for $path" unless( exists $configPaths{$path} );
  }

  die "Schema/type mismatch between oparse.json and oparse.type.json:\n"
    . join( "\n", @errors ) . "\n" if( scalar @errors );
}


# Enforce constraints that can be checked for one parsed array at a time.
# This runs after legacy compatibility fixups, so legacy plot ranges have
# already had the old leading "points" value removed before len=2 is checked.
sub describeArrayLengthRule
{
  my ($key, $type, $constraints, $arrayRef) = @_;
  my ($baseType) = split /:/, $type, 2;
  my $typeDescription = describeType( $baseType, 0 );

  if( exists $constraints->{'len'} )
  {
    my $expected = $constraints->{'len'};
    return "expected exactly $expected " . valueWord( $expected ) . " ($typeDescription)";
  }

  if( exists $constraints->{'signlen'} )
  {
    my $positiveLength = $constraints->{'signlen'};
    if( scalar @$arrayRef == 0 )
    {
      return "expected either 1 value with a negative first value, or "
        . "$positiveLength " . valueWord( $positiveLength )
        . " with a positive first value ($typeDescription)";
    }
    elsif( $arrayRef->[0] < 0 )
    {
      return "expected 1 value because first value is negative ($typeDescription)";
    }
    elsif( $arrayRef->[0] > 0 )
    {
      return "expected $positiveLength " . valueWord( $positiveLength )
        . " because first value is positive ($typeDescription)";
    }
    else
    {
      return "first value is 0, expected a negative value for automatic sizing "
        . "or a positive value for explicit sizing ($typeDescription)";
    }
  }

  return "invalid array length ($typeDescription)";
}


sub validateArrayLength
{
  my ($key, $sourceRef, $type, $constraints, $arrayRef) = @_;
  my $keyContext = formatInputKeyContext( $key, $sourceRef );
  my $length = scalar @$arrayRef;

  if( exists $constraints->{'len'} )
  {
    my $expected = $constraints->{'len'};
    die "Invalid array length for $keyContext: got $length " . valueWord( $length )
      . ", " . describeArrayLengthRule( $key, $type, $constraints, $arrayRef ) . "\n"
      unless( $length == $expected );
  }

  if( exists $constraints->{'signlen'} )
  {
    my $positiveLength = $constraints->{'signlen'};
    die "Invalid array length for $keyContext: got 0 values, "
      . describeArrayLengthRule( $key, $type, $constraints, $arrayRef ) . "\n"
      if( $length == 0 );

    if( $arrayRef->[0] < 0 )
    {
      die "Invalid array length for $keyContext: got $length " . valueWord( $length )
        . ", " . describeArrayLengthRule( $key, $type, $constraints, $arrayRef ) . "\n"
        unless( $length == 1 );
    }
    elsif( $arrayRef->[0] > 0 )
    {
      die "Invalid array length for $keyContext: got $length " . valueWord( $length )
        . ", " . describeArrayLengthRule( $key, $type, $constraints, $arrayRef ) . "\n"
        unless( $length == $positiveLength );
    }
    else
    {
      die "Invalid array length for $keyContext: "
        . describeArrayLengthRule( $key, $type, $constraints, $arrayRef ) . "\n";
    }
  }
}


# Fetch a dot-separated path from the final config tree.
sub getConfigValue
{
  my ($config, $key) = @_;
  my @path = split /\./, $key;
  my $ref = $config;
  foreach my $part ( @path )
  {
    return undef unless( ref( $ref ) eq 'HASH' && exists $ref->{$part} );
    $ref = $ref->{$part};
  }
  return $ref;
}


# Gather constraints that need final config values instead of one raw input array.
sub collectRelationalArrayConstraints
{
  my ($typeRef, $prefix, $oneOfRef, $checksRef) = @_;

  if( ref( $typeRef ) eq 'HASH' )
  {
    foreach my $key ( sort keys %$typeRef )
    {
      my $newPrefix = length $prefix ? "$prefix.$key" : $key;
      collectRelationalArrayConstraints( $typeRef->{$key}, $newPrefix, $oneOfRef, $checksRef );
    }
  }
  else
  {
    my (undef, $constraints) = parseTypeSpec( $typeRef, $prefix );
    push @{ $oneOfRef->{ $constraints->{'oneof'} } }, $prefix
      if( exists $constraints->{'oneof'} );
    push @$checksRef, { key => $prefix, type => $typeRef, constraints => $constraints }
      if( exists $constraints->{'lenmul'} || exists $constraints->{'lenmatch'} );
  }
}


# Enforce relational constraints after all supplied values have been assigned.
# This keeps checks independent of the order in the user's input file.
sub describeRelationalLengthRule
{
  my ($type, $constraintName, $rule, $refKey, $expected) = @_;
  my ($baseType) = split /:/, $type, 2;
  my $typeDescription = describeType( $baseType, 0 );

  if( $constraintName eq 'lenmul' )
  {
    my ($multiplier) = $rule =~ m/^([1-9]\d*)\*/;
    return "expected $expected " . valueWord( $expected ) . ": "
      . "$multiplier " . valueWord( $multiplier )
      . " for each value in $refKey ($typeDescription)";
  }

  if( $constraintName eq 'lenmatch' )
  {
    return "expected $expected " . valueWord( $expected )
      . ": one for each value in $refKey ($typeDescription)";
  }

  return "expected $expected " . valueWord( $expected ) . " ($typeDescription)";
}


sub validateRelationalArrayLengths
{
  my ($config, $typeDef, $suppliedRef, $sourceMapRef) = @_;
  my %oneOf;
  my @checks;

  collectRelationalArrayConstraints( $typeDef, '', \%oneOf, \@checks );

  foreach my $group ( sort keys %oneOf )
  {
    my @members = @{ $oneOf{$group} };
    my @supplied = grep { exists $suppliedRef->{$_} } @members;
    my @suppliedContext = map {
      formatInputKeyContext( $_, defined $sourceMapRef ? $sourceMapRef->{$_} : undef )
    } @supplied;
    die "Invalid array selection for $group: got " . scalar @supplied
      . " supplied inputs"
      . ( scalar @supplied ? " (" . join( ', ', @suppliedContext ) . ")" : "" )
      . ", expected exactly one of " . join( ', ', @members ) . "\n"
      unless( scalar @supplied == 1 );
  }

  foreach my $check ( @checks )
  {
    my $key = $check->{'key'};
    next unless( exists $suppliedRef->{$key} );
    my $keyContext = formatInputKeyContext( $key,
      defined $sourceMapRef ? $sourceMapRef->{$key} : undef );

    my $value = getConfigValue( $config, $key );
    my ($baseType) = split /:/, $check->{'type'}, 2;
    my $typeDescription = describeType( $baseType, 0 );
    die "Invalid array length for $keyContext: final value is missing or is not an array ($typeDescription)\n"
      unless( ref( $value ) eq 'ARRAY' );
    my $length = scalar @$value;

    if( exists $check->{'constraints'}->{'lenmul'} )
    {
      my $rule = $check->{'constraints'}->{'lenmul'};
      $rule =~ m/^([1-9]\d*)\*(.+)$/;
      my $multiplier = $1;
      my $refKey = $2;
      my $refValue = getConfigValue( $config, $refKey );
      die "Invalid array length for $keyContext: reference $refKey is missing or is not an array ($typeDescription)\n"
        unless( ref( $refValue ) eq 'ARRAY' );
      my $expected = $multiplier * scalar @$refValue;
      die "Invalid array length for $keyContext: got $length " . valueWord( $length )
        . ", " . describeRelationalLengthRule( $check->{'type'}, 'lenmul', $rule, $refKey, $expected ) . "\n"
        unless( $length == $expected );
    }

    if( exists $check->{'constraints'}->{'lenmatch'} )
    {
      my $refKey = $check->{'constraints'}->{'lenmatch'};
      my $refValue = getConfigValue( $config, $refKey );
      die "Invalid array length for $keyContext: reference $refKey is missing or is not an array ($typeDescription)\n"
        unless( ref( $refValue ) eq 'ARRAY' );
      my $expected = scalar @$refValue;
      die "Invalid array length for $keyContext: got $length " . valueWord( $length )
        . ", " . describeRelationalLengthRule( $check->{'type'}, 'lenmatch', undef, $refKey, $expected ) . "\n"
        unless( $length == $expected );
    }
  }
}
