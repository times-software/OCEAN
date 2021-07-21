#!/usr/bin/env perl
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

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/legacyBuildJSON\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }


my $json = JSON::PP->new;
my $input;

if( open( my $in, "oparse.json" ))
{
  local $/ = undef;
  $input = $json->decode(<$in>);
  close($in);
} 
else
{
  die "Failed to open oparse.json which should have been copied in from the OCEAN install directory!\n$!";
}

my $typeDef;
if( open( my $in, "oparse.type.json" ))
{
  local $/ = undef;
  $typeDef = $json->decode(<$in>);
  close($in);
}
else
{
  die "Failed to open oparse.type.json\n$!";
}

my %decoder = (
  'calc' => 'calc.mode',
  'dft' => 'dft.program',
  'nkpt' => 'bse.kmesh',
  'ngkpt' => 'dft.den.kmesh',
  'ngkpt.auto' => 'dft.den.auto',
  'qinunitsofbvectors.ipt' => 'calc.photon_q',
  'dft.split' => 'dft.bse.split',
  'dft.qe_redirect' => 'dft.redirect',
  'nbands' => 'bse.nbands',
  'dft_energy_range.ipt' => 'bse.dft_energy_range',
  'obf_energy_range' => 'nope.obf_energy_range',
  'obkpt.ipt' => 'nope.obkpt',
  'obf.nbands' => 'nope.obf_nbands',
  'trace_tol' => 'nope.trace_tol',
  'acc_level.ipt' => 'nope.acc_level',
  'k0.ipt' => 'bse.kshift',
  'fband' => 'dft.fband',
  'occopt' => 'dft.occopt',
  'mixing' => 'dft.mixing',
  'rscale' => 'structure.rscale',
  'rprim' => 'structure.rprim',
  'ntype' => 'nope.ntype',
  'typat' => 'structure.typat',
  'znucl' => 'structure.znucl',
  'zsymb' => 'structure.zsymb',
  'pplist' => 'psp.pp_list',
  'ppdatabase' => 'psp.pp_database',
  'ecut.qualtiy' => 'psp.ecut_quality',
  'natoms' => 'nope.natoms',
  'coord' => 'structure.coord',
  'taulist' => 'structure.taulist',
  'ecut' => 'dft.ecut', 
  'epsilon' => 'structure.epsilon',
  'etol' => 'dft.toldfe',
  'wftol' => 'dft.tolwfr',
  'nrun' => 'dft.nstep',
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
  'nedges.ipt' => 'nope.nedges',
  'edges.ipt' => 'calc.edges',
  'nbuse.ipt' => 'nope.nbuse',
  'xmesh.ipt' => 'bse.xmesh',
  'cnbse.rad' => 'bse.core.screen_radius',
  'metal' => 'structure.metal',
  'cksshift' => 'nope.cksshift',
  'cksstretch' => 'nope.cksstretch',
  'cnbse.niter' => 'bse.core.haydock.niter',
  'haydockconv.ipt' => 'bse.core.haydock.converge.thresh',
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
  'gw_control' => 'bse.core.gw.control',
  'gwcstr' => 'nope.bse_gw_cstr',
  'gwvstr' => 'nope.bse_gw_vstr',
  'gwgap' => 'nope.gwgap',
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
  'core_offset' => 'screen.core_offset',
  'ham_kpoints' => 'nope.ham_kpoints',
  'niter' => 'bse.val.haydock.niter',
  'backf' => 'bse.val.backf',
  'aldaf' => 'bse.val.aldaf',
  'qpflg' => 'bse.val.qpflg',
  'bwflg' => 'bse.val.bwflg',
  'bande' => 'bse.val.bande',
  'bflag' => 'bse.val.bflag',
  'lflag' => 'bse.val.lflag',
  'convergence' => 'nope.convergence',
  'decut' => 'bse.val.decut',
  'se_rs' => 'nope.se.rs',
  'se_metal' => 'nope.se.metal',
  'se_niter' => 'nope.se.niter',
  'spect.h' => 'bse.val.plot.range',
  'tot_charge' => 'dft.tot_charge',
  'nspin' => 'dft.nspin',
  'smag' => 'dft.smag',
  'ldau' => 'dft.ldau.Hubbard_U',
  'qe_scissor' => 'nope.qe_scissor',
  'nphoton' => 'nope.nphoton',
  'serbse' => 'nope.ser_bse',
  'spin_orbit' => 'bse.core.spin_orbit',
  'screen_energy_range.ipt' => 'screen.dft_energy_range',
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
  'echamp.inp' => 'bse.core.gmres.echamp',
  'bshift' => 'nope.bshift' );
  

# TEST MODE -- make sure all the input flags are correctly here
if( 0 )
{
  open IN, "oparse.h" or die;
  while( my $line = <IN> )
  { 
    if( $line =~ m/[\{\}]/ ) { die "Mix up in oparse.h: expected input\n$line";}
    $line = <IN>;
    die "$line. Mix up in oparse.h\n" unless( $line =~ m/\{/ );
    $line = <IN>;
    print "$line";
    $line = <IN> unless( $line =~ m/\}/ );
    $line =~ m/\}\s*(\S+)/ or die "$line";
    my $tag = $1;
    chomp($tag);
    die "Failed to find tag: $tag\n" unless( exists $decoder{ $tag } );
  }
}

if( 1 )
{
  my $enable = 1;
  $json->canonical([$enable]);
  $json->pretty([$enable]);

  open OUT, ">", "json.test" or die;
  print OUT $json->encode($input);
  close OUT;
}


print $typeDef->{'calc'} . "\n";
print $typeDef->{'calc'}->{'mode'} . "\n";

foreach my $key (keys %decoder)
{
  print "$key: ";
  die "Failed to find $key\n" unless -e $key;
  my $value ="";
  if( open( my $in, $key ))
  {
    local $/ = undef;
    $value .= <$in>;
    $value =~ s/\n/ /g;
    close($in);
  }
  else
  {
    die "Failed to open $key!\n$!";
  }
  my $newKey = $decoder{ $key };
  my @newKey = split /\./, $newKey;
  foreach my $i (@newKey )
  {
    print "$i ";
  }
  print "\n";

  
  my $type;
  my $hashref;
  if( $newKey[0] =~ m/nope/ )
  {
    $type = 's';
  }
  else
  {
    $type = $typeDef;
    for( my $i = 0; $i < scalar @newKey; $i++ )
    {
      $type = $type->{$newKey[$i]} ;
    }
    print "   $type ";
    $hashref = $input;
    for( my $i = 0; $i < scalar @newKey - 1; $i++ )
    {
      $hashref = $hashref->{$newKey[$i]};
    }
  }
  die "Failed! >$type<\n" unless( $type =~ m/\As|S|i|f|b|as|ai|af\Z/ );


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
      print "$i ";
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
    

    # Fix incompatible!
    if( $key =~ m/spect\.h|cnbse\.spect_range/ )
    {
      shift @rawArray;
    }
    # end fix
    $hashref->{$newKey[-1]} = [@rawArray];
  }
  else
  {
    unless( $value eq ' ' )
    {
      $value =~ s/^\s+//;
      $value =~ s/\s+$//;
    }

    # fix incompatible!
    if( $key =~ m/epsilon/ )
    {
      if( $value =~ m/dfpt/ )
      {
         $value = 0;
         $input->{'dft'}->{'epsilon'}->{'method'} = 'dfpt'
      }
      else
      {
         $input->{'dft'}->{'epsilon'}->{'method'} = 'input'
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
    print "$value ";
    $hashref->{$newKey[-1]} = $value;
  }
  print "\n";
}
  

  my $enable = 1;
  $json->canonical([$enable]);
  $json->pretty([$enable]);

  open OUT, ">", "jsonOut.test" or die;
  print OUT $json->encode($input);
  close OUT;
