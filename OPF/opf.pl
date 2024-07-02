#!/usr/bin/perl
# Copyright (C) 2014, 2016-2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#

use strict;
use File::Copy;
use File::Compare;
use File::Spec::Functions;
use Cwd;
use Cwd 'abs_path';
require JSON::PP;
JSON::PP->import;
use Storable qw(dclone);
use Time::HiRes qw( gettimeofday tv_interval );
use Digest::MD5 qw(md5_hex);


###########################
if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/opf\.pl/;
#  $ENV{"OCEAN_BIN"} = $1;
  $ENV{"OCEAN_BIN"} = abs_path( $1 );
  print "OCEAN_BIN not set. Setting it to $ENV{'OCEAN_BIN'}\n";
}

my $dir = getcwd;
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
###########################

my $dataFile = "../Common/postDefaultsOceanDatafile";
if( -e $dataFile )
{
  my $t0 = [gettimeofday];
  my $didAnything = 0;
  my $json = JSON::PP->new;
  my $commonOceanData;
  if( open( my $in, "<", $dataFile ))
  {
    local $/ = undef;
    $commonOceanData = $json->decode(<$in>);
    close($in);
  }
  else
  {
    die "Failed to open config file $dataFile\n$!";
  }


  my %todoHash;
  my %todoPsp;
  my %todoZNL;
  my %todoInput;
  if( $commonOceanData->{'calc'}->{'mode'} =~ m/val/i )
  {
    print "Valence run! Valence screening augmentation is not automatic!\n";
    exit 0;
    #TODO Would need to add every unique z to a list, and maybe try and disable calculating the core response?
  }
  else
  {
    foreach my $edge ( @{$commonOceanData->{'calc'}->{'edges'}} )
    {
      $edge =~ m/(\d+)\s(\d)\s(\d)/ or die;
      my $site = $1;
      my $N = $2;
      my $L = $3;

      # change from 1-index to 0-index
      my $type = @{$commonOceanData->{'structure'}->{'typat'}}[$site-1];
      my $Z = @{$commonOceanData->{'structure'}->{'znucl'}}[$type-1];
      my $psp_hash = @{$commonOceanData->{'psp'}->{'pphash'}}[$type-1];
      my $psp_name = @{$commonOceanData->{'psp'}->{'pp_list'}}[$type-1];

      print "$Z $N $L $psp_hash\n";
      
      $todoHash{ $Z } = $psp_hash;
      $todoPsp{ $Z } = $psp_name;
      
      my $nl = "$N $L";
      if( exists $todoZNL{ $Z } )
      {
        $todoZNL{ $Z }{ $nl } = 1;
      }
      else
      {
        $todoZNL{ $Z } = {};
        $todoZNL{ $Z }{ $nl } = 1;
      }
    }
    foreach my $Z (keys %todoHash )
    {
      print "$todoHash{ $Z } ";
      foreach (keys %{$todoZNL{ $Z }} )
      {
        print "$_ ";
      }
      print "\n";

      # input hash
      if( $commonOceanData->{'opf'}->{'program'} eq 'shirley' ) {
        $todoInput{ $Z } = &inputHashShirley( $Z, $commonOceanData );
      } elsif( $commonOceanData->{'opf'}->{'program'} eq 'hamann' ) {
        $todoInput{ $Z } = &inputHashHamann( $Z, $todoPsp{ $Z }, $commonOceanData );
      } else { 
        die "Un-supported opf->program\n";
      }
    }
  }

  my $runClean = 1;
  my $opfData;
  $dataFile = "opf.json";
  if( -e $dataFile && open( my $in, "<", $dataFile ))
  {
    local $/ = undef;
    $opfData = $json->decode(<$in>);
    close($in);

    #TODO checks for other opf options as they are added in
    goto CLEAN if( $commonOceanData->{'opf'}->{'program'} ne $opfData->{'program'} );
    if( $commonOceanData->{'opf'}->{'program'} eq "shirley" )
    {
      if( ( $commonOceanData->{'opf'}->{'shirley'}->{'caution'} && ! $opfData->{'shirley'}->{'caution'} ) 
       || ( ! $commonOceanData->{'opf'}->{'shirley'}->{'caution'} && $opfData->{'shirley'}->{'caution'} ) )
      {
        goto CLEAN;
      }
      if( $commonOceanData->{'opf'}->{'shirley'}->{'hfkgrid'}[0] != $opfData->{'shirley'}->{'hfkgrid'}[0] 
       || $commonOceanData->{'opf'}->{'shirley'}->{'hfkgrid'}[1] != $opfData->{'shirley'}->{'hfkgrid'}[1] )
      {
        goto CLEAN;
      }
    }

    goto CLEAN if( abs( $commonOceanData->{'opf'}->{'radius'} - $opfData->{'radius'} ) > 0.00000001 );

    $runClean = 0;
  }
  CLEAN:
    
  if( $runClean )
  {
    print "Clean run is necessary!!\n";
    $opfData = dclone $commonOceanData->{'opf'};
    my $enable = 1;
    $json->canonical([$enable]);
    $json->pretty([$enable]);
    open OUT, ">", "opf.json" or die;
    print OUT $json->encode($opfData);
    close OUT;
  }


  if( exists $opfData->{'completed'} )
  {
    print "Previous run detected\n";

    foreach my $Z (keys %{$opfData->{'completed'}} )
    { 
      next unless( exists $todoHash{ $Z } );
      my $delete = 1;
      unless( $todoHash{ $Z } eq $opfData->{'completed'}->{$Z}->{'psp_hash'} ) {
        print "Psp hash has changed: " . $todoHash{ $Z } . " " . $opfData->{'completed'}->{$Z}->{'psp_hash'} . "\n";
        $delete = 0;
      }
      if( exists $opfData->{'completed'}->{$Z}->{'input_hash'} ) { 
        unless( $todoInput{$Z} eq $opfData->{'completed'}->{$Z}->{'input_hash'} ) {
          print "Psp inputs have changed: " . $todoInput{$Z} . " " . $opfData->{'completed'}->{$Z}->{'input_hash'} . "\n";
        $delete = 0;
        }
      } else { 
        print "No previous input hash\n";
        $delete = 0;
      }
      
#      print "$todoHash{ $Z }  ---  $opfData->{'completed'}->{$Z}->{'psp_hash'}\n";
#      if( $todoHash{ $Z } eq $opfData->{'completed'}->{$Z}->{'psp_hash'} && 
#          exists $opfData->{'completed'}->{$Z}->{'input_hash'} &&
#          $todoHash{$Z} eq $opfData->{'completed'}->{$Z}->{'input_hash'} )
#      {
      if( $delete ) {
        foreach my $nl ( @{$opfData->{'completed'}->{$Z}->{'NL'}} )
        {
#          print "   $nl\n";
 #         foreach my $i ( keys %{$todoZNL{ $Z }} )
 #         {  print "   $i\n";
 #         }
          delete $todoZNL{ $Z }{ $nl } if( exists $todoZNL{ $Z }{ $nl } );
        }
      }
      else
      {
        delete $opfData->{'completed'}->{ $Z };
      }
    }
  }
  
  foreach my $Z ( keys %todoZNL )
  {
    if( scalar keys %{$todoZNL{ $Z }} == 0 )
    {
      delete $todoZNL{ $Z };
      delete $todoHash{ $Z };
    }
    else
    {
      my $j = scalar keys %{$todoZNL{ $Z }};
      print "### $j\n";
      print "$Z: ";
      foreach my $i ( keys %{$todoZNL{ $Z }} )
      { print "$i ";
      }
      print "\n";
    }
  }
   
  # write out new opf.json here!
  # only needed if we didn't already clean everything
  unless( $runClean )
  {
    my $enable = 1;
    $json->canonical([$enable]);
    $json->pretty([$enable]);
    open OUT, ">", "opf.json" or die;
    print OUT $json->encode($opfData);
    close OUT;
  }

  $opfData->{'completed'} = {}  unless( exists $opfData->{'completed'} );
  # then fake the doing and write out result for testing
  foreach my $Z (keys %todoHash )
  {
    next unless( scalar keys %{$todoZNL{ $Z }} > 0 );
    $didAnything = 1;

    print "Need to run Z = $Z with hash $todoHash{ $Z }\n" if( scalar keys %{$todoZNL{ $Z }} > 0 );
    $opfData->{'completed'}->{$Z} = {} unless( exists $opfData->{'completed'}->{$Z} );
    # FAKE
    if( $opfData->{'program'} =~ m/shirley/ )
    {
      runShirley( $Z, $todoPsp{ $Z }, $commonOceanData ); 
    }
    elsif( $opfData->{'program'} =~ m/ham/ )
    { 
      runONCV( $Z, $todoPsp{ $Z }, $commonOceanData );
    }
    else 
    {   
      die "Unsupported opf aux program\n";
    }
    # Can just write over since we already checked that they're the same
    $opfData->{'completed'}->{$Z}->{'psp_hash'} = $todoHash{ $Z };
    $opfData->{'completed'}->{$Z}->{'input_hash'} = $todoInput{ $Z };
    $opfData->{'completed'}->{$Z}->{'NL'} = [] unless( exists $opfData->{'completed'}->{'NL'} );
    foreach my $nl (keys %{$todoZNL{ $Z }} )
    {
      print "  Run for NL: $nl\n";
      push @{$opfData->{'completed'}->{$Z}->{'NL'}}, $nl;
    }
  }

  if( $didAnything ) {
    $opfData->{'time'} = tv_interval( $t0 );
  }
  my $enable = 1;
  $json->canonical([$enable]);
  $json->pretty([$enable]);
  open OUT, ">", "opf.json" or die;
  print OUT $json->encode($opfData);
  close OUT;
  exit 0;
}


my @CommonFiles = ("znucl", "opf.hfkgrid", "opf.fill", "opf.opts", "pplist", "ppdir", 
                   "ntype", "natoms", "typat", "taulist", "nedges", "edges", "caution", 
                   "scfac", "opf.program", "coord", "avecsinbohr.ipt" );
my @ExtraFiles = ("calc");


my $runOPF = 1;
if (-e "done" ) {
  $runOPF = 0;
  foreach (@CommonFiles) {
#    if (`diff -q $_ ../Common/$_`) {
    if( compare("$_", "../Common/$_" ) != 0 )   # should get diff or non-exist
    {
      print "$_ differs\n";
      $runOPF = 1;
      last;
    }
  }
}

# Two very different ways of doing things depending on which wether we use ONCVPSP or Shirley codes
my $program = 'shirley';
open IN, "../Common/opf.program" or die "Failed to open opf.program: $!\n";
if( <IN> =~ m/hamann/i )
{
  $program = 'hamann';
  print "Will run OPF using the oncvpsp.x code\n";
}
else
{
  print "Will run OPF using hfk.x\n";
}
close IN;

# test additional files
if( $runOPF == 0 )
{
  if( $program eq 'hamann' )
  {
    print "hamann\n";

  }
  else
  {
    # This might grab too much when we fix multi-element runs
    open IN, "opf.fill" or die "Failed to open opf.fill\n$!";
    while( my $line = <IN> )
    {
      if( $line =~ m/\d+\s+(\S+)/ )
      {
        my $file = $1;
        chomp $file;
        if( compare("$file", "../$file" ) != 0 )
        {
          print "$file differs\n";
          $runOPF = 1;
          last;
        }
      }
      else
      { print $line; }
    }

    if( $runOPF == 0 )
    {
      open IN, "opf.opts" or die "Failed to open opf.opts\n$!";
      while( my $line = <IN> )
      {
        if( $line =~ m/\d+\s+(\S+)/ )
        {
          my $file = $1;
          chomp $file;
          if( compare("$file", "../$file" ) != 0 )
          {
            print "$file differs\n";
            $runOPF = 1;
            last;
          }
        }
      }
    }
  }
}

if ($runOPF == 0 ) {
  print "Nothing new needed for OPF stage\n";
  open OUT, ">", "old" or die;
  print OUT "1\n";
  close OUT;
  exit 0;
}


unlink "done";
unlink "old";


foreach(@ExtraFiles)
{
  copy( "../Common/$_", "$_" ) == 1 or die "Failed to get $_ from Common/\n";
}

foreach (@CommonFiles) {
  copy( "../Common/$_", "$_") == 1 or die "Failed to get $_ from Common/\n";
}

my $allaug = 0;
if( -e "../Common/screen.allaug" ) {
  copy( "../Common/screen.allaug", "screen.allaug" )  == 1 or die "Failed to get screen.allaug from Common/\n";
  if( open IN, "screen.allaug" ) {
    if( <IN> =~ m/T/i )
    { $allaug = 1; }
    close IN;
  }
}

if( open CALC, "calc" )
{
  if( <CALC> =~ m/VAL/i )
  {
    if( $allaug == 0 ) {
      print "No OPF calc for valence run\n";
      close CALC;
      exit 0;
    }
    close CALC;
    ## need OPF run for all sites for all aug

  }
}


###################################
# Quickcheck
if( $program =~ m/shirley/ )
{
  open IN, "opf.opts" or die "Failed to open opf.opts: $!\n";
  if( <IN> =~ m/^!/ )
  {
    close IN;
    die "The input parameter opf.opts was not specified!\n";
  }
  close IN;

  open IN, "opf.fill" or die "Failed to open opf.fill: $!\n";
  if( <IN> =~ m/^!/ )
  {
    close IN;
    die "The input parameter opf.fill was not specified!\n";
  }
  close IN;
}

# Setup
###################################
print "Running OPF Setup\n";
system("$ENV{'OCEAN_BIN'}/pawsetup.x") == 0 or die "Failed to run pawsetup.x\n";

unless( -d "zpawinfo" )
{
  mkdir "zpawinfo" or die "$!";
}
###################################


if( $program =~ m/shirley/ )
{
  my $ppfilename;
  my $ppmodname;

  # Load up pspoptions
  ###################################
  open PSPO, "pspoptions" or die "Failed to open pspoptions\n";
  my %pspopts;
  my %pspfill;
  my %psplist;
  while (<PSPO>) {
    $_ =~ m/\s*(\d+)\s+(\w+)\s+(\S+)\s+(\S+)\s+(\S+)/ or die "Malformed pspoptions\n";
    $psplist{"$1"} = $3;
    $pspopts{"$1"} = $4;
    $pspfill{"$1"} = $5;
  #  copy( "../$3", "$3" ) == 1 or die "Failed to copy $3\n$1";
    copy( "../$4", "$4" ) == 1 or die "Failed to copy $4\n$!";
    copy( "../$5", "$5" ) == 1 or die "Failed to copy $5\n$!";
  }
  close PSPO;

  my %skipValidation;
  # PP mod
  ###################################
  foreach  my $znucl (keys %psplist )
  {
    $ppfilename = $psplist{$znucl};
    $ppmodname = $ppfilename . ".mod";
    print "Need $ppmodname\n\n";
    if( -e "../$ppmodname" )
    {
      print "Using found modified psps\n";
      copy( "../$ppmodname", $ppmodname ) == 1 or die "Failed to copy $ppmodname\n$!";
      $skipValidation{$znucl} = 1;
    }
    else
    {
      copy( "../$ppfilename", $ppfilename ) == 1 or die "Failed to copy $ppfilename\n$!";
      print "Converting $ppfilename\n";
      system("echo '$ppfilename\n$ppmodname' | $ENV{'OCEAN_BIN'}/fhi2eric.x") == 0
          or die "Failed to convert psp file $ppfilename\n";
    }
  }

  # shells
  ##################################
  if( 0 ) {
  open SHELLS, "screen.shells" or die "Failed to open screen.shells\n";
  my $numshells = 0;
  my $allshells = '';
  while (<SHELLS>) {
    chomp;
    $allshells .= $_ ." ";
  }
  close SHELLS;
  my @rads = split(/\s+/, $allshells);
  $numshells = $#rads + 1;
  open SHELLS, ">shells" or die "Failed to open shells for writing\n";
  print SHELLS "$numshells\n";
  print SHELLS "$allshells\n";
  close SHELLS;
  }

  # HFK
  ###################################
  print "Starting HFK section\n";

  open GRID, "opf.hfkgrid" or die;
  my $grid = <GRID>;
  chomp $grid;
  close GRID;

  my $optionfilename;

  foreach my $znucl (keys %psplist )
  {
    $ppfilename = $psplist{$znucl};
    $optionfilename = $pspopts{"$znucl"};
  #  print "$optionfilename\n";
    die "No option file for $ppfilename\n" unless (-e $optionfilename );
    copy( $optionfilename, "atomoptions" ) == 1 or die "Failed to copy $optionfilename to atomoptions\n$!";
    copy( "${ppfilename}.mod", "ppot" ) == 1 or die "Failed to copy ${ppfilename}.mod to ppot\n$!";

    if( $skipValidation{$znucl} == 1 )
    {
      print "  Skipping validation for $znucl. Pre-modified psp was provided\n";
    }
    else
    {
      system( "$ENV{'OCEAN_BIN'}/validate_opts.pl ${ppfilename} $optionfilename" ) == 0 
        or die "Failed to validate options file\nCheck $optionfilename\n";
    }

    open HFIN, ">hfin1" or die;
    print HFIN "initgrid\n";
    print HFIN "$znucl $grid\n";
    print HFIN "ppload\nmkcorcon\nscreencore\nquit\n";
    close HFIN;

    open HFIN, ">hfin2" or die;
    print HFIN "initgrid\n";
    print HFIN "$znucl $grid\n";
    print HFIN "ppload\nmkcorcon\ncalcso\nquit\n";
    close HFIN;

    open HFIN, ">hfin3" or die;
    print HFIN "initgrid\n";
    print HFIN "$znucl $grid\n";
    print HFIN "ppload\nmkcorcon\nspartanfip\n";
    print HFIN `cat "$pspfill{"$znucl"}"`;
    print HFIN "quit\n";
    close HFIN;


    print "Running hfk.x\n";


    system("$ENV{'OCEAN_BIN'}/hfk.x < hfin1 > hfk.${znucl}.1.log") == 0 or die;
    # Check the end of the log to see if we are ok
    my $hfk_status = `tail -n 1 hfk.${znucl}.1.log`;
    unless( $hfk_status =~ m/terminus/ )
    {
      die "The program hfk.x has exited incorrectly for hfin1.\nExiting ...\n";
    }

    system("$ENV{'OCEAN_BIN'}/hfk.x < hfin2 > hfk.${znucl}.2.log") == 0 or die;
    $hfk_status = `tail -n 1 hfk.${znucl}.2.log`;
    unless( $hfk_status =~ m/terminus/ )
    {
      die "The program hfk.x has exited incorrectly for hfin2.\nExiting ...\n";
    }
    my $corezfile = sprintf("corezetaz%03i",$znucl);
    move("xifile","$corezfile");

    system("$ENV{'OCEAN_BIN'}/hfk.x < hfin3 > hfk.${znucl}.3.log") == 0 or die;
    $hfk_status = `tail -n 1 hfk.${znucl}.3.log`;
    unless( $hfk_status =~ m/terminus/ )
    {
      die "The program hfk.x has exited incorrectly for hfin3.\nExiting ...\n";
    }
    print "Done running hfk.x\n";


    unless( -d "zdiag_${znucl}" )
    {
      ( mkdir "zdiag_${znucl}" ) or die "$!";
    }

  # Wander through the directory making some changes
    opendir DIR, $dir;
    while( my $file = readdir(DIR) )
    {
      if( ( $file =~ m/^[fg]k/ ) or ( $file =~ m/^(ae|ft|ps).z/ ) or
          ( $file =~ m/^(deflin|melfile|corez)/ ) or ( $file =~ m/^(rad|prj)filez/ ) or
          ( $file =~ m/^(vcallel|vvallel|vc_bare|vpseud|valence)/ ) or 
          ( $file =~ m/^coreorb/ ) )
      {
        move( $file, "zpawinfo" );
      }
      elsif( $file =~ m/^phr/ )
      {
        my $dest = sprintf( "${file}z%03i", $znucl );
        move( $file, "zpawinfo/$dest" );
      }
      elsif( ( $file =~ m/^melfilez\w$/ ) or ( $file =~ m/^(sm|am|di|pr|psft|aeft)\w$/ ) or
             ( $file =~ m/^(mt|dif)\w\w$/ ) or ( $file =~ m/^(map|ex)/ ) or ( $file =~ /hfin\d/ ) or
             ( $file =~ m/hfk.+log/ ) or ( $file =~ m/aetotal/ ) or ( $file =~ m/radf/ ) )
      {
        move( $file, "zdiag_${znucl}" );
      }
      elsif( ( $file =~ m/^angs\d$/ ) or ( $file =~ m/^ldep\d$/ ) or ( $file =~ m/^shellR\d\.\d\d$/ ) )
      {
        unlink $file;
      }
    }

    my @paw_files_to_kill = ( 'radpot', 'hapot', 'hfpot', 'config', 'corcon', 'valcon', 'chgsum',
                              'atmf99', 'skip', 'psld', 'rslt', 'tmp', 'aprog', 'vxcofr', 'xifile',
                              'actual', 'ppot' );
    foreach (@paw_files_to_kill)
    {
      unlink $_;
    }


  }


}
######################################
else  # oncvpsp method
{

  # Load up pspoptions
  ###################################
  open PSPO, "pspoptions" or die "Failed to open pspoptions\n";
  my %psplist;
  while (<PSPO>) {
    $_ =~ m/\s*(\d+)\s+(\w+)\s+(\S+)\s+(\S+)\s+(\S+)/ or die "Malformed pspoptions\n";
    $psplist{"$1"} = $3;
  }
  close PSPO;

  open GRID, "opf.hfkgrid" or die;
  my $grid = <GRID>;
  chomp $grid;
  close GRID;

  open IN, "ppdir" or die "Failed to open ppdir\n$!";
  my $ppdir = <IN>;
  chomp $ppdir;
  close IN;

  ###################################
  foreach  my $znucl (keys %psplist )
  {
    my $oncvpspInputFile = $psplist{"$znucl"};
    $oncvpspInputFile =~ s/.upf$//i;
    $oncvpspInputFile .= ".in";
    my $targ = catdir( "$ppdir", "$oncvpspInputFile" );
    unless( -e $targ )
    {
      print "Was trying to find input for $psplist{$znucl}:\t$oncvpspInputFile\n";
      print "$targ\n";
      die;
    }
    # ppdir
    copy( $targ, $oncvpspInputFile);
    my $targRad = -1;
    if( -e "overrideRadius" )
    {
      open IN, "overrideRadius" or die "$!";
      $targRad = <IN>;
      chomp $targRad;
      close IN;
    }
    open IN, "scfac" or die "$!";
    my $scfac = <IN>;
    chomp $scfac;
    close IN;
    open IN, ">>", $oncvpspInputFile or die "Failed to open $oncvpspInputFile\n$!";
    print IN ".true.\n$targRad   $scfac\n";
    close IN;

    my @oncvpsp;
    open IN, $oncvpspInputFile or die "Failed to open $oncvpspInputFile\n$!";
    while( my $line = <IN> )
    {
      unless( $line =~ m/^\s*#/ )
      {
        chomp $line;
        push @oncvpsp, $line;
      }
    }
    close IN;

    open OUT, ">atomoptions" or die "Failed to open atomoptions for writing\n";

    $oncvpsp[0] =~ m/(\w+)\s+(\d+\.?\d*)\s+(\d+)\s+(\d+)\s+(-?\d+)\s+\w+/ 
      or die "Failed reading $oncvpspInputFile\t$oncvpsp[0]\n";
    my $zee = $2;
    my $nc = $3;
    my $nv = $4;
    if( $zee != $znucl )
    {
      print "Mismatch between Z in psp.in and Z in input file\n$zee\t$znucl\n";
      die;
    }
    print OUT "$zee\n";
    
    my $lmax = 0;
    my @spdf;
    $spdf[0] = 0; $spdf[1] = 0; $spdf[2] = 0; $spdf[3] = 0;
    for( my $i = 1; $i <= $nc; $i++ )
    {
      print "CORE: $oncvpsp[$i]\n";
      $oncvpsp[$i] =~ m/^\s*(\d+)\s+(\d+)/ 
          or die "Failed reading core levels of $oncvpspInputFile\t$oncvpsp[$i]\n";
      $spdf[$2]++;
      $lmax = $2 if( $2 > $lmax );
    }

    print OUT "$spdf[0]   $spdf[1]   $spdf[2]   $spdf[3]\n";

    $spdf[0] = 0; $spdf[1] = 0; $spdf[2] = 0; $spdf[3] = 0;
    print OUT "scalar\nrel\nlda\n";
    for( my $i = $nc+1; $i <= $nc+$nv; $i++ )
    {
      print "    : $oncvpsp[$i]\n";
      $oncvpsp[$i] =~ m/^\s*(\d+)\s+(\d+)\s+(\S+)/ 
          or die "Failed reading core levels of $oncvpspInputFile\t$oncvpsp[$i]\n";
      $spdf[$2] += $3;
      if( $3 > 0.000001 ) 
      {
        $lmax = $2 if( $2 > $lmax );
      }
    }   
    close IN;

    # reset to max per principal
    $spdf[0] = 2  if( $spdf[0] > 2 );
    $spdf[1] = 6  if( $spdf[1] > 6 );
    $spdf[2] = 10 if( $spdf[2] > 10 );
    $spdf[3] = 14 if( $spdf[3] > 14 );

    # This hack is because, currently, hfk.x is used to screen the core
    #  and for some neutral atoms it doesn't work. 
    for( my $i = 0; $i <=3; $i++ )
    {
      $spdf[$i] *= 0.95;
    } 

    print OUT "$spdf[0]   $spdf[1]   $spdf[2]   $spdf[3]\n";
    print OUT "$spdf[0]   $spdf[1]   $spdf[2]   $spdf[3]\n";
    close OUT;


    my $log = $psplist{"$znucl"} . ".log";
    system("$ENV{'OCEAN_BIN'}/oncvpsp.x < $oncvpspInputFile > $log") == 0 
        or die "Failed to run oncvpsp.x. Check logfile OPF/$log\n\$!";
    unless( move( "ocean.mod", "ppot" ) ) 
    {
      open IN, "<", $log or die "It appears oncvpsp.x has not run correctly\n";
      my $correctONCVPSP = 0;
      while( my $line = <IN> ) 
      {
        if( $line =~ m/OCEAN/ ) 
        {
          $correctONCVPSP = 1;
          last;
        }
      }
      close IN;
      die "oncvpsp.x did not finish correctly. Check the logfile OPF/$log\n" if( $correctONCVPSP );
      die "Incorrect version of oncvpsp.x installed. Make sure you have the OCEAN-modified version\n";
    }
#    move( "ocean.mod", "ppot" ) or die;

    open HFIN, ">hfin1" or die;
    print HFIN "initgrid\n";
    print HFIN "$znucl $grid\n";
    print HFIN "ppload\n";
    print HFIN "mkcorcon\nscreencore\nmkcorcon\ncalcso\nquit\n";
    close HFIN;

    print "Running hfk.x\n";


    system("$ENV{'OCEAN_BIN'}/hfk.x < hfin1 > hfk.${znucl}.1.log") == 0 or die;
    # Check the end of the log to see if we are ok
    my $hfk_status = `tail -n 1 hfk.${znucl}.1.log`;
    unless( $hfk_status =~ m/terminus/ )
    {
      die "The program hfk.x has exited incorrectly for hfin1.\nExiting ...\n";
    }

    my $corezfile = sprintf("corezetaz%03i",$znucl);
    move("xifile","$corezfile");


  # Diagnostics
    my $zdiag = sprintf "zdiagz%3.3i", $znucl;
    print $zdiag . "\n";

    unless( -d "$zdiag" )
    {
      ( mkdir $zdiag ) or die "$!";
    }

    my $zee = sprintf "z%3.3i", $znucl;

  # Wander through the directory making some changes
    opendir DIR, $dir;
    while( my $file = readdir(DIR) )
    {
      if( ( $file =~ m/^[fg]k/ ) or ( $file =~ m/^(ae|ft|ps).z/ ) or
          ( $file =~ m/^(deflin|melfile|corez)/ ) or ( $file =~ m/^(rad|prj)filez/ ) or
          ( $file =~ m/^(vcallel|vvallel|vc_bare|vpseud|valence)/ ) or
          ( $file =~ m/^coreorb/ ) )
      {
        move( $file, "zpawinfo" );
      }
      elsif( $file =~ m/^phr/ )
      {
        my $dest = sprintf( "${file}z%03i", $znucl );
        move( $file, "zpawinfo/$dest" );
      }
      elsif( ( $file =~ m/^melfilez\w$/ ) or ( $file =~ m/^(sm|am|di|pr|psft|aeft)\w$/ ) or
             ( $file =~ m/^(mt|dif)\w\w$/ ) or ( $file =~ m/^(map|ex)/ ) or ( $file =~ /hfin\d/ ) or
             ( $file =~ m/hfk.+log/ ) or ( $file =~ m/aetotal/ ) or ( $file =~ m/radf/ ) )
      {
        move( $file, "$zdiag" );
      }
      elsif( ( $file =~ m/^angs\d$/ ) or ( $file =~ m/^ldep\d$/ ) or ( $file =~ m/^shellR\d\.\d\d$/ ) )
      {
        unlink $file;
      }
      elsif( $file =~ m/${zee}/ )
      {
        move( $file, "$zdiag/" );
      }
    }

    my @paw_files_to_kill = ( 'radpot', 'hapot', 'hfpot', 'config', 'corcon', 'valcon', 'chgsum',
                              'atmf99', 'skip', 'psld', 'rslt', 'tmp', 'aprog', 'vxcofr', 'xifile',
                              'actual', 'ppot' );
    foreach (@paw_files_to_kill)
    {
      unlink $_;
    }

  }


}
######################################
print "OPF section done\n";

open DONE, ">done" or exit 0;
print DONE "1\n";
close DONE;

exit 0;


sub runShirley
{
  my $Z = $_[0];
  my $pspFile = $_[1];
  my $hashRef = $_[2];

#  $pspFile =~ s/\.upf$//i;

  my $targ = catdir( $hashRef->{'psp'}->{'ppdir'}, $pspFile );
  copy( $targ, $pspFile ) == 1 or die "Failed to copy $pspFile\n   $targ\n$!";;
  my $ppmodname = $pspFile . ".mod";

  if( $pspFile =~ m/\.upf$/i ) {
    print "Converting $pspFile with upf2shirley.pl\n";
    system( "$ENV{'OCEAN_BIN'}/upf2shirley.pl $pspFile $ppmodname > upf2shirley.log" ) == 0
      or die "Failed to convert psp file $pspFile\n$!";
  } else {
    print "Converting $pspFile with fhi2eric.x\n";
    system("echo '$pspFile\n$ppmodname' | $ENV{'OCEAN_BIN'}/fhi2eric.x") == 0
        or die "Failed to convert psp file $pspFile\n";
  }

  print "Starting HFK section\n";

  # DELETE in future
  open OUT, ">", "scfac" or die;
  print OUT "1.0\n";
  close OUT;
  
  ### multi-opt hack
  my $optionfilename;
  my @opts = split ' ', $hashRef->{'opf'}->{'shirley'}->{'opts'};
  for( my $i = 0; $i<scalar @opts; $i+=2 )
  {
    print "$opts[$i] $opts[$i+1]\n";
    if( $opts[$i] == $Z )
    {
      $optionfilename = catdir( updir(), $opts[$i+1] );
      last;
    }
  }
  die "Failed to find opts file\n" unless( defined $optionfilename );
  copy( $optionfilename, "atomoptions" ) == 1 or die "Failed to copy $optionfilename to atomoptions\n$!";
  copy( "${pspFile}.mod", "ppot" ) == 1 or die "Failed to copy ${pspFile}.mod to ppot\n$!";

  system( "$ENV{'OCEAN_BIN'}/validate_opts.pl ${pspFile} atomoptions" ) == 0
        or die "Failed to validate options file\nCheck $optionfilename\n"; 

  my $fillFileName;
  my @fill = split ' ', $hashRef->{'opf'}->{'shirley'}->{'fill'};
  for( my $i = 0; $i<scalar @fill; $i+=2 )
  { 
#    print "$fill[$i] $fill[$i+1]\n";
    if( $fill[$i] == $Z )
    { 
      $fillFileName = catdir( updir(), $fill[$i+1] );
      last;
    }
  }
  die "Failed to find fill file\n" unless( defined $fillFileName );
  my $fill;
  {
    open IN, "<", "$fillFileName" or die "Failed to open fill file: $fillFileName\n$!";
    local $/ = undef;
    $fill = <IN>; 
    close IN;
  }
  
  my $grid = $hashRef->{'opf'}->{'shirley'}->{'hfkgrid'}[0] . " " . $hashRef->{'opf'}->{'shirley'}->{'hfkgrid'}[1];


  open HFIN, ">", "hfin1" or die;
  print HFIN "initgrid\n";
  print HFIN "$Z $grid\n";
  print HFIN "ppload\nmkcorcon\nscreencore\nquit\n";
  close HFIN;    

  open HFIN, ">hfin2" or die;
  print HFIN "initgrid\n";
  print HFIN "$Z $grid\n";
  print HFIN "ppload\nmkcorcon\ncalcso\nquit\n";
  close HFIN;

  open HFIN, ">hfin3" or die;
  print HFIN "initgrid\n";
  print HFIN "$Z $grid\n";
  print HFIN "ppload\nmkcorcon\nspartanfip\n";
  print HFIN $fill;
  print HFIN "quit\n";
  close HFIN;

  open OUT, ">", "caution" or die "Failed to open caution\n$!";
  if( $hashRef->{'opf'}->{'shirley'}->{'caution'} )
  {
    print OUT "true\n";
  }
  else
  {
    print OUT "false\n";
  }
#  print OUT $hashRef->{'opf'}->{'shirley'}->{'caution'} . "\n";
  close OUT;

  print "Running hfk.x: $Z\n";

  #change Z -> 0-padded znucl in what follows
  my $znucl = sprintf "%03i", $Z;

  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin1 > hfk.${znucl}.1.log") == 0 or die;
  # Check the end of the log to see if we are ok
  my $hfk_status = `tail -n 1 hfk.${znucl}.1.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly for hfin1.\nExiting ...\n";
  }

  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin2 > hfk.${znucl}.2.log") == 0 or die;
  my $hfk_status = `tail -n 1 hfk.${znucl}.2.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly for hfin2.\nExiting ...\n";
  }

  my $corezfile = sprintf("corezetaz%03i",$Z);
  move("xifile","$corezfile");

  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin3 > hfk.${znucl}.3.log") == 0 or die;
  my $hfk_status = `tail -n 1 hfk.${znucl}.3.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly for hfin3.\nExiting ...\n";
  }

  unless( -d "zdiag_${znucl}" )
  {
    ( mkdir "zdiag_${znucl}" ) or die "$!";
  }
  unless( -d "zpawinfo" )
  {
    mkdir "zpawinfo" or die "$!";
  }

  my $dir = getcwd;

# Wander through the directory making some changes
  opendir DIR, $dir;
  while( my $file = readdir(DIR) )
  {
    if( ( $file =~ m/^[fg]k/ ) or ( $file =~ m/^(ae|ft|ps).z/ ) or
        ( $file =~ m/^(deflin|melfile|corez)/ ) or ( $file =~ m/^(rad|prj)filez/ ) or
        ( $file =~ m/^(vcallel|vvallel|vc_bare|vpseud|valence)/ ) or
        ( $file =~ m/^coreorb/ ) or ( $file =~m/^c2cz/ ) )
    {
      move( $file, "zpawinfo" );
    }
    elsif( $file =~ m/^phr/ )
    {
      my $dest = sprintf( "${file}z%03i", $Z );
      move( $file, "zpawinfo/$dest" );
    }
    elsif( ( $file =~ m/^melfilez\w$/ ) or ( $file =~ m/^(sm|am|di|pr|psft|aeft)\w$/ ) or
           ( $file =~ m/^(mt|dif)\w\w$/ ) or ( $file =~ m/^(map|ex)/ ) or ( $file =~ /hfin\d/ ) or
           ( $file =~ m/hfk.+log/ ) or ( $file =~ m/aetotal/ ) or ( $file =~ m/radf/ ) or 
           ( $file =~ m/z${znucl}l/ ) or ( $file =~m/(scfac|caution|atomoptions)/  ) )
    {
      move( $file, "zdiag_${znucl}" );
    }
    elsif( ( $file =~ m/^angs\d$/ ) or ( $file =~ m/^ldep\d$/ ) or ( $file =~ m/^shellR\d\.\d\d$/ ) )
    {
      unlink $file;
    }
  }

  my @paw_files_to_kill = ( 'radpot', 'hapot', 'hfpot', 'config', 'corcon', 'valcon', 'chgsum',
                            'atmf99', 'skip', 'psld', 'rslt', 'tmp', 'aprog', 'vxcofr', 'xifile',
                            'actual', 'ppot' );
  foreach (@paw_files_to_kill)
  {
    unlink $_;
  }
}

sub runONCV
{
  my $Z = $_[0];
  my $pspFile = $_[1];
  my $hashRef = $_[2];

  
  my $oncvpspInputFile = $pspFile;
  $oncvpspInputFile =~ s/.upf$//i;
  $oncvpspInputFile .= ".in";
  my $targ = catdir( $hashRef->{'psp'}->{'ppdir'}, "$oncvpspInputFile" );

  unless( -e $targ )
  {
    print "Was trying to find input for $pspFile:\t$oncvpspInputFile\n";
    print "$targ\n";
    die;
  }

  copy( $targ, $oncvpspInputFile);

  my $targRad = -1;
  if( exists $hashRef->{'opf'}->{'radius'} ) {
    $targRad = $hashRef->{'opf'}->{'radius'};
  }
  my $scfac = 1.;
  open IN, ">>", $oncvpspInputFile or die "Failed to open $oncvpspInputFile\n$!";
  print IN ".true.\n$targRad   $scfac\n";
  close IN;

  my @oncvpsp;
  open IN, $oncvpspInputFile or die "Failed to open $oncvpspInputFile\n$!";
  while( my $line = <IN> )
  {
    unless( $line =~ m/^\s*#/ )
    {
      chomp $line;
      push @oncvpsp, $line;
    }
  }
  close IN;

  open OUT, ">atomoptions" or die "Failed to open atomoptions for writing\n";

  $oncvpsp[0] =~ m/(\w+)\s+(\d+\.?\d*)\s+(\d+)\s+(\d+)\s+(-?\d+)\s+\w+/
      or die "Failed reading $oncvpspInputFile\t$oncvpsp[0]\n";
  my $zee = $2;
  my $nc = $3;
  my $nv = $4;
  if( $zee != $Z )
  { 
    print "Mismatch between Z in psp.in and Z in input file\n$zee\t$Z\n";
    die;
  }
  print OUT "$zee\n";

  my $lmax = 0;
  my @spdf;
  $spdf[0] = 0; $spdf[1] = 0; $spdf[2] = 0; $spdf[3] = 0;
  for( my $i = 1; $i <= $nc; $i++ )
  { 
    print "CORE: $oncvpsp[$i]\n";
    $oncvpsp[$i] =~ m/^\s*(\d+)\s+(\d+)/ 
        or die "Failed reading core levels of $oncvpspInputFile\t$oncvpsp[$i]\n";
    $spdf[$2]++;
    $lmax = $2 if( $2 > $lmax );
  }

  print OUT "$spdf[0]   $spdf[1]   $spdf[2]   $spdf[3]\n";

  $spdf[0] = 0; $spdf[1] = 0; $spdf[2] = 0; $spdf[3] = 0;
  print OUT "scalar\nrel\nlda\n";
  for( my $i = $nc+1; $i <= $nc+$nv; $i++ )
  {
    print "    : $oncvpsp[$i]\n";
    $oncvpsp[$i] =~ m/^\s*(\d+)\s+(\d+)\s+(\S+)/
        or die "Failed reading core levels of $oncvpspInputFile\t$oncvpsp[$i]\n";
    $spdf[$2] += $3;
    if( $3 > 0.000001 )
    {
      $lmax = $2 if( $2 > $lmax );
    }
  }
  close IN;

  # reset to max per principal
  $spdf[0] = 2  if( $spdf[0] > 2 );
  $spdf[1] = 6  if( $spdf[1] > 6 );
  $spdf[2] = 10 if( $spdf[2] > 10 );
  $spdf[3] = 14 if( $spdf[3] > 14 );


  # This hack is because, currently, hfk.x is used to screen the core
  #  and for some neutral atoms it doesn't work. 
  for( my $i = 0; $i <=3; $i++ )
  {
    $spdf[$i] *= 0.95;
  }

  print OUT "$spdf[0]   $spdf[1]   $spdf[2]   $spdf[3]\n";
  print OUT "$spdf[0]   $spdf[1]   $spdf[2]   $spdf[3]\n";
  close OUT;

  my $log = $pspFile . ".log";
  system("$ENV{'OCEAN_BIN'}/oncvpsp.x < $oncvpspInputFile > $log") == 0 or die;
  unless( move( "ocean.mod", "ppot" ) )
  {
    open IN, "<", $log or die "It appears oncvpsp.x has not run correctly\n";
    my $correctONCVPSP = 0;
    while( my $line = <IN> )
    {
      if( $line =~ m/OCEAN/ )
      {
        $correctONCVPSP = 1;
        last;
      }
    }
    close IN;
    die "oncvpsp.x did not finish correctly. Check the logfile OPF/$log\n" if( $correctONCVPSP );
    die "Incorrect version of oncvpsp.x installed. Make sure you have the OCEAN-modified version\n";
  }
#  move( "ocean.mod", "ppot" ) or die;

  open HFIN, ">hfin1" or die;
  print HFIN "initgrid\n";
  print HFIN "$Z $hashRef->{'opf'}->{'shirley'}->{'hfkgrid'}[0] $hashRef->{'opf'}->{'shirley'}->{'hfkgrid'}[1]\n";
#    print HFIN "fakel\n0 $lmax\n";
  print HFIN "ppload\n";
  print HFIN "mkcorcon\nscreencore\nmkcorcon\ncalcso\nquit\n";
  close HFIN;

  open OUT, ">", "caution" or die "Failed to open caution\n$!";
  if( $hashRef->{'opf'}->{'shirley'}->{'caution'} )
  {
    print OUT "true\n";
  }
  else
  {
    print OUT "false\n";
  }
#  print OUT $hashRef->{'opf'}->{'shirley'}->{'caution'} . "\n";
  close OUT;

  print "Running hfk.x\n";
  system("$ENV{'OCEAN_BIN'}/hfk.x < hfin1 > hfk.${Z}.1.log") == 0 or die;
  # Check the end of the log to see if we are ok
  my $hfk_status = `tail -n 1 hfk.${Z}.1.log`;
  unless( $hfk_status =~ m/terminus/ )
  {
    die "The program hfk.x has exited incorrectly for hfin1.\nExiting ...\n";
  }

  my $corezfile = sprintf("corezetaz%03i",$Z);
  move("xifile","$corezfile");


# Diagnostics
  my $zdiag = sprintf "zdiagz%3.3i", $Z;
  print $zdiag . "\n";

  unless( -d "$zdiag" )
  {
    ( mkdir $zdiag ) or die "$!";
  }
  unless( -d "zpawinfo" )
  {
    mkdir "zpawinfo" or die "$!";
  }

  my $zee = sprintf "z%3.3i", $Z;

  # Wander through the directory making some changes
  opendir DIR, $dir;
  while( my $file = readdir(DIR) )
  { 
    if( ( $file =~ m/^[fg]k/ ) or ( $file =~ m/^(ae|ft|ps).z/ ) or
        ( $file =~ m/^(deflin|melfile|corez)/ ) or ( $file =~ m/^(rad|prj)filez/ ) or
        ( $file =~ m/^(vcallel|vvallel|vc_bare|vpseud|valence)/ ) or
        ( $file =~ m/^coreorb/ ) )
    { 
      move( $file, "zpawinfo" );
    }
    elsif( $file =~ m/^phr/ )
    { 
      my $dest = sprintf( "${file}z%03i", $Z );
      move( $file, "zpawinfo/$dest" );
    }
    elsif( ( $file =~ m/^melfilez\w$/ ) or ( $file =~ m/^(sm|am|di|pr|psft|aeft)\w$/ ) or
           ( $file =~ m/^(mt|dif)\w\w$/ ) or ( $file =~ m/^(map|ex)/ ) or ( $file =~ /hfin\d/ ) or
           ( $file =~ m/hfk.+log/ ) or ( $file =~ m/aetotal/ ) or ( $file =~ m/radf/ ) )
    { 
      move( $file, "$zdiag" );
    }
    elsif( ( $file =~ m/^angs\d$/ ) or ( $file =~ m/^ldep\d$/ ) or ( $file =~ m/^shellR\d\.\d\d$/ ) )
    { 
      unlink $file;
    }
    elsif( $file =~ m/${zee}/ )
    { 
      move( $file, "$zdiag/" );
    }
  }

  my @paw_files_to_kill = ( 'radpot', 'hapot', 'hfpot', 'config', 'corcon', 'valcon', 'chgsum',
                            'atmf99', 'skip', 'psld', 'rslt', 'tmp', 'aprog', 'vxcofr', 'xifile',
                            'actual', 'ppot' );
  foreach (@paw_files_to_kill)
  {
    unlink $_;
  }


}


sub inputHashShirley
{
  my ( $Z, $hashRef ) = @_;

  my $optfill = '';
  ### multi-opt hack
  my $filename;
  my @opts = split ' ', $hashRef->{'opf'}->{'shirley'}->{'opts'};
  for( my $i = 0; $i<scalar @opts; $i+=2 )
  {
    print "$opts[$i] $opts[$i+1]\n";
    if( $opts[$i] == $Z )
    {
      $filename = catdir( updir(), $opts[$i+1] );
      last;
    }
  }
  {
     open my $fh, '<', $filename or die;
     local $/ = undef;
     $optfill .= <$fh>;
     close $fh;
  }
  my @fill = split ' ', $hashRef->{'opf'}->{'shirley'}->{'fill'};
  for( my $i = 0; $i<scalar @fill; $i+=2 )
  {
    print "$fill[$i] $fill[$i+1]\n";
    if( $fill[$i] == $Z )
    {
      $filename = catdir( updir(), $fill[$i+1] );
      last;
    }
  }
  {
     open my $fh, '<', $filename or die;
     local $/ = undef;
     $optfill .= <$fh>;
     close $fh;
  }

  return md5_hex( $optfill );
}


sub inputHashHamann
{    
  my ( $Z, $pspFile, $hashRef ) = @_;

  my $oncvpspInputFile = $pspFile;
  $oncvpspInputFile =~ s/.upf$//i;
  $oncvpspInputFile .= ".in";
  my $filename = catdir( $hashRef->{'psp'}->{'ppdir'}, "$oncvpspInputFile" );

  my $optfill = '';
  {
     open my $fh, '<', $filename or die;
     local $/ = undef;
     $optfill .= <$fh>;
     close $fh;
  }

  return md5_hex( $optfill );
}
