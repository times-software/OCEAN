#!/usr/bin/perl

use strict;

if (! $ENV{"OCEAN_BIN"} ) {
  $0 =~ m/(.*)\/NewAbinitDriver\.pl/;
  $ENV{"OCEAN_BIN"} = $1;
  print "OCEAN_BIN not set. Setting it to $1\n";
}
if (! $ENV{"OCEAN_WORKDIR"}){ $ENV{"OCEAN_WORKDIR"} = `pwd` . "../" ; }
#if (!$ENV{"OCEAN_VERSION"}) {$ENV{"OCEAN_VERSION"} = `cat $ENV{"OCEAN_BIN"}/Version`; }
if (! $ENV{"OCEAN_ABINIT"} ) {$ENV{"OCEAN_ABINIT"} = $ENV{"OCEAN_BIN"} . "/abinis"; }

my @GeneralFiles = ("core", "scratch" );#, "getden");

my @KgenFiles = ("nkpt", "k0.ipt", "qinunitsofbvectors.ipt", "paw.nkpt");
my @BandFiles = ("nbands", "paw.nbands");
my @AbinitFiles = ( "rscale", "rprim", "ntype", "natoms", "typat",
    "verbatim", "coord", "taulist", "ecut", "etol", "nrun", "wftol",
    "fband", "occopt", "ngkpt", "abpad");
my @PPFiles = ("pplist", "znucl");
my @OtherFiles = ("epsilon");

foreach (@GeneralFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@KgenFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}
foreach (@BandFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

foreach (@PPFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

foreach (@AbinitFiles, @OtherFiles) {
  system("cp ../Common/$_ .") == 0 or die;
}

system("$ENV{'OCEAN_BIN'}/pp.pl znucl pplist finalpplist") == 0
    or die "Failed to run pp.pl\n";

my $core = `cat core`;
chomp($core);
`echo 1 > core`;

my $paral_kgb = 1;

open NKPT, "paw.nkpt" or die "Failed to open paw.nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. paw.nkpt\n";
my @pawnkpt = ($1, $2, $3);
close NKPT;
open NKPT, "nkpt" or die "Failed to open nkpt\n";
<NKPT> =~ m/(\d+)\s+(\d+)\s+(\d+)/ or die "Failed to parse. nkpt\n";
my @nkpt = ($1, $2, $3);
close NKPT;
my $pawnbands;
my $nbands;
open NBANDS, "paw.nbands" or die "Failed to open paw.nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse paw.nbands\n";
$pawnbands = $1;
close NBANDS;
open NBANDS, "nbands" or die "Failed to open nbands\n";
<NBANDS> =~ m/(\d+)/ or die "Failed to parse nbands\n";
$nbands = $1;
close NBANDS;

my $scratch = `cat scratch`;
chomp($scratch);

system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0 
  or die "Failed to run kgen2.x\n";
`cp kpts.0001 kpts.bse`;
`cp nkpt bse.nkpt`;
`cp paw.nkpt nkpt`;
system("$ENV{'OCEAN_BIN'}/kgen2.x") == 0
  or die "Failed to run kgen2.x\n";
`cp kpts.0001 kpts.paw`;
`cp bse.nkpt nkpt`;


if ( $pawnkpt[0] + $pawnkpt[1] + $pawnkpt[2] == 0 ) {
  `cp nkpt paw.nkpt`;
  @pawnkpt = [ @nkpt ];
  if ( $pawnbands == 0 ) {
    `cp nbands paw.nbands`;
    $pawnbands = $nbands;
  }
  elsif ( $nbands > $pawnbands ) {
    die "paw.nbands must be larger than nbands\b";
  }
}

my $ndtset = 3;
if ( $nkpt[0] == $pawnkpt[0] && $nkpt[1] == $pawnkpt[1] && $nkpt[2] == $pawnkpt[2] ) {
  print "ndtset = 2\n";
  $ndtset = 2;
  $nbands = $pawnbands;
}

my $kpt_tot = $nkpt[0] * $nkpt[1] * $nkpt[2];
if( $kpt_tot > $core && ( ( $kpt_tot % $core ) != 0 ) ) {
  
  print "$core is not a multiple of nkpts ($kpt_tot)\n";
}
my $paw_kpt_tot = $pawnkpt[0] * $pawnkpt[1] * $pawnkpt[2];
if( $ndtset == 3 && $paw_kpt_tot > $core && ( ( $paw_kpt_tot % $core ) != 0 ) ) {
  print "$core is not a multiple of paw nkpts ($paw_kpt_tot)\n";
}


if( $paral_kgb == 0 ) {
  `echo "$ndtset" > ndtset`;
}

  `echo symmorphi 0 > abfile`;
  `echo -n 'acell ' >> abfile`;
  `cat rscale >> abfile`;
  `echo rprim >> abfile`;
  `cat rprim >> abfile`;
  `echo -n 'ntypat ' >> abfile`;
  `cat ntype >> abfile`;
  `echo -n 'znucl ' >> abfile`;
  `cat znucl >> abfile`;
  `echo -n 'natom ' >> abfile`;
  `cat natoms >> abfile`;
  `echo -n 'typat ' >> abfile`;
  `cat typat >> abfile`;
  `cat coord >> abfile`;
  `cat taulist >> abfile`;
  `echo -n 'ecut ' >> abfile`;
  `cat ecut >> abfile`;
#  `echo -n 'toldfe ' >> abfile`;
#  `cat etol >> abfile`;
  `echo -n 'nstep ' >> abfile`;
  `cat nrun >> abfile`;
  `echo -n 'diemac ' >> abfile`;
  `cat epsilon >> abfile`;
  `cat verbatim >> abfile`;
  `echo -n 'occopt ' >> abfile`;
  `cat occopt >> abfile`;
###
  `cp abfile test_par`;
#  `echo "ndtset 2" >> test_par`;
  
  `echo -n 'fband ' >> test_par`;
  `cat fband >> test_par`;
  `echo prtden 1 >> test_par`;
  `echo kptopt 1 >> test_par`;
  `echo -n 'ngkpt ' >> test_par`;
  `cat ngkpt >> test_par`;
 
  open FILES, ">par.files";
  print FILES "test_par\n"
            . "par.out\n"
            . "SCx\n"
            . "PAR\n"
            . "$scratch\n";
  close FILES;
  `cat finalpplist >> par.files`;
  #
  `echo "paral_kgb -$core" >> test_par`;
  `echo 'wfoptalg 4' >> test_par`;
  `echo 'nloalg 4' >> test_par`;
  `echo 'fftalg 401' >> test_par`;
  `echo 'intxc 0' >> test_par`;
  `echo 'fft_opt_lob 2' >> test_par`;
  `echo 'istwfk *1' >> test_par`;
   
  print "Testing par options\n";
  system("$ENV{'OCEAN_BIN'}/abinis < par.files >& par.log");

  open LOG, "par.log" or die;
  my $logline;
  my $trigger = 0;
  my $nproc = 0;
  my $npkpt;
  my $npband;
  my $npfft;
  my $bandpp;
  my @gs_procs;
  my @temp;
  LOG: while ($logline = <LOG>) {
    if ($trigger == 0 ) {
      if( $logline =~ m/nproc\s+npkpt\s+npband\s+npfft\s+bandpp/ ) {
        $trigger = 1;
      }
    }
    else {
      if ( $logline =~ m/(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
        @temp = ("$1","$2","$3","$4","$5");
        push @gs_procs, [@temp];
        print "$temp[0]\n";
#        next LOG if( $1 < $nproc );
#        $nproc = $1;
#        $npkpt = $2;
#        $npband = $3;
#        $npfft = $4;
#        $bandpp = $5;
      }
      else {
        last LOG;
      }
    }
  }
  close( LOG );
  
  my @paw_procs;

  if( $ndtset == 3 && $core > $paw_kpt_tot) {
    print "ndtset = 3\n";
  `cp abfile test_par`; 
  `echo "nband $pawnbands" >> test_par`;
#  `echo 'iscf -2' >> test_par`;
#  `echo -n 'tolwfr ' >> test_par`;
#  `cat wftol >> test_par`;
#  `echo getden 1 >> test_par`;
  `echo kptopt 0 >> test_par`;
  `sed "s/kpt/kpt/" kpts.paw >> test_par`;
  `echo "paral_kgb -$core" >> test_par`;
  `echo 'wfoptalg 4' >> test_par`;
  `echo 'nloalg 4' >> test_par`;
  `echo 'fftalg 401' >> test_par`;
  `echo 'intxc 0' >> test_par`;
  `echo 'fft_opt_lob 2' >> test_par`;
  `echo 'istwfk *1' >> test_par`;

  system("$ENV{'OCEAN_BIN'}/abinis < par.files >& par3.log");
  open LOG, "par3.log" or die;
#  my @paw_procs;
  $trigger = 0;
  LOG: while ($logline = <LOG>) {
    if ($trigger == 0 ) {
      if( $logline =~ m/nproc\s+npkpt\s+npband\s+npfft\s+bandpp/ ) {
        $trigger = 1;
      }
    }
    else {
      if ( $logline =~ m/(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
        push @paw_procs, [($1,$2,$3,$4,$5)];
#        print "$1,$2,$3,$4,$5\n";
      }
      else {
        last LOG;
      }
    }
  }
  close( LOG );
  } # paw
  else {
    push @paw_procs, [( $core,$core,1,1,1 )];
  }

  my @bse_procs;
  if( $core > $kpt_tot ) {
  `cp abfile test_par`;
  `echo "nband $nbands" >> test_par`;
#  `echo 'iscf -2' >> test_par`;
#  `echo -n 'tolwfr ' >> test_par`;
#  `cat wftol >> test_par`;
#  `echo getden 1 >> test_par`;
  `echo kptopt 0 >> test_par`;
  `cat kpts.bse >> test_par`;
  `echo "paral_kgb -$core" >> test_par`;
  `echo 'wfoptalg 4' >> test_par`;
  `echo 'nloalg 4' >> test_par`;
  `echo 'fftalg 401' >> test_par`;
  `echo 'intxc 0' >> test_par`;
  `echo 'fft_opt_lob 2' >> test_par`;
  `echo 'istwfk *1' >> test_par`;

  system("$ENV{'OCEAN_BIN'}/abinis < par.files >& par2.log");
  open LOG, "par2.log" or die;
#  my @bse_procs;
  $trigger = 0;
  LOG: while ($logline = <LOG>) {
    if ($trigger == 0 ) {
      if( $logline =~ m/nproc\s+npkpt\s+npband\s+npfft\s+bandpp/ ) {
        $trigger = 1;
      }
    }
    else {
      if ( $logline =~ m/(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/ ) {
        push @bse_procs, [($1,$2,$3,$4,$5)];
      }
      else {
        last LOG;
      }
    }
  }
  close( LOG );
  }
  else {
    my $new_core = $core;
    my $cur_try = int( $kpt_tot / $core ) + 1;
    $cur_try = 2 if( $kpt_tot / $core < 2 );
    $cur_try = 1 if( $kpt_tot / $core < 1 ); 
#    my $cur_try = int( $kpt_tot / $core );
    for( ; $cur_try < sqrt($kpt_tot+1); $cur_try ++ ) {
      last if( ( $kpt_tot % $cur_try ) == 0 );
    }
    $new_core = $kpt_tot / $cur_try;
    print "Greatest multiple less than core: $new_core\n";
    push @bse_procs, [( $new_core,$new_core,1,1,1 )];
  }



  my $iter = 0;
  my $iter2;
  my $iter3;
#  print "$#gs_procs\t$#bse_procs\t$#paw_procs\n";
  PROC: for ( $iter = $#gs_procs; $iter >= 0; $iter-- ) {
    for( $iter2 = $#bse_procs; $iter2 >= 0; $iter2-- ) {
      if( $ndtset == 3 ) {
        for( $iter3 = $#paw_procs; $iter3 >= 0; $iter3-- ) {
#          print "$gs_procs[$iter][0]\t$paw_procs[$iter3][0]\t$bse_procs[$iter2][0]\n";
          last PROC if( $gs_procs[$iter][0] == $bse_procs[$iter2][0] && $gs_procs[$iter][0] == $paw_procs[$iter3][0]);
        }
      }
      else {
#        print "$gs_procs[$iter][0]   $bse_procs[$iter2][0]\n";
        last PROC if( $gs_procs[$iter][0] == $bse_procs[$iter2][0] );
      }
    }
  }

#  my @possible_procs;
#  my $good_iter = -1;
  if( $ndtset == 3 ) {
  for( my $match = 0; $match <= 0.25; $match += .05 ) {
    for( $iter3 = $#paw_procs; $iter3 >= 0; $iter3-- ) {
      for( $iter = $#gs_procs; $iter >= 0; $iter-- ) {
        if( ( $gs_procs[$iter][0] <= $paw_procs[$iter3][0] ) &&
            ( $gs_procs[$iter][3] == $paw_procs[$iter3][3] ) )
        {
#          print "$gs_procs[$iter][0] $paw_procs[$iter3][0]\n";
          for( $iter2 = $#bse_procs; $iter2 >= 0; $iter2-- ) {
            if( ( $gs_procs[$iter][3] == $bse_procs[$iter2][3] ) && 
#                ( abs( $bse_procs[$iter2][0] - $paw_procs[$iter3][0] ) < 20 ) ) {
                ( $bse_procs[$iter2][0] >= (1-$match)* $paw_procs[$iter3][0] ) && 
                ( $bse_procs[$iter2][0] <= (1+$match)* $paw_procs[$iter3][0] ) ) {
               # ( abs( $paw_procs[$iter3][0] - $bse_procs[$good_iter][0] ) >= 
               #   abs( $paw_procs[$iter3][0] - $bse_procs[$iter2][0] ) ) ) {
              goto ITER;
              #$good_iter = $iter2;
            }
          }
        }
      } 
    }
  }
  }
  else {
    print "!\n";
    for( $iter = $#gs_procs; $iter >= 0; $iter-- ) {
      for( $iter2 = $#bse_procs; $iter2 >= 0; $iter2-- ) {
        if( ( $gs_procs[$iter][0] <= $bse_procs[$iter3][0] ) &&
            ( $gs_procs[$iter][3] == $bse_procs[$iter3][3] ) )
        {
          goto ITER;
        }
      }
    }
  }
#  if( $good_iter == -1 ) {
#    die "Failed to match procs\n";
#  }
  ITER:
  print "nproc\tnpkpt\tnpband\tnpfft\n";
  print "$gs_procs[$iter][0]\t$gs_procs[$iter][1]\t$gs_procs[$iter][2]\t$gs_procs[$iter][3]\n";
  print "$paw_procs[$iter3][0]\t$paw_procs[$iter3][1]\t$paw_procs[$iter3][2]\t$paw_procs[$iter3][3]\n";
  print "$bse_procs[$iter2][0]\t$bse_procs[$iter2][1]\t$bse_procs[$iter2][2]\t$bse_procs[$iter2][3]\n";
#  die "Failed to find matching proc number\n" if ( $iter == 0 && $iter2 == 0 && $iter3 == 0 &&
#      ( $gs_procs[0][0] != $bse_procs[0][0] ) )
#  print "Use nproc = $gs_procs[$iter][0]\n";
  #
#  $iter = $#gs_procs;
#  $iter2 =  $#bse_procs;
  #
  #
  `cp abfile gs.in`;
  `cp abfile bse.in`;
  `cp abfile paw.in` if ($ndtset == 3);
  #
  if ( $gs_procs[$iter][0] != $gs_procs[$iter][1] ) {
    print "Using KGB Paral methods\n";
    `echo "npkpt $gs_procs[$iter][1]" >> gs.in`;
    `echo "npband $gs_procs[$iter][2]" >> gs.in`;
    `echo "npfft $gs_procs[$iter][3]" >> gs.in`;
#    `echo "bandpp $gs_procs[$iter][4]" >> gs.in`;
    `echo 'paral_kgb 1' >> gs.in`;
    `echo 'wfoptalg 4' >> gs.in`;
    `echo 'nloalg 4' >> gs.in`;
    `echo 'fftalg 401' >> gs.in`;
    `echo 'intxc 0' >> gs.in`;
    `echo 'fft_opt_lob 2' >> gs.in`;
    `echo 'istwfk *1' >> gs.in`;
  }
  else {
    print "Using kpt paral only\n";
  }
  `echo -n 'toldfe ' >> gs.in`;
  `cat etol >> gs.in`;
  `echo -n 'fband ' >> gs.in`;
  `cat fband >> gs.in`;
  `echo prtden 1 >> gs.in`;
  `echo prtwf 0 >> gs.in`;
  `echo kptopt 1 >> gs.in`;
  `echo -n 'ngkpt ' >> gs.in`;
  `cat ngkpt >> gs.in`;
  #  

  if ( $bse_procs[$iter2][0] != $bse_procs[$iter2][1] ) {
    print "Using KGB Paral methods\n";
    `echo "npkpt $bse_procs[$iter2][1]" >> bse.in`;
    `echo "npband $bse_procs[$iter2][2]" >> bse.in`;
    `echo "npfft $bse_procs[$iter2][3]" >> bse.in`;
#    `echo "bandpp $bse_procs[$iter2][4]" >> bse.in`;
    `echo 'paral_kgb 1' >> bse.in`;
    `echo 'wfoptalg 4' >> bse.in`;
    `echo 'nloalg 4' >> bse.in`;
    `echo 'fftalg 401' >> bse.in`;
    `echo 'intxc 0' >> bse.in`;
    `echo 'fft_opt_lob 2' >> bse.in`;
    `echo 'istwfk *1' >> bse.in`;
  }
  else {
    print "Using kpt paral only\n";
  }  
  `echo "nband $nbands" >> bse.in`;
  `echo 'iscf -2' >> bse.in`;
  `echo -n 'tolwfr ' >> bse.in`;
  `cat wftol >> bse.in`;
  `echo getden -1 >> bse.in`;
  `echo prtden 0 >> bse.in`;
  `echo kptopt 0 >> bse.in`;
  `cat kpts.bse >> bse.in`;

  if( $ndtset == 3 ) {
    if ( $paw_procs[$iter3][0] != $paw_procs[$iter3][1] ) {
      print "Using KGB Paral methods\n";
      `echo "npkpt $paw_procs[$iter3][1]" >> paw.in`;
      `echo "npband $paw_procs[$iter3][2]" >> paw.in`;
      `echo "npfft $paw_procs[$iter3][3]" >> paw.in`;
#      `echo "bandpp $paw_procs[$iter2][4]" >> paw.in`;
      `echo 'paral_kgb 1' >> paw.in`;
      `echo 'wfoptalg 4' >> paw.in`;
      `echo 'nloalg 4' >> paw.in`;
      `echo 'fftalg 401' >> paw.in`;
      `echo 'intxc 0' >> paw.in`;
      `echo 'fft_opt_lob 2' >> paw.in`;
      `echo 'istwfk *1' >> paw.in`;
    }
    else {
      print "Using kpt paral only\n";
    }  
    `echo "nband $pawnbands" >> paw.in`;
    `echo 'iscf -2' >> paw.in`;
    `echo -n 'tolwfr ' >> paw.in`;
    `cat wftol >> paw.in`;
    `echo getden -1 >> paw.in`;
    `echo prtden 0 >> paw.in`;
    `echo kptopt 0 >> paw.in`;
    `cat kpts.paw >> paw.in`;
  }
########
  
  `echo "ndtset $ndtset" >> abfile`;

  `echo -n 'fband1 ' >> abfile`;
  `cat fband >> abfile`;
  `echo prtden1 1 >> abfile`;
  `echo kptopt1 1 >> abfile`;
  `echo -n 'ngkpt1 ' >> abfile`;
  `cat ngkpt >> abfile`;

if( $ndtset == 3 ) {
  `echo "nband2 $pawnbands" >> abfile`;
  `echo 'iscf2 -2' >> abfile`;
  `echo -n 'tolwfr2 ' >> abfile`;
  `cat wftol >> abfile`;
  `echo getden2 1 >> abfile`;
  `echo kptopt2 0 >> abfile`;
  `sed "s/kpt/kpt2/" kpts.paw >> abfile`;
}

  `echo "nband$ndtset $nbands" >> abfile`;
  `echo 'iscf$ndtset -2' >> abfile`;
  `echo -n 'tolwfr$ndtset ' >> abfile`;
  `cat wftol >> abfile`;
  `echo getden$ndtset 1 >> abfile`;
  `echo kptopt$ndtset 0 >> abfile`;
  `sed "s/kpt/kpt$ndtset/" kpts.bse >> abfile`;

#  `cat kpts.bse >> abfile`;

open FILES, ">gs.files";
  print FILES "gs.in\n"
            . "gs.out\n"
            . "SCxx\n"
            . "SCx\n"
            . "$scratch\n";
close FILES;
`cat finalpplist >> gs.files`;

open FILES, ">bse.files";
  print FILES "bse.in\n"
            . "bse.out\n"
            . "SCx\n"
            . "RUN0001\n"
            . "$scratch\n";
close FILES;
`cat finalpplist >> bse.files`;

if( $ndtset == 3 ) {
  open FILES, ">paw.files";
    print FILES "paw.in\n"
              . "paw.out\n"
              . "SCx\n"
              . "RUN0002\n"
              . "$scratch\n";
  close FILES;
  `cat finalpplist >> paw.files`;
}

open RUNFILE, ">ab_runfile" or die;
  print RUNFILE "mpirun -np $gs_procs[$iter][0] --hostfile hf " 
              . "$ENV{'OCEAN_BIN'}/abinis < gs.files >& gs.log\n";
if( $ndtset == 3 ) {
  print RUNFILE "mpirun -np $paw_procs[$iter3][0] --hostfile hf " 
              . "$ENV{'OCEAN_BIN'}/abinis < paw.files >& paw.log\n";
}
  print RUNFILE "mpirun -np $bse_procs[$iter2][0] --hostfile hf " 
              . "$ENV{'OCEAN_BIN'}/abinis < bse.files >& bse.log\n";
close RUNFILE;

######
print `pwd`;
open RSCALE, "rscale" or die;
open RPRIM, "rprim" or die;
<RSCALE> =~  m/(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)\s+(\d+\.?\d+([eEfF][+-]?\d+)?)/ or die;
my @rscale = ($1, $3, $5);
print "$1\t$3\t$5\n";
close RSCALE;
  
open AVECS, ">avecsinbohr.ipt" or die;
for (my $i = 0; $i < 3; $i++ ) {
  <RPRIM> =~  m/([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)\s+([+-]?\d?\.?\d+([eEfF][+-]?\d+)?)/ or die "$_";
  print AVECS $1*$rscale[0] . "  " . $3*$rscale[1] .  "  " . $5*$rscale[2] . "\n";
  print "$1\t$3\t$5\n";
} 
close RPRIM; 
close AVECS;

