use strict;

my @avec;
my @alen;

my $loud = 1;
my $kden = 0;
my $xden = 0;

my $b2a = 1/.529177210903;

if( scalar @ARGV > 0 ) {

  my $input_filename = $ARGV[0];

  my $loud = 1;

  open IN, $input_filename or die "Failed to open $input_filename\n$!\n";

  my $file;
  while (my $line = <IN> )
  {
    chomp($line);
    # if there are comment characters -- #, *, or ! --
        #   remove them and everything following
    $line =~ s/[#\*\\!].*/ /;

    $file .= ' ' . $line;
  }
  close IN;

  $file =~ m/acell\s+\{\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/i or die "Failed to grab acell from the input $ARGV[0]\n";
  my @acell = ( $1, $2, $3 );

  $file =~ m/rprim\s+\{\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/i or die "Failed to parse rprim\n";
  $avec[0][0] = $1; $avec[0][1] = $2; $avec[0][2] = $3;
  $avec[1][0] = $4; $avec[1][1] = $5; $avec[1][2] = $6;
  $avec[2][0] = $7; $avec[2][1] = $8; $avec[2][2] = $9;

  for( my $i = 0; $i<3; $i++ ) {
    $avec[$i][0] *= $acell[$i];
    $avec[$i][1] *= $acell[$i];
    $avec[$i][2] *= $acell[$i];

    $alen[$i] = sqrt( $avec[$i][0]**2 + $avec[$i][1]**2 + $avec[$i][2]**2 )
  }
}
else {
  open AVEC, "avecsinbohr.ipt" or die "Failed to open avecsinbohr.ipt\n$!";
  for( my $i = 0; $i < 3; $i ++ )
  {
    <AVEC> =~ m/(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)/ or die "Couldn't parse avecsinbohr.ipt$_\n";
    $avec[$i][0] = $1;
    $avec[$i][1] = $2;
    $avec[$i][2] = $3;
    $alen[$i] = sqrt( $avec[$i][0]**2 + $avec[$i][1]**2 + $avec[$i][2]**2 )
  }
  close AVEC;
}

if( $loud == 1 ) {
  print "A-vector lengths\n";
  printf "%10.4f %10.4f %10.4f\n", $alen[0], $alen[1], $alen[2];
}


my $vol = 0;
for( my $i = 0; $i < 3; $i++ ) {
  $vol += $avec[0][$i] * ( $avec[1][$i-2] * $avec[2][$i-1] - $avec[2][$i-2] * $avec[1][$i-1]);
}

$vol = abs( $vol );
if( $loud == 1 ) {
#  printf "Vol: %10.4f\n", $vol;
}

my @bvec;
my @blen;

my $pref = 2*4*atan2(1,1)/$vol;

#print "B-vectors\n";
for (my $i = 0; $i < 3; $i++ ) {
  for (my $j = 0; $j < 3; $j++ ) {
    $bvec[$i][$j] = $pref * ($avec[$i-2][$j-2] * $avec[$i-1][$j-1] - $avec[$i-2][$j-1] * $avec[$i-1][$j-2]) ;
#    printf "%10.4f ", $bvec[$i][$j] if ( $loud == 1 );
  }
#  print "\n" if( $loud == 1 );
  $blen[$i] = sqrt( $bvec[$i][0]**2 + $bvec[$i][1]**2 + $bvec[$i][2]**2 );
}

if( $loud == 1 ) {
  print "B-vector lengths\n";
  printf "%10.4f %10.4f %10.4f\n", $blen[0], $blen[1], $blen[2];
  print "#######################\n";
}

if( $loud == 1 ) {
  my $klen = 1/$blen[0];
  #for (my $i = 1; $i <= 62; $i++ )
  my @kpt = (1,1,1);
  for( my $j = 1; $j<3; $j++ ) {
    my $t = $kpt[$j]/$blen[$j];
    $klen = $t if( $t<$klen);
  }
  my $klenPrevious = $klen;
  print " kx   ky   kz          kden   (inv Bohr)      (inv A)  radius (Bohr)\n";
  printf "%3i  %3i  %3i  %12.5f %12.5f %12.5f %12.5f\n", $kpt[0], $kpt[1], $kpt[2], $klen, 1/$klen,  $b2a/$klen, $klen*2*4*atan2(1,1);
  while( $klen < 20 )
  {
    my $t = $kpt[0]/$blen[0] * 1.01;
#    print "0 $t\n";
    for( my $j=1; $j<3; $j++ ) {
      $t = $kpt[$j]/$blen[$j] * 1.01 if ( $t > $kpt[$j]/$blen[$j] * 1.01 );
#      print "$j $t\n";
    }
    
  if( 0 ) {
    my $t = 10*$klen;
    for (my $j = 0; $j < 3; $j++ )
    {
      my $i = $kpt[$j]+1;
      while( isAllowed($i) == 0 ) {
        $i++;
      }
      $t = $i/$blen[$j] if( $t>$i/$blen[$j] );
#      print "$t $j $i $blen[$j]  " . ($i/$blen[$j]) . "\n";
    }
    $t *= 0.9999;
    }
    $klen *= 10;
    for (my $j = 0; $j < 3; $j++ )
    {
      $kpt[$j] = int( $t * $blen[$j] )+1;
      while( isAllowed( $kpt[$j] ) == 0 )
      {
        $kpt[$j]++;
      }
      $klen = $kpt[$j]/$blen[$j] if( $kpt[$j]/$blen[$j]<$klen);
    }
#    if( $klen - $klenPrevious > 0.001 ) {
      printf "%3i  %3i  %3i  %12.5f %12.5f %12.5f %12.5f\n", $kpt[0], $kpt[1], $kpt[2], $klen, 1/$klen,  $b2a/$klen, $klen*2*4*atan2(1,1);
#    }
    $klenPrevious = $klen;
  }

  print "#######################\n";

  my $xlen = 0;
  my @xpt = (1,1,1);
  for( my $j = 0; $j<3; $j++ ) {
    my $t = $xpt[$j]/$alen[$j];
    $xlen = $t if( $t>$xlen);
  }
  my $xlenPrevious = $xlen;
  print "  x    y    z          xden       (Bohr)          (A)\n";
  printf "%3i  %3i  %3i  %12.5f %12.5f %12.5f\n", $xpt[0], $xpt[1], $xpt[2], $xlen, 1/$xlen,  1/($b2a*$xlen);
  while( $xlen < 5.5 )
  {
    my $t = 10*$xlen;
    for (my $j = 0; $j < 3; $j++ )
    {
      my $i = $xpt[$j]+1;
      while( isAllowed($i) == 0 ) {
        $i++;
      }
      $t = $i/$alen[$j] if( $t>$i/$alen[$j] );
    }
    $t *= 0.9999;
    for (my $j = 0; $j < 3; $j++ )
    {
      $xpt[$j] = int( $t * $alen[$j] )+1;
      while( isAllowed( $xpt[$j] ) == 0 )
      {
        $xpt[$j]++;
      }
      $xlen = $xpt[$j]/$alen[$j] if( $xpt[$j]/$alen[$j]>$xlen);
    }
    if( $xlen - $xlenPrevious > 0.001 ) {
      printf "%3i  %3i  %3i  %12.5f %12.5f %12.5f\n", $xpt[0], $xpt[1], $xpt[2], $xlen, 1/$xlen,  1/($b2a*$xlen);
    }
    $xlenPrevious = $xlen;
  }
}

exit 0;


sub isAllowed {
  my $n = $_[0];

  for( my $i=2; $i <= 7; $i++ )
  {
    while(  $n % $i == 0 )
    {
      $n /= $i;
    }
  }
  if( $n == 1 || $n == 11 || $n == 13 )
  {
    return 1;
  }
  return 0;
}


