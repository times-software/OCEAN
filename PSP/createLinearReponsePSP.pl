use strict;

#die "Usage [PSP] [Core hole] [Strength]\n" if( scalar @ARGV <3 );
if( scalar @ARGV < 3 ) {
  print "Usage [PSP] [Core hole] [Strength]\n" .
        " For example:   [PSP] = C.UPF\n" .
        "                [Core hole] = vc_barez006n01l00\n" .
        "                [Strength] = 0.01 (for 1%)\n";
  die;
}
my $pspIn = $ARGV[0];
my $coreHole = $ARGV[1];
my $strength = $ARGV[2];

open PSPIN, "<", $pspIn or die "Failed to open $pspIn\n$!";
open CH, "<", $coreHole or die "Failed to open $coreHole\n$!";

my $pspOut = $pspIn;
my $UPF = 0;
if( $pspOut =~ m/upf$/i ) {
  $UPF = 1;
  $pspOut =~ s/\.upf//i;
  $pspOut .= '-mod.UPF';
} else {
  die "Only UPF supported currently\n";
}
open PSPOUT, ">", $pspOut or die "Failed to open $pspOut\n$!";

my @CHPot;
my @CHRad;

while( my $line = <CH> ) {
  $line =~ m/(\S+)\s+(\S+)/ or die "Bad line in $coreHole\n";
  push @CHRad, $1;
  push @CHPot, $2*$strength*2.0;  # Ha to Ryd
}
close CH;
my $CHnum = scalar @CHRad;
$CHnum -= 1;
#print "$CHRad[0] $CHRad[-1]  $CHnum\n";

my @PPRad;
while( my $line = <PSPIN> ) {
  if( $line =~ m/<PP_R / ) {
    my $radii;
    print PSPOUT $line;
    do {
      $line = <PSPIN>;
      $radii .= $line unless( $line =~ m/PP_R/);
      print PSPOUT $line;
    } until( $line =~ m/PP_R/ );
    @PPRad = split ' ', $radii;
#    foreach ( @PPRad ) { print $_ . "\n"; };
#    print scalar @PPRad . "\n";
    
  } elsif( $line =~ m/<PP_LOCAL/ ) {
    my $i = 0;
    my $j = 0;
    my $col = 4;
    if( $line =~ m/columns="(\d+)"/ ) {
      $col = $1;
    }
    print PSPOUT $line;
    do {
      $line = <PSPIN>;
      unless( $line =~ m/PP_LOCAL/ ) {
        my @TPot = split ' ', $line;
        foreach my $Pot (@TPot) {
          my $V;
          if( $PPRad[$i] <= $CHRad[0] ) {
            $V = $CHPot[0];
          } elsif( $PPRad[$i] > $CHRad[$CHnum] ) {
            print "$PPRad[$i] $CHRad[$CHnum]\n";
            $V = $CHPot[$CHnum];
          } else {
            while( $PPRad[$i] > $CHRad[$j] ) { $j++; }
            if( $j == 0 ) {
              $V = $CHPot[0];
            } else {
              my $b = ( $CHPot[$j] - $CHPot[$j-1] ) / ($CHRad[$j] - $CHRad[$j-1] );
              $V = ($PPRad[$i] - $CHRad[$j] ) * $b + $CHPot[$j];
            }
          }
#          printf "%.10E  %.10E  %.10E  $i  $j\n", $PPRad[$i], $Pot, $Pot+$V;
          printf  PSPOUT "%20.10E", $Pot+$V;
          $i++;
        }
        print PSPOUT "\n";
      }
    } until ( $line =~ m/PP_LOCAL/ );
    print PSPOUT $line;

  } elsif( $line =~ m/z_valence/i || $line =~ m/z valence/i ) {
    my $n = 0;
    if( $line =~ m/(\d+.?\d*[Ee]?\+?\d*)/ )
    {
      $n = $1;
    }
    if( $n == 0 ) { die; }
    $n = $n + $strength;
    printf PSPOUT "       z_valence=\" %.4f\"\n", $n;
  } else {
    print PSPOUT $line;
  }
}
close PSPIN;
close PSPOUT;

print "Done\n";
