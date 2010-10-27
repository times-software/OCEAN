#!/usr/bin/perl

use strict;
#use Findbin '$Bin';

$| = 1;
my $OCEAN_ABINIT;
unless ($ENV{'OCEAN_ABINIT'}) {
#  $OCEAN_ABINIT = "$Bin/abinit";
 die unless ($0 =~ m/(\S+)par_ab2.pl/ );
 $OCEAN_ABINIT = "$1/abinit";
}
else {
  $OCEAN_ABINIT = $ENV{'OCEAN_ABINIT'};
}
print "$OCEAN_ABINIT\n";

my $nfiles;
open(NFILES,"Nfiles") or die "Failed to open Nfiles\n";
if (<NFILES> =~ m/(\d+)/ ) {
  $nfiles = $1;
}
else {
  die "Nfiles fail\n";
}
close NFILES;

open(CORE,"core") or die "Failed to open core\n";
my $OCEAN_CORE;
if (<CORE> =~ m/(\d+)/ ) {
  $OCEAN_CORE = $1;
}
else {
  die "Failed to parse core\n";
}

my @proc;
my $file_counter = 0;
while ($file_counter < $nfiles ) {
  foreach (@proc) {
    waitpid($_,0);
    if ($? != 0) { die "Failed\n"; }
  }
  @proc =();
  for (my $core = 1; $core <= $OCEAN_CORE; $core++ ) {
   $file_counter++;
   if ($file_counter > $nfiles) {
      exit 0;
   }
   my $pid = fork();
   if ( $pid ) {
    push(@proc, $pid);
   }
   elsif ($pid == 0 ) {
     print "Starting ABINIT run $file_counter out of $nfiles\n";
     my $denin = sprintf("denin.files.%04i", $file_counter);

     if (system("$OCEAN_ABINIT < $denin >& ABINIT.$file_counter.log") != 0 ) {
       print STDERR  "Failed to run Abinit $file_counter\n";
       exit 1;
     }
     exit 0;
   }
   else {
     die "Failed to fork process for abinit\n";
   }
  }
}
  foreach (@proc) {
    waitpid($_,0);
    if ($? != 0) { die "Failed\n"; }
  }
  @proc =();

exit 0;
