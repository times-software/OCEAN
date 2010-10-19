#!/usr/bin/perl
# This program is passed the pp_file names from the common script
# and if it can't find them then it looks for them using the old, 
# dumb way of searching for the appropriate nucleus guys.

# The purpose of this is to create the pp list for the abinit runfile
# based either on a list given in the input file (and passed to this
# program or on the list of nuclei picking one that matches and throwing an
# error to stderr or the log.

use strict;
use FindBin qw($Bin);
use Cwd;
use Cwd 'abs_path';
use File::Basename;

my $WORKDIR = getcwd;
my $AI2NBSE_BIN = $Bin;
my $AI2NBSE_PP = $AI2NBSE_BIN . "/../PseudoPots";
$AI2NBSE_PP = abs_path($AI2NBSE_PP);
print "$AI2NBSE_PP\n";
$WORKDIR = abs_path("$WORKDIR/..");
#print STDERR "$WORKDIR\t$AI2NBSE_BIN\t$AI2NBSE_PP\n";

if ($#ARGV != 2) {
  die "Expecting 3 inputs. List of Z, PP names, PP output file.";
}

my @Z;
my @PP;
my $Zn_file = "";
my $PP_file = "";

open(Zn,<$ARGV[0]>) or die "File open failed\n";
open(PP,<$ARGV[1]>) or die "File open failed\n";
while (<Zn>) {
  chomp();
  $Zn_file .= $_ . " ";
}
while (<PP>) {
  chomp();
  $PP_file .= $_ . " ";
}
close(Zn);
close(PP);

@Z  = split(' ', $Zn_file);
@PP = split(' ', $PP_file); 

# at this point can be integrated as subroutine with two arrays
# or better as two different subs each taking one array.
open(OUTFILE,">$ARGV[2]");
my $pp_file;
if ($PP[0] ne "NULL") {
  foreach $pp_file (@PP) {
#    next if ($pp_file =~ m/^\s*$/ );
#    chomp($pp_file);
    if (-e "$WORKDIR/$pp_file") {
      print OUTFILE "$WORKDIR/$pp_file\n";
    }
    else {
#      print "Didn't find $WORKDIR/$pp_file\n";
      if (-e "$AI2NBSE_PP/$pp_file") {
        print OUTFILE "$AI2NBSE_PP/$pp_file\n";
      }
      else {
        die "FAILED to find PP file $pp_file. Quitting...\n";
      }
    }
  } 
}
else {
  my @files;
  my $zfile;
  my @files2;
  my $curZ;
  @files  = <$WORKDIR/*>;
  @files2 = <$AI2NBSE_PP/*>;
  print STDERR "Warning: PP files unspecified! I may choose unexpected ones\n";
  foreach $curZ (@Z) {
    foreach $zfile (@files) {
      $zfile = basename($zfile);
#      print  "$zfile\n";
      if ($zfile =~ m/^(\d+)/) {
#       print $1 . "\n";
       if ($1 == $curZ ) {
        $pp_file = $WORKDIR . "/" . $zfile;
        goto ZFOUND;
       }
      }
    }
    foreach $zfile (@files2) {
      $zfile = basename($zfile);
      if ($zfile =~ m/^(\d+)/) {
       if ($1 == $curZ) {
#       if ($_ =~ m/($curZ\S+\.pspnc)/ ) {
        $pp_file = $AI2NBSE_PP . "/" . $zfile;
        goto ZFOUND;
       }
      }
    }
    die "PP file for $curZ not found. Quiting...\n";
    ZFOUND : {
     print OUTFILE "$pp_file\n";
     print "$pp_file\n";
    }
  }
}
close(OUTFILE);
exit 0;
