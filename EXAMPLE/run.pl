#!/usr/bin/perl/

use strict;
use POSIX;

my $PPN = 32;
my $qsub = "qsub";
my $queuename = "";#-q vinson";
my $node_extra = ":vinson";
my $env_load = "source ~/.bashrc\nmodule load intel\nmodule load openmpi/2.0.0\nmodule load fftw"; #source ~/bin/intel/15";
my $sleep = 1;
my $instdir = "/users/jtv1/cluster/bin/ocean/ocean2x";
my $scratchdir = "/flash/jtv1/ocean_test/ocean2x/";

my @exclude_list = ("AN", "Cu/HighQuality", "diamond");

#my @DFT_flavors = ("abi", "qe", "obf" );
my %allowed_dft_flavors =  (
  "abi" => '1',
  "qe"  => '1',
  "obf" => '1' 
);

if( scalar(@ARGV) < 1 )
{
  print "List DFT flavors after script\n";
  exit 0;
}
my @DFT_flavors;
foreach my $test_dft (@ARGV)
{
  push @DFT_flavors, lc($test_dft) if exists $allowed_dft_flavors{ lc($test_dft) };
}

my @runlist;
my $pwd = `pwd`;
chomp($pwd);
print "$pwd\n";

#Open the current directory
my $fh;
open( $fh, "-|", "find", ".", "-type", "d" );

# Run through each directory 
while( my $example_dir = <$fh> )
{
  # skip .
  next if $example_dir =~ m/^\.\n/;
  chomp($example_dir);
  my $skip = 0;
  foreach my $exclude (@exclude_list)
  {
    if( $example_dir =~ m/^\.\/$exclude$/ )
    {
      $skip = 1;
      last;
    }
  }
  next if ( $skip == 1 );

  print "$example_dir\n";

  # Run through each DFT flavor
  foreach my $DFT_flavor (@DFT_flavors)
  {
    my $run_dir = $example_dir . "_" . $DFT_flavor;

    # make a copy at our scratch dir
    `cp -r $example_dir $scratchdir/$run_dir`;
    chdir "$scratchdir/$run_dir";

    my $file_name = "$example_dir" . ".in";
    print "$file_name\n";

    # add the DFT flavor to the input
    `echo "dft $DFT_flavor" >> $file_name`;

    #figure out how many nodes to ask for
    open IN, "$file_name" or die "$!";
    my $line;
    while( $line = <IN>) 
    {
      last if $line =~ m/para_prefix/;
    }
    close IN;
    $line =~ m/(\d+)/ or die "$line\n$!";
    my $nprocs = $1;

    my $nnodes = floor( ( $nprocs - 1  ) / $PPN ) + 1;
    $nnodes = sprintf("%.0f", $nnodes );

    # build pbs file
    open OUT, ">runit" or die;
    print OUT "#PBS -l nodes=$nnodes:ppn=${PPN}${node_extra} $queuename -N $run_dir\n";
    print OUT "$env_load\ncd \$PBS_O_WORKDIR\n\n";
    print OUT "time $instdir/ocean.pl $file_name\n";
    close OUT;

    #start run and store number to array runlist
    push @runlist, `$qsub runit`;

    chdir "$pwd";


  }
  sleep $sleep;
}







