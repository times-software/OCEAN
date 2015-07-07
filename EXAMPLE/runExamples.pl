#!/usr/bin/perl/

use strict;
use POSIX;

my $PPN = 12;
my $qsub = "qsub";
my $queuename = "-q quick";
my $node_extra = ":sandy";
my $env_load = "source ~/bin/intel/15";
my $sleep = 5;
my $instdir = "/home/jtv1/bin/intel/TEST/";
my $scratchdir = "/wrk/jtv1/TEST/";

my @DFT_flavors = ("qe"); #"abi", "qe", "obf" );

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
#exit 0;
  sleep 60;
}







