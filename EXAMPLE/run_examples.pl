#!/usr/bin/perl

use strict;
use POSIX;

my $PPN = 12;
my $qsub = "qsub";
my $batch_script = "pbs";
#my $batch_script = "sbatch";
my $queuename = "";
my $env_load = "";
my $node_extra = "";
my $accounting = "";
my $sleep = 1;
my $instdir = "/path/to/ocean";
my $scratchdir = "/scratch/OCEAN_TEST";

my @exclude_list = ("AN");

print "$scratchdir\n";
`mkdir -p $scratchdir`;

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

my $DFT_flavor;
if( exists $allowed_dft_flavors{ lc($ARGV[0]) } )
{
  $DFT_flavor = lc($ARGV[0]);
}
else
{
  die "Requested DFT type $ARGV[0] does not match known types\n";
}

my @runlist;
my $pwd = `pwd`;
chomp($pwd);
print "$pwd\n";

my @example_dir;

if( scalar @ARGV == 2 )
{
  print "$ARGV[1]\n";
  die "Requested example not found\n" unless -d $ARGV[1];
  push @example_dir, $ARGV[1];
}
else
{
  # Run through each directory 
 
  opendir DIR, "." ;
  while( my $example_dir = readdir DIR )
  {
  # skip .
#  next if $example_dir =~ m/^\.\n/;
    next unless -d $example_dir;
    next if $example_dir eq '.' or $example_dir eq '..';
  
    chomp($example_dir);
    my $skip = 0;
    foreach my $exclude (@exclude_list)
    {
      if( $example_dir =~ m/^$exclude$/ )
      {
        $skip = 1;
        last;
      }
    }
    next if ( $skip == 1 );
    push @example_dir, $example_dir;

    print "$example_dir\n";
  }
}
sleep $sleep;

foreach my $example_dir (@example_dir)
{

# Run through each DFT flavor
  my $run_dir = $example_dir . "_" . $DFT_flavor;

  # make a copy at our scratch dir
  `cp -r $example_dir $scratchdir/$run_dir`;
  chdir "$scratchdir/$run_dir" or die "Failed to chdir: $!\n";

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
  open OUT, ">runit" or die "$!";
  if( $batch_script eq 'sbatch' )
  {
    print OUT "#!/bin/bash\n";
    print OUT "#SBATCH -N $nnodes ${node_extra} $queuename -J $run_dir $accounting\n";
  }
  else
  {
    print OUT "#PBS -l nodes=$nnodes:ppn=${PPN}${node_extra} $queuename -N $run_dir $accounting\n";
  }
  print OUT "$env_load\n\n";
  print OUT "time $instdir/ocean.pl $file_name\n";
  close OUT;

  #start run and store number to array runlist
  push @runlist, `$qsub runit`;
  print $runlist[0] . "\n";

  chdir "$pwd";


  sleep $sleep;
}







