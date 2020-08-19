use strict;

use JSON::PP;
use File::Spec::Functions;

# line returns and indents
my $json = JSON::PP->new->pretty;
# Sort the keys before dumping the file
my $enable = 1;
$json = $json->canonical([$enable]);

my $outfile;
#my $oncvEXE = "/Users/jtv1/Documents/Code/oncvpsp/src/oncvpsp.x";
my $oncvEXE = $ARGV[0];
my $pspData->{"pseudopotentials"};

for( my $i = 1; $i < scalar @ARGV; $i++ )
{
  my $filename = $ARGV[$i];
  my $data;


  if( open( my $json_stream, $filename ))
  {
        local $/ = undef;
        $data = $json->decode(<$json_stream>);
        close($json_stream);
  }

  #Maken this into a consistency check
  if( defined $outfile ) {
    my $testName = $data->{ "dojo_info" }{"dojo_dir"};
    die "Name convention mismatch between input jsons!\n  $testName  $outfile\n" 
      unless( $testName eq $outfile );
  }
  else {
    $outfile = $data->{ "dojo_info" }{"dojo_dir"};
  }

  my @elements = keys $data->{ "pseudos_metadata" };
  foreach my $element (keys $data->{ "pseudos_metadata" })
  {
    # only do new things
    my $basename = $data->{ "pseudos_metadata" }{ $element }{ "basename" };
    next if(exists $pspData->{"pseudopotentials"}{ "$basename" } );
    

    open OUT, ">", "oncv.inp" or die;
    print OUT $data->{ "pseudos_metadata" }{ $element }{ "input" };
    close OUT;
    my $oncvOut = `$oncvEXE < oncv.inp`;
    print "$oncvOut";
  #  exit 0;

    $pspData->{"pseudopotentials"}{ "$basename" }{ "input" } = $data->{ "pseudos_metadata" }{ $element }{ "input" };
    $oncvOut =~ m/Begin PSPCODE8(.*)END_PSP/s;
    $pspData->{"pseudopotentials"}{ "$basename" }{ "psp8" } = $1;
    $oncvOut =~ m/Begin PSP_UPF(.*)END_PSP/s;
    $pspData->{"pseudopotentials"}{ "$basename" }{ "upf" } = $1;

    last;
  }
}




my $filename = "ONCVinfo";
my $headerData;
if( open( my $json_stream, $filename ))
{
      local $/ = undef;
      $headerData = $json->decode(<$json_stream>);
      close($json_stream);
}
$pspData->{ "psp_info" } = $headerData;

open OUT, ">", "test2.djson";
print OUT $json->encode($pspData);
close OUT;

$outfile .= ".json";
print "$outfile\n";

exit 0;
#my $filename = $data->{ "dojo_info" }{"dojo_dir"} . ".json";
open OUT, ">", "$outfile";
print OUT $json->encode($pspData);
close OUT;
