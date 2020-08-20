#!/usr/bin/perl
# Copyright (C) 2020 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
use strict;

use JSON::PP;
#use File::Spec::Functions;
use Compress::Zlib;
use MIME::Base64;

# line returns and indents
my $json = JSON::PP->new->pretty;
# Sort the keys before dumping the file
my $enable = 1;
$json = $json->canonical([$enable]);

my $outroot;
#my $oncvEXE = "/Users/jtv1/Documents/Code/oncvpsp/src/oncvpsp.x";
my $oncvEXE = $ARGV[0];
my $upfData->{"pseudopotentials"};
my $psp8Data->{"pseudopotentials"};

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
  if( defined $outroot ) {
    my $testName = $data->{ "dojo_info" }{"dojo_dir"};
    die "Name convention mismatch between input jsons!\n  $testName  $outroot\n" 
      unless( $testName eq $outroot );
  }
  else {
    $outroot = $data->{ "dojo_info" }{"dojo_dir"};
  }

  my @elements = keys $data->{ "pseudos_metadata" };
  foreach my $element (keys $data->{ "pseudos_metadata" })
  {
    # only do new things
    my $basename = $data->{ "pseudos_metadata" }{ $element }{ "basename" };
    next if(exists $upfData->{"pseudopotentials"}{ "$basename" } );
    

    open OUT, ">", "oncv.inp" or die;
    print OUT $data->{ "pseudos_metadata" }{ $element }{ "input" };
    close OUT;
    my $oncvOut = `$oncvEXE < oncv.inp`;
    print "$oncvOut";
  #  exit 0;

    $upfData->{"pseudopotentials"}{ "$basename" }{ "input" } = $data->{ "pseudos_metadata" }{ $element }{ "input" };
    $psp8Data->{"pseudopotentials"}{ "$basename" }{ "input" } = $data->{ "pseudos_metadata" }{ $element }{ "input" };

    $oncvOut =~ m/Begin PSPCODE8\s*\n*(.*)END_PSP.*Begin PSP_UPF/s;
    my $temp = $1;
    chomp( $temp );
    my $compressed = compress( $temp );
    $psp8Data->{"pseudopotentials"}{ "$basename" }{ "psp8" } = encode_base64( $compressed );
    $oncvOut =~ m/Begin PSP_UPF\s*\n*(.*)END_PSP/s;
    my $temp = $1;
    chomp( $temp );
    my $compressed = compress( $temp );
    $upfData->{"pseudopotentials"}{ "$basename" }{ "upf" } = encode_base64( $compressed );

#    last;
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
$upfData->{ "psp_info" } = $headerData;
$psp8Data->{ "psp_info" } = $headerData;

my $outfile;
$outfile = $outroot . ".upf.json";
print "$outfile\n";
open OUT, ">", "$outfile";
print OUT $json->encode($upfData);
close OUT;

$outfile = $outroot . ".psp8.json";
print "$outfile\n";
open OUT, ">", "$outfile";
print OUT $json->encode($psp8Data);
close OUT;

