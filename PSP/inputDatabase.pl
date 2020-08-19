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
use File::Spec::Functions;

my $filename = $ARGV[0];
#print "$ARGV[1]\n" if( scalar @ARGV > 1 ) ;
my $bonusTag = '';
$bonusTag = '-' . $ARGV[1] if( scalar @ARGV > 1 ) ;

#my $filename = "stringent.djson";

my $data;

# line returns and indents
my $json = JSON::PP->new->pretty;

# Sort the keys before dumping the file
my $enable = 1;
$json = $json->canonical([$enable]);

if( open( my $json_stream, $filename ))
{
      local $/ = undef;
      $data = $json->decode(<$json_stream>);
      close($json_stream);
}

print $data->{ "dojo_info" }{ "pp_type" } . "\n";


#open OUT, ">", "test.djson";
#print OUT $json->encode($data);
#close OUT;

my @elements = keys $data->{ "pseudos_metadata" };
foreach my $element (keys $data->{ "pseudos_metadata" })
{
#  print "$element ";
  my $baseName = $data->{ "pseudos_metadata" }{ "$element" }{ "basename" };
  $baseName =~ s/\.psp8/.in/;
  
  my $inputFile = catfile( "$element", "$baseName" );
#  print $inputFile . "\n";;

  my $psp;
  if( open IN, "$inputFile" )
  {
    local $/ = undef;
    $psp = <IN>;
    close( IN )
  }
  
  #Change psp8 or upf to both
  $psp =~ s/psp8/both/;
  $psp =~ s/upf/both/;

  $data->{ "pseudos_metadata" }{ "$element" }{ "input" } = "$psp";
  $data->{ "pseudos_metadata" }{ "$element" }{ "basename" } = $baseName;
  foreach my $i ( "md5", "dfact_meV", "dfactprime_meV", "tags" )
  {
    delete $data->{ "pseudos_metadata" }{ "$element" }{ "$i" };
  }
}


$filename = "Header";
my $headerData;
if( open( my $json_stream, $filename ))
{
      local $/ = undef;
      $headerData = $json->decode(<$json_stream>);
      close($json_stream);
}
$data->{ "dojo_info" } = $headerData->{ "dojo_info" };

my $filename = $data->{ "dojo_info" }{"dojo_dir"} . "$bonusTag" . ".json";
open OUT, ">", "$filename";
print OUT $json->encode($data);
close OUT;

if( 0 )
{
$filename = "ONCVinfo";
my $headerData;
if( open( my $json_stream, $filename ))
{
      local $/ = undef;
      $headerData = $json->decode(<$json_stream>);
      close($json_stream);
}
$data->{ "psp_info" } = $headerData;

open OUT, ">", "test2.djson";
print OUT $json->encode($data);
close OUT;
}
