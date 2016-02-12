#!/usr/bin/perl
# Copyright (C) 2015 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
# This perl script attempts to append the hash of the current git commit to 
#  the version number. If .git doesn't exist, or if the git command returns
#  something that is not a 40 char hash then nothing is appended.

use strict;


open IN, "VersionNumber" or die "Failed to open VersionNumber!\n$!";
my $version = <IN>;
close IN;

open OUT, ">Version" or die "Failed to open Version!\n$!";
my $stub = '';

if( -d ".git" )
{
  $stub = `git rev-parse HEAD`;
  print $stub;
  chomp $stub;
  # Check that stub is a hash (lowercase alpha & numeric only!
  unless( $stub =~ m/^[a-z0-9]+$/ )
  {
    $stub = '';
  }
}

#print (length $stub) ;
if( length $stub == 40 )
{
  $version .= $stub . "\n";
}

print OUT $version;

close OUT;
