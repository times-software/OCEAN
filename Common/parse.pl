# Copyright (C) 2013 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
#!/usr/bin/perl

use strict;

my $option = $ARGV[0];
my $input_filename = $ARGV[1];
unless(  -e $input_filename )
{
	die "Could not find file $input_filename!\n";
}


my $config_filename = $ARGV[2];
unless( -e $config_filename )
{
	die "Could not find configuration file $config_filename!\n";
}

open IN, $config_filename or die "Failed to open config file: $config_filename\n$!\n";

my @required_list = ();
my %optional_hash = ();
my %filename_hash = ();

while ( my $line = <IN> )
{
	next unless( $line =~ m/\S/ );
	chomp($line);
	my $keyd = $line;
	$line = <IN>;
	$line = <IN>;
	if( $line =~ m/\}\s+(\S+)/ )
	{
		$line =~ m/\}\s+(\S+)/;
		$filename_hash{ lc $keyd } = "$1";
		push @required_list, $keyd;
	}
	else
	{
		my $tmp_buffer = $line;
		while( $line = <IN> )
		{
			last if( $line =~ m/\}\s*(\S+)/ );
			$tmp_buffer .= $line;
		}
#		$tmp_buffer .= $_ until ( <IN> =~ m/\}\s*(\S+)/ );
			
		$optional_hash{ lc $keyd } = "$tmp_buffer";
#		$line = <IN>;
		$line =~ m/\}\s*(\S+)/ or die;
		$filename_hash{ lc $keyd } = "$1";
#		print "$keyd, $1, $tmp_buffer";
	}
}

open DERP, ">optional.txt";
for( keys %optional_hash )
{
	print DERP "$_ >>> " . $optional_hash{ $_ } . "\n";
}
close DERP;


close IN;

print "Config read in\n";



open IN, $input_filename or die "Failed to open $input_filename\n$!\n";

my $tmp_buffer;
my $have_brace = 0;
my $have_continue = 0;

my $case_fix;
my %input_hash = ();

open TEST, ">test.txt";

while (my $line = <IN> )
{
	# remove endline
	chomp($line);
	
	# if there are comment characters -- #, *, or ! --
        #   remove them and everything following
	$line =~ s/[#\*\\!].*/ /;

	# remove blank lines
	next unless( $line =~ m/\S/ );

	if( $have_continue )
	{
		if( $have_brace )
		{
			if( $line =~ m/\{/ )
                        {
				die "Found bonus {\n";
			}
			elsif( $line =~ m/\}/ )
			{
				$have_continue = 0;
				$have_brace = 0;
				$line =~ s/\}/ /;
				$tmp_buffer .= " " . $line;
				$tmp_buffer =~ s/\s+/ /g;
				chomp($tmp_buffer);
#				print TEST "$tmp_buffer\n";
				$tmp_buffer =~ m/^\s*(\S+)\s+(.+)/ or die "$tmp_buffer\nFailed to parse.\n";
				$input_hash{ lc $1 }  = "$2";
			}
			else
			{
				$tmp_buffer .= " " . $line  . "# ";
			}
		}
		elsif( $line =~ m/^\s*\S+\s*\{/ )
		{
			die "Mangled input\n$line\n";
		}
		elsif( $line =~ m/\{.+\}/ )
		{
			$line =~ s/\{/ /;
			$line =~ s/\}/ /;
			$have_continue = 0;
			$tmp_buffer .= " " . $line;
			$tmp_buffer =~ s/\s+/ /g;
#			print TEST "$tmp_buffer\n";
			chomp($tmp_buffer);
			$tmp_buffer =~ m/^\s*(\S+)\s+(.+)/ or die "$tmp_buffer\nFailed to parse.\n";
			$input_hash{ lc $1 }  = "$2";
		}
		elsif( $line =~ m/\{/ )
		{
			$have_brace = 1;
			$line =~ s/\{/ /;
			$tmp_buffer .= " " . $line . "# ";
		}
		elsif( $line =~ m/\S+/ )
		{
			$have_continue = 0;
			$tmp_buffer .= " " . $line;
			$tmp_buffer =~ s/\s+/ /g;
#			print TEST "$tmp_buffer\n";
			chomp($tmp_buffer);
			$tmp_buffer =~ m/^\s*(\S+)\s+(.+)/ or die "$tmp_buffer\nFailed to parse.\n";
			$input_hash{ lc $1 }  = "$2";
		}
		else{ die "Unexpected blank line\n"; }
	}
	elsif( $line =~ m/\{.+\}/ )
	{
		$line =~ s/\{/ /;
		$line =~ s/\}/ /;
		$line =~ s/\s+/ /g;
#		print TEST "$line\n";
		$tmp_buffer = $line;
		chomp($tmp_buffer);
		$tmp_buffer =~ m/^\s*(\S+)\s+(.+)/ or die "$tmp_buffer\nFailed to parse.\n";
		$input_hash{ lc $1 }  = "$2";
	}
	elsif( $line =~ m/^\S+\s*\{\s*\S+/ )
	{
		$have_continue = 1;
		$have_brace = 1;
		$line =~ s/\{/ /;
		$tmp_buffer = $line . "# ";
	}
	elsif( $line =~ m/^\S+\s*\{/ )
	{
		$have_continue = 1;
		$have_brace = 1;
		$line =~ s/\{/ /;
		$tmp_buffer = $line;
	}
	elsif( $line =~ m/\S+\s+\S+/ )
	{
#		print TEST "$line\n";
		$tmp_buffer = $line;
		chomp($tmp_buffer); 
		$tmp_buffer =~ m/^\s*(\S+)\s+(.+)/ or die "$tmp_buffer\nFailed to parse.\n";
		$input_hash{ lc $1 }  = "$2";
	}
	elsif( $line =~ m/\S+/ )
	{
		$have_continue = 1;
		$tmp_buffer = $line;
	}
	else
	{
		die "Unexpected blank line\n";
	}
#	$input_string .= $line . " " ;
}

for( keys %input_hash )
{
	my $derp = $input_hash{$_};
	$derp =~ s/\#\s+/\n/g;
	chomp($derp);
	print TEST "$_ >>>>\n" .  $derp . "\n";
#	print TEST "$_ " .  $input_hash{$_} . "\n";
}

close TEST;
close IN;

print "Input read in\n";

my $filename;
my $value;

foreach my $key (@required_list)
{
	$filename = $filename_hash{ $key };
	open OUT, ">$filename" or die "$!\n";
	unless( $value = $input_hash{ $key } )
	{
		print "Required input $key  not found\n";
		close OUT;
		die;
	}
	$value =~  s/\#\s+/\n/g;
	chomp($value);
	print OUT "$value\n";
	close OUT;
}

my $key;
for( keys %optional_hash )
{
	$key = $_;
	$filename = $filename_hash{ $key };
	open OUT, ">$filename" or die "$!\n";
	unless( $value = $input_hash{ $key } )
	{
		$value  = $optional_hash{ $key };
	}
	$value =~  s/\#\s+/\n/g;
	chomp($value);
	print OUT "$value\n";
	close OUT;
}




exit;
