use strict;

use JSON::PP;

#use JSON::SL;
#use JSON::Streaming::Reader;

#open IN, "oparse.description.json" or die;
#my $text;
#while( <IN> ) {
#  $text .= $_;
#}
#my $sl = JSON::SL->new();
#$sl->object_drip(1);
#my @res = $sl->feed($text);
#
#for( my $i=0; $i<10; $i++ ) {
#  print $res[$i]->{'Value'} . "  " . $res[$i]->{'Path'} . "\n";
#}
#
#if ( 0 ) {
#my $fh;
#my $jsonr = JSON::Streaming::Reader->for_stream($fh);
#
#my $level = 0;
#while( my $token = $jsonr->get_token) {
#    my ($key, $value) = @$token;
#    if ($key eq 'start_property') {
#        $level++;
#    } elsif ( $key eq 'end_property' ) {
#        $level--;
#    } elsif ( $key eq 'start_object' ) {
#      $token = $jsonr->get_token;
#      ($key, $value) = @$token;
#      print "$level  $key  $value\n";
#    }
#    
#}   
#
#
#}




my $json = JSON::PP->new;
my $description;
if( open my $in, "oparse.description.json" ) {
  local $/ = undef;
  $description = $json->decode(<$in>);
  close($in);
} else {
  die "Failed to open oparse.description.json\n$!";
}

my $defaults;
if( open my $in, "oparse.json" ) {
  local $/ = undef;
  $defaults = $json->decode(<$in>);
  close($in);
} else {
  die;
  $defaults = {};
}

my $type;
if( open my $in, "oparse.type.json" ) {
  local $/ = undef;
  $type = $json->decode(<$in>);
  close($in);
} else {
  $type = {};
}


my $level = 0;
my $tree = '';
print "# OCEAN inputs\n";
recursiveJSON( $description, $defaults, $type, $level, $tree );


# add in oparse.json and oparse.type.json files
# sub to pull apart "bse.etc." into perl hash look up, and if exists print out default and type info

sub recursiveJSON  {
  my ($hashRef, $defRef, $typRef, $level, $tree ) = @_;

  my $localLevel = $level + 1;
  my $localTree;# = $tree;
#  my @tmp = keys %{$hashRef};
#  my @keys = sort { $a <=> $b }@tmp; 
  foreach my $key (sort keys %{$hashRef})  {


    unless ( $key eq '.' ) {
      for( my $i =0; $i< $localLevel; $i++ ) {
        print "#"
      }
      if( $tree eq '' ) {
        $localTree = $key
      } else {
        $localTree = "$tree.$key";
      }
      print " <a name=\"$localTree\"></a>$key\n" ;
    }
    # Special thing for the references, which atm are at the top
    if( $localTree eq 'references' ) {
      my $refCount = 0;
      foreach my $ref (sort keys %{$hashRef->{"$key"}})  {
        $refCount++;
        print "$refCount. [$ref]: " . $hashRef->{"$key"}->{"$ref"} . "\n";
      }
      foreach my $ref (sort keys %{$hashRef->{"$key"}})  {
        $refCount++;
        print "[$ref]: " . $hashRef->{"$key"}->{"$ref"} . "\n";
      }
      print "\n\n";
      next;
    }
    if( ref $hashRef->{"$key"} eq 'HASH' ) {
      print "\n";
#      my $localTree;
#      $localTree = $tree . ".$key";
      recursiveJSON( $hashRef->{"$key"}, $defRef, $typRef, $localLevel, $localTree );
    } else {
      print $hashRef->{"$key"} . "\n\n";
      unless( $key eq '.' ) {
#        print "$localTree  ";
        my $t = extract( $localTree, $defRef );
        print "Default : $t\n" if( $t ne '' );
        my $t = extract( $localTree, $typRef );
        print "Type : $t\n" if( $t ne '' );
        print "\n";
      }
#      print "$t\n" unless( $t eq '' );
    }
  }
}
  


sub extract {
  my ($key, $hashRef ) = @_;

  my $outString = '';
  my @key = split /\./, $key;
#  print ">>>>>> $key $key[0]  $hashRef->{'bse'}->{'core'}->{'broaden'} \n";
  for( my $i = 0; $i < scalar @key - 1; $i++ ) {
    return '' unless( exists $hashRef->{$key[$i]} );
    $hashRef = $hashRef->{$key[$i]};
#    print "$key[$i] ";
  }
#  print "Default value: ";
  if( ref($hashRef->{$key[-1]}) eq 'ARRAY' ) {
#    print "[ ";
    $outString .= "[ ";
    foreach (@{$hashRef->{$key[-1]}}) {
#      print $_ . " ";
      $outString .= "$_ ";
    }
    $outString .= "]";
#    print "]\n";
  } else {
#    print "$hashRef->{$key[-1]}\n";
    if(  $hashRef->{$key[-1]} =~ /a*[ifbs]/ ) {
      $outString = $hashRef->{$key[-1]};
    } elsif( $hashRef->{$key[-1]} == JSON::PP::false ) {
      $outString = "False";
    } elsif ( $hashRef->{$key[-1]} == JSON::PP::true ) {
      $outString = "True";
    } else {
      $outString = $hashRef->{$key[-1]};
    }
  }

  return $outString;
#  return $hashRef;
}
