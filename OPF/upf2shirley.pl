#!/usr/bin/perl
# Copyright (C) 2021 OCEAN collaboration
#
# This file is part of the OCEAN project and distributed under the terms 
# of the University of Illinois/NCSA Open Source License. See the file 
# `License' in the root directory of the present distribution.
#
#
use strict;
use File::Basename;

die "Expected Pseudo file as input" unless( scalar @ARGV > 0 );
my $fileName = $ARGV[0];
my $outFileName = $fileName . ".mod";
if( scalar @ARGV > 1 ) {
  $outFileName = $ARGV[1];
}

open IN, $fileName or die "Failed to open $fileName";

my @rad;
my @local;
my @l0;
my @wfc0;
my @l1;
my @wfc1;
my @l2;
my @wfc2;
my @l3;
my @wfc3;

while( my $line = <IN> )
{
#  if( $line =~ m/</ )
#  {
#    print "### $line";
#  }
  my $val;
  my $tag;
  if( $line =~ m/</ )
  {
    next if( $line =~ m/<\// );
    $tag = $line;
    $tag =~ s/\<//;
    until( $tag =~ m/>/ )
    {
      $line = <IN>;
      $tag .= $line;
    }
    $tag =~ s/>//;
  }
  if( 0 ){
  if( $line =~ m/<([^\/].*)>/ )
  {
    $tag = $1;
    print "$tag\n";
  }
  elsif( $line =~ m/>/ )
  {  print '#' . $line; }
  elsif( $line =~ m/</ )
  {
    $tag = chomp($line);
    print "$tag\n";
    while( 1 ) {
      $line = <IN>;
      $tag .= $line;
      last if( $line =~ m/>/ );
    }
    print "$tag\n";
  }
  }
  if( $tag =~ m/PP_R$|PP_R\s|PP_LOCAL|PP_BETA|PP_PSWFC/ )
  {
    
    print "$tag\n";
    my $line2;
    do {
      $line2 = <IN>;
      if( $line2 =~ m/\s*</ ) {
        until( $line2 =~ m/>/) {
          chomp $line2;
          $line2 .= ' ' . <IN>;
        }
      }
      $val .= $line2;
#      print "$tag  $line2" if( $tag =~ m/PP_PSWFC/ ) ;
    } until( $line2 =~ m/<(.*)>/ );
#    print "------- $line2";

    ### PP_CHI
    if( $line2 =~ m/PP_CHI/ ) {
#      print "!!!!!!!!!!!!!!!!\n";
      my $line3;
      do {
        $tag = "PP_CHI";
        if( $line2 =~ m/\s*</ ) {
          until( $line2 =~ m/>/) {
            chomp $line2;
            $line2 .= ' ' . <IN>;
          }
        }
        $line2 =~ m/l="(\d+)"/ or die "parse failed\n$line2";
        my $el = $1;
#        print $line2;
        $val = '';
        do {
          $line3 = <IN>;
          $val .= $line3;
        } until( $line3 =~  m/<(.*)>/ );
#        print "$el -------- $line3";

        my @temp = split ' ', $val;
        my $i = 0;
        if( $el == 0 && scalar @wfc0 == 0 ) { 
#          print "   l = 0\n";
          while( $temp[$i] =~ m/\d\.\d+/ )
          {
            push @wfc0, $temp[$i];
            $i++;
          }
        } elsif( $el == 1 && scalar @wfc1 == 0 ) {
#          print "   l = 1\n" . scalar @temp . "\n$temp[0] $temp[10] 10\n";
#          print $val;
          while( $temp[$i] =~ m/\d\.\d+/ )
          {
            push @wfc1, $temp[$i];
            $i++;
          }
        } elsif( $el == 2 && scalar @wfc2 == 0 ) {
#          print "   l = 2\n";
          while( $temp[$i] =~ m/\d\.\d+/ )
          {
            push @wfc2, $temp[$i];
            $i++;
          }
        } elsif( $el == 3 && scalar @wfc3 == 0 ) {
#          print "   l = 3\n";
          while( $temp[$i] =~ m/\d\.\d+/ )
          {
            push @wfc3, $temp[$i];
            $i++;
          }
        } 
        $line2 = <IN>;

      } while ( $line2 =~ m/PP_CHI/);

    }
  

    if ( $tag =~ m/PP_R$|PP_R\s/ )
    {
      @rad = split ' ', $val;
      pop @rad;
    } 
    elsif ( $tag =~ m/PP_LOCAL/ )
    {
      @local = split ' ', $val;
      pop @local;
    }
    elsif( $tag =~ m/PP_BETA/ )
    {
      my @temp = split ' ', $val;
      my $l;
      if( $tag =~ m/angular_momentum="(\d)"/ )
      {
        $l = $1;
        print "### angular $l\n";
      } 
      else
      {
        my $beta = shift @temp;
        $l = shift @temp;
        shift @temp; shift @temp; shift @temp;
        print "### angular2 $l\n";
      }
      if( $l == 0 )
      {
        die "Multiple projectors not supported for conversion!" if( scalar @l0 > 0 );
        print "$temp[0] $temp[1] $temp[-1]\n";
        @l0 = @temp;
        print "$l0[0]\n";
        pop @l0;
      }
      elsif( $l == 1 )
      {
        die "Multiple projectors not supported for conversion!" if( scalar @l1 > 0 );
        @l1 = @temp;
        print "$l1[0]\n";
        pop @l1;
      }
      elsif( $l == 2 )
      {
        die "Multiple projectors not supported for conversion!" if( scalar @l2 > 0 );
        @l2 = @temp;
        print "$l2[0]\n";
        pop @l2;
      }
      elsif( $l == 3 )
      {
        die "Multiple projectors not supported for conversion!" if( scalar @l3 > 0 );
        @l3 = @temp;
        print "$l3[0]\n";
        pop @l3;
      }
      else{ die "l = $l" }
    }
    elsif( $tag =~ m/PP_PSWFC/ )
    {
      print "####### PP_PSWFC\n";
      my @temp = split ' ', $val;
      print scalar @temp . "\n";
      my $i = 0;
      do
      {
        if( $temp[$i] =~ m/^\d*s/i && scalar @wfc0 == 0 )
        {
          print "Gabbed s pswfn\n";
          $i++; $i++; $i++; $i++;
          while( $temp[$i] =~ /\d\.\d+/ )
          {
            push @wfc0, $temp[$i];
            $i++;
          }
        }
        elsif( $temp[$i] =~ m/^\d*p/i && scalar @wfc1 == 0 )
        {
          print "Gabbed p pswfn\n";
          $i++; $i++; $i++; $i++;
          while( $temp[$i] =~ /\d\.\d+/ )
          {
            push @wfc1, $temp[$i];
            $i++;
          }
        }
        elsif($temp[$i] =~ m/^\d*d/i && scalar @wfc2 == 0 )
        {
          print "Gabbed d pswfn\n";
          $i++; $i++; $i++; $i++;
          while( $temp[$i] =~ /\d\.\d+/ )
          {
            push @wfc2, $temp[$i];
            $i++;
          }
        }
        elsif($temp[$i] =~ m/^\d*f/i && scalar @wfc3 == 0 )
        {
          $i++; $i++; $i++; $i++;
          while( $temp[$i] =~ /\d\.\d+/ )
          {
            push @wfc3, $temp[$i];
            $i++;
          }
        }
        else
        {
          $i++;
        }
      }while( $i < scalar @temp );
    }
  }
} 
close IN;

my $lmax = 1;
$lmax++ if( scalar @l0 > 0 );
$lmax++ if( scalar @l1 > 0 );
$lmax++ if( scalar @l2 > 0 );
$lmax++ if( scalar @l3 > 0 );
$lmax = 4 if ( $lmax > 4 );


print scalar @wfc1 . "\n";

open OUT, ">", $outFileName;
print OUT "  "  . $lmax . "   " . (scalar @rad ) ."\n";
print OUT "0\n";
my $sum = 0;
for( my $i = 0; $i < scalar @l0; $i++ )
{
  if( $wfc0[$i] > 0.0000000000001 ) {
    print OUT "$rad[$i]  " . ($l0[$i]/(2*$wfc0[$i])+$local[$i]/2) . "\n";
  } else {
    print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
  }
}
for( my $i = scalar @l0; $i < scalar @local; $i++ )
{
  print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
}
#print OUT "\n\n";
print OUT "1\n";
for( my $i = 0; $i < scalar @l1; $i++ )
{
  if( $wfc1[$i] > 0.0000000000001 ) {
    print OUT "$rad[$i]  " . ($l1[$i]/(2*$wfc1[$i])+$local[$i]/2) . "\n";
  } else {
    print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
  }
}
for( my $i = scalar @l1; $i < scalar @local; $i++ )
{
  print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
}
#print OUT "\n\n";
print OUT "2\n";
for( my $i = 0; $i < scalar @l2; $i++ )
{
  if( $wfc2[$i] > 0.0000000000001 ) {
    print OUT "$rad[$i]  " . ($l2[$i]/(2*$wfc2[$i])+$local[$i]/2) . "\n";
  } else {
    print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
  }
}
for( my $i = scalar @l2; $i < scalar @local; $i++ )
{
  print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
}
#print OUT "\n\n";
print OUT "3\n";
for( my $i = 0; $i < scalar @l3; $i++ )
{
  if( $wfc3[$i] > 0.0000000000001 ) {
    print OUT "$rad[$i]  " . ($l3[$i]/(2*$wfc3[$i])+$local[$i]/2) . "\n";
  } else {
    print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
  }
}
for( my $i = scalar @l3; $i < scalar @local; $i++ )
{
  print OUT "$rad[$i]  " . $local[$i]/2 . "\n";
}
close OUT;



