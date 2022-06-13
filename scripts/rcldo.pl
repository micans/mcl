#!/usr/bin/perl

use strict;
use warnings;

my $mode = shift;


cumgra() if $mode eq 'cumgra';
distwrangle() if $mode eq 'distwrangle';

if ($mode eq 'clstag') {
  my $tag = get_tag(@ARGV);
  my $scale = $ENV{RCLPLOT_PARAM_SCALE} || 0;
  $tag /= 10 ** $scale;
  print "$tag\n";
}


sub get_tag {

  my $fname = shift;
  die "No file name supplied" unless $fname;
  $fname =~ /[rI](\d{2,})\b/ || die "Cannot find resolution parameter from file [$fname]\n";
  return $1;
}


sub distwrangle {

    my $T = shift @ARGV;
    my $scale = $ENV{RCLPLOT_PARAM_SCALE} || 0;

    while (<>) {
      chomp;
      /d1=(\d+)\s+d2=(\d+).*?[rI](\d{2,})\b.*?[rI](\d{2,})/ || die "No match on line [$_]\n";
      my $d1=$1;      # Leading zeroes do not lead to octal
      my $d2=$2;      # ^
      my $x=$3;
      my $y=$4;
      if ($scale) {
        $x /= 10**$scale;
        $y /= 10**$scale;
      }
      print $d1/$T, "\t$x\t$y\n";
      print $d2/$T, "\t$y\t$x\n";
   }
}


sub cumgra {

  my @x = 0..80;

  my @y = grep { $_ <= 10000 } map { 10 ** ($_/10) } @x;
  # print "@y\n";

  my $runningtotal = 0;
  my $y = 0;
  my $size = 0;

  my $fname = shift @ARGV;
  my $totalsize = shift @ARGV;

  my $tag = get_tag($fname);

  my $scale = $ENV{RCLPLOT_PARAM_SCALE} || 0;
  $tag /= 10**$ENV{RCLPLOT_PARAM_SCALE};

  while (<>) {
    chomp;
    $size = $_;
    while ($size > $y[$y]) {
      my $fraction = $runningtotal/$totalsize;
      print "$tag\t$y\t$fraction\t$y[$y]\n";
      $y++;
    }
    $runningtotal += $size;
  }
  while ($y < @y) {
    my $fraction = $runningtotal/$totalsize;
    print "$tag\t$y\t$fraction\t$y[$y]\n";
    $y++;
  }

}


#while ($size > $y[$y]) {
#  print "$y\t$runningtotal\n";
#  $y++;
#}

