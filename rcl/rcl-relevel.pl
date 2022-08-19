#!/usr/bin/perl

use strict;
use warnings;

my $LEVEL = shift || die "Please supply maximum number of levels for new hierarchy\n";
die "Number of levels should be positive\n" unless $LEVEL > 0;

my $N = $LEVEL - 1;

my $header = <>;
chomp $header;

die "Header not recognised\n" unless $header eq "level\ttree\ttype\tjoinval\tN1\tN2\tnesting\tid\tnodes";
                                                #  0      1     2     3      4   5      6      7   8

my $curnest = "";
my @output  = ();

while (<>) {

  chomp;
  my @F = split "\t";
  push @F, "" if $F[5] == 0 && @F < 9;      # N1 is 0, split dropped last column.
  my $nesting = $F[6];
  my @nesting = split "_", $F[6];
  my $mynest = @nesting > $N ? join '_', @nesting[0..$N] : $nesting;
  $F[6] = $mynest;
  $F[8] = [ split " ", $F[8] ];

  if ($mynest eq $curnest) {
    push @{$output[-1][8]}, @{$F[8]};
    $output[-1][5] +=  $F[5];
    # print STDERR "Join $mynest $nesting\n";
  }
  else {
    push @output, [ @F ];
    $curnest = $mynest;
    # print STDERR "Push $mynest\n";
  }
}


my $newid = 1;

print "$header\n";

for my $rec (@output) {
  my @items = sort { $a <=> $b } @{$rec->[8]};
  local $" = "\t";
  print "@$rec[0..6]";
  local $" = ' ';
  print "\t$newid\t@items\n";
  $newid++;
}

my $n = @output;
print STDERR "Found $n clusters at level/depth $LEVEL\n";

