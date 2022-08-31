#!/usr/bin/perl

use strict;
use warnings;

my $LEVEL = shift || die "Please supply maximum number of levels for new hierarchy\n";
die "Number of levels should be positive\n" unless $LEVEL > 0;

my $N = $LEVEL - 1;

my $header = <>;
chomp $header;

die "Header not recognised\n" unless $header eq "level\ttree\ttype\tjoinval\tN1\tN2\tnesting\telements";
                                                #  0      1     2     3      4   5      6      7

my $curnest = "";
my @output  = ();

my $current_parent_size = 0;

while (<>) {

  chomp;
  my @F = split "\t";
  push @F, "" if $F[5] == 0 && @F < 8;      # N2 is 0, split dropped last column.
  my $nesting = $F[6];
  my @nesting = split "_", $nesting;
  my $mynest = @nesting > $N ? join '_', @nesting[0..$N] : $nesting;
  $F[6] = $mynest;
  $F[7] = [ split " ", $F[7] ];

  $current_parent_size = $F[4] if $mynest =~ /_A$/ || $mynest eq 'A';
  $F[4] = $current_parent_size;

  if ($mynest eq $curnest) {
    push @{$output[-1][7]}, @{$F[7]};
    $output[-1][5] +=  $F[5];
    # print STDERR "Join $mynest $nesting\n";
  }
  else {
    push @output, [ @F ];
    $curnest = $mynest;
    # print STDERR "Push $mynest\n";
  }
}


print "$header\n";

for my $rec (@output) {
  my @items = sort { $a <=> $b } @{$rec->[7]};
  local $" = "\t";
  print "@$rec[0..6]";
  local $" = ' ';
  print "\t@items\n";
}

my $n = @output;
print STDERR "Found $n clusters at level/depth $LEVEL\n";

