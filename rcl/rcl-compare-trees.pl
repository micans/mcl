#!/usr/bin/perl

use strict;
use warnings;

sub read_tree {
  my ($fname, $min, $max) = @_;
  my @cls = ();
  open(F, "<$fname") || die "Cannot open $fname\n";
  my $header = <F>;
  die unless $header eq "level\tsize\tjoinval\tnesting\tnodes\n";
  while(<F>) {
    chomp;
    my ($l, $s, $j, $n, $set) = split "\t";
    my @set = split /\s+/, $set;
    die "$l count mismatch\n" unless $s == @set;
    push @cls,
      { size => $s, label => $l, joinval => $j, nodes => { map { ( $_, 1 ) } @set } }
      if $s >= $min && (!$max || $s <= $max);
  }
  return \@cls;
}

sub best_jaccard {
  my ($set1, $multiviewlist) = @_;
  my @best = (0, 0, 0, 0);
  my $nodes1 = $set1->{nodes};

  for my $set2 (@$multiviewlist) {
    my $nodes2 = $set2->{nodes};
    my $d2 = grep { !defined($nodes1->{$_}) } keys %$nodes2;
    my $d1 = grep { !defined($nodes2->{$_}) } keys %$nodes1;
    my $sh = (scalar keys %$nodes1) - $d1;
    my $jc = $sh/($sh+$d1+$d2);
    my $c1 = $sh/($sh+$d1);
    my $c2 = $sh/($sh+$d2);
    if ($jc > $best[0]) {
      @best = ($jc, $c1, $c2, $set2);
    }
  }
  return \@best;
}

my $f1 = shift || die "Need two files\n";
my $f2 = shift || die "Need two files\n";
my $minclsize = shift || die "Need two files and minimum cluster size\n";

my $list1 = read_tree($f1, $minclsize, 0);
my $list2 = read_tree($f2, $minclsize, 0);

my $lc1 = @$list1;
my $lc2 = @$list2;

local $" = "\t";
my @header = qw(jc size1 size2 jv1 jv2 jc2 s1 s2);
print "@header\n";
for my $set1 (@$list1) {
  my $best = best_jaccard($set1, $list2);
  my $set2 = $best->[3];
  print "$best->[0]\t$set1->{size}\t$set2->{size}\t$set1->{joinval}\t$set2->{joinval}\t@$best[0..2]\n";
}


