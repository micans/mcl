#!/usr/bin/perl -an

# Reads the output of clm close in --sl mode.
#
# Then proceeds to descend each node in the tree, as long as it finds two
# independent components below the node that are both of size >= script-argument
#
# Use e.g. rcl-mix.pl 200 sl.join-order  | sort -nr | cut -f 2 | mcxload -235-ai - -o res200.cls
# This will reverse sort on cluster size (larger clusters first).


use strict;
use warnings;
use List::Util qw(min max);


BEGIN {
  $::min = shift || die "Need min max";
  %::nodes = ();
  %::cid2node = ();
  $::L=1;
}
next if $. == 1;
my ($i, $x, $y, $xid, $yid, $val, $xcid, $ycid, $xcsz, $ycsz, $xycsz, $nedge, $ctr, $lss, $nsg) = @F;


my ($n1, $n2) = ($xcid, $ycid);
($n1, $n2)    = ($ycid, $xcid) if $ycid < $xcid;

my $upname    = "L$::L.$n1";

if ($xcsz == 1) {
  $::cid2node{$xcid} = "L0.$xcid";
  $::nodes{"L0.$xcid"} =
  {    name => "L0.$xcid"
  ,    size =>  1
  ,   items => [ $xid ]
  ,     ann => ""
  ,     bob => ""
  ,  csizes => []
  , lss => 0
  , nsg => 0
  } ;
}
if ($ycsz == 1) {
  $::cid2node{$ycid} = "L0.$ycid";
  $::nodes{"L0.$ycid"} =
  {    name => "L0.$ycid"
  ,    size =>  1
  ,   items => [ $yid ]
  ,     ann => ""
  ,     bob => ""
  ,  csizes => []
  , lss => 0
  , nsg => 0
  } ;
}

my $name1 = $::cid2node{$n1};
my $name2 = $::cid2node{$n2};


# Keep track of the maximum size of the smaller of any pair of nodes below the current node that are
# not related by descendancy.

# Given a node N; what is the max min size of two non-nesting nodes below it.
# Its max(mms(desc1), mms(desc2), min(|desc1|, |desc2|))

   $::nodes{$upname} =
   {   name  => $upname
   ,   parent => undef
   ,   size  => $::nodes{$name1}{size} + $::nodes{$name2}{size}
   ,   items => [ @{$::nodes{$name1}{items}}, @{$::nodes{$name2}{items}} ]
   ,    ann  => $name1
   ,    bob  => $name2
   ,  csizes => [ $::nodes{$name1}{size}, $::nodes{$name2}{size}]
   , lss => max( $::nodes{$name1}{lss}, $::nodes{$name2}{lss}, min($::nodes{$name1}{size}, $::nodes{$name2}{size}))
   , nsg => $nsg
   } ;

$::cid2node{$n1} = $upname;
$::cid2node{$n2} = $upname;

$::L++;
$::upname = $upname;

local $" = ' ';


END {
  my @stack = $::upname;
  while (@stack) {
    my $name = pop @stack;

    my $size = $::nodes{$name}{size};
    my @desc = @{$::nodes{$name}{csizes}};
    my $ann  = $::nodes{$name}{ann};
    my $bob  = $::nodes{$name}{bob};
    my $nsg  = sprintf("%.3f", $::nodes{$name}{nsg} / $::nodes{$name}{size});

    local $" = ' ';
    if ($::nodes{$name}{lss} >= $::min) {
      push @stack, $ann;
      push @stack, $bob;
    }
    else {
      print "$size\t$nsg\t@{$::nodes{$name}{items}}\n";
    }
  }
}


