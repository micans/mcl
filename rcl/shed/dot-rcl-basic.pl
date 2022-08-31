#!/usr/bin/perl

use strict;
use warnings;

print <<EOH;
digraph g {
  node [shape="circle", width=0.20, fixedsize=true, label="" ];
    edge [arrowhead="none"];
    L142_30_47 [shape=invtriangle, color=blue]
    L124_x55 [style=filled, fillcolor=red]
EOH

while (<>) {
  chomp;
  my ($p, $c, $v, $Nc) = split "\t";
  print "$p -> $c;\n";
}
print "}\n";

