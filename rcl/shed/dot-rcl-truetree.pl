#!/usr/bin/perl

# reads the .joindot file that can optionally be made rcl-res.pl
#

use strict;
use warnings;

my %joins = ();

while (<>) {
  chomp;
  my ($p, $c, $v, $Nc) = split "\t";
  push @{$joins{$v}}, [ $p, $c ];
}

my @levels = map { $_ + 0.01 } sort { $b <=> $a } keys %joins;
push @levels, 0;
my @printlevels = map { "lev$_" } 0..(@levels-1);
push @printlevels, map { my $i=$_;my $j=$i+1; "lev$j -> lev$i" } 0..(@levels-2);

                          # L142_30_47 [shape=invtriangle, color=blue]
                          # L124_x55 [style=filled, fillcolor=red]
local $" = ";\n";
print <<EOH;
digraph g {
  node [shape="circle", width=0.20, fixedsize=true, label="" ];
    edge [arrowhead="none"];
    ranksep = 0.1;
    subgraph levels {
label="levels";
@printlevels
    }
    subgraph tree {
EOH


for my $v (sort {$b <=> $a} keys %joins) {

  my @edges = @{$joins{$v}};
  my %c2p = ();
  my %isparent = ();

  my $level = 0;
  $level++ while $levels[$level] > $v;

  for my $e (@edges) {
    $c2p{$e->[1]} = $e->[0];
    $isparent{$e->[0]} = 1;
  }

  my @sinks = ();
  my @leafnodes = grep { !defined($isparent{$_}) } keys %c2p;

  for my $c (@leafnodes) {
    my $p = $c2p{$c};
    $p = $c2p{$p} while defined($c2p{$p});
    print "$p -> $c;\n";
    push @sinks, $p;
  }
  local $" = ' ';
  print "{ rank = same; lev$level @sinks; }\n";
}

print "}}\n";


__DATA__

make connected components
  find the top of each (grandparent); this will be the representative.
    link all leaves in the component to the representative, remove the intermediates.
    everything keeps its name.

