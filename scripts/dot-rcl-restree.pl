#!/usr/bin/perl

# Reads the .hi.RESOLUTION.resdot file that describes the resolution tree,
# with node information for levels and link for edges.
# Creates a GraphViz definition for a plot.
# If RCL_DOT_HI_MISSING is set, extra labels are displayed, indicating
# how many cells of a parent are not accounted for by the children displayed,
# i.e. cells that were split off into smaller clusters.

use strict;
use warnings;

my @levelinfo = ();
my @links  = ();
my %is_parent = ();

while (<>) { chomp;
  my @F = split "\t", $_;
  my $type = shift @F;
  if ($type eq 'node') {    push @levelinfo, [$F[0], $F[1], $F[2], $F[3]]; }   # name value, size, missing
  elsif ($type eq 'link') { push @links, [$F[0], $F[1]]; $is_parent{$F[0]}++; }
}
@levelinfo = sort { $b->[1] <=> $a->[1] } @levelinfo;
push @levelinfo, ['dummy', 0, 0, 0];


  # RCL scores live in 0-1000. By /50 the range becomes 0-20.
  # Graphviz widths/lengths and scaling ways are a bit dark and dangerous,
  # but 0-20 seems to help. Maybe we try /100 someday.
  #
my @rulerlevels = map { "lev$_" } 0..(@levelinfo-2);
push @rulerlevels, map
  { my $i=$_;my $j=$i+1; my $d=sprintf("%.2f",(($levelinfo[$i][1]-$levelinfo[$j][1])/50));
    "lev$j -> lev$i [minlen=$d]"
  } 0..(@levelinfo-2);

my @treelevels = map { qq{ { rank = same; lev$_ $levelinfo[$_][0]; }} } 0..(@levelinfo-2);
pop @levelinfo;

# Outcommented lines has xlable for number of nodes 'peeling off', but only relative to displayed
# clusters. Clutters the display and needs too much explaining. Here in case
# it's useful as a debug/inspection thing at some point.
my @treenames  = defined($ENV{RCL_DOT_HI_MISSING}) ? map
  { my $name = $_->[0];
    defined($is_parent{$name})
      ? qq{$name [label=$_->[2], xlabel=< <font point-size="40">[$_->[3]]</font> >];}
      : qq{$name [label=$_->[2]];}
  } @levelinfo
: map { qq{$_->[0] [label=$_->[2]];} } @levelinfo;


{ local $" = "\n";
  print <<EOH;
  digraph g {
    node [shape="circle", width=1.0, fixedsize=true, label="" ];
      edge [arrowhead=none]
      ranksep = 0.2;
      subgraph levels {
        label="levels";
        node [style="invis", shape=point, width=0.01];
        edge [style="invis"];
  @rulerlevels
      }
      subgraph tree {
      node [width=2, fontsize=40];
EOH

  print "@treelevels\n@treenames\n";
}

for my $it (@links) {
  print "$it->[0] -> $it->[1]\n";
}
print "}}\n";


