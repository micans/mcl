#!/usr/bin/perl

# Reads the .hi.RESOLUTION.resdot file that describes the resolution tree,
# with node information for levels and link for edges.
# Creates a GraphViz definition for a plot.
# If RCL_DOT_HI_MISSING is set and true, extra labels are displayed, indicating
# how many cells of a parent are not accounted for by the children displayed,
# i.e. cells that were split off into smaller clusters.

# example:
# RCL_DOT_HI_MISSING=0 dot-rcl-reshi.pl < def.hi.200-500-1250-3000.resdot > r8
# dot -Tpdf -Gsize=10,10\! < r8 > r8.pdf


use strict;
use warnings;

my %nodeinfo = ();
my @links  = ();
my $maxrung = 0;

$ENV{RCL_DOT_HI_MISSING} = 0 unless defined($ENV{RCL_DOT_HI_MISSING});

while (<>) { chomp;
  my @F = split "\t", $_;
  my $type = shift @F;
  if ($type eq 'node') {
    my $val = $F[1];                    # float in range 0-1000.0
    my $rung = int(($val-0.01)/20);     # integer in range 0-199.
    $nodeinfo{$F[0]} =
    { name => $F[0]
    , val  => $val
    , ival => int($val)
    , rung => $rung
    , size => $F[2]
    , miss => $F[3]
    , is_parent => 0
    };
    $maxrung = $rung if $rung > $maxrung;
  } elsif ($type eq 'link') { push @links,
    { parent => $F[0]
    , child  => $F[1]
    };
    $nodeinfo{$F[0]}{is_parent}++;        # dangersign; depends on links following nodes.
  }
}

print STDERR "highest rung at $maxrung\n";
my @rulerlevels =  map { "lev$_" } 0..($maxrung+1);
push @rulerlevels, map { my $j=$_+1; "lev$_ -> lev$j [minlen=1.0]" } 0..$maxrung;
my @treelevels = ();

# check if any parent+child are on the same rung. Moan if so.
# then ... push child down/up if possible (todo).
# down the tree/drawing, but leaves in tree are high in rung.
#
while (1) {
  my $amends = 0;
  for my $link (@links) {
    my ($p, $c) = ($link->{parent}, $link->{child});
    if ($nodeinfo{$p}{rung} == $nodeinfo{$c}{rung}) {
      print STDERR "Woe betides $p $c $nodeinfo{$p}{rung}\n";
      $nodeinfo{$c}{rung} += 1;
      $amends++;
    }
  }
  last unless $amends;
}

# now pin the nodes to the $maxrung-rung ladder.
for my $rung (0..$maxrung) {
  my @besties = grep { $nodeinfo{$_}{rung} == $rung } keys %nodeinfo;
  push @treelevels, qq{{ rank = same; lev$rung @besties; }};
}

my @treenames  = $ENV{RCL_DOT_HI_MISSING} ? map
  { $nodeinfo{$_}{is_parent}
      ? qq{$_ [label=$nodeinfo{$_}{size}, xlabel=< <font point-size="40">[$nodeinfo{$_}{miss}]</font> >];}
      : qq{$_ [label=$nodeinfo{$_}{size}];}
  } keys %nodeinfo
: map { qq{$_ [label=$nodeinfo{$_}{size}];} } keys %nodeinfo;


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
  print "$it->{parent} -> $it->{child}\n";
}
print "}}\n";


