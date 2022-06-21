#!/usr/bin/perl

# Reads the .hi.RESOLUTION.resdot file that describes the resolution tree,
# with node information for levels and link for edges.
# Creates a GraphViz definition for a plot.
# With --xlabel-remainder, extra labels are displayed, indicating
# how many cells of a parent are not accounted for by the children displayed,
# i.e. cells that were split off into smaller clusters.

use Getopt::Long;
use warnings;
use strict;

my @ARGV_COPY = @ARGV;
my $n_args    = @ARGV;

$::debug  = 0;
$::test   = 0;
my $help  = 0;
my $verbose = 0;
my $progname = 'rcl-dot-resmap.pl';

my $plot_xlabeltype     = "none";
my $plot_minres         = 0;
my $plot_labeltype      = 'size';

sub help {
   print <<EOH;
Usage:
   $progname [options]
Options:
--help              this
                    from a cluster as it splits into its subclusters.
--minres=<num>      do not display nodes of size below <num>
--label=<string>    size|ival|leaf|none
--xlabel=<string>   size|ival|peel currently only shown for parents
EOH
}

if
(! GetOptions
   (  "help"            =>   \$help
   ,  "test"            =>   \$::test
   ,  "debug=i"         =>   \$::debug
   ,  "minres=i"        =>   \$plot_minres
   ,  "label=s"         =>   \$plot_labeltype
   ,  "xlabel=s"        =>   \$plot_xlabeltype
   ,  "verbose"         =>   \$verbose
   )
)
   {  print STDERR "option processing failed\n";
      exit(1);
   }

if ($help) {
   help();
   exit 0;
}

die "Need size|ival|leaf\n" unless $plot_labeltype =~ /^(size|ival|leaf|none)$/;


my %nodeinfo = ();
my @links  = ();
my $maxrung = 0;
my $n_skip_edges = 0;
my $n_skip_nodes = 0;


while (<>) { chomp;
  my @F = split "\t", $_;
  my $type = shift @F;
  if ($type eq 'node') {
    die "Expect 4 data fields for node at line $.\n" unless @F == 4;
    my ($name, $val, $size, $peel) = @F;
    if ($size < $plot_minres) {
      $n_skip_nodes++;
      next;
    }
    my $rung = int(($val-0.01)/20);     # integer in range 0-50.
    $nodeinfo{$name} =
    { name => $name
    , val  => $val
    , ival => int($val)
    , rung => $rung
    , size => $size
    , peel => $peel
    , is_parent => 0
    , leaf => ""
    };
    if ($name =~ /_x(\S+)/) {
      $nodeinfo{$name}{leaf} = $1;
    }
                          # dangersign dependency on rcl-res naming convention
                          #
    $maxrung = $rung if $rung > $maxrung;
  } elsif ($type eq 'link') {
    die "Expect 2 data fields for node at line $.\n" unless @F == 2;
    my ($parent, $child) = @F;
    if (!defined($nodeinfo{$parent}) || !defined($nodeinfo{$child})) {
      $n_skip_edges++;
      next;
    }
    push @links,
    { parent => $F[0]
    , child  => $F[1]
    };
    $nodeinfo{$parent}{is_parent}++;    # dangersign; depends on links following nodes.
  }
}

print STDERR "Highest rung at $maxrung\n";
print STDERR "Blanked $n_skip_nodes nodes and $n_skip_edges edges\n";

my @rulerlevels =  map { "lev$_" } 0..($maxrung+1);
push @rulerlevels, map { my $j=$_+1; "lev$_ -> lev$j [minlen=1.0]" } 0..$maxrung;
my @treelevels = ();

# check if any parent+child are on the same rung. Moan if so.
# then ... push child down/up if possible (todo).
# down the tree/drawing, but leaves in tree are high in rung.
#
my $amends = 0;
while (1) {
  $amends = 0;
  for my $link (@links) {
    my ($p, $c) = ($link->{parent}, $link->{child});
    if ($nodeinfo{$p}{rung} >= $nodeinfo{$c}{rung}) {
      print STDERR "Equal rank $nodeinfo{$p}{rung} for $p and $c (amending)\n" if $verbose;
      $nodeinfo{$c}{rung} += 1;
      $amends++;
    }
  }
  last unless $amends;
}
print STDERR "Adjusted position for $amends nodes\n" if $amends;

# now pin the nodes to the $maxrung-rung ladder.
for my $rung (0..$maxrung) {
  my @besties = grep { $nodeinfo{$_}{rung} == $rung } keys %nodeinfo;
  push @treelevels, qq{{ rank = same; lev$rung @besties; }};
}

my @nodenames = ();
if ($plot_labeltype eq 'none') {
}
elsif ($plot_labeltype eq 'leaf') {
  @nodenames  = map { qq{$_ [label=$nodeinfo{$_}{leaf}];} } grep { length($nodeinfo{$_}{leaf}); } keys %nodeinfo;
}
elsif ($plot_labeltype =~ /^(size|ival)$/) {
  @nodenames  = map { qq{$_ [label=$nodeinfo{$_}{$plot_labeltype}];} } keys %nodeinfo;
}
if ($plot_xlabeltype =~ /^(size|ival|peel)$/) {
  push @nodenames,
    map {qq{$_ [xlabel=< <font point-size="40">[$nodeinfo{$_}{$plot_xlabeltype}]</font> >];} }
    grep { $nodeinfo{$_}{is_parent} } keys %nodeinfo;
}


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

  print "@treelevels\n@nodenames\n";
}

for my $it (@links) {
  print "$it->{parent} -> $it->{child}\n";
}
print "}}\n";


