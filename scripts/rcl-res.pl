#!/usr/bin/perl

# Only reads STDIN, which should be the output of clm close in --sl mode.  That
# output encodes the single-linkage join order of a tree.  The script further
# requires a prefix for file output and a list of resolution sizes.
#
# A cluster for a given resolution size R corresponds to a node in the tree
# where the cluster consists of all associated leaf nodes below. This cluster
# is either at least of size R and has no sub-split in the tree into two
# components each of size at least R, or it is smaller than R and was split off
# in order to allow another such split to happen elsewhere.
#
# For decreasing resolution sizes, the code descends each node in the tree, as
# long as it finds two independent components below the node that are both of
# size >= resolution.  For each resolution size the internal nodes that encode
# the clustering for that resolution are marked.  After this stage, the
# clusterings for the different resolutions are output.

# rcl.sh incorporates rcl-res.pl, see there for comprehensive usage example.
# Use e.g.
#     rcl-res.pl pfx 50 100 200 < sl.join-order
#     mcxload -235-ai pfx50.clusters -o pfx50.cls


use strict;
use warnings;
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);

$::prefix = shift || die "Need prefix for file names";
die "Need at least one resolution parameter\n" unless @ARGV;
for my $r (@ARGV) {
  die "Resolution check: strange number $r\n" unless looks_like_number($r);
}
@::resolution = sort { $a <=> $b } @ARGV;
$::reslimit = $::resolution[0];
$::resdisplaylimit = defined($ENV{RCL_RES_PLOT_LIMIT}) ? $ENV{RCL_RES_PLOT_LIMIT} : $::resolution[0];

@ARGV = ();
%::nodes = ();
%::cid2node = ();
$::L=1;
%::topoftree = ();
print STDERR "-- constructing tree:\n";

my $header = <>;
chomp $header;

die "Join order header line not recognised" unless $header =~ /^link.*nsg$/;

while (<>) {

   chomp;
   my @F = split "\t";

   die "Expect 15 elements (have \"@F\")\n" unless @F == 15;
   my ($i, $x, $y, $xid, $yid, $val, $xcid, $ycid, $xcsz, $ycsz, $xycsz, $nedge, $ctr, $lss, $nsg) = @F;
   die "Checks failed on line $.\n" unless
         looks_like_number($xid) && looks_like_number($yid)
      && looks_like_number($lss) && looks_like_number($nsg);
   print STDERR '.' if $. % 1000 == 1;

   my ($id1, $id2) = ($xcid, $ycid);
   ($id1, $id2)    = ($ycid, $xcid) if $ycid < $xcid;

   my $upname    = "L$::L" . "_$id1" . "_$xycsz";

                      # singletons have to be introduced into our tree/node listing
   if ($xcsz == 1) {
     $::cid2node{$xcid} = "L0_$xcid";
     $::nodes{"L0_$xcid"} =
     {    name => "L0_$xcid"
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
     $::cid2node{$ycid} = "L0_$ycid";
     $::nodes{"L0_$ycid"} =
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

   my $node1 = $::cid2node{$id1};     # id1 id2 are xcid ycid; these were assigned and updated
   my $node2 = $::cid2node{$id2};     # along join order for descendant nodes by clm close.


   # Keep track of the maximum size of the smaller of any pair of nodes below the current node that are
   # not related by descendancy.

   # Given a node N; what is the max min size of two non-nesting nodes below it.
   # Its max(mms(desc1), mms(desc2), min(|desc1|, |desc2|))

      $::nodes{$upname} =
      {   name  => $upname
      ,   parent => undef
      ,   size  => $::nodes{$node1}{size} + $::nodes{$node2}{size}
      ,    ann  => $node1
      ,    bob  => $node2
      ,  csizes => [ $::nodes{$node1}{size}, $::nodes{$node2}{size}]
      , lss => max( $::nodes{$node1}{lss}, $::nodes{$node2}{lss}, min($::nodes{$node1}{size}, $::nodes{$node2}{size}))
      , nsg => $nsg
      } ;

   # clm close outputs a line with id1 == id2 for all singleton nodes.
   print STDERR "LSS error check failed ($id1 $id2)\n" if $::nodes{$upname}{lss} != $lss && $id1 ne $id2;

   $::cid2node{$id1} = $upname;
   $::cid2node{$id2} = $upname;

   delete($::topoftree{$node1});
   delete($::topoftree{$node2});

   $::topoftree{$upname} = 1;
   $::L++;

}

print STDERR "\n" if $. >= 1000;
my @inputstack = ( sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } keys %::topoftree );
my @clustering = ();
my %resolutionstack = ();

print STDERR "-- computing tree nodes for resolution";
 # Start from top of tree(s), so we find the larger-size-resolution nodes
 # first.  inputstack is a set of nodes for which we know that they are
 # higher (or equal) in the tree relative to the nodes that answer our
 # resolution request.  At a resolution step, we can use the answer obtained
 # for the previous step as the new inputstack.
 #
for my $res (sort { $b <=> $a } @::resolution) { print STDERR " .. $res";

  while (@inputstack) {

    my $name = pop @inputstack;
    my $ann  = $::nodes{$name}{ann};
    my $bob  = $::nodes{$name}{bob};
 
    if ($::nodes{$name}{lss} >= $res) {
      push @inputstack, $ann;
      push @inputstack, $bob;
    }
    else {
      push @clustering, $name;
    }
  }

   # make copy, as we re-use clustering as inputstack.
   #
  $resolutionstack{$res} = [ @clustering ];
  @inputstack = @clustering;
  @clustering = ();
}

print STDERR "\n-- collecting clusters for resolution";
 # when collecting items, proceed from fine-grained to coarser clusterings,
 # so with low resolution first.
 #
my %digraph = ();   # collect links for dot plot
my %digraph_printname = ();

for my $res (sort { $a <=> $b } @::resolution) { print STDERR " .. $res";

  my $clsstack = $resolutionstack{$res};

  local $" = ' ';
  my $fname = "$::prefix.res$res.info";
  open(OUT, ">$fname") || die "Cannot write to $fname";

  for my $name ( sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } @$clsstack ) {

    my $size = $::nodes{$name}{size};
    my @nodestack = ( $name );
    my @items = ();

    while (@nodestack) {

      my $nodename = pop(@nodestack);
          # Below items are either cached from a previous more fine-grained clustering
          # or they are leaf nodes
      if (defined($::nodes{$nodename}{items})) {

        push @items, @{$::nodes{$nodename}{items}};

        if ($nodename ne $name && $::nodes{$nodename}{size} >= $::reslimit) {
          $digraph{$name}{$nodename} = 1;
        }
      }
      else {
        push @nodestack, ($::nodes{$nodename}{ann}, $::nodes{$nodename}{bob});
      }
    }
    @items = sort { $a <=> $b } @items;
    $::nodes{$name}{items} = \@items unless defined($::nodes{$name}{items});
    $digraph_printname{$name} = "$size" if $size >= $::resdisplaylimit;
 
    my $nitems = @items;
    print STDERR "Error res $res size difference $size / $nitems\n" unless $nitems == $size;
 
    my $nsg  = sprintf("%.3f", $::nodes{$name}{nsg} / $::nodes{$name}{size});
    print OUT "$name\t$size\t$nsg\t@items\n";
  }
  close(OUT);
}

open(DIGRAPH, ">$::prefix.digraph") || die "Cannot open $::prefix.digraph for writing";
print DIGRAPH <<EOT;
digraph g {
  forcelabels = true;
EOT

for my $n (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } keys %digraph_printname ) {
  my $size = $::nodes{$n}{size};
  my $sum  = 0;
  my $missing = "0";
  if (defined($digraph{$n})) {
    $sum += $::nodes{$_}{size} for keys %{$digraph{$n}};
    $missing = sprintf("%d", 100 * ($size - $sum) / $size);
  }
  print DIGRAPH qq{  $n [label="$digraph_printname{$n}", xlabel="[$missing]"];\n};
}
for my $n1 (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } keys %digraph_printname ) {
  for my $n2 (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } keys %{$digraph{$n1}} ) {
    print DIGRAPH "  $n1 -> $n2;\n" if $::nodes{$n2}{size} >= $::resdisplaylimit;
  }
}
print DIGRAPH "}\n";
close(DIGRAPH);

print STDERR "\n";


