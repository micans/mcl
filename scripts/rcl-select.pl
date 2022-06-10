#!/usr/bin/perl

# Copyright 2022 Stijn van Dongen

# This program is free software; you can redistribute it and/or modify it
# under the terms of version 3 of the GNU General Public License as published
# by the Free Software Foundation. It should be shipped with MCL in the top
# level directory as the file COPYING.

# rcl-select.pl :
# (1) From a binary tree pick sets of internal nodes that represent balanced flat clusterings
# (2) From a binary tree pick internal nodes that represent significant merge and pre-merge states
# Only reads STDIN, which should be the output of clm close in --sl mode.  That
# output encodes the single-linkage join order of a tree.  The script further
# requires a prefix for file output and a list of resolution sizes.

# --- The first mode of output ---
# The output is a list of flat clusterings, one for each resolution size.
# These clusterings usually share clusters between them (i.e. clusters do not
# always split at each resolution level), but do form a (not strictly) nesting
# set of clusterings.  Also output is the 'dot' specification of a plot that
# shows the structure of the hierarchy (ignoring clusters below the smallest
# resolution size).  This file can be put through GraphViz dot to obtain the
# plot.
#
# A cluster corresponds to a tree node. The cluster consists of all associated
# leaf nodes below this node. For a given resolution size R each cluster C must
# either be of size at least R without a sub-split below C's tree node into two
# other clusters of size at least R, or C is smaller than R and was split off
# in order to allow another such split to happen elsewhere. In the last case
# C will not have been split any further.
#
# For decreasing resolution sizes, the code descends each node in the tree, as
# long as it finds two independent components below the node that are both of
# size >= resolution.  For each resolution size the internal nodes that encode
# the clustering for that resolution are marked.  After this stage, the
# clusterings for the different resolutions are output, going back up the tree
# from small resolution / fine-grained clusters to larger resolution /
# coarse-grained clusters, and merging or copying clusters from the previous
# stage.

# --- The second mode of output ---
# The tree is descended, and any split where the smallest subtree is at least
# size reslimit (the smallest specified resolution) is taken.
# Additionally, any clusters found in the first mode are added if not found
# by this method.


# rcl.sh incorporates rcl-select.pl, see there for comprehensive usage example.
# Use e.g.
#     rcl-select.pl pfx 50 100 200 < sl.join-order
#     mcxload -235-ai pfx50.info -o pfx50.cls

# TODO:
# equijoin not called if one branch has lower ival. Document reasoning or reconsider.
# (detect circular/nonDAG input to prevent memory/forever issues (defensive))

use strict;
use warnings;
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);


# Globals yes, too lazy for now to package into a state object.

$::jiggery = defined($ENV{RCL_JIGGERY})? $ENV{RCL_JIGGERY} : 1;
die "Only jiggery 1 is currently available\n" unless $::jiggery == 1;

$::prefix = shift || die "Need prefix for file names";
die "Need at least one resolution parameter\n" unless @ARGV;
for my $r (@ARGV) {
  die "Resolution check: strange number $r\n" unless looks_like_number($r);
}
@::resolution = sort { $a <=> $b } @ARGV;
$::reslimit = $::resolution[0];
$::reslimithi = $::resolution[-1];
$::resolutiontag = join '-', @::resolution;

@ARGV = ();
%::nodes = ();
$::nodes{dummy}{items} = [];     # used for singletons; see below
$::nodes{dummy}{size}  = 0;      #
$::nodes{dummy}{lss}   = 0;      #

$::L=1;
%::topoftree = ();


sub pick_hierarchy {

  my $pick = {};
  my ($toplevelstack, $lowlimit,$hilimit, $pick_by_level) = @_;

  for my $node (@$toplevelstack)
  {  my $root   = {  name => 'ROOT'
                  , pickdepth => 0
                  , treedepth => 0
                  , longname  => ''
                  , bigsplit  => 0
                 };
    my $state  = {  name => $node
                  , visit => 0
                  , fraycount => 0  # number of fray nodes between this and resolution node.
                  , treedepth => 1
                  , pickdepth => 0
                  , longname  => ''
                  , debug     => 0
                  , lastpick  => ""
                  , bigsplit  => 0};
                                
    my @stack = ($root, $state);

    my $huh = 0;

    while (@stack > 1) {

      my $state   = $stack[-1];
      my $pstate  = $stack[-2];       # parent state

      my $name    = $state->{name};
      my $pname   = $pstate->{name};
      my $size    = $::nodes{$name}{size};
      my $ival    = $::nodes{$name}{ival};

      my ($ann, $bob) = map { $::nodes{$name}{$_} } qw(ann bob);

      my $mrk = '';

      if ($state->{visit} == 0) {

        $state->{bigsplit} = 1 * ($::nodes{$ann}{size} >= $lowlimit && $::nodes{$bob}{size} >= $lowlimit);

        my $bothlss = 1 * ($::nodes{$ann}{lss} >= $lowlimit && $::nodes{$bob}{lss} >= $lowlimit);
        my $somelss = 1 * ($::nodes{$ann}{lss} >= $lowlimit || $::nodes{$bob}{lss} >= $lowlimit);
        my $onelss  = 1 * ($somelss && !$bothlss);
        my $nonelss = 1 * (!$somelss);
        my $pbigsplit = $pstate->{bigsplit};

        $mrk .= 'o' if ($pstate->{name} ne 'ROOT' && $pbigsplit && $::nodes{$pstate->{name}}{size} > $hilimit && $size <= $hilimit);
        $mrk .= 's' if $state->{bigsplit};
        $mrk .= '=' if $bothlss;
        $mrk .= '-' if $onelss;
        $mrk .= '|' if $nonelss;
        if ($pbigsplit && $nonelss) {;      # these are close merges generally. 
          $mrk .= $ival >= 1.05 * $::nodes{$pname}{ival} ? 'y' : 'x';
        }
        $mrk .= 'p' if $pbigsplit;

        my $wanted =  $state->{bigsplit};
        $wanted    =  1 if $pbigsplit && $nonelss;     # this is a close merge generanecdotally
        if (!$wanted && defined($pick_by_level->{$name})) {
          $mrk .= 'f';
          $wanted = 1;
        }

        # print STDERR "(potential equijoin) $name $::nodes{$name}{lss} $::nodes{$ann}{size} $::nodes{$bob}{size}\n" if $state->{bigsplit} && $::nodes{$name}{lss} < $lowlimit;

        if ($wanted) {
          my $ival = $::nodes{$name}{ival};
          my $longname = $pstate->{longname} ? $pstate->{longname} . '::' . $name : $name;

          if ($size <= $hilimit) {
            $state->{pickdepth}++;
            $state->{longname}  = $longname;      # $state2->{longname} = $longname;      # fixme, this is a bit unhappy. hv. what if successive nodes are chosen - same longname?
            $state->{lastpick}  = $name;
            if ($pstate->{lastpick}) {
              push @{$pick->{$pstate->{lastpick}}{children}}, $name;
            }

            $pick->{$name}{longname} = $state->{longname};
            $pick->{$name}{mark} = $mrk;
            $pick->{$name}{pname}= $pstate->{lastpick};
            $pick->{$name}{level}= $state->{pickdepth};
            $pick->{$name}{children}  = [];
          }
        }
      }

      my $state2  = { %$state };      # need to inherit a few things ..
      $state2->{visit} = 0;           # building this for ann, then for bob anew.
      $state2->{treedepth}++;

      $state->{visit} += 1;

      if ($state->{visit} == 3) {
        pop @stack;
      }
      else {
        if ($state->{bigsplit}) {
          $state2->{fraycount} = 0;
        }
        if ($state->{visit} == 1) {
          if ($state->{bigsplit} || 2 * $::nodes{$ann}{lss} >= $lowlimit) {
            $state2->{fraycount} += $::nodes{$bob}{size} unless $state->{bigsplit};
            $state2->{name} = $ann;
            push @stack, $state2;
          }
        }
        elsif ($state->{visit} == 2) {
          if ($state->{bigsplit} || 2 * $::nodes{$bob}{lss} >= $lowlimit) {
            $state2->{fraycount} += $::nodes{$ann}{size} unless $state->{bigsplit};
            $state2->{name} = $bob;
            push @stack, $state2;
          }
        }
      }
    }
  }
  return $pick;
}


sub flat_pick_levels {

  my $ips = shift;
  my @inputstack = @$ips;       # typically the top of trees, representing connected components in network.

  my @clustering = ();
  my %resolutionstack = ();
  my %pick_by_level = ();
  
  for my $res (sort { $b <=> $a } @::resolution) { print STDERR " .. $res";

    while (@inputstack) {

      my $name = pop @inputstack;
      my $ann  = $::nodes{$name}{ann};
      my $bob  = $::nodes{$name}{bob};
   
      if ($::jiggery == 1) {

         my $la = $ann eq 'null' ? '-' : $::nodes{$ann}{lss};
         my $lb = $bob eq 'null' ? '-' : $::nodes{$bob}{lss};
         my $action = 'retain';

         if ($::nodes{$name}{size} == 1 || (2 * $::nodes{$ann}{lss} <= $res && 2 * $::nodes{$bob}{lss} <= $res)) {
           push @clustering, $name;
           $pick_by_level{$name} = $res;
         }
         else {                 # there is a merge node of size >= $res at/below either ann or bob.
                                # so we descend down to such merge nodes (discarding e.g. volatile nodes)
                                # OTOH once such a merge node does not exist we stop descending (so
                                # keeping volatile nodes). TBC.
           $action = 'descend';
           push @inputstack, $ann;
           push @inputstack, $bob;
         }
         # print STDERR "\n$res $action $name (lss $la $lb $ann $bob)\n" if $name eq 'L33169_445';
         # print STDERR "\n$res $action $name (lss $la $lb $ann $bob)\n" if $ann eq 'L33169_445' || $bob eq 'L33169_445';
      }
      else {
         die "No other jiggery is available right now\n";
      }
    }

     # make copy, as we re-use clustering as inputstack.
     #
    $resolutionstack{$res} = [ @clustering ];
    @inputstack = @clustering;
    @clustering = ();
  }

  return (\%resolutionstack, \%pick_by_level);
  print STDERR "---\n";
}

sub dot_full_tree {
  if (defined($ENV{RCL_RES_DOT_TREE}) && $ENV{RCL_RES_DOT_TREE} == $::L) {
    open(DOTTREE, ">$::prefix.joindot") || die "Cannot open $::prefix.joindot";
    for my $node
    ( grep { $::nodes{$_}{size} > 1 }
      sort { $::nodes{$b}{val} <=> $::nodes{$a}{val} }
      keys %::nodes
    ) {
      my $val = $::nodes{$node}{val};
      my $ann = $::nodes{$node}{ann};
      my $bob = $::nodes{$node}{bob};
      print DOTTREE "$node\t$::nodes{$node}{ann}\t$val\t$::nodes{$ann}{size}\n";
      print DOTTREE "$node\t$::nodes{$node}{bob}\t$val\t$::nodes{$bob}{size}\n";
    }
    close(DOTTREE);
  }
}


sub dot_flat_tree {

  my $flatpick = shift;
  my $dotname = "$::prefix.hi.$::resolutiontag.resdot";
  open(RESDOT, ">$dotname") || die "Cannot open $dotname for writing";

  for my $n (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } keys %$flatpick) {
    my $size = $::nodes{$n}{size};
    my $sum  = 0;
    my $pctloss = "0";
    $sum += $::nodes{$_}{size} for @{$flatpick->{$n}{children}};
    $pctloss = sprintf("%d", 100 * ($size - $sum) / $size) if $sum;
    print RESDOT "node\t$n\t$::nodes{$n}{val}\t$size\t$pctloss\n";
  }
  for my $n1 (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } keys %$flatpick) {
    for my $n2 (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } @{$flatpick->{$n1}{children}} ) {
      print RESDOT "link\t$n1\t$n2\n"; # Could implement filter here
                                       # to avoid printing out the very smallest stuff.
    }
  }
  close(RESDOT);
}


sub read_full_tree {

  my $epsilon = 0.2;
  my $header = <>;
  chomp $header;
  my $header_expect = "link\tval\tNID\tANN\tBOB\txcsz\tycsz\txycsz\tiss\tlss\tannid\tbobid";

  die "Join order header line not recognised (expect [$header_expect])" unless $header eq $header_expect;
  print STDERR "-- constructing tree:\n";

  my $equijoin_report = defined($ENV{RCL_EQUIJOIN}) ? $ENV{RCL_EQUIJOIN} : 0;

  while (<>) {

     chomp;
     my @F = split "\t";

     die "Expect 12 elements (have \"@F\")\n" unless @F == 12;
     my ($i, $val, $upname, $ann, $bob, $xcsz, $ycsz, $xycsz, $iss, $lss, $annid, $bobid) = @F;
     die "Checks failed on line $.\n" unless
           looks_like_number($xcsz) && looks_like_number($ycsz)
        && looks_like_number($iss) && looks_like_number($lss);
     print STDERR '.' if $. % 1000 == 1;

                        # leaves have to be introduced into our tree/node listing
     if ($xcsz == 1) {
        $ann =~ /leaf_(\d+)/ || die "Missing leaf (Ann) on line $.\n";
        my $leafid = $1;
        $::nodes{$ann} =
        {  name => $ann
        ,  size =>  1
        ,  items => [ $leafid ]
        ,  ann => "null"
        ,  bob => "null"
        ,  csizes => []
        ,  lss => 0
        ,  iss => 0
        ,  val => 1000
        ,  ival => 1000
        } ;
     }
     if ($ycsz == 1) {
        $bob =~ /leaf_(\d+)/ || die "Missing leaf (Bob) on line $.\n";
        my $leafid = $1;
        $::nodes{$bob} =
        {  name => $bob
        ,  size =>  1
        ,  items => [ $leafid ]
        ,  ann => "null"
        ,  bob => "null"
        ,  csizes => []
        ,  lss => 0
        ,  iss => 0
        ,  val => 1000
        ,  ival => 1000
        } ;
     }

     # LSS: largest sub split. keep track of the maximum size of the smaller of
     # any pair of nodes below the current node that are not related by
     # descendancy.  Given a node N the max min size of two non-nesting
     # nodes below it is max(mms(desc1), mms(desc2), min(|desc1|, |desc2|)).
     # clm close and rcl-select.pl both compute it - a bit pointless but lets just
     # call it a sanity check.

     # ISS: immediate sub split.

     # $ann eq $bob is how clm close denotes a singleton in the network - the
     # only type of node that does not participate in a join.
     # A dummy node exists (see above) that has only size, items and lss with none
     # of the other fields set. Currently that node is only accessed when
     # items are picked up in the cluster aggregation step. If code is added
     # and pokes at other attributes they will be undefined and we will know.

     $bob = 'dummy' if $ann eq $bob;
     die "Parent node $upname already exists\n" if defined($::nodes{$upname});

     my $equijoin = 0;
     $equijoin = $val + $epsilon >= $::nodes{$ann}{val} && $val + $epsilon >= $::nodes{$bob}{val} unless $upname =~ /^sgl_/;

     my $properjoin = $equijoin ? 0 : 1;
     if ($equijoin && $equijoin_report && $::nodes{$ann}{size} >= $equijoin_report && $::nodes{$bob}{size} >= $equijoin_report) {
        print STDERR "-- equijoin $upname($val) $ann($::nodes{$ann}{val}) $bob($::nodes{$bob}{val}) $::nodes{$ann}{size} $::nodes{$bob}{size}\n";
     }

     $iss =  $properjoin * min($::nodes{$ann}{size}, $::nodes{$bob}{size});

     $::nodes{$upname} =
     {  name  => $upname
     ,  parent => undef
     ,  size  => $::nodes{$ann}{size} + $::nodes{$bob}{size}
     ,  ann   => $ann
     ,  bob   => $bob
     ,  csizes => [ $::nodes{$ann}{size}, $::nodes{$bob}{size}]
     ,  iss   => $iss
     ,  lss   => max( $::nodes{$ann}{lss}, $::nodes{$bob}{lss}, $iss )
     ,  val   => $val
     ,  ival  => int($val + 0.5)
     ,  parent => ''
     } ;

     $::nodes{$ann}{parent} = $upname;
     $::nodes{$bob}{parent} = $upname;

  #  print STDERR "LSS error check failed ($ann $bob)\n" if $::nodes{$upname}{lss} != $lss && $ann ne $bob;
  #  above check no longer ok because of equijoin.
  #  these comments as a reminder of the fact that clm close computes lss and equijoin is important -
  #  this program's lss overrides clm-close-lss.

     delete($::topoftree{$ann});
     delete($::topoftree{$bob});

     $::topoftree{$upname} = 1;
     $::L++;
  }
print STDERR "\n" if $. >= 1000;
  return [ grep { $::nodes{$_}{size} > 1 } sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } keys %::topoftree ];
}


sub set_flat_level {
  my ($fp, $name, $level) = @_;
  $fp->{$name}{level} = $level;
  for (@{$fp->{$name}{children}}) {
    set_flat_level($fp, $_, $level+1);
  }
}

sub flat_cls_collect {

  my $resolutionstack = shift;

  print STDERR "\n-- collecting clusters for resolution\n";
   # when collecting items, proceed from fine-grained to coarser clusterings,
   # so with low resolution first.
   #
  my %flatpick  = ();
  my @datasizes = ();

  my $res_index = 0;

  for my $res (sort { $a <=> $b } @::resolution) { print STDERR " .. $res";

    $res_index++;

    my $clsstack = $resolutionstack->{$res};
    my $datasize = 0;

    local $" = ' ';
    my $fname = "$::prefix.res$res.info";
    open(OUT, ">$fname") || die "Cannot write to $fname";

    for my $name ( sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } @$clsstack ) {

      my $size = $::nodes{$name}{size};
      my $ival = int(0.5 + $::nodes{$name}{val});
      my @nodestack = ( $name );
      my @items = ();

      if (!defined($flatpick{$name})) {
        $flatpick{$name} = {};
        $flatpick{$name}{children} = [];
      }

      while (@nodestack) {

        my $nodename = pop(@nodestack);
            # Below items are either cached from a previous more fine-grained clustering
            # or they are leaf nodes
        if (defined($::nodes{$nodename}{items})) {

          push @items, @{$::nodes{$nodename}{items}};

          if ($nodename ne $name) {
          # if (defined($flatpick{$name})) 
            push @{$flatpick{$name}{children}}, $nodename;
               # the above depends on the fact that anytime we find items,
               # we are guaranteed that nodename is an immediate subclustering
               # of name (there is nothing inbetween);
               # this depends on how we find clusters iteratively
               # with the resolution parameter descreasing in steps.
               # It's not entirely neat, this dependency.
          }
        }
        else {
          push @nodestack, ($::nodes{$nodename}{ann}, $::nodes{$nodename}{bob});
        }
      }
      @items = sort { $a <=> $b } @items;
      $::nodes{$name}{items}    = \@items unless defined($::nodes{$name}{items});

      $flatpick{$name}{children} = [] unless defined($flatpick{$name}{children});
   
      my $nitems = @items;
      print STDERR "Error res $res size difference $size / $nitems\n" unless $nitems == $size;
   
      print OUT "$name\t$ival\t$size\t@items\n";
      $datasize += $size;
      $flatpick{$name}{tag} .= $res_index;
    }
    push @datasizes, $datasize;
    close(OUT);
  }
  local $" = ' ';
  print STDERR " cls node counts: @datasizes\n";

  for my $name (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } @{$resolutionstack->{$::reslimithi}} ) {
    set_flat_level(\%flatpick, $name, 1);
  }

  return (\%flatpick, \@datasizes);
}


# note the flatpick hierarchy has clusters with size < $::reslimit,
# so below needs checking and picking.

sub printlistnode {

  my ($pick, $fh, $level, $nodelist, $ni, $parent) = @_;

  my $size = $::nodes{$ni}{size};
  my $ival = $::nodes{$ni}{ival};
  # return unless $size >= $::reslimit;       # perhaps argumentise.

  my $sumofchildren = 0;
  my @children = grep { $::nodes{$_}{size} >= $::reslimit } @{$pick->{$ni}{children}};
  my $up    = $parent ? $ival - $::nodes{$parent}{ival} : $ival;
  my $down  =   @children
              ? (sort { $a <=> $b } map { $::nodes{$_}{ival} - $ival } @children)[0]
              : '-';
  for (@children) {
    $sumofchildren += $::nodes{$_}{size};
  }

  my $nloss = $sumofchildren ? sprintf("%d", 100 * ($size - $sumofchildren) / $size) : "-";
  my $tag = join('::', (@{$nodelist}, $ni));
  local $" = ' ';
die "suprisingly no items for [$ni]\n" if !defined($::nodes{$ni}{items});
  print $fh "$level\t$size\t$ival\t$nloss\t$up\t$down\t$tag\t@{$::nodes{$ni}{items}}\n";
  for my $nj (sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} } @children ) {
    printlistnode($pick, $fh, $level+1, [ @$nodelist, $ni ], $nj, $ni);
  }
}

# fixme datasize: embed check/definition properly.

sub print_hierarchy {
  my ($type, $listname, $pick, $datasize) = @_;
     # This output encodes the top-level hierarchy of the RCL clustering,
     # with explicit levels, descendancy encoded in concatenated labels,
     # and all the nodes contained within each cluster.

  my $toplevelsum = 0;

  my @huh = grep { !defined($pick->{$_}{level}) } keys %$pick;
  die "No level for @huh\n" if @huh;

  my @toplevelnames = sort { $::nodes{$b}{size} <=> $::nodes{$a}{size} }
                      grep { $pick->{$_}{level} == 1 && $::nodes{$_}{size} >= $::reslimit } keys %$pick;
  $toplevelsum += $::nodes{$_}{size} for @toplevelnames;

  my $nloss = $toplevelsum ? sprintf("%d", 100 * ($datasize - $toplevelsum) / $datasize) : "-";
  print "$type: $toplevelsum $datasize ($nloss loss at first level)\n";

  my $down = (sort { $a <=> $b } map { $::nodes{$_}{ival} } @toplevelnames)[0];

  open(RESLIST, ">$listname") || die "Cannot open $listname for writing";
  print RESLIST "level\tsize\tjoinval\tpctloss\tup\tdown\tnesting\tnodes\n";
  print RESLIST "0\t$datasize\t0\t$nloss\t0\t$down\troot\t-\n";

  for my $n (@toplevelnames)
  {   printlistnode($pick, \*RESLIST, 1, [], $n, '');
  }
  close(RESLIST);
}

sub get_sibling {

  my $name = shift;
  my $sib = "";

  if ($::nodes{$name}{parent}) {
    my $ann = $::nodes{$::nodes{$name}{parent}}{ann};
    my $bob = $::nodes{$::nodes{$name}{parent}}{bob};
    if ($ann eq $name) {
      $sib = $bob;
    } elsif ($bob eq $name) {
      $sib = $ann;
    }
    else { die "sibbobannhuh $name\n"; }
  }
  return $sib;
}


sub get_node_items {

  my $name = shift;
  my @nodestack = ( $name );
  return $::nodes{$name}{items} if defined($::nodes{$name}{items});

  my @items = ();

  while (@nodestack) {

    my $nodename = pop(@nodestack);
        # Below items are either cached from a previous more fine-grained clustering
        # or they are leaf nodes
    if (defined($::nodes{$nodename}{items})) {
      push @items, @{$::nodes{$nodename}{items}};
    }
    else {
      push @nodestack, ($::nodes{$nodename}{ann}, $::nodes{$nodename}{bob});
    }
  }
  @items = sort { $a <=> $b } @items;
  $::nodes{$name}{items} = \@items;
  return \@items;
}


sub hich_cls_collect {
  my $hichpick = shift;
  #
  # get longer names / deeper nesting clusters before shorter names / top level clusters.
  # Anything found is cached, this helps efficiency.
  #
  for (sort { $hichpick->{$b}{longname} cmp $hichpick->{$a}{longname} } keys %$hichpick) {
    get_node_items($_);
  }
}


sub hich_dump {

  my ($hichpick, $flatpick) = @_;

  my $treename1 = "$::prefix.subtree.$::resolution[0]-$::resolution[-1].txt";
  my $treename2 = "$::prefix.subtreecls.$::resolution[0]-$::resolution[-1].txt";
  open(SUBTREE1, ">$treename1") || die "Cannot open $treename1 for writing";
  open(SUBTREE2, ">$treename2") || die "Cannot open $treename2 for writing";

  for (sort { $hichpick->{$b}{longname} cmp $hichpick->{$a}{longname} } keys %$hichpick) {
    # print "$_\t$hichpick->{$_}{longname}\n";
    die "Expected items for $_\n" unless defined($::nodes{$_}{items});
    my $items = $::nodes{$_}{items};
    local $"    = "\t";
    my $level   = $hichpick->{$_}{level};
    my $parent  = $level == 1 ? '-' : $hichpick->{$_}{pname};
    print SUBTREE2 "$level\t$_\t$parent\t$::nodes{$_}{ival}\t$::nodes{$_}{size}\t$hichpick->{$_}{longname}\t@$items\n";
  }
  close(SUBTREE2);

  my $prefix = "";

  for my $name (sort {$hichpick->{$a}{longname} cmp $hichpick->{$b}{longname}} keys %$hichpick) {
    my $longname = $hichpick->{$name}{longname};
    if ($longname !~ /^$prefix/) {
      print SUBTREE1 "\n";
    }
    $prefix=$longname;
    my $ival = $::nodes{$name}{ival};
    my $size = $::nodes{$name}{size};
                      # need two checks to avoid creating $flatpick->{$name}.
    my $fp  = defined($flatpick->{$name}) && defined($flatpick->{$name}{tag}) ? $flatpick->{$name}{tag} : '-';
    my $mrk = $hichpick->{$name}{mark};
    my $lss = $::nodes{$name}{lss};
    my $sib = get_sibling($name);
    my $sibinfo = $sib ? "s=$::nodes{$sib}{size}/$::nodes{$sib}{lss}" : "";
    my @children = ();
    my @children_info = ();
    @children = @{$hichpick->{$name}{children}};
    @children_info = map { "$_\[v=$::nodes{$_}{ival}][i=$::nodes{$_}{iss}][l=$::nodes{$_}{lss}]" } @children;
    my $gap = 2000;
    for my $c (@children) {
      my $g = $::nodes{$c}{ival} - $ival;
      $gap = $g if $g < $gap;
    }
    $gap = '-' unless @children;
    my $nchild = @children;
    my $ann = $::nodes{$name}{ann};
    my $bob = $::nodes{$name}{bob};
    my $chdinfo = "c=$::nodes{$ann}{lss}/$::nodes{$bob}{lss}";

    local $" = ', ';
    print SUBTREE1 "i=$ival\tl=$lss\tg=$gap\t$fp\t$mrk\t$longname    ($hichpick->{$name}{pname} <$lss> @children_info)\n";

  } close(SUBTREE1);
}

sub dump_subtree {

  my ($root, $res) = @_;
  my @stack = ([$root, "", 1]);
  print "Searching $root with $res\n";

  while (@stack) {
    my $item = $stack[-1];
    my ($name, $longname, $pbigsplit) = @$item;
    my $ann = $::nodes{$name}{ann};
    my $bob = $::nodes{$name}{bob};
    $::nodes{$name}{visit} = 0 unless defined($::nodes{$name}{visit});
    my $bigsplit = 1 * ($::nodes{$ann}{size} >= $res && $::nodes{$bob}{size} >= $res);
    if ($::nodes{$name}{visit} == 0) {
      if ($pbigsplit || $bigsplit) {
        $item->[1] .= "::$name";          # fixme hv check other stack code for longname.
        my $lss = $::nodes{$name}{lss};
        my $ival = $::nodes{$name}{ival};
        print "$ival\t$lss\t$item->[1]\n";
      }
    }
    $::nodes{$name}{visit}++;
    if ($::nodes{$name}{visit} == 1 && ($::nodes{$ann}{lss} >= $res || $bigsplit)) {
      push @stack, [$ann, $item->[1], $bigsplit];
    }
    elsif ($::nodes{$name}{visit} == 2 && ($::nodes{$bob}{lss} >= $res || $bigsplit)) {
      push @stack, [$bob, $item->[1], $bigsplit];
    }
    elsif ($::nodes{$name}{visit} == 3) {
      pop @stack;
    }
  }
}


if ($::prefix =~ s/^\+//) {
  my $toplevelstack = read_full_tree();
  dump_subtree($::prefix, $::reslimit);
  exit 0;
}


my $toplevelstack = read_full_tree();
dot_full_tree();


      print STDERR "-- computing resolution hierarchy (toplevelstack @$toplevelstack)\n";
my ($resolutionstack, $pick_by_level) = flat_pick_levels($toplevelstack);
my ($flatpick, $datasizes) = flat_cls_collect($resolutionstack);
      print STDERR "---\n";

      print STDERR "-- computing tree descend hierarchy\n";
my $hichpick = pick_hierarchy($toplevelstack, $::resolution[0], $::resolution[-1]*2, $pick_by_level);
hich_cls_collect($hichpick);
      print STDERR "---\n";

dot_flat_tree($flatpick);
hich_dump($hichpick, $flatpick);


my $flatname = "$::prefix.hi.$::resolutiontag.txt";
print_hierarchy('flat', $flatname, $flatpick, $datasizes->[0]);


my $hichname = "$::prefix.hic.$::resolutiontag.txt";
print_hierarchy('hier', $hichname, $hichpick, $datasizes->[0]);

