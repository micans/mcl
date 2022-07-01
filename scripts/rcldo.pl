#!/usr/bin/perl

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my $mode = shift || die "Need mode";

  # Out of pure laziness and spite these are global for now so I don't need to pass them into
  # newick(). TODO revisit when code stabilises.
my %hm_tree  = ();
my %hm_nodes = ( root => { level => 0, type => 'cls', ival => 0, print => 0 } );
my $hm_order  =  1;


if ($mode eq 'cumgra') {
  cumgra();
}
elsif ($mode eq 'distwrangle') {
  distwrangle();
}
elsif ($mode eq 'label') {
  die "Need <tabfile>\n" unless @ARGV == 1;
  my $tabfile = shift;
  labelwrangle($tabfile);
}
elsif ($mode eq 'clstag') {
  die "Need file name\n" unless @ARGV == 1;
  my $tag = get_tag(@ARGV);
  my $scale = $ENV{RCLPLOT_PARAM_SCALE} || 0;
  $tag /= 10 ** $scale;
  print "$tag\n";
}
elsif ($mode eq 'heatannot' || $mode eq 'heatannotcls') {
  die "Need <annotationfile> <partitionhierarchyfile> <tabfile> <outputbase>\n" unless @ARGV == 4;
  my ($fnannot, $fnhier, $fntab, $fnbase) = @ARGV;
  my ($dfannot, $termlist) = read_node_annotation_table($fnannot);
  my $tab = read_tab($fntab);
  my $dispatch = $mode eq 'heatannot' ? 'rclhm' : 'cls';
  read_partition_hierarchy($dispatch, $fnhier, $dfannot, $termlist, $tab, $fnbase);
}
else {
  die "Unknown mode $mode\n";
}

$::median_warning = 0;

sub get_tag {

  my $fname = shift;
  die "No file name supplied" unless $fname;
  $fname =~ /[rI](\d{2,})\b/ || die "Cannot find resolution parameter from file [$fname]\n";
  return $1;
}


sub labelwrangle {
  my $tabfile = shift;
  my $tab = read_tab($tabfile);
  
  my $header = <>;
  chomp $header;
  my @header = split "\t", $header;
  my $N = @header;
  my %mapcount = ();

  my $node_index = 0;
  for (@header) {
    last if $_ eq 'nodes';
    $node_index++;
  }
  die "No nodes column found\n" unless $node_index < $N;
  while (<>) {
    chomp;
            # -1 below is to retain empty fields.
    my @F = split "\t", $_, -1;
    my $NF = @F;
    die "Column count mismatch ($N/$NF) on line $.\n" unless $N == $NF;
    my $nodelist = $F[$node_index];
    my @labels = ();
    for my $n (split " ", $nodelist, -1) {
      die "Weird entry [$n] on line $.\n" unless $n =~ /^[0-9]+$/;
      die "[$n] not present in tab on line $.\n" unless defined($tab->{$n});
      push @labels, $tab->{$n};
      $mapcount{$n}++;
    }
    $F[$node_index] = join " ", @labels;
    local $" = "\t";
    print "@F\n";
  }
  my %hist = ( 0 => 0 );
  for (values %mapcount) {
    $hist{$_}++;
  }
  for (sort { $a <=> $b } keys %hist) {
    printf STDERR "%6d keys were mapped %d times\n", $hist{$_}, $_;
  }
}

sub distwrangle {

    my $T = shift @ARGV;
    my $scale = $ENV{RCLPLOT_PARAM_SCALE} || 0;

    while (<>) {
      chomp;
      /d1=(\d+)\s+d2=(\d+).*?[rI](\d{2,})\b.*?[rI](\d{2,})/ || die "No match on line [$_]\n";
      my $d1=$1;      # Leading zeroes do not lead to octal
      my $d2=$2;      # ^
      my $x=$3;
      my $y=$4;
      if ($scale) {
        $x /= 10**$scale;
        $y /= 10**$scale;
      }
      print $d1/$T, "\t$x\t$y\n";
      print $d2/$T, "\t$y\t$x\n";
   }
}


sub cumgra {

  my @x = 0..80;

        # TODO make below flexible.
  my @y = grep { $_ <= 10000 } map { 10 ** ($_/10) } @x;

  my $runningtotal = 0;
  my $y = 0;
  my $size = 0;

  my $fname = shift @ARGV;
  my $totalsize = shift @ARGV;

  my $tag = get_tag($fname);

  my $scale = $ENV{RCLPLOT_PARAM_SCALE} || 0;
  $tag /= 10**$ENV{RCLPLOT_PARAM_SCALE};

  while (<>) {
    chomp;
    $size = $_;
    while ($size > $y[$y]) {
      my $fraction = $runningtotal/$totalsize;
      print "$tag\t$y\t$fraction\t$y[$y]\n";
      $y++;
    }
    $runningtotal += $size;
  }
  while ($y < @y) {
    my $fraction = $runningtotal/$totalsize;
    print "$tag\t$y\t$fraction\t$y[$y]\n";
    $y++;
  }

}


sub attach_annotation {
  my ($h, $v) = @_;
  my %h = ();
  for (my $i=0; $i<@$h; $i++) {
    $h{$h->[$i]} = $v->[$i];
  }
  return \%h;
}

sub read_node_annotation_table {
  my $fntable = shift;
  open(CT, "<$fntable") || die "Cannot open annotation table $fntable\n";
  my $header = <CT>;
  chomp $header;
  my @header = split "\t", $header;
  die "No leading tab $header" unless $header[0] eq "";
  shift @header;
  my %dfannot = ();
  while (<CT>) {
    chomp;
    my @F = split "\t";
    my $cell = shift @F;
    die "Cardinal error $.\n" unless @F == @header;
    $dfannot{$cell} = attach_annotation(\@header, [@F]);
  } close(CT);
  return (\%dfannot, \@header);
}

sub read_tab {
  my $tabfile = shift;
  open(TAB, "<$tabfile") || die "No tab\n";
  my %tab = ();
  while (<TAB>) {
    chomp;
    my @F = split "\t";
    die "Tab error [$_]\n" unless @F == 2;
    $tab{$F[0]} = $F[1];
  }
  return \%tab;
}

sub get_median {
  my $ref = shift;
  my $N = @$ref;
  if (!$N) {
    print STDERR "## warning empty array for median\n" if $::median_warning++ < 10;
    return 0;
  }
  my $i = int(($N+1)/2);
  $i = $N-1 if $i >= $N;
  return $ref->[$i];
}


  # Used for sorting nodes in dendrogram to ensure same sorting as in input.
  # e.g. A   B_A   B_B_A   B_B_C  ..... ....   AB_A
  # TODO revisit this make it less fuzzy.
sub bycluskey {
  return length($a) <=> length($b) || $a cmp $b;
}


sub hm_newick {
  my ($lim, $depth, $key, $branchlength) = @_;
          # Dangersign NewickSelection.
          # Below is brittle. The idea is that we never stop at
          # internal nodes, but only filter on size once at a leaf node. This
          # ensures (with clearest/simplest/tightest coupling) that our
          # dendrogram matches the table filtering.  However, the
          # implementation needs to be more explicit about this (e.g.
          # is_internal).  The hm_tree, hm_nodes data structures need more
          # scrutitinising.
  my @children = grep { !defined($hm_nodes{$_}{size}) || $hm_nodes{$_}{size} >= $lim } sort bycluskey keys %{$hm_tree{$key}};
  if (!@children) {
    die "No info for childless node $key\n" unless defined($hm_nodes{$key});
    die "Leaf key error for $key $hm_order\n" if defined($hm_nodes{$key}) && !defined($hm_nodes{$key}{ival});
    my ($ival, $level, $type) = map { $hm_nodes{$key}{$_} } qw(ival level type);
    my $up = $ival - $branchlength;
# print STDERR "l\t$depth\t$key\t$branchlength\t$ival\t$up\n";
    return 'x' . sprintf("%04d", $hm_order++) . ".$key" . ':' . $up;
  }
  else {
              # TODO: check ival/pival, definition of up.
    my $ival = $branchlength;
    if (defined($hm_nodes{$key}) && !defined($hm_nodes{$key}{ival})) {
      my @huh = keys %{$hm_nodes{$key}};
      die "Internal key error for $key $hm_order [@huh]\n"
    }
    $ival = $hm_nodes{$key}{ival} if defined($hm_nodes{$key});
    my $up = $ival - $branchlength;
# print STDERR "i\t$depth\t$key\t$branchlength\t$ival\t$up\n";
    my @newick = ();
    for my $child (@children) {
      my $nwk = hm_newick($lim, $depth+1, $child, $branchlength + $up);
      push @newick, $nwk;
    }
    return "(" . (join ",", @newick) . "):$up";
  }
}


sub read_partition_hierarchy {

  my $inputmode = shift;
  my $datafile = shift;
  my $dfannot  = shift;
  my $termlist = shift;
  my $tab = shift;
  my $fnbase   = shift;

  my $lim = $ENV{RCLPLOT_HEAT_LIMIT} || 80;
  my $restful = defined($ENV{RCLPLOT_HEAT_NOREST}) ? 0 : 1;

  open(CLS, "<$datafile") || die "No partition input file\n";
  my %cls = ();
  my $toplevel = 'root';

  my %dfuniverse = ();
  my $NU = keys %$dfannot;
  print STDERR "-- computing term universe for $NU terms\n";

  for my $term (@$termlist) {
    for (sort keys %$dfannot ) {
      die "No $term score for $_\n" unless defined($dfannot->{$_}{$term});
      $dfuniverse{$term} += $dfannot->{$_}{$term};
    }
  }

  open (GLORIOUS, ">$fnbase.sum.txt") || die "Cannot open frequency output table $fnbase.fq.txt\n";
  open (NWK, ">$fnbase.nwk")   || die "Cannot open Newick output table $fnbase.nwk\n";
  print STDERR "-- computing cluster/term aggregate scores \n";

  local $" = "\t";
  print GLORIOUS "Type\tJoinval\tSize\t@$termlist\n";

  my @bglevel = map { log($dfuniverse{$_}) / log(10) } @$termlist;
  my @bgsum   = map { $dfuniverse{$_} } @$termlist;

  print GLORIOUS "NA\tNA\t$NU\t@bgsum\n";

    # level, type(rest/cls), Ncls, Nmiss, elements
  while (<CLS>) {
    chomp;
    if ($. == 1 && $inputmode eq 'rclhm') {
      die "Header [$_] not matched\n" unless $_ eq "level\ttype\tjoinval\tN1\tN2\tnesting\tnodes";
      next;
    }
    my ($level, $type, $joinval, $N1, $N2, $nesting, $elems) = (1, 'cls', 0, 0, "A", "");
    if ($inputmode eq 'rclhm') {
      ($level, $type, $joinval, $N1, $N2, $nesting, $elems) = split "\t";
      $hm_nodes{$nesting}{ival}  = $joinval;
      $hm_nodes{$nesting}{level} = $level;
      $hm_nodes{$nesting}{type}  = $type;
      $hm_nodes{$nesting}{size}  = $N2;
      $hm_nodes{$nesting}{print} = 1;
      $hm_tree{$nesting} = {};
                                                # _A are the residual clusters; their joinval is the
                                                # joinval for the parent cluster (which we need).
                                                # Hence we need a _A residual cluster for
                                                # every non-leaf RCL cluster, even if it is empty.
      while ($nesting =~ s/^(.*)(_.+?)$/$1/) {
        $hm_tree{$1}{"$1$2"} = 1;               # Dangersign format dependency.
        if ($2 eq '_A') {                       # Dangersign even more so.
          $hm_nodes{$1}{ival}  = $joinval;
        }
        $nesting = $1;
      }
      $hm_tree{'root'}{$nesting} = 1;
    }
    elsif ($inputmode eq 'cls') {
      $elems = $_;
    }
    else {
      die "Pity\n";
    }
    my @elems = length($elems) ? map { $tab->{$_} } split /\s+/, $elems : ();
    die "Error cls input $.\n" unless $inputmode eq 'cls' || $N2 == @elems;

    $N2 = $N1 = @elems unless $N2;    # This is for mode 'cls'.
    next unless $N2 >= $lim;          # Dangersign; NewickSelection dependency.

    if (@elems <= 1) {
      print STDERR "## small cluster [$level $type $joinval $N1 $N2 $elems]\n";
    }

    my %df = ();
    for my $t (@$termlist) {
      my $sum = 0;
      for my $e (@elems) {
        next if $e eq 'dummy';
        die "No $e $t\n" unless defined($dfannot->{$e}{$t});
        $sum += $dfannot->{$e}{$t};
      }
      $df{$t}{sum} = $sum;
    }
    my @values_gl = map { $df{$_}{sum}} @$termlist;
    print GLORIOUS "$type\t$joinval\t$N2\t@values_gl\n";

  } close(CLS);

  if ($inputmode eq 'rclhm') {
    my $nwk = hm_newick($lim, 0, 'root', 0, 0);
    print NWK "($nwk)\n";
  }

  close(NWK);
  close(GLORIOUS);
}


sub read_cls {
  my $clsfile = shift;
  my $dfannot  = shift;
  my $termlist = shift;
  my $tab = shift;

  open(CLS, "<$clsfile") || die "No cls\n";
  my $header = <CLS>;
  chomp $header;
  my %cls = ();
  my $toplevel = 'root';
  die "Header error\n" unless $header eq "level\tsize\tjoinval\tpctloss\tup\tdown\tnesting\tnodes";
  while (<CLS>) {
    chomp;
    my @F = split "\t";
    my $nesting = $F[6];
    next if $nesting eq 'root';
    if ($nesting !~ /^$toplevel/) {
      print "\n";
      $toplevel = $nesting;
    }
    my @elems = map { $tab->{$_} } split /\s+/, $F[7];
    my %df = ();
    for my $t (@$termlist) {
      for my $e (@elems) {
        $df{$t} += $dfannot->{$e}{$t};
      }
    }
    my $N = @elems;
    my @o = sort { $df{$b} <=> $df{$a} } keys %df;
    my $result = join "\t", map { my $t = $o[$_]; "$t\t$df{$t}" } 0..1;
    my $score  = sprintf("%.2f", $df{$o[0]} / $df{$o[1]});
    $nesting =~ s/L\d+_//g;
    $nesting =~ s/::/:/g;
    my $margin = defined($ENV{RCL_MARGIN}) ? $ENV{RCL_MARGIN} : 40;
    printf "%-*s\t%s\t%s\n", $margin, $nesting, $score, $result;
  } close(CLS);
}


