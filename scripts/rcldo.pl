#!/usr/bin/perl

use strict;
use warnings;

my $mode = shift || die "Need mode";


if ($mode eq 'cumgra') {
  cumgra();
}
elsif ($mode eq 'distwrangle') {
  distwrangle();
}
elsif ($mode eq 'clstag') {
  die "Need file name\n" unless @ARGV == 1;
  my $tag = get_tag(@ARGV);
  my $scale = $ENV{RCLPLOT_PARAM_SCALE} || 0;
  $tag /= 10 ** $scale;
  print "$tag\n";
}
elsif ($mode eq 'heatannot' || $mode eq 'heatannotcls') {
  die "Need <annotationfile> <partitionhierarchyfile> <tabfile>\n" unless @ARGV == 3;
  my ($ct, $celltypes) = read_node_annotation_table($ARGV[0]);
  my $tab = read_tab($ARGV[2]);
  my $dispatch = $mode eq 'heatannot' ? 'rclph' : 'cls';
  read_partition_hierarchy($dispatch, $ARGV[1], $ct, $celltypes, $tab);
}
else {
  die "Unknown mode $mode\n";
}


sub get_tag {

  my $fname = shift;
  die "No file name supplied" unless $fname;
  $fname =~ /[rI](\d{2,})\b/ || die "Cannot find resolution parameter from file [$fname]\n";
  return $1;
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
  my %ct = ();
  while (<CT>) {
    chomp;
    my @F = split "\t";
    my $cell = shift @F;
    die "Cardinal error $.\n" unless @F == @header;
    $ct{$cell} = attach_annotation(\@header, [@F]);
  } close(CT);
  return (\%ct, \@header);
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


sub read_partition_hierarchy {

  my $inputmode = shift;
  my $datafile = shift;
  my $ct  = shift;
  my $celltypes = shift;
  my $tab = shift;
  my $lim = $ENV{RCLPLOT_HEAT_LIMIT} || 10;

  open(CLS, "<$datafile") || die "No partition input file\n";
  my %cls = ();
  my $toplevel = 'root';

  local $" = "\t";
  print "Type\tSize\t@$celltypes\n";

    # level, type(rest/cls), Ncls, Nmiss, elements
  while (<CLS>) {
    chomp;
    my ($level, $type, $N1, $N2, $elems) = (1, 'cls', 0, 0, "");
    if ($inputmode eq 'rclph') {
      ($level, $type, $N1, $N2, $elems) = split "\t";
    }
    elsif ($inputmode eq 'cls') {
      $elems = $_;
    }
    else {
      die "Pity\n";
    }
    my @elems = $elems ? map { $tab->{$_} } split /\s+/, $elems : ();
    $N2 = $N1 = @elems unless $N2;    # This is for mode 'cls'.
    my %df = ();
    for my $t (@$celltypes) {
      $df{$t} = 0;
      for my $e (@elems) {
        next if $e eq 'dummy';
die "No $e $t\n" unless defined($ct->{$e}{$t});
        $df{$t} += $ct->{$e}{$t};
      }
    }
    my @values = map { $df{$_} } @$celltypes;
    print "$type\t$N2\t@values\n" if $N2 >= $lim;
  } close(CLS);
}


sub read_cls {
  my $clsfile = shift;
  my $ct  = shift;
  my $celltypes = shift;
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
    for my $t (@$celltypes) {
      for my $e (@elems) {
        $df{$t} += $ct->{$e}{$t};
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


