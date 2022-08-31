#!/usr/bin/perl -an

# This prints the links in the input graph between two clusters that are merged,
# and nothing else (e.g. no intra-cluster links).
# only use is by rcl-mkimg.sh
  
use strict;
use warnings;

BEGIN { $::rclfile = shift; %::sets = (); }
next if $. == 1;
my ($i, $x, $y, $xid, $yid, $val, $xcid, $ycid, $xcsz, $ycsz, $ignore) = @F;

my $xnew = 0;
my $ynew = 0;

if ($xcsz == 1) {
  push @{$::sets{$xcid}}, $xcid;
  $xnew = 1;
}
if ($ycsz == 1) {
  push @{$::sets{$ycid}}, $ycid;
  $ynew = 1;
}

my $of = sprintf "link/out.%03d", $i;

if ($xid == 87 || $yid == 87) {
  local $" = ' ';
  print STDERR "$i xset @{$::sets{$xcid}}\n";
  print STDERR "$i yset @{$::sets{$ycid}}\n";
  local $" = ',';
  my @xset = sort { $a <=> $b } @{$::sets{$xcid}};
  my @yset = sort { $a <=> $b } @{$::sets{$ycid}};
  print STDERR qq{mcxsubs -imx $::rclfile "dom(c, i(@xset)), dom(r, i(@yset))" | mcxdump\n};
}
local $" = ',';
system qq{echo "$x\t$y\t$val\t$xnew\t$ynew" > $of};
# print "\n";
system qq{mcxsubs -imx $::rclfile "dom(c, i(@{$::sets{$xcid}})), dom(r, i(@{$::sets{$ycid}}))" | mcxdump >> $of};
# print "\n";
# print "$i\t$val\t(@{$::sets{$xcid}})\t(@{$::sets{$ycid}})\n";

if ($xcid < $ycid) {
   push @{$::sets{$xcid}}, @{$::sets{$ycid}};
}
else {
   push @{$::sets{$ycid}}, @{$::sets{$xcid}};
}

