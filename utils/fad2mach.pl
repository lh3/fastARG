#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
die("Usage: fad2mach.pl <in.fad>\n") if (@ARGV == 0 && -t STDIN);

my (@pos, @seq, $id);
$id = 1;
while (<>) {
  my @t = split;
  push(@pos, "M $t[0]");
  for (0 .. length($t[1]) - 1) {
	$seq[$_] .= substr($t[1], $_, 1);
  }
}

my @conv = ("\t1/1", "\t2/2", "\t1/2");
print STDERR join("\n", @pos), "\n";
for (0 .. $#seq) {
  print "S_$_\tS_$_\t0\t0\tF";
  $_ = $seq[$_];
  s/(.)/$conv[$1]/eg;
  print $_, "\n";
}
