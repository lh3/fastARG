#!/usr/bin/perl -w

use strict;
use warnings;

die("Usage: fad2bgl.pl <in.fad>\n") if (-t STDIN && @ARGV == 0);
my ($m, $n) = (0, 0);
my @conv = (" A A", " C C", " A C", " ? ?");
while (<>) {
  if (/(\d+)\s(\S+)/) {
	my @t = split;
	if ($n == 0) { # the first line
	  $n = length($t[1]);
	  print join(" ", "I", "id");
	  print " $_ $_" for (1 .. $n);
	  print "\n";
	}
	$_ = $t[1];
	s/(.)/$conv[$1]/eg;
	print join(" ", "M", $t[0]), $_, "\n";
	++$m;
  }
}
