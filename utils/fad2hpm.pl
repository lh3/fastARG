#!/usr/bin/perl -w

use strict;
use warnings;

die("Usage: fad2hpm.pl <in.fad>\n") if (-t STDIN && @ARGV == 0);
my ($m, $n) = (0, 0);
while (<>) {
  if (/(\d+)\s(\S+)/) {
	my @t = split;
	if ($n == 0) { # the first line
	  $n = length($t[1]);
	  print "Id ", join(" ", (1..$n)), "\n";
	}
	$_ = $t[1]; tr/0123/1210/; print join(" ", $t[0], split("", $_)), "\n";
	$_ = $t[1]; tr/0123/1220/; print join(" ", $t[0], split("", $_)), "\n";
	++$m;
  }
}
