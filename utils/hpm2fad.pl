#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('d', \%opts);
die("Usage: hpm2fad.pl [-d] <in.hpm>\n") if (-t STDIN && @ARGV == 0);
my ($n, $m) = (0, 0);
while (<>) {
  my @t = split;
  if ($n == 0) {
	$n = @t - 1;
  } else {
	$_ = <>;
	my @s = split;
	my $str = '';
	if (!defined($opts{d})) {
	  for (1 .. $#t) {
		$str .= ($t[$_]-1) . ($s[$_]-1);
	  }
	} else {
	  for (1 .. $#t) {
		$str .= ($t[$_] == $s[$_])? $t[$_]-1 : 2;
	  }
	}
	print "$t[0]\t$str\n";
  }
}
