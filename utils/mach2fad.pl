#!/usr/bin/perl -w

use strict;
use warnings;

my $n = 0;
my @site;
while (<>) {
  my @t = split;
  $_ = <>;
  my @s = split;
  my $str = '';
  my $l = length($t[2]);
  $t[2] =~ tr/12/01/; $s[2] =~ tr/12/01/;
  if ($n == 0) {
	$site[$_] = '' for (0 .. $l-1);
  }
  for (0 .. $l-1) {
	$site[$_] .= substr($t[2], $_, 1) . substr($s[2], $_, 1);
  }
  ++$n;
}

for (0 .. $#site) {
  print "$_\t", $site[$_], "\n";
}
