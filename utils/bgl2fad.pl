#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('d', \%opts);
die("Usage: bgl2fad.pl [-d] <in.bgl>\n") if (-t STDIN && @ARGV == 0);
my ($n, $m) = (0, 0);
while (<>) {
  my @t = split;
  if (/^M/) {
	shift(@t);
	my $pos = shift(@t);
	my $str = '';
	if (defined($opts{d})) {
	  for (my $i = 0; $i < @t; $i += 2) {
		if ($t[$i] eq $t[$i+1]) {
		  $str .= ($t[$i] eq 'A')? 0 : 1;
		} else {
		  $str .= 2;
		}
	  }
	} else {
	  $str = join("", @t);
	  $str =~ tr/AC/01/;
	}
	print "$pos $str\n";
  }
}
