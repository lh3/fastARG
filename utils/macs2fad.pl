#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('d', \%opts);
die("Usage: macs2fad.pl <in.macs>\n") if (@ARGV == 0 && -t STDIN);

my @conv = (0, 2, 2, 1);
my ($len, $flag, $pos, $is_d, $last_coor);
$is_d = defined($opts{d})? 1 : 0;
$last_coor = -1;
while (<>) {
  if (/^COMMAND.*macs\s+\d+\s+(\d+)/) {
	$len = $1; $flag = 0; $pos = 0;
  } elsif (/^BEGIN_SEGREGATING_SITES/) {
	$flag = 1;
  } elsif ($flag && /^\d+\s+(\S+)\s+(\S+)/) {
	my $coor = int($1 * $len + 0.5) + 1;
	++$coor if ($coor == $last_coor);
	if ($is_d) {
	  $_ = $2;
	  s/(\d)(\d)/$conv[$1<<1|$2]/ge;
	  print "$coor\t", $_, "\n";
	} else {
	  print "$coor\t", $2, "\n";
	}
	$last_coor = $coor;
  } elsif (/^END_SEGREGATING_SITES/) {
	$flag = 0;
  }
}
