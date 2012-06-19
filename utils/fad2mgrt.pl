#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts('d', \%opts);
die("Usage: fad2mgrt.pl [-d] <in.fad>\n") if (@ARGV == 0 && -t STDIN);

my (@pos, @seq, $id);
$id = 1;
while (<>) {
  my @t = split;
  push(@pos, $t[0]);
  for (0 .. length($t[1]) - 1) {
	$seq[$_] .= substr($t[1], $_, 1);
  }
}

my $n_cases = int(@seq/2);
my $n_ctrls = @seq - $n_cases;
$n_cases *=2, $n_ctrls *= 2 if (defined($opts{d}));
print "$n_cases $n_ctrls ", scalar(@pos), "\n";
print join("\n", @pos), "\n";
for (0 .. $#seq) {
  if (defined($opts{d})) {
	$seq[$_] =~ tr/2/U/;
	print $seq[$_], "\n";
  }
  print $seq[$_], "\n";
}
