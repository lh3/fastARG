#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

&main;

sub main{
  my $version = '0.1.3';
  &usage($version) if (@ARGV < 1);
  my $command = shift(@ARGV);
  my %func = (macs2fad=>\&macs2fad, arg2dot=>\&arg2dot, fad2hit=>\&fad2hit, hit2fad=>\&hit2fad,
			  bgl2fad=>\&bgl2fad, fad2bgl=>\&fad2bgl, fad2mach=>\&fad2mach, mach2fad=>\&mach2fad);
  die("Unknown command \"$command\".\n") if (!defined($func{$command}));
  &{$func{$command}}();
}

sub arg2dot {
  my %opts = (c=>11, p=>-1);
  getopts('c:ip:', \%opts);
  die(qq/Usage: fastARG.pl arg2dot [-i] [-p -1] <in.arg>\n/) if (-t STDIN && @ARGV == 0);
  my ($m, $n);
  print qq/
digraph {
  graph [ranksep=0.1, nodesep=0.2];
  node [margin=\"0.03,0.02\", style=\"setlinewidth(0.33),filled\", width=0.1, height=0.1, fontname="Helvetica-bold", fontsize=7, fillcolor="white", color="white"];
  edge [arrowhead=normal, arrowsize=0.33, style="setlinewidth(0.33)", fontsize=6];
/;
  my (@node, %leaves);
  my $pos = $opts{p};
  while (<>) {
	if (/^N\s(\d+)\s(\d+)/) {
	  $n = $1; $m = $2;
	} elsif (/^[CR]/) {
	  my @t = split;
	  if ($pos >= 0 && $pos < $m) {
		next unless ($pos >= $t[3] && $pos < $t[4]);
	  }
	  $leaves{$t[2]} = 1 if ($t[2] < $n);
	  if (!defined($node[$t[1]])) {
		my $attr = '';
		$attr .= 'label="" ' if ($t[1] >= $n && !defined($opts{i}));
		$attr .= 'fontcolor="white" fillcolor=' . (/^R/? qq/"red" / : qq/"black" /);
		chop($attr);
		print qq/  $t[1] [$attr];/ if ($attr);
		$node[$t[1]] = 1;
	  }
	  my $width = $t[5] / ($t[4] - $t[3]) * $opts{c} + 0.33;
	  $width = 10 if ($width > 10);
	  my $minlen = int($width + 0.67 + 0.499);
	  print qq/  $t[1] -> $t[2] [label="$t[3],$t[4]", style="setlinewidth($width)", minlen=$minlen];\n/;
	}
  }
  print qq/{ rank = same; /, join("; ", keys(%leaves)), "; };\n";
  print qq/}\n/;
}

sub macs2fad {
  my %opts = ();
  getopts('d', \%opts);
  die("Usage: fastARG.pl macs2fad <in.macs>\n") if (@ARGV == 0 && -t STDIN);
  my @conv = ('00', '22', '22', '11');
  my ($n, $len, $flag, $pos, $is_d, $last_coor);
  $is_d = defined($opts{d})? 1 : 0;
  $last_coor = -1;
  while (<>) {
	if (/^COMMAND.*macs\s+(\d+)\s+(\d+)/) {
	  $n = $1; $len = $2; $flag = 0; $pos = 0;
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
}

#############################
# format conversion for HIT #
#############################

sub fad2hit {
  die("Usage: fastARG.pl fad2hit <in.fad>\n") if (-t STDIN && @ARGV == 0);
  while (<>) {
	if (/(\d+)\s(\S+)/) {
	  my @t = split;
	  my ($s1, $s2) = ('', '');
	  for (my $i = 0; $i < length($t[1]); $i += 2) {
		$s1 .= substr($t[1], $i, 1);
		$s2 .= substr($t[1], $i+1, 1);
	  }
	  $s1 =~ tr/0123/1210/; $s2 =~ tr/0123/1220/;
	  print $s1, "\n", $s2, "\n";
	}
  }
}

sub hit2fad {
  my %opts = ();
  getopts('d', \%opts);
  die("Usage: fastARG.pl [-d] hit2fad <in.hit>\n") if (-t STDIN && @ARGV == 0);
  my ($n, $m) = (0, 0);
  while (<>) {
	my $h1 = $_;
	my $h2 = <>;
	chop($h1); chop($h2);
	$h1 =~ tr/12/01/; $h2 =~ tr/12/01/;
	my $str = '';
	if (!defined($opts{d})) {
	  for (0 .. length($h1)) {
		$str .= substr($h1, $_, 1) . substr($h2, $_, 1);
	  }
	} else {
	}
	print "$n\t$str\n";
	++$n;
  }
}

################################
# format conversion for BEAGLE #
################################

sub fad2bgl {
  die("Usage: fastARG fad2bgl <in.fad>\n") if (-t STDIN && @ARGV == 0);
  my ($m, $n) = (0, 0);
  my %conv = (00=>" A A", 11=>" C C", 22=>" A C", 33=>" ? ?");
  while (<>) {
	if (/(\d+)\s(\S+)/) {
	  my @t = split;
	  if ($n == 0) { # the first line
		$n = length($t[1]);
		print join(" ", "I", "id");
		print " $_ $_" for (1 .. $n/2);
		print "\n";
	  }
	  $_ = $t[1];
	  s/00/ A A/g; s/11/ C C/g; s/22/ A C/g; s/33/ ? ?/g;
	  print join(" ", "M", $t[0]), $_, "\n";
	  ++$m;
	}
  }
}

sub bgl2fad {
  my %opts = ();
  getopts('d', \%opts);
  die("Usage: fastARG bgl2fad [-d] <in.bgl>\n") if (-t STDIN && @ARGV == 0);
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
}

##############################
# format conversion for MACH #
##############################

sub fad2mach {
  my %opts = ();
  die("Usage: fastARG.pl fad2mach <in.fad>\n") if (@ARGV == 0 && -t STDIN);
  my (@pos, @seq);
  while (<>) {
	my @t = split;
	push(@pos, "M $t[0]");
	for (my $i = 0; $i < length($t[1]); $i += 2) {
	  $seq[$i/2] .= substr($t[1], $i, 1);
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
}

sub mach2fad {
  die("Usage: fastARG.pl mach2fad <in.fad>\n") if (@ARGV == 0 && -t STDIN);
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
}

#########
# Usage #
#########

sub usage
{
  my $version = shift;
  die qq(
Program: fastARG.pl (helper script for fastARG)
Version: $version
Contact: Heng Li <lh3\@sanger.ac.uk>

Usage:   fastARG.pl <command> <arguments> [...]

Command: macs2fad     convert MaCS output to fastARG input
         arg2dot      convert the arg format to the dot format for graphviz
         fad2hit      convert fad to HIT
         hit2fad      convert HIT to fad
         fad2bgl      convert fad to BEAGLE
         bgl2fad      convert BEAGLE to fad
         fad2mach     convert fad to MACH
         mach2fad     convert MACH to fad
\n);
}
