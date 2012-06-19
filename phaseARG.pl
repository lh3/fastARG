#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Cwd qw/getcwd abs_path/;

&main;
exit;

sub main {
  my %opts = (r=>20, n=>40);
  getopts('r:n:', \%opts);
  die("Usage: phaseARG.pl <in.fad> <out.prefix>\n") if (@ARGV != 2);
  my ($n_rounds, $n_splits, %data);
  my ($fn, $prefix) = @ARGV;
  my $fastARG = &gwhich('fastARG');
  system("rm -f $prefix.*.fad");
  &arg_data_read($fn, \%data);
  $n_rounds = $opts{r};
  $opts{n} = $data{n} if ($opts{n} > $data{n});
  $n_splits = int($data{n} / $opts{n} + .499);
  &init_tmp_file(">$prefix.pileup", $data{m}, $data{n}, $n_rounds);
  my $fh = &gopen("+<$prefix.pileup");
  my @ida;
  $ida[$_] = 0 for (0 .. $data{n});
  my %globals = (PRE=>$prefix, DATA=>\%data, PROG=>$fastARG, M=>$data{m}, N=>$data{n}, IDA=>\@ida, NR=>$n_rounds,
				 FH=>$fh);
  for (my $r = 0; $r < $n_rounds; ++$r) {
	my @list = &gen_sub($data{n}, $n_splits);
	for my $x (@list) {
	  &arg_phase1(\%globals, $x, $r%2==1);
	}
  }
  close($fh);
  unlink("$prefix.tmp.fad");
  &pileup2fads($prefix, \%data, $n_rounds);
}

sub init_tmp_file {
  my ($fn, $m, $n, $n_rounds) = @_;
  my $fh = &gopen($fn);
  for (my $i = 0; $i < $n; ++$i) {
	for (my $j = 0; $j < $n_rounds; ++$j) {
	  print $fh "?"x$m, "\n";
	}
	print $fh "-"x$m, "\n";
  }
  close($fh) unless (ref($fn) eq 'GLOB');
}

sub pileup2fads {
  my ($prefix, $data, $n_rounds) = @_;
  my ($m, $n) = ($data->{m}, $data->{n});
  my @r;
  my @conv = ('00', '01', '10', '11');
  my $fh = &gopen("$prefix.pileup");
  for (my $j = 0; $j < $n_rounds; ++$j) {
	my $x = \%{$r[$j]};
	$x->{m} = $data->{m}; $x->{n} = $data->{n}; @{$x->{p}} = @{$data->{p}};
  }
  for (my $i = 0; $i < $n; ++$i) {
	for (my $j = 0; $j < $n_rounds; ++$j) {
	  $_ = <$fh>; chomp;
	  s/(\d)/$conv[$1]/eg;
	  $r[$j]{s}[$i] = $_;
	}
	<$fh>;
  }
  close($fh);
  for (my $j = 0; $j < $n_rounds; ++$j) {
	&arg_data_print($r[$j], sprintf(">$prefix.%.4d.fad", $j));
  }
}

sub arg_phase1 {
  my ($g, $list, $is_b) = @_;
  my @inds = split(',', $list);
  my (%sub, %rst);
  my $fn = "$g->{PRE}.tmp.fad";
  my $ida = $g->{IDA};
  &arg_data_sub($g->{DATA}, \@inds, \%sub);
  &arg_data_print(\%sub, ">$g->{PRE}.tmp.fad");
  my $b = $is_b? '-b' : '';
  &arg_data_read("$g->{PROG} build -Ls0 $b $fn 2>/dev/null | $g->{PROG} leafseq - |", \%rst);
  for (my $i = 0; $i < @inds; ++$i) {
	my $x = $inds[$i];
	my $offset = ($x * ($g->{NR} + 1) + $ida->[$x]) * ($g->{M} + 1);
	seek($g->{FH}, $offset, 0);
	my $y = $rst{s}[$i];
	$y =~ s/(\d)(\d)/$1<<1|$2/eg;
	print {$g->{FH}} $y;
  }
  ++$ida->[$_] for (@inds);
}

sub arg_data_read {
  my ($fn, $data) = @_;
  my $fh = &gopen($fn);
  %$data = ();
  my ($d, $pos) = (\@{$data->{s}}, \@{$data->{p}});
  while (<$fh>) {
	my @t = split;
	die("[arg_read_data] number of haplotypes is not an even number\n") if (length($t[1])%2 != 0);
	push(@$pos, $t[0]);
	my $l = length($t[1]);
	for (my $i = 0; $i < $l; $i += 2) {
	  my $id = int($i/2);
	  if (defined $d->[$id]) {
		$d->[$id] .= substr($t[1], $i, 2);
	  } else {
		$d->[$id] = substr($t[1], $i, 2);
	  }
	}
  }
  close($fh) unless (ref($fn) eq 'GLOB');
  $data->{m} = length($d->[0])/2;
  $data->{n} = @$d;
}

sub arg_data_print {
  my ($data, $fn) = @_;
  my $fh = &gopen($fn);
  my ($m, $n) = ($data->{m}, $data->{n});
  my ($d, $pos) = (\@{$data->{s}}, \@{$data->{p}});
  for my $x (0 .. $m-1) {
	my $str = '';
	$str .= substr($_, $x*2, 2) for (@$d);
	print $fh "$pos->[$x]\t", $str, "\n";
  }
  close($fh) unless (ref($fn) eq 'GLOB');
}

sub arg_data_sub {
  my ($data, $inds, $sub) = @_;
  %$sub = ();
  my ($d, $pos) = (\@{$data->{s}}, \@{$data->{p}});
  my $sd = \@{$sub->{s}};
  @{$sub->{p}} = @{$data->{p}};
  $sub->{m} = $data->{m};
  $sub->{n} = @$inds;
  $sd->[$_] = $d->[$inds->[$_]] for (0 .. @$inds-1);
}

sub gen_sub {
  my ($m, $n) = @_;
  die if ($m < $n);
  my $x = int($m/$n); # floor
  my $y = $m - $n * $x;
  my @a = &shuffle($m);
  my @b;
  if ($y) {
	push(@b, join(',', @a[$_*($x+1) .. ($_+1)*($x+1)-1])) for (0 .. $y-1);
  }
  my $z = $y * ($x + 1);
  push(@b, join(',', @a[$z+$_*$x .. $z+($_+1)*$x-1])) for (0 .. $n-$y-1);
  return @b;
}

sub shuffle {
  my $n = shift;
  my @a = (0 .. $n-1);
  for (my $i = $n - 1; $i > 0; --$i) {
	my $x = int(rand() * $i);
	$a[$i] ^= $a[$x]; $a[$x] ^= $a[$i]; $a[$i] ^= $a[$x];
  }
  return @a;
}

sub dirname {
  my $prog = shift;
  my $cwd = getcwd;
  return $cwd if ($prog !~ /\//);
  $prog =~ s/\/[^\s\/]+$//g;
  return $prog;
}

sub which {
  my $file = shift;
  my $path = (@_)? shift : $ENV{PATH};
  return if (!defined($path));
  foreach my $x (split(":", $path)) {
	$x =~ s/\/$//;
	return "$x/$file" if (-x "$x/$file" && -f "$x/$file");
  }
  return;
}

sub gwhich {
  my $progname = shift;
  my $addtional_path = shift if (@_);
  my $dirname = &dirname($0);
  my $tmp;
  chomp($dirname);
  if (-x $progname && -f $progname) {
	return abs_path($progname);
  } elsif (defined($dirname) && (-x "$dirname/$progname" && -f "$dirname/$progname")) {
	return abs_path("$dirname/$progname");
  } elsif (($tmp = &which($progname))) { # on the $PATH
	return $tmp;
  } else {
	warn("[gwhich] fail to find executable $progname anywhere.");
	return;
  }
}

sub gopen {
  my $f = shift;
  return $f if (ref($f) eq 'GLOB');
  if (!ref($f)) {
	my $fh;
	unless (open($fh, $f)) {
	  die("[gopen] fail to open file $f.\n");
	}
	return $fh;
  }
}
