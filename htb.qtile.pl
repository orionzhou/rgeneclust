#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  htb.qtile.pl - tile HMM hits along each qry-seq (hmmscan)

=head1 SYNOPSIS
  
  htb.qtile.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Location;
use Htb;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
#pod2usage(2) if !$fi || !$fo;

if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $h = {};
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [split "\t"];
  next unless @$ps == 15;
  my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $hId, $hBeg, $hEnd, $hSrd, $hSize, 
    $e, $score, $qls, $hls) = @$ps;
  $qSrd eq "+" || die "$id qSrd: $qSrd\n";
  $hSrd eq "+" || die "$id hSrd: $hSrd\n";
  my ($rql, $rhl) = (locStr2Ary($qls), locStr2Ary($hls));
  my $ql = [ map {[$qBeg + $_->[0] - 1, $qBeg + $_->[1] - 1]} @$rql ]; 
  my $hl = [ map {[$hBeg + $_->[0] - 1, $hBeg + $_->[1] - 1]} @$rhl ]; 
  $h->{$qId} ||= [];
  push @{$h->{$qId}}, [$id, $qBeg, $qEnd, $qSize, 
    $hId, $hBeg, $hEnd, $hSize, $e, $score, $ql, $hl];
}
close $fhi;

print $fho join("\t", @HEAD_HTB)."\n";
for my $qId (sort(keys(%$h))) {
  my @ary = @{$h->{$qId}};
  my $nid = join("_", map {$_->[0]} @ary);

  my @es = map {$_->[8]} @ary;
  my @qlocs = map {[$_->[1], $_->[2]]} @ary;
  my $ref = tiling(\@qlocs, \@es, 1);
  for (@$ref) {
    my ($qb, $qe, $idx) = @$_;
    my ($id, $qBeg, $qEnd, $qSize, $hId, $hBeg, $hEnd, $hSize, 
      $e, $score, $ql, $hl) = @{$ary[$idx]};
    my ($hb, $he) = map {coordTransform($_, $ql, "+", $hl, "+")} ($qb, $qe);
    ($hb, $he) = ($he, $hb) if $hb > $he;
    print $fho join("\t", $id, $qId, $qb, $qe, "+", $qSize,
      $hId, $hb, $he, "+", $hSize, 
      $e, $score, '', '')."\n";
  }
}
#print STDERR scalar(keys(%$h))." rows picked\n";
close $fho;


__END__
