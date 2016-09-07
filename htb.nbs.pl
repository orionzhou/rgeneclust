#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  htb.nbs.pl - 

=head1 SYNOPSIS
  
  htb.nbs.pl [-help] [-in input-file] [-out output-file]

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
use Common;
use Htb;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fo;

my $cary = [ 
  ['P-loop', 26, 34],
  ['Kin-2' , 105, 112],
  ['RNBS-B', 132, 138],
  ['GLPL'  , 193, 199]
];
open(my $fho, ">$fo") || die "cannot write $fo\n";
print $fho join("\t", qw/id size doms e domstr tag beg end/)."\n";

my $ti = readTable(-in=>$fi, -header=>1);
for my $i (1..$ti->nofRow) {
  my ($id, $size, $doms, $e, $qb, $qe, $hb, $he, $qlocs, $hlocs) = $ti->row($i-1);
  my ($rql, $rhl) = (locStr2Ary($qlocs), locStr2Ary($hlocs));
  my $ql = [ map {[$qb + $_->[0] - 1, $qb + $_->[1] - 1]} @$rql ]; 
  my $hl = [ map {[$hb + $_->[0] - 1, $hb + $_->[1] - 1]} @$rhl ]; 
  my $h = {};
  my @ary;
  for (@$cary) {
    my ($dom, $db, $de) = @$_;
    my $shl = [[$db, $de]];
    my ($lo, $olen) = posOvlp($shl, $hl);
    next if $olen / ($de-$db+1) < 0.6;
    ($db, $de) = @{$lo->[0]};
    my ($pb, $pe) = map {coordTransform($_, $hl, "+", $ql, "+")} ($db, $de);
    $h->{$dom} = [$pb, $pe, $db, $de];
    push @ary, [$dom, $pb, $pe, $db, $de];
  }
  my ($tag, $tb, $te) = ("") x 3;
  if(exists($h->{"P-loop"}) && exists($h->{"GLPL"})) {
    $tag = 1;
    ($tb, $te) = ($h->{"P-loop"}->[0], $h->{"GLPL"}->[1]);
    if($tb >= $te) {
      print "$id suck\n".Dumper($ql, $hl, $h);
    }
  }
  my $str = join(" ", map {"$_->[0]\[$_->[1]-$_->[2]\]"} @ary);
  print $fho join("\t", $id, $size, $doms, $e, $str, $tag, $tb, $te)."\n"; 
}
close $fho;


__END__
