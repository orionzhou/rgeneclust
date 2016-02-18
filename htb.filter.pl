#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  htb.filter.pl - filter Htb records 

=head1 SYNOPSIS
  
  htb.filter.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -e (--evalue) max evalue (default: 1) 
    -s (--score)  min score (default: 0) 
    -l (--len)    min qry/tgt length (default: 1)
    -q (--qcov)   min qry coverage (default: 0) 
    -t (--tcov)   min tgt coverage (default: 0) 

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

my ($fi, $fo) = ('') x 2;
my ($max_e, $min_score, $min_len, $min_tcov, $min_qcov) = (1, 0, 1, 0, 0);
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "evalue|e=f"  => \$max_e,
  "score|s=f"   => \$min_score,
  "len|l=i"     => \$min_len,
  "tcov|t=f"    => \$min_tcov,
  "qcov|q=f"    => \$min_qcov,
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

print $fho join("\t", @HEAD_HTB)."\n";

my $cnt = 0;
while( <$fhi> ) {
  chomp;
  next if /(^id)|(^\#)|(^\s*$)/;
  my $ps = [split("\t", $_, -1)];
  next unless @$ps == 15;
  my ($id, $qId, $qBeg, $qEnd, $qSrd, $qSize,
    $hId, $hBeg, $hEnd, $hSrd, $hSize, 
    $e, $score, $qlS, $hls) = @$ps;
  $score >= $min_score || next;
  $e < $max_e || next;
  ($hEnd - $hBeg + 1) >= $min_len || next;
  ($qEnd - $qBeg + 1) >= $min_len || next;
  ($hEnd - $hBeg + 1) / $hSize >= $min_tcov || next;
  ($qEnd - $qBeg + 1) / $qSize >= $min_qcov || next;
  print $fho join("\t", @$ps)."\n";
  $cnt ++;
}
print STDERR "$cnt rows passed filter\n";
close $fhi;
close $fho;


__END__
