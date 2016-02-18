#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  hmm.stat.pl - 

=head1 SYNOPSIS
  
  hmm.stat.pl [-help] [-i input-file] [-o output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file (HMM)
    -o (--out)    output file (Tbl)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin";
use Common;
use Data::Dumper;
use File::Path qw/make_path remove_tree/;
use File::Basename;
use List::Util qw/min max sum/;

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

my $fho;
if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

sub check_gap {
  my ($fi) = @_;
  my $ai = Bio::AlignIO->new(-file=>"<$fi");
  my $gap = 0;
  while(my $aln = $ai->next_aln()) {
    for my $seq ($aln->each_seq()) {
      if($seq->seq =~ /[\-\.]/) {
        $gap = 1;
        last;
      }
    }
  }
  return $gap;
}
  
print $fho join("\t", qw/id acc accl nseq len gap conseq/)."\n";
my $lines = runCmd("hmmstat $fi", 2);
for (@$lines) {
  /(^\#)|(^\s*$)/ && next;
  my @ps = split " ";
  my ($id, $accl, $nseq, $len) = @ps[1,2,3,5];
  $accl = "" if $accl eq "-";
  my $acc = $accl;
  $acc =~ s/\.\d+$//;
  
#  my $seqlines = runCmd("hmmemit -c $fi", 2);
#  my $seq = join("", @$seqlines[1..@$seqlines-1]);

  print $fho join("\t", $id, $acc, $accl, $nseq, $len, '', '')."\n";
}
close $fho;
