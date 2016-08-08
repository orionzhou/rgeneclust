#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  dna2pro.pl - translate a set of DNA sequences to protein (amino acid sequences)

=head1 SYNOPSIS
  
  dna2pro.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input (fasta) file (default: stdin)
    -o (--out)    output (fasta) file (default: stdout)

=cut
  
#### END of POD documentation.
#------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;

my ($fi, $fo) = ('') x 2;
my $help_flag;

#------------------------------ MAIN -----------------------------------#
GetOptions(
  "help|h"  => \$help_flag,
  "in|i=s"  => \$fi,
  "out|o=s" => \$fo,
) or pod2usage(2);
pod2usage(1) if $help_flag;

my ($fhi, $fho);
if ($fi eq "" || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, "<$fi") || die "cannot read $fi\n";
}

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

my $seqHI = Bio::SeqIO->new(-fh => $fhi, -format => 'fasta');
my $seqHO = Bio::SeqIO->new(-fh => $fho, -format => 'fasta');
while(my $dna = $seqHI->next_seq()) {
    my $seqpro = $dna->translate()->seq;
    $seqpro =~ s/\*$//;
    $seqpro =~ s/\*/X/;
    my $pro = Bio::Seq->new(-id => $dna->id, -seq => $seqpro);
    $seqHO->write_seq($pro);
}
$seqHI->close();
$seqHO->close();
close $fho;

exit 0;
