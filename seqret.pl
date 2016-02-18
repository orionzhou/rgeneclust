#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  seqret.pl - retrieve sequences with given intervals from a seq-db

=head1 SYNOPSIS
  
  seqret.pl [-help] [-db fasta-db] [-out output-fasta]
                    [-bed BED-file] [-str interval-string]

  Options:
    -h (--help)   brief help message
    -d (--db)     sequence database (fasta)
    -o (--out)    output file (can be 'stdout')
    -b (--bed)    BED file with genomic intervals
    -s (--str)    interval string (e.g.: chr5:21-50,chr8:40-60)

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------

use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Bio::DB::Fasta;
use Bio::SeqIO;

my ($fd, $fo, $fb, $str) = ('') x 4;
my ($fho);
my $help_flag;

#------------------------------ MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "db|d=s"   => \$fd,
  "out|o=s"  => \$fo,
  "bed|b=s"  => \$fb,
  "str|s=s"  => \$str,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fd;
pod2usage(2) if !$fb && !$str;

my $db = Bio::DB::Fasta->new($fd);

if ($fo eq "" || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "Can't open file $fo for writing: $!\n";
}

my $seqHO = Bio::SeqIO->new(-fh=>$fho, -format=>'fasta');
my $cnt = 0;

if($fb && -s $fb) {
  open(my $fhb, "<$fb") or die "cannot read $fb\n";
  while(<$fhb>) {
    chomp;
    my ($seqid, $beg, $end) = split "\t";
    my $id;
    if(!defined($beg) || !defined($end)) {
      $beg = 1;
      $end = $db->length($seqid);
      $id = $seqid;
      defined $end || die "$id not in db\n";
    } else {
      $beg += 1; # 0-based coordinate
      $id = join("-", $seqid, $beg, $end);
    } 
    $beg <= $end || die "loc error in $fb\n$seqid:$beg-$end\n";
    
    my $seq = $db->seq($seqid, $beg, $end);
    defined $seq || die "$id not in db\n";
    $seqHO->write_seq( Bio::Seq->new(-id=>$id, -seq=>$seq) );
    $cnt ++;
  }
  close $fhb;
}

if($str) {
  my @ps = split(",", $str);
  for (@ps) {
    my ($seqid, $beg, $end);
    if(/^([\w\-]+)\:(\d+)\-(\d+)$/) {
      ($seqid, $beg, $end) = ($1, $2, $3);
    } elsif(/^([\w\-]+)\:(\d+)$/) {
      ($seqid, $beg) = ($1, $2);
      $end = $db->length($seqid);
    } elsif(/^([\w\-]+)$/) {
      $seqid = $1;
      ($beg, $end) = (1, $db->length($seqid));
    } else {
      die "unknown locstring: $str\n";
    }
    my $id = join("-", $seqid, $beg, $end);
    my $seq = $db->seq($seqid, $beg, $end);
    $seqHO->write_seq( Bio::Seq->new(-id=>$id, -seq=>$seq) );
    $cnt ++;
  }
}
$seqHO->close();
printf "  %4d sequences extracted\n", $cnt;

exit 0;



