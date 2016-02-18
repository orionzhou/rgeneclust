#!/usr/bin/perl -w
#
# POD documentation
#---------------------------------------------------------------------------
=pod BEGIN
  
=head1 NAME
  
  hmmc2htb.pl - convert hmmscan output to Htb format

=head1 SYNOPSIS
  
  hmmscan2htb.pl [-help] [-in input-file] [-out output-file]

  Options:
    -h (--help)   brief help message
    -i (--in)     input file
    -o (--out)    output file
    -s (--seq)    protein fasta file
    -m (--hmm)    hmm file

=cut
  
#### END of POD documentation.
#---------------------------------------------------------------------------


use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Common;
use Location;
use Htb;
use List::Util qw/min max sum/;
use List::MoreUtils qw/first_index first_value insert_after apply indexes pairwise zip uniq/;

my ($fi, $fo) = ('', '');
my ($fs, $fm) = ('', '');
my ($fhi, $fho);
my $help_flag;

#--------------------------------- MAIN -----------------------------------#
GetOptions(
  "help|h"   => \$help_flag,
  "in|i=s"   => \$fi,
  "out|o=s"  => \$fo,
  "seq|s=s"  => \$fs,
  "hmm|m=s"  => \$fm,
) or pod2usage(2);
pod2usage(1) if $help_flag;
pod2usage(2) if !$fi || !$fs || !$fm;

if ($fi eq '' || $fi eq "stdin" || $fi eq "-") {
  $fhi = \*STDIN;
} else {
  open ($fhi, $fi) || die "cannot read $fi\n";
}

if ($fo eq '' || $fo eq "stdout" || $fo eq "-") {
  $fho = \*STDOUT;
} else {
  open ($fho, ">$fo") || die "cannot write $fo\n";
}

if(-s "$fs.index" && (-M "$fs.index") > (-M $fs)) {
  runCmd("rm -rf $fs.index");
}
my $db = Bio::DB::Fasta->new($fs);

-s "$fm.info" || runCmd("hmm.stat.pl -i $fm -o $fm.info");
my $t = readTable(-in => "$fm.info", -header => 1);
my $hm = { map {$t->elm($_, "id") => $t->elm($_, "len")} 0..$t->lastRow };

my $in = Bio::SearchIO->new(-fh => $fhi, -format => "hmmsearch3");
print $fho join("\t", @HEAD_HTB)."\n";

my $id = 1;
while( my $res = $in->next_result ) {
  while( my $hit = $res->next_hit ) {
    while( my $hsp = $hit->next_hsp ) {
      my $qid = $res->query_name;
      my $hid = $hit->name;

      exists $hm->{$hid} || die "no size for $hid\n";
      my $hsize = $hm->{$hid};
      my $qsize = $db->length($qid);
      my ($qbeg, $qend) = ($hsp->start("hit"), $hsp->end("hit"));
      my ($hbeg, $hend) = ($hsp->start("query"), $hsp->end("query"));
      my ($alnH, $alnQ) = ($hsp->query_string, $hsp->hit_string);
      my $alnP = $hsp->{"PP_SEQ"};
      my ($e, $score) = ($hsp->{"IEVALUE"}, $hsp->score());
      my ($alnQ_n, $alnH_n, $locQ, $locH) = 
        check_parse_aln($qid, $qbeg, $qend, $hid, $hbeg, $hend, $alnQ, $alnH);
      my ($locQs, $locHs) = map {locAry2Str($_)} ($locQ, $locH);
      print $fho join("\t", $id++, $qid, $qbeg, $qend, "+", $qsize, 
        $hid, $hbeg, $hend, "+", $hsize, 
        $e, $score, $locQs, $locHs)."\n";
    }
  }
}
close $fhi;
close $fho;

sub parse_aln {
  my ($alnQ, $alnT) = @_;
  
  my @poss;
  my @locmap;
  my ($posA, $posQ, $posT) = (0, 0, 0);
  for my $i (1..length($alnQ)) {
    my $chQ = substr($alnQ, $i-1, 1);
    my $chT = substr($alnT, $i-1, 1);
    if($chQ !~ /[\.\-]/ || $chT !~ /[\.\-]/) {
      push @poss, $i;
      if($chQ =~ /[\.\-]/ && $chT !~ /[\.\-]/) {
        push @locmap, [++$posA, '', ++$posT];
      } elsif($chQ !~ /[\.\-]/ && $chT =~ /[\.\-]/) {
        push @locmap, [++$posA, ++$posQ, ''];
      } else {
        push @locmap, [++$posA, ++$posQ, ++$posT];
      }
    }
  }

  my ($lenA, $lenQ, $lenT) = ($posA, $posQ, $posT);
  my $gapQ = grep {$_->[1] eq ''} @locmap;
  my $gapT = grep {$_->[2] eq ''} @locmap;

  my $alnQ_n = join("", map {substr($alnQ, $_-1, 1)} @poss);
  my $alnT_n = join("", map {substr($alnT, $_-1, 1)} @poss);

  my @idxs_ins = indexes {$_->[1] eq '' || $_->[2] eq ''} @locmap;
  my ($idxs_nogap) = posDiff([[0, @locmap-1]], [map {[$_, $_]} @idxs_ins]);
  my ($locA, $locQ, $locT) = ([], [], []);
  for (@$idxs_nogap) {
      my ($idxB, $idxE) = @$_;
      push @$locA, [$locmap[$idxB]->[0], $locmap[$idxE]->[0]];
      push @$locQ, [$locmap[$idxB]->[1], $locmap[$idxE]->[1]];
      push @$locT, [$locmap[$idxB]->[2], $locmap[$idxE]->[2]];
  }

  return ($alnQ_n, $alnT_n, $lenA, $lenQ, $lenT, $gapQ, $gapT, $locA, $locQ, $locT);
}
sub check_parse_aln {
  my ($idQ, $begQ, $endQ, $idH, $begH, $endH, $alnQ, $alnH) = @_;
  
  my ($lenAQ, $lenAH) = map {length($_)} ($alnQ, $alnH);
  $lenAQ == $lenAH || die("$idQ $idH:$begH-$endH length not equal [alnQ:$lenAQ <> alnH:$lenAH");

  my ($alnQ_n, $alnH_n, $lenA, $lenQ, $lenH, $gapQ, $gapH, $locA, $locQ, $locH) = parse_aln($alnQ, $alnH);
  $lenA == $lenQ + $gapQ || die("$idQ $idH:$begH-$endH lenA[$lenA] != lenQ[$lenQ] + gapQ[$gapQ]\n$alnQ_n\n$alnH_n");
  $lenA == $lenH + $gapH || die("$idQ $idH:$begH-$endH lenA[$lenA] != lenH[$lenH] + gapH[$gapH]\n$alnQ_n\n$alnH_n");
  $lenQ == locAryLen($locQ) + $gapH || die("$idQ $idH:$begH-$endH lenQ[$lenQ] ".locAry2Str($locQ)." gapH[$gapH]\n$alnQ_n\n$alnH_n");
  $lenQ == $endQ - $begQ + 1 || die("$idQ $idH:$begH-$endH lenQ[$lenQ] != $begQ-$endQ");
  $lenH == locAryLen($locH) + $gapQ || die("$idQ $idH:$begH-$endH lenH[$lenH] ".locAry2Str($locH)." gapQ[$gapQ]\n$alnQ_n\n$alnH_n");
  $lenH == $endH - $begH + 1 || die("$idQ $idH:$begH-$endH lenH[$lenH] != $begH-$endH");

# $locQ = [ map {[$_->[0] + $begQ - 1, $_->[1] + $begQ - 1]} @$locQ ];
#  $locH = [ map {[$_->[0] + $begH - 1, $_->[1] + $begH - 1]} @$locH ];

  return ($alnQ_n, $alnH_n, $locQ, $locH);
}

__END__
