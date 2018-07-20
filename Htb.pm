package Htb;
use strict;
use Common;
use Location; 
use Data::Dumper;
use List::Util qw/min max sum/;
use POSIX qw/ceil floor/;
use List::MoreUtils qw/first_index last_index insert_after apply indexes pairwise zip uniq/;
use vars qw/$VERSION @ISA @EXPORT @EXPORT_OK/;
require Exporter;
@ISA = qw/Exporter AutoLoader/;
@EXPORT = qw/@HEAD_HTB/;
@EXPORT_OK = qw//;

our @HEAD_HTB  = qw/id qid qbeg qend qsrd qsize hid hbeg hend hsrd hsize e score qloc hloc/;

my ($id, $qid, $qbeg, $qend, $qsrd, $qsize,
  $hid, $hbeg, $hend, $hsrd, $hsize, $e, $score, $qls, $hls) = ();




1;
__END__
