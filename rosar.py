#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math 
import os.path as op
import numpy as np
import argparse

import multiprocessing
nproc = multiprocessing.cpu_count()

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
def merge_seqs(fis, orgs, fo):
    print "merging input files to %s" % fo
    seqs = []
    for i in range(0,len(orgs)):
        seq_it = SeqIO.parse(open(fis[i], "rU"), "fasta")
        seqs1 = [SeqRecord(rcd.seq, id = orgs[i] + "|" + rcd.id,
            description = '') for rcd in seq_it]
        seqs += seqs1
    fho = open(fo, "w")
    SeqIO.write(seqs, fho, "fasta")
    fho.close()

def pfam_scan(fm, fi, fo, nproc):
    print "running hmmscan, using %d processors" % nproc
    cmd = "time hmmscan --cpu %d -o %s.1.txt %s %s" % (nproc, fo, fm, fi) 
    print cmd
    os.system(cmd)
    cmd = "hmmc2htb.pl -i %s.1.txt -o %s.2.htb -m %s -s %s" % (fo, fo, fm, fi)
    os.system(cmd)
    cmd = "htb.qtile.pl -i %s.2.htb -o %s.3.htb" % (fo, fo)
    os.system(cmd)
    cmd = "htb.filter.pl -i %s.3.htb -l 10 -e 0.01 -o %s.4.htb" % (fo, fo)
    os.system(cmd)
    cmd = "cut -f2-4,6,7-9,11-13 %s.4.htb > %s" % (fo, fo)
    os.system(cmd)
from itertools import groupby
from operator import itemgetter
def extract_nbs(fi, fo):
    ary = np.genfromtxt(fi, names = True, dtype = None)
    ary_sorted = sorted(ary, key = lambda x: (x[0], x[1]))

    fho = open(fo, "w")
    print >>fho, "\t".join(["org", "gid", "doms", "conf"])
    for key, valuesiter in groupby(ary_sorted, key = itemgetter(0)):
        res = [v[4] for v in valuesiter if v[8] < 1e-5 and v[2]-v[1]+1 >= 20]
        (org, gid) = key.split("|", 1)
        if 'NB-ARC' in res:
            print >>fho, "%s\t%s\t%s\t%s" % (org, gid, ",".join(res), '')
    fho.close()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Identify, cluster and characterizie plant NBS-LRR genes'
    )
 #   parser.add_argument('program', type=str, help='progname')
    parser.add_argument(
        'seqfile', nargs = '+', help = 'protein sequence file(s)'
    )
    parser.add_argument(
        '--out', dest = "output", default = "test", help = 'output directory (default: "test")'
    )
    parser.add_argument(
        '--cpu', dest = "ncpu", default = nproc, help = 'number processors to use (default: all/%d)' % nproc)
    args = parser.parse_args()

    (fis, dirw) = (args.seqfile, args.output)
    if not op.exists(dirw): os.makedirs(dirw)
    orgs = []
    for fi in fis:
        if not os.access(fi, os.R_OK):
            print "no access to input file: %s" % fi
            sys.exit(1)
        orgs.append(op.splitext(op.basename(fi))[0])
    
    f_pfam = '/home/youngn/zhoup/Data/db/pfam_v29/Pfam-A.hmm'
    if not os.access(f_pfam, os.R_OK):
        print "no access to the Pfam-A.hmm file"
        sys.exit(1)
    
    print "%d input files detected" % len(fis)
    print "species to work on: %s" % "  ".join(orgs)
    print "output directory: %s" % dirw
    
    f01 = op.join(dirw, "01.fas")
#    merge_seqs(fis, orgs, f01)
    f11 = op.join(dirw, "11")
#    pfam_scan(f_pfam, f01, f11, nproc)
    extract_nbs('11.htb', '21.tbl')
