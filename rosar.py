#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math 
import os.path as op
import numpy as np
import argparse
import usearch

import multiprocessing
nproc = multiprocessing.cpu_count()

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
def read_cfg(fc):
    (orgs, fis) = ([], [])
    fhc = open(fc, "r")
    for line in fhc:
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        (org, fi) = line.split(",")
        if not os.access(fi, os.R_OK):
            print "no access to input file: %s" % fi
            print os.access(fi, os.F_OK)
            sys.exit(1)
        orgs.append(org)
        fis.append(fi)
    fhc.close()
    return (orgs, fis)
def merge_seqs(fis, orgs, fo):
    print "merging input files to %s" % fo
    seqs = []
    for i in range(0,len(orgs)):
        handle = 0
        if (fis[i].endswith(".gz")):
            handle = gzip.open(fis[i], "rb")
        else:
            handle = open(fis[i], "rU")
        seq_it = SeqIO.parse(handle, "fasta")
        handle.close

        seqs1 = [SeqRecord(rcd.seq, id = orgs[i] + "|" + rcd.id,
            description = '') for rcd in seq_it]
        seqs += seqs1
    fho = open(fo, "w")
    SeqIO.write(seqs, fho, "fasta")
    fho.close()
def cds2pro(fi, fo):
    print "translating CDSs to proteins"
    seq_it = SeqIO.parse(open(fi, "rU"), "fasta")
    seqs = [SeqRecord(rcd.seq.translate(stop_symbol = ""), id = rcd.id,
        description = '') for rcd in seq_it]
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
def extract_nbs(fi, fo, fd_cds, fd_pro, fs_cds, fs_pro):
    ary = np.genfromtxt(fi, names = True, dtype = None)
    ary_sorted = sorted(ary, key = lambda x: (x[0], x[1]))

    cds_dict = SeqIO.index(fd_cds, "fasta")
    pro_dict = SeqIO.index(fd_pro, "fasta")

    fho = open(fo, "w")
    print >>fho, "\t".join(["id", "doms", "conf"])
    cds_rcds = []
    pro_rcds = []
    for key, valuesiter in groupby(ary_sorted, key = itemgetter(0)):
        res = [[v[1], v[2], v[4]] for v in valuesiter if v[8] < 1e-5 and v[2]-v[1]+1 >= 20]
        doms = [v[2] for v in res]
        if 'NB-ARC' in doms:
            print >>fho, "%s\t%s\t%s" % (key, ",".join(doms), '')
            cdss = [str(cds_dict[key].seq[(v[0]-1)*3:v[1]*3]) for v in res]
            cds_rcd = SeqRecord(Seq("".join(cdss)), id = key, description = "")
            cds_rcds.append(cds_rcd)
            pro_rcds.append(pro_dict[key])
    fho.close()
    SeqIO.write(cds_rcds, fs_cds, "fasta")
    SeqIO.write(pro_rcds, fs_pro, "fasta")
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Identify, cluster and characterizie plant NBS-LRR genes'
    )
 #   parser.add_argument('program', type=str, help='progname')
    parser.add_argument(
        'cfgfile', help = 'config file (a text file with species identifier followed by the absolute path of CDS fasta in each line)'
    )
    parser.add_argument(
        '--out', dest = "output", default = "test", help = 'output directory (default: "test")'
    )
    parser.add_argument(
        '--cpu', dest = "ncpu", default = nproc, help = 'number processors to use (default: all/%d)' % nproc)
    args = parser.parse_args()

    (fc, dirw) = (args.cfgfile, args.output)
    (orgs, fis) = read_cfg(fc)
    if not op.exists(dirw): os.makedirs(dirw)
    
    f_pfam = '/home/youngn/zhoup/Data/db/pfam_v29/Pfam-A.hmm'
    if not os.access(f_pfam, os.R_OK):
        print "no access to the Pfam-A.hmm file"
        sys.exit(1)
    
    cdir = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] = os.environ['PATH']+':'+cdir
    cwd = os.getcwd()
    os.chdir(dirw)
    
    print "%d input files detected" % len(fis)
    print "species to work on: %s" % "  ".join(orgs)
    print "output directory: %s" % dirw
    
#    merge_seqs(fis, orgs, "01.cds.fas")
#    cds2pro("01.cds.fas", "05.pro.fas")
#    pfam_scan(f_pfam, '05.pro.fas', '11', nproc)
#    os.system("ncoils.py 05.pro.fas")
#    extract_nbs("11.htb", "21.tbl", "01.cds.fas", "05.pro.fas", "23.nbs.cds.fas", "22.nbs.pro.fas")
    cmd = "usearch -cluster_fast %s -sort length -id %g -uc %s" % ('23.nbs.cds.fas', 0.8, '31.uc')
    os.system(cmd)
    usearch.usearch2tbl('31.uc', '32.tbl')
    sys.exit(1);
