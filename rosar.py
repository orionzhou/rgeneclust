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
def read_tsv(fi, names = False, delimiter = "\t"):
    ary = np.genfromtxt(fi, names = True, dtype = None)
    aryo = []
    for i in range(len(ary)):
        row = []
        for j in range(len(ary[i])):
            val = ary[i][j]
            if type(ary[i][j]) is np.bytes_:
                val = str(ary[i][j], 'utf-8')
            row.append(val)
        aryo.append(row)
    return aryo
def read_cfg(fc):
    (orgs, fis, fgs) = ([], [], [])
    fhc = open(fc, "r")
    for line in fhc:
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        ps = line.split(",")
        (org, fi) = ps[0:2]
        if len(ps) > 2:
            fg = ps[2]
        else:
            fg = ''
        if not os.access(fi, os.R_OK):
            print("no access to input file: %s" % fi)
            print(os.access(fi, os.F_OK))
            sys.exit(1)
        orgs.append(org)
        fis.append(fi)
        fgs.append(fg)
    fhc.close()
    return (orgs, fis, fgs)
def merge_seqs(fis, orgs, fo):
    print("merging input files to %s" % fo)
    seqs = []
    for i in range(0,len(orgs)):
        handle = 0
        if (fis[i].endswith(".gz")):
            handle = gzip.open(fis[i], "rb")
        else:
            handle = open(fis[i], "rU")
        seq_it = SeqIO.parse(handle, "fasta")
        handle.close
        
        seqs1 = [SeqRecord(rcd.seq, id = orgs[i] + "|" + rcd.id, description = '') for rcd in seq_it]
        seqs += seqs1
    fho = open(fo, "w")
    SeqIO.write(seqs, fho, "fasta")
    fho.close()

def getDigits(num): 
    digit = 1
    while int(num / 10) >= 1:
        num = int(num/10)
        digit += 1
    return digit
def pfam_scan(fm, fi, fo, nproc):
    print("running hmmscan, using %d processors" % nproc)
    cmd = "hmmscan --cpu %s -o %s.1.txt %s %s" % (nproc, fo, fm, fi)
    os.system(cmd)
    cmd = "hmmc2htb.pl -i %s.1.txt -o %s.2.htb -m %s -s %s" % (fo, fo, fm, fi)
    os.system(cmd)
    cmd = "htb.qtile.pl -i %s.2.htb -o %s.3.htb" % (fo, fo)
    os.system(cmd)
    cmd = "htb.filter.pl -i %s.3.htb -l 10 -e 0.01 -o %s.4.htb" % (fo, fo)
    os.system(cmd)
    cmd = "cut -f2-4,6,7-9,11-13,14-15 %s.4.htb > %s.tsv" % (fo, fo)
    os.system(cmd)
from itertools import groupby
from operator import itemgetter
def filter_nbs(fi, fo):
    ary = read_tsv(fi, names = True)
    ary_sorted = sorted(ary, key = lambda x: (x[0], x[1]))

    fho = open(fo, "w")
    fho.write("\t".join(["id", "size", "doms", "nbs_e", "nbs_qb", "nbs_qe", "nbs_hb", "nbs_he", "nbs_qloc", "nbs_hloc"]) + "\n")
    for key, valuesiter in groupby(ary_sorted, key = itemgetter(0)):
        res = [v for v in valuesiter if v[8] < 1e-5 and v[2]-v[1]+1 >= 20]
        nbsdoms = [v for v in res if v[4] == 'NB-ARC']
        if len(nbsdoms) == 0:
            continue
        nbsdoms = sorted(nbsdoms, key = lambda x: x[8], reverse = True)
        nbsdom = nbsdoms[0]

        doms = [v[4] for v in res]
        (org, gid) = str(key).split("|", 1)
        size = res[0][2]
        aryo = (key, size, ",".join(doms))
        aryo += tuple([nbsdom[i] for i in [8,1,2,5,6,10,11]]) 
        fho.write("%s\t%s\t%s\t%g\t%d\t%d\t%d\t%d\t%s\t%s\n" % aryo)
    fho.close()
def extract_nbs(fi, fo, fd_cds, fd_pro, fo_pro, fn_cds, fn_pro):
    cds_dict = SeqIO.index(fd_cds, "fasta")
    pro_dict = SeqIO.index(fd_pro, "fasta")
    fhi = open(fi, "r")
    fho = open(fo, "w")
    header = fhi.readline().strip("\n")
    fho.write(header + "\n")
    (nbs_cds_rcds, nbs_pro_rcds, pro_rcds) = ([], [], [])
    for line in fhi:
        line = line.strip("\n")
        gid, size, doms, e, domstr, tag, beg, end = line.split("\t")
        if tag == '': continue
        beg, end = int(beg), int(end)
        cds = str(cds_dict[gid].seq[(beg-1)*3:end*3])
        pro = str(pro_dict[gid].seq[(beg-1):end])
        cds_rcd = SeqRecord(Seq(cds), id = gid, description = "")
        pro_rcd = SeqRecord(Seq(pro), id = gid, description = "")
        nbs_cds_rcds.append(cds_rcd)
        nbs_pro_rcds.append(pro_rcd)
        pro_rcds.append(pro_dict[gid])
        fho.write(line + "\n")
    fho.close()
    SeqIO.write(nbs_cds_rcds, fn_cds, "fasta")
    SeqIO.write(nbs_pro_rcds, fn_pro, "fasta")
    SeqIO.write(pro_rcds, fo_pro, "fasta")
from Bio import AlignIO
from Bio.Align import AlignInfo
def make_cluster_consensus(fi, fs, fo, diro):
    fhi = open(fi, "r")
    dc = {}
    for line in fhi:
        line = line.strip("\n")
        (grp, sid) = line.split("\t")
        if grp == 'grp': continue
        if grp in dc:
            dc[grp].append(sid)
        else:
            dc[grp] = [sid]
    fhi.close()
    
    fhs = open(fs, "rU")
    ds = SeqIO.to_dict(SeqIO.parse(fhs, "fasta"))
    fhs.close()

    if not op.exists(diro):
        os.makedirs(diro)
    else:
        os.system("rm -rf %s/*" % diro)
    cnt = 0
    con_rcds = [SeqRecord(
        Seq('dfedlVGieahlkkleslLcldsddeVrmiGIwGmgGIGKTTLAraLynqlssqFdlscFvenskefsveqrpigldelgmqeqlLskilnqkdieienhlgvlkerLkdkKVLlVLDDVdkleqLdaLaketpWfGpGSRIIiTTrDkslLkahginhiYeVkcLseeeAlqLFcryAFgqksppdgfeeLskvevvklcgGLPLALkVlGssLrgKgskeeWedalprLetsldgenIesvLkvSYDgLdeedKdcFLyiAcFFniekedlVkylLaeggeldgrvglkvLvdksLitisddgrveMHdLlremgreiaseegcrpgkrqfLvdapeicdvltnrtgtsVlGIsLDslsinkieliisekafkrMrnLrfLkiyssh'), 
        id = 'Arabidopsis', description='')]
    for grp in dc.keys():
        sids = dc[grp]
        if len(sids) == 1:
            sid = sids[0]
            new_rcd = SeqRecord(ds[sid].seq, id=grp, description=sid)
            con_rcds.append(new_rcd)
            continue
        rcds = [ds[sid] for sid in sids]
        rcds.sort(key=lambda x: len(x), reverse=True)
        pre = "%s/%s" % (diro, grp)
        fht = open("%s.1.fas" % pre, "w")
        SeqIO.write(rcds, fht, "fasta")
        fht.close()
        os.system("muscle -quiet -in %s.1.fas -out %s.2.fas" % (pre, pre))
        os.system("hmmbuild --amino %s.3.hmm %s.2.fas" % (pre, pre))
        os.system("hmmemit -C %s.3.hmm > %s.4.fas" % (pre, pre))

        #aln = AlignIO.read("%s/%s.fas" % (diro, grp), 'fasta')
        #aln_info = AlignInfo.SummaryInfo(aln)
        #con_seq = aln_info.dumb_consensus(threshold=0.5)
        fhs = open("%s.4.fas" % pre, "r")
        rcds = list(SeqIO.parse(fhs, "fasta"))
        fhs.close()
        con_rcd = SeqRecord(rcds[0].seq, id=grp, description="%d" % len(sids))
        con_rcds.append(con_rcd)
        cnt += 1
    #os.system("rm %s/*.in.fas" % diro)
    
    fho = open(fo, 'w')
    SeqIO.write(con_rcds, fho, "fasta")
    fho.close()
    print("%d non-singleton clusters written" % cnt)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Identify, cluster and characterizie plant NBS-LRR genes'
    )
 #   parser.add_argument('program', type=str, help='progname')
    parser.add_argument(
        'cfgfile', help = 'config file (a text file with species identifier followed by the absolute path of CDS fasta in each line)'
    )
    parser.add_argument(
        'outdir', help = 'output directory'
    )
    parser.add_argument(
        '--cpu', dest = "ncpu", default = nproc, type = int, help = 'number processors to use (default: all/%d)' % nproc)
    args = parser.parse_args()

    (fc, dirw) = (args.cfgfile, args.outdir)
    (orgs, fis, fgs) = read_cfg(fc)
    if not op.exists(dirw): os.makedirs(dirw)
    
    cdir = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] = os.environ['PATH']+':'+cdir
    cwd = os.getcwd()
    os.chdir(dirw)
    
    f_hmm = op.join(cdir, 'rdoms.hmm')
    if not os.access(f_hmm, os.R_OK):
        print("no access to %s" % f_hmm)
        sys.exit(1)
    f_pfm = '/home/youngn/zhoux379/data/db/pfam_v29/Pfam-A.hmm'
    if not os.access(f_pfm, os.R_OK):
        print("no access to %s" % f_pfm)
        sys.exit(1)
    
    print("%d input files detected" % len(fis))
    print("species to work on: %s" % "  ".join(orgs))
    print("output directory: %s" % dirw)
   
    merge_seqs(fis, orgs, "00.cds.raw.fas")
    os.system("seq.check.pl -i 00.cds.raw.fas -o 01.cds.fas")
    os.system("dna2pro.pl -i 01.cds.fas -o 05.pro.fas")
    pfam_scan(f_hmm, '05.pro.fas', '11', args.ncpu)
    filter_nbs("11.tsv", "13.tsv")
    os.system("htb.nbs.pl -i 13.tsv -o 14.tsv") 
    os.system("ncoils.py 05.pro.fas")
    extract_nbs("14.tsv", "21.tsv", "01.cds.fas", "05.pro.fas", "22.pro.fas", "23.cds.nbs.fas", "23.pro.nbs.fas")
    os.system("seqlen.py 22.pro.fas 22.pro.tbl")
    pfam_scan(f_pfm, '22.pro.fas', '27.pfam', args.ncpu)
   
    idty = 0.8
    cmd = "usearch -cluster_fast %s -sort length -id 0.8 -uc 31.uc" % '23.cds.nbs.fas'
    cmd = "usearch -cluster_agg %s -linkage min -clusterout 31.cluster -id %.02f" % ('23.pro.nbs.fas', idty)
    print(cmd)
    os.system(cmd)
    ###usearch.uc2tbl('31.uc', '32.tbl')
    usearch.cluster2tbl('31.cluster', '32.tbl')

    make_cluster_consensus('32.tbl', '23.pro.nbs.fas', '42.con.fas', '41_alns')
    os.system("muscle -in 42.con.fas -out 43.aln -clw")
    os.system("aln.conv.py --fmt 1 43.aln 43.fas")
    os.system("aln.conv.py --fmt 2 43.aln 43.phy")
    ##os.system("FastTree -wag -boot 1000 43.fas > 44.ft.nwk")
    os.system('phy.py 42.con.fas 43 --seqtype aa')
    sys.exit(1);
