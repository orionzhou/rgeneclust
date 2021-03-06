#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import math 
import re 
import os.path as op
import numpy as np

def uc2tbl(fi, fo):
    fhi = open(fi, "r")
    fho = open(fo, "w")
    fho.write("grp\tid\n")

    for line in fhi:
        line = line.strip("\n")
        (tag, grp, len_size, idty, srd, ign1, ign2, cigar, id, tgt) = line.split("\t")
        grp = int(grp) + 1
        if tag == "H" or tag == "S":
            fho.write("%d\t%s\n" % (grp, id))
    fhi.close()
    fho.close()
def cluster2tbl(fi, fo):
    #cmd = "echo -e \"grp\tid\" | cat - %s > /tmp/out && mv /tmp/out %s" % (fi, fo)
    #os.system(cmd)
    fhi = open(fi, "r")
    fho = open(fo, "w")
    fho.write("grp\tid\n")
    for line in fhi:
        line = line.strip("\n")
        (grp, gid) = line.split("\t")
        grp = int(grp) + 1
        fho.write("%d\t%s\n" % (grp, gid))
    fhi.close()
    fho.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'cluster a set of protein sequences'
    )
    parser.add_argument(
        'fi', help = 'input file (fasta)'
    )
    parser.add_argument(
        'dirw', help = 'output directory'
    )
    parser.add_argument(
        'idty', nargs = '?', help = 'clustering similarity (0.7)'
    )
    args = parser.parse_args()
   
    idty = 0.7
    (fi, dirw) = (args.fi, args.dirw)
    if args.idty: idty = float(args.idty)
    if not op.exists(dirw): os.makedirs(dirw)
    
    f31 = op.join(dirw, "31.uc")
    f32 = op.join(dirw, "32.tbl")
    cmd = "usearch -cluster_fast %s -sort length -id %g -uc %s" % (fi, idty, f31)
    os.system(cmd)
    usearch2tbl(f31, f32)
