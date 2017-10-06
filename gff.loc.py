#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Extract gene locations from GFF file'
    )
    parser.add_argument(
        'fi', help = 'input file (GFF)'
    )
    parser.add_argument(
        'fo', help = 'output file (tsv)'
    )
    parser.add_argument(
        '--opt', dest='opt', default='mRNA', help='Type (3rd column): mRNA'
    )
    parser.add_argument(
        '--tag', dest='tag', default='Name', help='tag: Name'
    )
    args = parser.parse_args()
    
    (fi, fo, opt, tag) = (args.fi, args.fo, args.opt, args.tag)
    fho = open(fo, "w")
    fho.write("id\tchr\tbeg\tend\n")
    
    fhi = open(fi, "r")
    for line in fhi:
        line = line.strip("\n")
        line = line.strip("\r")
        if line == "":
            break
        ps = line.split("\t")
        if len(ps) < 9:
            continue
        (chrom, beg, end, opt2, desc) = [ps[i] for i in [0,3,4,2,8]]
        if opt2 == opt:
            gid = ''
            pairs = desc.split(";")
            for pair in pairs:
                x = pair.split("=", 1)
                if len(x) < 2:
                    continue
                (tag2, value) = x
                if tag2 == tag:
                    gid = value
                    break
            if gid == '':
                print("no gid for:\n%s" % line)
            fho.write("%s\t%s\t%s\t%s\n" % (gid, chrom, beg, end))
    fho.close()


