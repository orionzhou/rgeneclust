#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import argparse

import multiprocessing
nproc = multiprocessing.cpu_count()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'identify sequences with coiled-coil domains'
    )
    parser.add_argument(
        'seqfile', help = 'protein sequence file'
    )
    parser.add_argument(
        '--out', dest = "output", default = "ncoil", help = 'output directory (default: "ncoil")'
    )
    args = parser.parse_args()

    (fi, dirw) = (args.seqfile, args.output)
    if not op.exists(dirw): os.makedirs(dirw)
    if not os.access(fi, os.R_OK):
        print "no access to input file: %s" % fi
        sys.exit(1)
   
    cdir = os.path.dirname(os.path.realpath(__file__))
    os.environ['COILSDIR'] = op.join(cdir, 'coils')
    f_ncoil = op.join(cdir, "coils/ncoils-osf")
    fi = op.realpath(fi)

    cwd = os.getcwd()
    os.chdir(dirw)
    os.system("rm *.out")
    while not op.exists('01.out') or op.getsize('01.out') == 0:
        os.system("%s -c -win 21 -min_P 0.75 <%s >%s" % (f_ncoil, fi, '01.out'))
    while not op.exists('02.out') or op.getsize('02.out') == 0:
        os.system("%s -c -w -win 21 -min_P 0.75 <%s >%s" % (f_ncoil, fi, '02.out'))
    while not op.exists('03.out') or op.getsize('03.out') == 0:
        os.system("%s -c -win 28 -min_P 0.75 <%s >%s" % (f_ncoil, fi, '03.out'))
    while not op.exists('04.out') or op.getsize('04.out') == 0:
        os.system("%s -c -w -win 28 -min_P 0.75 <%s >%s" % (f_ncoil, fi, '04.out'))
    os.system("perl -ne '/:\s+(\S+)/; print $1, \"\\n\";' %s |sort >%s" % ('01.out', '11.txt'))
    os.system("perl -ne '/:\s+(\S+)/; print $1, \"\\n\";' %s |sort >%s" % ('02.out', '12.txt'))
    os.system("perl -ne '/:\s+(\S+)/; print $1, \"\\n\";' %s |sort >%s" % ('03.out', '13.txt'))
    os.system("perl -ne '/:\s+(\S+)/; print $1, \"\\n\";' %s |sort >%s" % ('04.out', '14.txt'))
    
    os.system('comm -1 -2 11.txt 12.txt > 21.txt')
    os.system('comm -1 -2 13.txt 14.txt > 22.txt')
    os.system("comm 21.txt 22.txt | perl -pe 's/^\\t+//' > 31.txt")
    os.chdir(cwd)
