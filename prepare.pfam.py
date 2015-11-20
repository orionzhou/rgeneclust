#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import os.path as op
import argparse

import urllib2
from urlparse import urlparse
from ftplib import FTP

def internet_on():
    try:
        resp = urllib2.urlopen("http://www.google.com/", timeout = 3)
        return True
    except urllib2.URLError as err: pass
    return False
def download_from_url(url, fo):
    if not internet_on():
        print "no internet connection"
        sys.exit(1)
    o = urlparse(url)
    fname = op.basename(o.path)
    if(o.scheme == 'ftp'):
        print "downloading " + url
        ftp = FTP(o.netloc)
        ftp.login()
        ftp.cwd(op.dirname(o.path))
        ftp.retrbinary('RETR '+fname, open(fo, 'wb').write)
    else:
        print "%s: scheme not supported"
    return fo
def uncompress(f_gzipped):
    fo, ext = op.splitext(f_gzipped)
    if(ext == '.gz'):
        print "uncompressing " + f_gzipped
        os.system("gunzip " + f_gzipped)
    else:
        print "not a valid gzip file: " + f_gzipped
        sys.exit(1)
    return fo
def hmmpress(f_hmm):
    fname, ext = op.splitext(f_hmm)
    if ext == '.hmm':
        print "prepare HMM database for scan..."
        os.system("hmmpress " + f_hmm)
    else:
        print "not a valide HMM file: " + f_hmm
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'download and prepare Pfam database for hmmscan'
    )
    parser.add_argument(
        'out', default = "test", help = 'output directory (default: "test")'
    )
    args = parser.parse_args()

    dirw = args.out
    if not op.exists(dirw): os.makedirs(dirw)

    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
    f_pfam_gz = op.join(dirw, op.basename(url))
    download_from_url(url, f_pfam_gz)
    f_pfam = op.splitext(f_pfam_gz)[0]
    uncompress(f_pfam_gz)
    hmmpress(f_pfam)
