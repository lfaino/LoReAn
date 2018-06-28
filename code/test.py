#! /usr/bin/env python3


import sys

from Bio import SeqIO


# OTHER SCRIPTS


def main(name_file):
    sizes = [rec.id for rec in SeqIO.parse(name_file, "fasta")]
    print (len(sizes))
if __name__ == '__main__':
    main(*sys.argv[1:])