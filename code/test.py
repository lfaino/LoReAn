#! /usr/bin/env python3


import gzip
import subprocess
import sys

import mapping
from Bio import SeqIO
from Bio.SeqUtils import GC


def main(reference_database, reads, threads, max_intron_length, working_dir, verbose):


    subset_fasta = reads + "subset.100000.fasta"
    with open(subset_fasta, "w") as fh:
        for rec in SeqIO.parse(reads, "fasta"):
            SeqIO.write(rec, fh, "fasta")

    bam = mapping.minimap(reference_database, subset_fasta, threads, max_intron_length, working_dir, verbose)
    fasta_gz = bam + ".fasta.gz"
    cmd = "extractSoftclipped %s > %s" % (bam, fasta_gz)
    extract_clip = subprocess.Popen(cmd, shell=True)
    extract_clip.communicate()

    list_short = []
    list_long = []
    dict_uniq = {}
    with gzip.open(fasta_gz, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            name_seq = str(rec.id)
            name = name_seq.split("_")[0]
            if name in dict_uniq:
                if len(dict_uniq[name].seq) > len(rec.seq):
                    list_long.append(dict_uniq[name])
                    list_short.append(rec)
                else:
                    list_short.append(dict_uniq[name])
                    list_long.append(rec)
            else:
                dict_uniq[name] = rec
    long_file = "long.fasta"
    short_file = "short.fasta"
    with open(long_file, "w") as fh:
        SeqIO.write(list_long, fh, "fasta")
    with open(short_file, "w") as fh:
        SeqIO.write(list_short, fh, "fasta")
    kmer_start = 21
    poly_AAA = True
    while poly_AAA and kmer_start < 70:
        cmd = "jellyfish count -s 100000 -m %s -o %s.kmer %s" % (kmer_start, kmer_start, long_file)
        jelly_count = subprocess.Popen(cmd, cwd = "./" , shell=True)
        jelly_count.communicate()

        cmd = "jellyfish dump -L 2 -ct %s.kmer | sort -k2n | tail -n 1" % kmer_start
        jelly_dump = subprocess.Popen(cmd, cwd = "./", stdout=subprocess.PIPE, shell=True)
        out_dump= jelly_dump.communicate()[0].decode('utf-8')
        mer = out_dump.split("\t")
        gc_mer = GC(mer[0])
        kmer_start += 5
        if gc_mer > 40:
            poly_AAA = False
        else:
            mer_poly = mer
    print (mer_poly)



if __name__ == '__main__':
    main(*sys.argv[1:])