#!/usr/bin/env python3


import datetime
import os
import subprocess
import sys
import tempfile
import warnings

from Bio import SeqIO
from Bio.Seq import Seq

#======================================================================================================================

GFFREAD = 'gffread -C -g %s -y %s %s'

IPRSCAN = 'interproscan.sh -i %s -cpu %s'

IPRSCAN_HELP = 'interproscan.sh -h'

#======================================================================================================================


def check_iprscan():
    com = IPRSCAN_HELP
    call = subprocess.Popen(com, shell=True)
    err, log = call.communicate()

def iprscan(ref, gff_file, wd, threads):

    warnings.filterwarnings("ignore")
    fmtdate = '%H:%M:%S %d-%m'
    now = datetime.datetime.now().strftime(fmtdate)
    fasta_file_outfile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".log")
    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".err")
    prot_file_out = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".fasta")
    prot_file_mod = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.mod.", suffix=".fasta")

    com = GFFREAD % (os.path.abspath(ref), prot_file_out.name, os.path.abspath(gff_file))
    call = subprocess.Popen(com, stdout=fasta_file_outfile, cwd = wd, stderr=errorFilefile, shell=True)
    call.communicate()
    input_file = open(prot_file_out.name)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    count = len(fasta_dict)

    bad_prot = []
    for id in fasta_dict:
        if "." in str(fasta_dict[id].seq):
            count += 1
            bad_prot.append(fasta_dict[id].description)
            prot = str(fasta_dict[id].seq)
            prot_mod = prot.replace(".","")
            fasta_dict[id].seq = Seq(prot_mod)
        SeqIO.write(fasta_dict[id], prot_file_mod, "fasta")

    bad_gene = wd + "/bad_gene.txt"
    with open(bad_gene, "w") as fh:
        for line in bad_prot:
            fh.write(line + "\n")

    sys.stdout.write(("\n###INTERPROSCAN ANALYSIS STARTED AT:\t" + now + "\t###\n###RUNNING ANALYSIS FOR \t\033[32m" + str(count) + "\033[0m\t mRNA\t###\n"))

    cmd = IPRSCAN %(prot_file_mod.name, threads)
    err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix=prot_file_mod.name, suffix=".err")
    log = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix=prot_file_mod.name, suffix=".log")
    iprscan = subprocess.Popen(cmd, cwd=wd, stderr = err, stdout = log, shell=True)
    iprscan.communicate()

    done_prot = {}
    tsv_file = prot_file_mod.name + ".tsv"
    with open(tsv_file, "r") as fh:
        for line in fh:
            mRNA = line.split("\t")[0]
            done_prot[mRNA] = mRNA
    sys.stdout.write(("\n###FINISHED TO RUN INTERPROSCAN ANALYSIS AT:\t" + now + "\t###\n###PROTEINS DOMAINS WERE FOUND FOR \t\033[32m" + str(len(done_prot)) + "\033[0m\t PROTEINS\t###\n"))
    final_annot = gff_file + ".tsv"
    os.rename(prot_file_mod.name + ".tsv", final_annot)
    return final_annot, bad_gene



if __name__ == '__main__':
    iprscan(*sys.argv[1:])


