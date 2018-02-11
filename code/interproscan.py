#!/usr/bin/env python3


import datetime
import os
import subprocess
import sys
import tempfile

import dirsAndFiles as logistic
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq

#======================================================================================================================

GFFREAD = 'gffread -g %s -y %s %s'

IPRSCAN = 'interproscan.sh -i %s -cpu %s'

#======================================================================================================================




def iprscan(ref, gff_file, root, threads):


    fmtdate = '%H:%M:%S %d-%m'
    now = datetime.datetime.now().strftime(fmtdate)
    wd = os.path.join(root,"tmp")
    logistic.check_create_dir(wd)
    fasta_file_outfile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".log")
    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".err")
    prot_file_out = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".fasta")
    prot_file_mod = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.mod.", suffix=".fasta")

    com = GFFREAD % (os.path.abspath(ref), prot_file_out.name, os.path.abspath(gff_file))
    call = subprocess.Popen(com, stdout=fasta_file_outfile, cwd = wd, stderr=errorFilefile, shell=True)
    call.communicate()
    count = 0
    input_file = open(prot_file_out.name)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    for id in fasta_dict:
        count += 1
        if "." in fasta_dict[id].seq:
            mutable_seq = MutableSeq(str(fasta_dict[id].seq), IUPAC.unambiguous_dna)
            mutable_seq.remove(".")
            fasta_dict[id].seq = mutable_seq
            SeqIO.write(fasta_dict[id], prot_file_mod, "fasta")
        else:
            SeqIO.write(fasta_dict[id], prot_file_mod, "fasta")
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

    return (prot_file_mod.name + ".tsv")



if __name__ == '__main__':
    iprscan(*sys.argv[1:])


