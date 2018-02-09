#!/usr/bin/env python3


import os
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq

#======================================================================================================================

GFFREAD = 'gffread -g %s -y %s %s'

IPRSCAN = 'interproscan.sh -i %s -cpu %s'

#======================================================================================================================




def iprscan(ref, gff_file, root, threads):

    wd = os.path.join(root,"tmp")
    os.mkdir(wd)
    fasta_file_outfile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".log")
    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".err")
    prot_file_out = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.", suffix=".fasta")
    prot_file_mod = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="prot_gffread.mod.", suffix=".fasta")

    com = GFFREAD % (os.path.abspath(ref), prot_file_out.name, os.path.abspath(gff_file))
    call = subprocess.Popen(com, stdout=fasta_file_outfile, cwd = wd, stderr=errorFilefile, shell=True)
    call.communicate()

    input_file = open(prot_file_out.name)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    for id in fasta_dict:
        if "." in fasta_dict[id].seq:
            mutable_seq = MutableSeq(str(fasta_dict[id].seq), IUPAC.unambiguous_dna)
            mutable_seq.remove(".")
            fasta_dict[id].seq = mutable_seq
            SeqIO.write(fasta_dict[id], prot_file_mod, "fasta")
        else:
            SeqIO.write(fasta_dict[id], prot_file_mod, "fasta")

    cmd = IPRSCAN %(prot_file_mod.name, threads)
    err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix=prot_file_mod.name, suffix=".err")
    log = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix=prot_file_mod.name, suffix=".log")
    iprscan = subprocess.Popen(cmd, cwd=wd, stderr = err, stdout = log, shell=True)
    iprscan.communicate()





if __name__ == '__main__':
    iprscan(*sys.argv[1:])


