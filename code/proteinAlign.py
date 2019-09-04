import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

import tqdm
from Bio import SeqIO

#bedtools getfasta -fo /tmp/pybedtools.423huf8a.tmp -fi /data/lfainoData/lorean/LoReAn_Example/JR2/LoReAn_annotation/run/chr8.fasta.masked.fasta.rename.fasta -bed /tmp/pybedtools.aeqgr5mv.tmp

EXONERATE = 'exonerate --model protein2genome --bestn 1 --refine region --showtargetgff TRUE --showquerygff TRUE --showalignment FALSE --showvulgar FALSE  --query %s --target %s'
TBLASTN = 'tblastn -query %s -subject %s -outfmt "6 sseqid sstart send slen evalue" -max_target_seqs 1'
CONVERT = 'exonerate_gff_to_alignment_gff3.pl /dev/stdin'
BEDTOOLS = 'bedtools getfasta -fo %s -fi %s -bed %s'

def protAlign(genome, fasta, nproc, wd, verbose):
    list_fasta = []
    for record in SeqIO.parse(fasta, "fasta"):
        list_fasta.append([record, genome])
    pool = Pool(processes=int(nproc))
    results_get = []
    sys.stdout.write("###RUNNING EXONERATE ###\n")
    for x in tqdm.tqdm(pool.imap_unordered(runExonerate, list_fasta), total=len(list_fasta)):
        results_get.append(x)
        pass
    with open(os.path.join(wd, "protein_evidence.gff3"), "w") as fh:
        for gene in results_get:
            if gene is not None:
                coords = gene.split("\n")
                for line in coords:
                    if line.strip() != "":
                        elem = line.split("\t")
                        location = elem[0].split(":")
                        start = location[1].split("-")[0]
                        elem[0] = location[0]
                        elem[3] = str(int(start) + int(elem[3]))
                        elem[4] = str(int(start) + int(elem[4]))
                        fh.write("\t".join(elem) + "\n")

    return


def runExonerate(sequence):
    outfile_fasta = tempfile.NamedTemporaryFile()
    with open(outfile_fasta.name, "w") as fp:
        SeqIO.write(sequence[0], fp, "fasta")
    com = TBLASTN % (outfile_fasta.name, sequence[1])
    call = subprocess.Popen(com, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out_b, err_b = call.communicate()
    list_match = out_b.decode().split("\n")
    count = 0
    if len(list_match) > 0:
        for coord in list_match:
            count = count + 1
            elem = coord.split("\t")
            if coord != '' and float(elem[4]) < 1e-5:
                begin = int(elem[1]) - 100000
                stop = int(elem[2]) + 100000
                if begin < 0 :
                    elem[1] = "0"
                else:
                    elem[1] = str(begin)
                if stop > int(elem[3]):
                    elem[2] = elem[3]
                else:
                    elem[2] = str(stop)
                del elem[-1]
                del elem[-1]
                new_coords = "\t".join(elem) + "\n"
                outfile_bed = tempfile.NamedTemporaryFile()
                with open(outfile_bed.name, "w") as fp:
                    fp.write(new_coords)
                outfile_fo_fasta = tempfile.NamedTemporaryFile()
                com_bedtools = BEDTOOLS % (outfile_fo_fasta.name, sequence[1], outfile_bed.name)
                call_bedtools= subprocess.Popen(com_bedtools, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                call_bedtools.communicate()
                com_exo = EXONERATE % (outfile_fasta.name, outfile_fo_fasta.name)
                com_conv = CONVERT
                call_exo = subprocess.Popen(com_exo, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                call_conv = subprocess.Popen(com_conv, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=call_exo.stdout)
                out, err = call_conv.communicate()
                output = out.decode()
                return (output)

if __name__ == '__main__':
    protAlign(*sys.argv[1:])


#Contig1 nap-nr_minus_rice.fasta nucleotide_to_protein_match     8208    8276    50.00   +       .       ID=match.nap.nr_minus_rice.fasta.120;Target=RF|XP_623193.1|66524404|XM_623190 1 23