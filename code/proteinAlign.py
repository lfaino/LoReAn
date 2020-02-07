import os
import re
import subprocess
import sys
import tempfile
import warnings
from multiprocessing import Pool

import tqdm
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#bedtools getfasta -fo /tmp/pybedtools.423huf8a.tmp -fi /data/lfainoData/lorean/LoReAn_Example/JR2/LoReAn_annotation/run/chr8.fasta.masked.fasta.rename.fasta -bed /tmp/pybedtools.aeqgr5mv.tmp

EXONERATE = 'exonerate --model protein2genome --bestn 1  --showvulgar no --showalignment no --showquerygff no --showtargetgff yes --query %s --target %s' #--refine region
BLASTP = 'diamond blastp -q %s --db %s -k 1 -p %s --out %s --evalue 1e-15'
CONVERT = 'exonerate_gff_to_alignment_gff3.pl /dev/stdin '
BEDTOOLS = 'bedtools getfasta -fo %s -fi %s -bed %s'
MAKEDB = 'diamond makedb --in %s -d %s -p %s'

warnings.filterwarnings("ignore")
gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'Z', 'TAG':'Z',
      'TGC':'C', 'TGT':'C', 'TGA':'Z', 'TGG':'W'}

basepairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def translate_frameshifted( sequence ):
      translate = ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
      return translate

def reverse_complement( sequence ):
      reversed_sequence = (sequence[::-1])
      rc = ''.join([basepairs.get(reversed_sequence[i], 'X') for i in range(len(sequence))])
      return rc

def transeq(data):
    dummy = int(data[1])
    record = data[0]
    if dummy == 0:
        prot = (translate_frameshifted(record.seq[0:]))
        prot_rec = (SeqRecord(Seq(prot, IUPAC.protein), id=record.id + "_strand0plus"))
    if dummy == 1:
        prot = (translate_frameshifted(record.seq[1:]))  # second frame
        prot_rec = (SeqRecord(Seq(prot, IUPAC.protein), id=record.id + "_strand1plus"))
    if dummy == 2:
        prot = (translate_frameshifted(record.seq[2:]))  # third frame
        prot_rec =(SeqRecord(Seq(prot, IUPAC.protein), id=record.id + "_strand2plus"))
    if dummy == 3:
        prot = (translate_frameshifted(reverse_complement(record.seq)))  # negative first frame
        prot_rec = (SeqRecord(Seq(prot, IUPAC.protein), id=record.id + "_strand0minus"))
    if dummy == 4:
        prot = (translate_frameshifted(reverse_complement(record.seq[:len(record.seq) - 1])))  # negative second frame
        prot_rec =(SeqRecord(Seq(prot, IUPAC.protein), id=record.id + "_strand1minus"))
    if dummy == 5:
        prot = (translate_frameshifted(reverse_complement(record.seq[:len(record.seq) - 2])))  # negative third frame
        prot_rec = (SeqRecord(Seq(prot, IUPAC.protein), id=record.id + "_strand2minus"))
    return(prot_rec)


def protAlign(genome, fasta, nproc, wd, verbose):
    output_dimaonds = os.path.join(wd, "output_diamonds.txt")
    output_dimaonds_done = os.path.join(wd, "output_diamonds.done.txt")
    genome_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
    if not os.path.exists(output_dimaonds) and not os.path.exists(output_dimaonds_done):
        translate_genome = genome + ".prot"#os.path.join(wd, "translatedProtGenome.fasta")
        results_get = []
        list_fasta = []
        for record in genome_dict:
            count = 0
            for strand in range(6):
                count += 1
                list_fasta.append([genome_dict[record], str(strand)])
        pool = Pool(processes=int(nproc))
        for x in tqdm.tqdm(pool.imap_unordered(transeq, list_fasta), total=len(list_fasta)):
            results_get.append(x)
            pass
        sys.stdout.write("\n###RUNNING DIAMOND MAKEDB###\n")
        with open(translate_genome, "w") as output_handle:
            SeqIO.write(results_get, output_handle, "fasta")
        com = MAKEDB % (translate_genome, translate_genome, str(nproc))
        if verbose:
            sys.stdout.write(com)
        call = subprocess.Popen(com, shell=True , cwd = wd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        err, out = call.communicate()
        if verbose:
            sys.stdout.write(err.decode())
            sys.stdout.write(out.decode())
        sys.stdout.write("\n###RUNNING DIAMOND ###\n")
        com = BLASTP % (fasta, translate_genome, str(nproc), output_dimaonds)
        if verbose:
            sys.stdout.write(com)
        call = subprocess.Popen(com, shell=True, stdout= subprocess.PIPE, stderr=subprocess.PIPE)
        err, out = call.communicate()
        if verbose:
            sys.stdout.write(err.decode())
            sys.stdout.write(out.decode())
        os.mknod(output_dimaonds_done)

    list_match = []
    with open(output_dimaonds, "r") as fh:
        for line in fh:
            list_match.append(line)

    record_dict = {}
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id not in record_dict:
            record_dict[record.id] = record
    list_fasta = []
    for align in list_match:
        if align != "":
            name_prot = align.split("\t")
            if len(name_prot) == 12:
                chr = name_prot[1].split("_")[0]
                #chr.pop()
                #chr_name = "_".join(chr)
                if verbose:
                    list_fasta.append([align, genome, record_dict[name_prot[0]], len(genome_dict[chr].seq), wd, "True"])
                else:
                    list_fasta.append([align, genome, record_dict[name_prot[0]], len(genome_dict[chr].seq), wd, "False"])
    pool = Pool(processes=int(nproc))
    results_get = []
    sys.stdout.write("###RUNNING EXONERATE ###\n")
    for x in tqdm.tqdm(pool.imap_unordered(runExonerate, list_fasta), total=len(list_fasta)):
        results_get.append(x)
        pass
    match_id = 0
    final_ouput = os.path.join(wd, "protein_evidence.gff3")
    #print(final_ouput)
    with open(final_ouput, "w") as fh:
        for match in results_get:
            match_id += 1
            id = "ID=match" + str(match_id)
            if match is not None:
                coords = match.split("\n")
                for line in coords:
                    if not line.startswith("#") and "exonerate:protein2genome:local" in line:
                        elem = line.split("\t")
                        if len(elem) == 9 and len(re.split(':|-', elem[0])) == 3 and elem[2] == "gene":
                            target = elem[8].split(" ; ")[1].split(" ")[1]
                        elif len(elem) == 9 and len(re.split(':|-', elem[0])) == 3 and elem[2] == "exon":
                            loc = elem[0].split(":")[0]
                            elem[3] = str(int(elem[3]) + int(elem[0].split(":")[1].split("-")[0]))
                            elem[4] = str(int(elem[4]) + int(elem[0].split(":")[1].split("-")[0]))
                            elem[0] = loc
                            elem[8] = id + ";Target=" + target
                            elem[2] = "nucleotide_to_protein_match"
                            elem[1] = "exonerate"
                            fh.write("\t".join(elem) + "\n")
    return(final_ouput)


def runExonerate(sequence):
    align, genome, prot, length, wd, verbose = sequence
    elem = align.split("\t")
    name_prot = os.path.join(wd,  elem[0] + ".fasta")
    #name_gff = os.path.join(wd,  elem[0] + ".gff")
    with open(name_prot, "w") as output_handle:
        SeqIO.write(prot , output_handle, "fasta")
    if len(elem) == 12:
        if float(elem[10]) < 1e-100:
            if elem[1].endswith("plus"):
                begin = (int(elem[8]) * 3) - 100000
                stop = (int(elem[9]) * 3)  + 100000
                if begin < 0 :
                    begin = "0"
                else:
                    begin = str(begin)
                if stop > int(elem[9] * 3):
                    stop = elem[3] * 3 - 3
                else:
                    stop = str(stop)
            else:
                stop = int(length) - (int(elem[8]) * 3)
                begin = int(length) -(int(elem[9]) * 3)
                begin = begin - 100000
                stop = stop  + 100000
                if begin < 0 :
                    begin = "0"
                else:
                    begin = str(begin)
                if stop > int(elem[9] * 3):
                    stop = elem[3] * 3 - 3
                else:
                    stop = str(stop)
            #chr = elem[1].split("_")[0]
            chr = elem[1].split("_")[0]
            #chr.pop()
            #chr_name = "_".join(chr)
            new_coords = "\t".join([chr, begin, stop]) + "\n"
            if verbose:
                outfile_bed = tempfile.NamedTemporaryFile(delete=False, suffix=".bed")
            else:
                outfile_bed = tempfile.NamedTemporaryFile(suffix=".bed")
            with open(outfile_bed.name, "w") as fp:
                fp.write(new_coords)
            if verbose:
                outfile_fo_fasta = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
            else:
                outfile_fo_fasta = tempfile.NamedTemporaryFile(suffix=".fasta")
            com_bedtools = BEDTOOLS % (outfile_fo_fasta.name, genome, outfile_bed.name)
            call_bedtools= subprocess.Popen(com_bedtools, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            call_bedtools.communicate()

            com_exo = EXONERATE % (name_prot, outfile_fo_fasta.name)
            if verbose:
                sys.stdout.write(com_exo)
            #com_conv = CONVERT
            call_exo = subprocess.Popen(com_exo, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #call_conv = subprocess.Popen(com_conv, shell=True, stdout=subprocess.PIPE, stdin=call_exo.stdout, stderr=subprocess.PIPE)
            output_final, err = call_exo.communicate()
            if output_final.decode != "":
                return (output_final.decode())


if __name__ == '__main__':
    protAlign(*sys.argv[1:])


#Contig1 nap-nr_minus_rice.fasta nucleotide_to_protein_match     8208    8276    50.00   +       .       ID=match.nap.nr_minus_rice.fasta.120;Target=RF|XP_623193.1|66524404|XM_623190 1 23
