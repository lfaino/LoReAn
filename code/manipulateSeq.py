#!/usr/bin/env python3
import ctypes as ct
import datetime
import os
import subprocess
import sys
import tempfile
from glob import glob
from pathlib import Path

import align as align
from Bio import SeqIO
from Bio.Seq import reverse_complement

#==========================================================================================================
# COMMANDS LIST

BEDTOOLS_SORT = 'bedtools sort -i %s'

BEDTOOLS_MERGE = 'bedtools merge -d 1'

BEDTOOLS_MASK = 'bedtools maskfasta -fi %s -bed %s -fo %s'

AWK = 'awk \'{if ($3 - $2 > %s) print $0} \''

BUILD_TABLE  = 'build_lmer_table -sequence  %s -freq %s'

REPEAT_SCOUT =  'RepeatScout -sequence %s -output %s -freq %s'

REPEAT_MASKER = 'RepeatMasker %s -e ncbi -lib %s -gff -pa %s -dir %s'


#==========================================================================================================

def to_int(seq, lEle, dEle2Int):
    """
    it is used for the alignment of the primers to the reads
    :param seq: 
    :param lEle: 
    :param dEle2Int: 
    :return: 
    """
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i,ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n
    return num


def align_one(ssw, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    """
    this function calculate the score and other results from the alignment
    """
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, int(nMaskLen))
    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)
    return nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar


def align_call(elem):
    """
    this function call the aligner software
    :param elem: 
    :return: 
    """
    record = elem[0]
    adapter = elem[1]
    lEle = []
    dRc = {} 
    dEle2Int = {}
    dInt2Ele = {}
    nMatch = 2
    nMismatch = 1
    nOpen = 1
    nExt = -2
    nFlag = 0
    #if not args.sMatrix:
    lEle = ['A', 'C', 'G', 'T', 'N']
    for i,ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    nEleNum = len(lEle)
    lScore = [0 for i in range(nEleNum**2)]
    for i in range(nEleNum-1):
        for j in range(nEleNum-1):
            if lEle[i] == lEle[j]:
                lScore[i*nEleNum+j] = nMatch
            else:
                lScore[i*nEleNum+j] = -nMismatch
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore
    
    ssw = ssw_lib.CSsw("./")
    sQSeq = record.seq
    sQId = record.id
    if len(sQSeq) > 30:
        nMaskLen = len(sQSeq) / 2
    else:
        nMaskLen = 15
    outputAlign = []
    qNum = to_int(sQSeq, lEle, dEle2Int)
    qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)
    sQRcSeq = reverse_complement(sQSeq)
    qRcNum = to_int(sQRcSeq, lEle, dEle2Int)
    qRcProfile = ssw.ssw_init(qRcNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)
    sRSeq = adapter.seq
    sRId = adapter.id
    rNum = to_int(sRSeq, lEle, dEle2Int)
    res = align_one(ssw, qProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)
    resRc = None
    resRc = align_one(ssw, qRcProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)
    strand = 0
    if res[0] == resRc[0]:
        next
    if res[0] > resRc[0]:
        res = res
        strand = 0
        outputAlign = [sRId , sQId, strand, res]
    elif res[0] < resRc[0]:
        res = resRc
        strand = 1
        outputAlign = [sRId , sQId, strand, res]
    ssw.init_destroy(qProfile)
    ssw.init_destroy(qRcProfile)
    return outputAlign


def filterLongReads(fastq_filename, min_length, max_length, wd, adapter, threads, a):
    """
    Filters out reads longer than length provided and it is used to call the alignemnt and parse the outputs
    """
    scoring = [3, -6, -5, -2]
    align_score_value = ""

    if a and not adapter:
        out_filename = wd + fastq_filename + '.long_reads.filtered.fasta'
    else:
        out_filename = fastq_filename + '.long_reads.filtered.fasta'

    filter_count = 0

    if not os.path.isfile(out_filename):
        with open(out_filename, "w") as output_handle:
            if fastq_filename.endswith('fastq') or fastq_filename.endswith('fq'):
                for record in SeqIO.parse(fastq_filename, "fastq"):
                    if len(str(record.seq)) > int(min_length) < int(max_length):
                        record.description= ""
                        record.name = ""
                        record.id = str(filter_count)
                        filter_count += 1
                        SeqIO.write(record, output_handle, "fasta")
            elif fastq_filename.endswith('fasta') or fastq_filename.endswith('fa'):
                for record in SeqIO.parse(fastq_filename, "fasta"):
                    if int(min_length) < len(str(record.seq)) < int(max_length):
                        record.description= ""
                        record.name = ""
                        record.id = str(filter_count)
                        filter_count += 1
                        SeqIO.write(record, output_handle, "fasta")
    else:
        sys.stdout.write(('Filtered FASTQ existed already: ' + out_filename + ' --- skipping\n'))

    if adapter:
        out_filename_oriented = wd + fastq_filename + '.longreads.filtered.oriented.fasta'
        filter_count = align.adapter_alignment(out_filename, adapter, scoring, align_score_value, out_filename_oriented, threads)
        fmtdate = '%H:%M:%S %d-%m'
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write("###FINISHED FILTERING AT:\t" + now + "###\n\n###LOREAN KEPT\t\033[32m" + str(filter_count) +
                         "\033[0m\tREADS AFTER LENGTH FILTERING AND ORIENTATION###\n")

        return out_filename_oriented, filter_count
    else:
        fmtdate = '%H:%M:%S %d-%m'
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write("###FINISHED FILTERING AT:\t" + now + "###\n\n###LOREAN KEPT\t\033[32m" + str(filter_count) +
                         "\033[0m\tREADS AFTER LENGTH FILTERING###\n")

        return out_filename



def maskedgenome(wd, ref, gff3, length, verbose):
    """
    this module is used to mask the genome when a gff or bed file is provided
    """

    outputmerge = tempfile.NamedTemporaryFile(delete=False, mode="w", prefix="genome.", suffix=".masked.gff3", dir=wd)
    cmd = BEDTOOLS_SORT % gff3
    cmd1 = BEDTOOLS_MERGE
    cmd2 = AWK % length
    if verbose:
        print(cmd)
        print(cmd1)
        print(cmd2)
    bedsort = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    bedmerge = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stdin=bedsort.stdout, shell=True)
    awk = subprocess.Popen(cmd2, stdin=bedmerge.stdout, stdout=outputmerge, shell=True)
    awk.communicate()

    masked = ref + ".masked.fasta"
    cmd = BEDTOOLS_MASK % (ref, outputmerge.name, masked)
    if verbose:
        print(cmd)
    maskfasta=subprocess.Popen(cmd, cwd=wd, shell=True)
    maskfasta.communicate()
    return masked


def repeatsfind(genome, working_dir, repeat_lenght, threads_use, verbose):

    name_gff = genome.split("/")[-1] + ".out.gff"
    gff_path = Path(working_dir + "/" + genome.split("/")[-1] + ".out.gff")

    if gff_path.is_file():
        gff = [y for x in os.walk(working_dir) for y in glob(os.path.join(x[0], name_gff))][0]
    else:
        freq_file = working_dir + genome.split("/")[-1] + ".freq"
        log = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir)
        err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir)

        cmd = BUILD_TABLE % (genome, freq_file)
        if verbose:
            print(cmd)
        build = subprocess.Popen(cmd, cwd=working_dir, stdout=log, stderr=err, shell=True)
        build.communicate()

        fasta_out = working_dir + genome.split("/")[-1] + ".repeats.fasta"
        log = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir)
        err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir)

        cmd = REPEAT_SCOUT % (genome, fasta_out, freq_file)
        if verbose:
            print(cmd)
        scout = subprocess.Popen(cmd, cwd=working_dir, stdout=log, stderr=err, shell=True)
        scout.communicate()

        log = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir)
        err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir)

        cmd = REPEAT_MASKER % (genome, fasta_out, str(threads_use), working_dir)
        if verbose:
            print(cmd)
        mask = subprocess.Popen(cmd, cwd=working_dir, stdout=log, stderr=err, shell=True)
        mask.communicate()
        gff = [y for x in os.walk(working_dir) for y in glob(os.path.join(x[0], name_gff))][0]

    genome_masked = maskedgenome(working_dir, genome, gff, repeat_lenght, verbose)
    return genome_masked


if __name__ == '__main__':
    #filterLongReads(fastq_filename, min_length, max_length, wd, adapter, threads, a, alignm_score_value)
    filterLongReads(*sys.argv[1:])