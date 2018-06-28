#!/usr/bin/env python3
import datetime
import os
import subprocess
import sys
import tempfile
from glob import glob
from pathlib import Path

import align as align
from Bio import SeqIO

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


def filterLongReads(fastq_filename, min_length, max_length, wd, adapter, threads, align_score_value, a):
    """
    Filters out reads longer than length provided and it is used to call the alignemnt and parse the outputs
    """
    scoring = [3, -6, -5, -2]

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

        return out_filename_oriented
    else:
        sizes = [rec.id for rec in SeqIO.parse(out_filename, "fasta")]
        fmtdate = '%H:%M:%S %d-%m'
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write("###FINISHED FILTERING AT:\t" + now + "###\n\n###LOREAN KEPT\t\033[32m" + str(len(sizes)) +
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