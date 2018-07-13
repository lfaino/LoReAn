#!/usr/bin/env python3

import os
import re
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool, Manager

import numpy as np
import progressbar
from Bio import SeqIO

#==========================================================================================================
# COMMANDS LIST

ASSEMBLY  = 'iAssembler.pl -i %s -h %s  -p %s -o %s_output  2> %s.log'

BEDTOOLS_GETFASTA = 'bedtools getfasta -fi %s -bed %s -fo %s -name -split'

BEDTOOLS_MERGE_ST = 'bedtools merge -s -c 4,4 -o count,distinct'

BEDTOOLS_MERGE = 'bedtools merge  -c 4,4 -o count,distinct'

CAT = 'cat %s'

BEDTOOLS_SORT = 'bedtools sort'

#==========================================================================================================

count_sequences = 0
length_cluster = 0

def gffread(gff3_file, reference, working_dir, verbose):
    """
    Runs gffread on a gff3 file to produce fasta files with the
    matching records
    """
    out_name = working_dir + 'getFasta.fasta'
    cmd = BEDTOOLS_GETFASTA % (reference, gff3_file, out_name)
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n\n' % cmd)
            log = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir, prefix="startUser.", suffix=".out")
            err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=working_dir, prefix="startUser.", suffix=".out")
        else:
            log = tempfile.NamedTemporaryFile(mode='w', dir=working_dir, prefix="startUser.", suffix=".out")
            err = tempfile.NamedTemporaryFile(mode='w', dir=working_dir, prefix="startUser.", suffix=".out")

        bedtools = subprocess.Popen(cmd, cwd=working_dir, stdout=log, stderr=err, shell=True)
        bedtools.communicate()
    except:
        raise NameError('')
    return out_name


def cluster_pipeline(gff3_file, strand, verbose):
    """
    here clusters of sequences from the same locus are prepared
    """

    cat = CAT % gff3_file
    btsort1 = BEDTOOLS_SORT

    if strand:
        btmerge1 = BEDTOOLS_MERGE_ST
        sys.stdout.write("\t ###CLUSTERING IN\033[32m STRANDED MODE\033[0m###\n")

    else:
        btmerge1 = BEDTOOLS_MERGE
        sys.stdout.write("\t###CLUSTERING IN\033[32m NON-STRANDED MODE\033[0m ###\n")

    btsort2 = BEDTOOLS_SORT
    # Sort the GFF3 file
    cat_call = subprocess.Popen(cat, stdout=subprocess.PIPE, shell=True)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % cat)
    btsort1_call = subprocess.Popen(btsort1, stdin=cat_call.stdout, stdout=subprocess.PIPE, shell=True)
    # Merge the BED entries, count number of reads on each merged entry
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % btsort1)
    btmerge1_call = subprocess.Popen(btmerge1, stdin=btsort1_call.stdout, stdout=subprocess.PIPE, shell=True)
    # NSort it again and returns
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % btmerge1)
    btsort2_call = subprocess.Popen(btsort2, stdin=btmerge1_call.stdout, stdout=subprocess.PIPE, shell=True)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % btsort2)
    outputBT = btsort2_call.communicate()[0]
    final_output = outputBT.splitlines()

    return final_output


def fasta2Dict(fasta_filename):
    """
    Prepare a dictionary of all the sequences that is used together with
    the fasta file to make single fasta files for the assembly
    """
    fasta_file = open(fasta_filename, 'r')
    fasta_dict2 = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    fasta_dict = {}
    for key, seq2 in list(fasta_dict2.items()):
        seq = str(seq2.seq).replace("N", "")
        fasta_dict[key] = seq
        del fasta_dict2[key]
    fasta_file.close()
    return fasta_dict


def write_fastas(count, bedline, fasta_dict, min_length, min_evidence, max_evidence, wd):
    """
    From the output list of the pipeline, recovers the ID and goes back to the
    fasta file to retrieve the sequence
    """


    line = (bedline.decode("utf-8")).split('\t')
    if len(line) == 6:
        chrm, start, end, strand, number, idents = (line[0], line[1], line[2], line[3], line[4], line[5])
    elif len(line) == 5:
        chrm, start, end, number, idents = (line[0], line[1], line[2], line[3], line[4])
    ids_short = []
    if len(idents.split(',')) > int(min_evidence) and len(idents.split(',')) < int(max_evidence):
        for element in re.split(',|;', idents):
            ids_short.append(element)
    else:
            return False


    unique_ids = list(set(ids_short))
    cluster_filename = wd + '_'.join([chrm, start, end]) + '_' + str(count) + '.fasta'
    cluster_file = open(cluster_filename, 'w')
    read_count = 1

    for iden in unique_ids:
        if iden in fasta_dict:
            if len(str(fasta_dict[iden])) > int(min_length):
                if len(iden) < 40:
                    cluster_file.write('>' + iden + '\n' +  str(fasta_dict[iden]) + '\n')
                else:
                    cluster_file.write('>' + str(read_count) + '\n' + str(fasta_dict[iden]) + '\n')
                    read_count += 1
                del fasta_dict[iden]
            else:
                del fasta_dict[iden]

    cluster_file.close()
    cluster_file = open(cluster_filename, 'r')
    nlines = 0
    for line in cluster_file:
        if line.startswith('>'):
            nlines += 1
    if nlines < int(min_evidence) or nlines > int(max_evidence):
        return False

    cluster_file.close()
    return cluster_filename


def generate_fasta(clusterList, fasta_dict, min_evidence, max_evidence, overlap_length, strand, wd):
    """write fasta clusters
    """

    if min_evidence == "" and strand:
        array_perc = [float(line.decode().split("\t")[4]) for line in clusterList]
        min_evidence = np.percentile(array_perc, 75)
        print("\n" + "\033[32m ### LOREAN SET THE MIN READS SUPPORT FOR A CLUSTER TO " + str(min_evidence) + " AUTOMATICALLY ### \n")
        print('\033[0m')
    elif min_evidence == "":
        array_perc = [float(line.decode().split("\t")[3]) for line in clusterList]
        min_evidence = np.percentile(array_perc, 75)
        print("\n" + "\033[32m ### LOREAN SET THE MIN READS SUPPORT FOR A CLUSTER TO " + str(min_evidence) + " AUTOMATICALLY ### \n")
        print('\033[0m')
    else:
        print("\n" + "\033[32m ### LOREAN SET THE MIN READS SUPPORT FOR A CLUSTER TO " + str(min_evidence) + " ### \n")
        print('\033[0m')

    cluster_counter = 1
    for record in clusterList:
        # Write fasta for each cluster
        write_fastas(cluster_counter, record, fasta_dict, overlap_length, min_evidence, max_evidence, wd)
        cluster_counter += 1


def assembly(overlap_length, percent_identity, threads, wd, verbose):
    """
    """
    manage = Manager()
    queue = manage.Queue()
    pool = Pool(processes=int(threads), maxtasksperchild=10)

    new_commands = []
    for root, dirs, file in os.walk(wd):
        for fasta_file in file:
            complete_data = (fasta_file, percent_identity, overlap_length, wd, verbose, queue)
            new_commands.append(complete_data)
    results = pool.map_async(iAssembler, new_commands, chunksize=1)
    with progressbar.ProgressBar(max_value=len(new_commands)) as bar:
        while not results.ready():
            size = queue.qsize()
            bar.update(size)
            time.sleep(1)

def iAssembler(new_commands):
    """
    Call iAssembler to assemble every cluster in fasta_list
    """

    cmd = ASSEMBLY %(new_commands[0], new_commands[2], new_commands[1], new_commands[0], new_commands[0])
    new_commands[5].put(cmd)
    outputDir = new_commands[3] + new_commands[0] + '_output/'  # whole path
    log_name = new_commands[3] + "Assembly.log"
    log = open(log_name, 'w')

    try:
        if new_commands[4]:
            sys.stderr.write('Executing: %s\n\n' % cmd)
        assembly = subprocess.Popen(cmd, cwd=new_commands[3], stderr = log, stdout = log, shell=True)
        assembly.communicate()
    except:
        return False
    log.close()
    return outputDir


if __name__ == '__main__':
    cluster_pipeline(*sys.argv[1:])