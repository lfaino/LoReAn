#!/usr/bin/env python3

import os
import re
import subprocess
import sys
import time
from multiprocessing import Pool, Manager

import progressbar
from Bio import SeqIO

#==========================================================================================================
# COMMANDS LIST

ASSEMBLY  = 'iAssembler.pl -i %s -h %s  -p %s -o %s_output  2> %s.log'

BEDTOOLS_GETFASTA =  'bedtools getfasta -fi %s -bed %s -fo %s -name -split'

BEDTOOLS_MERGE_ST = 'bedtools merge -s -d %s -c 4,4 -o count,distinct'

BEDTOOLS_MERGE = 'bedtools merge -d %s -c 4,4 -o count,distinct'

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
        bedtools = subprocess.Popen(cmd, shell=1)
        bedtools.communicate()
    except:
        raise NameError('')
    return out_name

def cluster_pipeline(gff3_file, merge_distance, strand):
    """
    here the clusters of sequence from the same locus are prepared
    """

    cat = CAT %(gff3_file)
    btsort1 = BEDTOOLS_SORT

    dist = '-' + str(merge_distance)
    if strand:
        btmerge1 = BEDTOOLS_MERGE_ST % (str(dist))
        sys.stdout.write("\t###CLUSTERING IN STRANDED MODE###\n")

    else:
        btmerge1 = BEDTOOLS_MERGE_ST % (str(dist))
        sys.stdout.write("\t###CLUSTERING IN NON-STRANDED MODE###\n")

    btsort2 = BEDTOOLS_SORT

    # Sort the GFF3 file
    cat_call = subprocess.Popen(cat, stdout=subprocess.PIPE, shell =1)
    btsort1_call = subprocess.Popen(btsort1, stdin=cat_call.stdout, stdout=subprocess.PIPE, shell =1)
    # Merge the BED entries, count number of reads on each merged entry
    btmerge1_call = subprocess.Popen(btmerge1, stdin=btsort1_call.stdout, stdout=subprocess.PIPE, shell =1)
    # NSort it again and returns
    btsort2_call = subprocess.Popen(btsort2, stdin=btmerge1_call.stdout, stdout=subprocess.PIPE, shell =1)
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

    global idents, end, start, chrm
    line = (bedline.decode("utf-8")).split('\t')
    if len(line) == 6:
        chrm, start, end, strand, number, idents = (
            line[0], line[1], line[2], line[3], line[4], line[5])
    elif len(line) == 5:
        chrm, start, end, number, idents = (
            line[0], line[1], line[2], line[3], line[4])

    ids_short = []
    if len(idents.split(',')) > int(min_evidence) and len(
            idents.split(',')) < int(max_evidence):
        for element in re.split(',|;', idents):
            # if 'ID=' in element:
            # We keep read.mrna and evm.model records
            ids_short.append(element)
    else:
        return False

    unique_ids = list(set(ids_short))
    clusterFilename = wd + \
                      '_'.join([chrm, start, end]) + '_' + str(count) + '.fasta'
    clusterFile = open(clusterFilename, 'w')
    read_count = 1

    for iden in unique_ids:
        if iden in fasta_dict:
            if len(str(fasta_dict[iden])) > int(min_length):
                if len(iden) < 40:
                    clusterFile.write('>' + iden + '\n' +
                                      str(fasta_dict[iden]) + '\n')
                else:
                    clusterFile.write('>' + str(read_count) +
                                      '\n' + str(fasta_dict[iden]) + '\n')
                    read_count += 1
                del fasta_dict[iden]
            else:
                del fasta_dict[iden]

    clusterFile.close()
    clusterFile = open(clusterFilename, 'r')
    nlines = 0
    for line in clusterFile:
        if line.startswith('>'):
            nlines += 1
    if nlines < int(min_evidence) or nlines > int(max_evidence):
        return False

    clusterFile.close()
    return clusterFilename

def generate_fasta(clusterList, fasta_dict, min_evidence, max_evidence, overlap_length, wd):
    """write fasta clusters
    """
    cluster_counter = 1
    for record in clusterList:
        # Write fasta for each cluster
        write_fastas(
            cluster_counter,
            record,
            fasta_dict,
            overlap_length,
            min_evidence,
            max_evidence,
            wd)
        cluster_counter += 1

def assembly(overlap_length, percent_identity, threads, wd, verbose):
    """
    """
    manage = Manager()
    queue = manage.Queue()
    pool = Pool(int(threads))

    new_commands = []
    for root, dirs, file in os.walk(wd):
        for fasta_file in file:
            complete_data = (fasta_file, percent_identity, overlap_length, wd, verbose, queue)
            new_commands.append(complete_data)
    #with Pool(int(threads)) as pool:
    results = pool.map_async(iAssembler, new_commands)
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
        assembly = subprocess.Popen(cmd, cwd=new_commands[3], stderr = log, stdout = log, shell=1)
        assembly.communicate()
    except:
        return False
    log.close()
    return outputDir
