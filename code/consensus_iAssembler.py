#!/usr/bin/env python3

import subprocess
import os
import re
from Bio import SeqIO
from queue import Queue
from threading import Thread

count_sequences = 0
length_cluster = 0


def gff_read(gff3_file, reference, working_dir):
    """
    Runs gffread on a gff3 file to produce fasta files with the
    matching records
    :param gff3_file:
    :param reference:
    :param working_dir:
    :return:
    """
    out_name = working_dir + 'getFasta.fasta'
    args = [
        'bedtools',
        'getfasta',
        '-fi',
        reference,
        '-bed',
        gff3_file,
        '-fo',
        out_name,
        '-name',
        '-split']

    if os.path.isfile(out_name):
        print((
            'gff3read file existed already: ' +
            out_name +
            ' --- skipping\n'))
        return out_name

    try:
        subprocess.check_call(args)
    except:
        raise NameError('')
    return out_name


def cluster_pipeline(gff3_file, merge_distance, strand):
    """
    here the cluseter of sequence from the same locus are prepared
    :param gff3_file:
    :param merge_distance:
    :param strand:
    :return:
    """

    btsort1 = ['bedtools', 'sort', '-i', gff3_file]
    dist = '-' + str(merge_distance)
    if strand:
        btmerge1 = [
            'bedtools',
            'merge',
            '-s',
            '-d',
            str(dist),
            '-c',
            '4,4',
            '-o',
            'count,distinct']

        print("\t###CLUSTERING IN STRANDED MODE###\n")
    else:
        btmerge1 = [
            'bedtools',
            'merge',
            '-d',
            str(dist),
            '-c',
            '4,4',
            '-o',
            'count,distinct']
        print("\t###CLUSTERING IN NON-STRANDED MODE###\n")
    btsort2 = ['bedtools', 'sort']

    # Sort the GFF3 file
    btsort1_call = subprocess.Popen(btsort1, stdout=subprocess.PIPE)
    # Merge the BED entries, count number of reads on each merged entry
    btmerge1_call = subprocess.Popen(btmerge1, stdin=btsort1_call.stdout, stdout=subprocess.PIPE)
    # NSort it again and returns
    btsort2_call = subprocess.Popen(btsort2, stdin=btmerge1_call.stdout, stdout=subprocess.PIPE)
    outputBT = btsort2_call.communicate()[0]
    final_output = outputBT.splitlines()
    return final_output


def fasta2Dict(fasta_filename):
    """
    Prepare a dictionary of all the sequences that is used together with
    the fasta file to make single fasta files for the assembly
    :param fasta_filename:
    :return:
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
    :param count:
    :param bedline:
    :param fasta_dict:
    :param min_length:
    :param min_evidence:
    :param max_evidence:
    :param wd:
    :return:
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


def iAssembler(fasta_file, overlap_length, percent_identity, wd):
    """
    Call iAssembler to assemble every cluster in fasta_list
    :param fasta_file:
    :param overlap_length:
    :param percent_identity:
    :param wd:
    :return:
    """
    newFasta = fasta_file.split('/')[-1]
    args = ['iAssembler.pl', '-i', newFasta, '-h', str(overlap_length),
            '-p', str(percent_identity), '-o', newFasta + '_output', '2> ']
    outputDir = wd + newFasta + '_output/'  # whole path
    try:
        subprocess.call(
            args,
            cwd=wd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        # TRY TO CHANGE TO SUBPROCESSPOPEN
    except:
        return False
    return outputDir


def assembleParse(queue, overlap_length, percent_identity, wd):
    """to join the assembly and the parsing process"""
    while True:
        try:
            clusterFasta = queue.get()
        except:
            break
        queue.task_done()
        global count_sequences
        count_sequences += 1
        global length_cluster
    return


def func_star(a_b):
    return assembleParse(*a_b)


def assembly(clusterList, fasta_dict, min_evidence, max_evidence, overlap_length,
             percent_identity, threads, wd):
    """
    handles the assembly process and parsing in a multithreaded way
    :param clusterList:
    :param fasta_dict:
    :param min_evidence:
    :param max_evidence:
    :param overlap_length:
    :param percent_identity:
    :param threads:
    :param wd:
    :return:
    """

    assembly_queue = Queue()
    cluster_counter = 1
    for record in clusterList:
        # Write fasta for each cluster
        clusterFasta = write_fastas(
            cluster_counter,
            record,
            fasta_dict,
            overlap_length,
            min_evidence,
            max_evidence,
            wd)
        cluster_counter += 1
        if clusterFasta:
            # If it is a cluster to assemble put it on queue
            assembly_queue.put(clusterFasta)
        else:
            continue
    global length_cluster
    length_cluster = str(assembly_queue.qsize())
    for i in range(int(threads)):
        t = Thread(
            target=assembleParse,
            args=(
                assembly_queue,
                overlap_length,
                percent_identity,
                wd))
        t.daemon = True
        t.start()
    assembly_queue.join()
