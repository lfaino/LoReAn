#!/usr/bin/python3

import os
import sys

from Bio import SeqIO

count_sequences = 0
length_cluster = 0



def parse_only(threshold_float, wd, verbose):
    """
    to join the assembly and the parsing process
    """
    evm_list = []
    if verbose:
        sys.stderr.write('Executing: Parse assembled consensus and EVM\n')
    for root, dirs, _ in os.walk(wd):
        for direc in dirs:
            if 'output' in direc:
                outputDir = wd + direc + '/'
                if outputDir:
                    parse_contigs(outputDir, threshold_float, verbose)
            try:
                t_filename = root + direc + '/contig_member'
                t_file = open(t_filename, 'r')
                for line in t_file:
                    line = line.strip()
                    elements = line.split('\t')
                    for el in elements:
                        if 'evm' in el:
                            evm_list.append(el)
                t_file.close()
            except IOError:
                continue
        global count_sequences
        count_sequences += 1
        global length_cluster
    return evm_list


def parse_contigs(output_assembly, threshold_float, verbose):
    """
    Parses the output from iAssembler, to a single FASTA file
    """

    if verbose:
        sys.stderr.write('Executing: Parse assembled consensus\n')
    fname = output_assembly + 'contig_member'
    fname_exists = os.path.isfile(fname)
    if fname_exists:
        # part all the files in the tmp assembly folder
        contigInfo = open(output_assembly + 'contig_member', 'r')
        contigs = {}
        total_reads = 0
        for line in contigInfo:
            line = line.strip()
            line = line.split('\t')
            # count the reads for each assembled Unitig
            read_number = len(line) - 1
            for element in line:
                if 'evm' in element:
                    contigs[line[0]] = [read_number, element]
            if line[0] not in contigs:
                contigs[line[0]] = [read_number]
            # sum all the reads whitin one cluster
            total_reads += read_number
        contigInfo.close()
        # calculate the number of reads considering the threshold
        threshold = total_reads * float(threshold_float)
        global count_sequences
        real_contigs = {}
        count_unitigs = 1
        for key, element in list(contigs.items()):
            # to retrieve only supported assembly
            if len(element) == 2:
                if element[0] >= threshold and 'evm' in element[1]:
                    real_contigs[key] = element[1] + ' ' + str(
                        element[0]) + '_above_threshold_' + str(threshold) + ' loc_' + str(count_sequences)
                elif element[0] < threshold and 'evm' in element[1]:
                    real_contigs[key] = element[1] + ' ' + str(
                        element[0]) + '_below_threshold_' + str(threshold) + ' loc_' + str(count_sequences)
            elif len(element) == 1:
                if element[0] >= threshold:
                    real_contigs[key] = 'Unitig' + str(count_sequences) + '_' + str(count_unitigs) + ' ' + str(
                        element[0]) + '_above_threshold_' + str(threshold) + ' loc_' + str(count_sequences)
                    count_unitigs += 1
        # writes the outputs
        fileAssembly = output_assembly + 'unigene_seq.new.fasta'
        contigSeq = open(fileAssembly, 'r')
        contigDict = SeqIO.to_dict(SeqIO.parse(contigSeq, 'fasta'))
        output_filename = output_assembly[:-1] + '_assembled.fasta'
        outputFile = open(output_filename, 'w')
        for iden, write_iden in list(real_contigs.items()):
            if iden in contigDict:
                outputFile.write('>' + write_iden + '\n' + str(contigDict[iden].seq) + '\n')
        contigSeq.close()
        outputFile.close()
    return


def cat_assembled(wd):
    """
    collect the assembled contigs and generate a multifasta file
    """
    sys.stdout.write('\t###GENERATE FASTA FILE FROM CONTIGS###\n')
    wd_tmp = wd
    fileName = wd_tmp + 'assembly.fasta'
    testFasta = open(fileName, 'w')
    for root, dirs, files in os.walk(wd_tmp):
        for name in files:
            wd_fasta = os.path.join(root, name)
            if 'assembled.fasta' in wd_fasta:
                t_file = open(wd_fasta, 'r')
                for line in t_file:
                    testFasta.write(line)
                t_file.close()

    testFasta.close()
    return fileName


def cat_assembled_all(wd):
    """
    collect the assembled contigs and generate a multifasta file
    """
    sys.stdout.write('\t###GENERATE FASTA FILE FROM CONTIGS###\n')
    wd_tmp = wd
    fileName = wd_tmp + 'assembly.all.fasta'
    out_file = open(fileName, "w")
    count_seq_ass = 0
    for root, dirs, files in os.walk(wd_tmp):
        for name in dirs:
            wd_dir = os.path.join(root, name)
            if 'fasta_output' in wd_dir:
                wd_fasta = os.path.join(wd_dir, "unigene_seq.new.fasta")
                fasta_sequences = SeqIO.parse(open(wd_fasta),'fasta')
                for fasta in fasta_sequences:
                    fasta_name = ">assembled_" + str(count_seq_ass)
                    count_seq_ass += 1
                    cdna = fasta_name + "\n" + str(fasta.seq) + "\n"
                    out_file.write(cdna)
    return fileName


def add_EVM(whole_fasta_name, output_filename, output_merged_fasta_name):
    """
    this module looks for genes that were not used in the consensus stage. usually are gene models without long reads
    support
    """
    sys.stdout.write('\t###APPEND EVM NOT USED FROM CONTIGS BUILDING###\n')
    '''Adds the EVM records that are not present in the final contig evidence'''
    whole_fasta = open(whole_fasta_name, 'r')
    out_fasta_file = open(output_filename, 'r')
    outputMerged = open(output_merged_fasta_name, 'w')
    wholeDict = SeqIO.to_dict(SeqIO.parse(whole_fasta, 'fasta'))
    count = 0
    dictOut = {}
    outFasta = SeqIO.parse(out_fasta_file, 'fasta')
    for record in outFasta:
        if record.id in dictOut:
            dictOut[str(record.id) + '_' + str(count)] = str(record.seq)
            count += 1
        else:
            dictOut[record.id] = str(record.seq)
    for key in list(wholeDict.keys()):
        if 'evm' in key and key not in dictOut:
            ident = '>Gene' + str(count) + '_' + key
            outputMerged.write(
                ident + '\n' + str(wholeDict[key].seq) + '\n')
            count += 1
    for key, element in list(dictOut.items()):
        ident = '>Gene' + str(count) + '_' + key
        outputMerged.write(ident + '\n' + str(element) + '\n')
        count += 1

    whole_fasta.close()
    outFasta.close()
    outputMerged.close()

