#!/usr/bin/python3

import os
from Bio import SeqIO

count_sequences = 0
length_cluster = 0


def parse_only(threshold_float, tmp_wd):
    """
    to join the assembly and the parsing process
    :param threshold_float:
    :param tmp_wd:
    :return:
    """
    evm_list = []
    wd = tmp_wd + 'consensus/tmp/'
    for root, dirs, _ in os.walk(wd):
        for direc in dirs:
            if 'output' in direc:
                outputDir = wd + direc + '/'
                if outputDir:
                    parse_contigs(outputDir, threshold_float)
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


def parse_contigs(outputAssembly, threshold_float):
    """
    Parses the output from iAssembler, to a single FASTA file
    :param outputAssembly: 
    :param threshold_float: 
    :return:
    """
    fname = outputAssembly + 'contig_member'
    fname_exists = os.path.isfile(fname)
    if fname_exists:
        # part all the files in the tmp assembly folder
        contigInfo = open(outputAssembly + 'contig_member', 'r')
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
        fileAssembly = outputAssembly + 'unigene_seq.new.fasta'
        contigSeq = open(fileAssembly, 'r')
        contigDict = SeqIO.to_dict(SeqIO.parse(contigSeq, 'fasta'))
        output_filename = outputAssembly[:-1] + '_assembled.fasta'
        outputFile = open(output_filename, 'w')
        for iden, write_iden in list(real_contigs.items()):
            if iden in contigDict:
                outputFile.write('>' + write_iden + '\n' +
                                 str(contigDict[iden].seq) + '\n')
        contigSeq.close()
        outputFile.close()
    return


def catAssembled(wd):
    """
    collect the assembled contigs and generate a multifasta file
    :param wd: 
    :return: 
    """
    print('\t###GENERATE FASTA FILE FROM CONTIGS###\n')
    '''C at all the assembled single fasta files in to a uniq file'''
    wd_tmp = wd + 'consensus/tmp/'
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


def addEVM(whole_fasta_name, output_filename, evm_nosupport, output_merged_fasta_name):
    """
    this module looks for genes that were not used in the consensus stage. usually are gene models without long reads
    support
    :param whole_fasta_name: 
    :param output_filename: 
    :param evm_nosupport:
    :param output_merged_fasta_name: 
    :return: 
    """
    print('\t###APPEND EVM NOT USED FROM CONTIGS BUILDING###\n')
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
    seq_count = 1
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


def getEVMnoUnitig(target_wd):
    """
    this module change the name of the assembled contigs in evm names and replace the unigene identifiers
    :param target_wd: 
    :return: 
    """
    print('\t###EXTRACT EVM NAME FROM ASSEMBLED CONTIGS###\n')
    '''Gets the name of evm prediction in the assembly that do not have support'''
    evm_list = []
    for root, dirs, _ in os.walk(target_wd):
        for direc in dirs:
            t_filename = root + direc + '/contig_member'
            try:
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
    return evm_list
