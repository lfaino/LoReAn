#!/usr/bin/python3

import os
import re
import subprocess
import sys
import tempfile

from Bio import SeqIO

count_sequences = 0
length_cluster = 0

#==========================================================================================================
# COMMANDS LIST

GFFREAD_W = 'gffread -g %s -W -w %s %s'

#==========================================================================================================


def parse_only(threshold_float, wd, verbose):
    """
    to join the assembly and the parsing process
    """
    evm_list = []
    if verbose:
        sys.stderr.write('Executing: Parse assembled consensus and EVM\n')
    for root, dirs, _ in os.walk(wd):
        for direc in dirs:
            if direc.endswith('output'):
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
        location_list = output_assembly.split("/")[-2].split(".")[0].split("_")[:-1]
        location = "_".join(location_list)
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
                    real_contigs[key] = element[1] + '_' + location + ' ' + str(
                        element[0]) + '_above_threshold_' + str(threshold) + ' loc_' + str(count_sequences)
                elif element[0] < threshold and 'evm' in element[1]:
                    real_contigs[key] = element[1] + '_' + location + ' ' + str(
                        element[0]) + '_below_threshold_' + str(threshold) + ' loc_' + str(count_sequences)
            elif len(element) == 1:
                if element[0] >= threshold:
                    real_contigs[key] = 'Unitig' + str(count_sequences) + '_' + str(count_unitigs) + '_' + location + ' ' + str(
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
            if wd_fasta.endswith('_assembled.fasta'):
                input_file = open(wd_fasta)
                fasta_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
                evm = [key for key in fasta_dict if "evm" in key]
                above = [key for key in fasta_dict if "above" in fasta_dict[key].description]
                if len(fasta_dict) > 1:
                    if len(evm) > 0 and len(above) > 0:
                        if evm[0] in above:
                            SeqIO.write(fasta_dict[evm[0]], testFasta, "fasta")
                        else:
                            SeqIO.write(fasta_dict[evm[0]], testFasta, "fasta")
                            SeqIO.write(fasta_dict[above[0]], testFasta, "fasta")
                    elif len(evm) > 1:
                        SeqIO.write(fasta_dict[evm[0]], testFasta, "fasta")
                    elif len(above) > 1:
                        SeqIO.write(fasta_dict[above[0]], testFasta, "fasta")

                elif len(evm) > 0:
                    SeqIO.write(fasta_dict[evm[0]], testFasta, "fasta")
                elif len(above) > 0:
                    SeqIO.write(fasta_dict[above[0]], testFasta, "fasta")
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


def add_EVM(tmp_assembly, final_update, ref, consensus_wd, verbose):

    """
    this module looks for genes that were not used in the consensus stage. usually are gene models without long reads
    support
    """

    out = tempfile.NamedTemporaryFile(delete=False, prefix="gffreads", dir=consensus_wd)
    err = tempfile.NamedTemporaryFile(delete=False, prefix="gffreads", dir=consensus_wd)
    log = tempfile.NamedTemporaryFile(delete=False, prefix="gffreads", dir=consensus_wd)

    merged_fasta_filename = consensus_wd + 'assembly.wEVM.fasta'

    sys.stdout.write('\t###APPEND EVM NOT USED FROM CONTIGS BUILDING###\n')


    com = GFFREAD_W % (ref, out.name, final_update)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=log, cwd = consensus_wd, stderr=err, shell=True)
    call.communicate()

    evm_common = []
    evm_all = []
    evm_fasta = open(out.name, 'r')
    evm_dict = SeqIO.to_dict(SeqIO.parse(evm_fasta, 'fasta'))
    for key in evm_dict:
        evm_all.append(key)
    new_dict = {}
    count = 0

    out_fasta_file = open(tmp_assembly, 'r')
    assembly_dict = SeqIO.to_dict(SeqIO.parse(out_fasta_file, 'fasta'))
    for record in assembly_dict:
        if record.startswith('evm'):
            original = record.split('_')[:-3]
            if original[0] in evm_dict:
                evm_common.append(original[0])
                new_dict[record] = assembly_dict[record]
        else:
            count += 1
            name = 'Gene' + str(count) + '_' + assembly_dict[record].id
            assembly_dict[record].id = name
            new_dict[record] = assembly_dict[record]


    diff = list(set(evm_common)^set(evm_all))

    with open(merged_fasta_filename, 'w') as fh:
        SeqIO.write(new_dict.values(), fh, 'fasta')
        for item in diff:
            desc = evm_dict[item].description.split(" ")[2]
            evm_dict[item].description = ""
            loc = re.split(":|\||-", desc)
            region = "_".join([loc[1], loc[2], loc[3]])
            id = evm_dict[item].id + "_" + region
            evm_dict[item].id = id
            SeqIO.write(evm_dict[item], fh, 'fasta')
    return merged_fasta_filename


if __name__ == '__main__':
    add_EVM(*sys.argv[1:])
