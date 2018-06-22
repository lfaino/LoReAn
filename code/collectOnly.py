#!/usr/bin/python3

import gffutils
import gffutils.gffwriter as gffwriter
import os
import subprocess
import sys
import tempfile
from Bio import SeqIO

count_sequences = 0
length_cluster = 0

#==========================================================================================================
# COMMANDS LIST

GFFREAD_W = 'gffread -g %s -W -w %s %s'

GFFREAD_M = 'gffread -M -F -o %s %s'

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
        contig_info = open(output_assembly + 'contig_member', 'r')
        contigs = {}
        total_reads = 0
        for line in contig_info:
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
        contig_info.close()
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
        file_assembly = output_assembly + 'unigene_seq.new.fasta'
        if os.path.exists(file_assembly):
            contig_seq = open(file_assembly, 'r')
            contig_dict = SeqIO.to_dict(SeqIO.parse(contig_seq, 'fasta'))
            output_filename = output_assembly[:-1] + '_assembled.fasta'
            output_file = open(output_filename, 'w')
            for iden, write_iden in list(real_contigs.items()):
                if iden in contig_dict:
                    output_file.write('>' + write_iden + '\n' + str(contig_dict[iden].seq) + '\n')
            contig_seq.close()
            output_file.close()
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
                if os.path.exists(wd_fasta):
                    fasta_sequences = SeqIO.parse(open(wd_fasta),'fasta')
                    for fasta in fasta_sequences:
                        fasta_name = ">assembled_" + str(count_seq_ass)
                        count_seq_ass += 1
                        cdna = fasta_name + "\n" + str(fasta.seq) + "\n"
                        out_file.write(cdna)
    return fileName


def add_EVM(final_update, wd, consensus_mapped_gff3):

    """
    """

    db_evm = gffutils.create_db(final_update, ':memory:', merge_strategy='create_unique', keep_order=True)
    ids_evm = [gene.attributes["ID"][0] for gene in db_evm.features_of_type("mRNA")]

    db_gmap = gffutils.create_db(consensus_mapped_gff3, ':memory:', merge_strategy='create_unique', keep_order=True)
    ids_gmap_full = [gene.attributes["ID"][0] for gene in db_gmap.features_of_type("gene")]
    ids_gmap = [gene.attributes["ID"][0].split("_")[0] for gene in db_gmap.features_of_type("gene")]

    uniq_evm = [evm for evm in ids_evm if not evm in ids_gmap]

    mRNA = []
    for evm in uniq_evm:
        for line in db_evm.parents(evm, order_by='start'):
            mRNA.append(line.attributes["ID"][0])
    mRNA_uniq = list(set(mRNA))
    outfile = tempfile.NamedTemporaryFile(delete=False, prefix="additional.1.", suffix=".gff3", dir=wd)
    gff_out_s = gffwriter.GFFWriter(outfile.name)

    for name in mRNA_uniq:
        for i in db_evm.children(name, order_by='start'):
            gff_out_s.write_rec(i)
        gff_out_s.write_rec(db_evm[name])
    for name in ids_gmap_full:
        for i in db_gmap.children(name, order_by='start'):
            gff_out_s.write_rec(i)
        gff_out_s.write_rec(db_gmap[name])
    gff_out_s.close()

    return outfile.name
    #return outfile_out.name


def transform_func(x):

     if 'locus' in x.featuretype:
         x.featuretype = "gene"
         x.attributes['ID'] = x.attributes['locus']
         x.source = "LoReAn"
         return x
     elif 'mRNA' in x.featuretype:
         x.attributes['Parent'] = x.attributes['locus']
         x.source = "LoReAn"
         return x
     else:
         x.source = "LoReAn"
         return x


if __name__ == '__main__':
    add_EVM(*sys.argv[1:])
