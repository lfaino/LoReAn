from __future__ import division
from multiprocessing import Pool
import sys
import subprocess
import argparse
import os
import re
import shutil
from dirs_and_files import check_create_dir
from Bio import SeqIO
from Queue import Queue
from threading import Thread, Lock
import itertools


def lengthSupport(tmp_wd):
    '''with this function, we try to reduce the UTRs if they have low support (at least 10% of the total reads)
    we do this by using the information from the iAssembly outputs'''
    wd = tmp_wd  # +'consensus/tmp/'
    print '\n\t###FILTERING ASSEMBLED GENE REGIONS WITH HIGH SUPPORT###\n'
    for root, dirs, _ in os.walk(wd):
        for direc in dirs:
            if 'output' in direc:
                outputDir = wd + '/' + direc + '/'
                member = outputDir + 'contig_member'
                fname_exists = os.path.isfile(member)
                if fname_exists:
                    contigs = open(member, 'r')
                    fasta_seq = outputDir + 'unigene_seq.fasta'
                    fasta_seq_new = outputDir + 'unigene_seq.new.fasta'
                    outFasta = SeqIO.parse(fasta_seq, 'fasta')
                    name_list = []
                    unitigDict = {}
                    origDict = {}
                    for line in contigs:
                        line = line.strip()
                        line = line.split('\t')
                        if len(line) > 10:
                            unitigDict[line[0]] = len(line) - 1
                            name_list.append(line[0])
                        else:
                            origDict[line[0]] = len(line) - 1
                    contigs.close()
                    genome = open(outputDir + 'genome.txt', 'w')
                    contig = open(outputDir + 'unigene.sam', 'r')
                    for line in contig:
                        line = line.strip()
                        if "@SQ" in line:
                            line = line.replace(':', '\t').split('\t')
                            if line[2] in name_list:
                                genome.write(line[2] + '\t' + line[4] + '\n')
                    genome.close()
                    support = open(outputDir + 'unigene_mp', 'r')
                    coverage = open(outputDir + 'coverage.bed', 'w')
                    for line in support:
                        line = line.strip()
                        line = line.split('\t')
                        if line[2] in name_list:
                            coverage.write(
                                line[2] + '\t' + line[6] + '\t' + line[7] + '\n')
                    coverage.close()
                    BTcov = [
                        'bedtools',
                        'genomecov',
                        '-bg',
                        '-i',
                        outputDir +
                        'coverage.bed',
                        '-g',
                        outputDir +
                        'genome.txt']
                    BTcov_call = subprocess.check_output(BTcov)
                    bedFile = BTcov_call.splitlines()
                    bed = open(outputDir + 'assembly.cov', 'w')
                    for el in bedFile:
                        bed.write(el + '\n')
                    bed.close()
                    bed = open(outputDir + 'assembly.cov', 'r')
                    rangeDict = {}
                    for line in bed:
                        listStartEnd = []
                        info = line.split('\t')
                        # print info[0]
                        if info[0] in unitigDict and int(info[3]) > int(
                                int(unitigDict[info[0]]) / 20) and info[0] in rangeDict:
                            rangeDict[info[0]].append(info[2])
                        elif info[0] in unitigDict and int(info[3]) > int(int(unitigDict[info[0]]) / 20):
                            listStartEnd.append(info[1])
                            listStartEnd.append(info[2])
                            rangeDict[info[0]] = listStartEnd
                            # print rangeDict[info[0]]
                    finalDict = {}
                    for el in name_list:
                        listFinal = []
                        if el in rangeDict:
                            start = rangeDict[el][0]
                            end = int(rangeDict[el][-1])
                            listFinal.append(start)
                            listFinal.append(end)
                            finalDict[el] = listFinal
                    sequence = open(fasta_seq_new, 'w')
                    for record in outFasta:
                        name_seq = record.id
                        if name_seq in finalDict:
                            location = finalDict[name_seq]
                            nucleotide = record.seq
                            newNucl = nucleotide[int(finalDict[name_seq][0]):int(
                                finalDict[name_seq][1])]
                            sequence.write(
                                '>' + str(name_seq) + '\n' + str(newNucl) + '\n')
                        else:
                            sequence.write('>' + str(record.id) +
                                           '\n' + str(record.seq) + '\n')
        contigs.close()
    return


################
####MAIN###
################

def main():

    print '\n\n\n###############\n###FINISHED###\n###############\n\n'


if __name__ == '__main__':
    main()
