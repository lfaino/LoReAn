#!/usr/bin/python



#######################################
######DOCUMENT THIS FILE################
#######################################





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


count_sequences = 0
length_cluster = 0

###############
###FUNCTIONS###
###############



       
def parseOnly (threshold_float, unitigs ,tmp_wd):
    '''to join the assembly and the parsing process'''
    evm_list = []
    wd = tmp_wd +'consensus/tmp/'
    for root , dirs, _ in os.walk(wd):
    	for direc in dirs:
            if 'output' in direc:
                outputDir = wd + direc + '/'
                if outputDir:
                    parse_contigs(outputDir, threshold_float, unitigs)
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
        count_sequences+= 1
        global length_cluster
    return evm_list


def parse_contigs(outputAssembly, threshold_float, unitigs):
    '''Parses the output from iAssembler, to a single FASTA file'''
    #PARSE CONTIG MEMBERS FILE
    fname = outputAssembly+'contig_member'
    fname_exists = os.path.isfile(fname)
    if fname_exists:
        ###part all the files in the tmp assembly folder
        contigInfo = open(outputAssembly+'contig_member', 'r')
        contigs = {}
        total_reads = 0    
        for line in contigInfo:
            line = line.strip() 
            line = line.split('\t')
            ###count the reads for each assembled Unitig
            read_number = len(line) - 1
            for element in line:
                if 'evm' in element:
                    contigs[line[0]] = [read_number, element]      
            if not contigs.has_key(line[0]):
                contigs[line[0]] = [read_number]
            ### sum all the reads whitin one cluster
            total_reads += read_number
        contigInfo.close()
        ### calculate the number of reads considering the threshold
        threshold = total_reads*float(threshold_float)
        global count_sequences
        real_contigs = {}
        count_unitigs = 1
        for key, element in contigs.items():
            ### to retrieve only supported assembly
            if unitigs:
                if len(element) == 2:
                    if element[0] >= threshold and 'evm' in element[1]:
                        real_contigs[key] = element[1]+' '+str(element[0])+'_above_threshold_'+str(threshold)+' loc_'+str(count_sequences)
                elif len(element) == 1:
                    if element[0] >= threshold:
                        real_contigs[key] = 'Unitig'+str(count_sequences) + '_' + str(count_unitigs)+' '+str(element[0])+'_above_threshold_'+str(threshold)+' loc_'+str(count_sequences)
                        count_unitigs += 1
            ### to retrieve all the evm even if not supported
            else:
                if len(element) == 2:
                    if element[0] >= threshold and 'evm' in element[1]:
                        real_contigs[key] = element[1]+' '+str(element[0])+'_above_threshold_'+str(threshold)+' loc_'+str(count_sequences)
                    elif element[0] < threshold and 'evm' in element[1]:
                        real_contigs[key] = element[1]+' '+str(element[0])+'_below_threshold_'+str(threshold)+' loc_'+str(count_sequences)       
                elif len(element) == 1:
                    if element[0] >= threshold:
                        real_contigs[key] = 'Unitig'+str(count_sequences) + '_' + str(count_unitigs)+' '+str(element[0])+'_above_threshold_'+str(threshold)+' loc_'+str(count_sequences)
                        count_unitigs += 1
        ### writes the outputs
        fileAssembly = outputAssembly+'unigene_seq.new.fasta'
        contigSeq = open(fileAssembly, 'r')   
        contigDict = SeqIO.to_dict(SeqIO.parse(contigSeq, 'fasta'))
        outputFilename = outputAssembly[:-1] + '_assembled.fasta'
        outputFile = open(outputFilename, 'w')
        for iden, write_iden in real_contigs.items():
            if contigDict.has_key(iden):            
                outputFile.write('>'+write_iden+'\n'+str(contigDict[iden].seq)+'\n')
        contigSeq.close()  
        outputFile.close()
    return

def catAssembled(wd):
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


def addEVM(wholeFastaName, outputFilename, unitigs, evm_nosupport, outputMergedFastaName):
    '''Adds the EVM records that are not present in the final contig evidence'''
    wholeFasta = open(wholeFastaName, 'r')
    outFastaFile = open(outputFilename, 'r')    
    outputMerged = open(outputMergedFastaName, 'w')
    wholeDict = SeqIO.to_dict(SeqIO.parse(wholeFasta, 'fasta'))
    count = 0
    if unitigs:
        outFasta = SeqIO.parse(outputMergedFastaName, 'fasta')
        for record in outFasta:
            count += 1
            ident = '>Gene'+str(count)
            outputMerged.write(ident+'\n'+str(record.seq)+'\n')
            exons = description.split(',')
            if len(exons) == 1 and  "evm" not in record.id:
                SeqIO.write(record, s_file , "fasta")
            else:
                SeqIO.write(record, m_file , "fasta")
        newDict = {}
        for key, values in wholeDict.items():
            if key in evm_nosupport:
                continue
            else:
                newDict[key] = values
        dictOut = {}
        outFasta = SeqIO.parse(outFastaFile, 'fasta') 
        for record in outFasta:
#            count = 0
            if dictOut.has_key(record.id):
                dictOut[str(record.id)+'_'+str(count)] = str(record.seq)
                count += 1                    
            else:
                dictOut[record.id] = str(record.seq)
        seq_count = 1
        for key in newDict.keys():        
            if 'evm' in key and not dictOut.has_key(key):   
                ident = '>Gene'+str(count)+'_'+key
                outputMerged.write(ident+'\n'+str(newDict[key].seq)+'\n')
                count += 1
        for key, element in dictOut.items():
            ident = '>Gene'+str(count)+'_'+key
            outputMerged.write(ident+'\n'+str(element)+'\n')
            count += 1
    else:
        dictOut = {}
        outFasta = SeqIO.parse(outFastaFile, 'fasta') 
        
        for record in outFasta:
#            count = 0
            if dictOut.has_key(record.id):
                dictOut[str(record.id)+'_'+str(count)] = str(record.seq)
                count += 1                    
            else:
                dictOut[record.id] = str(record.seq)
        seq_count = 1
        for key in wholeDict.keys():        
            if 'evm' in key and not dictOut.has_key(key):   
                ident = '>Gene'+str(count)+'_'+key
                outputMerged.write(ident+'\n'+str(wholeDict[key].seq)+'\n')
                count += 1
        for key, element in dictOut.items():
            ident = '>Gene'+str(count)+'_'+key
            outputMerged.write(ident+'\n'+str(element)+'\n')
            count += 1
            
    wholeFasta.close()
    outFasta.close()
    outputMerged.close()
    
def getEVMnoUnitig(target_wd):
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


###############
###MAIN###
###############

def main():
    '''Main body of the function'''
#    wd = 'consensus/tmp/'   
#    check_create_dir(wd)
#    wd = os.path.abspath(wd) + '/'
#    print wd
        
#    gff3FileName = os.path.abspath(sys.argv[1])
#    reference = os.path.abspath(sys.argv[1])
#    fastaFilename = wd + 'gffread.fasta'
#    threshold_float = 0.3
#    unitigs = False
        
    evm_list = parseOnly(threshold_float, unitigs, wd)
    outputFilename = catAssembled(wd)       
    mergedFastaFilename = wd+'assembly.wEVM.fasta'
    addEVM(fastaFilename, outputFilename, unitigs , evm_list, mergedFastaFilename)
    
    print '\n\n\n###############\n###FINISHED###\n###############\n\n'

    
if __name__ == '__main__':
    main()

