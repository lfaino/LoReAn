#!/usr/bin/env python



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


def gffread(gff3File, reference, working_dir):
    '''Runs gffread on a gff3 file to produce fasta files with the 
    matching records'''
    out_name = working_dir + 'getFasta.fasta'
    args = ['bedtools','getfasta', '-fi', reference, '-bed', gff3File, '-fo', out_name, '-name', '-split']
      
    if os.path.isfile(out_name): 
        print ('gff3read file existed already: ' + out_name + ' --- skipping\n')
        return out_name
    
    
    try:
        subprocess.check_call(args)
        #print '>gff3read worked. Output is: ' + out_name +'\n'
    except:
        #print 'gff3read did not work properly\n'
        raise NameError('')
    return out_name

                    
def cluster_pipeline(gff3File, mergeDistance, strand ,wd): 
    '''here the cluseter of sequence from the same locus are prepared'''    
    tmpFile = wd + 'clusters.tmp.bed'
    BTsort1 = ['bedtools', 'sort', '-i', gff3File]
    dist = '-' + str(mergeDistance)
    if strand:
        BTmerge1 = ['bedtools', 'merge','-s', '-i', gff3File + '.sorted.bed', '-d', str(dist), '-c', '4,4', '-o', 'count,distinct']
        
        print "\n\t###CLUSTERING IN STRANDED MODE###\n"
    else:
        BTmerge1 = ['bedtools', 'merge', '-i', gff3File + '.sorted.bed', '-d', str(dist), '-c', '4,4', '-o', 'count,distinct']
        print "\n\t###CLUSTERING IN NON-STRANDED MODE###\n"
    BTsort2 = ['bedtools', 'sort', '-i', gff3File + '.sorted.merged.bed']
    
    #Sort the GFF3 file
    BTsort1_call = subprocess.Popen(BTsort1, stdout = file(gff3File + ".sorted.bed", "w"))
    BTsort1_call.communicate()
    
    #Merge the BED entries, count number of reads on each merged entry    
    BTmerge1_call = subprocess.Popen(BTmerge1,stdout = file(gff3File + ".sorted.merged.bed", "w"))
    BTmerge1_call.communicate()
    #NSort it again and returns
    BTsort2_call = subprocess.check_output(BTsort2)        
    final_output = BTsort2_call.splitlines()
    return final_output
    
    
    
    
def fasta2Dict(fastaFilename):
    '''Prepare a dictionary of all the sequences that is used together with 
    the fasta file to make single fasta files for the assembly'''
    fastaFile = open(fastaFilename, 'r')
    fastaDict2 = SeqIO.to_dict(SeqIO.parse(fastaFile, 'fasta'))
    fastaDict = {}
    for key, seq2 in fastaDict2.items():
        seq = str(seq2.seq).replace("N", "")
        fastaDict[key] = seq
    fastaFile.close()
    return fastaDict
    
    
    
    
    
def write_fastas(count, bedline, fastaDict, min_length, min_evidence, max_evidence, wd):
    '''From the output list of the pipeline, recovers the ID and goes back to the 
    fasta file to retrieve the sequence'''
    
    line = bedline.split('\t')
    if len(line) == 6: 
        chrm, start, end, strand, number, idents = (line[0], line[1], line[2], line[3], line[4], line[5])
    elif len(line) == 5:
        chrm, start, end, number, idents = (line[0], line[1], line[2], line[3], line[4])

    ids_short = []
    if len(idents.split(',')) > int(min_evidence) and len(idents.split(',')) < int(max_evidence):
        for element in re.split(',|;', idents):
            #if 'ID=' in element:
            #We keep read.mrna and evm.model records                
            ids_short.append(element)
    else:
        return False
    
    unique_ids = list(set(ids_short))
    clusterFilename = wd + '_'.join([chrm, start, end]) + '_' + str(count) + '.fasta'
    clusterFile = open(clusterFilename, 'w')
    read_count = 1
    
    
    for iden in unique_ids:
        if fastaDict.has_key(iden):
            if len(str(fastaDict[iden])) > int(min_length):
                if len(iden) < 40:
                    clusterFile.write('>'+iden+'\n'+str(fastaDict[iden])+'\n')
                else:
                    clusterFile.write('>'+str(read_count)+'\n'+str(fastaDict[iden])+'\n')
                    read_count += 1
                del fastaDict[iden]
            else:
                del fastaDict[iden]
          
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
 
  
def iAssembler(fastaFile, overlapLength, percent_identity, overhang, wd):    
    '''Call iAssembler to assemble every cluster in fasta_list'''
    newFasta = fastaFile.split('/')[-1]
    args = ['iAssembler.pl', '-i', newFasta, '-h', str(overlapLength), 
            '-p', str(percent_identity), '-o', newFasta+'_output' , '2> ']    
    outputDir = wd + newFasta +'_output/' #whole path
    try:
        subprocess.call(args, cwd = wd, stdout = subprocess.PIPE,stderr = subprocess.PIPE)
        ##TRY TO CHANGE TO SUBPROCESSPOPEN
    except:
        return False
    return outputDir
  

def assembleParse (queue, overlapLength, percent_identity, overhang, wd):
    '''to join the assembly and the parsing process'''
    while True:    
        try:
            clusterFasta = queue.get()     
        except:
            break
        outputDir = iAssembler(clusterFasta, overlapLength, percent_identity, overhang, wd)
        queue.task_done()   
        global count_sequences
        count_sequences += 1
        global length_cluster
    return
def func_star(a_b):
    
    return assembleParse(*a_b)        

def assembly(clusterList, fastaDict, min_evidence, max_evidence, overlapLength,
             percent_identity, overhang, threads, wd):
    '''handles the assembly process and parsing in a multithreaded way'''
    assembly_queue = Queue()
    cluster_counter = 1
    for record in clusterList:
        #Write fasta for each cluster
        clusterFasta = write_fastas(cluster_counter, record, fastaDict, overlapLength, min_evidence, max_evidence, wd) 
        cluster_counter += 1
        if clusterFasta:
            #If it is a cluster to assemble put it on queue
            assembly_queue.put(clusterFasta)
        else: 
            continue
    global length_cluster
    length_cluster = str(assembly_queue.qsize())
    for i in range(int(threads)):
        t = Thread(target=assembleParse, args=(assembly_queue, overlapLength, percent_identity, overhang, wd))
        t.daemon = True
        t.start()    
    assembly_queue.join()
    

###############
###MAIN###
###############

def main():
    '''Main body of the function'''
    wd = './consensus/'   
    check_create_dir(wd)
    #wd = os.path.abspath(wd) + '/'
        
    gff3FileName = os.path.abspath(sys.argv[1])
    reference = os.path.abspath(sys.argv[2])
    fastaFilename = gffread(gff3FileName, reference, wd)
    gffreadDict = fasta2Dict(fastaFilename)
    #print gffreadDict
    cluster_list = cluster_pipeline(gff3FileName, 200, False ,wd)    
     
    
    assembly(cluster_list, gffreadDict, 3, 20, 200, 97, 2, 20, wd)
#    outputFilename = catAssembled(wd)
            
#    mergedFastaFilename = wd+'assembly.wEVM.fasta'
#    addEVM(fastaFilename, outputFilename, unitigs ,mergedFastaFilename)
    
    print '\n\n\n###############\n###FINISHED###\n###############\n\n'

    
    
    


if __name__ == '__main__':
    main()
