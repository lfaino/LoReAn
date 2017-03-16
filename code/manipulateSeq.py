from Bio import pairwise2
from Bio import SeqIO
import sys
from sys import argv
import numpy as np
import os
from Bio.Seq import reverse_complement
from multiprocessing import Pool
import itertools


def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return align_call(*a_b)

def filterLongReads(fastqFilename, min_length, max_length, wd, a):
    '''Filters out reads longer than length provided'''
    finalSeq= []
    if a:
        outFilename = wd + fastqFilename + '.longreads.filtered.fasta'
    else:
        outFilename = fastqFilename + '.longreads.filtered.fasta'
    
    if fastqFilename.endswith('fastq') or fastqFilename.endswith('fq'):
        fastqFile = open(fastqFilename, 'r')
        fastq = SeqIO.parse(fastqFile, 'fastq')
        if os.path.isfile(outFilename):
            print(('Filtered FASTQ existed already: ' +
                outFilename + ' --- skipping\n'))
            return outFilename, 0
        outFile = open(outFilename, 'w')
        filter_count = 0
        for record in fastq:
            if len(str(record.seq)) > int(min_length) and len(
                    str(record.seq)) < int(max_length):
                record.id = str(filter_count)
                finalSeq.append(record)
                filter_count += 1
    elif  fastqFilename.endswith('fasta') or fastqFilename.endswith('fa'):
        fastqFile = open(fastqFilename, 'r')
        fastq = SeqIO.parse(fastqFile, 'fasta')
        if os.path.isfile(outFilename):
            print(('Filtered FASTA existed already: ' +
                outFilename + ' --- skipping\n'))
            return outFilename, 0
        outFile = open(outFilename, 'w')
        filter_count = 0
        for record in fastq:
            if len(str(record.seq)) > int(min_length) and len(
                    str(record.seq)) < int(max_length):
                record.id = str(filter_count)
                finalSeq.append(record)
                filter_count += 1
    
    SeqIO.write(finalSeq, outFilename, "fasta")
    fastqFile.close()


    return (outFilename, filter_count)

def findOrientation(fastqFilename, min_length, max_length, wd, fastaAdapt, threads):
    '''Filters out reads longer than length provided'''
    outFilename = wd + fastqFilename + '.longreads.filtered.oriented.fasta'
    seqDict = {}
    scoreDict = {}
    listScore=[]
    listSeqGood = []
    listAdapter = []
    finalSeq = []
    listSeqAdap = []
    finalDNA = []
    pool = Pool(int(threads))
    for adpt in SeqIO.parse(fastaAdapt, "fasta"):
        listAdapter.append(adpt.id)
        listSeqAdap.append(adpt)
    if fastqFilename.endswith('fasta') or fastqFilename.endswith('fa'):
        fastqFile = open(fastqFilename, 'r')
        fastq = SeqIO.parse(fastqFile, 'fasta')
        if os.path.isfile(outFilename):
            print(('Filtered FASTA existed already: ' +
                outFilename + ' --- skipping\n'))
            return outFilename, 0
        outFile = open(outFilename, 'w')
        filter_count = 0
        for record in fastq:
            if len(str(record.seq)) > int(min_length) and len(str(record.seq)) < int(max_length):
                record.id = "seq" + str(filter_count)
                filter_count += 1
                listSeqGood.append(record)
        for adapter in listSeqAdap:
            results = pool.map(func_star, zip(listSeqGood, itertools.repeat(adapter)))
            finalDNA = finalDNA + results
        for record in finalDNA:
            if record[0].id in seqDict:
                otherString = [record[0], record[1]]
                seqDict[record[0].id] = seqDict[record[0].id] + otherString
                scoreDict[record[0].id] = scoreDict[record[0].id] + [record[2]]
            else:
                seqDict[record[0].id] = [record[0], record[1]]
                scoreDict[record[0].id] =  [record[2]]
        for key in scoreDict:
            if scoreDict[key][0] > scoreDict[key][1]:
                listScore.append(float(scoreDict[key][0]))
            elif scoreDict[key][0] < scoreDict[key][1]:
                listScore.append(float(scoreDict[key][1]))
        a = np.array(listScore)
        mean = np.mean(a)
        stderrS = np.std(a)
        valueOptimal = mean + stderrS
        count = 0
        revcom = 0
        same = 0
        lost = 0
        for key, score in list(scoreDict.items()):
            if scoreDict[key][0] > valueOptimal and scoreDict[key][1] > valueOptimal:
                lost += 1
                next
            else:
                if scoreDict[key][0] > scoreDict[key][1] and scoreDict[key][0] > valueOptimal and listAdapter[0] in seqDict[key][1]:
                    filter_count += 1
                    finalSeq.append(seqDict[key][0])
                    same += 1
                elif scoreDict[key][0] < scoreDict[key][1] and scoreDict[key][1] > valueOptimal and listAdapter[1] in seqDict[key][3]:
                    sequenze = reverse_complement(seqDict[key][2].seq)
                    seqDict[key][2].seq = sequenze
                    finalSeq.append(seqDict[key][2])
                    filter_count += 1
                    revcom += 1
    elif fastqFilename.endswith('fastq') or fastqFilename.endswith('fq'):
        fastqFile = open(fastqFilename, 'r')
        fastq = SeqIO.parse(fastqFile, 'fasta')
        if os.path.isfile(outFilename):
            print(('Filtered FASTA existed already: ' +
                outFilename + ' --- skipping\n'))
            return outFilename, 0
        outFile = open(outFilename, 'w')
        filter_count = 0
        for record in fastq:
            if len(str(record.seq)) > int(min_length) and len(str(record.seq)) < int(max_length):
                record.id = "seq" + str(filter_count)
                filter_count += 1
                listSeqGood.append(record)
        for adapter in listSeqAdap:
            results = pool.map(func_star, zip(listSeqGood, itertools.repeat(adapter)))
            finalDNA = finalDNA + results
        for record in finalDNA:
            if record[0].id in seqDict:
                otherString = [record[0], record[1]]
                seqDict[record[0].id] = seqDict[record[0].id] + otherString
                scoreDict[record[0].id] = scoreDict[record[0].id] + [record[2]]
            else:
                seqDict[record[0].id] = [record[0], record[1]]
                scoreDict[record[0].id] =  [record[2]]
        for key in scoreDict:
            if scoreDict[key][0] > scoreDict[key][1]:
                listScore.append(float(scoreDict[key][0]))
            elif scoreDict[key][0] < scoreDict[key][1]:
                listScore.append(float(scoreDict[key][1]))
        a = np.array(listScore)
        mean = np.mean(a)
        stderrS = np.std(a)
        valueOptimal = mean + stderrS
        count = 0
        revcom = 0
        same = 0
        lost = 0
        for key, score in list(scoreDict.items()):
            if scoreDict[key][0] > valueOptimal and scoreDict[key][1] > valueOptimal:
                lost += 1
                next
            else:
                if scoreDict[key][0] > scoreDict[key][1] and scoreDict[key][0] > valueOptimal and listAdapter[0] in seqDict[key][1]:
                    filter_count += 1
                    finalSeq.append(seqDict[key][0])
                    same += 1
                elif scoreDict[key][0] < scoreDict[key][1] and scoreDict[key][1] > valueOptimal and listAdapter[1] in seqDict[key][3]:
                    sequenze = reverse_complement(seqDict[key][2].seq)
                    seqDict[key][2].seq = sequenze
                    finalSeq.append(seqDict[key][2])
                    filter_count += 1
                    revcom += 1
    else:
        print("Can not recognize file type")
        
    outFile.close()
    SeqIO.write(finalSeq, outFilename, "fasta")
    fastqFile.close()
    return (outFilename, filter_count, lost)
            
    #return()
   
def align_call(record, adapter):
    #for record in listSeq:
    score = pairwise2.align.localms(record.seq, adapter.seq, 1, -10, -5, -4.9, one_alignment_only=1, score_only=1)
    record.name = score
    listValue = [record, adapter.id, score]
    return listValue


if __name__ == '__main__':
    
    fastaSeq = argv[1]
    min_length = argv[2]
    max_length = argv[3]
    wd = argv[4]
    fastaAdapt= argv[5]
    threads = argv[6]
    
    #filterLongReads(fastaSeq, min_length, max_length, wd, a = True)
    findOrientation(fastaSeq, min_length, max_length, wd, fastaAdapt, threads)
