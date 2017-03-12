from Bio import pairwise2
from Bio import SeqIO
import sys
from sys import argv
from Bio.pairwise2 import format_alignment
import numpy as np
import os
from Bio.Seq import reverse_complement


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
            print ('Filtered FASTQ existed already: ' +
                outFilename + ' --- skipping\n')
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
            print ('Filtered FASTA existed already: ' +
                outFilename + ' --- skipping\n')
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


    return outFilename, filter_count

def findOrientation(fastqFilename, min_length, max_length, wd, adapterSeq):
    '''Filters out reads longer than length provided'''
    outFilename = wd + fastqFilename + '.longreads.filtered.oriented.fasta'
    seqDict = {}
    scoreDict = {}
    listScore=[]
    listAdapter = []
    finalSeq = []

    for adpt in SeqIO.parse(fastaAdapt, "fasta"):
        listAdapter.append(adpt.id)
    if fastqFilename.endswith('fasta') or fastqFilename.endswith('fa'):
        fastqFile = open(fastqFilename, 'r')
        fastq = SeqIO.parse(fastqFile, 'fasta')
        if os.path.isfile(outFilename):
            print ('Filtered FASTA existed already: ' +
                outFilename + ' --- skipping\n')
            return outFilename, 0
        outFile = open(outFilename, 'w')
        filter_count = 0
        for record in fastq:
            if len(str(record.seq)) > int(min_length) and len(str(record.seq)) < int(max_length):
                record.id = "seq" + str(filter_count)
                filter_count += 1
                #print filter_count
                for adtp in SeqIO.parse(fastaAdapt, "fasta"):
                    a = pairwise2.align.localms(record.seq, adtp.seq, 1, -3, -1, -0.99, one_alignment_only=1, score_only=1)
                    if seqDict.has_key(record.id):
                        otherString = [record, adtp.id]
                        seqDict[record.id] = seqDict[record.id] + otherString
                        scoreDict[record.id] = scoreDict[record.id] + [a]
                    else:
                        seqDict[record.id] = [record, adtp.id]
                        scoreDict[record.id] = [a]
        for key in scoreDict:
            if scoreDict[key][0] > scoreDict[key][1]:
                listScore.append(float(scoreDict[key][0]))
            elif scoreDict[key][0] < scoreDict[key][1]:
                listScore.append(float(scoreDict[key][1]))
        a = np.array(listScore)
        mean = np.mean(a)
        stderrS = np.std(a)
        valueOptimal = mean + stderrS
        count =0
        for key, score in scoreDict.items():
            if scoreDict[key][0] > valueOptimal and scoreDict[key][1] > valueOptimal:
                next
            else:
                if scoreDict[key][0] > scoreDict[key][1] and scoreDict[key][0] > valueOptimal:
                    if listAdapter[0] in seqDict[key][1]:
                        filter_count += 1
                        finalSeq.append(seqDict[key][0])
                    elif listAdapter[1] in seqDict[key][1]:
                        sequenze = reverse_complement(seqDict[key][0].seq)
                        seqDict[key][0].seq = sequenze
                        finalSeq.append(seqDict[key][0])
                        filter_count += 1
                elif scoreDict[key][0] < scoreDict[key][1] and scoreDict[key][1] > valueOptimal:
                    if listAdapter[1] in seqDict[key][3]:
                        sequenze = reverse_complement(seqDict[key][0].seq)
                        seqDict[key][0].seq = sequenze
                        finalSeq.append(seqDict[key][0])
                        #finalSeq.append(seqDict[key][0].reverse_complement())
                        filter_count += 1
                    elif listAdapter[1] in seqDict[key][3]:
                        finalSeq.append(seqDict[key][0])
                        filter_count += 1

    elif fastqFilename.endswith('fastq') or fastqFilename.endswith('fq'):
        fastqFile = open(fastqFilename, 'r')
        fastq = SeqIO.parse(fastqFile, 'fastq')
        if os.path.isfile(outFilename):
            print ('Filtered FASTQ existed already: ' +
                outFilename + ' --- skipping\n')
            return outFilename, 0
        outFile = open(outFilename, 'w')
        filter_count = 0
        for record in fastq:
            if len(str(record.seq)) > int(min_length) and len(str(record.seq)) < int(max_length):
                record.id = "seq" + str(filter_count)
                filter_count += 1
                for adtp in SeqIO.parse(fastaAdapt, "fasta"):
                    a = pairwise2.align.localms(record.seq, adtp.seq, 1, -3, -1, -0.99, one_alignment_only=1, score_only=1)
                    if seqDict.has_key(record.id):
                        otherString = [record, adtp.id]
                        seqDict[record.id] = seqDict[record.id] + otherString
                        scoreDict[record.id] = scoreDict[record.id] + [a]
                    else:
                        seqDict[record.id] = [record, adtp.id]
                        scoreDict[record.id] = [a]
        for key in scoreDict:
            if scoreDict[key][0] > scoreDict[key][1]:
                listScore.append(float(scoreDict[key][0]))
            elif scoreDict[key][0] < scoreDict[key][1]:
                listScore.append(float(scoreDict[key][1]))
        a = np.array(listScore)
        mean = np.mean(a)
        stderrS = np.std(a)
        valueOptimal = mean + stderrS
        for key, score in scoreDict.items():
            if scoreDict[key][0] > valueOptimal and scoreDict[key][1] > valueOptimal:
                next
            else:
                if scoreDict[key][0] > scoreDict[key][1] and scoreDict[key][0] > valueOptimal:
                    if listAdapter[0] in seqDict[key][1]:
                        filter_count += 1
                        finalSeq.append(seqDict[key][0])
                    elif listAdapter[1] in seqDict[key][1]:
                        sequenze = reverse_complement(seqDict[key][0].seq)
                        seqDict[key][0].seq = sequenze
                        finalSeq.append(seqDict[key][0])
                        filter_count += 1
                elif scoreDict[key][0] < scoreDict[key][1] and scoreDict[key][1] > valueOptimal:
                    if listAdapter[1] in seqDict[key][3]:
                        sequenze = reverse_complement(seqDict[key][0].seq)
                        seqDict[key][0].seq = sequenze
                        finalSeq.append(seqDict[key][0])
                        #finalSeq.append(seqDict[key][0].reverse_complement())
                        filter_count += 1
                    elif listAdapter[1] in seqDict[key][3]:
                        finalSeq.append(seqDict[key][0])
                        filter_count += 1
    else:
        print "Can not recognize file type"
        
    
    outFile.close()
    SeqIO.write(finalSeq, outFilename, "fasta")
    fastqFile.close()

            
    return()

if __name__ == '__main__':
    
    fastaSeq = argv[1]
    min_length = argv[2]
    max_length = argv[3]
    wd = argv[4]
    fastaAdapt= argv[5]
    
    filterLongReads(fastaSeq, min_length, max_length, wd, a = True)
    findOrientation(fastaSeq, min_length, max_length, wd, fastaAdapt)
