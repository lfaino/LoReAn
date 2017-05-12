#!/usr/bin/env python3
from Bio import pairwise2
from Bio import SeqIO
import sys
from sys import argv
import numpy as np
import os
from Bio.Seq import reverse_complement
import subprocess
import itertools
from Bio.SeqUtils import GC
import os.path as op
import ctypes as ct
import timeit as ti
import math
from Bio import SeqIO
import ssw_lib
from Bio.Seq import reverse_complement
import io
from multiprocessing import Pool

def to_int(seq, lEle, dEle2Int):
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i,ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n
    return num

def align_one(ssw, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, int(nMaskLen))
    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)
    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)

def align_call(elem):
    record = elem[0]
    adapter = elem[1]
    lEle = []
    dRc = {} 
    dEle2Int = {}
    dInt2Ele = {}
    nMatch = 2
    nMismatch = 1
    nOpen = 1
    nExt = -2
    nFlag = 0
    #if not args.sMatrix:
    lEle = ['A', 'C', 'G', 'T', 'N']
    for i,ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    nEleNum = len(lEle)
    lScore = [0 for i in range(nEleNum**2)]
    for i in range(nEleNum-1):
        for j in range(nEleNum-1):
            if lEle[i] == lEle[j]:
                lScore[i*nEleNum+j] = nMatch
            else:
                lScore[i*nEleNum+j] = -nMismatch
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore
    
    ssw = ssw_lib.CSsw("./")
    sQSeq = record.seq
    sQId = record.id
    if len(sQSeq) > 30:
        nMaskLen = len(sQSeq) / 2
    else:
        nMaskLen = 15
    outputAlign = []
    qNum = to_int(sQSeq, lEle, dEle2Int)
    qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)
    sQRcSeq = reverse_complement(sQSeq)
    qRcNum = to_int(sQRcSeq, lEle, dEle2Int)
    qRcProfile = ssw.ssw_init(qRcNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)
    sRSeq = adapter.seq
    sRId = adapter.id
    rNum = to_int(sRSeq, lEle, dEle2Int)
    res = align_one(ssw, qProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)
    resRc = None
    resRc = align_one(ssw, qRcProfile, rNum, len(sRSeq), nOpen, nExt, nFlag, nMaskLen)
    strand = 0
    if res[0] == resRc[0]:
        next
    if res[0] > resRc[0]:
        resPrint = res
        strand = 0
        outputAlign = [sRId , sQId, strand, resPrint[0]]
    elif res[0] < resRc[0]:
        resPrint = resRc
        strand = 1
        outputAlign = [sRId , sQId, strand, resPrint[0]]
    ssw.init_destroy(qProfile)
    ssw.init_destroy(qRcProfile)
    return outputAlign

def filterLongReads(fastqFilename, min_length, max_length, wd, adapter, threads, a):
    '''Filters out reads longer than length provided'''
    finalSeq= []
    seqDict = {}
    firstDictScore = {}
    firstDictSeq = {}
    scoreDict = {}
    listScore=[]
    listSeqGood = []
    listAdapter = []
    finalSeq = []
    listSeqAdap = []
    record_dict = {}
    maxScore = 0
    if a and not adapter:
        outFilename = wd + fastqFilename + '.longreads.filtered.fasta'
    elif a and adapter:
        outFilename = wd + fastqFilename + '.longreads.filtered.oriented.fasta'
    else:
        outFilename = fastqFilename + '.longreads.filtered.fasta'
    filter_count = 0
    if os.path.isfile(outFilename):
            print(('Filtered FASTQ existed already: ' +
                outFilename + ' --- skipping\n'))
            return outFilename, 0
    if fastqFilename.endswith('fastq') or fastqFilename.endswith('fq'):
        for record in SeqIO.parse(fastqFilename, "fastq"):
            if len(str(record.seq)) > int(min_length) and len(str(record.seq)) < int(max_length):
                record.description= ""
                record.name = ""
                record.id = str(filter_count)
                filter_count += 1
                record_dict[record.id] = record
    elif fastqFilename.endswith('fasta') or fastqFilename.endswith('fa'):
        for record in SeqIO.parse(fastqFilename, "fasta"):
            if len(str(record.seq)) > int(min_length) and len(str(record.seq)) < int(max_length):
                record.description= ""
                record.name = ""
                record.id = str(filter_count)
                filter_count += 1
                record_dict[record.id] = record
    if adapter:
        for adpt in SeqIO.parse(adapter, "fasta"):
            listAdapter.append(adpt.id)
            listSeqAdap.append(adpt)
    outFile = open(outFilename, 'w')
    
    if len(listAdapter) == 1:
        for key in record_dict:
            if len(str(record_dict[key].seq)) > int(min_length) and len(str(record_dict[key].seq)) < int(max_length):
                for adpter in listSeqAdap:
                    alingRes = align_call(record_dict[key], adpter)
                    if len(alingRes) == 0:
                        next
                    else:
                        seqDict[record_dict[key].id] = [record_dict[key], alingRes[2]]
                        scoreDict[record_dict[key].id] =  alingRes[3]
        for key in scoreDict:
            #for score in scoreDict[key]:
            if scoreDict[key] > maxScore:
                maxScore = scoreDict[key]
        valueOptimal = maxScore - (maxScore/10)
        for key in scoreDict:
            if scoreDict[key]  > valueOptimal and seqDict[key][1] == 0:
                filter_count += 1
                finalSeq.append(seqDict[key][0])
            elif scoreDict[key]  > valueOptimal and seqDict[key][1] == 1:
                filter_count += 1
                sequenze = reverse_complement(seqDict[key][0].seq)
                seqDict[key][0].seq = sequenze
                finalSeq.append(seqDict[key][0])

    elif len(listAdapter) == 2:
        listCommand = []
        for key in record_dict:
            for adpter in listSeqAdap:
                listCommand.append([record_dict[key], adpter])
        with Pool(int(threads)) as p:
            alignResul = p.map(align_call, listCommand)
        for alingRes in alignResul:
            if len(alingRes) == 0:
                next
            elif alingRes[1] in firstDictSeq:
                seqDict[alingRes[1]] =  firstDictSeq[alingRes[1]]  +  [ alingRes[2], alingRes[0]]
                scoreDict[alingRes[1]] = firstDictScore[alingRes[1]]  +  [ alingRes[0], alingRes[3]]
            else:
                firstDictSeq[alingRes[1]] = [record_dict[alingRes[1]], alingRes[2], alingRes[0]]
                firstDictScore[alingRes[1]] =  [ alingRes[0], alingRes[3]]
        maxScoreFirst = 0
        maxScoreSecond = 0
        for key in scoreDict:
            score = scoreDict[key]
            if score[0] == listAdapter[0] and score[1] > maxScoreFirst:
                maxScoreFirst = score[1]
            if score[2] == listAdapter[0] and score[3] > maxScoreFirst:
                maxScoreFirst = score[3]
            if score[0] == listAdapter[1] and score[1] > maxScoreSecond:
                maxScoreSecond = score[1]
            if score[2] == listAdapter[1] and score[3] > maxScoreSecond:
                maxScoreSecond = score[3]
        valueOptimalFirst = maxScoreFirst - (maxScoreFirst/20)
        valueOptimalSecond = maxScoreSecond - (maxScoreSecond/20)
        listReadsOverLimit = []
        for key in scoreDict:
            score = scoreDict[key]
            if (score[0] == listAdapter[0] and score[1] > valueOptimalFirst) and (score[2] == listAdapter[1] and score[3] > valueOptimalSecond):
                listReadsOverLimit.append(key)
            elif (score[2] == listAdapter[0] and score[3] > valueOptimalFirst) and (score[0] == listAdapter[1] and score[1] > valueOptimalSecond):
                listReadsOverLimit.append(key)

        for key in listReadsOverLimit:
            if seqDict[key][1] == 1 and seqDict[key][3] == 0:
                if seqDict[key][2] == listAdapter[0] and seqDict[key][4] == listAdapter[1]:
                    finalSeq.append(seqDict[key][0])
                elif seqDict[key][2] in listAdapter[1] and seqDict[key][4] in listAdapter[0]:
                    sequenze = reverse_complement(seqDict[key][0].seq)
                    seqDict[key][0].seq = sequenze
                    finalSeq.append(seqDict[key][0])
            elif seqDict[key][1] == 0 and seqDict[key][3] == 1:
                if seqDict[key][2] == listAdapter[0] and seqDict[key][4] == listAdapter[1]:
                    sequenze = reverse_complement(seqDict[key][0].seq)
                    seqDict[key][0].seq = sequenze
                    finalSeq.append(seqDict[key][0])
                elif seqDict[key][2] in listAdapter[1] and seqDict[key][4] in listAdapter[0]:
                    finalSeq.append(seqDict[key][0])
    elif len(listAdapter) == 0:
        for key in record_dict:
            if len(str(record_dict[key].seq)) > int(min_length) and len(str(record_dict[key].seq)) < int(max_length):
                finalSeq.append(record_dict[key])
    SeqIO.write(finalSeq, outFilename, "fasta")
    return (outFilename, filter_count)

def maskedgenome(wd, ref , gff3):
    if '/' in ref:
        out_name = wd + "/" +  ref.split('/')[-1] + '.masked.fasta'
    else:
        out_name = wd + "/" +  ref + '.masked.fasta'
    #out_name = gmap_ref
    outmerged = wd + "/" + gff3 + '.masked.gff3'
    outputmerge = open(outmerged, 'w')
    cat = subprocess.Popen(['cat', gff3], stdout = subprocess.PIPE)
    bedsort = subprocess.Popen(['bedtools', 'sort'], stdin = cat.stdout, stdout = subprocess.PIPE)
    bedmerge = subprocess.Popen(['bedtools', 'merge'], stdout = outputmerge, stdin = bedsort.stdout)
    out = bedmerge.communicate()
    outputmerge.close()
    maskfasta = subprocess.Popen(['bedtools', 'maskfasta', '-fi', ref , '-bed', outmerged, '-fo', out_name])
    maskfasta.communicate()
    return out_name

if __name__ == '__main__':

    fastqFilename = argv[1]
    min_length = argv[2]
    max_length = argv[3]
    wd = argv[4]
    adapter = argv[5] 
    threads = argv[6]
    a = True
    filterLongReads(fastqFilename, min_length, max_length, wd, adapter, threads, a)
