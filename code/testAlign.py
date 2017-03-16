from Bio import pairwise2
from Bio import SeqIO
import sys
from sys import argv
from Bio.pairwise2 import format_alignment
import numpy as np
import os
from Bio.Seq import reverse_complement
from Queue import Queue
from threading import Thread, Lock




def align_multi(threads, adapter, fasta_list):
    '''handles the assembly process and parsing in a multithreaded way'''

    align_queue = Queue()
    counter = 1
    for record in fasta_list:
        align_queue.put(record)
    global length_cluster
    length_cluster = str(align_queue.qsize())

    for i in range(int(threads)):
        t[i] = Thread(target=alignParse, args=(align_queue, adapter))
        t.daemon = True
        t[i].start()
    align_queue.join()
    #parseAugustus(wd, single_fasta_list)


def alignParse(queue, adapter):
    '''to join the assembly and the parsing process'''
    while True:
        try:
            record = queue.get()
        except:
            break
        listValue  = align_call(record, adapter)
        queue.task_done()
        global count_sequences
        count_sequences += 1
        global length_cluster
    return listValue


def func_star(a_b):
    return augustusParse(*a_b)


def align_call(record, adapter):
    score = pairwise2.align.localms(record.seq, adapter.seq, 1, -3, -1, -0.99, one_alignment_only=1, score_only=1)
    listValue = [record, adtp.id, score]
    return listValue

if __name__ == '__main__':
    
    fastaSeq = argv[1]
    min_length = argv[2]
    max_length = argv[3]
    wd = argv[4]
    fastaAdapt= argv[5]
    
    filterLongReads(fastaSeq, min_length, max_length, wd, a = True)
    findOrientation(fastaSeq, min_length, max_length, wd, fastaAdapt)
