from Bio import pairwise2
from Bio import SeqIO
import sys


def findOrientation(fastaSeq, fastaAdapt):
        
    for t in SeqIO.parse(fastaSeq, "fasta"):
        for q in SeqIO.parse(fastaAdapt, "fasta"):
            a = pairwise2.align.localms(t.seq, q.seq, 1, -3, -7, -2, one_alignment_only=1, score_only=1)
            print t.id, "\t", q.id, "\t", a
            
    return()
if __name__ == '__main__':
    
    
    fastaSeq = argv[1]
    fastaAdapt= argv[2]
    #wd = argv[3]
    #prefix = argv[4]
    out1 = findOrientation(fastaSeq, fastaAdapt)
    #out2 = appendID(out1)
    #out3 = removeOverlap(out2)
    #out4 = removeDiscrepancy(out3, evmgff)
    #newNames(out4)
    #removeOverlap(out2)
    #genename(out2, prefix)
