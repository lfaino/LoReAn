#! /usr/bin/env python

#
#    Helper script for aligning PacBio reads: cleans reads using reference sequence
#
#    Copyright (C) 2015  Michael Hamilton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import pysam
from argparse import ArgumentParser, ArgumentTypeError
import os, sys, numpy
from SpliceGrapher.formats.FastaLoader import FastaLoader                
from SpliceGrapher.shared.utils import *    
parser = ArgumentParser(description='Fix aligned reads using reference genome')                                                                  
parser.add_argument('-v', '--verbose', dest="verbose",
                    action='store_true', default=False,
                    help='Verbose mode')
parser.add_argument('-a', '--anchorIgnoreLim', dest="anchorAdjustLim",
                    action='store', type=int, default=0, 
                    help='Amount of anchor to ignore when correcting aligned reads.')

parser.add_argument('-f', '--fasta_output', dest="fasta_output",
                    action='store', type=str, default='fixed.fa', 
                    help='File name for writing fixed FASTA records')
parser.add_argument('-u', '--unaligned', dest="unaligned_output",
                    action='store', type=str, default='unaligned.txt', 
                    help='File name for writing unaligned reads')
parser.add_argument('-r', '--filtered', dest="filtered_output",
                    action='store', type=str, default='filtered.txt', 
                    help='File name for writing rejected FASTA records')
parser.add_argument('-j', '--junctions', dest="junctions_output",
                    action='store', type=str, default='junctions.txt', 
                    help='File name for writing junctions summary')

parser.add_argument('-s', '--sam_output', dest="sam_output",
                    action='store', type=str, default='fixed.sam', 
                    help='Name for writing fixed SAM records')

parser.add_argument('-t', '--thresh', dest="thresh",
                    action='store', type=int, default=-1,
                    help='Mapping quality threshold')

parser.add_argument('-e', '--edr', dest="edr",
                    action='store', type=float, default=0.10,
                    help='Edit distance ratio threshold')

parser.add_argument('reference', action='store', 
                    type=str, help='Reference sequence')
parser.add_argument('bam_input', action='store', 
                    type=str, help='BAM file input file to fix')
args = parser.parse_args()
MAPQTHRESH = args.thresh
MATCH    = 0   #M
INSERT   = 1   #I
DELETE   = 2   #D
GAP      = 3   #N
SOFT_CLIP= 4   #4
HARD_CLIP= 5   #H
PAD      = 6   #P
EQUAL    = 7   #=
DIFF     = 8   #X
AMATCHL  = 101 #match for keeping left anchors untouched
AMATCHR  = 102 #match for keeping right anchors untouched
AINSERT  = 11  #insert for keeping anchors untouched
ADELETE  = 12  #delete for keeping anchors untouched

loader = FastaLoader(args.reference, verbose=args.verbose)
bamfile = pysam.Samfile(args.bam_input , 'rb' )
outfile = pysam.Samfile( args.sam_output, 'wb', header=bamfile.header )
fastafile = open(args.fasta_output, 'w')
unalignedfile = open(args.unaligned_output, 'w')
filteredfile = open(args.filtered_output, 'w')
badSJs =  set()
filtered  = 0
unaligned = 0
accepted  = 0
mismatches = 0
deletions = [ ]
insertions = [ ]

mismatchesSJ = 0
deletionsSJ = [ ]
insertionsSJ = [ ]
totN = 0 
totSJ = 0
sSJ = 0
pstrand = 0
nstrand = 0
astrand = 0
spliced = 0
indicator = ProgressIndicator(100000, verbose=args.verbose)
canonicalSS = ['GTAG', 'GCAG']
SJdict = {}
revCDict = {'G':'C','C':'G','T':'A','A':'T','N':'N'}
def revComp(s):
    return ''.join([revCDict[x.upper()] for x in s])[::-1]
    
def junctionAdjust( cigarList, aLim ):
    """
    Adjust cigar so that we do not fix up to k nucleotides in anchors
    of a gap in the alignment.
     -------GT^AG-------
     |--k--|     |--k--| 
    """
    gapIdxs = [i for i in xrange(len(cigarList)) if cigarList[i][0] == 3]
    # fix left anchor
    for ridx,gidx in enumerate(gapIdxs):
        # do not iterate past a previous splice junction
        if ridx == 0:
            stop = -1
        else:
            stop = max(gapIdxs[ridx-1], -1)
        
        aRem = aLim
        iIdx = gidx - 1
        while iIdx > stop and aRem > 0:
            t,l = cigarList[iIdx]
            # MATCH
            # create a match record with amount to 
            if t == MATCH:
                if aRem < l:
                    keep = aRem
                    adj = l - aRem
                else:
                    keep = l
                    adj = 0
                cigarList[iIdx] = (AMATCHL, '%d,%d' % (adj, keep))
                aRem -= keep
            elif t == INSERT:
                cigarList[iIdx] = (AINSERT, l)
                aRem -= l
            elif t == DELETE:
                cigarList[iIdx] = (ADELETE, l)
            iIdx -= 1
    # fix right anchor
    for ridx,gidx in enumerate(gapIdxs):
        # do not iterate past a following splice junction
        if ridx == len(gapIdxs) - 1:
            stop = len(cigarList)
        else:
            stop = min(gapIdxs[ridx+1], len(cigarList))
        
        aRem = aLim
        iIdx = gidx + 1
        while iIdx < stop and aRem > 0:
            t,l = cigarList[iIdx]
            # MATCH
            if t == MATCH:
                if aRem < l:
                    keep = aRem
                    adj = l - aRem
                else:
                    adj = 0
                    keep = l
                cigarList[iIdx] = (AMATCHR, '%d,%d' % (adj, keep))
                aRem -= keep
            elif t == INSERT:
                cigarList[iIdx] = (AINSERT, l)
                aRem -= l
            elif t == DELETE:
                cigarList[iIdx] = (ADELETE, l)
            iIdx += 1

if args.verbose:
    sys.stderr.write('Processing SAM records\n')

for read in bamfile:
    ed = 0 
    rmismatches = 0
    rdeletions = [ ]
    rinsertions = [ ]

    rmismatchesSJ = 0
    rdeletionsSJ = [ ]
    rinsertionsSJ = [ ]

    ambSJ = False
    ambseq = False
    rtotSJ = 0
    rsSJ = 0
    rSJdict = {}
    indicator.update()
    rpos = 0
    pos = read.pos
    cleanSeq = ''
    cleanCigar = [ ]
    matchLen = 0
    # check if unmapped
    if read.is_unmapped:
        unalignedfile.write('>%s\n' % (read.qname))
        unaligned += 1
        continue
    # check if alignment exceeds edit distance ratio threshold
    tagD = dict(read.tags)
    strand = tagD.get('XS', None)
    qname = read.qname.split('@')[0]
    try:
        XL, XR = [int(x) for x in read.qname.split('@')[-1].split('_') ]
    except:
        XL, XR = 0,0
    if read.cigar[0][0] == SOFT_CLIP:
        XL += read.cigar[0][1]
    if read.cigar[-1][0] == SOFT_CLIP:
        XR += read.cigar[-1][1]
        
    junctions = [ ]
    cigarL = read.cigar[:] 

    if args.anchorAdjustLim > 0:
        junctionAdjust(cigarL, args.anchorAdjustLim)
    for i,rec in enumerate(cigarL):
        t,l = rec
        if t == MATCH:
            refseq =  loader.subsequence( bamfile.getrname(read.tid), pos, pos+l-1 )
            pos = pos+ l
            readseq  = read.query[rpos:(rpos+l)]
            for sidx in xrange(len(refseq)):
                if refseq[sidx] != readseq[sidx]: rmismatches += 1
            if 'N' in refseq:
                cleanSeq = ''.join([cleanSeq,readseq] )
                ambseq = True
            else:
                cleanSeq = ''.join([cleanSeq,refseq] )
            rpos = rpos + l

            matchLen += l
        # adjusted match left anchor
        elif t == AMATCHL:
            adj, keep = [int(x) for x in l.split(',')]
            # adjust
            if adj > 0:
                refseq =  loader.subsequence( bamfile.getrname(read.tid), pos, pos+adj-1 )
                pos = pos+ adj
                readseq  = read.query[rpos:(rpos+adj)]
                for sidx in xrange(len(refseq)):
                    if refseq[sidx] != readseq[sidx]: rmismatches += 1

                if 'N' in refseq:
                    cleanSeq = ''.join([cleanSeq,readseq] )
                    ambseq = True
                else:
                    cleanSeq = ''.join([cleanSeq,refseq] )
                rpos = rpos + adj
                matchLen += adj
            #keep
            readseq  = read.query[rpos:(rpos+keep)]
            refseq =  loader.subsequence( bamfile.getrname(read.tid), pos, pos+keep-1 )
            for sidx in xrange(len(refseq)):
                if refseq[sidx] != readseq[sidx]: 
                    rmismatchesSJ += 1
                    ed += 1
            cleanSeq = ''.join([cleanSeq,readseq] )
            pos += keep
            rpos += keep
            matchLen += keep
        # adjusted match right anchor
        elif t == AMATCHR:
            adj, keep = [int(x) for x in l.split(',')]
            #keep
            readseq  = read.query[rpos:(rpos+keep)]
            refseq =  loader.subsequence( bamfile.getrname(read.tid), pos, pos+keep-1 )
            cleanSeq = ''.join([cleanSeq,readseq] )
            pos += keep
            rpos += keep
            matchLen += keep
            for sidx in xrange(len(refseq)):
                if refseq[sidx] != readseq[sidx]: 
                    rmismatchesSJ += 1
                    ed += 1 

            # adjust
            if adj > 0:
                refseq =  loader.subsequence( bamfile.getrname(read.tid), pos, pos+adj-1 )
                pos = pos+ adj
                readseq  = read.query[rpos:(rpos+adj)]
                for sidx in xrange(len(refseq)):
                    if refseq[sidx] != readseq[sidx]: 
                        rmismatches += 1
                        ed += 1
                if 'N' in refseq:
                    cleanSeq = ''.join([cleanSeq,readseq] )
                    ambseq = True
                else:
                    cleanSeq = ''.join([cleanSeq,refseq] )
                rpos = rpos + adj
                matchLen += adj
        #insertion
        elif t == INSERT:
            rpos += l
            rinsertions.append(l)
            ed += l
        #adjusted insertion
        elif t == AINSERT:
            readseq  = read.query[rpos:(rpos+l)]
            cleanSeq = ''.join([cleanSeq,readseq] )
            rpos += l
            rinsertionsSJ.append(l)
            ed += l 
        #deletion
        elif t == DELETE:
            refseq = loader.subsequence( bamfile.getrname(read.tid), pos, pos+l-1 )
            cleanSeq = ''.join([cleanSeq, refseq])
            pos += l
            matchLen += l
            rdeletions.append(l)
            ed += l 
        #adjusted deletion
        elif t == ADELETE:
            pos += l
            rdeletionsSJ.append(l)
            ed += 1
        elif t == GAP:
            rtotSJ += 1
            if matchLen == 0:
                print read.qname
            cleanCigar.append( (0,matchLen) )
            matchLen = 0    
            cleanCigar.append( (t,l) )
            p5 = loader.subsequence( bamfile.getrname(read.tid), pos, pos+1 )
            p3 = loader.subsequence( bamfile.getrname(read.tid), pos+l-2, pos+l-1 )
            if 'N' in p5 or 'N' in p3: 
                ambSJ = True
            elif strand == '+':
                pair = '%s-%s' % (p5,p3)
                rsSJ += 1
                rSJdict[pair] = rSJdict.get(pair,0) + 1
                junctions.append(pair)
            elif strand == '-':
                pair = '%s-%s' % (revComp(p3),revComp(p5))
                rsSJ += 1
                rSJdict[pair] = rSJdict.get(pair,0) + 1
                junctions.append(pair)
            pos += l 
        else:
            if i > 0 and i+1 < len(read.cigar):
                print i, len(read.cigar)
                print t
                print read.cigar
    totN += read.qlen
    if ambSJ:
        filtered += 1
        filteredfile.write('%s_ambSJ\n' % (read.qname))
        continue
    if ambseq:
        filtered += 1
        filteredfile.write('%s_ambseq\n' % (read.qname))
        continue
        

    #ed = rmismatches + sum(rdeletions) + sum(rinsertions) + rmismatchesSJ + sum(rdeletionsSJ) + sum( rinsertionsSJ)

    if float(ed) / len(read.query) > args.edr:
        filtered += 1
        filteredfile.write('>%s_edr\n' % (read.qname))
        continue

    mismatches += rmismatches
    deletions.extend(rdeletions)
    insertions.extend(rinsertions)

    mismatchesSJ += rmismatchesSJ
    deletionsSJ.extend(rdeletionsSJ)
    insertionsSJ.extend(rinsertionsSJ)
    totSJ += rtotSJ
    if rtotSJ > 0:
        spliced += 1
    sSJ += rsSJ
    for p in rSJdict:
        SJdict[p] = SJdict.get(p,0) + rSJdict[p]
    cleanCigar.append( (0,matchLen) )
    if strand == '+': pstrand += 1
    elif strand == '-': nstrand += 1
    else: astrand += 1
    a = pysam.AlignedRead()
    a.qname = '%s@%d_%d' % (qname, XL, XR)
    a.seq = cleanSeq if len(cleanSeq) > 0 else read.seq
    a.flag = read.flag
    a.rname = read.rname
    a.pos = read.pos
    a.mapq = read.mapq
    a.cigar = cleanCigar if len(cleanCigar) > 0 else '*'
    a.mrnm = read.mrnm
    a.mpos = read.mpos
    a.isize = read.isize
    a.tags = read.tags
    if junctions:
        a.setTag( 'XD', ','.join(junctions))
    a.setTag( 'XL', XL )
    a.setTag( 'XR', XR )
    accepted += 1
    # no more errors
    if ed == 0:
        outfile.write(a)
    
    elif len(cleanSeq) > 0:
        fastafile.write('>%s\n%s\n' % (a.qname, cleanSeq))
    else:
        fastafile.write('>%s\n%s\n' % (a.qname, read.seq))
outfile.close()
indicator.finish()
if args.verbose:
    mismatchesPerc = float(mismatches)/totN * 100
    if totSJ == 0:
        mismatchesSJPerc = 0
        deletionsSJPerc = 0
        insertionsSJPerc =0
    else:
        mismatchesSJPerc = float(mismatchesSJ)/totSJ * 100
        deletionsSJPerc = float(len(deletionsSJ))/totSJ* 100
        insertionsSJPerc = float(len(insertionsSJ))/totSJ* 100
    deletionsPerc = float(len(deletions))/totN * 100
    if len(deletions) == 0:
        deletionsMean = 0
    else:    deletionsMean = float(numpy.mean(deletions))
    if len(deletionsSJ) == 0:
        deletionsSJMean = 0
    else:    deletionsSJMean = numpy.mean(deletionsSJ)

    insertionsPerc = float(len(insertions))/totN* 100

    if len(insertions) == 0: insertionsMean = 0
    else:     insertionsMean = numpy.mean(insertions)
    if len(insertionsSJ) == 0: insertionsSJMean = 0
    else:     insertionsSJMean = numpy.mean(insertionsSJ)



    tot = accepted + unaligned + filtered
    acceptPerc = float(accepted)/tot * 100
    unalignedPerc = float(unaligned)/tot * 100
    filteredPerc = float(filtered)/tot * 100
    sys.stderr.write('accepted: %d(%0.2f%%) filtered:%d(%0.2f%%) unaligned:%d(%0.2f%%)\n' % (accepted, 
                                                                                        acceptPerc, 
                                                                                        filtered, 
                                                                                        filteredPerc, 
                                                                                        unaligned,
                                                                                        unalignedPerc,
                                                                                        ))

    sys.stderr.write('Processed %d reads: %d (%0.2f%%) spliced alignments, %d (%0.2f%%) contiguous alignments\n' % (tot, spliced, 
                                                                                                                    float(spliced)/accepted*100,
                                                                                                                    accepted-spliced,
                                                                                                                    (1-float(spliced)/accepted)*100))
    sys.stderr.write('\t%d (%0.2f%%) + strand, %d (%0.2f%%) - strand, %d (%0.2f%%) ambiguous\n' % ( pstrand,
                                                                                                    float(pstrand)/accepted*100,
                                                                                                    nstrand,
                                                                                                    float(nstrand)/accepted*100,
                                                                                                    astrand,
                                                                                                    float(astrand)/accepted*100))
    sys.stderr.write('accepted: %d stranded junctions, %d total junctions\n' % (sSJ, totSJ) ) 
    sjitems = SJdict.items()
    sjitems.sort(key=lambda x: x[1], reverse=1)
    junction_ofile = open( args.junctions_output, 'w' )
    for pair in sjitems[:3]:
        if sSJ == 0:
            sys.stderr.write('No stranded splice junction\n')
        else:
            sys.stderr.write('%s: %d (%0.2f%%)\n' % ( pair[0], pair[1], pair[1]/float(sSJ)*100))
            junction_ofile.write('%s: %d (%0.2f%%)\n' % ( pair[0], pair[1], pair[1]/float(sSJ)*100))
    rem = 0
    for pair in sjitems[3:]:
        rem += pair[1]
        if sSJ == 0:
            pass
        else:
            sys.stderr.write('%s: %d (%0.2f%%)\n' % ( 'others', rem, rem/float(sSJ)*100))
            junction_ofile.write('%s: %d (%0.2f%%)\n' % ( 'others', rem, rem/float(sSJ)*100))
    junction_ofile.close()
    sys.stderr.write('fixed: %d (%0.2f%%) mismatches %d (%0.2f%%) deletions (%0.2f BPs/del) %d (%0.2f%%) insertions (%0.2f BPs/ins)\n' % (mismatches, mismatchesPerc, len(deletions), deletionsPerc, deletionsMean, len(insertions), insertionsPerc, insertionsMean))
    if args.anchorAdjustLim > 0:
        sys.stderr.write('SJs: %d mismatches %d deletions (%0.2f BPs/del) %d insertions (%0.2f BPs/ins)\n' % (mismatchesSJ, len(deletionsSJ), deletionsSJMean, len(insertionsSJ), insertionsSJMean))
