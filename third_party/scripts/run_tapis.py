#! /usr/bin/env python
import pysam,os,numpy,collections,subprocess
from argparse import ArgumentParser, ArgumentTypeError
from SpliceGrapher.formats.loader import loadGeneModels
from SpliceGrapher.formats.GeneModel import *        
from SpliceGrapher.formats import fasta
from SpliceGrapher.SpliceGraph       import *
from SpliceGrapher.shared.GeneModelConverter import *
from bx.intervals.cluster import ClusterTree
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from SpliceGrapher.shared.utils       import *
from SpliceGrapher.plot.PlotterConfig import *
from SpliceGrapher.view.ViewerUtils   import *
from SpliceGrapher.plot.PlotUtils import *
from sys import maxint as MAXINT
from itertools import chain
import os,sys
import warnings  
import numpy, ConfigParser
from matplotlib.patches import Rectangle
#warnings.filterwarnings('ignore')
from matplotlib import rc
rc('text', usetex=True)
import numpy as np
DEFAULT_FONT   = 12
DEFAULT_HEIGHT = 11.0
DEFAULT_WIDTH  = 8.5


parser = ArgumentParser(description='Assemble transcripts from PacBio alignments')
parser.add_argument('-v', '--verbose', dest="verbose",
                    action='store_true', default=False,
                    help='Verbose mode')
parser.add_argument('-p', '--plot', dest="plot",
                    action='store_true', default=False,
                    help='Plot novel gene graphs and poly(A) figures, default is no plotting')
parser.add_argument('-o', '--outdir', dest="outdir",
                    action='store', type=str, default='tapis_out',
                    help='Output directory for TAPIS results, default=tapis_out')
parser.add_argument('-t', '--trimMax', dest="trimMax",
                    action='store', type=int, default=5,
                    help='Maximum length of read trimming to tolerate on 3\' end of reads, default=5')
parser.add_argument('-w', '--w', dest="w",
                    action='store', type=int, default=5,
                    help='Width of peaks when searching for poly(A) sites, default=5')
parser.add_argument('-m', '--minDist', dest="minDist",
                    action='store', type=int, default=20,
                    help='Minimum distance between any two poly(A) sites, default=20')
parser.add_argument('-s', '--minSupport', dest="minSupport",
                    action='store', type=int, default=2,
                    help='Minimum number of trusted reads supporting a poly-A site, default=2')
parser.add_argument('geneModel', action='store', 
                    type=str, help='Gene models annotation file (GFF/GTF)')
parser.add_argument('bamfile', action='store', 
                    type=str, help='Aligned reads file (sorted and indexed)')


            
args = parser.parse_args()
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
if args.plot and not os.path.exists(os.path.join(args.outdir, 'novelGraphs')):
    os.makedirs(os.path.join(args.outdir, 'novelGraphs'))
if args.plot and not os.path.exists(os.path.join(args.outdir, 'polyAFigures')):
    os.makedirs(os.path.join(args.outdir, 'polyAFigures'))

###########################################
geneModel = loadGeneModels( args.geneModel, verbose=args.verbose)
genes = geneModel.getAllGenes()
genes.sort()

def junctionItr( read ):
    """
    Junction generator for spliced reads
    """
    for bidx in xrange( len(read.blocks) - 1 ):
        jct = ( read.blocks[bidx][1]+1, read.blocks[bidx+1][0])
        if abs( jct[0] - jct[1]) > 10:
            yield jct
    
def clusterToGraphP( cluster, chromosome, name='noneThusFar'):
    cluster.sort( key=lambda x: len(x.blocks), reverse=1)
    unresolved = set()
    graph = SpliceGraph.SpliceGraph(name, chromosome, '+')
    pid = None
    i = 0
    gapless = [ ]
    for read in cluster:
        if len(read.blocks) == 1:
            gapless.append(read)
            continue
        graphT = SpliceGraph.SpliceGraph(name+'_temp', chromosome, '+')
        pid = None
        for block in read.blocks:
            graphT.addNode( str(i), block[0], block[1] )
            if pid != None:
                graphT.addEdge( str(pid), str(i) )
            pid = i
            i += 1

        updateRoot(graphT, graph)
        updateRoot(graph, graphT)
        updateLeaf(graph, graphT)
        updateLeaf(graphT, graph)
        graph = graph.union(graphT)
    # merge gapless reads
    gapless.sort(key=lambda x: abs(x.blocks[0][0] - x.blocks[0][1]), reverse=1)
    merged = [ ]
    for read in gapless:
        min1, max1 = read.blocks[0]
        for i,r in enumerate(merged):
            min2, max2 = r
            if min(max1, max2) - max(min1,min2) > 0.75 * min(abs(max1-min1), abs(max2-min2)):
                merged[i] = (min(min1,min2), max(max1,max2))
                break
        else:
            merged.append( (min1,max1))
    gcount=0
    for minpos,maxpos in merged:
        for node in graph.resolvedNodes():
            if minpos >= node.minpos and maxpos <= node.maxpos:
                break
        else:
            graph.addNode('g_%d'% gcount, minpos,maxpos)
            gcount += 1
    return graph

def clusterToGraphN( cluster, chromosome, name='noneThusFar'):
    cluster.sort( key=lambda x: len(x.blocks), reverse=1)
    unresolved = set()
    graph = SpliceGraph.SpliceGraph(name, chromosome, '-')

    pid = None
    i = 0
    for block in cluster[0].blocks[::-1]:
        graph.addNode( str(i), block[1], block[0] )

        if pid != None:
            graph.addEdge( str(pid), str(i) )
        pid = i
        i += 1
    gapless = [ ]
    for read in cluster[1:]:
        if len(read.blocks) == 1:
            gapless.append(read)
            continue
        graphT = SpliceGraph.SpliceGraph(name+'_temp', chromosome, '-')
        pid = None
        for block in read.blocks[::-1]:
            graphT.addNode( str(i), block[1], block[0] )
            if pid != None:
                graphT.addEdge( str(pid), str(i) )
            pid = i
            i += 1

        updateRoot(graph, graphT)
        updateRoot(graphT, graph)
        updateLeaf(graph, graphT)
        updateLeaf(graphT, graph)

        graph = graph.union(graphT)

    # merge gapless reads
    gapless.sort(key=lambda x: abs(x.blocks[0][0] - x.blocks[0][1]), reverse=1)
    merged = [ ]
    for read in gapless:
        min1, max1 = read.blocks[0]
        for i,r in enumerate(merged):
            min2, max2 = r
            if min(max1, max2) - max(min1, min2)  > 0.75 * min(abs(max1-min1), abs(max2-min2)):
                merged[i] = (min(min1,min2), max(max1,max2))
                break
        else:
            merged.append( (min1,max1))
    gcount=0
    for minpos,maxpos in merged:
        for node in graph.resolvedNodes():
            if minpos >= node.minpos and maxpos <= node.maxpos:
                break
        else:
            graph.addNode('g_%d'% gcount, maxpos, minpos)
            gcount += 1

    return graph

def clusterToTranscripts( cluster, strand ):
    """
    Predict transcripts from full-length and non-full length CDNA reads
    on a given strand.
    """
    cluster.sort( key=lambda x: len(x.blocks), reverse=True)
    fl = [ ]
    for read in cluster:
        tDict = dict(read.tags)
        if 'XR' not in tDict or 'XL' not in tDict:
            sys.stderr.write('Warning: trimmed ends info not available, recommend running alignPacBio.py first!\n')
        else:
            if (strand == '+' and tDict['XR'] < 5) or (strand == '-' and \
                                                       tDict['XL'] < 5):
                fl.append(read)
    transcripts = [ ]
    # trusted reads
    for read in fl:
        t = [jct for jct in junctionItr(read)]
        for transcript in transcripts:
            if transcript[len(transcript)-len(t):] == t:
                break
        else:
            transcripts.append( t )
    return transcripts

def clusterReads(bamfile, cluster_treesP, cluster_treesN, readDict):
    bamfile = pysam.Samfile(bamfile, 'rb')
    if args.verbose:
        sys.stderr.write('Clustering reads\n')
    indicator = ProgressIndicator(100000, verbose=args.verbose)
    for read in bamfile:
        chrom = bamfile.getrname(read.tid)
        start = read.blocks[0][0]
        end   = read.blocks[-1][1]
        tDict = dict(read.tags)
        # filter reads with too much trimming on 3' end
        if 'XS' not in tDict:
            sys.stderr.write('Warning: skipping read %s (strand orientation not available in alignment record)\n')
            continue
        if 'XR' not in tDict or 'XL' not in tDict:
            sys.stderr.write('Warning: trimmed ends info not available, recommend running alignPacBio.py first!\n')
        else:
            if tDict['XS'] == '+' and tDict['XR'] > args.trimMax or \
               tDict['XS'] == '-' and tDict['XL'] > args.trimMax :
                continue
        # no rjcts
        #read.setTag('XC', cond)
        readDict[clusterReads.c] = read
        if tDict['XS'] == '+':
            cluster_treesP[chrom].insert(start,end,clusterReads.c)
            clusterReads.c += 1
        else:
            cluster_treesN[chrom].insert(start,end,clusterReads.c)
            clusterReads.c += 1
        indicator.update()
    indicator.finish()

def plotNovel(graph, cluster, outname, shrink_introns=False):
    """
    Plot pacBio reads for a gene
    """
    tsize = 16
    verbose=False
    rc('text', usetex=False)
    minPos = graph.minpos
    maxPos = graph.maxpos
    plt.figure(frameon=False)
    
    titlePadding = getTitlePadding(16)
    topLine      = 0.99-titlePadding
    c = 1
    height         = topLine * (1.0/(1+c) - titlePadding)
    patchDict = {}
    # Plot gene model 
    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    topLine        = topLine - height - titlePadding 

    SpliceGraphView(graph, curAxes,
                    xLimits=(minPos-40, maxPos+40)).plot(xLabels=True)
    curAxes.set_title('PacBio Model', size=tsize)
    curAxes.set_yticks([])
    xlim = (minPos-40, maxPos+40)
    if graph.strand == '-':
        xlim = xlim[::-1]
    curAxes.set_xlim(xlim)
    curAxes.set_xticklabels([])
    
    #plot reads
    if graph.strand == '+':
        cluster.sort( key=lambda x: x.blocks[0][0], reverse=1 )
    else:
        cluster.sort( key=lambda x: x.blocks[-1][1] )
    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])

    curAxes.get_yaxis().set_visible(False)
    for i, r in enumerate(cluster):
        tdict = dict(r.tags)
        if 1:#tdict['XF']:
            curAxes.plot( [r.blocks[0][0], r.blocks[-1][1] ], 
                          [(i+1),(i+1)], 'k', ls='-')
            for b in r.blocks:
                curAxes.add_patch(Rectangle( (b[0]+1, i+.8), b[1] - b[0], .4,
                                              alpha=1, facecolor='k'))

        
    curAxes.set_xlim(xlim)
    curAxes.get_xaxis().get_major_formatter().set_useOffset(False)
    xticks = curAxes.get_xticks()
    curAxes.set_xticklabels([str(int(x)) for x in xticks], size=12)
    #curAxes.set_xticks(xticks)
    plt.ylim((0, 2+i))
    

    topLine        = topLine - height - titlePadding
    curAxes.set_title('Reads', size=tsize)
    #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    
    #plotLegend(patchDict)                    
    #plt.show()
    plt.savefig(outname, dpi=400)
    plt.close()

def plotCluster(referenceGenes, graph, cluster, start, end, geneModel, outname, graph2=None, shrink_introns=False):
    """
    Plot pacBio reads for a gene
    """
    tsize = 16
    refGraph = makeSpliceGraph( referenceGenes[0] )
    # union together multiple reference genes if necessary
    for g in referenceGenes[1:]:
        refGraph = refGraph.union(makeSpliceGraph(g))
    refGraph.annotate()
    verbose=False
    rc('text', usetex=False)
    minPos = min( refGraph.minpos, graph.minpos)
    maxPos = max( refGraph.maxpos, graph.maxpos)
    plt.figure(frameon=False)
    
    titlePadding = getTitlePadding(16)
    topLine      = 0.99-titlePadding
    c = 2
    if graph2: c += 1
    height         = topLine * (1.0/(1+c) - titlePadding)
    patchDict = {}
    # Plot gene model 
    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    topLine        = topLine - height - titlePadding 
    GeneView(referenceGenes, curAxes).plot()
    SpliceGraphView(refGraph, curAxes, xLimits=(minPos-40, maxPos+40)).plot()
    if len(referenceGenes) == 1:
        curAxes.set_title('Gene model for %s' % refGraph.name, size=tsize)
    else:
        curAxes.set_title('Gene models for %s' % refGraph.name, size=tsize)
    #patchDict.update(patches)
    curAxes.set_yticks([])
    xlim = (minPos-40, maxPos+40)
    if graph.strand == '-':
        xlim = xlim[::-1]
    curAxes.set_xlim(xlim)
    #xticks = setXticks(int(min(xlim)), int(max(xlim)))

    curAxes.set_xticklabels([])
    
    # plot prediction
    if 1:
        curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
        topLine        = topLine - height - titlePadding         
        GeneView(referenceGenes, curAxes).plot()
        SpliceGraphView(graph, curAxes, xLimits=(minPos, maxPos)).plot()
        curAxes.set_title('PacBio Model', size=tsize)
        curAxes.set_yticks([])
        curAxes.set_xlim(xlim)
        curAxes.set_xticklabels([])
    #plot reads
    if graph.strand == '+':
        cluster.sort( key=lambda x: x.blocks[0][0], reverse=1 )
    else:
        cluster.sort( key=lambda x: x.blocks[-1][1] )
    # plot isoforms
    if graph2:
        curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
        topLine        = topLine - height - titlePadding
        GeneView(referenceGenes, curAxes).plot()
        IsoformView(graph2, curAxes).plot(isoformLabels=True, 
                                          sortByName=True)
        curAxes.set_title('Isoforms', size=tsize)
        curAxes.set_yticks([])
        curAxes.set_xlim(xlim)
        curAxes.set_xticklabels([])

    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    GeneView(referenceGenes, curAxes).plot()
    curAxes.get_yaxis().set_visible(False)
    #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    for i, r in enumerate(cluster):
        tdict = dict(r.tags)
        if 1:#tdict['XF']:
            curAxes.plot( [r.blocks[0][0], r.blocks[-1][1] ], 
                          [(i+1),(i+1)], 'k', ls='-')
            for b in r.blocks:
                curAxes.add_patch(Rectangle( (b[0]+1, i+.8), b[1] - b[0], .4,
                                              alpha=1, facecolor='k'))

        
    curAxes.set_xlim(xlim)
    curAxes.get_xaxis().get_major_formatter().set_useOffset(False)
    xticks = curAxes.get_xticks()
    curAxes.set_xticklabels([str(int(x)) for x in xticks])
    #curAxes.set_xticks(xticks)
    plt.ylim((0, 2+i))
    

    topLine        = topLine - height - titlePadding

    curAxes.set_title('Reads', size=tsize)
    #plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
    
    #plotLegend(patchDict)                    
    plt.savefig(outname, dpi=400)
    plt.close()

def summarizeClusters(treesP, treesN):
    misAnnot = 0
    totClusters = 0
    novelGenesP = collections.defaultdict(list)
    novelGenesN = collections.defaultdict(list)
    gCounts = collections.defaultdict(int)
    for key in cluster_treesP:
        regions = cluster_treesP[key].getregions()
        for region in regions:
            start, end = region[:2]
            cluster = filter(lambda x: dict(x.tags)['XF'] and dict(x.tags)['XR'] <= 5,
                             [ readDict[rid] for rid in region[-1] ])
            if len(cluster) == 0:
                continue
            totClusters += 1
            referenceGenes = geneModel.getGenesInRange(key.lower(), 
                                                       start, end,
                                                       strand='+')

            if len(referenceGenes) == 0:
                novelGenesP[key].append( (start,end))
            elif len(referenceGenes) > 1:
                misAnnot += len(referenceGenes)
            else:
                gCounts[referenceGenes[0].id] += 1

    for key in cluster_treesN:
        regions = cluster_treesN[key].getregions()
        for region in regions:
            start, end = region[:2]
            cluster = filter(lambda x: dict(x.tags)['XF'] and dict(x.tags)['XL']<=5,
                             [ readDict[rid] for rid in region[-1] ])
            if len(cluster) == 0:
                continue
            totClusters += 1
            referenceGenes = geneModel.getGenesInRange(key.lower(), 
                                                       start, end,
                                                       strand='-')

            if len(referenceGenes) == 0:
                novelGenesN[key].append( (start,end))
            elif len(referenceGenes) > 1:
                misAnnot += len(referenceGenes)
            else:
                gCounts[referenceGenes[0].id] += 1
    
def resolveMultiCluster(cluster, refgenes, strand):
    """
    Try and determine if cluster can be paritioned:
    1) merges two or more misannotated genes
    2) can be separated into multiple genes to to weakly-
       supported gene fusion introns
    """
    # get Poly-A sites

    uniqueOverlap = numpy.zeros(len(refgenes))
    polyALocs = [ list() for _ in xrange(len(refgenes))]
    graphs = [makeSpliceGraph(gene) for gene in refgenes]
    
    for read in cluster:
        minPos, maxPos = read.blocks[0][0], read.blocks[-1][1]
        overLap = numpy.zeros(len(refgenes))
        for i,gene in enumerate(refgenes):
            if max(minPos, refgenes[i].minpos) < min(maxPos, refgenes[i].maxpos):
                overLap[i] += 1
        if strand == '-':
            closest = numpy.argmin( [abs(minPos - refgenes[i].minpos) for i in xrange(len(refgenes))])
            polyALocs[closest].append( minPos )
        else:
            closest = numpy.argmin( [abs(maxPos - refgenes[i].maxpos) for i in xrange(len(refgenes))])
            
            polyALocs[closest].append( maxPos )
        if sum(overLap) == 1:
            uniqueOverlap += overLap
    cluster = [ readDict[rid] for rid in region[-1] ]
    fusionReads = [ ]    
    if numpy.sum(numpy.array( [len(x) for x in polyALocs]) > 0 ) > 1:
        clusters = [ list() for _ in xrange(len(refgenes))]
        if strand == '-':
            pends = [ min(x) if x else numpy.inf for x in polyALocs]
            for read in cluster:
                minPos, maxPos = read.blocks[0][0], read.blocks[-1][1]
                closest = numpy.argmin( [abs(minPos -
                                             refgenes[i].minpos) \
                                         for i in xrange(len(refgenes))])
                # see if read crosses another gene's poly-A end
                for i in xrange(len(refgenes)):
                    if i == closest: continue
                    if maxPos > pends[i] and minPos < pends[i]:
                        fusionReads.append(read)
                        break
                else:
                    clusters[closest].append(read)
                
        else:
            pends = [ max(x) if x else -numpy.inf for x in polyALocs]
            for read in cluster:
                minPos, maxPos = read.blocks[0][0], read.blocks[-1][1]
                closest = numpy.argmin( [abs(maxPos -
                                             refgenes[i].maxpos) \
                                         for i in xrange(len(refgenes))])
                # see if read crosses another gene's poly-A end
                for i in xrange(len(refgenes)):
                    if i == closest: continue
                    if minPos < pends[i] and maxPos > pends[i]:
                        fusionReads.append(read)
                        break
                else:
                    clusters[closest].append(read)
        # return if genes are misannotated
        return (len(fusionReads) == 0,\
            [ (clusters[i], refgenes[i]) for i in xrange(len(refgenes))])

def subsumedIso( iso1, iso2, strand, verbose=False ):
    """
    Checks if iso1 is subsumed by iso2
    """
    if len(iso1) > len(iso2):
        return False
    if len(iso1) == 1 and len(iso2) == 1:
        min1,max1 = iso1[0]
        min2,max2 = iso2[0]
        return min(max1, max2) >= max(min1,min2) 
    
    if strand == '+':
        iso1 = iso1[::-1]
        iso2 = iso2[::-1]
        
    for i,exon in enumerate(iso1):

        min1,max1 = exon
        min2,max2 = iso2[i]
        
        # for last exons, all we care is if
        # exons start at same place
        if i == 0:
            # if it's a single exon isoform,
            # then just need to check if exon starts at or downstream
            # of exon 2 start postion
            if len(iso1) == 1:
                if strand == '+' and min1 < min2 or \
                   strand == '-' and max1 > max2:
                    break
            elif strand == '+' and min1 != min2 or \
               strand == '-' and max1 != max2:
                break

        # check first exon of iso1
        elif i == len(iso1) - 1:
            # if both exons are first exons,
            # then all we care is if
            # both exons end the same place
            if strand == '+' and max1 != max2 or \
               strand == '-' and min1 != min2:
                break
            # otherwise we need to check that iso1 first exon does
            # not start before other exon
            if i < len(iso2) - 1:
                if strand == '+' and min1 < min2 or \
                   strand == '-' and max1 > max2:
                    break
        # internal exon
        else:
            if min1 != min2 or max1 != max2:
                break
    else:
        return True
    return False

def clusterToIsoforms(cluster, strand):
    if strand == '+':
        # sort by reads' smallest positions
        cluster.sort(key=lambda x: x.blocks[0][0])
    else:
        # sort by reads' largest positions
        cluster.sort(key=lambda x: x.blocks[-1][1], reverse=1)

    isos = []
    gapless = [ ]
    # handle spliced reads
    for read in cluster:
        readBlocks = [ (block[0]+1, block[1]) for block in read.blocks]
        # trim some from start exon
        if strand == '+':
            if readBlocks[0][1] - readBlocks[0][0] <= 10:
                readBlocks.pop(0)
            else:
                readBlocks[0] = (readBlocks[0][0]+10, readBlocks[0][1])
        else:
            if readBlocks[-1][1] - readBlocks[-1][0] <= 10:
                readBlocks.pop()
            else:
                readBlocks[-1] = (readBlocks[-1][0], readBlocks[-1][1]-10)

        if len(read.blocks) == 1:
            gapless.append(read)
            continue
        for iso in isos:

            if subsumedIso( readBlocks, iso, strand):
                break
        # no isoform subsumed the read so add
        else:
            isos.append(readBlocks)
                        
    gapless.sort(key=lambda x: abs(x.blocks[0][0] - x.blocks[0][1]), reverse=1)
    merged = [ ]
    for read in gapless:
        min1, max1 = read.blocks[0]
        if strand == '+':
            min1 += 10
        else:
            max1 -= 10
        for i,r in enumerate(merged):
            min2, max2 = r
            if min(max1, max2) - max(min1,min2) > 1:
                merged[i] = (min(min1,min2), max(max1,max2))
                break
        else:
            merged.append( (min1,max1))

    for minpos,maxpos in merged:
        for iso in isos:
            if strand == '+' and minpos >= iso[-1][0]:
                break
            elif strand == '-' and maxpos <= iso[0][1]:
                break
        else:
            isos.append([(minpos,maxpos)])
    return isos

 
def processGene(isos, cluster, gene):
    """
    Summarize isos for gene
    """
    strand = gene.strand
    if strand == '+':
        # sort by isos' smallest positions
        isos.sort(key=lambda x: x[0][0])
    else:
        # sort by reads' largest positions
        isos.sort(key=lambda x: x[-1][1], reverse=1)
    # annotated splice forms
    geneGraph = makeSpliceGraph(gene)
    starts = sorted([(node.minpos, node.maxpos) \
                     for node in geneGraph.getRoots()])
    graph = SpliceGraph.SpliceGraph(gene.id+'_pred', gene.chromosome,
                                    gene.strand)
    for isoNum, iso in enumerate(isos):
        graphT = SpliceGraph.SpliceGraph(gene.id+'_isoG_'+str(isoNum),
                                         gene.chromosome,
                                         gene.strand)
        pid = None
        i = 0
        isoName = '%s_%d' % (gene.id, isoNum+1)
        if strand == '-':
            for block in iso[::-1]:
                graphT.addNode( str(i), block[1], block[0] )
                graphT.getNode(block[1], block[0] ).addIsoform(isoName)
                if pid != None:
                    graphT.addEdge( str(pid), str(i) )
                pid = i
                i += 1
        else:
            for block in iso:
                graphT.addNode( str(i), block[0], block[1] )
                graphT.getNode(block[0], block[1] ).addIsoform(isoName)
                if pid != None:
                    graphT.addEdge( str(pid), str(i) )
                pid = i
                i += 1

        graph = graph.union(graphT, mergeEnds=True)
    graph.annotate()
    knownIntrons = set([(edge.minpos,edge.maxpos) for edge in edgeSet(geneGraph)])
    predIntrons = set( [(edge.minpos, edge.maxpos) for edge in edgeSet(graph)])
    novel = [ ]
    fullLength = [ ]
    
    # check if iso is full-length
    # check if iso is novel    
    for isoNum,iso in enumerate(isos):
        minPos,maxPos = iso[0] if strand == '+' else iso[-1]
        for minG, maxG in starts:
            if strand == '+':
                if len(iso) == 1 and minPos < maxG:
                    fullLength.append(True)
                    break
                elif maxG == maxPos:
                    fullLength.append(True)
                    break
            else:
                if len(iso) == 1 and maxPos > minG:
                    fullLength.append(True)
                    break
                elif maxG == maxPos:
                    fullLength.append(True)
                    break
        else:
            fullLength.append(False)
            
        isoName = '%s_%d' % (gene.id, isoNum+1)
        for node in graph.isoformDict()[isoName]:
            if node.isRoot():
                root = node
            if node.isLeaf():
                leaf = node
        for gisoR in geneGraph.isoformDict().values():
            giso = sorted( [ (node.minpos, node.maxpos) for node in gisoR])
            if subsumedIso( iso, giso, strand):
                # make sure removing isoform does not remove any AS
                novel.append(False)
                for node in graph.resolvedNodes():
                    if isoName in node.isoformSet:
                        node.isoformSet.remove(isoName)
                break
        # no isoform subsumed the predicted isoform
        else:
            novel.append(True)
    #graph2 = graph.union(geneGraph, mergeEnds=True )
    #graph2.annotate()
    #plotCluster([gene], graph, cluster, graph.minpos,
    #            graph.maxpos, geneModel,
    #            '%s.png' % gene.id, graph2=graph2)

    return graph, novel, fullLength

def writeNovelGenes(novelClustersP, novelClustersN):
    if args.verbose:
        sys.stderr.write('Processing novel clusters...\n')
    fout = open(os.path.join(args.outdir,'novelGenes.csv'), 'w')
    fout.write('chromosome\tminPos\tmaxPos\tstrand\texpression\tASforms\n')
    fastaout = open(os.path.join(args.outdir,'novelGenes.fa'), 'w')
    novelReads = {}
    indicator = ProgressIndicator(5000, verbose=args.verbose)
    for key in novelClustersP:
        for cluster in novelClustersP[key]:
            indicator.update()
            longest = sorted( cluster, key=lambda x: len(x.query), reverse=1)[0]
            minpos  = min(x.blocks[0][0] for x in cluster)
            maxpos  = max(x.blocks[-1][1] for x in cluster) 
            nid = '%s_%d_%d_+' % (key, minpos,maxpos)
            novelReads[nid] = len(cluster)
            graph = clusterToGraphP(cluster, key,'none')
            graph.annotate()

            if args.plot:
                plotNovel(graph, cluster, os.path.join(args.outdir, 'novelGraphs', '%s.pdf'%nid))
            fout.write('%s\t%d\t%d\t+\t%d\t%s\n' % (key,
                                                    minpos,
                                                    maxpos,
                                                    len(cluster),
                                                    ','.join(graph.altForms())))
            fastaout.write('>%s\n%s\n' % (nid, longest.query))
    for key in novelClustersN:
        for cluster in novelClustersN[key]:
            indicator.update()
            longest = sorted( cluster, key=lambda x: len(x.query), reverse=1)[0]
            minpos  = min(x.blocks[0][0] for x in cluster)
            maxpos  = max(x.blocks[-1][1] for x in cluster) 
            nid = '%s_%d_%d_-' % (key, minpos,maxpos)
            novelReads[nid] = len(cluster)
            graph = clusterToGraphP(cluster, key,'none')
            graph.annotate()

            if args.plot:
                plotNovel(graph, cluster, os.path.join(args.outdir, 'novelGraphs', '%s.pdf'%nid))
            fout.write('%s\t%d\t%d\t-\t%d\t%s\n' % (key,
                                                    minpos,
                                                    maxpos,
                                                    len(cluster),
                                                    ','.join(graph.altForms())))

            fastaout.write('>%s\n%s\n' % (nid, longest.query))
    fastaout.close()
    fout.close()
    indicator.finish()

def getPeaks(depths):
    """
    Get polyA peaks for given set of depths
    """
    w = args.w
    minDist = args.minDist
    N = len(depths)
    keepGoing = True
    peaks = []
    counts = [ ]
    while True:
        currPeaks = numpy.zeros(len(depths))
        for i in xrange(N):
            for c in peaks:
                if i <= c and c - i + 1 < minDist or \
                   i >= c and i - c + 1 < minDist:
                    break
            else:
                if numpy.sum(depths[ max(0, i-w-1):min(N,i+w)]) >= args.minSupport:
                    currPeaks[i] = depths[i]*2+numpy.median(depths[ max(0, i-w-1):min(N,i+w)])
        if numpy.max(currPeaks) ==0:
            break
        cp = numpy.argmax(currPeaks)
        if cp not in peaks:
            peaks.append(cp)
            counts.append(currPeaks[cp])
    return peaks, counts

def polyA_analysis(polyAMap):
    ofile = open(os.path.join(args.outdir,'polyA_summary.csv'), 'w') 
    ofile.write('gene\tstrand\taligned reads\tnum sites\tlocations\n')
    if args.verbose:
        sys.stderr.write('Performing poly(A) analysis\n')
    polyAcounts = collections.defaultdict(int)    
    polyAGene = collections.defaultdict(list)
    indicator = ProgressIndicator(5000, verbose=args.verbose)
    for gene in polyAMap:
        indicator.update()
        locs = polyAMap[gene]
        minpos = min(locs)
        maxpos = max(locs)
        depth = numpy.zeros( maxpos-minpos+1, int)
        offset = min(locs)
        for loc in locs:
            depth[loc-offset] += 1
        peaks,counts = getPeaks(depth)
        polyAGene[gene] = [x + offset for x in peaks]
        polyAcounts[len(peaks)] += 1
        X = numpy.arange(min(locs), max(locs)+1)


        ofile.write('%s\t%s\t%d\t%d\t%s\n' % (gene.id, gene.strand, len(locs),
                                              len(peaks), ','.join([str(x+offset) for x in peaks])))

        if args.plot and len(peaks) > 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            remove_border(ax)
            ax.bar(X,depth,1, color='.33')
            if gene.strand == '+':
                ax.set_xlim((minpos-10,maxpos+10))
            else:
                ax.set_xlim((maxpos+10,minpos-10))
            xticks = ax.get_xticks()
            ax.set_xticklabels([str(int(x)) for x in xticks], rotation=45, size=15)
            yticks = ax.get_yticks()
            ax.set_yticklabels([str(int(x)) for x in yticks], size=15)
            ax.set_title(gene.id, size=20)
            ax.set_ylabel('Depth', size=18)
            ax.set_xlabel('Coordinates', size=18)
            plt.tight_layout()
            plt.savefig(os.path.join(args.outdir,'polyAFigures', '%s.png' % gene.id))
            plt.close()    
    ofile.close()
    indicator.finish()
    
def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
    """
    Minimize chartjunk by stripping out unnecesasry plot borders and axis ticks
    
    The top/right/left/bottom keywords toggle whether the corresponding plot border is drawn
    """
    ax = axes or plt.gca()
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)
    
    #turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    
    #now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()
    if right:
        ax.yaxis.tick_right()

def writeGtf(geneIsos):
    with open(os.path.join(args.outdir,'assembled.gtf'), 'w') as fout:
        for gene in geneIsos:
            for isoNum, iso in enumerate(geneIsos[gene]):
                trans_id = '%s_%d' % (gene.id, isoNum+1)
                for exonNum, exon in enumerate(iso):
                    fout.write('%s\tprotein_coding\texon\t%d\t%d\t.\t%s\t.\t gene_id "%s"; transcript_id "%s"; exon_number "%d"; gene_name "%s"; transcript_name "%s";\n' %
                    (gene.chromosome,
                     exon[0],
                     exon[1],
                     gene.strand,
                     gene.id,
                     trans_id,
                     exonNum+1,
                     gene.id,
                     gene.id
                 ))
                    fout.write('%s\tprotein_coding\tCDS\t%d\t%d\t.\t%s\t.\t gene_id "%s"; transcript_id "%s"; exon_number "%d"; gene_name "%s"; transcript_name "%s"; protein_id "%s";\n' %
                    (gene.chromosome,
                     exon[0],
                     exon[1],
                     gene.strand,
                     gene.id,
                     trans_id,
                     exonNum+1,
                     gene.id,
                     gene.id,
                     trans_id
                 ))

################################################################################
if __name__ == '__main__':
    # build clusters
    clusterReads.c = 0
    readDict = { }
    clusterDist = 50
    clusterMembers = 1
    cluster_treesP = collections.defaultdict(lambda:ClusterTree(clusterDist, 
                                                                clusterMembers))
    cluster_treesN = collections.defaultdict(lambda:ClusterTree(clusterDist,
                                                                clusterMembers))
    clusterReads(args.bamfile, cluster_treesP, cluster_treesN, readDict)
    keys = list(set.union(*[set(cluster_treesN.keys()), 
                            set(cluster_treesP.keys())]))

    # Transcript assembly
    geneIsos = collections.defaultdict(list)
    geneReads = collections.defaultdict(list)
    novelClustersN = collections.defaultdict(list)
    novelClustersP = collections.defaultdict(list)
    allIsos = []
    allClusters = [ ]
    allReads = 0
    isosP = collections.defaultdict(list)
    isosN = collections.defaultdict(list)
    allNovel = 0
    allMulti = 0
    fusionReads = 0
    # reduce reads to unique isoforms
    # handle multi-gene overlappers
    if args.verbose:
        sys.stderr.write('Processing read clusters\n')
    indicator = ProgressIndicator(5000, verbose=args.verbose)
    allgenes = 0
    for key in keys:
        if key in cluster_treesN:
            regions = cluster_treesN[key].getregions()
            for region in regions:
                indicator.update()
                cluster = [ readDict[rid] for rid in region[-1] ]
                isos = clusterToIsoforms(cluster, '-')
                allIsos.extend(isos)
                allReads += len(cluster)
                start, end = region[:2]
                refGenes = geneModel.getGenesInRange(key,start, end, strand='-')
                allgenes += len(refGenes)
                if len(refGenes) == 1:
                    geneReads[refGenes[0]].extend(cluster)
                    geneIsos[refGenes[0]].extend(isos)
                    isosN[key].extend(isos)
                # fusion reads--ignore
                elif len(refGenes) > 1:
                    allMulti += len(isos)
                    fusionReads += len(cluster)
                # novel cluster
                else:
                    allNovel += len(isos)
                    novelClustersN[key].append(cluster)

        if key in cluster_treesP:
            regions = cluster_treesP[key].getregions()
            for region in regions:
                indicator.update()
                cluster = [ readDict[rid] for rid in region[-1] ]
                isos = clusterToIsoforms(cluster, '+')
                allIsos.extend(isos)
                allReads += len(cluster)
                start, end = region[:2]
                refGenes = geneModel.getGenesInRange(key,start,end,strand='+')
                allgenes += len(refGenes)
                if len(refGenes) == 1:
                    geneReads[refGenes[0]].extend(cluster)
                    geneIsos[refGenes[0]].extend(isos)
                    isosP[key].extend(isos)
                # fusion reads--ignore
                elif len(refGenes) > 1:
                    allMulti += len(isos)
                    fusionReads += len(cluster)
                # novel cluster
                else:
                    allNovel += len(isos)
                    novelClustersP[key].append(cluster)
    indicator.finish()
    if args.verbose:
        sys.stderr.write('Assembled %d transcripts\n' % len(allIsos))
        sys.stderr.write('%d (%0.2f%%) from novel genes\n' % (allNovel, float(allNovel)/len(allIsos)))
        sys.stderr.write('%d (%0.2f%%) spanning multiple genes\n' % (allMulti, float(allMulti)/len(allIsos)))
    writeGtf(geneIsos)
    writeNovelGenes(novelClustersP, novelClustersN)
    
    polyAMap = {}
    for gene in geneReads:
        if gene.strand == '+':
            polyAMap[gene] = [x.blocks[-1][1] for x in geneReads[gene]]
        else:
            polyAMap[gene] = [x.blocks[0][0] for x in geneReads[gene]]
    polyA_analysis(polyAMap)


    #compute AS statistics
    if args.verbose:
        sys.stderr.write('Computing alternative splicing statistics\n')

    graphsdir = os.path.join(args.outdir, 'graphs')
    if not os.path.exists(graphsdir):
        os.makedirs(graphsdir)

    cmd = 'gene_model_to_splicegraph.py -a -A -d %s -m %s' % (graphsdir, os.path.join(args.outdir,'assembled.gtf'))
    if args.verbose:
        cmd += ' -v'
        sys.stderr.write('Executing %s\n' % (cmd) )
    
    subprocess.call(cmd,shell=1)

    cmd = 'splicegraph_statistics.py -a %s -o %s' % (graphsdir, os.path.join(args.outdir, 'AS_statistics.txt'))
    if args.verbose:
        cmd += ' -v'
        sys.stderr.write('Executing %s\n' % (cmd) )
    subprocess.call(cmd,shell=1)
