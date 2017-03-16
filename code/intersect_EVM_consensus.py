#!/usr/bin/env python

'''
MASTER THESIS PROJECT
Author: Jose A. Espejo
Date: September 2015 - March 2016


USAGE: intersect_EVM_consensus.py consensus.gff3 evm.gff3
'''

###############
###IMPORTS###
###############


import sys
import subprocess
#import argparse
#import os
#import re
#import shutil
#from dirs_and_files import check_create_dir
#from Bio import SeqIO
from BCBio import GFF
from operator import itemgetter
import dirs_and_files as logistic


###############
###FUNCTIONS###
###############

def parseConsensus(consensusFilename):
    print("\t###PARSE THE CONSENSUS AFTER ASSEMBLY###\n")
    '''Reads the consensus GFF3 file and returns 2 dictionaries,
    one with features and other with locations'''
    consensusFile = open(consensusFilename, 'r')
    consensusDict = {}
    locationDict = {}

    for record in GFF.parse(consensusFile, target_lines=1):
        split = record.features[0].id.split('.')
        iden = split[0] + '.' + split[1][-1]
        if record.features[0].type == 'gene':
            consensusDict[iden] = [record.features[0].id]
            hit = record.id
            start = int(record.features[0].location.start) + 1
            end = int(record.features[0].location.end)
            strand = record.features[0].location.strand
            if strand == 1:
                strand = '+'
            elif strand == -1:
                strand = '-'
            locationDict[iden] = [hit, start, end, strand]
        elif record.features[0].type == 'mRNA':
            consensusDict[iden].append(record.features[0].id)
        else:
            parent = record.features[0].qualifiers['Parent'][0]
            split = parent.split('.')
            iden = split[0] + '.' + split[1][-1]
            consensusDict[iden].append(record.features[0].id)
    consensusFile.close()
    return consensusDict, locationDict


def location2Bed(locationDict):
    print("\t###location2Bed###\n")
    '''From the location dictionary, returns an array with elements
    as if they were sorted bed files:
    Chr \t start \t end \t identifier'''
    bed_records = []
    for key, element in list(locationDict.items()):
        record = element + [key]
        bed_records.append(record)  # Get records
    sorted_bed_records = sorted(bed_records, key=itemgetter(0, 1))
    outBedName = 'consensus.bedrecords.tmp'
    outBed = open(outBedName, 'w')
    for record in sorted_bed_records:
        # Locations are not str, and if I do it before I cannot sort properly
        outBed.write('\n' + ('\t'.join(map(str, record))))
    outBed.close()
    return outBedName


def geneGFF3(gff3Filename):
    print("\t###geneGFF3###\n")
    '''From a gff3 file, writes a "tmp" file with only genes'''
    gff3File = open(gff3Filename, 'r')
    tmpFilename = gff3Filename + '.gene.tmp'
    tmpFile = open(tmpFilename, 'w')
    GFF.write(GFF.parse(gff3File, limit_info=dict(gff_type=['gene'])), tmpFile)
    gff3File.close()
    tmpFile.close()
    return tmpFilename


def intersectConsensusEVM(consensusBed, evmGeneGFF):
    print("\t###intersectConsensusEVM###\n")
    '''Uses bedtools intersect to get the consensus genes that are
    not in the EVM file, returning only the identifiers'''
    args = ['bedtools', 'intersect', '-a', consensusBed, '-b', evmGeneGFF,
            '-v']
    intersect = subprocess.check_output(args)
    geneRecords = intersect.splitlines()

    geneIDs = []
    for record in geneRecords:
        geneID = record.split('\t')[-1]
        geneIDs.append(geneID)
    return geneIDs


def getFeatures(geneIDs, consensusDict):
    print("\t###getFeatures###\n")
    '''Returns an array with the features to write, when the geneID matches'''
    featureIDlist = []
    for gene in geneIDs:
        featureIDlist += consensusDict[gene]
    return featureIDlist


def featuresGFF3(features, originalGFF3Filename, outFilename):
    print("\t###featuresGFF3###\n")
    '''Writes the features from the list that are in the original GFF3
    to the file in filename'''
    inFile = open(originalGFF3Filename, 'r')
    outFile = open(outFilename, 'w')
    recordsOut = []
    for record in GFF.parse(inFile, target_lines=1):
        if record.features[0].id in features:
            recordsOut.append(record)
    GFF.write(recordsOut, outFile)

    inFile.close()
    outFile.close()


def parseGFF3_out(gff3Filename):
    print("\t###parseGFF3_out###\n")
    '''Parses the GFF3 file resulting from featuresGFF3 as it writes a lot of
    annoying ##'''
    outFilename = '.'.join(gff3Filename.split('.')[:-1]) + '.gff3'
    inFile = open(gff3Filename, 'r')
    outFile = open(outFilename, 'w')
    for i, line in enumerate(inFile):
        if i == 0:
            outFile.write(line)
        else:
            if line.startswith('##'):
                continue
            else:
                outFile.write(line)
    inFile.close()
    outFile.close()
    return outFilename


###############
###MAIN###
###############

def main():
    '''Main body of the function'''
    consensusFilename = sys.argv[1]
    modelFilename = sys.argv[2]
    consensusDict, locationDict = parseConsensus(consensusFilename)
    consensusBed = location2Bed(locationDict)
    geneModelFilename = geneGFF3(modelFilename)
    geneIDs = intersectConsensusEVM(consensusBed, geneModelFilename)

    features = getFeatures(geneIDs, consensusDict)
    outTmp = 'consensus.intersect.tmp'
    featuresGFF3(features, consensusFilename, outTmp)
    consensusOut = parseGFF3_out(outTmp)
    logistic.catTwoFiles(
        consensusOut,
        modelFilename,
        'consensus.model.all.gff3')

    print('\n\n\n################\n####FINISHED####\n################\n\n')


if __name__ == '__main__':
    main()
