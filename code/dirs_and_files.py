#!/usr/bin/env python3

import errno
import os
import subprocess

#==========================================================================================================
# COMMANDS LIST

GFFREAD = 'gffread -o- -T %s'

GTF2BED = 'gtf2bed.py transcript %s'

BEDTOOLS = 'bedtools bamtobed -split -bed12 -i %s'

#==========================================================================================================


def check_dir(path):
    '''Checks if a directory exists'''
    if not os.path.isdir(path):  # Fast5 dir needs to exist.
        raise IOError(path + ' directory not found')


def check_file(path_file):
    '''Check if a file exists'''
    if not os.path.isfile(path_file):
        raise IOError(path_file + ' file not found')


def check_create_dir(path):
    """
    Checks a directory and creates it if it does not exist
    """

    try:
        os.makedirs(path)  # Try to make it, otherwise don't do anything
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def copy_file(in_file, directory):
    '''Copies a file into a directory'''
    args = ['cp', '-f', in_file, directory]
    try:
        subprocess.check_call(args)
    except:
        raise NameError('')


def catTwoBeds(gmap, evm, outFilename):
    '''convert in to bed12 and concatenates the two bed12 files'''
    gtf = evm + ".gtf"
    bed12_evm = evm + ".bed12"
    bed12file = open(bed12_evm, "w")
    gtffile = open(gtf, "w")
    gffread_con = GFFREAD % (evm)
    gffread_call = subprocess.Popen(gffread_con, stdout=gtffile, shell=1)
    gffread_call.communicate()
    gft2bed = GTF2BED % (gtf)
    evm_call = subprocess.Popen(gft2bed, stdout=bed12file, shell=1)
    evm_call.communicate()
    bed12_gmap = gmap + ".bed12"
    bed12gmapfile = open(bed12_gmap, "w")
    bedtools = BEDTOOLS %(gmap)
    bedtools_call = subprocess.Popen(bedtools, stdout=bed12gmapfile, shell=1)
    bedtools_call.communicate()
    inFile1 = open(bed12_gmap, 'r')
    inFile2 = open(bed12_evm, 'r')
    outFile = open(outFilename, 'w')
    for File in [inFile1, inFile2]:
        for line in File:
            outFile.write(line)
    inFile1.close()
    inFile2.close()
    outFile.close()
    outNameNew = outFilename + 'new.bed'
    lastFile = open(outFilename, 'r')
    o = open(outNameNew, 'w')
    countLine = 0
    for line in lastFile:
        countLine += 1
        linenew = line.split('\t')
        if "evm" in linenew[3]:
            o.write(line)
        else:
            linenew[3] = str(countLine)
            aline = '\t'.join(linenew)
            o.write(aline)
    o.close()
    bed12gmapfile.close()
    bed12file.close()
    gtffile.close()
    return outNameNew


def catTwoFasta(trinity, consens, allSeq, wd):
    '''Concatenates the two fasta file into one output'''
    outFileFasta = wd + "/allFasta.fasta.clean"
    if os.path.isfile(outFileFasta):

        allOutFasta = outFileFasta + ".long.clean"
        cat_con = ['cat', trinity, allSeq, consens]
        allOutFastafile = open(allOutFasta, "w")
        # sys.stdout.write (cat_con)
        cat_call = subprocess.Popen(cat_con, stdout=allOutFastafile)
        cat_call.communicate()
        outFileFasta = allOutFasta
        allOutFastafile.close()
    else:
        allOutFastafile = open(outFileFasta, "w")
        cat_con = ['cat', trinity, allSeq, consens]
        # sys.stdout.write (cat_con)
        # sys.stdout.write ("in")
        cat_call = subprocess.Popen(cat_con, stdout=allOutFastafile)
        cat_call.communicate()
        allOutFastafile.close()

    return outFileFasta
