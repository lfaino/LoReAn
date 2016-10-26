#!/usr/bin/env python

''' 
MASTER THESIS PROJECT
Author: Jose A. Espejo
Date: September 2015 - March 2016

Directory and file handling
'''
import os
import sys
import errno
import subprocess

def check_dir(path):
    '''Checks if a directory exists'''
    if not os.path.isdir(path): #Fast5 dir needs to exist.
        raise IOError(path + ' directory not found')
  
def check_file(path_file):
    '''Check if a file exists'''
    if not os.path.isfile(path_file):
        raise IOError(path_file + ' file not found')
        
def check_create_dir(path):
    '''Checks a directory and creates it if it does not exist'''
    try:
        os.makedirs(path) #Try to make it, otherwise don't do anything
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def copy_file(in_file, directory):
    '''Copies a file into a directory'''
    args = ['cp', '-f', in_file, directory]   
    print '\nCMD: ' + ' '.join(args) + '\n'
    try:
        subprocess.check_call(args)
    except:
        print 'Could not move file: '+ in_file + '\n'
        raise NameError('')
    

def catTwoBeds(gmap, evm, outFilename):
    '''convert in to bed12 and concatenates the two bed12 files'''
    gtf = evm + ".gtf"
    bed12_evm = evm + ".bed12"
    gffread_con = ['gffread', '-o-', '-T', evm ]
    gffread_call = subprocess.Popen(gffread_con, stdout = file(gtf, "w"))
    gffread_call.communicate()
    gft2bed = ['gtf2bed.py',  'transcript', gtf ]
    evm_call = subprocess.Popen(gft2bed, stdout = file(bed12_evm, "w"))
    evm_call.communicate()
    bed12_gmap = gmap + ".bed12"
    gmap_con = ['bedtools','bamtobed', '-split', '-bed12', '-i',  gmap ]
    gmap_call = subprocess.Popen(gmap_con, stdout = file(bed12_gmap , "w"))
    gmap_call.communicate()
    '''Concatenates the two inFiles into the outFile'''
    inFile1 = open(bed12_gmap , 'r')
    inFile2 = open(bed12_evm, 'r')
    outFile = open(outFilename, 'w')
    for File in [inFile1, inFile2]:
        for line in File:
            outFile.write(line)
    inFile1.close()
    inFile2.close()
    outFile.close()
    return outFilename

def catTwoFasta(trinity, consens, wd):
    '''Concatenates the two fasta file into one output'''
    outFileFasta = wd + "/allFasta.fasta.clean"
    cat_con = ['cat', trinity, consens ]
    cat_call = subprocess.Popen(cat_con, stdout = file(outFileFasta , "w"))
    cat_call.communicate()
    return outFileFasta

def catTwoFiles(inFilename1, inFilename2, outFilename):
    '''Concatenates the two inFiles into the outFile'''
    inFile1 = open(inFilename1, 'r')
    inFile2 = open(inFilename2, 'r')
    outFile = open(outFilename, 'w')
    for File in [inFile1, inFile2]:
        for line in File:
            outFile.write(line)
    inFile1.close()
    inFile2.close()
    outFile.close()
    return outFilename


    
    
    
    

###############
###MAIN###
###############

def main():
    '''Main body of the function'''
    gmap = os.path.abspath(sys.argv[1])
    evm = os.path.abspath(sys.argv[2])
    outFilename = os.path.abspath(sys.argv[3])
    catTwoFiles(gmap, evm, outFilename)
       
    print '\n\n\n###############\n###FINISHED###\n###############\n\n'

    
    
    


if __name__ == '__main__':
    main()
