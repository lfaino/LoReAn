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
from subprocess import Popen, PIPE
from pybedtools import BedTool

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

def catTwoFiles(trinityGFF3, evm_gff3, wd):
    '''Concatenates the two inFiles into the outFile'''
    
    skipped_lines = 0
    first_skipped_line = 0
    
    p1 = subprocess.Popen(['gffread' , '-o-' , '-T', evm_gff3 ], stdout=PIPE)
    p2 = subprocess.Popen(['gffread' , '-o-' , '-T', trinityGFF3 ], stdout=PIPE)
    
    p3 = subprocess.Popen(['gtf2bed.py', '/dev/stdin'], stdin=p1.stdout, stdout=file(evm_gff3 + ".bed12", "w"))
    out1 = p3.communicate()
    p4 = subprocess.Popen(['gtf2bed.py', '/dev/stdin'], stdin=p2.stdout, stdout=file(trinityGFF3 + ".bed12", "w"))
    out2 = p4.communicate()
    
    file1 = BedTool(evm_gff3 + ".bed12")
    file2 = BedTool(trinityGFF3 + ".bed12")
    
    intergenic_snps = file2.intersect(file1, v = True)
    intergenic_snps.saveas(wd + '/merged.bed12')
    
    outFile = wd + '/combinedGMAP_EVM.bed12'
    
    o = open(outFile,'w')
    a = open(wd + "/merged.bed12", 'r')
    for line in a:
        o.write(line)
       
    f = open(evm_gff3 + ".bed12", 'r')
    for line in f:
        o.write(line)
    
    output_name = wd + '/combinedGMAP_EVM.gff3'
    out = open( output_name, 'w' )
    i = 0
    for i, line in enumerate( file( wd + '/combinedGMAP_EVM.bed12' ) ):
        complete_bed = False
        line = line.rstrip( '\r\n' )
        if line and not line.startswith( '#' ) and not line.startswith( 'track' ) and not line.startswith( 'browser' ):
            try:
                elems = line.split( '\t' )
                if len( elems ) == 12:
                    complete_bed = True
                chrom = elems[0]
                if complete_bed:
                    feature = "mRNA"
                else:
                    try:
                        feature = elems[3]
                    except:
                        feature = 'feature%d' % ( i + 1 )
                start = int( elems[1] ) + 1
                end = int( elems[2] )
                try:
                    score = elems[4]
                except:
                    score = '0'
                try:
                    strand = elems[5]
                except:
                    strand = '+'
                try:
                    group = elems[3]
                except:
                    group = 'group%d' % ( i + 1 )
                if complete_bed:
                    out.write( '%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\tID=%s;\n' % ( chrom, feature, start, end, score, strand,  group  ) )
                else:
                    out.write( '%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\t%s;\n' % ( chrom, feature, start, end, score, strand, group  ) )
                if complete_bed:
                    # We have all the info necessary to annotate exons for genes and mRNAs
                    block_count = int( elems[9] )
                    block_sizes = elems[10].split( ',' )
                    block_starts = elems[11].split( ',' )
                    for j in range( block_count ):
                        exon_start = int( start ) + int( block_starts[j] )
                        exon_end = exon_start + int( block_sizes[j] ) - 1
                        out.write( '%s\tbed2gff\texon\t%d\t%d\t%s\t%s\t.\tParent=%s\n' % ( chrom, exon_start, exon_end, score, strand, group ) )
            except:
                skipped_lines += 1
                if not first_skipped_line:
                    first_skipped_line = i + 1
    out.close()
    
    
    
    f.close()
    o.close() 
    
    
    return output_name
    

    
    
    
    

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
