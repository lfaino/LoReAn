#!/usr/bin/env python

''' 
MASTER THESIS PROJECT
Author: Jose A. Espejo
Date: September 2015 - March 2016

Handles all the mapping. Needs the mapper (STAR/GMAP) in the path.
'''
import os
import math
import subprocess
import dirs_and_files as logistic
from Bio import SeqIO

def star_build(reference, genome_dir, threads, wd):
    '''Builds star reference index'''
    check_file = genome_dir + 'SAindex'
    if os.path.isfile(check_file):# If the ref is there do not build it again
        print ('STAR index existed already: '+
            check_file + ' --- skipping')
        return None
    else:
        ##CHECK IF THE INDEX IS TOO BIG
        fastaFile = open(reference, 'r')
        fastaParsed = SeqIO.parse(fastaFile, 'fasta')
        filter_count = 0
        genome_size = 0
        for record in fastaParsed:
            genome_size += len(str(record.seq))
            filter_count += 1
        fastaFile.close()
         #IF INDEX IS TOO BIG NEEDS TO BE FRAGMENTED
        if filter_count > 5000 and genome_size > 1000000000:
            fragmented_genome = int(math.log((int(genome_size)/int(filter_count)),2)) ##WHAT IS THIS???
            args = ['STAR', '--runThreadN', str(threads), '--runMode', 'genomeGenerate', '--genomeChrBinNbits' , '=',  'min(', str(fragmented_genome) , ')', '--genomeDir', genome_dir, '--genomeFastaFiles', reference]
        else:  #Genome small enough
            args = ['STAR', '--runThreadN', str(threads), '--runMode', 'genomeGenerate', '--genomeDir', genome_dir, '--genomeFastaFiles', reference]
            
        log_name = wd + 'star_build.log'
        log = open(log_name, 'w')
        
        print '\nCMD: ' + ' '.join(args) + '\n'   
        #Call
        try:
            subprocess.check_call(args, stderr = log, cwd = wd)
            print 'STAR database built'
        except:
            print 'STAR build failed'
            raise NameError('')
        
        log.close()
        
    return None

def star_map(reference, reads, threads, genome_dir, max_intron_length, wd):
    '''Calls STAR to map reads to reference'''
       
    prefix = 'STAR_shortreads'
    if isinstance(reads, str): #Only one file
        args = ['STAR', '--runThreadN', threads, '--genomeDir', genome_dir, 
            '--outSAMtype', 'BAM', 'SortedByCoordinate',  '--alignIntronMax', 
            max_intron_length, '--alignMatesGapMax', '1000', '--outFilterMismatchNmax', 
            '15', '--outFileNamePrefix', prefix, '--outSAMstrandField', 'intronMotif',
            '--outFilterIntronMotifs', 'RemoveNoncanonical', '--outSAMattrIHstart', '1',
            '--readFilesIn', reads] #### '--limitGenomeGenerateRAM', '1100000000',
    elif isinstance(reads, list): #Paired end reads
        args = ['STAR', '--runThreadN', threads, '--genomeDir', genome_dir, 
            '--outSAMtype', 'BAM', 'SortedByCoordinate', '--alignIntronMax',
            max_intron_length, '--alignMatesGapMax', '1000', '--outFilterMismatchNmax', 
            '15', '--outFileNamePrefix', prefix, '--outSAMstrandField', 'intronMotif',
            '--outFilterIntronMotifs', 'RemoveNoncanonical', '--outSAMattrIHstart', '1',
            '--readFilesIn', reads[0], reads[1]]
    
    
    filename = wd+prefix+'Aligned.sortedByCoord.out.bam'
   
    if os.path.isfile(filename): 
        print ('STAR alignment file existed already: ' + filename + ' --- skipping\n')
        return filename
    

    
    log_name = wd + 'star_map.log'
    log = open(log_name, 'w')
    
    print '\nCMD: ' + ' '.join(args) + '\n'

    try:
        subprocess.check_call(args, stderr = log, cwd = wd)
        print '>STAR worked. Output is: ' + filename +'\n'
    except:
        print 'STAR did not work properly\n'
        raise NameError('')
    
    
    log.close()
    return filename

def star(reference, fastq_reads, threads, max_intron_length, wd):
    '''Calls the mapper STAR to map fastq_reads to reference.
    First builds the reference index, then maps'''
    #Create dir for genome index
    refer = reference.split('/')[-1]
    genome_dir = wd + refer + '_STARindex/'
    logistic.check_create_dir(genome_dir)
    #Build the reference
    print '\n#BUILD INDEX\n'
    star_build(reference, genome_dir, threads, wd)
        
    #Mapping
    print '\n#MAP\n'
    out_file = star_map(reference, fastq_reads, threads, genome_dir, max_intron_length, wd)
    return out_file
    

def filterLongReads(fastqFilename, min_length, max_length, wd):
    '''Filters out reads longer than length provided'''
    fastqFile = open(fastqFilename, 'r')
    fastq = SeqIO.parse(fastqFile, 'fastq')
    outFilename = wd + 'longreads.filtered.fastq'
    if os.path.isfile(outFilename): 
       print ('Filtered FASTQ existed already: ' + 
              outFilename + ' --- skipping\n')
       return outFilename, 0
    outFile = open(outFilename, 'w')
    filter_count = 0
    for record in fastq:
        if len(str(record.seq)) > int(min_length) and len(str(record.seq)) < int(max_length):
            record.id = str(filter_count)
            SeqIO.write(record, outFile, 'fastq')
            filter_count += 1
    fastqFile.close()
    outFile.close()
    return outFilename, filter_count



def gmap_build(reference, working_dir):
    '''Build the GMAP indexed reference from the fasta reference file    '''
    #  gmap_build -d <genome> -D path/to/DB[-k <kmer size>] <fasta_files...>
    refer = reference.split('/')[-1] + '_GMAPindex'
    
    args = ['gmap_build', '-d', refer, '-D', 
            working_dir, reference]
    # If the ref is there do not build it again
    refer_path = working_dir + refer
    
    if os.path.isdir(refer_path): 
        print ('GMAP database existed already: '+
               refer_path + ' --- skipping')
        return refer
    
    log_name = working_dir + 'gmap_build.log'
    log = open(log_name, 'w')
    
    print '\nCMD: ' + ' '.join(args) + '\n'   
    #Call
    try:
        subprocess.check_call(args, stderr = log)
        print 'GMAP database built'
    except:
        print 'GMAP build failed'
        raise NameError('')
    
    log.close()
    return refer

def gmap(type_out, reference, fastq_reads, threads, out_format, min_intron_length, max_intron_length, exon_length, wd, Fflag):
    '''Calls the mapper GMAP to map fastq_reads to reference.
    First builds the gmap database index with gmap_build(),
    then uses gmap_map() to map'''
    
    #Build the reference
    print '\n#BUILD INDEX\n'
    reference_db = gmap_build(reference, wd)
        
    #Mapping
    print '\n#MAP\n'
    out_file = gmap_map(reference_db, fastq_reads, threads, out_format, min_intron_length, max_intron_length, exon_length, wd, Fflag, type_out)
    return out_file

def gmap_map(reference_database, reads, threads, out_format, min_intron_length, max_intron_length, exon_length, working_dir, Fflag, type_out):
    '''Calls gmap to map reads to reference
    Out_format can be samse of gff3 (2)'''
    if out_format == 'samse':
        filename = working_dir+'gmap.long_reads.sam'
    elif out_format == '2' or out_format == 'gff3_gene':
        if type_out == 'cons':
            filename = working_dir + 'gmap.cluster_consensus.gff3'
        elif type_out == 'trin':
            filename = working_dir + 'gmap.trinity.gff3'
    else:
        raise NameError('Unknown format: '+out_format+'for GMAP. Accepted are samse or 2 (gff3_gene)')
    if os.path.isfile(filename) and os.path.getsize(filename) > 1 :# If the ref is there do not build it again
        print ('STAR index existed already: '+
            filename + ' --- skipping')
    else:
        out_f = open(filename, 'w')
        log_name = working_dir + 'gmap_map.log'
        log = open(log_name, 'w')
        if not Fflag:
            args = ['gmap', '-D', str(working_dir), '-d', str(reference_database),  '-H', str(exon_length), '--cross-species', '--expand-offsets', '1', '-B',  '5', '--min-intronlength', str(min_intron_length), '-n' , '1' , '--microexon-spliceprob' , '1' , '-K', str(max_intron_length), '-f', str(out_format), '-t', str(threads), reads]
            print '\nCMD: ' + ' '.join(args) + '\n'
            try:
                subprocess.check_call(args, stdout = out_f, stderr = log)
                print '>GMAP worked. Output is: ' + filename +'\n'
            except:
                print 'GMAP did not work properly\n'
                raise NameError('')
        else:
            args = ['gmap', '-D', str(working_dir), '-d', str(reference_database),  '-H', str(exon_length), '--cross-species', '--expand-offsets', '1', '-B',  '5', 
                '--min-intronlength', str(min_intron_length), '-F','-Y', '--microexon-spliceprob' , '1' ,'-n', '1', '-K', str(max_intron_length),
                '-f', str(out_format), '-t', str(threads), reads]
            print '\nCMD: ' + ' '.join(args) + '\n'
            try:
                subprocess.check_call(args, stdout = out_f, stderr = log)
                print '>GMAP worked. Output is: ' + filename +'\n'
            except:
                print 'GMAP did not work properly\n'
                raise NameError('')
        out_f.close()
        log.close()
    return filename
      
def samtools_view(sam_file, wd):
    '''SAM to BAM'''
    bam_filename = sam_file + '.bam'
    args = ['samtools', 'view', '-bS', '-o', bam_filename, sam_file]
    if os.path.isfile(bam_filename): 
        print ('BAM file existed already: ' + bam_filename + ' --- skipping\n')
        return bam_filename
    log_name = wd + 'samtools_view.log'
    log = open(log_name, 'w')   
    print '\nCMD: ' + ' '.join(args) + '\n'      
    try:
        subprocess.check_call(args, stderr = log)
        print '> SAM converted to BAM in ' + bam_filename
    except:
        print 'Samtools view failed'
        raise NameError('')
    
    log.close()
    return bam_filename

def samtools_sort(bam_file, threads, wd):
    '''BAM sorting'''
    s_bam_filename = bam_file + '.sorted'
    args = ['samtools', 'sort', '-@', str(threads), bam_file , s_bam_filename]
    s_bam_filename = s_bam_filename
    if os.path.isfile(s_bam_filename): 
        print ('Sorted BAM file existed already: ' + s_bam_filename + 
               '.bam --- skipping\n')
        return s_bam_filename
    log_name = wd + 'samtools_sort.log'
    log = open(log_name, 'w')  
    print '\nCMD: ' + ' '.join(args) + '\n'      
    try:
        subprocess.check_call(args, stderr = log)
        
        print '> BAM sorted in ' + s_bam_filename + '\n'
    except:
        print 'Samtools sort failed'
        raise NameError('')
    log.close()
    sor_bam_filename = s_bam_filename + ".bam"
    return sor_bam_filename
   
def sam_to_sorted_bam (sam_file, threads, wd):
    print '#SAM to BAM'
    bam_filename = samtools_view(sam_file, wd)
    
    print'\n#Sorting BAM'
    s_bam_filename = samtools_sort(bam_filename, threads, wd)
    return s_bam_filename
