#!/usr/bin/env python3

import sys
import os
import math
import subprocess
import dirs_and_files as logistic
from Bio import SeqIO

#==========================================================================================================
# COMMANDS LIST

GMAP_BULD  = 'gmap_build -k 13 -d %s -D %s %s'

GMAP_GFF = 'gmap -D %s  -d %s --trim-end-exons %s --cross-species --expand-offsets 1 -B 5 --min-intronlength %s -n  1 \
--microexon-spliceprob 1 -F -K  %s -t %s -f %s %s'

GMAP_SAM = 'gmap -D %s  -d %s --trim-end-exons %s --cross-species --expand-offsets 1 -B 5 --min-intronlength %s -n  1 \
--microexon-spliceprob 1 -K  %s -t %s -f %s %s'

STAR = 'STAR --runThreadN %s --genomeDir %s --outSAMtype BAM Unsorted --alignIntronMax %s --alignMatesGapMax %s \
--outFilterMismatchNmax 15 --outFileNamePrefix %s --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical \
--outSAMattrIHstart 1 --outFilterMultimapNmax 3 --readFilesIn %s'

STAR_BUILD = 'STAR --runThreadN %s --runMode genomeGenerate --genomeDir %s --genomeSAindexNbases 6 --genomeFastaFiles %s'

SAMTOOLS_VIEW = 'samtools view -bS -o %s %s'

SAMTOOLS_SORT = 'samtools sort -@ %s %s %s'

#==========================================================================================================

def gmap_map(reference_database, reads, threads, out_format, min_intron_length, max_intron_length, exon_length, working_dir, Fflag, type_out, verbose):
    '''Calls gmap to map reads to reference
    Out_format can be samse of gff3 (2)'''
    global filename
    if out_format == 'samse':
        filename = working_dir + 'gmap.long_reads.sam'
    elif out_format == '2' or out_format == 'gff3_gene':
        if type_out == 'cons':
            filename = working_dir + 'gmap.cluster_consensus.gff3'
        elif type_out == 'trin':
            filename = working_dir + 'gmap.trinity.gff3'
    else:
        raise NameError(
            'Unknown format: ' +
            out_format +
            'for GMAP. Accepted are samse or 2 (gff3_gene)')
    if os.path.isfile(filename) and os.path.getsize(filename) > 1:  # If the ref is there do not build it again
        print(('STAR index existed already: ' +
               filename + ' --- skipping'))
    else:
        out_f = open(filename, 'w')
        log_name = working_dir + 'gmap_map.log'
        log = open(log_name, 'w')
        if not Fflag:
            cmd = GMAP_SAM % (working_dir, reference_database, exon_length, min_intron_length, max_intron_length, threads, out_format, reads)
            try:
                if verbose:
                    sys.stderr.write('Executing: %s\n\n' % cmd)
                gmapmap = subprocess.Popen(cmd, stdout=out_f, stderr=log, shell=1)
                gmapmap.communicate()
                # print '>GMAP worked. Output is: ' + filename +'\n'
            except:
                # print 'GMAP did not work properly\n'
                raise NameError('')
        else:
            cmd = GMAP_GFF % (working_dir, reference_database, exon_length, min_intron_length, max_intron_length, threads, out_format, reads)
            try:
                if verbose:
                    sys.stderr.write('Executing: %s\n\n' % cmd)
                gmapmap = subprocess.Popen(cmd, stdout=out_f, stderr=log, shell=1)
                gmapmap.communicate()
                # print '>GMAP worked. Output is: ' + filename +'\n'
            except:
                # print 'GMAP did not work properly\n'
                raise NameError('')
        out_f.close()
        log.close()
    return filename

def star_build(reference, genome_dir, threads, wd, verbose):
    '''Builds star reference index'''
    check_file = genome_dir + 'SAindex'
    if os.path.isfile(check_file):  # If the ref is there do not build it again
        print(('STAR index existed already: ' +
               check_file + ' --- skipping'))
        return None
    else:
        # CHECK IF THE INDEX IS TOO BIG
        fastaFile = open(reference, 'r')
        fastaParsed = SeqIO.parse(fastaFile, 'fasta')
        filter_count = 0
        genome_size = 0

        for record in fastaParsed:
            genome_size += len(str(record.seq))
            filter_count += 1
        fastaFile.close()
        # IF INDEX IS TOO BIG NEEDS TO BE FRAGMENTED

        if filter_count > 5000 and genome_size > 1000000000:
            fragmented_genome = int(
                math.log(
                    (int(genome_size) / int(filter_count)),
                    2))  # WHAT IS THIS???
            STAR_BUILD_NEW = STAR_BUILD + ' --genomeChrBinNbits =, min 18, %s'
            cmd = STAR_BUILD_NEW % (threads, genome_dir, reference, fragmented_genome)
        else:  # Genome small enough
            cmd = STAR_BUILD % (threads, genome_dir, reference)

        log_name_err = wd + 'star_build.err.log'
        log_err = open(log_name_err, 'w')
        log_name = wd + 'star_build.log'
        log = open(log_name, 'w')

        try:
            if verbose:
                sys.stderr.write('Executing: %s\n\n' % cmd)
            star_build = subprocess.Popen(cmd, stdout=log, stderr=log_err, shell=1, cwd=wd)
            star_build.communicate()
        except:
            print('STAR build failed')
            raise NameError('')

        log.close()
        log_err.close()

    return None

def star_map(reads, threads, genome_dir, max_intron_length, wd, verbose):
    '''
    mapping short reads using STAR
    :param reads:
    :param threads:
    :param genome_dir:
    :param max_intron_length:
    :param wd:
    :return:
    '''

    global args
    prefix = 'STAR_shortreads'
    if isinstance(reads, str):  # Only one file
        cmd = STAR % (threads, genome_dir, max_intron_length, max_intron_length, prefix, reads) # '--limitGenomeGenerateRAM', '1100000000',
    elif isinstance(reads, list):  # Paired end reads
        cmd = STAR % (threads, genome_dir, max_intron_length, max_intron_length, prefix, reads[0], reads[1])
    filename = wd + prefix + 'Aligned.out.bam'

    if os.path.isfile(filename):
        print((
            'STAR alignment file existed already: ' +
            filename +
            ' --- skipping\n'))
        return filename

    log_name_err = wd + 'star.err.log'
    log_err = open(log_name_err, 'w')
    log_name = wd + 'star.log'
    log = open(log_name, 'w')

    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        star = subprocess.Popen(cmd, stdout=log, stderr=log,  shell=1, cwd=wd)
        star.communicate()
        # print '>STAR worked. Output is: ' + filename +'\n'
    except:
        raise NameError('')
    log.close()
    log_err.close()

    return filename


def star(reference, fastq_reads, threads, max_intron_length, wd, verbose):
    '''Calls the mapper STAR to map fastq_reads to reference.
    First builds the reference index, then maps'''
    # Create dir for genome index
    refer = reference.split('/')[-1]
    genome_dir = wd + refer + '_STARindex/'
    logistic.check_create_dir(genome_dir)
    # Build the reference
    print('\t###BUILD INDEX###\n')
    star_build(reference, genome_dir, threads, wd, verbose)

    # Mapping
    print('\t###MAP###\n')
    out_file = star_map(
        fastq_reads,
        threads,
        genome_dir,
        max_intron_length,
        wd, verbose)
    return out_file

def gmap_build(reference, working_dir, verbose):
    """
    Build the GMAP indexed reference from the fasta reference file
    :param reference:
    :param working_dir:
    :return:
    """
    if '/' in reference:
        refer = reference.split('/')[-1] + '_GMAPindex'
    else:
        refer = reference + '_GMAPindex'

    cmd = GMAP_BULD % (refer, working_dir, reference)
    # If the ref is there do not build it again
    refer_path = working_dir + refer
    if os.path.isdir(refer_path):
        print(('GMAP database existed already: ' +
               refer_path + ' --- skipping'))
        return refer
    log_name = working_dir + 'gmap_build.log'
    log = open(log_name, 'w')

    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        gmap_build = subprocess.Popen(cmd, stdout=log, stderr=log, shell =1)
        gmap_build.communicate()
    except:
        print('GMAP build failed')
        raise NameError('')
    log.close()
    return refer


def gmap(type_out, reference, fastq_reads, threads, out_format, min_intron_length, max_intron_length, exon_length, wd, verbose, Fflag):
    '''Calls the mapper GMAP to map fastq_reads to reference.
    First builds the gmap database index with gmap_build(),
    then uses gmap_map() to map'''
    # Build the reference
    print('\t###BUILD INDEX###\n')
    reference_db = gmap_build(reference, wd, verbose)
    # Mapping
    print('\t###MAP###\n')
    out_file = gmap_map(
        reference_db,
        fastq_reads,
        threads,
        out_format,
        min_intron_length,
        max_intron_length,
        exon_length,
        wd,
        Fflag,
        type_out, verbose)
    return out_file

def samtools_view(sam_file, wd, verbose):
    '''SAM to BAM'''
    bam_filename = sam_file + '.bam'
    cmd = SAMTOOLS_VIEW % (bam_filename, sam_file)
    if os.path.isfile(bam_filename):
        print(('BAM file existed already: ' + bam_filename + ' --- skipping\n'))
        return bam_filename
    log_name = wd + 'samtools_view.log'
    log = open(log_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        samtools = subprocess.Popen(cmd, stderr=log, shell=1)
        samtools.communicate()
    except:
        raise NameError('')

    log.close()
    return bam_filename


def samtools_sort(bam_file, threads, wd, verbose):
    '''
    run a sorting of a bam file
    :param bam_file:
    :param threads:
    :param wd:
    :return:
    '''
    s_bam_filename = bam_file + '.sorted'

    if not os.path.isfile(s_bam_filename + '.bam'):
        cmd = SAMTOOLS_SORT % (threads, bam_file, s_bam_filename)
        s_bam_filename = s_bam_filename
        if os.path.isfile(s_bam_filename):
            print(('Sorted BAM file existed already: ' + s_bam_filename +
                   '.bam --- skipping\n'))
            return s_bam_filename
        log_name = wd + 'samtools_sort.log'
        log = open(log_name, 'w')
        try:
            if verbose:
                sys.stderr.write('Executing: %s\n' % cmd)
            samtools = subprocess.Popen(cmd, stderr=log, shell=1)
            samtools.communicate()
        except:
            raise NameError('')
        log.close()
    sor_bam_filename = s_bam_filename + ".bam"
    return sor_bam_filename

def sam_to_sorted_bam(sam_file, threads, wd, verbose):
    print('\t###SAM to BAM###\n')
    bam_filename = samtools_view(sam_file, wd, verbose)

    print('\t###SORTING BAM###\n')
    s_bam_filename = samtools_sort(bam_filename, threads, wd, verbose)
    return s_bam_filename
