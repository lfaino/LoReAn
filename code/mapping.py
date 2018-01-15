#!/usr/bin/env python3

import datetime
import math
import os
import subprocess
import sys

import dirsAndFiles as logistic
from Bio import SeqIO

#==========================================================================================================
# COMMANDS LIST

GMAP_BULD  = 'gmap_build -k 13 -d %s -D %s %s'

GMAP_GFF = 'gmap -D %s  -d %s --trim-end-exons %s --cross-species --expand-offsets 1 -B 5 --min-intronlength %s -n  1 \
--microexon-spliceprob 1 -F -K  %s -t %s -f %s %s'

GMAP_SAM = 'gmap -D %s  -d %s --trim-end-exons %s --cross-species --expand-offsets 1 -B 5 --min-intronlength %s -n  1 \
--microexon-spliceprob 1 -K  %s -t %s -f %s %s'

STAR_SINGLE = 'STAR --runThreadN %s --genomeDir %s --outSAMtype BAM Unsorted --alignIntronMax %s --alignMatesGapMax %s ' \
       '--outFilterMismatchNmax 15 --outFileNamePrefix %s --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical ' \
       '--outSAMattrIHstart 1 --outFilterMultimapNmax 3 --readFilesIn %s'

STAR_PAIRED = 'STAR --runThreadN %s --genomeDir %s --outSAMtype BAM Unsorted --alignIntronMax %s --alignMatesGapMax %s ' \
       '--outFilterMismatchNmax 15 --outFileNamePrefix %s --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical ' \
       '--outSAMattrIHstart 1 --outFilterMultimapNmax 3 --readFilesIn %s %s'

STAR_BUILD = 'STAR --runThreadN %s --runMode genomeGenerate --genomeDir %s --genomeSAindexNbases 6 --genomeFastaFiles %s'

SAMTOOLS_VIEW = 'samtools view -bS -o %s %s'

SAMTOOLS_SORT = 'samtools sort -@ %s %s %s'

#==========================================================================================================


def gmap_map(reference_database, reads, threads, out_format, min_intron_length, max_intron_length, exon_length, working_dir, Fflag, type_out, verbose):
    '''Calls gmap to map reads to reference
    Out_format can be samse of gff3 (2)'''

    if out_format == 'samse':
        filename = working_dir + 'gmap.long_reads.sam'
    elif out_format == '2' or out_format == 'gff3_gene':
        if type_out == 'cons':
            filename = working_dir + 'gmap.cluster_consensus.gff3'
        elif type_out == 'trin':
            filename = working_dir + 'gmap.trinity.gff3'
        elif type_out == 'ext':
            filename = working_dir + 'external.gff3'

    else:
        raise NameError(
            'Unknown format: ' + out_format + 'for GMAP. Accepted are samse or 2 (gff3_gene)')
    if os.path.isfile(filename) and os.path.getsize(filename) > 1:  # If the ref is there do not build it again
        sys.stdout.write(('GMAP index existed already: ' + filename + ' --- skipping'))
    else:
        out_f = open(filename, 'w')
        log_name = working_dir + 'gmap_map.log'
        log = open(log_name, 'w')
        if not Fflag:
            cmd = GMAP_SAM % (working_dir, reference_database, exon_length, min_intron_length, max_intron_length, threads, out_format, reads)
            try:
                if verbose:
                    sys.stderr.write('Executing: %s\n\n' % cmd)
                gmapmap = subprocess.Popen(cmd, stdout=out_f, stderr=log, shell=True)
                gmapmap.communicate()
            except:
                raise NameError('')
        else:
            cmd = GMAP_GFF % (working_dir, reference_database, exon_length, min_intron_length, max_intron_length, threads, out_format, reads)
            try:
                if verbose:
                    sys.stderr.write('Executing: %s\n\n' % cmd)
                gmapmap = subprocess.Popen(cmd, stdout=out_f, stderr=log, shell=True)
                gmapmap.communicate()
            except:
                raise NameError('')
        out_f.close()
        log.close()
    return filename


def star_build(reference, genome_dir, threads, wd, verbose):
    #TODO STOP if the build of the index fails
    '''Builds star reference index'''
    check_file = genome_dir + 'SAindex'
    if os.path.isfile(check_file):  # If the ref is there do not build it again
        sys.stdout.write(('STAR index existed already: ' +
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
            star_build = subprocess.Popen(cmd, stdout=log, stderr=log_err, shell=True, cwd=wd)
            star_build.communicate()
        except:
            sys.stdout.write('STAR build failed')
            raise NameError('')

        log.close()
        log_err.close()

    return None


def star_map(reads, threads, genome_dir, max_intron_length, wd, verbose):
    '''
    mapping short reads using STAR
    '''
    #TODO STOP if the mapping gices an empty file
    prefix = 'STAR_shortreads'
    if isinstance(reads, str):  # Only one file
        cmd = STAR_SINGLE % (threads, genome_dir, max_intron_length, max_intron_length, prefix, reads) # '--limitGenomeGenerateRAM', '1100000000',
    elif isinstance(reads, list):  # Paired end reads
        cmd = STAR_PAIRED % (threads, genome_dir, max_intron_length, max_intron_length, prefix, reads[0], reads[1])
    filename = wd + prefix + 'Aligned.out.bam'

    if os.path.isfile(filename):
        sys.stdout.write((
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
        star = subprocess.Popen(cmd, stdout=log, stderr=log,  shell=True, cwd=wd)
        star.communicate()
        # sys.stdout.write '>STAR worked. Output is: ' + filename +'\n'
    except:
        raise NameError('')


    log.close()
    log_err.close()

    if os.path.exists(filename) and os.path.getsize(filename) > 0:
        return filename
    else:
        sys.exit("##### STAR DID NOT WORK. PLEASE, CHECK THE FASTQ FILES #####\n")


def star(reference, fastq_reads, threads, max_intron_length, wd, verbose):
    '''Calls the mapper STAR to map fastq_reads to reference.
    First builds the reference index, then maps'''
    # Create dir for genome index
    refer = reference.split('/')[-1]
    genome_dir = wd + refer + '_STARindex/'
    logistic.check_create_dir(genome_dir)
    # Build the reference
    sys.stdout.write('\t###BUILD INDEX###\n')
    star_build(reference, genome_dir, threads, wd, verbose)

    # Mapping
    sys.stdout.write('\t###MAP###\n')
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
    """

    if '/' in reference:
        refer = reference.split('/')[-1] + '_GMAPindex'
    else:
        refer = reference + '_GMAPindex'

    cmd = GMAP_BULD % (refer, working_dir, reference)
    # If the ref is there do not build it again
    refer_path = working_dir + refer
    if os.path.isdir(refer_path):
        sys.stdout.write(('GMAP database existed already: ' + refer_path + ' --- skipping'))
        return refer
    log_name = working_dir + 'gmap_build.log'
    log = open(log_name, 'w')

    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        gmap_build = subprocess.Popen(cmd, stdout=log, stderr=log, shell =1)
        gmap_build.communicate()
    except:
        sys.stdout.write('GMAP build failed')
        raise NameError('')
    log.close()
    return refer


def gmap(type_out, reference, fastq_reads, threads, out_format, min_intron_length, max_intron_length, exon_length, wd, verbose, Fflag):
    '''Calls the mapper GMAP to map fastq_reads to reference.
    First builds the gmap database index with gmap_build(),
    then uses gmap_map() to map'''
    # Build the reference
    fmtdate = '%H:%M:%S %d-%m'
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write('\n###GMAP MAPPING  STARTED AT:\t' + now + '\t###\n')
    sys.stdout.write('\t###BUILD INDEX###\n')
    reference_db = gmap_build(reference, wd, verbose)
    # Mapping
    sys.stdout.write('\t###MAP###\n')
    out_file = gmap_map(reference_db, fastq_reads, threads, out_format, min_intron_length, max_intron_length, exon_length,
                        wd, Fflag, type_out, verbose)
    return out_file


def samtools_view(sam_file, wd, verbose):
    '''SAM to BAM'''
    bam_filename = sam_file + '.bam'
    cmd = SAMTOOLS_VIEW % (bam_filename, sam_file)
    if os.path.isfile(bam_filename):
        sys.stdout.write(('BAM file existed already: ' + bam_filename + ' --- skipping\n'))
        return bam_filename
    log_name = wd + 'samtools_view.log'
    log = open(log_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        samtools = subprocess.Popen(cmd, stderr=log, shell=True)
        samtools.communicate()
    except:
        raise NameError('')

    log.close()
    return bam_filename


def samtools_sort(bam_file, threads, wd, verbose):
    '''
    run a sorting of a bam file
    '''
    s_bam_filename = bam_file + '.sorted'

    if not os.path.isfile(s_bam_filename + '.bam'):
        cmd = SAMTOOLS_SORT % (threads, bam_file, s_bam_filename)
        s_bam_filename = s_bam_filename
        if os.path.isfile(s_bam_filename):
            sys.stdout.write(('Sorted BAM file existed already: ' + s_bam_filename + '.bam --- skipping\n'))
            return s_bam_filename
        log_name = wd + 'samtools_sort.log'
        log = open(log_name, 'w')
        try:
            if verbose:
                sys.stderr.write('Executing: %s\n' % cmd)
            samtools = subprocess.Popen(cmd, stderr=log, shell=True)
            samtools.communicate()
        except:
            raise NameError('')
        log.close()
    sor_bam_filename = s_bam_filename + ".bam"
    return sor_bam_filename


def sam_to_sorted_bam(sam_file, threads, wd, verbose):
    sys.stdout.write('\t###SAM to BAM###\n')
    bam_filename = samtools_view(sam_file, wd, verbose)

    sys.stdout.write('\t###SORTING BAM###\n')
    s_bam_filename = samtools_sort(bam_filename, threads, wd, verbose)
    return s_bam_filename
