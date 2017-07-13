#!/usr/bin/env python3

import os
import math
import subprocess
import dirs_and_files as logistic
from Bio import SeqIO


def star_build(reference, genome_dir, threads, wd):
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
            args = [
                'STAR',
                '--runThreadN',
                str(threads),
                '--runMode',
                'genomeGenerate',
                '--genomeChrBinNbits',
                '=',
                'min(',
                str(fragmented_genome),
                ')', '--genomeSAindexNbases', '6',
                '--genomeDir',
                genome_dir,
                '--genomeFastaFiles',
                reference]
        else:  # Genome small enough
            args = [
                'STAR',
                '--runThreadN',
                str(threads),
                '--runMode',
                'genomeGenerate',
                '--genomeDir',
                genome_dir, '--genomeSAindexNbases', '6',
                '--genomeFastaFiles',
                reference]

        log_name_err = wd + 'star_build.err.log'
        log_err = open(log_name_err, 'w')
        log_name = wd + 'star_build.log'
        log = open(log_name, 'w')

        # Call
        try:
            subprocess.check_call(args, stdout=log, stderr=log_err, cwd=wd)

        except:
            print('STAR build failed')
            raise NameError('')

        log.close()
        log_err.close()

    return None


def star_map(reference, reads, threads, genome_dir, max_intron_length, wd):
    '''Calls STAR to map reads to reference'''

    global args
    prefix = 'STAR_shortreads'
    if isinstance(reads, str):  # Only one file
        args = [
            'STAR',
            '--runThreadN',
            threads,
            '--genomeDir',
            genome_dir,
            '--outSAMtype',
            'BAM',
            'Unsorted',
            '--alignIntronMax',
            max_intron_length,
            '--alignMatesGapMax',
            '1000',
            '--outFilterMismatchNmax',
            '15',
            '--outFileNamePrefix',
            prefix,
            '--outSAMstrandField',
            'intronMotif',
            '--outFilterIntronMotifs',
            'RemoveNoncanonical',
            '--outSAMattrIHstart',
            '1',
            '--outFilterMultimapNmax',
            '3',
            '--readFilesIn',
            reads]  # '--limitGenomeGenerateRAM', '1100000000',
    elif isinstance(reads, list):  # Paired end reads
        args = [
            'STAR',
            '--runThreadN',
            threads,
            '--genomeDir',
            genome_dir,
            '--outSAMtype',
            'BAM',
            'Unsorted',
            '--alignIntronMax',
            max_intron_length,
            '--alignMatesGapMax',
            '1000',
            '--outFilterMismatchNmax',
            '15',
            '--outFileNamePrefix',
            prefix,
            '--outSAMstrandField',
            'intronMotif',
            '--outFilterIntronMotifs',
            'RemoveNoncanonical',
            '--outSAMattrIHstart',
            '1',
            '--outFilterMultimapNmax',
            '3',
            '--readFilesIn',
            reads[0],
            reads[1]]

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
        subprocess.check_call(args, stdout=log, stderr=log, cwd=wd)
        # print '>STAR worked. Output is: ' + filename +'\n'
    except:
        # print 'STAR did not work properly\n'
        raise NameError('')

    log.close()
    log_err.close()

    return filename


def star(reference, fastq_reads, threads, max_intron_length, wd):
    '''Calls the mapper STAR to map fastq_reads to reference.
    First builds the reference index, then maps'''
    # Create dir for genome index
    refer = reference.split('/')[-1]
    genome_dir = wd + refer + '_STARindex/'
    logistic.check_create_dir(genome_dir)
    # Build the reference
    print('\t###BUILD INDEX###\n')
    star_build(reference, genome_dir, threads, wd)

    # Mapping
    print('\t###MAP###\n')
    out_file = star_map(
        reference,
        fastq_reads,
        threads,
        genome_dir,
        max_intron_length,
        wd)
    return out_file

def gmap_build(reference, working_dir):
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

    args = ['gmap_build', '-k' , '13' , '-d', refer, '-D',
            working_dir, reference]
    # If the ref is there do not build it again
    refer_path = working_dir + refer

    if os.path.isdir(refer_path):
        print(('GMAP database existed already: ' +
               refer_path + ' --- skipping'))
        return refer

    log_name = working_dir + 'gmap_build.log'
    log = open(log_name, 'w')

    # Call
    try:
        subprocess.check_call(args, stdout=log, stderr=log)
    except:
        print('GMAP build failed')
        raise NameError('')

    log.close()
    return refer


def gmap(
        type_out,
        reference,
        fastq_reads,
        threads,
        out_format,
        min_intron_length,
        max_intron_length,
        exon_length,
        wd,
        Fflag):
    '''Calls the mapper GMAP to map fastq_reads to reference.
    First builds the gmap database index with gmap_build(),
    then uses gmap_map() to map'''
    # Build the reference
    print('\t###BUILD INDEX###\n')
    reference_db = gmap_build(reference, wd)
    # Mapping
    print('\t###MAP###\n')
    if Fflag:
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
            type_out)
    else:
        out_file = alignLong(reference_db, fastq_reads, threads, max_intron_length, reference, wd)
    return out_file


def alignLong(reference_db, fastq_reads, threads, max_intron_length, reference_seq, wd):
    """
    :param reference_db:
    :param fastq_reads:
    :param threads:
    :param max_intron_length:
    :param wd:
    :return:
    """

    args = [ 'alignPacBio.py', '-p', str(threads), '-K', str(max_intron_length), wd, reference_db, reference_seq, fastq_reads]
    map = subprocess.Popen(args, cwd = wd)
    map.communicate()
    output_file = wd + "/filtered/aligned.bam"
    return (output_file)


def gmap_map(
        reference_database,
        reads,
        threads,
        out_format,
        min_intron_length,
        max_intron_length,
        exon_length,
        working_dir,
        Fflag,
        type_out):
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
    if os.path.isfile(filename) and os.path.getsize(
            filename) > 1:  # If the ref is there do not build it again
        print(('STAR index existed already: ' +
               filename + ' --- skipping'))
    else:
        out_f = open(filename, 'w')
        log_name = working_dir + 'gmap_map.log'
        log = open(log_name, 'w')
        args = [ 'gmap',  '-D', str(working_dir),  '-d', str(reference_database), '--trim-end-exons',  str(exon_length),  '--cross-species',  '--expand-offsets', '1',  '-B',  '5',   '--min-intronlength',
str(min_intron_length),  '-n',  '1',  '--microexon-spliceprob', '1', '--max-intronlength-middle',  str(max_intron_length), '--max-intronlength-ends', str(max_intron_length), '-f', str(out_format),  '-t', str(threads)]
        if not Fflag:
            argsready = args + [reads]
            try:
                gmapmap = subprocess.Popen(argsready, stdout=out_f, stderr=log)
                gmapmap.communicate()
                # print '>GMAP worked. Output is: ' + filename +'\n'
            except:
                # print 'GMAP did not work properly\n'
                raise NameError('')
        else:
            argsready = args + ['F', reads]
            try:
                gmapmap = subprocess.Popen(argsready, stdout=out_f, stderr=log)
                gmapmap.communicate()
                # print '>GMAP worked. Output is: ' + filename +'\n'
            except:
                # print 'GMAP did not work properly\n'
                raise NameError('')
        out_f.close()
        log.close()
    return filename

def samtools_view(sam_file, wd):
    '''SAM to BAM'''
    bam_filename = sam_file + '.bam'
    args = ['samtools', 'view', '-bS', '-o', bam_filename, sam_file]
    if os.path.isfile(bam_filename):
        print(('BAM file existed already: ' + bam_filename + ' --- skipping\n'))
        return bam_filename
    log_name = wd + 'samtools_view.log'
    log = open(log_name, 'w')
    try:
        subprocess.check_call(args, stderr=log)
        # print '> SAM converted to BAM in ' + bam_filename
    except:
        # print 'Samtools view failed'
        raise NameError('')

    log.close()
    return bam_filename


def samtools_sort(bam_file, threads, wd):
    '''
    run a sorting of a bam file
    :param bam_file:
    :param threads:
    :param wd:
    :return:
    '''
    s_bam_filename = bam_file + '.sorted'

    if not os.path.isfile(s_bam_filename + '.bam'):
        args = [
            'samtools',
            'sort',
            '-@',
            str(threads),
            bam_file,
            s_bam_filename]
        s_bam_filename = s_bam_filename
        if os.path.isfile(s_bam_filename):
            print(('Sorted BAM file existed already: ' + s_bam_filename +
                   '.bam --- skipping\n'))
            return s_bam_filename
        log_name = wd + 'samtools_sort.log'
        log = open(log_name, 'w')
        try:
            subprocess.check_call(args, stderr=log)

            # print '> BAM sorted in ' + s_bam_filename + '\n'
        except:
            # print 'Samtools sort failed'
            raise NameError('')
        log.close()
    sor_bam_filename = s_bam_filename + ".bam"
    return sor_bam_filename


def sam_to_sorted_bam(sam_file, threads, wd):
    print('\t###SAM to BAM###\n')
    bam_filename = samtools_view(sam_file, wd)

    print('\t###SORTING BAM###\n')
    s_bam_filename = samtools_sort(bam_filename, threads, wd)
    return s_bam_filename
