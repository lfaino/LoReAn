#!/usr/bin/env python3
import os
import subprocess
import sys



#==========================================================================================================
# COMMANDS LIST

TRINITY = 'Trinity --genome_guided_bam %s --genome_guided_max_intron %s --max_memory 10G --output %s --CPU %s --full_cleanup'

LAUNCH_PASA =  'Launch_PASA_pipeline.pl -c %s -C -r -R -g %s -t %s --ALIGNERS gmap --TRANSDECODER -I %s --CPU %s'

GMES_FU = 'gmes_petap.pl --ES --fungus --core %s --sequence %s'

GMES = 'gmes_petap.pl --ES --core %s --sequence %s'

BRAKER_FU = 'braker.pl --cores=%s  --useexisting --species=%s --workingdir=%s --genome=%s --fungus --bam=%s'

BRAKER = 'braker.pl --cores=%s  --useexisting --species=%s --workingdir=%s --genome=%s --bam=%s'

#==========================================================================================================


def trinity(bam_file, wd, max_intron_length, threads, verbose):
    """Calls genome guided trinity on the BAM file to generate
    assembled transcripts"""
    out_dir = wd + 'trinity_out_dir/'
    cmd = TRINITY % (bam_file, max_intron_length,  out_dir, threads)
    out_name = out_dir + 'Trinity-GG.fasta'
    if os.path.isfile(out_name):
        print((
            'Trinity-GG file existed already: ' +
            out_name +
            ' --- skipping\n'))
        return out_name
    log_name_err = wd + 'trinity.err.log'
    log_err = open(log_name_err, 'w')
    log_name = wd + 'trinity.log'
    log = open(log_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        trinity_call = subprocess.Popen(cmd, stdout=log, stderr=log_err, shell=1)
        trinity_call.communicate()
    except:
        print('Trinity did not work properly\n')
        raise NameError('')
    log_err.close()
    log.close()
    return out_name

def pasa_configuration(pasa_dir, pasa_db):
    '''Creates a PASA configuration file. Database name will be the reference name'''
    conf_file = pasa_dir + 'alignAssembly.config'
    if os.path.isfile(conf_file):
        print((
            'PASA configuration file existed already: ' +
            conf_file +
            ' --- skipping\n'))
        return conf_file
    conf = open(conf_file, 'w')
    lines = [
        'MYSQLDB=' + pasa_db,
        'validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>',
        'validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>',
        'subcluster_builder.dbi:-m=50']
    for line in lines:
        conf.write(line + '\n')
    conf.close()
    return conf_file

def pasa_call(pasa_dir, conf_file, pasa_db, reference, transcripts, max_intron_length, threads, verbose):
    '''PASA to construct a database of transcripts. It will overwrite any
    database with the same name -the one of the reference-.'''
    cmd = LAUNCH_PASA % (conf_file, reference, transcripts, max_intron_length, threads)
    out_file = pasa_dir + pasa_db + '.pasa_assemblies.gff3'
    # print out_file, os.path.isfile(out_file)
    if os.path.isfile(out_file):
        print(('PASA output existed already: ' + out_file + ' --- skipping\n'))
        return out_file
    log_name = pasa_dir + 'pasa.err.log'
    log_out_name = pasa_dir + 'pasa.out.log'
    log = open(log_name, 'w')
    out_log = open(log_out_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        pasa = subprocess.Popen(cmd, stdout=out_log, stderr=log, cwd=pasa_dir, shell=1)
        pasa.communicate()
    except:
        print('PASA failed')
        raise NameError('')
    log.close()
    out_log.close()
    return out_file

def braker_call(wd, reference, bam_file, species_name, threads, fungus, verbose):
    '''Calls braker, may take a while'''
    # perl ~/bin/BRAKER1/braker.pl --cores=3 --workingdir=/home/jose/mapper_testing/braker1_output/ --species=gmap_gff3
    #--genome=/home/jose/Reference/JR2_Chr8/Verticillium_dahliaejr2.GCA_000400815.2.29.dna.chromosome.8.fa
    #--bam=/home/jose/mapper_testing/gmap/gmap_Chr8_2Dall.sorted.bam
    print ("\n###RUNNING BRAKER1 ###\n")

    if fungus:
        cmd = BRAKER_FU % (threads, species_name, wd, reference, bam_file)
    else:
        cmd = BRAKER % (threads, species_name, wd, reference, bam_file)

    log_name = wd + 'braker.log'
    log = open(log_name, 'w')
    log_name_err = wd + 'braker.error.log'
    log_err = open(log_name_err, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        braker_ex = subprocess.Popen(cmd, stdout=log, stderr=log_err, shell=1)
        braker_ex.communicate()
    except:
        raise NameError('')
    log.close()
    log_err.close()
    return

def gmes_call(wd, ref, fungus, threads, verbose):
    log_name = wd + 'gm_es.gff'
    log = open(log_name, 'w')
    log_name_err = wd + 'gm_es.err.log'
    log_e = open(log_name_err, 'w')
    print ("\n###RUNNING GENEMARK ###\n")

    if fungus:
        cmd = GMES_FU % (threads, ref)
    else:
        cmd = GMES % (threads, ref)

    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        genemarks = subprocess.Popen(cmd, stderr=log_e, stdout=log, cwd=wd, shell=1)
        genemarks.communicate()
    except:
        raise NameError('')
    log.close()
    log_e.close()

    return wd
