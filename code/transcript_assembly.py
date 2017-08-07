#!/usr/bin/env python3
import os
import subprocess
import sys



#==========================================================================================================
# COMMANDS LIST

TRINITY = 'Trinity --genome_guided_bam %s --genome_guided_max_intron %s --max_memory 5G --output %s --CPU %s --full_cleanup'

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
        sys.stdout.write((
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
        trinity_call = subprocess.Popen(cmd, stdout=log, stderr=log_err, shell=1, cwd = wd)
        trinity_call.communicate()
    except:
        sys.stdout.write('Trinity did not work properly\n')
        raise NameError('')
    log_err.close()
    log.close()
    return out_name

def braker_call(wd, reference, bam_file, species_name, threads, fungus, verbose):
    '''Calls braker, may take a while'''
    # perl ~/bin/BRAKER1/braker.pl --cores=3 --workingdir=/home/jose/mapper_testing/braker1_output/ --species=gmap_gff3
    #--genome=/home/jose/Reference/JR2_Chr8/Verticillium_dahliaejr2.GCA_000400815.2.29.dna.chromosome.8.fa
    #--bam=/home/jose/mapper_testing/gmap/gmap_Chr8_2Dall.sorted.bam
    sys.stdout.write ("\n###RUNNING BRAKER1 ###\n")

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


    if int(threads) < 1:
        threads = 1

    if fungus:
        cmd = GMES_FU % (threads, ref)
    else:
        cmd = GMES % (threads, ref)

    log_name = wd + 'gm_es.gff'

    if os.path.exists(log_name) and os.path.getsize(log_name) > 0:
        sys.stderr.write('Already executed: %s\n' % cmd)
        pass
    else:
        log = open(log_name, 'w')
        log_name_err = wd + 'gm_es.err.log'
        log_e = open(log_name_err, 'w')
        sys.stdout.write("\n###RUNNING GENEMARK ###\n")
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
