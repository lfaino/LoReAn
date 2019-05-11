#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
from multiprocessing import Pool

import dirsAndFiles as logistic
import math

#==========================================================================================================
# COMMANDS LIST

TRINITY = 'Trinity --genome_guided_bam %s --genome_guided_max_intron %s --grid_node_max_memory %sG --max_memory %sG --output %s --CPU %s --full_cleanup --no_path_merging'
#TODO more option can be changed in trinity run
#TODO test if trinity can speed up

GMES_FU = 'gmes_petap.pl --ES --fungus --core %s --sequence %s'

GMES = 'gmes_petap.pl --ES --core %s --sequence %s'

BRAKER_FU = 'braker.pl --cores=%s  --useexisting --species=%s --workingdir=%s --genome=%s --fungus --bam=%s'

BRAKER = 'braker.pl --cores=%s  --useexisting --species=%s --workingdir=%s --genome=%s --bam=%s'

SAMTOOLS_VIEW = 'samtools view -@ %s -b -o %s %s %s'

SAMTOOLS_INDEX = 'samtools index %s'

#==========================================================================================================


def trinity(bam_file, wd, max_intron_length, threads, list_fasta_names, verbose):
    list_command = []
    threads_single_trinity = int(math.modf(int(threads)/len(list_fasta_names))[1])

    if threads_single_trinity > 1:
        pool_threads = len(list_fasta_names)
    else:
        pool_threads = 1
        threads_single_trinity = 1

    cmdI = SAMTOOLS_INDEX % (bam_file)
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmdI)
        samtools_call = subprocess.Popen(cmdI, shell=True, cwd=wd)
        samtools_call.communicate()
    except:
        sys.stdout.write('SAMTOOLS split did not work properly\n')
    for line in list_fasta_names:
        chrPath = os.path.join(wd, "trinity_out_dir_" + line.split("/")[-1])
        logistic.check_create_dir(chrPath)
        chrName = line.split("/")[-1].split(".")[0]
        bamName = os.path.join(chrPath, chrName + ".bam")
        cmdS = SAMTOOLS_VIEW % (threads, bamName, bam_file, chrName)
        try:
            if verbose:
                sys.stderr.write('Executing: %s\n' % cmdS)
            samtools_call = subprocess.Popen(cmdS, shell=True, cwd=wd)
            samtools_call.communicate()
            if os.path.getsize(bamName) > 0:
                list_command.append([bamName, chrPath, max_intron_length, threads_single_trinity, verbose])
        except:
            sys.stdout.write('SAMTOOLS split did not work properly\n')
            raise NameError('')

    with Pool(processes=pool_threads, maxtasksperchild=1000) as pool:
        results = pool.map(trinity_single, list_command)
    out_name = os.path.join(wd, 'Trinity-GG.fasta')
    with open(out_name,'wb') as wfd:
        for f in results:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
    return out_name, results



def trinity_single(run):
    """Calls genome guided trinity on the BAM file to generate
    assembled transcripts"""
    bam_file, out_dir, max_intron_length, threads, verbose = run
    #out_dir = wd + 'trinity_out_dir/'
    cmd = TRINITY % (bam_file, max_intron_length, '3','10', out_dir, threads)
    out_name = os.path.join(out_dir, 'Trinity-GG.fasta')
    if os.path.isfile(out_name):
        sys.stdout.write(('Trinity-GG file existed already: ' + out_name + ' --- skipping\n'))
        return out_name
    log_name_err = out_dir + 'trinity.err.log'
    log_err = open(log_name_err, 'w')
    log_name = out_dir + 'trinity.log'
    log = open(log_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        trinity_call = subprocess.Popen(cmd, stdout=log, stderr=log_err, shell=True, cwd = out_dir)
        trinity_call.communicate()
    except:
        sys.stdout.write('Trinity did not work properly\n')
        raise NameError('')
    log_err.close()
    log.close()
    return out_name


def braker_call(wd, reference, bam_file, species_name, threads, fungus, verbose):
    '''Calls braker, may take a while'''

    sys.stdout.write ("###RUNNING BRAKER1 ###\n")

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
        braker_ex = subprocess.Popen(cmd, cwd=wd, stdout=log, stderr=log_err, shell=True)
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
        sys.stdout.write("###RUNNING GENEMARK ###\n")
        try:
            if verbose:
                sys.stderr.write('Executing: %s\n' % cmd)
            genemarks = subprocess.Popen(cmd, stderr=log_e, stdout=log, cwd=wd, shell=True)
            genemarks.communicate()
        except:
            raise NameError('')
        log.close()
        log_e.close()

    return wd


def find_species(home):
    augustus_species_cmd = 'augustus --species=help'
    process = subprocess.Popen(augustus_species_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out_augustus, err_augustus = process.communicate()
    list_file = [os.path.join(home, folder) for folder in os.listdir(home) if os.path.isfile(os.path.join(home, folder))
                 and ".bashrc" == folder]
    with open(list_file[0]) as bashrc:
        for path in bashrc:
            if "AUGUSTUS_CONFIG_PATH" in path:
                augustus_specie_dir = path.split("=~")[1].rsplit()[0]
                augustus_species = [species for species in os.listdir(home + augustus_specie_dir + "species")]

    return augustus_species, err_augustus
