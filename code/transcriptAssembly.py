#!/usr/bin/env python3

import os
import subprocess
import sys
from collections import namedtuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from simplesam import Reader

#==========================================================================================================
# COMMANDS LIST

TRINITY = 'Trinity --genome_guided_bam %s --genome_guided_max_intron %s --grid_node_max_memory %sG --max_memory %sG --output %s --CPU %s --full_cleanup --no_path_merging --max_reads_per_graph 100000'
#TODO more option can be changed in trinity run
#TODO test if trinity can speed up

GMES_FU = 'gmes_petap.pl --ES --fungus --core %s --sequence %s'

GMES = 'gmes_petap.pl --ES --core %s --sequence %s'

BRAKER_FU = 'braker.pl --cores=%s  --useexisting --species=%s --workingdir=%s --genome=%s --fungus --bam=%s'

BRAKER = 'braker.pl --cores=%s  --useexisting --species=%s --workingdir=%s --genome=%s --bam=%s'

SAMTOOLS_VIEW = 'samtools view -@ %s -b -o %s %s %s'

SAMTOOLS_INDEX = 'samtools index %s'

BAMTOFASTQ = 'bedtools bamtofastq -i %s -fq /dev/stdout'

#==========================================================================================================


def bamtofastq(bam, verbose):
    fasta = bam + ".fasta"
    in_file = open(bam, 'r')
    in_sam = Reader(in_file)
    with open(fasta, "w") as output_handle:
        for line in in_sam:
            if line.mapped:
                record = SeqRecord(Seq(str(line.seq)),name = str(line.qname))
                SeqIO.write(record, output_handle, "fasta")
    return fasta

def trinity(bam_file, wd, max_intron_length, threads, verbose):
    """Calls genome guided trinity on the BAM file to generate
    assembled transcripts"""
    #bam_file, out_dir, max_intron_length, threads, verbose = run


    MemInfoEntry = namedtuple('MemInfoEntry', ['value', 'unit'])
    meminfo = {}
    with open('/proc/meminfo') as file:
        for line in file:
            key, value, *unit = line.strip().split()
            meminfo[key.rstrip(':')] = MemInfoEntry(value, unit)
    memtotal = round((int(meminfo['MemTotal'].value)/10000000)-2)

    if memtotal < threads:
        threads = str(memtotal)

    out_dir = wd + 'trinity_out_dir/'
    cmd = TRINITY % (bam_file, max_intron_length, '3','10', out_dir, threads)
    out_name = os.path.join(out_dir, 'Trinity-GG.fasta')
    if os.path.isfile(out_name):
        sys.stdout.write(('Trinity-GG file existed already: ' + out_name + ' --- skipping\n'))
        return out_name
    log_name_err = wd + 'trinity.err.log'
    log_err = open(log_name_err, 'w')
    log_name = wd + 'trinity.log'
    log = open(log_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        trinity_call = subprocess.Popen(cmd, stdout=log, stderr=log_err, shell=True, cwd=wd)
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
