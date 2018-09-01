#!/usr/bin/env python3

import os
import re
import shutil
import subprocess
import sys

#==========================================================================================================
# COMMANDS LIST

LAUNCH_PASA =  'Launch_PASA_pipeline.pl -c %s -C -r -R -g %s -t %s --ALIGNERS gmap --TRANSDECODER -I %s --CPU %s'

LOAD_ANOT = 'Load_Current_Gene_Annotations.dbi -c %s -g %s -P %s'

COMP_ANNOT = 'Launch_PASA_pipeline.pl --ALT_SPLICE --TRANSDECODER --CPU %s -c %s -A -g %s -t %s'

#COMP_ANNOT_CREATE = 'Launch_PASA_pipeline.pl --ALT_SPLICE --TRANSDECODER --CPU %s -C -c %s -A -g %s -t %s'

LOCATION_CONF_ORIGINAL = "/opt/LoReAn/third_party/software/PASApipeline/pasa_conf/conf.txt"

LOCATION_CONF_NEW = "/opt/LoReAn/third_party/software/PASApipeline/pasa_conf/conf.txt"

CREATE_DATABASE = 'create_mysql_cdnaassembly_db.dbi -c %s -S /opt/LoReAn/third_party/software/PASApipeline/schema/cdna_alignment_mysqlschema'

#==========================================================================================================


def pasa_annot_configuration(pasa_dir, pasa_db):
    '''Creates a PASA annotation configuration file'''
    conf_file = pasa_dir + 'annotCompare.config'
    conf = open(conf_file, 'w')
    lines = [
        'DATABASE=' + pasa_db,
        'cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=<__MIN_PERCENT_OVERLAP__>',
        'cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=<__MIN_PERCENT_PROT_CODING__>',
        'cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=<__MIN_PERID_PROT_COMPARE__>',
        'cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=<__MIN_PERCENT_LENGTH_FL_COMPARE__>',
        'cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=<__MIN_PERCENT_LENGTH_NONFL_COMPARE__>',
        'cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE=<__MIN_FL_ORF_SIZE__>',
        'cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=<__MIN_PERCENT_ALIGN_LENGTH__>',
        'cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=<__MIN_PERCENT_OVERLAP_GENE_REPLACE__>',
        'cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE=<__STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE__>',
        'cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS=<__TRUST_FL_STATUS__>',
        'cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS=<__MAX_UTR_EXONS__>',
        'cDNA_annotation_comparer.dbi:--GENETIC_CODE=<__GENETIC_CODE__>']
    for line in lines:
        conf.write(line + '\n')

    conf.close()

    return conf_file

def pasa_mysql_configuration(pasa_data):
    '''Creates a PASA conf file'''

    user_name, user_pws = pasa_data.split(",")
    with open(LOCATION_CONF_ORIGINAL, "r") as cfh:
        with open(LOCATION_CONF_NEW, "w") as cnfh:
            for line in cfh:
                if line.startswith("MYSQL_RW_USER") or line.startswith("MYSQL_RO_USER"):
                    cnfh.write(("=".join([line.split("=")[0], user_name])) + "\n")
                elif line.startswith("MYSQL_RW_PASSWORD") or line.startswith("MYSQL_RO_PASSWORD"):
                    cnfh.write(("=".join([line.split("=")[0], user_pws])) + "\n")
                else:
                    cnfh.write(line)
    sys.stdout.write('\n### PASA WILL RUN USING MYSQL  HAVING %s AS USER AND %s AS USER PASSWORD ###\n' % (user_name, user_pws))



def load_gff3_pasa(pasa_dir, align_conf_file, reference, gff3_file, verbose):
    '''Loads a gff3 file into a PASA database '''
    cmd = LOAD_ANOT % (align_conf_file, reference, gff3_file)
    log_name = pasa_dir + 'load_gff3.log'
    stdout_file = pasa_dir + 'load_gff3.stdout'
    log = open(log_name, 'w')
    stdout_f = open(stdout_file, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        load = subprocess.Popen(cmd, stderr=log, stdout=stdout_f, cwd=pasa_dir, shell=True)
        load.communicate()
        processID = load.pid
        # sys.stdout.write '> GFF3 loaded to PASA DB \n'
    except:
        # sys.stdout.write 'Could not load GFF3 to PASA DB\n'
        raise NameError('')

    log.close()
    stdout_f.close()
    return processID


def annot_comparison(pasa_dir, annot_conf_file, reference, transcripts_file, n_cpu, round_n, verbose):
    '''Loads a gff3 file into a PASA database '''

    cmd = COMP_ANNOT % (n_cpu, annot_conf_file, reference, transcripts_file)
    log_name = pasa_dir + 'update_gff3.log'
    log = open(log_name, 'w')
    log_out_name = pasa_dir + 'pasa.out.log'
    out_log = open(log_out_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        pasa_call = subprocess.Popen(cmd, stdout=out_log, stderr=log, cwd=pasa_dir, shell=True)
        pasa_call.communicate()
    except:
        raise NameError('')
    out_log.close()
    log.close()

def parse_pasa_update(round_n, pasa_dir, pasa_db, verbose):
    '''Parses through the files in the PASA directory, finds the update file and
    renames it and returns it'''
    pasa_files = os.listdir(pasa_dir)

    pattern_build = '^' + pasa_db + '.gene_structures_post_PASA_updates.[0-9]+.gff3$'

    pasa_pattern = re.compile(pattern_build)
    for filename in pasa_files:
        match = re.match(pasa_pattern, filename)
        if match:
            update_file = filename

    new_filename = pasa_dir + 'FinalAnnotationLorean' + round_n + '.gff3'
    root = os.path.join(pasa_dir, update_file)
    shutil.move(root, new_filename)

    return new_filename


def update_database(n_cpu, round_n, pasa_dir, pasa_db, reference, transcripts_file, gff3_file, verbose):
    '''Updates the gff3 file with the PASA database'''
    sys.stdout.write('\t###CREATING CONFIGURATION FILE###\n')
    annot_conf_file = pasa_annot_configuration(pasa_dir, pasa_db)
    create_pasa_database(annot_conf_file, pasa_dir, verbose)
    align_conf_file = pasa_configuration(pasa_dir, pasa_db, verbose)
    sys.stdout.write('\t###LOADING GFF3 FILE INTO DATABASE###\n')
    load_gff3_pasa(pasa_dir, align_conf_file, reference, gff3_file, verbose)
    sys.stdout.write('\t###UPDATING GFF3 FILE###\n')
    annot_comparison(pasa_dir, annot_conf_file, reference, transcripts_file, n_cpu, round_n, verbose)
    sys.stdout.write('\t###PARSING OUTPUT###\n')
    gff3_out = parse_pasa_update(round_n, pasa_dir, pasa_db, verbose)

    return gff3_out


def pasa_configuration(pasa_dir, pasa_db, verbose):
    '''Creates a PASA configuration file. Database name will be the reference name'''
    conf_file = pasa_dir + 'alignAssembly.config'
    if os.path.isfile(conf_file):
        sys.stdout.write((
            'PASA configuration file existed already: ' +
            conf_file +
            ' --- skipping\n'))
        return conf_file
    conf = open(conf_file, 'w')
    lines = [
        'DATABASE=' + pasa_db,
        'validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>',
        'validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>',
        'subcluster_builder.dbi:-m=50']
    for line in lines:
        conf.write(line + '\n')
    conf.close()
    return conf_file


def pasa_call(pasa_dir, pasa_db, reference, transcripts, max_intron_length, threads, verbose):
    '''PASA to construct a database of transcripts. It will overwrite any
    database with the same name -the one of the reference-.'''

    align_pasa_conf = pasa_configuration(pasa_dir, pasa_db, verbose)
    cmd = LAUNCH_PASA % (align_pasa_conf, reference, transcripts, max_intron_length, threads)
    out_file = pasa_dir + pasa_db + '.pasa_assemblies.gff3'
    # sys.stdout.write out_file, os.path.isfile(out_file)
    if os.path.isfile(out_file):
        sys.stdout.write(('PASA output existed already: ' + out_file + ' --- skipping\n'))
        return out_file
    log_name = pasa_dir + 'pasa.err.log'
    log_out_name = pasa_dir + 'pasa.out.log'
    log = open(log_name, 'w')
    out_log = open(log_out_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        pasa = subprocess.Popen(cmd, stdout=out_log, stderr=log, cwd=pasa_dir, shell=True)
        pasa.communicate()
    except:
        sys.stdout.write('PASA failed')
        raise NameError('')
    log.close()
    out_log.close()
    return out_file


def create_pasa_database(conf_file, pasa_dir, verbose):

    cmd = CREATE_DATABASE % conf_file
    log_name = pasa_dir + 'create_databae.log'
    log = open(log_name, 'w')
    log_out_name = pasa_dir + 'create_databae.out.log'
    out_log = open(log_out_name, 'w')
    try:
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        pasa_call = subprocess.Popen(cmd, stdout=out_log, stderr=log, cwd=pasa_dir, shell=True)
        pasa_call.communicate()
    except:
        raise NameError('')