#!/usr/bin/env python

''' 
MASTER THESIS PROJECT
Author: Jose A. Espejo
Date: September 2015 - March 2016

EVM pipeline including PASA database update

'''
import os
import subprocess
from dirs_and_files import check_create_dir
import re

def evm_partitions(evm_output, reference, gene_preds, transcripts, proteins, segmentSize, overlapSize):
    '''Genome sequences and  gff3 files are partitioned based on individual contigs, 
    and large contigs are segmented into smaller overlapping chunks.'''
    
    if proteins == '' and transcripts != '':
        args = ['partition_EVM_inputs.pl', '--genome', reference,
            '--gene_predictions', gene_preds, '--transcript_alignments', 
            transcripts, '--segmentSize', segmentSize, '--overlapSize', overlapSize, 
            '--partition_listing', 'partitions_list.out']
    elif proteins == '' and transcripts == '':
        args = ['partition_EVM_inputs.pl', '--genome', reference,
            '--gene_predictions', gene_preds, '--segmentSize', segmentSize,    
            '--overlapSize', overlapSize, 
            '--partition_listing', 'partitions_list.out']
    elif transcripts == '' and proteins != '':
        args = ['partition_EVM_inputs.pl', '--genome', reference,
            '--gene_predictions', gene_preds, '--protein_alignments', proteins, 
            '--segmentSize', segmentSize, '--overlapSize', overlapSize, 
            '--partition_listing', 'partitions_list.out']
    else:
        args = ['partition_EVM_inputs.pl', '--genome', reference,
            '--gene_predictions', gene_preds, '--protein_alignments', proteins,
            '--transcript_alignments', transcripts, '--segmentSize', segmentSize, 
            '--overlapSize', overlapSize, '--partition_listing', 'partitions_list.out']
    
    partitions = evm_output + 'partitions_list.out'
    if os.path.isfile(partitions): 
        print ('\nPartitions file file existed already: ' + partitions + ' --- skipping\n')
        return partitions
    
    stdout_file = evm_output + 'partitions.stdout'
    log_name = evm_output + 'partitions.log'
    log = open(log_name, 'w')
    stdout_f = open(stdout_file, 'w')
    
    try:
        subprocess.check_call(args, stdout = stdout_f, stderr = log, cwd = evm_output)
        #print '> Partitions created\n'
    except:
        #print 'Partitions could not be created\n'
        raise NameError('')
    
   
    log.close()
    stdout_f.close()
    return partitions

def evm_write_commands(evm_output, reference, weights, gene_preds, transcripts, proteins, partitions):
    '''
    Writes a file with the necessary commnd for the partitions
    '''
    if proteins != '':
        args = ['write_EVM_commands.pl', '--genome', reference,
            '--weights', weights, '--gene_predictions', gene_preds, 
            '--transcript_alignments', transcripts, '--protein_alignments', 
            proteins,'--output_file_name',
            'evm.out', '--partitions', 'partitions_list.out']
    else:
        args = ['write_EVM_commands.pl', '--genome', reference,
            '--weights', weights, '--gene_predictions', gene_preds, 
            '--transcript_alignments', transcripts, '--output_file_name',
            'evm.out', '--partitions', 'partitions_list.out']
    
    command_file = evm_output + 'commands.list'
    if os.path.isfile(command_file): 
        print ('\nCommand file existed already: ' + command_file + ' --- skipping\n')
        return command_file
    
    log_name = evm_output + 'write_commands.log'
    command = open(command_file, 'w')
    log = open(log_name, 'w')
    
    try:
        subprocess.check_call(args, stdout = command, stderr = log, cwd = evm_output)
        #print '> Command list created. Output is: ' + command_file +'\n'
    except:
        #print 'Command file could not be created\n'
        raise NameError('')
    
    command.close()
    log.close()
    return command_file

def evm_run(evm_output, command_list, threads):
    '''
    Runs all the commands in commands.list in parallel
    '''
    args1 = ['cat', command_list]
    args2 = ['parallel' , '-j', str(threads), '--']
    
    
    out_file = evm_output + 'evm.out.combined.gff3' #THIS OUTPUT FROM THE WHOLE PIPELINE
    if os.path.isfile(out_file): 
        print ('\nEVM output existed already: ' + out_file + ' --- skipping\n')
        return ''
        
    log_name = evm_output + 'run.log'
    log = open(log_name, 'w')
    
    try:
        cat = subprocess.Popen(args1, stdout = subprocess.PIPE, cwd = evm_output)
        parallel = subprocess.check_call(args2, stdin = cat.stdout, cwd = evm_output, stderr = log)
        
        #print '> EVM finished.\n'
    except:
        #print 'Error in EVM run\n'
        raise NameError('')
    
    log.close()

def evm_combine(evm_output, partitions):
    '''
    perl /home/jose/bin/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl 
    --partitions partitions_list.out --output_file_name evm.out
    '''
    args = ['recombine_EVM_partial_outputs.pl', '--partitions', partitions,
            '--output_file_name', 'evm.out']
    
    out_file = evm_output + 'evm.out.combined.gff3' #THIS OUTPUT FROM THE WHOLE PIPELINE
    if os.path.isfile(out_file): 
        print ('\nEVM output existed already: ' + out_file + ' --- skipping\n')
        return ''
    st_file = evm_output + 'combine_partitions.stdout'
    log_name = evm_output + 'combine_partitions.log'
    st_out = open(st_file, 'w')
    log = open(log_name, 'w')
    

    try:
        subprocess.check_call(args, stdout = st_out, stderr = log, cwd = evm_output)
        #print '> Partitions combined\n'
    except:
        #print 'Partitions could not be combined\n'
        raise NameError('')
    
    st_out.close()
    log.close()

def evm_to_gff3(evm_output, partitions, reference):
    '''
    perl /home/jose/bin/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl 
    --partitions partitions_list.out --output evm.out 
    --genome ~/Reference/JR2_Chr8/Verticillium_dahliaejr2.GCA_000400815.2.29.dna.chromosome.8.fa
    '''
    args = ['convert_EVM_outputs_to_GFF3.pl', '--partitions', partitions,
            '--output', 'evm.out', '--genome', reference]
    
    out_file = evm_output + 'evm.out.combined.gff3' #THIS OUTPUT FROM THE WHOLE PIPELINE
    if os.path.isfile(out_file): 
        print ('\nEVM output existed already: ' + out_file + ' --- skipping\n')
        return ''
    st_file = evm_output + 'evm_to_gff3.stdout'
    log_name = evm_output + 'evm_to_gff3.log'
    st_out = open(st_file, 'w')
    log = open(log_name, 'w')
    
    try:
        subprocess.check_call(args, stdout = st_out, stderr = log, cwd = evm_output)
        #print '> Converted to GFF3 \n'
    except:
        #print 'Could not convert to GFF3\n'
        raise NameError('')
    
    st_out.close()
    log.close()

def combine_gff3(evm_output):
    '''combine the individual evm file in one file'''
    fileName = evm_output + 'evm.out.combined.gff3'
    testFasta = open(fileName, 'w')
    for root, dirs, files in os.walk(evm_output):
        for name in files:
            wd_fasta = os.path.join(root, name)
            if 'evm.out.gff3' in wd_fasta:
                t_file = open(wd_fasta, 'r')
                for line in t_file:
                    testFasta.write(line)
                t_file.close()
                    
    testFasta.close()
    return fileName
    
def evm_pipeline(working_dir, threads, reference, weights, gene_preds, transcripts, proteins, segmentSize, overlapSize):
    ''' Groups the five different calls to run the EVM_pipeline, and concatenates and convert the proper output. 
    It will spit out evm.out.gff3 '''
    
    # Creates output directory
    print '\t###CREATE AN OUTPUT DIRECTORY###\n'
    evm_output = working_dir + 'evm_output/'
    check_create_dir(evm_output)
    #print '> Output directory is: ' + evm_output + '\n'
    
    # Partitions
    print '\t###PARTITIONING THE INPUTS###\n'
    partitions = evm_partitions(evm_output, reference, gene_preds, transcripts, proteins, 
                                segmentSize, overlapSize)
    
    # Write Commands
    print '\t###GROUPING COMMANDS###\n'
    command_list = evm_write_commands(evm_output, reference, weights, gene_preds, transcripts, 
                                     proteins, partitions)
    
    # Run
    print '\t###RUNNING EVM###\n'
    evm_run(evm_output, command_list, threads)
    
    # Combine partitions
    print '\t###COMBINING PARTITIONS###\n'
    evm_combine(evm_output, partitions)
    
    # Convert to GFF3
    print '\t###CONVERTING TO GFF3###\n'
    evm_to_gff3(evm_output, partitions, reference)
    
    # Combine the different chromosomes
    evm_gff3 = combine_gff3(evm_output)
    
    return evm_gff3

def pasa_annot_configuration(pasa_dir, pasa_db):
    '''Creates a PASA annotation configuration file'''
    conf_file = pasa_dir + 'annotCompare.config'
    if os.path.isfile(conf_file): 
        print ('PASA annotation configuration file existed already: ' + conf_file + ' --- skipping\n')
        return conf_file
    conf = open(conf_file, 'w')
    
    lines = ['MYSQLDB='+pasa_db, 
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
             'cDNA_annotation_comparer.dbi:--GENETIC_CODE=<__GENETIC_CODE__>'
             ]
    for line in lines:
        conf.write(line + '\n')   
       
    conf.close()
    #print '> PASA annotation configuration file created in: ' + conf_file + '\n'
    return conf_file

def load_gff3_pasa(pasa_dir, align_conf_file, reference, gff3_file):
    '''Loads a gff3 file into a PASA database '''
    args = ['Load_Current_Gene_Annotations.dbi', '-c', align_conf_file, '-g',
            reference, '-P', gff3_file]
    log_name = pasa_dir + 'load_gff3.log'
    stdout_file = pasa_dir + 'load_gff3.stdout'
    log = open(log_name, 'w')
    stdout_f = open(stdout_file, 'w')
    try:
        load = subprocess.Popen(args, stderr = log, stdout = stdout_f, cwd = pasa_dir)
        load.communicate()
        processID = load.pid
        #print '> GFF3 loaded to PASA DB \n'
    except:
        #print 'Could not load GFF3 to PASA DB\n'
        raise NameError('')
        
    log.close()
    stdout_f.close()
    return processID

def annot_comparison(processID, pasa_dir, pasa_db, annot_conf_file, reference, transcripts_file, n_cpu, valueS):
    '''Loads a gff3 file into a PASA database '''
    '''Loads a gff3 file into a PASA database '''
    
    if valueS is "a":
	print "\t UPDATING WITH ALT-SPLICE"
        args = ['Launch_PASA_pipeline.pl', '--ALT_SPLICE','--TRANSDECODER','--CPU', str(n_cpu) , '-c', annot_conf_file, '-A', '-g',
            reference, '-t', transcripts_file]
    else:
	print "\t UPDATING WITHOUT ALT-SPLICE"
        args = ['Launch_PASA_pipeline.pl', '--TRANSDECODER','--CPU', str(n_cpu) , '-c', annot_conf_file, '-A', '-g',
            reference, '-t', transcripts_file]
    log_name = pasa_dir + 'update_gff3.log'
    log = open(log_name, 'w')
    log_out_name = pasa_dir + 'pasa.out.log'
    out_log = open(log_out_name, 'w')

   

    try:
        update_process = subprocess.check_call(args, stdout = out_log ,stderr = log, cwd = pasa_dir)
    except:
        raise NameError('')

    out_log.close()
    log.close()
    
def parse_pasa_update(round_n, pasa_dir, pasa_db):
    '''Parses through the files in the PASA directory, finds the update file and 
    renames it and returns it'''
    pasa_files = os.listdir(pasa_dir)
    
    pattern_build = '^'+pasa_db+'.gene_structures_post_PASA_updates.[0-9]+.gff3$'
    
    pasa_pattern = re.compile(pattern_build)
    for filename in pasa_files:
        match = re.match(pasa_pattern, filename)
        
        if match:
            update_file = filename
        
    new_filename = pasa_dir + 'annotation.PASAupdated.round'+str(round_n)+'.gff3'
    args = ['mv', update_file, new_filename]


    try:
        subprocess.check_call(args, cwd = pasa_dir)
    except:
        raise NameError('')
    return new_filename

def update_database(n_cpu , round_n,  pasa_dir, pasa_db, align_conf_file, reference, transcripts_file, gff3_file, valueS):
    '''Updates the gff3 file with the PASA database'''       
    print '\t###CREATING CONFIGURATION FILE###\n'
    annot_conf_file = pasa_annot_configuration(pasa_dir, pasa_db)
       
    print '\t###LOADING GFF3 FILE INTO DATABASE###\n'
    processID = load_gff3_pasa(pasa_dir, align_conf_file, reference, gff3_file)
    print '\t###UPDATING GFF3 FILE###\n'
    annot_comparison(processID, pasa_dir, pasa_db, annot_conf_file, reference, 
                     transcripts_file, n_cpu, valueS)
    print '\t###PARSING OUTPUT###\n'
    gff3_out = parse_pasa_update(round_n, pasa_dir, pasa_db)
    return gff3_out

