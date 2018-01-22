#!/usr/bin/env python3

import os
import subprocess
import sys
import tempfile

#=======================================================================================================================

PART = 'partition_EVM_inputs.pl --genome %s --gene_predictions %s --protein_alignments %s --segmentSize %s ' \
             '--overlapSize %s --partition_listing %s'

PART_TRANS = 'partition_EVM_inputs.pl --genome %s --gene_predictions %s --protein_alignments %s --segmentSize %s ' \
             '--overlapSize %s --partition_listing %s --transcript_alignments %s'


EVM_WRITE_COMM = 'write_EVM_commands.pl --genome %s --weights %s --gene_predictions %s --transcript_alignments %s ' \
                 '--protein_alignments %s --output_file_name %s --partitions  %s'

GT_RETAINID = 'gt gff3 -sort -tidy -addintrons -retainids %s'

GT_STATS = 'gt stat -addintrons'
#=======================================================================================================================



def evm_pipeline(working_dir, threads, reference, weights, gene_preds, transcripts, proteins, segmentSize, overlapSize, verbose):
    ''' Groups the five different calls to run the EVM_pipeline, and concatenates and convert the proper output.
    It will spit out evm.out.gff3 '''

    # Partitions
    sys.stdout.write('\t###PARTITIONING THE INPUTS###\n')
    partitions = evm_partitions(working_dir, reference, gene_preds, transcripts, proteins, segmentSize, overlapSize, verbose)

    # Write Commands
    sys.stdout.write('\t###GROUPING COMMANDS###\n')
    command_list = evm_write_commands(working_dir, reference, weights, gene_preds, transcripts, proteins, partitions,verbose)

    # Run
    sys.stdout.write('\t###RUNNING EVM###\n')
    evm_run(working_dir, command_list, threads)

    # Combine partitions
    sys.stdout.write('\t###COMBINING PARTITIONS###\n')
    evm_combine(working_dir, partitions)

    # Convert to GFF3
    sys.stdout.write('\t###CONVERTING TO GFF3###\n')
    evm_to_gff3(working_dir, partitions, reference)

    # Combine the different chromosomes
    evm_gff3 = combine_gff3(working_dir)
    gff3_stat_file = gff3_stats(evm_gff3, working_dir)

    return evm_gff3, gff3_stat_file


def gff3_stats(gff3_file, working_dir):
    '''
    generates stats for gff3 file
    '''

    gt_com = GT_RETAINID % gff3_file
    gt_stat = GT_STATS
    file_out = gff3_file + ".stats"
    file1 = open(file_out, "w")
    err1 = tempfile.NamedTemporaryFile(delete=False, mode="w", prefix="grs", dir=working_dir)
    gt_call = subprocess.Popen(gt_com, stdout=subprocess.PIPE, stderr=err1, shell=True)
    gt_call_stat = subprocess.Popen(gt_stat, stdin=gt_call.stdout, stdout=file1, stderr=err1, shell=True)
    gt_call_stat.communicate()
    file1.close()

    if "evm" in gff3_file:
        sys.stdout.write("\033[31m ### EVM GFF3 STATS ### \n\033[0m")
    else:
        sys.stdout.write("\033[31m ### LOREAN GFF3 STATS ### \n\033[0m")

    with open(file_out) as stats:
        for line in stats:
            sys.stdout.write("\033[32m" + line + "\033[0m")

    return file_out



def evm_partitions(evm_output, reference, gene_preds, transcripts, proteins, segmentSize, overlapSize, verbose):
    '''Genome sequences and  gff3 files are partitioned based on individual contigs,
    and large contigs are segmented into smaller overlapping chunks.'''


    partitions = evm_output + 'partitions_list.out'

    if transcripts == '' and proteins != '':
        cmd = PART % (reference, gene_preds, proteins, segmentSize, overlapSize, partitions)
    else:
        cmd = PART_TRANS % (reference, gene_preds, proteins, segmentSize, overlapSize, partitions, transcripts)

    stdout_file = evm_output + 'partitions.stdout'
    log_name = evm_output + 'partitions.log'
    log = open(log_name, 'w')
    stdout_f = open(stdout_file, 'w')

    try:
        if verbose:
            sys.stderr.write('Executing: %s\n\n' % cmd)
        evm_call = subprocess.Popen(cmd, stdout=stdout_f, stderr=log, cwd=evm_output, shell=True)
        evm_call.communicate()
    except:
        raise NameError( cmd + 'not working')
    log.close()
    stdout_f.close()
    return partitions


def evm_write_commands(evm_output, reference, weights, gene_preds, transcripts, proteins, partitions, verbose):
    '''
    Writes a file with the necessary commnd for the partitions
    '''
    evm_output_file = 'evm.out'
    cmd = EVM_WRITE_COMM % (reference, weights, gene_preds, transcripts, proteins, evm_output_file, partitions)
    command_file = evm_output + 'commands.list'

    log_name = evm_output + 'write_commands.log'
    command = open(command_file, 'w')
    log = open(log_name, 'w')

    try:
        if verbose:
            sys.stderr.write('Executing: %s\n\n' % cmd)
        evm_call = subprocess.Popen(cmd, stdout=command, stderr=log, cwd=evm_output, shell=True)
        evm_call.communicate()
    except:
        raise NameError('')

    command.close()
    log.close()
    return command_file


def evm_run(evm_output, command_list, threads):
    '''
    Runs all the commands in commands.list in parallel
    '''
    args1 = ['cat', command_list]
    args2 = ['parallel', '-j', str(threads), '--']

    # THIS OUTPUT FROM THE WHOLE PIPELINE
    out_file = evm_output + 'evm.out.combined.gff3'
    if os.path.isfile(out_file):
        sys.stdout.write(('\nEVM output existed already: ' + out_file + ' --- skipping\n'))
        return ''

    log_name = evm_output + 'run.log'
    log = open(log_name, 'w')

    try:
        cat = subprocess.Popen(args1, stdout=subprocess.PIPE, cwd=evm_output)
        # evm_call.communicate()
        parallel = subprocess.Popen(args2, stdin=cat.stdout, cwd=evm_output, stderr=log)
        parallel.communicate()
    except:
        raise NameError('')
    log.close()


def evm_combine(evm_output, partitions):
    '''
    perl /home/jose/bin/EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl
    --partitions partitions_list.out --output_file_name evm.out
    '''
    args = ['recombine_EVM_partial_outputs.pl', '--partitions', partitions, '--output_file_name', 'evm.out']

    # THIS OUTPUT FROM THE WHOLE PIPELINE
    out_file = evm_output + 'evm.out.combined.gff3'
    if os.path.isfile(out_file):
        sys.stdout.write(('\nEVM output existed already: ' + out_file + ' --- skipping\n'))
        return ''
    st_file = evm_output + 'combine_partitions.stdout'
    log_name = evm_output + 'combine_partitions.log'
    st_out = open(st_file, 'w')
    log = open(log_name, 'w')

    try:
        evm_call = subprocess.Popen(args, stdout=st_out, stderr=log, cwd=evm_output)
        evm_call.communicate()
    except:
        raise NameError('')

    st_out.close()
    log.close()


def evm_to_gff3(evm_output, partitions, reference):
    '''
    perl /home/jose/bin/EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl
    --partitions partitions_list.out --output evm.out
    --genome ~/Reference/JR2_Chr8/Verticillium_dahliaejr2.GCA_000400815.2.29.dna.chromosome.8.fa
    '''
    args = ['convert_EVM_outputs_to_GFF3.pl', '--partitions', partitions, '--output', 'evm.out', '--genome', reference]

    # THIS OUTPUT FROM THE WHOLE PIPELINE
    out_file = evm_output + 'evm.out.combined.gff3'
    if os.path.isfile(out_file):
        sys.stdout.write(('\nEVM output existed already: ' + out_file + ' --- skipping\n'))
        return ''
    st_file = evm_output + 'evm_to_gff3.stdout'
    log_name = evm_output + 'evm_to_gff3.log'
    st_out = open(st_file, 'w')
    log = open(log_name, 'w')

    try:
        evm_call = subprocess.Popen(args, stdout=st_out, stderr=log, cwd=evm_output)
        evm_call.communicate()
    except:
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


