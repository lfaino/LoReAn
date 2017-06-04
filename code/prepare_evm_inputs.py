#!/usr/bin/env python3

import os
import subprocess


def convert_augustus(aug_file, wd):
    '''Converts augustus.gff to augustus.gff3 (from BRAKER1) using the EVM
    script EVMUtils/misc/augustus_GTF_to_EVM_GFF3.pl which needs to be in PATH
    '''
    print('\t###CONVERTING AUGUSTUS TO GFF3###\n')
    args = ['augustus_GTF_to_EVM_GFF3.pl', aug_file]
    #COMMANDS.append(' '.join(args))
    out_file = aug_file + '3'
    if os.path.isfile(out_file):
        print((
            'Augustus GFF3 file existed already: ' +
            out_file +
            ' --- skipping\n'))
        return out_file
    log_name = wd + '.augustus_GTF_to_EVM_GFF3.pl.log'
    log = open(log_name, 'w')
    out_f = open(out_file, 'w')

    try:
        subprocess.check_call(args, stdout=out_f, stderr=log)

        # print '> Augustus to GFF3 completed: ' + out_file
    except:
        # print ' Augustus to GFF3 failed'
        raise NameError('')

    log.close()
    out_f.close()

    return out_file


def convert_genemark(genemark_file, wd):
    '''Converts genemark.gtf to genemark.gff3 (from BRAKER1) using gtf2gff3.pl,
    which needs to be in PATH'''

    print('\t###CONVERTING GENEMARK TO GFF3###\n')
    args = ['gtf2gff3.pl', genemark_file]
    #COMMANDS.append(' '.join(args))

    out_file = genemark_file + '.gff3'

    if os.path.isfile(out_file):
        print((
            'GeneMark GFF3 file existed already: ' +
            out_file +
            ' --- skipping\n'))
        return out_file

    log_name = wd + '.genemark_GTF_to_EVM_GFF3.pl.log'
    log = open(log_name, 'w')
    out_f = open(out_file, 'w')

    try:
        subprocess.check_call(args, stdout=out_f, stderr=log)

        # print '> Genemark to GFF3 completed: ' + out_file + '\n'
    except:
        # print ' Genemark to GFF3 failed'
        raise NameError('')

    log.close()
    out_f.close()

    return out_file


def move_single_file(filename, key, evm_dir, new_file_d):
    '''Moves a single file into the directory and appends the new path to the dictionary'''
    args = ['cp', filename, evm_dir]
    #COMMANDS.append(' '.join(args))

    true_filename = filename.split('/')[-1]
    out_file = evm_dir + true_filename
    if os.path.isfile(out_file):
        print(('File in EVM_dir already: ' + out_file + ' --- skipping\n'))
        new_file_d[key] = out_file
        return new_file_d

    try:
        subprocess.check_call(args)
        new_file_d[key] = out_file
        return new_file_d
    except:
        # print 'Could not move ' + filename
        raise NameError('')


def move_cat_files(file_list, key, evm_dir, new_file_d):
    '''Moves and concatenate files to evm dir (case of GFF3 when using long
    and short reads)'''
    args = ['cat'] + file_list
    #COMMANDS.append(' '.join(args))

    out_file = evm_dir + key + '.gff3'
    if os.path.isfile(out_file):
        print(('File in EVM_dir already: ' + out_file + ' --- skipping\n'))
        new_file_d[key] = out_file
        return new_file_d

    file_ = open(out_file, 'w')
    try:
        subprocess.check_call(args, stdout=file_)
        new_file_d[key] = out_file
        file_.close()
        return new_file_d
    except:
        print('Could not move ' + filename)
        raise NameError('')


def move_EVM_inputs(evm_dir, inputs):
    '''Takes a dictionary with files that are inputs for EVM and groups them in
    the same directory'''
    print('\t###MOVING IMPORTANT FILES###\n')
    new_files = {}
    for key, filename in list(inputs.items()):
        if isinstance(
                filename,
                list):  # FOR THE GFF3 alignment files in case of short & long reads
            new_files = move_cat_files(filename, key, evm_dir, new_files)
        else:
            new_files = move_single_file(filename, key, evm_dir, new_files)

    # print '> EVM input dir full of files: ' + evm_dir
    return new_files


def cat_EVM_inputs(evm_dir):  # , inputs):
    '''Takes the files in EVM input directory and concatenates the needed
    files to prepare the EVM command. Augustus, Genemark and Transdecoder
    go into gene_predictions.gff3 and pasa asemblies and transcript
    alignments go into transcripts.gff3'''
    # GENE PREDICTIONS

    print('\t###CONCATENATING FILES###\n')

    # GENE PREDICTION
    file_list = []

    ab_initio_list = ['cat']
    protein_list = []
    transcript_list = []
    list_soft = []
    transcript_file = ''
    protein_file = ''

    for root, dirs, files in os.walk(evm_dir):
        for name in files:
            if 'augustus' in name:
                ab_initio_list.append(os.path.join(root, name))
                list_soft.append('augustus')
            elif 'genemark' in name:
                ab_initio_list.append(os.path.join(root, name))
                list_soft.append('genemark')
            elif 'PASA' in name or 'pasa' in name:
                transcript_file = os.path.join(root, name)
                transcript_list.append(os.path.join(root, name))
                list_soft.append('pasa')
            elif 'protein' in name:
                protein_file = os.path.join(root, name)
                protein_list.append(os.path.join(root, name))
                list_soft.append('aat')
            elif 'trinity' in name:
                ab_initio_list.append(os.path.join(root, name))
                list_soft.append('gmap')

    pred_filename = evm_dir + 'gene_predictions.gff3'

    if os.path.isfile(pred_filename):
        print(('Gene predictions GFF3 file existed already: ' +
               pred_filename + ' --- skipping\n'))
    else:

        # print '\nCMD: ' + ' '.join(ab_initio_list) + '\n'
        pred_file = open(pred_filename, 'w')
        try:
            subprocess.check_call(
                ab_initio_list,
                stdout=pred_file,
                cwd=evm_dir)
            # print '> Gene prediction concatenation completed'
        except:
            # print 'Gene prediction concatenation failed'
            raise NameError('')
        pred_file.close()

    return (list_soft, pred_filename, transcript_file, protein_file)


def group_EVM_inputs(evm_dir, inputs):
    '''Moves all the inputs to EVM directory and concatenates them
    in the same file'''
    # Move
    new_inputs = move_EVM_inputs(evm_dir, inputs)
    # Concatenate
    list_soft, pred_file, transcript_file, protein_file = cat_EVM_inputs(
        evm_dir)
    return (list_soft, pred_file, transcript_file, protein_file)


def evm_weight(evm_dir, weights_dic, evidences, pasa_name, gmap_name):
    '''Writes a weight file "weights.txt" on evm_dir'''
    w_filename = evm_dir + 'weights.txt'
    list_match = []

    evidence_dic = {
        'GeneMark.hmm': 'ABINITIO_PREDICTION',
        'Augustus': 'ABINITIO_PREDICTION',
        'AAT': 'PROTEIN',
        pasa_name: 'TRANSCRIPT',
        gmap_name: 'ABINITIO_PREDICTION'}
    software_links = {
        'genemark': 'GeneMark.hmm',
        'augustus': 'Augustus',
        'aat': 'AAT',
        'pasa': pasa_name,
        'gmap': gmap_name}

    for software in software_links:
        if software in evidences:
            list_match.append(software_links[software])
    w_file = open(w_filename, 'w')
    for present_soft in list_match:
        if present_soft in evidence_dic:
            #            print present_soft
            w_file.write('\t'.join(
                [evidence_dic[present_soft], present_soft, weights_dic[present_soft]]))
            w_file.write('\n')

    w_file.close()

    return w_filename
