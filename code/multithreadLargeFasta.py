#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool

count_sequences = 0
count_sequences_aat = 0
length_cluster = 0
length_cluster_aat = 0

#==========================================================================================================
# COMMANDS LIST

AUGUSTUS = 'augustus --species=%s %s'

AAT = 'AAT.pl -P -b -q %s -s %s r"--dps" r" \'-f 100 -i 30 -a 200\'" r"--filter" r"\'-c 10\'", r"--nap", r"\'-x 10\'"'

#==========================================================================================================


def single_fasta(ref, wd_split):
    """
    From a fasta file make single files with each sequence
    """

    fasta_file = open(ref, 'r')
    single_fasta_list = []
    count = 0
    dict_ref_name = {}
    ref_rename = ref + ".rename.fasta"
    with open(ref_rename, "w") as fh:
        for record in SeqIO.parse(fasta_file, "fasta"):
            count += 1
            new_name = "seq" + str(count)
            dict_ref_name[new_name] = record.id
            new_rec = SeqRecord(record.seq, new_name, '', '')
            fasta_name = wd_split + '/' + new_name + '.fasta'
            single_fasta_list.append(fasta_name)
            output_handle = open(fasta_name, "w")
            SeqIO.write(new_rec, output_handle, "fasta")
            SeqIO.write(new_rec, fh, "fasta")
            output_handle.close()
    return single_fasta_list, dict_ref_name, ref_rename


def augustus_multi(threads, species, single_fasta_list, wd, verbose):
    '''handles the assembly process and parsing in a multithreaded way'''
    if int(threads) < 1:
        threads = 1
    all_augustus = []
    augustus = [wd, species, verbose]
    for record in single_fasta_list:
        single_command = augustus + [record]
        all_augustus.append(single_command)
    sys.stdout.write ("###RUNNING AUGUSTUS ###\n")
    with Pool(processes=int(threads), maxtasksperchild=1000) as pool:
        pool.map(augustus_call, all_augustus)

    parseAugustus(wd)
    return


def augustus_call(all_augustus):
    '''
    augustus call
    '''
    cmd = AUGUSTUS % (all_augustus[1], all_augustus[3])
    chromo = all_augustus[3].split('/')[-1]
    wd_augu = all_augustus[0] + '/' + chromo + '.augustus.gff'
    if os.path.exists(wd_augu) and os.path.getsize(wd_augu) > 0:
        sys.stderr.write('Already executed: %s\n' % cmd)
        pass
    else:
        log_name = wd_augu
        log = open(log_name, 'w')
        log_name_err = all_augustus[0] + 'augustus.err.log'
        log_e = open(log_name_err, 'w')
        try:
            if all_augustus[2]:
                sys.stderr.write('Executing: %s\n' % cmd)
            augustus = subprocess.Popen(cmd, stderr=log_e, stdout=log, cwd=all_augustus[0], shell=True)
            augustus.communicate()
            if all_augustus[2]:
                sys.stderr.write('Done: %s\n' % cmd)
        except:
            raise NameError('Augustus Failed')
        log.close()
        log_e.close()

    return all_augustus[0]


def parseAugustus(wd):
    '''From all the augustus output after the multithread generate a single gff file'''

    fileName = wd + '/augustus.gff'
    with open(fileName, 'wb') as outfile:
        for root, dirs, files in os.walk(wd):
            for name in files:
                if name.endswith("fasta.augustus.gff"):
                    filename = root + '/' + name
                    with open(filename, 'rb') as fd:
                        shutil.copyfileobj(fd, outfile, 1024*1024*10)

    return fileName


def aat_multi(threads, protein_evidence, single_fasta_list, wd, verbose):
    '''handles the assembly process and parsing in a multithreaded way'''

    if int(threads) < 1:
        threads = 1
    all_aat = []
    aat = [wd, protein_evidence, verbose]
    for record in single_fasta_list:
        single_command = aat + [record]
        all_aat.append(single_command)
    sys.stdout.write ("###RUNNING AAT ###\n")
    with Pool(processes=int(threads), maxtasksperchild=1000) as p:
        p.map(aat_call, all_aat)
    parseAAT(wd)
    return


def aat_call(all_aat):
    '''Calls genome guided trinity on the BAM file to generate assembled transcripts'''


    cmd = AAT % (all_aat[3], all_aat[1] )

    chr_fasta_name = all_aat[3].split("/")[-1]
    chr_name = chr_fasta_name.split(".")[0]
    prot_fasta_name = all_aat[1].split("/")[-1]
    file_created = all_aat[0] + "/" + chr_name + "." + prot_fasta_name + ".nap.btab"

    if os.path.exists(file_created) and os.path.getsize(file_created) > 0:
        sys.stderr.write('Already executed: %s\n' % cmd)
        pass
    else:
        log_name = all_aat[3] + 'AAT.log'
        log = open(log_name, 'w')
        stdout_f = open(all_aat[3] + 'AAT.stdout', 'w')
        if all_aat[2]:
            sys.stderr.write('Executing: %s\n' % cmd)
        aat_process = subprocess.Popen(cmd, stderr=log, stdout=stdout_f, cwd=all_aat[0], shell=True)
        aat_process.communicate()
        if all_aat[2]:
            sys.stderr.write('Done: %s\n' % cmd)
        log.close()
        stdout_f.close()

    return


def parseAAT(wd):
    '''From all the augustus output after the multithread generate a single gff file'''

    fileName = wd + '/AAT.gff'

    with open(fileName, 'wb') as outfile:
        for root, dirs, files in os.walk(wd):
            for name in files:
                if name.endswith("nap.btab"):
                    filename = root + '/' + name
                    with open(filename, 'rb') as fd:
                        shutil.copyfileobj(fd, outfile, 1024*1024*10)

    args_btab = ['AAT_btab_to_gff3.pl', fileName, 'P', ]
    outFilenameGff = wd + '/protein_evidence.gff3'
    stdout_file = open(outFilenameGff, 'w')
    outFilenameGff_err = wd + '/protein_evidence.err'
    stderr_file = open(outFilenameGff_err, 'w')
    aat_call = subprocess.Popen(args_btab, stdout=stdout_file, stderr=stderr_file, cwd=wd)
    aat_call.communicate()
    stdout_file.close()
    stderr_file.close()
    return outFilenameGff


if __name__ == '__main__':
    single_fasta(*sys.argv[1:])