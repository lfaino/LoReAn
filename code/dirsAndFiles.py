#!/usr/bin/env python3

import errno
import os
import shutil
import subprocess
import sys
import tempfile

import mapping

#==========================================================================================================
# COMMANDS LIST

GFFREAD = 'gffread -o- -T %s'

GTF2BED = 'gtf2bed.py transcript %s'

BEDTOOLS = 'bedtools bamtobed -split -bed12 -i %s'

GT_CHANGE_NAME = 'gt gff3 -setsource update -sort -tidy %s'

#==========================================================================================================


def check_dir(path):
    '''Checks if a directory exists'''
    if not os.path.isdir(path):  # Fast5 dir needs to exist.
        raise IOError(path + ' directory not found')


def check_file(path_file):
    '''Check if a file exists'''
    if not os.path.isfile(path_file):
        raise IOError(path_file + ' file not found')


def check_create_dir(path):
    """
    Checks a directory and creates it if it does not exist
    """

    try:
        os.makedirs(path)  # Try to make it, otherwise don't do anything
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def copy_file(in_file, directory):
    '''Copies a file into a directory'''
    args = ['cp', '-f', in_file, directory]
    try:
        subprocess.check_call(args)
    except:
        raise NameError('')


def change_ids(update, wd, verbose):
    new_name_update = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="update.", suffix=".gff")
    new_name_update_err = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="update.", suffix=".gff.err")
    gt_con = GT_CHANGE_NAME % update
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_con)
    gffread_call = subprocess.Popen(gt_con, stdout=new_name_update,stderr=new_name_update_err, shell=True)
    gffread_call.communicate()
    return (new_name_update.name)


def catTwoBeds(gmap_bam, evm_orig, verbose, wd):
    '''convert in to bed12 and concatenates the two bed12 files'''

    err = tempfile.NamedTemporaryFile()
    outFilename = wd + 'mergedGmapEvm.beforeAssembly.gff3'
    evm = evm_orig
    gtf = evm + ".gtf"
    bed12_evm = evm + ".bed12"
    bed12file = open(bed12_evm, "w")
    gtffile = open(gtf, "w")
    gffread_con = GFFREAD % evm
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gffread_con)
    gffread_call = subprocess.Popen(gffread_con, stdout=gtffile, stderr=err, shell=True)
    gffread_call.communicate()
    err = tempfile.NamedTemporaryFile()
    gft2bed = GTF2BED % gtf
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gft2bed)
    evm_call = subprocess.Popen(gft2bed, stdout=bed12file, stderr=err, shell=True)
    evm_call.communicate()
    #bed12_evm = tmp.name

    bed12_gmap = gmap_bam + ".bed12"
    bed12gmapfile = open(bed12_gmap, "w")
    bedtools = BEDTOOLS % gmap_bam
    err = tempfile.NamedTemporaryFile()
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % bedtools)
    bedtools_call = subprocess.Popen(bedtools, stdout=bed12gmapfile, stderr=err, shell=True)
    bedtools_call.communicate()
    inFile1 = open(bed12_gmap, 'r')
    inFile2 = open(bed12_evm, 'r')

    outFile = open(outFilename, 'w')
    for File in [inFile1, inFile2]:
        for line in File:
            outFile.write(line)
    inFile1.close()
    inFile2.close()
    outFile.close()
    outNameNew = outFilename + 'new.bed'
    lastFile = open(outFilename, 'r')
    o = open(outNameNew, 'w')
    countLine = 0
    for line in lastFile:
        countLine += 1
        linenew = line.split('\t')
        if "evm" in linenew[3]:
            o.write(line)
        else:
            linenew[3] = str(countLine)
            aline = '\t'.join(linenew)
            o.write(aline)
    o.close()
    bed12gmapfile.close()
    bed12file.close()
    gtffile.close()
    return outNameNew


def cat_two_fasta(trinity, consens, long_fasta, wd):
    '''Concatenates the two fasta file into one output'''
    fastas = [trinity, consens, long_fasta]
    outFileFasta = wd + "/allFasta.fasta.clean"
    with open(outFileFasta, 'wb') as outfile:
        for fasta in fastas:
            with open(fasta, 'rb') as fd:
                shutil.copyfileobj(fd, outfile, 1024*1024*10)


    return outFileFasta


def check_gmap(threads_use, type, min_intron_length, max_intron_length, end_exon, gmap_wd, verbose):

    print("\n ### LOREAN IS CHECKING THAT GMAP IS CORRECTLY BUILD ### \n")

    genome_gmap = "/opt/LoReAn/third_party/check_gmap/chr8.testGMAP.fasta"
    long_fasta = "/opt/LoReAn/third_party/check_gmap/exons.testGMAP.fasta"

    long_sam = mapping.gmap('test', genome_gmap, long_fasta, threads_use, type, min_intron_length, max_intron_length,
                            end_exon, gmap_wd, verbose, Fflag=False)

    if os.path.exists(long_sam) and os.path.getsize(long_sam) > 0:
        print("\n" + "\033[32m ### GMAP IS CORRECTLY BUILD ### \n")
        print('\033[0m')

    else:
        print("\n" + "\033[32m ### GMAP REQUIRES A RE-BUILD; THIS WILL TAKE TIME ### \n")
        print('\033[0m')
        log_name = gmap_wd + '/gmap_compile.log'
        err_name = gmap_wd + '/gmap_compile.err'
        log = open(log_name, 'w')
        err = open(err_name, 'w')

        gmap_installation_dir = '/opt/LoReAn/third_party/software/gmap'
        cmd = 'sudo make clean; sudo ./configure ; sudo make ; sudo make install'
        if verbose:
            sys.stderr.write('Executing: %s\n' % cmd)
        gmap_build = subprocess.Popen(cmd, stdout=log, stderr=err, cwd = gmap_installation_dir, shell = True)
        gmap_build.communicate()

        sys.stdout.write("\n### LOREAN FINISHED TO COMPILE ### \n")

        log.close()
        err.close()

        print("\n ### LOREAN IS CHECKING THAT GMAP IS CORRECTLY BUILD (2)### \n")

        long_sam = mapping.gmap('test', genome_gmap, long_fasta, threads_use, type, min_intron_length, max_intron_length,
                                end_exon, gmap_wd, verbose, Fflag=False)

        if os.path.exists(long_sam) and os.path.getsize(long_sam) > 0:
            print("\n" + "\033[32m ### GMAP WORKS CORRECTLY ### \n")
            print('\033[0m')
        else:
            # print("\n" + "\033[31m ### GMAP DID NOT COMPILE CORRECTLY ### \n")
            # print('\033[0m')
            sys.exit("\n" + "\033[31m ### GMAP DID NOT COMPILE CORRECTLY ### \n")


def augustus_species_func(home):
    #print(home)
    check_species = 'augustus --species=help'
    process = subprocess.Popen(check_species, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out_augustus, err_augustus = process.communicate()
    #list_file = [os.path.join(home, o) for o in os.listdir(home) if os.path.isfile(os.path.join(home, o)) and ".bashrc.lorean" == o]
    with open("/opt/LoReAn/third_party/conf_files/environment", "r") as bashrc:
        for path in bashrc:
            if path.startswith("AUGUSTUS_CONFIG_PATH"):
                augustus_specie_dir = path.split("=")[1].rsplit("\"")[1]
                augustus_species = [d for d in os.listdir(augustus_specie_dir)]
    return augustus_species, err_augustus
