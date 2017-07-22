#!/usr/bin/env python3

import os
import protein_alignment
import transcript_assembly as transcripts
import dirs_and_files as logistic
import multithread_large_fasta as multiple


def BrakerAAT(queue,ref, bamFile, species_name, proteinEvidence, threads, fungus, list_fasta_names, wd):
    '''Handles Braker and AAT so that we can run them in parallel'''
    # DIVIDE THREADS BY 2
    use = (round(int(threads) / 2) - 1)
    while True:
        dummy = queue.get()
        if dummy == 0:
            print ("\n###RUNNING BRAKER1 ###\n")
            transcripts.braker_call(wd, ref, bamFile, species_name, use, fungus)
        if dummy == 1:
            aat_wd = wd + 'AAT/'
            logistic.check_create_dir(aat_wd)
            proteinFastaFile = os.path.abspath(proteinEvidence)
            print ("\n### RUNNING AAT ###\n")
            AATflag = multiple.aat_multi(ref, use, proteinFastaFile, list_fasta_names, aat_wd)
            if AATflag:
                protein_alignment.parseAAT(aat_wd)
        queue.task_done()
    return


def AugustGmesAAT(queue, ref, species, protein_evidence, threads, fungus, list_fasta_names, wd, verbose):
    use = (round(int(threads) / 3)-1)
    use_gmes = str(use)
    augustus_wd = wd + 'augustus/'
    logistic.check_create_dir(augustus_wd)
    gmes_wd = wd + 'gmes/'
    logistic.check_create_dir(gmes_wd)
    aat_wd = wd + 'AAT/'
    logistic.check_create_dir(aat_wd)
    while True:
        dummy = queue.get()
        if dummy == 0:
            print ("\n###RUNNING AUGUSTUS ###\n")
            multiple.augustus_multi(use, species, list_fasta_names, augustus_wd, verbose)
        if dummy == 1:
            print ("\n###RUNNING AAT ###\n")
            multiple.aat_multi(ref, use, protein_evidence, list_fasta_names, aat_wd,)
        if dummy == 2:
            print ("\n###RUNNING GENEMARK ###\n")
            transcripts.gmes_call(gmes_wd, ref, fungus, use_gmes)
        queue.task_done()
    return


def AugustGmes(queue, ref, species, fungus, threads, list_fasta_names, wd):
    '''Handles Braker and AAT so that we can run them in parallel'''
    use = (round(int(threads) / 2)-1)
    use_gmes = str(use)
    # DIVIDE THREADS BY 2

    augustus_wd = wd + 'augustus'
    logistic.check_create_dir(augustus_wd)
    gmes_wd = wd + 'gmes/'
    logistic.check_create_dir(gmes_wd)

    while True:
        dummy = queue.get()
        if dummy == 0:
            multiple.augustus_multi(
                ref, use, species, list_fasta_names, augustus_wd)
        if dummy == 1:
            transcripts.gmes_call(gmes_wd, ref, fungus, use_gmes)

        queue.task_done()
    return


