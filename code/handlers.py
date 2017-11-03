#!/usr/bin/env python3

import dirs_and_files as logistic
import multithread_large_fasta as multiple
import transcript_assembly as transcripts


def BrakerAAT(queue,ref, bamFile, species_name, protein_evidence, threads, fungus, list_fasta_names, wd, verbose):
    '''Handles Braker and AAT so that we can run them in parallel'''
    # DIVIDE THREADS BY 2
    use = (round(int(threads) / 2) - 1)
    use_gmes = str(use)
    aat_wd = wd + 'AAT/'
    logistic.check_create_dir(aat_wd)
    while True:
        dummy = queue.get()
        if dummy == 0:
            transcripts.braker_call(wd, ref, bamFile, species_name, use, fungus, verbose)
        if dummy == 1:
            multiple.aat_multi(use, protein_evidence, list_fasta_names, aat_wd, verbose)
        queue.task_done()


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
            multiple.augustus_multi(use, species, list_fasta_names, augustus_wd, verbose)
        if dummy == 1:
            multiple.aat_multi(use, protein_evidence, list_fasta_names, aat_wd, verbose)
        if dummy == 2:
            transcripts.gmes_call(gmes_wd, ref, fungus, use_gmes, verbose)
        queue.task_done()

