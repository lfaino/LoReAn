#!/usr/bin/env python

''' 
MASTER THESIS PROJECT
Author: Jose A. Espejo
Date: September 2015 - March 2016

handle programs that can run in parallel
'''


import os 
import subprocess
import protein_alignment
import transcript_assembly as transcripts
import dirs_and_files as logistic
import multithread_large_fasta as multiple


def BrakerAAT(queue, ref, bamFile, species_name, proteinEvidence, threads, fungus, list_fasta_names, wd):
    '''Handles Braker and AAT so that we can run them in parallel'''
    ###DIVIDE THREADS BY 2
    use = (int(threads)/2)
    while True:
        dummy = queue.get()
        if dummy == 0:
            transcripts.braker_call(wd, ref, bamFile, species_name, use, fungus)
        if dummy == 1:
            aat_wd = wd+'AAT/'
            logistic.check_create_dir(aat_wd)
            proteinFastaFile = os.path.abspath(proteinEvidence)
            AATflag = multiple.aat_multi(ref, use, proteinFastaFile, list_fasta_names, aat_wd)
            if AATflag:
                protein_alignment.parseAAT(aat_wd)
        queue.task_done()
    return 
        
       
def AugustGmesAAT(queue, ref, species, protein_evidence, threads, fungus, list_fasta_names, wd):
    '''Handles Braker and AAT so that we can run them in parallel'''
    ###DIVIDE THREADS BY 3
    use = (int(threads)/3)
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
            multiple.augustus_multi(ref, use, species, list_fasta_names, augustus_wd)
        if dummy == 1:
            multiple.aat_multi(ref, use, protein_evidence, list_fasta_names, aat_wd)
        if dummy == 2:
            transcripts.gmes_call(gmes_wd, ref, fungus, use_gmes)
        queue.task_done()
    return 



def AugustGmes(queue, ref, species, fungus, threads, list_fasta_names, wd):
    '''Handles Braker and AAT so that we can run them in parallel'''
    use = (int(threads)/2)
    use_gmes = str(use)
    ###DIVIDE THREADS BY 2
    
    augustus_wd = wd + 'augustus'
    logistic.check_create_dir(augustus_wd)
    gmes_wd = wd + 'gmes/'
    logistic.check_create_dir(gmes_wd)
    
    while True:
        dummy = queue.get()
        if dummy == 0:
            multiple.augustus_multi(ref, use, species, list_fasta_names, augustus_wd)
        if dummy == 1:
            transcripts.gmes_call(gmes_wd, ref, fungus, use_gmes)
            
        queue.task_done()
    return 



def main():
    '''Main body of the function'''
    ref = os.path.abspath(sys.argv[1])
    threads = sys.argv[2]
    species = sys.argv[3]
    wd = os.path.abspath(sys.argv[4])
    AugustGmes(ref, threads, species, wd)
    
    print '\n\n\n###############\n###FINISHED###\n###############\n\n'

    
    
    


if __name__ == '__main__':
    main()
 

