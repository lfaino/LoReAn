#! /usr/bin/env python

###############
###TO DO AGENDA###
###############

# DISCARD EVIDENCE BEFORE WRITING FASTA

# ALSO IMPLEMENT INTERRUPTIONS IF FILES IN CONSENSUS PIPELINE

# Implement interruptions at certain points if wished


''' 
MASTER THESIS PROJECT
Author: Jose A. Espejo
Date: September 2015 - March 2016

Pipeline to go from Oxford Nanopore's long reads to all necessary for EVM modeler
in Verticillium dahlie.

''' 

###############
###IMPORTS###
###############

##LIBRARIES

from __future__ import division
from multiprocessing import Pool
import sys
import subprocess
import argparse
import os
import errno
import re
import shutil
from Bio import SeqIO
from BCBio import GFF
from operator import itemgetter
from Queue import Queue
from threading import Thread
import itertools
import shutil


##OTHER SCRIPTS
import dirs_and_files as logistic
import mapping
import transcript_assembly as transcripts
import prepare_evm_inputs as inputEvm
import evm_pipeline
import consensus_iAssembler as consensus
import intersect_EVM_consensus as intersectEvm
import protein_alignment
import collect_only as collect
import handlers as handler
import get_right_strand as grs
import multithread_large_fasta as multiple
import reduceUTRs as utrs

####################################
### CHEKS BEFORE START LOREAN ######
####################################


os.system('sudo usermod -d /var/lib/mysql/ mysql')
os.system('sudo /etc/init.d/mysql start')
os.system('mysql --user="root" --password="lorean" --execute="set global sql_mode=\'STRICT_TRANS_TABLES,ERROR_FOR_DIVISION_BY_ZERO,NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION\';"')

if not os.path.isfile("/home/lorean/.gm_key"):
    if os.path.isfile("/data/gm_key"):
        os.system('cp /data/gm_key /home/lorean/.gm_key')
    else:
        print 'Key for GeneMark-ES not found
        LOREAN STOPS HERE. Please, place the GeneMark-ES key in the folder where you have your data.'
else:



###############
###FUNCTIONS###
###############


    def arguments():
        '''Parses the arguments from the program invocation'''
        
        #Call the argument parse
        parser = argparse.ArgumentParser(prog='lorean', 
                                        usage='%(prog)s [options] reference species_name',
                                        description = 'LoReAn - Automated genome annotation pipeline that integrates long reads',
                                        epilog = 'Jose Espejo Valle-Inclan - April 2016') 
        
        #Specify arguments
        #parser.add_argument('--example', nargs='?', const=1, type=int, default=1)

        parser.add_argument("ref", 
                            help="Path to reference file")
        parser.add_argument("species", 
                            help="Species name for AUGUSTUS training. No re-training if species already present in AUGUSTUS config folder")
        parser.add_argument("--stranded", 
                            help="Run LoReAn on stranded mode [FALSE]",
                            action = 'store_true')
        parser.add_argument("--fungus", 
                            help="Run LoReAn on fungal species [FALSE]",
                            action = 'store_true')
        parser.add_argument("--collect_only", 
                            help="Collect only assebmled transcripts [FALSE]",
                            action = 'store_true')
        parser.add_argument("--only_unitigs", 
                            help="Removes gene models that are not supported by long reads [FALSE]",
                            action = 'store_true') 
        ###TO CHECK WITH JOSE HERE
        parser.add_argument("--no_braker", 
                            help="Do not run braker if short or long reads are supplied [FALSE]",
                            action = 'store_true')
        parser.add_argument("--keep_tmp", 
                            help="Keep temporary files [FALSE]",
                            action = 'store_true')
        #parser.add_argument("--PacBio_non_CCS", 
                            #help="Long reads are not from PacBio CCS (Circular consensus sequences) [FALSE]",
                            #action = 'store_true')
        parser.add_argument('--short_reads', nargs="?", default="",
                            help="Path to short reads FASTQ. If paired end, comma-separated (1-1.fq,1-2.fq). BAM sorted files are allowed; the extension of the file should be filename.sorted.bam []",
                            metavar = 'FASTQ_file')
        parser.add_argument('--long_reads', nargs="?", default="",
                            help="Path to long reads FASTQ []",
                            metavar = 'FASTQ_file')
        parser.add_argument('--protein_evidence', nargs="?", default="",
                            help="Path to protein sequences FASTA file []",
                            metavar = 'FASTA_prot')
        #parser.add_argument('--min_long_read', nargs="?", default=100,
                            #help="Filter out long reads shorter than this value (shorter reads may affect mapping and assembling) [100]", type=int)
        parser.add_argument('--max_long_read', nargs="?", default=20000,
                            help="Filter out long reads longer than this value (longer reads may affect mapping and assembling) [20000]", type=int)
        parser.add_argument("--pasa_db", nargs="?", default="pipeline_run",
                            help="PASA database name [pipeline_run]")
        parser.add_argument('-wd', "--working_dir", nargs="?", default="./",
                            help="Working directory (will create if not present) [./]")
        parser.add_argument('-t', "--threads", nargs="?", default="1",
                            help="Number of threads [1]",
                            metavar = 'N')
        parser.add_argument("--overhang", nargs="?", default="20",
                            help="CAP3 max overhang percent length; this value should be > 3 [20]",
                            metavar = 'N')
        parser.add_argument("--augustus", nargs="?", default="1",
                            help="Weight assigned to AUGUSTUS evidence for EVM [1]",
                            metavar = 'N')    
        parser.add_argument("--genemark", nargs="?", default="1",
                            help="Weight assigned to GENEMARK evidence for EVM [1]",
                            metavar = 'N')
        parser.add_argument("--trinity", nargs="?", default="1",
                            help="Weight assigned to Trinity mapped with GMAP evidence for EVM [1]",
                            metavar = 'N')    
        #parser.add_argument("--alignment", nargs="?", default="10",
                            #help="Weight assigned to short read ALIGNMENT evidence for EVM [10]",
                            #metavar = 'N')
        parser.add_argument("--pasa", nargs="?", default="5",
                            help="Weight assigned to PASA evidence for EVM [5]",
                            metavar = 'N')
        parser.add_argument("--AAT", nargs="?", default="1",
                            help="Weight assigned to AAT protein evidence for EVM [1]",
                            metavar = 'N')
        parser.add_argument("--segmentSize", nargs="?", default="100000",
                            help="Segment size for EVM partitions [100000]",
                            metavar = 'N')
        parser.add_argument("--overlapSize", nargs="?", default="10000",
                            help="Overlap size for EVM partitions [10000]",
                            metavar = 'N')    
        parser.add_argument("--min_intron_length", nargs="?", default="9",
                            help="Minimal intron length for GMAP [9]",
                            metavar = 'N')
        parser.add_argument("--max_intron_length", nargs="?", default="1000",
                            help="Maximal intron length for GMAP, STAR and TRINITY [1000]",
                            metavar = 'N')
        parser.add_argument("-H", nargs="?", default="20",
                            help="Minimal length for end exon with GMAP [20]",
                            metavar = 'N')       
        #parser.add_argument("--cluster_min_read_length", nargs="?", default="100",
                            #help="Minimal read length to be considered for CONSENSUS pipeline [100]",
                            #metavar = 'N')   
        parser.add_argument("--cluster_min_evidence", nargs="?", default="5",
                            help="Minimal evidence needed to form a cluster [5]",
                            metavar = 'N') 
        parser.add_argument("--cluster_max_evidence", nargs="?", default="5000",
                            help="Maximal evidence to form a cluster.Prevents the clustering or rRNA genes i.e. [5000]",
                            metavar = 'N') 
        #parser.add_argument("--mergeDistance", nargs="?", default="1",
                            #help="Minimal distance for reads to be merged [1]",
                            #metavar = 'N') 
        parser.add_argument("--assembly_overlapLength", nargs="?", default="200",
                            help="Minimal length (in nt) of overlap for ASSEMBLY [200]",
                            metavar = 'N') 
        parser.add_argument("--assembly_percentIdentity", nargs="?", default="97",
                            help="Minimal identity for the ASSEMBLY (95-100) [97]",
                            metavar = 'N')
        parser.add_argument("--assembly_readThreshold", nargs="?", default="0.3",
                            help="Fraction of reads supporting an assembled UNITIG to keep on the ASSEMBLY (0.1-1) [0.3]",
                            metavar = 'F')    
        parser.add_argument("--no_EVM", 
                            help="Run until the preparation of EVM inputs [FALSE]",
                            action = 'store_true')       
        parser.add_argument("--no_consensus", 
                            help="Do not run the long reads consensus pipeline [FALSE]",
                            action = 'store_true')
        parser.add_argument("--no_update", 
                            help="Do not run the PASA update[FALSE]",
                            action = 'store_true')

        
        #Parse and set defaults for optional arguments
        args = parser.parse_args()       
        return args
    ###############
    ###MAIN###
    ###############


    def main():
        '''Core of the program'''
        
        #Parse the arguments
        args = arguments()

        #Useful variables for later
        wd = os.path.abspath(args.working_dir) + '/'
        wd_gmap = wd
        ref = os.path.abspath(args.ref)

        FinalFiles = [] #STORE THE IMPORTANT OUTPUT FILES
        
        #Arrange stuff before start
        logistic.check_create_dir(wd)
        logistic.check_file(ref)
                
        ###COLLECT ONLY ONLY RUNS PART OF THE CONSENSUS PIPELINE

        if not args.collect_only:
            list_fasta_names = multiple.single_fasta(ref, wd)
            if args.short_reads != '' or args.long_reads != '':
                print '='*70
                print '\n###MAPPING###\n'
                    
                ##SHORT READS
                if 'fastq' in args.short_reads or 'fq' in args.short_reads :
                    if ',' in args.short_reads:
                        pairedEndFiles = args.short_reads.split(',')
                        short_1 = c
                        short_2 = os.path.abspath(pairedEndFiles[1])
                        short_reads_file = [short_1, short_2]
                        print "Paired short reads detected: " + ' '.join(short_reads_file)
                    else:
                        short_reads_file = os.path.abspath(args.short_reads)
                        print "Short reads detected: " + short_reads_file

                    #Map with STAR
                    short_sorted_bam = mapping.star(ref, short_reads_file, args.threads, args.max_intron_length, wd) 
                    
                    #Keep the output
                    FinalFiles.append(short_sorted_bam)
                ##BAM SORTED FILES GET IN HERE
                elif 'sorted.bam' in args.short_reads:
                    short_sorted_bam = os.path.abspath(args.short_reads)
                    
                
                
                else:
                    short_sorted_bam = False
                    print 'No short reads file'
                    
                ##LONG READS
                if 'fastq' in args.long_reads or 'fq' in args.long_reads :
                    ### with this operation, reads are filtered for their length. Nanopore reads can be chimaras or sequencing artefacts. filtering on length reduces the amount of sequencing artefacts
                    print "##FILTERING OUT LONG READS"
                    long_fastq, filter_count = mapping.filterLongReads(args.long_reads, args.assembly_overlapLength, args.max_long_read, wd)
                    
                        
                        
                    if filter_count != 0:
                        print "\n> " + str(filter_count) + " reads kept\n"
                    if not short_sorted_bam:
                        ##If short reads have been mapped dont do it
                        print '##GMAP'
                        long_sam = mapping.gmap('sam', ref, long_fastq, args.threads, 'samse', args.min_intron_length, args.max_intron_length, args.H, wd, Fflag = False)
                        split_option = 'split'        
                        
                        #Convert to sorted BAM
                        long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, args.threads, wd) 
                        
                        #Keep the output
                        FinalFiles.append(long_sorted_bam)
                    else:
                        long_sorted_bam = False

                else:
                    print 'No long reads file'
                    long_sorted_bam = False
                print '='*70
                
                
                if short_sorted_bam: #If there are short reads, these will serve to the transcript assembly pipeline
                    default_bam = short_sorted_bam
                else:
                    default_bam = long_sorted_bam

                ##TRANSCRIPT ASSEMBLY

                #TRINITY
                print '='*70
                print '\n###TRINITY###\n'
                trinity_out = transcripts.trinity(default_bam, wd, args.max_intron_length, args.threads)
                print '='*70
                trinityGFF3 = mapping.gmap('trin', ref, trinity_out, args.threads, 'gff3_gene', args.min_intron_length, args.max_intron_length, args.H, wd, Fflag = True )
                
                    
                #PASA Pipeline
                print '='*70
                print '\n###PASA###\n'
                #Create PASA folder and configuration file
                print '##Configuration\n'
                pasa_dir = wd + 'PASA/'
                logistic.check_create_dir(pasa_dir)
                align_pasa_conf = transcripts.pasa_configuration(pasa_dir, args.pasa_db)
                #Launch PASA
                print '\n##LAUNCH PASA\n'
                pasa_gff3 = transcripts.pasa_call(pasa_dir, align_pasa_conf, args.pasa_db, ref, trinity_out, args.max_intron_length, args.threads)
                print '='*70
                #TransDecoder 
                #print '='*70
                #print '\n###TRANSDECODER###\n'
                #Two parts, LongORFs and Predict
                #Create a directory
                #td_dir = wd + 'transdecoder/'
                #logistic.check_create_dir(td_dir)
                #transD_gff3 = transcripts.transdecoder(td_dir, trinity_out, args.threads)
                print '='*70
                ##HERE WE PARALLELIZE PROCESSES WHEN MULTIPLE THREADS ARE USED
                if not args.no_braker:
                    if (args.protein_evidence != '' and args.short_reads) or (args.protein_evidence != '' and args.long_reads): ### USING PROTEINS AND READS 
                        queue = Queue()
                        print '='*70
                        print '\n###RUNNING BRAKER1 AND AAT###\n'
                        print '='*70
                        for i in range(2):
                            queue.put(i) #QUEUE WITH A ZERO AND A ONE
                        for i in range(2):
                            t = Thread(target=handler.BrakerAAT, args =(queue, ref, default_bam, args.species, args.protein_evidence, args.threads, args.fungus, list_fasta_names, wd))
                            t.daemon = True
                            t.start()
                        queue.join()
                    elif args.short_reads or args.long_reads:  ### USING ONLY READS WITH BRAKER
                        print '='*70
                        print '\n###BRAKER1###\n'
                        print '='*70
                        braker_out = transcripts.braker_call(wd, ref, default_bam, args.species, args.threads, arg.fungus)
                        print '='*70
                
                elif args.protein_evidence != '' and args.short_reads != '':
                    check_species = ['augustus', '--species=help']
                    process = subprocess.Popen(check_species, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
                    out, err = process.communicate()
                    protein_loc = os.path.abspath(args.protein_evidence)
                    if args.species in err:
                        print '='*70
                        print '\n###RUNNING AUGUSTUS, GENEMARK-ES AND AAT###\n'
                        print '='*70
                        queue = Queue()
                        for i in range(3):
                            queue.put(i) #QUEUE WITH A ZERO AND A ONE
                            for i in range(3):
                                t = Thread(target=handler.AugustGmesAAT, args =(queue, ref, args.species, protein_loc, args.threads, args.fungus, list_fasta_names, wd))
                                #AugustGmesAAT(queue, ref, species, protein_evidence, threads, fungus, list_fasta_names, wd):
                                t.daemon = True
                                t.start()
                        queue.join()
                    else:
                        print '='*70
                        print '='*70
                        sys.exit("#####UNRECOGNIZED SPECIES FOR AUGUSTUS#####\n")

                
                    
            elif args.protein_evidence != '' and args.short_reads == ''  and args.long_reads == '':  ### USING PROTEINS AND AUGUSTUS AND GMES_PETAP (TOM IMPLEMENT THE LAST) 
                check_species = ['augustus', '--species=help']
                process = subprocess.Popen(check_species, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
                out, err = process.communicate()
                protein_loc = os.path.abspath(args.protein_evidence)
                if args.species in err:
                    print '='*70
                    print '\n###RUNNING AUGUSTUS, GENEMARK-ES AND AAT###\n'
                    print '='*70
                    queue = Queue()
                    for i in range(3):
                        queue.put(i) #QUEUE WITH A ZERO AND A ONE
                        for i in range(3):
                            t = Thread(target=handler.AugustGmesAAT, args =(queue, ref, args.species, protein_loc, args.threads, args.fungus, list_fasta_names, wd))
                            #AugustGmesAAT(queue, ref, species, protein_evidence, threads, fungus, list_fasta_names, wd):
                            t.daemon = True
                            t.start()
                    queue.join()
                else:
                    print '='*70
                    print '='*70
                    sys.exit("#####UNRECOGNIZED SPECIES FOR AUGUSTUS#####\n")



            elif args.protein_evidence == '' and args.short_reads == '' and args.long_reads == '': ### USING AUGUSTUS AND GMES_PETAP
                check_species = ['augustus', '--species=help']
                process = subprocess.Popen(check_species, stdout=subprocess.PIPE,  stderr=subprocess.PIPE)
                out, err = process.communicate()
                if args.species in err:
                    print '='*70
                    print '\n###RUNNING AUGUSTUS AND GENEMARK-ES###\n'
                    print '='*70
                    queue = Queue()
                    for i in range(2):
                        queue.put(i) #QUEUE WITH A ZERO AND A ONE
                        for i in range(2):
                            t = Thread(target=handler.AugustGmes, args =(queue, ref, args.species, args.fungus, args.threads, list_fasta_names, wd))
                            t.daemon = True
                            t.start()
                    queue.join()
                else:
                    print '='*70
                    print '='*70
                    sys.exit("#####UNRECOGNIZED SPECIES FOR AUGUSTUS#####\n")

            ##Prepare EVM input files
            print '='*70
            print '\n###CONVERTING FILES###\n'
            ##HERE WE CONVERT FILES FOR EVM AND PLACE THEM IN INPUT FOLDER
            gmap_name = args.ref + '_GMAPindex'

            if args.short_reads != '' and args.protein_evidence != '': ### WE HAVE SHORT READS AND PROTEINS
                if args.no_braker:
                    augustus_file = wd + 'augustus/augustus.gff'
                    augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                    genemark_file =  wd +  'gmes/genemark.gtf'
                    genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)                
                else:
                    braker_out=wd+'braker/'+args.species+'/'
                    augustus_file = braker_out +'augustus.gff'
                    augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                    genemark_file = braker_out+'GeneMark-ET/genemark.gtf'
                    genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                mergedProtGFF3 = wd+'AAT/protein_evidence.gff3'
                trinity_path= wd + 'gmap.trinity.gff3'
                weights_dic = {'Augustus':args.augustus, 'assembler-'+args.pasa_db:args.pasa, 
                    'GeneMark.hmm':args.genemark, 'AAT':args.AAT, gmap_name :args.trinity}
                
            elif args.short_reads == '' and args.protein_evidence != '': ### WE HAVE PROTEINS BUT NOT SHORT READS
                augustus_file = wd + 'augustus/augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file =  wd +  'gmes/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                mergedProtGFF3 = wd+'AAT/protein_evidence.gff3'
                weights_dic = {'Augustus':args.augustus, 'GeneMark.hmm':args.genemark, 'AAT':args.AAT}
                
            elif args.short_reads != '' and args.protein_evidence == '': ### WE HAVE SHORT READS BUT NOT PROTEINS
                braker_out=wd+'braker/'+args.species+'/'
                augustus_file = braker_out +'augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file = braker_out+'GeneMark-ET/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                trinity_path= wd + 'gmap.trinity.gff3'
                weights_dic = {'Augustus':args.augustus, 'assembler-'+args.pasa_db:args.pasa, 
                    'GeneMark.hmm':args.genemark, gmap_name :args.trinity}
                
            elif args.short_reads == '' and args.protein_evidence == '':  ### WE NOT USE EITHER PROTEINS OR SHORT READS
                augustus_file = wd + 'augustus/augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file =  wd +  'gmes/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                weights_dic = {'Augustus':args.augustus, 'GeneMark.hmm':args.genemark}

            
            print '='*70
            print '='*70
            ## HERE WE RUN EVM; WE PREPARE FILES THAT ARE REQUIRED BY EVM LIKE WEIGTH TABLE
            if not args.no_EVM:
                print'\n###PREPARE EVM INPUTS###\n'
                
                print'##CREATING EVM DIRECTORY\n'

                evm_dir = wd+'evm_inputs/'
                logistic.check_create_dir(evm_dir)
                #print '> EVM input directory created in ' + evm_dir
                print '\n##PREPARING THE FILES\n'
                if args.protein_evidence != '' and args.short_reads != '' :  ### WE HAVE SHORT READS AND PROTEINS
                    evm_inputs = {'pasa':pasa_gff3, 'augustus':augustus_gff3, 'genemark':genemark_gff3, 'AAT':mergedProtGFF3, 'gmap': trinity_path}
                elif  args.protein_evidence == '' and args.short_reads != '' :  ### WE HAVE SHORT READS BUT NOT PROTEINS
                    evm_inputs = {'pasa':pasa_gff3, 
                            'augustus':augustus_gff3, 'genemark':genemark_gff3,'gmap': trinity_path}
                elif  args.protein_evidence != '' and args.short_reads == '' : ### WE HAVE PROTEINS BUT NOT SHORT READS
                    evm_inputs = { 'augustus':augustus_gff3, 'genemark':genemark_gff3, 'AAT':mergedProtGFF3}
                elif  args.protein_evidence == '' and args.short_reads == '' : ### WE NOT USE EITHER PROTEINS OR SHORT READS
                    evm_inputs = { 'augustus':augustus_gff3, 'genemark':genemark_gff3}
                
                list_soft, pred_file, transcript_file, protein_file= inputEvm.group_EVM_inputs(evm_dir, evm_inputs)
                
                print'##PREPARING WEIGHTS FILE\n'
                pasa_name = 'assembler-'+args.pasa_db
                weight_file = inputEvm.evm_weight(evm_dir, weights_dic, list_soft, pasa_name, gmap_name)
                print '='*70
                
                ###EVM PIPELINE
                print '='*70
                
                print'\n###RUNNING EVM PIPELINE###\n'
                if args.protein_evidence != '' and args.short_reads != '' :  ### WE HAVE SHORT READS AND PROTEINS
                    evm_gff3 = evm_pipeline.evm_pipeline(wd, args.threads, ref, weight_file, pred_file, transcript_file, 
                                                    protein_file, args.segmentSize, args.overlapSize)
                elif  args.protein_evidence == '' and args.short_reads != '' :  ### WE HAVE SHORT READS BUT NOT PROTEINS
                    protein_file = ''
                    evm_gff3 = evm_pipeline.evm_pipeline(wd, args.threads, ref, weight_file, pred_file, transcript_file, 
                                                    protein_file, args.segmentSize, args.overlapSize)
                elif  args.protein_evidence != '' and args.short_reads == '' : ### WE HAVE PROTEINS BUT NOT SHORT READS
                    transcript_file = ''
                    evm_gff3 = evm_pipeline.evm_pipeline(wd, args.threads, ref, weight_file, pred_file, transcript_file, 
                                                    protein_file, args.segmentSize, args.overlapSize)
                elif  args.protein_evidence == '' and args.short_reads == '' : ### WE NOT USE EITHER PROTEINS OR SHORT READS
                    protein_file = ''
                    transcript_file = ''
                    evm_gff3 = evm_pipeline.evm_pipeline(wd, args.threads, ref, weight_file, pred_file, transcript_file, 
                                                    protein_file, args.segmentSize, args.overlapSize)
                
                
                #print '\n>>> EVM pipeline finished!!\n'
                
                ##KEEP THIS OUTPUT
                FinalFiles.append(evm_gff3)
                print '='*70 
                
                if args.short_reads == ''  and args.long_reads == '':
                    sys.exit("##### EVM IS FINISHED #####\n")
        
                
            
                
            else:
                print "NOT RUNNING EVM PIPELINE"
                evm_gff3 = wd + '/evm_output/evm.out.combined.gff3'
                print '='*70        
                print '='*70  

                ###RE-RUN PASA PIPELINE
                
            print '='*70 + '\n'
            ##HERE WE CAN EXCLUDE TO RUN AGAIN PASA TO UPDATE THE DATABASE AFTER EVM; #We only want to update if it ran with short reads
            
            if args.short_reads != "" and not args.no_update:
                
                    
                print'\n###UPDATE WITH PASA DATABASE###\n'
                
                #for round_n in range(1,3): #Two rounds, 1 & 2
                round_n = "1"    
                firstRound = pasa_dir + 'annotation.PASAupdated.round1.gff3'
                if os.path.isfile(firstRound): 
                    print ('Update already performed --- skipping')
                    updatedGff3 = firstRound
                elif args.long_reads == "":
                    trinity_evm =  wd + 'trinity_evm_combined.gff3'
                    trinity_evm = logistic.catTwoFiles(trinityGFF3, evm_gff3, trinity_evm)
                    updatedGff3 = evm_pipeline.update_database(args.threads ,  round_n, pasa_dir, args.pasa_db, align_pasa_conf, ref, 
                                trinity_out, trinity_evm)
                else:
                    print '\n##UPDATE ROUND ###\n'
                    updatedGff3 = evm_pipeline.update_database(args.threads ,  round_n, pasa_dir, args.pasa_db, align_pasa_conf, ref, 
                                trinity_out, evm_gff3)
                    
                ##Keep this output
                FinalFiles.append(firstRound)
                print '='*70
            else:
                updatedGff3 = evm_gff3
                
            #updatedGff3 = wd+'PASA/annotation.PASAupdated.round1.gff3'      
            ##HERE WE CHECK IF WE HAVE LONG READS; IF LONG READS ARE NOT PROVIDED, THE SOFTWARE STOPS
            if args.long_reads == '':
                sys.exit("#####ANNOTATION FINISHED WITHOUT USING LONG READS#####\n")
            ## HERE WE START WITH LONG READS
            else:
                print '='*70
                print '\n###CONSENSUS###\n'
                if args.long_reads and not args.no_consensus: 
                    #Means there are long reads to map and user wants to run this pipeline
                    print "\n###CONSENSUS PIPELINE###\n"
                    consensus_wd = (wd+'consensus/')   
                    logistic.check_create_dir(consensus_wd)
                    print "\n#MAPPING TO GFF3\n"
            ##HERE WE MAP THE READS ON THE GENOME USING GMAP
                    long_sam = mapping.gmap('sam', ref, long_fastq, args.threads, 'samse', args.min_intron_length, args.max_intron_length, args.H, wd, Fflag = False ) ## change in 1 and 2
                    long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, args.threads, wd)
                    
            ##HERE WE MERGE THE GMAP OUTPUT WITH THE EVM OUTPUT TO HAVE ONE FILE
                    mergedmapGFF3 = consensus_wd + 'mergedGmapEvm.beforeAssembly.gff3'
            ##HERE WE CHECK IF WE HAVE THE PASA UPDATED FILE OR THE EVM ORIGINAL FILE        
                    if os.path.isfile(updatedGff3):
                        ##HERE WE MERGE THE TWO FILES
                        logistic.catTwoBeds(long_sorted_bam, updatedGff3, mergedmapGFF3)
                    else:
                        logistic.catTwoBeds(long_sorted_bam, evm_gff3, mergedmapGFF3)
                    print "\n#GFFREAD\n"
                    ##HERE WE TRANSFORM THE COODINATES INTO SEQUENCES USING THE REFERENCE
                    gffreadFastaFile = consensus.gffread(mergedmapGFF3, ref, consensus_wd)
                    ##HERE WE STORE THE SEQUENCE IN A DICTIONARY
                    gffreadDict = consensus.fasta2Dict(gffreadFastaFile)
                    print "\n#CLUSTERING\n"
                    ##HERE WE CLUSTER THE SEQUENCES BASED ON THE GENOME POSITION
                    cluster_list = consensus.cluster_pipeline(mergedmapGFF3, args.assembly_overlapLength, args.stranded, consensus_wd)
                    print "\n#CONSENSUS FOR EACH CLUSTER\n"
                    ##HERE WE MAKE CONSENSUS FOR EACH CLUSTER
                    tmp_wd = consensus_wd+'tmp/'
                    logistic.check_create_dir(tmp_wd)
                    tmp_assembly_file = tmp_wd + 'assembly.fasta'
                    if os.path.isfile(tmp_assembly_file):
                        print 'No assembly'
                    else:
                        consensus.assembly(cluster_list, gffreadDict, args.cluster_min_evidence, args.cluster_max_evidence, 
                        args.assembly_overlapLength, args.assembly_percentIdentity, args.overhang , args.threads, tmp_wd)
                        utrs.lengthSupport(tmp_wd)
                    
        ## WITH THE ELSE, WE ALLOW THE USER TO DECIDE TO CHANGE THE ASSEMBLY PARAMETERS AND COLLECT DIFFERENT ASSEMBLED SEQUENCES WITHOT RUNNING THE FULL PIPELINE
        else:
            print '='*70     
        
            print "\n###COLLECT ONLY SEQUENCES###\n"
            
            ##PLACE WHERE THE EVM ANNOTATION IS LOCATED
            updatedGff3 = wd+'PASA/annotation.PASAupdated.round1.gff3'
            consensus_wd = (wd+'consensus/')
            tmp_wd = consensus_wd+'tmp/'
            if not os.path.isfile(updatedGff3):
                updatedGff3 =  wd + '/evm_output/evm.out.combined.gff3'
            mergedmapGFF3 = consensus_wd + 'mergedGmapEvm.beforeAssembly.gff3'
            gffreadFastaFile = consensus.gffread(mergedmapGFF3, ref, consensus_wd)
        ##HERE WE COLLECT THE ASSEMBLED SEQUENCES. WE COLLWCT ONLY SEQUENCE THAT PASS THE FILTER
        evm_nosupport =  collect.parseOnly(args.assembly_readThreshold, args.only_unitigs, wd)
        tmp_assembly = collect.catAssembled(wd)
        ##HERE WE COLLECT THE NEW ASSEMBLED SEQUENCES AND WE COLLECT THE OLD EVM DATA
        mergedFastaFilename = consensus_wd+'assembly.wEVM.fasta'
        collect.addEVM(gffreadFastaFile, tmp_assembly, args.only_unitigs , evm_nosupport, mergedFastaFilename)

        FinalFiles.append(mergedFastaFilename)
        #shutil.rmtree(tmp_wd)
        
        print '='*70     
        print '='*70        

        print "\n###MAPPING CONSENSUS ASSEMBLY###\n"
        ##HERE WE MAP ALL THE FASTA FILES TO THE GENOME USING GMAP
        consensusMappedGFF3 = mapping.gmap('cons', ref, mergedFastaFilename, args.threads, 'gff3_gene', args.min_intron_length, args.max_intron_length, args.H, wd, Fflag = True )
        FinalFiles.append(consensusMappedGFF3)
        
        print '='*70     
        print '='*70        

        print "\n###GETTING THE STRAND RIGHT###\n"
        ##IN THIS STEP WE CORRECT FOR STRAND. GMAP CAN NOT DECIDE THE STRAND FOR SINGLE EXONS GENE MODELS. WE USE THE ORIENTATION FROM EVM IF GMAP INVERT THE ORIGINAL STRAND
        outputList_gmap_multi = grs.gffread_multiexons(consensusMappedGFF3, multiExonFlag = True)
        outputList_gmap_all = grs.gffread_multiexons(consensusMappedGFF3)
        if os.path.isfile(updatedGff3):
            outputList_pasa_all = grs.gffread_multiexons(updatedGff3)
        else:
            outputList_pasa_all = grs.gffread_multiexons(evm_gff3)
        gmap_dict_multi = grs.gffread_parser(outputList_gmap_multi)
        gmap_dict_all = grs.gffread_parser(outputList_gmap_all)
        pasa_dict_all = grs.gffread_parser(outputList_pasa_all)
        
        gmapOut, pasaOut = grs.compare_dicts(gmap_dict_multi, gmap_dict_all, pasa_dict_all)
        finalOutput = grs.combineGff3(gmapOut, pasaOut, outputList_gmap_all, outputList_pasa_all, wd)
        ##HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN PASA TO CORRECT SMALL ERRORS
        fastaAll = logistic.catTwoFasta(trinity_out, mergedFastaFilename, pasa_dir)
        final = evm_pipeline.update_database(args.threads, "2", pasa_dir, args.pasa_db, align_pasa_conf, ref, fastaAll, finalOutput)
        FinalFiles.append(final) 
        FinalFiles.append(finalOutput)
        
        print '='*70        

            
        print'\n###CREATING OUTPUT DIRECTORY###\n'
        final_output_dir = wd+'output/'
        logistic.check_create_dir(final_output_dir)
        print "\n###DONE###\n "
        
        
        print "\n##PLACING OUTPUT FILES IN OUTPUT DIRECTORY"
        for filename in FinalFiles:
            if filename != '':
                logistic.copy_file(filename, final_output_dir)
        if args.keep_tmp:
            dirs_list = ['/PASA/', 'augustus/', 'gmes/', 'AAT/',  'split/']
            for dirs in dirs_list:
                dest = wd + dirs
                shutil.rmtree(dest, ignore_errors=True)
        
        print '\n\n\n###############\n###FINISHED###\n###############\n\n'
            
            
    
    
    
 
if __name__ == '__main__':
    main()
