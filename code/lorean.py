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

# LIBRARIES

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


# OTHER SCRIPTS
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
import parseGff3 as parsegff3

####################################
### CHEKS BEFORE START LOREAN ######
####################################


os.system('sudo usermod -d /var/lib/mysql/ mysql')
os.system('sudo /etc/init.d/mysql start')
os.system('mysql --user="root" --password="lorean" --execute="set global sql_mode=\'STRICT_TRANS_TABLES,ERROR_FOR_DIVISION_BY_ZERO,NO_AUTO_CREATE_USER,NO_ENGINE_SUBSTITUTION\';"')


###############
###FUNCTIONS###
###############


def arguments():
    '''Parses the arguments from the program invocation'''

    # Call the argument parse
    parser = argparse.ArgumentParser(
        prog='lorean',
        usage='%(prog)s [options] protein_sequences reference species_name',
        description='LoReAn - Automated genome annotation pipeline that integrates long reads',
        epilog='Luigi Faino - March 2017')

    parser.add_argument('protein_evidence',
                        help="Path to protein sequences FASTA file []")
    parser.add_argument("ref",
                        help="Path to reference file")
    parser.add_argument(
        "species",
        help="Species name for AUGUSTUS training. No re-training if species already present in AUGUSTUS config folder")
    parser.add_argument("--stranded",
                        help="Run LoReAn on stranded mode [FALSE]",
                        action='store_true')
    parser.add_argument("--fungus",
                        help="Use this option for fungal species (used in Gene Mark-ES)  [FALSE]",
                        action='store_true')
    parser.add_argument(
        "--only_unitigs",
        help="Removes gene models that are not supported by long reads [FALSE]",
        action='store_true')
    # TO CHECK WITH JOSE HERE
    parser.add_argument("--keep_tmp",
                        help="Keep temporary files [FALSE]",
                        action='store_true')
    parser.add_argument(
        '--short_reads',
        nargs="?",
        default="",
        help="Path to short reads FASTQ. If paired end, comma-separated (1-1.fq,1-2.fq). BAM sorted files are allowed; the extension of the file should be filename.sorted.bam []",
        metavar='FASTQ_file')
    parser.add_argument('--long_reads', nargs="?", default="",
                        help="Path to long reads FASTQ []",
                        metavar='FASTQ_file')
    parser.add_argument(
        '--max_long_read',
        nargs="?",
        default=20000,
        help="Filter out long reads longer than this value (longer reads may affect mapping and assembling) [20000]",
        type=int)
    parser.add_argument("--pasa_db", nargs="?", default="annotation",
                        help="PASA database name [pipeline_run]")
    parser.add_argument(
        "--prefix_gene",
        nargs="?",
        default="species",
        help="Prefix to add to the final Gff3 gene name [specie]")
    parser.add_argument(
        '-wd',
        "--working_dir",
        nargs="?",
        default="annotation",
        help="Working directory (will create if not present) [./]")
    parser.add_argument('-t', "--threads", nargs="?", default="1",
                        help="Number of threads [1]",
                        metavar='N')
    parser.add_argument(
        "--overhang",
        nargs="?",
        default="20",
        help="CAP3 max overhang percent length; this value should be > 3 [20]",
        metavar='N')
    parser.add_argument(
        "--augustus",
        nargs="?",
        default="1",
        help="Weight assigned to AUGUSTUS evidence for EVM [1]",
        metavar='N')
    parser.add_argument(
        "--genemark",
        nargs="?",
        default="1",
        help="Weight assigned to GENEMARK evidence for EVM [1]",
        metavar='N')
    parser.add_argument(
        "--trinity",
        nargs="?",
        default="1",
        help="Weight assigned to Trinity mapped with GMAP evidence for EVM [1]",
        metavar='N')
    parser.add_argument("--pasa", nargs="?", default="5",
                        help="Weight assigned to PASA evidence for EVM [5]",
                        metavar='N')
    parser.add_argument(
        "--AAT",
        nargs="?",
        default="1",
        help="Weight assigned to AAT protein evidence for EVM [1]",
        metavar='N')
    parser.add_argument("--segmentSize", nargs="?", default="100000",
                        help="Segment size for EVM partitions [100000]",
                        metavar='N')
    parser.add_argument("--overlapSize", nargs="?", default="10000",
                        help="Overlap size for EVM partitions [10000]",
                        metavar='N')
    parser.add_argument("--min_intron_length", nargs="?", default="9",
                        help="Minimal intron length for GMAP [9]",
                        metavar='N')
    parser.add_argument(
        "--max_intron_length",
        nargs="?",
        default="1000",
        help="Maximal intron length for GMAP, STAR and TRINITY [1000]",
        metavar='N')
    parser.add_argument("-H", nargs="?", default="20",
                        help="Minimal length for end exon with GMAP [20]",
                        metavar='N')
    parser.add_argument("--cluster_min_evidence", nargs="?", default="5",
                        help="Minimal evidence needed to form a cluster [5]",
                        metavar='N')
    parser.add_argument(
        "--cluster_max_evidence",
        nargs="?",
        default="5000",
        help="Maximal evidence to form a cluster.Prevents the clustering or rRNA genes i.e. [5000]",
        metavar='N')
    parser.add_argument(
        "--assembly_overlapLength",
        nargs="?",
        default="200",
        help="Minimal length (in nt) of overlap for ASSEMBLY [200]",
        metavar='N')
    parser.add_argument("--assembly_percentIdentity", nargs="?", default="97",
                        help="Minimal identity for the ASSEMBLY (95-100) [97]",
                        metavar='N')
    parser.add_argument(
        "--assembly_readThreshold",
        nargs="?",
        default="0.3",
        help="Fraction of reads supporting an assembled UNITIG to keep on the ASSEMBLY (0.1-1) [0.3]",
        metavar='F')
    parser.add_argument("--no_EVM",
                        help="Run until the preparation of EVM inputs [FALSE]",
                        action='store_true')
    parser.add_argument(
        "--no_consensus",
        help="Do not run the long reads consensus pipeline [FALSE]",
        action='store_true')
    parser.add_argument("--no_update",
                        help="Do not run the PASA update[FALSE]",
                        action='store_true')
    parser.add_argument("--collect_only",
                        help="Collect only assebmled transcripts [FALSE]",
                        action='store_true')
    #parser.add_argument("--SS_lib_type",
                        #help="For strand specific short reads RNA-Seq set the orientation (see Trinity manual):
                        #if paired options are RF or FR; 
                        #if single options are F or R (dUTP method = RF)",
                        #nargs="?")


    # Parse and set defaults for optional arguments
    args = parser.parse_args()
    return args
###############
###MAIN###
###############


def main():

    if os.path.isfile("/data/gm_key"):
        os.system('cp /data/gm_key /home/lorean/.gm_key')
        '''Core of the program'''

        # Parse the arguments
        args = arguments()

        # Useful variables for later
        wd_base = os.path.abspath(args.working_dir) + '/'
        wd = wd_base + 'run/'
        wd_gmap = wd
        ref = os.path.abspath(args.ref)

        FinalFiles = []  # STORE THE IMPORTANT OUTPUT FILES

        # Arrange stuff before start
        logistic.check_create_dir(wd)
        logistic.check_file(ref)
        gmap_wd = wd + '/gmap_output/'
        logistic.check_create_dir(gmap_wd)

        # COLLECT ONLY ONLY RUNS PART OF THE CONSENSUS PIPELINE
        if not args.collect_only:
            list_fasta_names = multiple.single_fasta(ref, wd)
            if args.short_reads != '' or args.long_reads != '':

                print '\n###STAR MAPPING###\n'
                # SHORT READS
                if 'fastq' in args.short_reads or 'fq' in args.short_reads:
                    star_out = wd + '/STAR/'
                    logistic.check_create_dir(star_out)
                    if ',' in args.short_reads:
                        pairedEndFiles = args.short_reads.split(',')
                        short_1 = os.path.abspath(pairedEndFiles[0])
                        short_2 = os.path.abspath(pairedEndFiles[1])
                        short_reads_file = [short_1, short_2]

                    else:
                        short_reads_file = os.path.abspath(args.short_reads)

                    # Map with STAR
                    short_bam = mapping.star(
                        ref,
                        short_reads_file,
                        args.threads,
                        args.max_intron_length,
                        star_out)
                    short_sorted_bam = mapping.samtools_sort(
                        short_bam, args.threads, wd)

                    # Keep the output
                    FinalFiles.append(short_sorted_bam)
                # BAM SORTED FILES GET IN HERE
                elif 'bam' in args.short_reads:
                    star_out = wd + '/STAR/'
                    logistic.check_create_dir(star_out)
                    short_sorted_bam = os.path.abspath(args.short_reads)
                    cmdstring = "mv %s %s" % (short_sorted_bam, star_out)
                    os.system(cmdstring)
                    bam_file = args.short_reads.split("/")
                    short_bam = star_out + "/" + bam_file[2]
                    short_sorted_bam = mapping.samtools_sort(
                        short_bam, args.threads, wd)
                    # print short_sorted_bam
                else:
                    short_sorted_bam = False
                    print 'No short reads file'
                # LONG READS
                if 'fastq' in args.long_reads or 'fq' in args.long_reads or 'fasta' in args.long_reads or 'fa' in args.long_reads  :
                    # with this operation, reads are filtered for their length.
                    # Nanopore reads can be chimaras or sequencing artefacts.
                    # filtering on length reduces the amount of sequencing
                    # artefacts

                    print "\n###FILTERING OUT LONG READS###\n"

                    long_fasta, filter_count = mapping.filterLongReads(
                        args.long_reads, args.assembly_overlapLength, args.max_long_read, gmap_wd , a = True)

                    if filter_count != 0:
                        print "LOREAN KEPT " + str(filter_count) + " READS AFTER LENGTH FILTERING\n"
                    if not short_sorted_bam:
                        # If short reads have been mapped dont do it

                        print '\n###GMAP###\n'

                        long_sam = mapping.gmap(
                            'sam',
                            ref,
                            long_fasta,
                            args.threads,
                            'samse',
                            args.min_intron_length,
                            args.max_intron_length,
                            args.H,
                            gmap_wd,
                            Fflag=False)
                        

                        # Convert to sorted BAM
                        long_sorted_bam = mapping.sam_to_sorted_bam(
                            long_sam, args.threads, wd)

                        # Keep the output
                        FinalFiles.append(long_sorted_bam)
                    else:
                        long_sorted_bam = False

                else:

                    print '\n###NO LONG READS FILE###\n'
                    long_sorted_bam = False

                if short_sorted_bam:  # If there are short reads, these will serve to the transcript assembly pipeline
                    default_bam = short_sorted_bam
                else:
                    default_bam = long_sorted_bam

                # TRANSCRIPT ASSEMBLY

                # TRINITY

                print '\n###TRINITY###\n'

                trin_dir = wd + 'Trinity/'
                logistic.check_create_dir(trin_dir)

                trinity_out = transcripts.trinity(
                    default_bam, trin_dir, args.max_intron_length, args.threads)

                # PASA Pipeline

                print '\n###PASA###\n'

                # Create PASA folder and configuration file
                pasa_dir = wd + 'PASA/'
                logistic.check_create_dir(pasa_dir)
                align_pasa_conf = transcripts.pasa_configuration(
                    pasa_dir, args.pasa_db)
                # Launch PASA
                pasa_gff3 = transcripts.pasa_call(
                    pasa_dir,
                    align_pasa_conf,
                    args.pasa_db,
                    ref,
                    trinity_out,
                    args.max_intron_length,
                    args.threads)
                # HERE WE PARALLELIZE PROCESSES WHEN MULTIPLE THREADS ARE USED

                if args.short_reads != '' or args.long_reads != '':
                    check_species = ['augustus', '--species=help']
                    process = subprocess.Popen(
                        check_species, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = process.communicate()
                    protein_loc = os.path.abspath(args.protein_evidence)
                    if args.species in err:

                        print '\n###RUNNING AUGUSTUS, GENEMARK-ES AND AAT###\n'

                        queue = Queue()
                        for i in range(3):
                            queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                            for i in range(3):
                                t = Thread(
                                    target=handler.AugustGmesAAT,
                                    args=(
                                        queue,
                                        ref,
                                        args.species,
                                        protein_loc,
                                        args.threads,
                                        args.fungus,
                                        list_fasta_names,
                                        wd))
                                # AugustGmesAAT(queue, ref, species,
                                # protein_evidence, threads, fungus,
                                # list_fasta_names, wd):
                                t.daemon = True
                                t.start()
                        queue.join()
                    elif args.short_reads:  # USING PROTEINS AND SHORT READS
                        queue = Queue()

                        print '\n###RUNNING BRAKER1 AND AAT###\n'

                        for i in range(2):
                            queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                        for i in range(2):
                            t = Thread(
                                target=handler.BrakerAAT,
                                args=(
                                    queue,
                                    ref,
                                    default_bam,
                                    args.species,
                                    args.protein_evidence,
                                    args.threads,
                                    args.fungus,
                                    list_fasta_names,
                                    wd))
                            t.daemon = True
                            t.start()
                        queue.join()

                    elif args.long_reads:  # USING PROTEINS AND LONG READS
                        queue = Queue()

                        print '\n###RUNNING BRAKER1 AND AAT###\n'

                        for i in range(2):
                            queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                        for i in range(2):
                            t = Thread(
                                target=handler.BrakerAAT,
                                args=(
                                    queue,
                                    ref,
                                    long_sorted_bam,
                                    args.species,
                                    args.protein_evidence,
                                    args.threads,
                                    args.fungus,
                                    list_fasta_names,
                                    wd))
                            t.daemon = True
                            t.start()
                        queue.join()

                print '\n###GMAP###\n'

                trinityGFF3 = mapping.gmap(
                    'trin',
                    ref,
                    trinity_out,
                    args.threads,
                    'gff3_gene',
                    args.min_intron_length,
                    args.max_intron_length,
                    args.H,
                    gmap_wd,
                    Fflag=True)

            # USING PROTEINS AND AUGUSTUS AND GMES_PETAP (TOM IMPLEMENT THE
            # LAST)
            elif args.short_reads == '' and args.long_reads == '':
                check_species = ['augustus', '--species=help']
                process = subprocess.Popen(
                    check_species,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
                out, err = process.communicate()
                protein_loc = os.path.abspath(args.protein_evidence)
                if args.species in err:

                    print '\n###RUNNING AUGUSTUS, GENEMARK-ES AND AAT###\n'

                    queue = Queue()
                    for i in range(3):
                        queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                        for i in range(3):
                            t = Thread(
                                target=handler.AugustGmesAAT,
                                args=(
                                    queue,
                                    ref,
                                    args.species,
                                    protein_loc,
                                    args.threads,
                                    args.fungus,
                                    list_fasta_names,
                                    wd))
                            # AugustGmesAAT(queue, ref, species,
                            # protein_evidence, threads, fungus,
                            # list_fasta_names, wd):
                            t.daemon = True
                            t.start()
                    queue.join()
                else:

                    sys.exit("#####UNRECOGNIZED SPECIES FOR AUGUSTUS#####\n")

            # Prepare EVM input files

            print '\n###RUNNING EVM###\n'

            # HERE WE CONVERT FILES FOR EVM AND PLACE THEM IN INPUT FOLDER
            gmap_name = args.ref + '_GMAPindex'
            if args.short_reads != '' or args.long_reads != '':  # WE HAVE SHORT READS AND PROTEINS
                braker_out = wd + 'braker/' + args.species + '/'
                if os.path.exists(braker_out):
                    augustus_file = braker_out + 'augustus.gff'
                    augustus_gff3 = inputEvm.convert_augustus(
                        augustus_file, wd)
                    genemark_file = braker_out + 'GeneMark-ET/genemark.gtf'
                    genemark_gff3 = inputEvm.convert_genemark(
                        genemark_file, wd)
                else:
                    augustus_file = wd + 'augustus/augustus.gff'
                    augustus_gff3 = inputEvm.convert_augustus(
                        augustus_file, wd)
                    genemark_file = wd + 'gmes/genemark.gtf'
                    genemark_gff3 = inputEvm.convert_genemark(
                        genemark_file, wd)

                mergedProtGFF3 = wd + 'AAT/protein_evidence.gff3'
                trinity_path = gmap_wd + 'gmap.trinity.gff3'
                weights_dic = {
                    'Augustus': args.augustus,
                    'assembler-' +
                    args.pasa_db: args.pasa,
                    'GeneMark.hmm': args.genemark,
                    'AAT': args.AAT,
                    gmap_name: args.trinity}

            elif args.short_reads == '' and args.long_reads == '':  # WE HAVE PROTEINS BUT NOT SHORT READS
                augustus_file = wd + 'augustus/augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file = wd + 'gmes/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                mergedProtGFF3 = wd + 'AAT/protein_evidence.gff3'
                weights_dic = {
                    'Augustus': args.augustus,
                    'GeneMark.hmm': args.genemark,
                    'AAT': args.AAT}

            # HERE WE RUN EVM; WE PREPARE FILES THAT ARE REQUIRED BY EVM LIKE
            # WEIGTH TABLE
            if not args.no_EVM:
                evm_dir = wd + 'evm_inputs/'
                logistic.check_create_dir(evm_dir)
                # print '> EVM input directory created in ' + evm_dir
                if args.short_reads != '' or args.long_reads != '':  # WE HAVE SHORT READS AND PROTEINS
                    evm_inputs = {
                        'pasa': pasa_gff3,
                        'augustus': augustus_gff3,
                        'genemark': genemark_gff3,
                        'AAT': mergedProtGFF3,
                        'gmap': trinity_path}
                elif args.short_reads == '' and args.long_reads == '':  # WE HAVE PROTEINS BUT NOT SHORT READS
                    evm_inputs = {
                        'augustus': augustus_gff3,
                        'genemark': genemark_gff3,
                        'AAT': mergedProtGFF3}

                list_soft, pred_file, transcript_file, protein_file = inputEvm.group_EVM_inputs(
                    evm_dir, evm_inputs)

                pasa_name = 'assembler-' + args.pasa_db
                weight_file = inputEvm.evm_weight(
                    evm_dir, weights_dic, list_soft, pasa_name, gmap_name)

                # EVM PIPELINE

                if args.short_reads != '' or args.long_reads != '':  # WE HAVE SHORT READS AND PROTEINS
                    evm_gff3 = evm_pipeline.evm_pipeline(
                        wd,
                        args.threads,
                        ref,
                        weight_file,
                        pred_file,
                        transcript_file,
                        protein_file,
                        args.segmentSize,
                        args.overlapSize)
                elif args.short_reads == '' and args.long_reads == '':  # WE HAVE PROTEINS BUT NOT SHORT READS
                    transcript_file = ''
                    evm_gff3 = evm_pipeline.evm_pipeline(
                        wd,
                        args.threads,
                        ref,
                        weight_file,
                        pred_file,
                        transcript_file,
                        protein_file,
                        args.segmentSize,
                        args.overlapSize)

                # print '\n>>> EVM pipeline finished!!\n'

                # KEEP THIS OUTPUT
                FinalFiles.append(evm_gff3)

                if args.short_reads == '' and args.long_reads == '':
                    sys.exit("##### EVM IS FINISHED #####\n")

            else:
                print "\n###NOT RUNNING EVM PIPELINE###\n"
                evm_gff3 = wd + '/evm_output/evm.out.combined.gff3'

                # RE-RUN PASA PIPELINE

            # HERE WE CAN EXCLUDE TO RUN AGAIN PASA TO UPDATE THE DATABASE
            # AFTER EVM; #We only want to update if it ran with short reads
            round_n = 0
            if (args.short_reads != "" and not args.no_update) and (args.long_reads == "" and not args.no_update):
                print'\n###UPDATE WITH PASA DATABASE###\n'

                # for round_n in range(1,3): #Two rounds, 1 & 2
                
                firstRound = pasa_dir + 'annotation.PASAupdated.round1.gff3'
                if os.path.isfile(firstRound):

                    print ('UPDATE ALREADY PERFORMED --- skipping')
                    updatedGff3 = firstRound
                    round_n = 1
                else:
                    round_n += 1
                    print '\n##UPDATE ROUND ###\n'
                    if args.long_reads == "":
                        finalOutput = evm_pipeline.update_database(
                            args.threads,
                            str(round_n),
                            pasa_dir,
                            args.pasa_db,
                            align_pasa_conf,
                            ref,
                            trinity_out,
                            evm_gff3,
                            "a")
                        final = parsegff3.genename(finalOutput, args.prefix_gene)
                        updatedGff3 = grs.newNames(final)

            else:
                updatedGff3 = evm_gff3

            #updatedGff3 = wd+'PASA/annotation.PASAupdated.round1.gff3'
            # HERE WE CHECK IF WE HAVE LONG READS; IF LONG READS ARE NOT
            # PROVIDED, THE SOFTWARE STOPS
            if args.long_reads == '':
                final_output_dir = wd + 'output/'
                logistic.check_create_dir(final_output_dir)
                for filename in FinalFiles:
                    if filename != '':
                        logistic.copy_file(filename, final_output_dir)
                cmdstring = "chmod -R 775 %s" % (wd)
                os.system(cmdstring)

                sys.exit(
                    "#####ANNOTATION FINISHED WITHOUT USING LONG READS#####\n")

            # HERE WE START WITH LONG READS
            else:

                print '\n###RUNNING iASSEMBLER###\n'

                if args.long_reads and not args.no_consensus:
                    # Means there are long reads to map and user wants to run
                    # this pipeline

                    consensus_wd = (wd + 'consensus/')
                    logistic.check_create_dir(consensus_wd)
                    # HERE WE MAP THE READS ON THE GENOME USING GMAP

                    if not long_sorted_bam:
                        long_sam = mapping.gmap(
                            'sam',
                            ref,
                            long_fasta,
                            args.threads,
                            'samse',
                            args.min_intron_length,
                            args.max_intron_length,
                            args.H,
                            gmap_wd,
                            Fflag=False)  # change in 1 and 2
                        long_sorted_bam = mapping.sam_to_sorted_bam(
                            long_sam, args.threads, wd)

            # HERE WE MERGE THE GMAP OUTPUT WITH THE EVM OUTPUT TO HAVE ONE
            # FILE
                    fileName = consensus_wd + 'mergedGmapEvm.beforeAssembly.gff3'
            # HERE WE CHECK IF WE HAVE THE PASA UPDATED FILE OR THE EVM
            # ORIGINAL FILE
                    if os.path.isfile(updatedGff3):
                        # HERE WE MERGE THE TWO FILES
                        mergedmapGFF3 = logistic.catTwoBeds(
                            long_sorted_bam, updatedGff3, fileName)
                    else:
                        mergedmapGFF3 = logistic.catTwoBeds(
                            long_sorted_bam, evm_gff3, fileName)

                    print "\n\t#GFFREAD\n"

                    # HERE WE TRANSFORM THE COODINATES INTO SEQUENCES USING THE
                    # REFERENCE
                    gffreadFastaFile = consensus.gffread(
                        mergedmapGFF3, ref, consensus_wd)
                    # HERE WE STORE THE SEQUENCE IN A DICTIONARY
                    
                    long_fasta, filter_count = mapping.filterLongReads(
                        gffreadFastaFile, args.assembly_overlapLength, args.max_long_read, consensus_wd, a = False)

                    gffreadDict = consensus.fasta2Dict(gffreadFastaFile)

                    print "\n\t#CLUSTERING\n"

                    # HERE WE CLUSTER THE SEQUENCES BASED ON THE GENOME
                    # POSITION
                    cluster_list = consensus.cluster_pipeline(
                        mergedmapGFF3, args.assembly_overlapLength, args.stranded, consensus_wd)

                    print "\n\t#CONSENSUS FOR EACH CLUSTER\n"

                    # HERE WE MAKE CONSENSUS FOR EACH CLUSTER
                    tmp_wd = consensus_wd + 'tmp/'
                    logistic.check_create_dir(tmp_wd)
                    tmp_assembly_file = tmp_wd + 'assembly.fasta'
                    if os.path.isfile(tmp_assembly_file):
                        print 'No assembly'
                    else:
                        consensus.assembly(
                            cluster_list,
                            gffreadDict,
                            args.cluster_min_evidence,
                            args.cluster_max_evidence,
                            args.assembly_overlapLength,
                            args.assembly_percentIdentity,
                            args.overhang,
                            args.threads,
                            tmp_wd)
                        utrs.lengthSupport(tmp_wd)

        # WITH THE ELSE, WE ALLOW THE USER TO DECIDE TO CHANGE THE ASSEMBLY
        # PARAMETERS AND COLLECT DIFFERENT ASSEMBLED SEQUENCES WITHOT RUNNING
        # THE FULL PIPELINE
        else:

            print "\n###COLLECT ONLY SEQUENCES###\n"

            # PLACE WHERE THE EVM ANNOTATION IS LOCATED
            updatedGff3 = wd + 'PASA/annotation.PASAupdated.round1.gff3'
            consensus_wd = (wd + 'consensus/')
            tmp_wd = consensus_wd + 'tmp/'
            if not os.path.isfile(updatedGff3):
                updatedGff3 = wd + '/evm_output/evm.out.combined.gff3'
            mergedmapGFF3 = consensus_wd + 'mergedGmapEvm.beforeAssembly.gff3'
            gffreadFastaFile = consensus.gffread(
                mergedmapGFF3, ref, consensus_wd)
        # HERE WE COLLECT THE ASSEMBLED SEQUENCES. WE COLLWCT ONLY SEQUENCE
        # THAT PASS THE FILTER
        evm_nosupport = collect.parseOnly(
            args.assembly_readThreshold, args.only_unitigs, wd)
        tmp_assembly = collect.catAssembled(wd)
        # HERE WE COLLECT THE NEW ASSEMBLED SEQUENCES AND WE COLLECT THE OLD
        # EVM DATA
        mergedFastaFilename = consensus_wd + 'assembly.wEVM.fasta'
        collect.addEVM(
            gffreadFastaFile,
            tmp_assembly,
            args.only_unitigs,
            evm_nosupport,
            mergedFastaFilename)

        # shutil.rmtree(tmp_wd)

        print "\n###MAPPING CONSENSUS ASSEMBLIES###\n"

        # HERE WE MAP ALL THE FASTA FILES TO THE GENOME USING GMAP
        consensusMappedGFF3 = mapping.gmap(
            'cons',
            ref,
            mergedFastaFilename,
            args.threads,
            'gff3_gene',
            args.min_intron_length,
            args.max_intron_length,
            args.H,
            gmap_wd,
            Fflag=True)

        print "\n###GETTING THE STRAND RIGHT###\n"

        # IN THIS STEP WE CORRECT FOR STRAND. GMAP CAN NOT DECIDE THE STRAND
        # FOR SINGLE EXONS GENE MODELS. WE USE THE ORIENTATION FROM EVM IF GMAP
        # INVERT THE ORIGINAL STRAND
        
        
        #outputList_gmap_multi = grs.gffread_multiexons(
            #consensusMappedGFF3, multiExonFlag=True)
        #outputList_gmap_all = grs.gffread_multiexons(consensusMappedGFF3)
        #if os.path.isfile(updatedGff3):
            #outputList_pasa_all = grs.gffread_multiexons(updatedGff3)
        #else:
            #outputList_pasa_all = grs.gffread_multiexons(evm_gff3)
        #gmap_dict_multi = grs.gffread_parser(outputList_gmap_multi)
        #gmap_dict_all = grs.gffread_parser(outputList_gmap_all)
        #pasa_dict_all = grs.gffread_parser(outputList_pasa_all)

        #gmapOut, pasaOut = grs.compare_dicts(
            #gmap_dict_multi, gmap_dict_all, pasa_dict_all)
        #finalOutput = grs.combineGff3(
            #gmapOut,
            #pasaOut,
            #outputList_gmap_all,
            #outputList_pasa_all,
            #gmap_wd)
            
        finalOutput = grs.strand(evm_gff3, consensusMappedGFF3, gmap_wd)
        gffPasa = grs.appendID(finalOutput)
        noOverl = grs.removeOverlap(gffPasa)
        #simplified = grs.parseGff(finalOutput)
        noDisc = grs.removeDiscrepancy(noOverl, evm_gff3)
        uniqGene = grs.newNames(noDisc)

    #out1 = strand(evmgff, gmapgff, wd)
    #out2 = appendID(out1)
    #out3 = removeOverlap(out2)
    #out4 = removeDiscrepancy(out3, evmgff)
    #newNames(out4)


        
        # HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN
        # PASA TO CORRECT SMALL ERRORS

        fastaAll = logistic.catTwoFasta(
            trinity_out, mergedFastaFilename, long_fasta, pasa_dir)
        round_n += 1
        finalupdate = evm_pipeline.update_database(
            args.threads,
            str(round_n),
            pasa_dir,
            args.pasa_db,
            align_pasa_conf,
            ref,
            fastaAll,
            uniqGene,
            "a")
        round_n += 1
        finalupdate2 = evm_pipeline.update_database(
            args.threads,
            str(round_n),
            pasa_dir,
            args.pasa_db,
            align_pasa_conf,
            ref,
            fastaAll,
            finalupdate,
            "a")
        finalupdate3 = grs.genename(finalupdate2, args.prefix_gene)


        FinalFiles.append(finalupdate3)


        print'\n###CREATING OUTPUT DIRECTORY###\n'

        final_output_dir = wd_base + 'output/'
        logistic.check_create_dir(final_output_dir)

        print "\n##PLACING OUTPUT FILES IN OUTPUT DIRECTORY"

        for filename in FinalFiles:
            if filename != '':
                logistic.copy_file(filename, final_output_dir)
                cmdstring = "chmod -R 775 %s" % (wd)
                os.system(cmdstring)
        if args.keep_tmp:
            dirs_list = ['/PASA/', 'augustus/', 'gmes/', 'AAT/', 'split/']
            for dirs in dirs_list:
                dest = wd + dirs
                shutil.rmtree(dest, ignore_errors=True)

        print '\n\n\n###############\n###FINISHED###\n###############\n\n'

    else:
        print 'Key for GeneMark-ES not found.  Please, place the GeneMark-ES key in the folder where you have your data.'
        sys.exit("#####LOREAN STOPS HERE.#####\n")

if __name__ == '__main__':
    main()
