#! /usr/bin/env python3

###############
###IMPORTS###
###############

# LIBRARIES
# TODO
# check for strand on single exon gene based on reads mapping


import datetime
import os
import sys

import collect_only as collect
import consensus_iAssembler as consensus
# OTHER SCRIPTS
import dirs_and_files as logistic
import get_right_strand as grs
import manipulateSeq as mseq
import mapping
import pasa as pasa
import reduceUTRs as utrs


###############
###MAIN###
###############


def update(args, consensus_wd, fmtdate, genome_gmap, gmap_wd, ref):

    updatedGff3 = args.update

    if 'fastq' in args.long_reads or 'fq' in args.long_reads or 'fasta' in args.long_reads or 'fa' in args.long_reads:
        # with this operation, reads are filtered for their length.
        # Nanopore reads can be chimaras or sequencing artefacts.
        # filtering on length reduces the amount of sequencing
        # artefacts
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n###FILTERING OUT LONG READS STARTED AT:\t" + now + "\t###\n"))
        long_fasta, filter_count = mseq.filterLongReads(args.long_reads, args.assembly_overlapLength,
                                                        args.max_long_read, gmap_wd, args.adapter, args.threads,
                                                        a=True)
        if filter_count != 0:
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(("###FINISHED FILTERING AT:\t" + now + "###\n\n###LOREAN KEPT\t" + str(
                filter_count) + "\tREADS AFTER LENGTH FILTERING###\n"))

            long_sam = mapping.gmap('sam', genome_gmap, long_fasta, args.threads, 'samse',
                                    args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                    args.verbose, Fflag=False)
            long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, args.threads, gmap_wd, args.verbose)

    # FILE
    fileName = consensus_wd + 'mergedGmapEvm.beforeAssembly.gff3'
    # HERE WE CHECK IF WE HAVE THE PASA UPDATED FILE OR THE EVM
    # ORIGINAL FILE
    if os.path.isfile(updatedGff3):
        # HERE WE MERGE THE TWO FILES
        mergedmapGFF3 = logistic.catTwoBeds(long_sorted_bam, updatedGff3, fileName, args.update)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("\n\t###GFFREAD\t" + now + "\t###\n"))

    # HERE WE TRANSFORM THE COODINATES INTO SEQUENCES USING THE
    # REFERENCE
    gffreadFastaFile = consensus.gffread(mergedmapGFF3, ref, consensus_wd, args.verbose)
    # HERE WE STORE THE SEQUENCE IN A DICTIONARY
    fake = []
    gffreadDict = consensus.fasta2Dict(gffreadFastaFile)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("\n\t#CLUSTERING\t" + now + "\t###\n"))

    # HERE WE CLUSTER THE SEQUENCES BASED ON THE GENOME
    # POSITION
    cluster_list = consensus.cluster_pipeline(mergedmapGFF3, args.assembly_overlapLength, args.stranded)
    now = datetime.datetime.now().strftime(fmtdate)

    sys.stdout.write(("\n\t#CONSENSUS FOR EACH CLUSTER\t" + now + "\t###\n"))

    # HERE WE MAKE CONSENSUS FOR EACH CLUSTER
    tmp_wd = consensus_wd + 'tmp/'
    logistic.check_create_dir(tmp_wd)
    tmp_assembly_file = tmp_wd + 'assembly.fasta'
    if os.path.isfile(tmp_assembly_file):
        sys.stdout.write('No assembly')
    else:
        consensus.generate_fasta(cluster_list, gffreadDict, args.cluster_min_evidence,
                                 args.cluster_max_evidence, args.assembly_overlapLength, tmp_wd)
        consensus.assembly(args.assembly_overlapLength, args.assembly_percentIdentity, args.threads, tmp_wd,
                           args.verbose)
        utrs.lengthSupport(tmp_wd, args.threads)

        # WITH THE ELSE, WE ALLOW THE USER TO DECIDE TO CHANGE THE ASSEMBLY
        # PARAMETERS AND COLLECT DIFFERENT ASSEMBLED SEQUENCES WITHOT RUNNING
        # THE FULL PIPELINE
        # HERE WE COLLECT THE ASSEMBLED SEQUENCES. WE COLLWCT ONLY SEQUENCE
        # THAT PASS THE FILTER
    tmp_consensus = os.path.join(consensus_wd, 'tmp/')
    collect.parse_only(args.assembly_readThreshold, tmp_consensus)
    tmp_assembly = collect.catAssembled(tmp_consensus)
    # HERE WE COLLECT THE NEW ASSEMBLED SEQUENCES AND WE COLLECT THE OLD
    # EVM DATA
    mergedFastaFilename = consensus_wd + 'assembly.wEVM.fasta'
    collect.addEVM(gffreadFastaFile, tmp_assembly, mergedFastaFilename)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("\n###MAPPING CONSENSUS ASSEMBLIES\t" + now + "\t###\n"))

    # HERE WE MAP ALL THE FASTA FILES TO THE GENOME USING GMAP
    consensusMappedGFF3 = mapping.gmap('cons', genome_gmap, mergedFastaFilename, args.threads, 'gff3_gene',
                                       args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                       args.verbose,
                                       Fflag=True)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("\n###GETTING THE STRAND RIGHT\t" + now + "\t###\n"))
    # IN THIS STEP WE CORRECT FOR STRAND. GMAP CAN NOT DECIDE THE STRAND
    # FOR SINGLE EXONS GENE MODELS. WE USE THE ORIENTATION FROM EVM IF GMAP
    # INVERT THE ORIGINAL STRAND

    strandMappedGFF3 = grs.strand(evm, consensusMappedGFF3, ref, args.threads, gmap_wd, args.verbose)
    gffPasa = grs.appendID(strandMappedGFF3)
    noOverl = grs.removeOverlap(gffPasa, args.verbose)
    noDisc = grs.removeDiscrepancy(noOverl, evm_gff3, args.verbose)
    uniqGene = grs.newNames(noDisc)

    finalupdate3 = grs.genename(uniqGene, args.prefix_gene, args.verbose)
    print(("\n###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))
    finalupdate4 = grs.exonerate(ref, finalupdate3, args.threads, exonerate_wd, args.verbose)
    finalupdate5 = grs.genename(finalupdate4, args.prefix_gene, args.verbose)

    # HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN
    # PASA TO CORRECT SMALL ERRORS

    sys.stdout.write(("\n###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))

    fastaAll = logistic.catTwoFasta(trinity_out, mergedFastaFilename, long_fasta, pasa_dir)
    round_n += 1

    finalupdate = pasa.update_database(args.threads, str(round_n), pasa_dir, args.pasa_db, align_pasa_conf, ref,
                                       fastaAll,
                                       finalupdate5, args.verbose)
    round_n += 1
    finalupdate2 = pasa.update_database(args.threads, str(round_n), pasa_dir, args.pasa_db, align_pasa_conf, ref,
                                        fastaAll,
                                        finalupdate, args.verbose)
    finalUpdate = grs.genename(finalupdate2, args.prefix_gene, args.verbose)

    FinalFiles.append(finalUpdate)

    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(('\n###CREATING OUTPUT DIRECTORY\t' + now + '\t###\n'))

    final_output_dir = os.path.join(root, 'output_annotation/')
    logistic.check_create_dir(final_output_dir)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("\n##PLACING OUTPUT FILES IN OUTPUT DIRECTORY\t" + now + "\t###\n"))

    for filename in FinalFiles:
        if filename != '':
            logistic.copy_file(filename, final_output_dir)
            cmdstring = "chmod -R 775 %s" % wd
            os.system(cmdstring)
    if not args.keep_tmp and args.working_dir == '':
        temp_dir.cleanup()
