#! /usr/bin/env python3

###############
###IMPORTS###
###############

# LIBRARIES
#TODO
#check for strand on single exon gene based on reads mapping


import datetime
import os
import shutil
import subprocess
import sys
import time
from queue import Queue
from threading import Thread

import collect_only as collect
import arguments as arguments
import consensus_iAssembler as consensus
# OTHER SCRIPTS
import dirs_and_files as logistic
import evm_pipeline
import get_right_strand as grs
import handlers as handler
import manipulateSeq as mseq
import mapping
import multithread_large_fasta as multiple
import prepare_evm_inputs as inputEvm
import reduceUTRs as utrs
import transcript_assembly as transcripts


###############
###MAIN###
###############

def main():
    if os.path.isfile("/data/gm_key"):
        '''Core of the program'''
        # Parse the arguments
        args = arguments.setting()
        fmtdate = '%H:%M:%S %d-%m'
        now = datetime.datetime.now().strftime(fmtdate)
        # Useful variables for later
        wd_base = os.path.abspath(args.working_dir) + '/'
        wd = wd_base + 'run/'
        ref = os.path.abspath(args.ref)

        gmap_name = args.ref + '_GMAPindex'
        pasa_name = 'assembler-' + args.pasa_db

        if args.short_reads == '' and args.long_reads == '':
            weights_dic = {'Augustus': args.augustus, 'GeneMark.hmm': args.genemark, 'AAT': args.AAT}

        elif args.short_reads != '' or args.long_reads != '':
            weights_dic = {'Augustus': args.augustus, pasa_name : args.pasa, 'GeneMark.hmm': args.genemark,
                           'AAT': args.AAT, gmap_name: args.trinity}
        FinalFiles = []  # STORE THE IMPORTANT OUTPUT FILES

        logistic.check_create_dir(wd)
        logistic.check_file(ref)
        gmap_wd = wd + '/gmap_output/'
        exonerate_wd = wd + '/exonerate/'
        pasa_dir = wd + 'PASA/'
        star_out = wd + '/STAR/'
        trin_dir = wd + 'Trinity/'
        evm_dir = wd + 'evm_inputs/'
        braker_out = wd + 'braker/' + args.species + '/'


        logistic.check_create_dir(evm_dir)
        logistic.check_create_dir(trin_dir)
        logistic.check_create_dir(star_out)
        logistic.check_create_dir(pasa_dir)
        logistic.check_create_dir(gmap_wd)
        logistic.check_create_dir(exonerate_wd)

        check_species = ['augustus', '--species=help']
        process = subprocess.Popen(check_species, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outAugustus, errAugustus = process.communicate()
        protein_loc = os.path.abspath(args.protein_evidence)

        if args.repeat_masked:
            genome_gmap = mseq.maskedgenome(gmap_wd, ref, args.repeat_masked)
        else:
            genome_gmap = ref

        # COLLECT ONLY ONLY RUNS PART OF THE CONSENSUS PIPELINE
        list_fasta_names = multiple.single_fasta(ref, wd)
        if args.short_reads or args.long_reads:
            now = datetime.datetime.now().strftime(fmtdate)
            print(('\n###STAR MAPPING  STARTED AT:\t'  + now + '\t###\n'))
            # SHORT READS
            if 'fastq' in args.short_reads or 'fq' in args.short_reads:
                if ',' in args.short_reads:
                    pairedEndFiles = args.short_reads.split(',')
                    short_1 = os.path.abspath(pairedEndFiles[0])
                    short_2 = os.path.abspath(pairedEndFiles[1])
                    short_reads_file = [short_1, short_2]
                else:
                    short_reads_file = os.path.abspath(args.short_reads)
                # Map with STAR
                short_bam = mapping.star(ref, short_reads_file, args.threads, args.max_intron_length, star_out, args.verbose)
                short_sorted_bam = mapping.samtools_sort(short_bam, args.threads, wd, args.verbose)
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
                short_bam = star_out + "/" + bam_file[-1]
                short_sorted_bam = mapping.samtools_sort(short_bam, args.threads, wd, args.verbose)
                # print short_sorted_bam
            else:
                short_sorted_bam = False
                print('No short reads file')
            # LONG READS
            if 'fastq' in args.long_reads or 'fq' in args.long_reads or 'fasta' in args.long_reads or 'fa' in args.long_reads  :
                # with this operation, reads are filtered for their length.
                # Nanopore reads can be chimaras or sequencing artefacts.
                # filtering on length reduces the amount of sequencing
                # artefacts
                now = datetime.datetime.now().strftime(fmtdate)
                print(("\n###FILTERING OUT LONG READS STARTED AT:\t"  +  now   + "\t###\n"))
                long_fasta, filter_count = mseq.filterLongReads(args.long_reads, args.assembly_overlapLength,
                                                                args.max_long_read, gmap_wd, args.adapter, args.threads,
                                                                a = True)
                if filter_count != 0:
                    now = datetime.datetime.now().strftime(fmtdate)
                    print(("###FINISHED FILTERING AT:\t" + now + "###\n\n###LOREAN KEPT\t" + str(filter_count) + "\tREADS AFTER LENGTH FILTERING###\n"))
                if not short_sorted_bam:
                    # If short reads have been mapped dont do it
                    now = datetime.datetime.now().strftime(fmtdate)
                    print(('\n###GMAP\t'  + now  + 't###\n'))
                    long_sam = mapping.gmap('sam', genome_gmap, long_fasta, args.threads, 'samse', args.min_intron_length,
                                            args.max_intron_length, args.end_exon, gmap_wd, args.verbose, Fflag=False)
                    # Convert to sorted BAM
                    long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, args.threads, wd, args.verbose)

                    # Keep the output
                    FinalFiles.append(long_sorted_bam)
                else:
                    long_sorted_bam = False

            else:
                now = datetime.datetime.now().strftime(fmtdate)
                print(('\n###NO LONG READS FILE\t'  + now  + '\t###\n'))
                long_sorted_bam = False
            if short_sorted_bam:  # If there are short reads, these will serve to the transcript assembly pipeline
                default_bam = short_sorted_bam
            else:
                default_bam = long_sorted_bam
            # TRANSCRIPT ASSEMBLY
            # TRINITY
            now = datetime.datetime.now().strftime(fmtdate)
            print(('\n###TRINITY STARTS AT:\t'  + now  + '\t###\n'))
            if int(args.threads) > 1:
                trinity_cpu = int(int(args.threads)/int(2))
            else:
                trinity_cpu = int(args.threads)
            trinity_out = transcripts.trinity(default_bam, trin_dir, args.max_intron_length, trinity_cpu, args.verbose)
            trinityGFF3 = mapping.gmap('trin', genome_gmap, trinity_out, args.threads, 'gff3_gene', args.min_intron_length,
                                       args.max_intron_length, args.end_exon, gmap_wd, args.verbose, Fflag=True)
            trinity_path = trinityGFF3

            # PASA Pipeline
            now = datetime.datetime.now().strftime(fmtdate)
            print(('\n###PASA STARTS AT:\t'  + now  + '\t###\n'))
            # Create PASA folder and configuration file
            align_pasa_conf = transcripts.pasa_configuration(pasa_dir, args.pasa_db)
            # Launch PASA
            pasa_gff3 = transcripts.pasa_call(pasa_dir, align_pasa_conf, args.pasa_db, ref, trinity_out, args.max_intron_length,
                                              args.threads, args.verbose)

            # HERE WE PARALLELIZE PROCESSES WHEN MULTIPLE THREADS ARE USED
            if args.species in (errAugustus.decode("utf-8")):
                now = datetime.datetime.now().strftime(fmtdate)
                print(('\n###AUGUSTUS, GENEMARK-ES AND AAT STARTED AT:'  + now  + '\t###\n'))
                queue = Queue()
                for i in range(3):
                    queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                    for i in range(3):
                        t = Thread(target=handler.AugustGmesAAT,args=(queue, ref, args.species, protein_loc,
                                                                      args.threads, args.fungus, list_fasta_names, wd, args.verbose))
                        t.daemon = True
                        t.start()
                queue.join()
                augustus_file = wd + 'augustus/augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file = wd + 'gmes/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                mergedProtGFF3 = wd + 'AAT/protein_evidence.gff3'

            elif args.short_reads:  # USING PROTEINS AND SHORT READS
                    now = datetime.datetime.now().strftime(fmtdate)
                    print(('\n###BRAKER1 (USING SHORT READS) AND AAT STARTED AT:\t'  + now  + '\t###\n'))
                    queue = Queue()
                    for i in range(2):
                        queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                        for i in range(2):
                            t = Thread(target=handler.BrakerAAT, args=(queue, ref, default_bam, args.species, protein_loc,
                                                                       args.threads, args.fungus, list_fasta_names, wd, args.verbose))
                            t.daemon = True
                            t.start()
                    queue.join()
                    augustus_file = braker_out + 'augustus.gff'
                    augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                    genemark_file = braker_out + 'GeneMark-ET/genemark.gtf'
                    genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                    mergedProtGFF3 = wd + 'AAT/protein_evidence.gff3'

            else:  # USING PROTEINS AND LONG READS
                queue = Queue()
                now = datetime.datetime.now().strftime(fmtdate)
                print(('\n###BRAKER1 (USING LONG READS) AND AAT STARTED AT: \t'  + now  + '\t###\n'))
                for i in range(2):
                    queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                for i in range(2):
                    t = Thread(target=handler.BrakerAAT, args=(queue, ref, long_sorted_bam, args.species, protein_loc,
                                                               args.threads, args.fungus, list_fasta_names, wd, args.verbose))
                    t.daemon = True
                    t.start()
                queue.join()
                augustus_file = braker_out + 'augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file = braker_out + 'GeneMark-ET/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                mergedProtGFF3 = wd + 'AAT/protein_evidence.gff3'
        elif args.species in (errAugustus.decode("utf-8")):
            now = datetime.datetime.now().strftime(fmtdate)
            print(('\n###AUGUSTUS, GENEMARK-ES AND AAT STARTED AT:'  + now  + '\t###\n'))
            queue = Queue()
            for i in range(3):
                queue.put(i)  # QUEUE WITH A ZERO AND A ONE
                for i in range(3):
                    t = Thread(target=handler.AugustGmesAAT,args=(queue, ref, args.species, protein_loc,
                                                                  args.threads, args.fungus, list_fasta_names, wd, args.verbose))
                    t.daemon = True
                    t.start()
            queue.join()
            augustus_file = wd + 'augustus/augustus.gff'
            augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
            genemark_file = wd + 'gmes/genemark.gtf'
            genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
            mergedProtGFF3 = wd + 'AAT/protein_evidence.gff3'
        else:
            now = datetime.datetime.now().strftime(fmtdate)
            sys.exit("#####UNRECOGNIZED SPECIES FOR AUGUSTUS AND NO READS\t"  + now + "\t#####\n")
        # Prepare EVM input files
        now = datetime.datetime.now().strftime(fmtdate)
        print(('\n###EVM STARTED AT:\t'  + now  + '\t###\n'))
        # HERE WE CONVERT FILES FOR EVM AND PLACE THEM IN INPUT FOLDER

        if not args.short_reads and not args.long_reads:
            evm_inputs = {'augustus': augustus_gff3, 'genemark': genemark_gff3, 'AAT': mergedProtGFF3}
        elif args.short_reads or args.long_reads:
            evm_inputs = {'pasa': pasa_gff3,'augustus': augustus_gff3, 'genemark': genemark_gff3, 'AAT': mergedProtGFF3,
                          'gmap': trinity_path}

        # HERE WE RUN EVM; WE PREPARE FILES THAT ARE REQUIRED BY EVM LIKE
        # WEIGTH TABLE

        # print '> EVM input directory created in ' + evm_dir

        list_soft, pred_file, transcript_file, protein_file = inputEvm.group_EVM_inputs(evm_dir, evm_inputs)
        weight_file = inputEvm.evm_weight(evm_dir, weights_dic, list_soft, pasa_name, gmap_name)
        # EVM PIPELINE


        if args.short_reads or args.long_reads:  # WE HAVE SHORT READS AND PROTEINS
            evm_gff3 = evm_pipeline.evm_pipeline(wd, args.threads, genome_gmap, weight_file, pred_file, transcript_file,
                                                 protein_file, args.segmentSize, args.overlapSize)
        elif not args.short_reads and not args.long_reads:  # WE HAVE PROTEINS BUT NOT SHORT READS
            transcript_file = ''
            evm_gff3 = evm_pipeline.evm_pipeline(wd, args.threads, genome_gmap, weight_file, pred_file, transcript_file,
                                                 protein_file, args.segmentSize, args.overlapSize)
        # KEEP THIS OUTPUT
        FinalFiles.append(evm_gff3)
        if not args.short_reads and not args.long_reads:
            now = datetime.datetime.now().strftime(fmtdate)
            sys.exit("##### EVM FINISHED AT:\t"  + now  + "\t#####\n")
            # RE-RUN PASA PIPELINE
        # HERE WE CAN EXCLUDE TO RUN AGAIN PASA TO UPDATE THE DATABASE
        # AFTER EVM; #We only want to update if it ran with short reads
        round_n = 0
        if args.short_reads and not args.long_reads:
            now = datetime.datetime.now().strftime(fmtdate)
            print(('\n###UPDATE WITH PASA DATABASE STARTED AT:\t ' +   now  + '\t###\n'))
            round_n += 1
            finalOutput = evm_pipeline.update_database(args.threads, str(round_n), pasa_dir, args.pasa_db,
                                                       align_pasa_conf, ref, trinity_out, evm_gff3, "a")
            finalUpdate = grs.genename(finalOutput, args.prefix_gene)
            updatedGff3 = grs.newNames(finalUpdate)
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
            cmdstring = "chmod -R 775 %s" % wd
            os.system(cmdstring)
            now = datetime.datetime.now().strftime(fmtdate)
            sys.exit("#####ANNOTATION FINISHED WITHOUT USING LONG READS\t"  + now  + "\t#####\n")

        # HERE WE START WITH LONG READS
        else:
            now = datetime.datetime.now().strftime(fmtdate)
            print(('\n###RUNNING iASSEMBLER\t'  + now  + '\t###\n'))

            if args.long_reads:
                # Means there are long reads to map and user wants to run
                # this pipeline

                consensus_wd = (wd + 'consensus/')
                logistic.check_create_dir(consensus_wd)
                # HERE WE MAP THE READS ON THE GENOME USING GMAP

                if not long_sorted_bam:
                    long_sam = mapping.gmap( 'sam', genome_gmap, long_fasta, args.threads, 'samse', args.min_intron_length,
                                             args.max_intron_length, args.end_exon, gmap_wd, args.verbose, Fflag=False)
                    long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, args.threads, wd, args.verbose)

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
                now = datetime.datetime.now().strftime(fmtdate)
                print(("\n\t###GFFREAD\t"  + now  + "\t###\n"))

                # HERE WE TRANSFORM THE COODINATES INTO SEQUENCES USING THE
                # REFERENCE
                gffreadFastaFile = consensus.gffread(
                    mergedmapGFF3, ref, consensus_wd, args.verbose)
                # HERE WE STORE THE SEQUENCE IN A DICTIONARY
                fake = []
                long_fasta, filter_count = mseq.filterLongReads(
                    gffreadFastaFile, args.assembly_overlapLength, args.max_long_read, consensus_wd,  fake, args.threads, a = False)

                gffreadDict = consensus.fasta2Dict(gffreadFastaFile)
                now = datetime.datetime.now().strftime(fmtdate)
                print(("\n\t#CLUSTERING\t"  + now  + "\t###\n"))

                # HERE WE CLUSTER THE SEQUENCES BASED ON THE GENOME
                # POSITION
                cluster_list = consensus.cluster_pipeline(
                    mergedmapGFF3, args.assembly_overlapLength, args.stranded)
                now = datetime.datetime.now().strftime(fmtdate)

                print(("\n\t#CONSENSUS FOR EACH CLUSTER\t"  + now  + "\t###\n"))

                # HERE WE MAKE CONSENSUS FOR EACH CLUSTER
                tmp_wd = consensus_wd + 'tmp/'
                logistic.check_create_dir(tmp_wd)
                tmp_assembly_file = tmp_wd + 'assembly.fasta'
                if os.path.isfile(tmp_assembly_file):
                    print('No assembly')
                else:
                    consensus.generate_fasta(
                        cluster_list,
                        gffreadDict,
                        args.cluster_min_evidence,
                        args.cluster_max_evidence,
                        args.assembly_overlapLength,
                        tmp_wd)
                    consensus.assembly(args.assembly_overlapLength, args.assembly_percentIdentity, args.threads, tmp_wd, args.verbose)
                    utrs.lengthSupport(tmp_wd, args.threads)

        # WITH THE ELSE, WE ALLOW THE USER TO DECIDE TO CHANGE THE ASSEMBLY
        # PARAMETERS AND COLLECT DIFFERENT ASSEMBLED SEQUENCES WITHOT RUNNING
        # THE FULL PIPELINE
        # HERE WE COLLECT THE ASSEMBLED SEQUENCES. WE COLLWCT ONLY SEQUENCE
        # THAT PASS THE FILTER
        evm_nosupport = collect.parse_only(args.assembly_readThreshold, wd)
        tmp_assembly = collect.catAssembled(wd)
        # HERE WE COLLECT THE NEW ASSEMBLED SEQUENCES AND WE COLLECT THE OLD
        # EVM DATA
        mergedFastaFilename = consensus_wd + 'assembly.wEVM.fasta'
        collect.addEVM(
            gffreadFastaFile,
            tmp_assembly,
            mergedFastaFilename)
        now = datetime.datetime.now().strftime(fmtdate)
        print(("\n###MAPPING CONSENSUS ASSEMBLIES\t"  + now + "\t###\n"))

        # HERE WE MAP ALL THE FASTA FILES TO THE GENOME USING GMAP
        consensusMappedGFF3 = mapping.gmap(
            'cons',
            genome_gmap,
            mergedFastaFilename,
            args.threads,
            'gff3_gene',
            args.min_intron_length,
            args.max_intron_length,
            args.end_exon,
            gmap_wd, args.verbose,
            Fflag=True)
        now = datetime.datetime.now().strftime(fmtdate)
        print(("\n###GETTING THE STRAND RIGHT\t"  + now  + "\t###\n"))
        # IN THIS STEP WE CORRECT FOR STRAND. GMAP CAN NOT DECIDE THE STRAND
        # FOR SINGLE EXONS GENE MODELS. WE USE THE ORIENTATION FROM EVM IF GMAP
        # INVERT THE ORIGINAL STRAND
        #finalOutput = grs.strand(evm_gff3, consensusMappedGFF3, gmap_wd)
        strandMappedGFF3 = grs.strand(evm_gff3, consensusMappedGFF3, ref, args.threads, gmap_wd)
        gffPasa = grs.appendID(strandMappedGFF3)
        noOverl = grs.removeOverlap(gffPasa)
        noDisc = grs.removeDiscrepancy(noOverl, evm_gff3)
        uniqGene = grs.newNames(noDisc)
        # HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN
        # PASA TO CORRECT SMALL ERRORS

        
        finalupdate3 = grs.genename(uniqGene, args.prefix_gene)
        print(("\n###FIXING GENES NON STARTING WITH MET\t"  + now  + "\t###\n"))
        finalupdate4 = grs.exonerate(ref, finalupdate3, args.threads, exonerate_wd)
        finalupdate5 = grs.genename(finalupdate4, args.prefix_gene)
        
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
            finalupdate5,
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
        finalUpdate = grs.genename(finalupdate2, args.prefix_gene)
        FinalFiles.append(finalUpdate)

        now = datetime.datetime.now().strftime(fmtdate)
        print(('\n###CREATING OUTPUT DIRECTORY\t'  + now  + '\t###\n'))

        final_output_dir = wd_base + 'output/'
        logistic.check_create_dir(final_output_dir)
        now = datetime.datetime.now().strftime(fmtdate)
        print(("\n##PLACING OUTPUT FILES IN OUTPUT DIRECTORY\t"  + now  + "\t###\n"))

        for filename in FinalFiles:
            if filename != '':
                logistic.copy_file(filename, final_output_dir)
                cmdstring = "chmod -R 775 %s" % wd
                os.system(cmdstring)
        if args.keep_tmp:
            dirs_list = ['/PASA/', 'augustus/', 'gmes/', 'AAT/', 'split/']
            for dirs in dirs_list:
                dest = wd + dirs
                shutil.rmtree(dest, ignore_errors=True)

    else:
        print('Key for GeneMark-ES not found.  Please, place the GeneMark-ES key in the folder where you have your data.')
        sys.exit("#####LOREAN STOPS HERE.#####\n")

if __name__ == '__main__':
    realstart = time.perf_counter()
    main()
    realend = time.perf_counter()
    realt = round(realend - realstart)
    m, s = divmod(realt, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    print("###LOREAN FINISHED WITHOUT ERRORS IN:" + " " + str(d) + " days " + str(h)  + " hours " + str(m) + " min " + str(s) + " sec\t###\n")
