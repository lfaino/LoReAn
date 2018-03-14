#! /usr/bin/env python3

###############
###IMPORTS###
###############

# LIBRARIES
import datetime
import multiprocessing
import os
import sys
import tempfile
import time
from queue import Queue
from threading import Thread

# OTHER SCRIPTS
import arguments as arguments
import collectOnly as collect
import consensusIAssembler as consensus
import dirsAndFiles as logistic
import evmPipeline
import getRightStrand as grs
import handlers as handler
import interproscan as iprscan
import manipulateSeq as mseq
import mapping
import multithreadLargeFasta as multiple
import pasa as pasa
import prepareEvmInputs as inputEvm
import reduceUTRs as utrs
import transcriptAssembly as transcripts
import update as update


###############
###MAIN###
###############

def main():

    fmtdate = '%H:%M:%S %d-%m'
    now = datetime.datetime.now().strftime(fmtdate)
    home = os.path.expanduser("~")
    args = arguments.setting()
    max_threads = multiprocessing.cpu_count()
    gmap_name = args.reference + '_GMAPindex'
    pasa_name = 'assembler-' + args.pasa_db
    protein_loc = os.path.abspath(args.proteins)
    iprscan_log = iprscan.check_iprscan()
    # Useful variables for later
    root = os.getcwd()
    
    output_dir = os.path.join(root, "LoReAn_" + args.working_dir)
    logistic.check_create_dir(output_dir)

    wd = os.path.join(output_dir, "run/")
    if args.keep_tmp:
        logistic.check_create_dir(wd)
    elif not os.path.exists(wd) and args.verbose:
        logistic.check_create_dir(wd)
    else:
        temp_dir = tempfile.TemporaryDirectory(prefix='run_', dir=output_dir, suffix="/", )
        wd = temp_dir.name

    ref_orig = os.path.abspath(args.reference)
    ref_link = os.path.join(wd, args.reference)
    if not os.path.exists(ref_link):
        os.link(ref_orig, ref_link)

    if args.upgrade:
        update.upgrade()

    elif os.path.isfile(home + "/.gm_key") and args.proteins != "":
        fasta = (".fasta", ".fa", ".fas", ".fsta")
        fastq = (".fastq", ".fq")
        '''Core of the program'''
        # Parse the arguments
        if int(args.threads) > max_threads:
            threads_use = str(max_threads)
            sys.stdout.write(('\n### MAX NUMBER OF USED THREADS IS ' + str(max_threads) + ' AND NOT ' + args.threads + ' AS SET ###\n'))
        else:
            threads_use = args.threads
        if args.external:
            external_file = args.external
        else:
            external_file = ''

        if args.short_reads == '' and args.long_reads == '':
            if external_file.endswith("gff3") or external_file.endswith(fasta):
                weights_dic = {'Augustus': args.augustus_weigth, 'GeneMark.hmm': args.genemark_weigth, 'AAT': args.AAT_weigth,
                               'external' : args.external_weigth}
            else:
                weights_dic = {'Augustus': args.augustus_weigth, 'GeneMark.hmm': args.genemark_weigth, 'AAT': args.AAT_weigth}
        elif args.short_reads != '' or args.long_reads != '':
            if external_file.endswith("gff3") or external_file.endswith(fasta):
                weights_dic = {'Augustus': args.augustus_weigth, pasa_name: args.pasa_weigth, 'GeneMark.hmm': args.genemark_weigth,
                               'AAT': args.AAT_weigth, gmap_name: args.trinity_weigth, 'external' : args.external_weigth}
            else:
                weights_dic = {'Augustus': args.augustus_weigth, pasa_name: args.pasa_weigth, 'GeneMark.hmm': args.genemark_weigth,
                           'AAT': args.AAT_weigth, gmap_name: args.trinity_weigth}

        final_files = []  # STORE THE IMPORTANT OUTPUT FILES

        logistic.check_create_dir(wd)
        logistic.check_file(ref_link)
        gmap_wd = wd + '/gmap_output/'
        exonerate_wd = wd + '/exonerate/'
        pasa_dir = wd + 'PASA/'
        star_out = wd + '/STAR/'
        trin_dir = wd + '/Trinity/'
        evm_inputs_dir = wd + '/evm_inputs/'
        braker_folder = wd + '/braker/'
        evm_output_dir = wd + '/evm_output/'
        interproscan_out_dir = wd + 'interproscan'
        logistic.check_create_dir(evm_inputs_dir)
        logistic.check_create_dir(evm_output_dir)
        logistic.check_create_dir(trin_dir)
        logistic.check_create_dir(star_out)
        logistic.check_create_dir(pasa_dir)
        logistic.check_create_dir(gmap_wd)
        logistic.check_create_dir(interproscan_out_dir)

        if args.long_reads:
            logistic.check_create_dir(exonerate_wd)
            consensus_wd = (wd + '/consensus/')
            logistic.check_create_dir(consensus_wd)

        logistic.check_gmap(threads_use, 'samse', args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                            args.verbose)

        augustus_species, err_augustus = logistic.augustus_species_func(home)

        if args.repeat_masked:
            masked_ref = mseq.maskedgenome(gmap_wd, ref_link, args.repeat_masked, args.repeat)
        else:
            masked_ref = ref_link

        # COLLECT ONLY ONLY RUNS PART OF THE CONSENSUS PIPELINE
        list_fasta_names, dict_ref_name = multiple.single_fasta(masked_ref, wd)
        if args.short_reads or args.long_reads:
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('\n###STAR MAPPING  STARTED AT:\t' + now + '\t###\n'))
            # SHORT READS
            if args.short_reads.endswith(fastq):
                if ',' in args.short_reads:
                    paired_end_files = args.short_reads.split(',')
                    short_1 = os.path.abspath(paired_end_files[0])
                    short_2 = os.path.abspath(paired_end_files[1])
                    short_reads_file = [short_1, short_2]
                else:
                    short_reads_file = os.path.abspath(args.short_reads)
                # Map with STAR
                short_bam = mapping.star(masked_ref, short_reads_file, threads_use, args.max_intron_length, star_out,
                                         args.verbose)
                short_sorted_bam = mapping.samtools_sort(short_bam, threads_use, wd, args.verbose)
                # Keep the output
                final_files.append(short_sorted_bam)
            # BAM SORTED FILES GET IN HERE
            elif args.short_reads.endswith("bam"):
                logistic.check_create_dir(star_out)
                short_sorted_bam = os.path.abspath(args.short_reads)
                bam_file = short_sorted_bam.split("/")
                short_bam = star_out + "/" + bam_file[-1]
                if not os.path.exists(ref_link):
                    os.link(short_sorted_bam, short_bam)

            else:
                short_sorted_bam = False
                sys.stdout.write("\n\033[31m ### NO SHORT READS ### \033[0m\n")

            # LONG READS
            if 'fastq' in args.long_reads or 'fq' in args.long_reads or 'fasta' in args.long_reads or 'fa' in args.long_reads:
                # with this operation, reads are filtered for their length.
                # Nanopore reads can be chimaras or sequencing artefacts.
                # filtering on length reduces the amount of sequencing
                # artefacts
                now = datetime.datetime.now().strftime(fmtdate)
                sys.stdout.write(("\n###FILTERING OUT LONG READS STARTED AT:\t" + now + "\t###\n"))
                long_fasta = mseq.filterLongReads(args.long_reads, args.assembly_overlap_length, args.max_long_read, gmap_wd,
                                                  args.adapter, threads_use,  a=True)
                if not short_sorted_bam:
                    # If short reads have been mapped dont do it
                    now = datetime.datetime.now().strftime(fmtdate)
                    sys.stdout.write(('\n###GMAP\t' + now + 't###\n'))
                    long_sam = mapping.gmap('sam', masked_ref, long_fasta, threads_use, 'samse',
                                            args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                            args.verbose, Fflag=False)
                    # Convert to sorted BAM
                    long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, threads_use, wd, args.verbose)

                    # Keep the output
                    final_files.append(long_sorted_bam)
                else:
                    long_sorted_bam = False

            else:
                now = datetime.datetime.now().strftime(fmtdate)
                sys.stdout.write(('\n###NO LONG READS FILE\t' + now + '\t###\n'))
                long_sorted_bam = False
            if short_sorted_bam:  # If there are short reads, these will serve to the transcript assembly pipeline
                default_bam = short_sorted_bam
            else:
                default_bam = long_sorted_bam
            # TRANSCRIPT ASSEMBLY
            # TRINITY
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('\n###TRINITY STARTS AT:\t' + now + '\t###\n'))
            if int(threads_use) > 1:
                trinity_cpu = int(int(threads_use) / int(2))
            else:
                trinity_cpu = int(threads_use)
            trinity_out = transcripts.trinity(default_bam, trin_dir, args.max_intron_length, trinity_cpu, args.verbose)
            trinity_gff3 = mapping.gmap('trin', masked_ref, trinity_out, threads_use, 'gff3_gene',
                                        args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                        args.verbose, Fflag=True)
            trinity_path = trinity_gff3
            # PASA Pipeline
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('\n###PASA STARTS AT:\t' + now + '\t###\n'))
            # Create PASA folder and configuration file
            #align_pasa_conf = pasa.pasa_configuration(pasa_dir, args.pasa_db, args.verbose)
            # Launch PASA
            pasa_gff3 = pasa.pasa_call(pasa_dir, args.pasa_db, masked_ref, trinity_out, args.max_intron_length, threads_use, args.verbose)

            # HERE WE PARALLELIZE PROCESSES WHEN MULTIPLE THREADS ARE USED
            if args.species in (err_augustus.decode("utf-8")) or args.species in augustus_species:
                now = datetime.datetime.now().strftime(fmtdate)
                sys.stdout.write(('\n###AUGUSTUS, GENEMARK-ES AND AAT STARTED AT:' + now + '\t###\n'))
                queue = Queue()
                for software in range(3):
                    queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                    for software in range(3):
                        t = Thread(target=handler.august_gmes_aat, args=(queue, masked_ref, args.species, protein_loc,
                                                                       threads_use, args.fungus, list_fasta_names, wd,
                                                                       args.verbose))
                        t.daemon = True
                        t.start()
                queue.join()
                augustus_file = wd + 'augustus/augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file = wd + 'gmes/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                merged_prot_gff3 = wd + 'AAT/protein_evidence.gff3'

            elif args.short_reads:  # USING PROTEINS AND SHORT READS
                logistic.check_create_dir(braker_folder)
                now = datetime.datetime.now().strftime(fmtdate)
                sys.stdout.write(('\n###BRAKER1 (USING SHORT READS) AND AAT STARTED AT:\t' + now + '\t###\n'))
                queue = Queue()
                for software in range(2):
                    queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                    for software in range(2):
                        t = Thread(target=handler.braker_aat, args=(queue, masked_ref, default_bam, args.species, protein_loc,
                                                                   threads_use, args.fungus, list_fasta_names, wd,
                                                                    braker_folder, args.verbose))
                        t.daemon = True
                        t.start()
                queue.join()
                augustus_file = braker_folder + 'augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file = braker_folder + 'GeneMark-ET/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                merged_prot_gff3 = wd + 'AAT/protein_evidence.gff3'

            else:  # USING PROTEINS AND LONG READS
                queue = Queue()
                now = datetime.datetime.now().strftime(fmtdate)
                sys.stdout.write(('\n###BRAKER1 (USING LONG READS) AND AAT STARTED AT: \t' + now + '\t###\n'))
                logistic.check_create_dir(braker_folder)
                for software in range(2):
                    queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                    for software in range(2):
                        t = Thread(target=handler.braker_aat,
                                   args=(queue, masked_ref, long_sorted_bam, args.species, protein_loc,
                                         threads_use, args.fungus, list_fasta_names, wd, braker_folder, args.verbose))
                        t.daemon = True
                        t.start()
                queue.join()
                augustus_file = braker_folder + 'augustus.gff'
                augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                genemark_file = braker_folder + 'GeneMark-ET/genemark.gtf'
                genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                merged_prot_gff3 = wd + 'AAT/protein_evidence.gff3'
        elif args.species in (err_augustus.decode("utf-8")) or args.species in augustus_species:
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('\n###AUGUSTUS, GENEMARK-ES AND AAT STARTED AT:' + now + '\t###\n'))
            queue = Queue()
            for software in range(3):
                queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                for software in range(3):
                    t = Thread(target=handler.august_gmes_aat, args=(queue, masked_ref, args.species, protein_loc,
                                                                   threads_use, args.fungus, list_fasta_names, wd,
                                                                   args.verbose))
                    t.daemon = True
                    t.start()
            queue.join()
            augustus_file = wd + 'augustus/augustus.gff'
            augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
            genemark_file = wd + 'gmes/genemark.gtf'
            genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
            merged_prot_gff3 = wd + 'AAT/protein_evidence.gff3'
        else:
            now = datetime.datetime.now().strftime(fmtdate)
            sys.exit("#####UNRECOGNIZED SPECIES FOR AUGUSTUS AND NO READS\t" + now + "\t#####\n")
        # Prepare EVM input files
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(('\n###EVM STARTED AT:\t' + now + '\t###\n'))
        # HERE WE CONVERT FILES FOR EVM AND PLACE THEM IN INPUT FOLDER

        if not args.short_reads and not args.long_reads:
            if external_file:
                if external_file.endswith(fasta):
                    external_file_gff3 = mapping.gmap('ext', masked_ref, external_file, threads_use, 'gff3_gene',
                                                      args.min_intron_length, args.max_intron_length, args.end_exon,
                                                      gmap_wd, args.verbose, Fflag=True)
                    external_file_changed = update.external(external_file_gff3, gmap_wd, args.verbose)
                elif external_file.endswith("gff3"):
                    external_file_changed = update.external(external_file, gmap_wd, args.verbose)
                evm_inputs = {'augustus': augustus_gff3, 'genemark': genemark_gff3, 'AAT': merged_prot_gff3,
                              'external': external_file_changed}
            else:
                evm_inputs = {'augustus': augustus_gff3, 'genemark': genemark_gff3, 'AAT': merged_prot_gff3}
        elif args.short_reads or args.long_reads:
            if args.external:
                external_file = args.external
                if external_file.endswith(fasta):
                    external_file_gff3 = mapping.gmap('ext', masked_ref, external_file, threads_use, 'gff3_gene',
                                                      args.min_intron_length, args.max_intron_length, args.end_exon,
                                                      gmap_wd, args.verbose, Fflag=True)
                    external_file_changed = update.external(external_file_gff3, gmap_wd, args.verbose)
                elif external_file.endswith("gff3"):
                    external_file_changed = update.external(external_file, gmap_wd, args.verbose)
                evm_inputs = {'pasa': pasa_gff3, 'augustus': augustus_gff3, 'genemark': genemark_gff3,
                                  'AAT': merged_prot_gff3, 'gmap': trinity_path,'external': external_file_changed}
            else:
                evm_inputs = {'pasa': pasa_gff3, 'augustus': augustus_gff3, 'genemark': genemark_gff3,
                              'AAT': merged_prot_gff3, 'gmap': trinity_path}
        # HERE WE RUN EVM; WE PREPARE FILES THAT ARE REQUIRED BY EVM LIKE
        # WEIGTH TABLE

        list_soft, pred_file, transcript_file, protein_file = inputEvm.group_EVM_inputs(evm_inputs_dir, evm_inputs)
        weight_file = inputEvm.evm_weight(evm_inputs_dir, weights_dic, list_soft, pasa_name, gmap_name)
        # EVM PIPELINE
        round_n = 1

        if args.short_reads or args.long_reads:  # WE HAVE SHORT READS AND PROTEINS
            evm_gff3 = evmPipeline.evm_pipeline(evm_output_dir, threads_use, masked_ref, weight_file, pred_file,
                                                transcript_file, protein_file, args.segmentSize, args.overlap_size,
                                                args.verbose)
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('\n###UPDATE WITH PASA DATABASE STARTED AT:\t ' + now + '\t###\n'))
            round_n += 1
            final_output = pasa.update_database(threads_use, str(round_n), pasa_dir, args.pasa_db, ref_link, trinity_out,
                                                evm_gff3, args.verbose)
            final_update = grs.genename(final_output, args.prefix_gene, args.verbose, pasa_dir)
            final_update_stats = evmPipeline.gff3_stats(final_update, pasa_dir)
            final_update = grs.genename_evm(final_update, args.verbose, pasa_dir)
            final_files.append(final_update)
            final_files.append(final_update_stats)
            if "command" not in (iprscan_log.decode("utf-8")):
                annot, bad_models = iprscan.iprscan(ref_link, final_update, interproscan_out_dir, args.threads)
                final_files.append(annot)
                final_files.append(bad_models)
            if args.long_reads == '':
                final_output_dir = wd + 'output/'
                logistic.check_create_dir(final_output_dir)
                for filename in final_files:
                    if filename != '':
                        logistic.copy_file(filename, final_output_dir)
                cmdstring = "chmod -R 775 %s" % wd
                os.system(cmdstring)
                now = datetime.datetime.now().strftime(fmtdate)
                sys.exit("#####LOREAN FINISHED WITHOUT USING LONG READS\t" + now + "\t. GOOD BYE.#####\n")

        elif not args.short_reads and not args.long_reads:  # WE HAVE PROTEINS BUT NOT SHORT READS
            transcript_file = ''
            evm_gff3 = evmPipeline.evm_pipeline(evm_output_dir, threads_use, masked_ref, weight_file, pred_file,
                                                transcript_file, protein_file, args.segmentSize, args.overlap_size,
                                                args.verbose)
            last_gff3 = grs.genename(evm_gff3, args.prefix_gene, args.verbose, pasa_dir)
            final_update_stats = evmPipeline.gff3_stats(last_gff3, pasa_dir)
            now = datetime.datetime.now().strftime(fmtdate)
            final_files.append(last_gff3)
            final_files.append(final_update_stats)
            if "command" not in (iprscan_log.decode("utf-8")):
                annot, bad_models = iprscan.iprscan(ref_link, last_gff3, interproscan_out_dir, args.threads)
                final_files.append(annot)
                final_files.append(bad_models)
            final_output_dir = wd + 'output/'
            logistic.check_create_dir(final_output_dir)
            for filename in final_files:
                if filename != '':
                    logistic.copy_file(filename, final_output_dir)
            cmdstring = "chmod -R 775 %s" % wd
            os.system(cmdstring)
            now = datetime.datetime.now().strftime(fmtdate)
            sys.exit("##### EVM FINISHED AT:\t" + now + "\t#####\n")

        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(('\n###RUNNING iASSEMBLER\t' + now + '\t###\n'))

        if not long_sorted_bam:
            long_sam = mapping.gmap('sam', masked_ref, long_fasta, threads_use, 'samse',
                                    args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                    args.verbose, Fflag=False)
            long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, threads_use, wd, args.verbose)
            final_files.append(long_sorted_bam)

        # HERE WE MERGE THE GMAP OUTPUT WITH THE EVM OUTPUT TO HAVE ONE            # FILE
        # HERE WE CHECK IF WE HAVE THE PASA UPDATED FILE OR THE EVM
        # ORIGINAL FILE

        mergedmap_gff3 = logistic.catTwoBeds(long_sorted_bam, final_update, args.verbose, consensus_wd)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n\t###GFFREAD\t" + now + "\t###\n"))

        # HERE WE TRANSFORM THE COODINATES INTO SEQUENCES USING THE
        # REFERENCE
        gffread_fasta_file = consensus.gffread(mergedmap_gff3, masked_ref, consensus_wd, args.verbose)
        # HERE WE STORE THE SEQUENCE IN A DICTIONARY
        fake = []
        long_fasta = mseq.filterLongReads(gffread_fasta_file, args.assembly_overlap_length, args.max_long_read,
                                          consensus_wd, fake, threads_use, a=False)
        gffread_dict = consensus.fasta2Dict(gffread_fasta_file)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n\t#CLUSTERING\t" + now + "\t###\n"))

        # HERE WE CLUSTER THE SEQUENCES BASED ON THE GENOME POSITION
        cluster_list = consensus.cluster_pipeline(mergedmap_gff3, args.assembly_overlap_length, args.stranded, args.verbose)
        now = datetime.datetime.now().strftime(fmtdate)

        sys.stdout.write(("\n\t#CONSENSUS FOR EACH CLUSTER\t" + now + "\t###\n"))

        # HERE WE MAKE CONSENSUS FOR EACH CLUSTER
        tmp_wd = consensus_wd + 'tmp/'
        logistic.check_create_dir(tmp_wd)
        tmp_assembly_file = tmp_wd + 'assembly.fasta'
        if os.path.isfile(tmp_assembly_file):
            sys.stdout.write('No assembly')
        else:
            consensus.generate_fasta(cluster_list, gffread_dict, args.cluster_min_evidence,
                                     args.cluster_max_evidence, args.assembly_overlap_length, tmp_wd)
            consensus.assembly(args.assembly_overlap_length, args.assembly_percent_identity, threads_use, tmp_wd,
                               args.verbose)
            utrs.lengthSupport(tmp_wd, threads_use)

        # WITH THE ELSE, WE ALLOW THE USER TO DECIDE TO CHANGE THE ASSEMBLY
        # PARAMETERS AND COLLECT DIFFERENT ASSEMBLED SEQUENCES WITHOT RUNNING
        # THE FULL PIPELINE
        # HERE WE COLLECT THE ASSEMBLED SEQUENCES. WE COLLECT ONLY SEQUENCE
        # THAT PASS THE FILTER
        tmp_consensus = os.path.join(consensus_wd , 'tmp/')
        collect.parse_only(args.assembly_read_threshold, tmp_consensus, args.verbose)
        tmp_assembly = collect.cat_assembled(tmp_consensus)
        tmp_assembly_all = collect.cat_assembled_all(tmp_consensus)
        # HERE WE COLLECT THE NEW ASSEMBLED SEQUENCES AND WE COLLECT THE OLD
        # EVM DATA
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n###MAPPING CONSENSUS ASSEMBLIES\t" + now + "\t###\n"))

        # HERE WE MAP ALL THE FASTA FILES TO THE GENOME USING GMAP
        consensus_mapped_gff3 = mapping.gmap('cons', masked_ref, tmp_assembly, threads_use, 'gff3_gene',
                                             args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                             args.verbose, Fflag=True)

        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n###GETTING THE STRAND RIGHT\t" + now + "\t###\n"))

        strand_mapped_gff3 = grs.strand(final_update, consensus_mapped_gff3, masked_ref, threads_use, gmap_wd, args.verbose)
        merged_gff3 = collect.add_EVM(final_update, gmap_wd, strand_mapped_gff3)
        gff_pasa = grs.appendID(merged_gff3)
        no_overl = grs.removeOverlap(gff_pasa, args.verbose)
        no_disc = grs.removeDiscrepancy(no_overl, final_update, args.verbose)
        uniq_gene = grs.genename(no_disc, args.prefix_gene, args.verbose, exonerate_wd)

        update1 = grs.genename(uniq_gene, args.prefix_gene, args.verbose, gmap_wd)
        print(("\n###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))
        update2 = grs.exonerate(masked_ref, update1, threads_use, exonerate_wd, args.verbose)
        update3 = grs.genename(update2, args.prefix_gene, args.verbose, exonerate_wd)

        update4 = grs.add_removed_evm(update3, final_update, exonerate_wd)
        update4_1 = grs.genename(update4, args.prefix_gene, args.verbose, exonerate_wd)
        # HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN
        # PASA TO CORRECT SMALL ERRORS

        sys.stdout.write(("\n###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))

        fasta_all = logistic.cat_two_fasta(trinity_out, tmp_assembly_all, long_fasta, pasa_dir)
        round_n += 1

        update5 = pasa.update_database(threads_use, str(round_n), pasa_dir, args.pasa_db,  ref_link, fasta_all, update4_1,
                                       args.verbose)
        round_n += 1
        update6 = pasa.update_database(threads_use, str(round_n), pasa_dir, args.pasa_db,  ref_link, fasta_all, update5,
                                       args.verbose)
        final_update_update = grs.genename(update6, args.prefix_gene, args.verbose, pasa_dir)

        final_files.append(final_update_update)

        final_update_stats = evmPipeline.gff3_stats(final_update_update, pasa_dir)
        final_files.append(final_update_stats)

        if "command" not in (iprscan_log.decode("utf-8")):
            annot, bad_models = iprscan.iprscan(ref_link, final_update_update, interproscan_out_dir, args.threads)
            final_files.append(annot)
            final_files.append(bad_models)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(('\n###CREATING OUTPUT DIRECTORY\t' + now + '\t###\n'))

        final_output_dir = os.path.join(output_dir,  args.species + '_output' )

        logistic.check_create_dir(final_output_dir)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n##PLACING OUTPUT FILES IN OUTPUT DIRECTORY\t" + now + "\t###\n"))

        for filename in final_files:
            if filename != '':
                logistic.copy_file(filename, final_output_dir)
                cmdstring = "chmod -R 775 %s" % wd
                os.system(cmdstring)
        if not args.keep_tmp:
            temp_dir.cleanup()
        sys.exit("##### LOREAN FINISHED HERE. GOOD BYE. #####\n")
    else:
        sys.exit("#####LOREAN STOPS HERE. CHECK THAT THE PROTEIN AND SPECIES OPTION HAVE BOTH AN ARGUMENT. CHECK THAT "
                 "THE gm_key IS IN THE FOLDER#####\n")


if __name__ == '__main__':
    realstart = time.perf_counter()
    main()
    realend = time.perf_counter()
    realt = round(realend - realstart)
    m, s = divmod(realt, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    sys.stdout.write(
        "###LOREAN FINISHED WITHOUT ERRORS IN:" + " " + str(d) + " days " + str(h) + " hours " + str(m) + " min " + str(
            s) + " sec\t###\n")
