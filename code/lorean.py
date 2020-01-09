#! /usr/bin/env python3


###############
###IMPORTS###
###############

# LIBRARIES
import datetime
import multiprocessing
import os
import random
import shutil
import string
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
    if args.pasa_db == "":
        pasadb = ''.join(random.sample(string.ascii_lowercase, 5))
    else:
        pasadb = args.pasa_db
    max_threads = multiprocessing.cpu_count()
    gmap_name = args.reference + '_GMAPindex'
    pasa_name = 'assembler-' + pasadb
    if args.upgrade == "":
        protein_loc = os.path.abspath(args.proteins)
    iprscan_log = iprscan.check_iprscan()
    # Useful variables for later
    root = os.getcwd()

    if args.out_dir != "" and args.out_dir.startswith("/"):
        output_dir = os.path.join(args.out_dir, "LoReAn" + args.working_dir)
    else:
        output_dir = os.path.join(root, "LoReAn_" + args.working_dir)

    logistic.check_create_dir(output_dir)


    if args.keep_tmp:
        wd = os.path.join(output_dir, "run/")
        logistic.check_create_dir(wd)
    elif args.verbose:
        wd = os.path.join(output_dir, "run/")
        logistic.check_create_dir(wd)
    else:
        temp_dir = tempfile.TemporaryDirectory(prefix='run_', dir=output_dir, suffix="/", )
        wd = temp_dir.name

    if args.upgrade == "":
        if not os.path.isfile(home + "/.gm_key"):
            sys.exit("#####LOREAN STOPS HERE. CHECK THAT THE gm_key IS IN THE HOME FOLDER#####\n")
        if args.proteins == "":
            sys.exit("#####LOREAN STOPS HERE. CHECK THAT THE PROTEIN OPTION IS SET#####\n")

    if args.stranded or args.adapter:
        if args.adapter == '':
            adapter_value = True
            sys.stdout.write('### RUNNING IN STRAND MODE AND FINDING ADAPTER AUTOMATICALLY ###\n')
            stranded_value = True
        else:
            adapter_value = args.adapter
            sys.stdout.write('### RUNNING IN STRAND MODE AND USING ADAPTER PROVIDED ###\n')
            stranded_value = True
    else:
        stranded_value = False
        sys.stdout.write('### RUNNING IN NON-STRAND MODE ###\n')
        adapter_value = False
    ref_orig = os.path.abspath(args.reference)
    ref_link = os.path.join(wd, args.reference)
    if not os.path.exists(ref_link):
        shutil.copyfile(ref_orig, ref_link)
    long_reads = args.long_reads
    fasta = (".fasta", ".fa", ".fas", ".fsta")
    fastq = (".fastq", ".fq")
    '''Core of the program'''
    # Parse the arguments
    if int(args.threads) > max_threads:
        threads_use = str(max_threads)
        sys.stdout.write(('### MAX NUMBER OF USED THREADS IS ' + str(max_threads) + ' AND NOT ' + args.threads + ' AS SET ###\n'))
    else:
        threads_use = args.threads
    if args.external:
        external_file = args.external
    else:
        external_file = ''
    if args.upgrade == "":
        if args.species == "":
            sys.exit("#####PLEASE DEFINE A SPECIES NAME\t" + now + "\t#####\n")
        else:
            if args.short_reads == '' and long_reads == '':
                if external_file.endswith("gff3") or external_file.endswith(fasta):
                    weights_dic = {'Augustus': args.augustus_weigth, 'GeneMark.hmm': args.genemark_weigth, 'exonerate': args.exonerate_weigth,
                                   'external' : args.external_weigth}
                else:
                    weights_dic = {'Augustus': args.augustus_weigth, 'GeneMark.hmm': args.genemark_weigth, 'exonerate': args.exonerate_weigth}
            elif args.short_reads != '' or long_reads != '':
                if external_file.endswith("gff3") or external_file.endswith(fasta):
                    weights_dic = {'Augustus': args.augustus_weigth, pasa_name: args.pasa_weigth, 'GeneMark.hmm': args.genemark_weigth,
                                   'exonerate': args.exonerate_weigth, gmap_name: args.trinity_weigth, 'external' : args.external_weigth}
                else:
                    weights_dic = {'Augustus': args.augustus_weigth, pasa_name: args.pasa_weigth, 'GeneMark.hmm': args.genemark_weigth,
                               'exonerate': args.exonerate_weigth, gmap_name: args.trinity_weigth}
    final_files = []  # STORE THE IMPORTANT OUTPUT FILES
    logistic.check_create_dir(wd)
    logistic.check_file(ref_link)
    gmap_wd = os.path.join(wd ,'gmap_output/')
    exonerate_wd = os.path.join(wd , 'exonerate')
    pasa_dir = os.path.join(wd , 'PASA/')
    star_out = os.path.join(wd , 'STAR/')
    trin_dir = os.path.join(wd , 'Trinity/')
    evm_inputs_dir = os.path.join(wd , 'evm_inputs/')
    braker_folder = os.path.join(wd , 'braker/')
    evm_output_dir = os.path.join(wd , 'evm_output/')
    interproscan_out_dir = os.path.join(wd , 'interproscan')
    wd_split = os.path.join(wd , 'split/')
    logistic.check_create_dir(wd_split)
    logistic.check_create_dir(evm_inputs_dir)
    logistic.check_create_dir(evm_output_dir)
    logistic.check_create_dir(trin_dir)
    logistic.check_create_dir(star_out)
    logistic.check_create_dir(pasa_dir)
    logistic.check_create_dir(gmap_wd)
    logistic.check_create_dir(exonerate_wd)
    if args.interproscan:
        logistic.check_create_dir(interproscan_out_dir)
    if long_reads:
        consensus_wd = os.path.join(wd , 'consensus/')
        logistic.check_create_dir(consensus_wd)
    if long_reads != "" or args.short_reads != "":
        logistic.check_gmap(threads_use, 'samse', args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                        args.verbose)
    augustus_species, err_augustus = logistic.augustus_species_func(home)
    if args.repeat_masked != "":
        sys.stdout.write(('###MASKING THE GENOME STARTED AT:\t' + now + '\t###\n'))
        masked_ref = mseq.maskedgenome(wd_split, ref_link, args.repeat_masked, args.repeat_lenght, args.verbose)
    elif args.mask_genome:
        sys.stdout.write(('###RUNNNG REPEATSCOUT AND REPEATMASK TO MASK THE GENOME STARTED AT:\t' + now + '\t###\n'))
        masked_ref = mseq.repeatsfind(ref_link, wd_split, args.repeat_lenght, threads_use, args.verbose)
    else:
        masked_ref = ref_link
    list_fasta_names, dict_ref_name, ref_rename = multiple.single_fasta(masked_ref, wd_split)
    if args.short_reads or long_reads:
        if int(threads_use) > 1:
            trinity_cpu = int(int(threads_use) / int(2))
        else:
            trinity_cpu = int(threads_use)
        now = datetime.datetime.now().strftime(fmtdate)
        # SHORT READS
        if args.short_reads.endswith(fastq):
            sys.stdout.write(('###STAR MAPPING  STARTED AT:\t' + now + '\t###\n'))
            if ',' in args.short_reads:
                paired_end_files = args.short_reads.split(',')
                short_1 = os.path.abspath(paired_end_files[0])
                short_2 = os.path.abspath(paired_end_files[1])
                short_reads_file = [short_1, short_2]
            else:
                short_reads_file = os.path.abspath(args.short_reads)
            # Map with STAR
            short_bam = mapping.star(ref_rename, short_reads_file, threads_use, args.max_intron_length, star_out,
                                     args.verbose)
            short_sorted_bam = mapping.samtools_sort(short_bam, threads_use, wd, args.verbose)
            final_mapping_star = mapping.change_chr(short_sorted_bam, dict_ref_name, star_out, threads_use, args.verbose, "short")
            default_bam = short_sorted_bam
            # Keep the output
            final_files.append(final_mapping_star)
            # TRANSCRIPT ASSEMBLY
            # TRINITY
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('###TRINITY STARTS AT:\t' + now + '\t###\n'))
            trinity_out = transcripts.trinity(short_sorted_bam, trin_dir, args.max_intron_length, trinity_cpu, args.verbose)
            if args.upgrade == "":
                trinity_gff3 = mapping.gmap('trin', ref_rename, trinity_out, threads_use, 'gff3_gene',
                                    args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                    args.verbose, Fflag=True)
                trinity_path = trinity_gff3
            long_sorted_bam = False
        # BAM SORTED FILES GET IN HERE
        elif args.short_reads.endswith("bam") or long_reads.endswith("bam"):
            logistic.check_create_dir(star_out)
            if args.short_reads.endswith("bam"):
                map_reads = os.path.abspath(args.short_reads)
                short_sorted_bam = mapping.change_chr_to_seq(map_reads, dict_ref_name, star_out, threads_use, args.verbose)
            else:
                map_reads = os.path.abspath(long_reads)
                short_sorted_bam = mapping.change_chr_to_seq(map_reads, dict_ref_name, star_out, threads_use, args.verbose)
                mapping.samtools_index(short_sorted_bam, star_out, args.verbose)
                long_reads = transcripts.bamtofastq(short_sorted_bam, args.verbose)
            #short_sorted_bam = os.path.abspath(args.short_reads)
            default_bam = short_sorted_bam
            # TRANSCRIPT ASSEMBLY
            # TRINITY
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('###TRINITY STARTS AT:\t' + now + '\t###\n'))
            trinity_out = transcripts.trinity(short_sorted_bam, trin_dir, args.max_intron_length, trinity_cpu, args.verbose)
            if args.upgrade == "":
                trinity_gff3 = mapping.gmap('trin', ref_rename, trinity_out, threads_use, 'gff3_gene',
                                        args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                        args.verbose, Fflag=True)
                trinity_path = trinity_gff3
            long_sorted_bam = False
        # LONG READS
        elif long_reads.endswith(fastq) or long_reads.endswith(fasta):
            # with this operation, reads are filtered for their length.
            # Nanopore reads can be chimeras or sequencing artefacts.
            # filtering on length reduces the amount of sequencing
            # artefacts
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(("###FILTERING OUT LONG READS STARTED AT:\t" + now + "\t###\n"))
            long_fasta, stranded_value = mseq.filterLongReads(long_reads, args.assembly_overlap_length, args.max_long_read, gmap_wd,
                                              adapter_value, threads_use, args.adapter_match_score, ref_rename,
                                              args.max_intron_length, args.verbose, stranded_value)


                # If short reads have been mapped dont do it
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('###GMAP\t' + now + 't###\n'))
            if args.minimap2:
                long_sam = mapping.minimap(ref_rename, long_fasta, threads_use, args.max_intron_length, gmap_wd, args.verbose)
            else:
                long_sam = mapping.gmap('sam', ref_rename, long_fasta, threads_use, 'samse',
                                    args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                    args.verbose, Fflag=False)

            # Convert to sorted BAM
            long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, threads_use, gmap_wd, args.verbose)
            sam_orig_id = mapping.change_chr(long_sorted_bam, dict_ref_name, gmap_wd, threads_use, args.verbose, "long")
            default_bam = long_sorted_bam
            # Keep the output

            final_files.append(sam_orig_id)
            # TRANSCRIPT ASSEMBLY
            # TRINITY
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('###TRINITY STARTS AT:\t' + now + '\t###\n'))
            trinity_out = transcripts.trinity(long_sorted_bam, trin_dir, args.max_intron_length, trinity_cpu, args.verbose)
            if args.upgrade == "":
                trinity_gff3 = mapping.gmap('trin', ref_rename, trinity_out, threads_use, 'gff3_gene',
                                        args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                        args.verbose, Fflag=True)
                trinity_path = trinity_gff3
        else:
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('###NO LONG READS FILE OR SHORT READS\t' + now + '\t###\n'))
        # PASA Pipeline
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(('###PASA STARTS AT:\t' + now + '\t###\n'))
        # Create PASA folder and configuration file
        #align_pasa_conf = pasa.pasa_configuration(pasa_dir, pasadb, args.verbose)
        # Launch PASA

        if args.upgrade == "":
            if os.path.isfile(home + "/.gm_key") and args.proteins != "":
                pasa_gff3 = pasa.pasa_call(pasa_dir, pasadb, ref_rename, trinity_out, args.max_intron_length,
                                           threads_use, args.verbose)
        # HERE WE PARALLELIZE PROCESSES WHEN MULTIPLE THREADS ARE USED
                if args.species in (err_augustus.decode("utf-8")) or args.species in augustus_species:
                    now = datetime.datetime.now().strftime(fmtdate)
                    sys.stdout.write(('###AUGUSTUS, GENEMARK-ES AND AAT STARTED AT:' + now + '\t###\n'))
                    queue = Queue()
                    for software in range(3):
                        queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                        for software in range(3):
                            t = Thread(target=handler.august_gmes_exonerate, args=(queue, ref_rename, args.species, protein_loc,
                                                                           threads_use, args.fungus, list_fasta_names, wd, exonerate_wd,
                                                                           args.verbose))
                            t.daemon = True
                            t.start()
                    queue.join()
                    augustus_file = wd + 'augustus/augustus.gff'
                    augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                    genemark_file = wd + 'gmes/genemark.gtf'
                    genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                    merged_prot_gff3 = wd + 'exonerate/protein_evidence.gff3'

                elif args.short_reads or long_reads:  # USING PROTEINS AND SHORT READS
                    logistic.check_create_dir(braker_folder)
                    now = datetime.datetime.now().strftime(fmtdate)
                    sys.stdout.write(('###BRAKER1 (USING SHORT READS) AND AAT STARTED AT:\t' + now + '\t###\n'))
                    queue = Queue()
                    for software in range(2):
                        queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                        for software in range(2):
                            t = Thread(target=handler.braker_exonerate, args=(queue, ref_rename, default_bam, args.species, protein_loc,
                                                                       threads_use, args.fungus, wd,
                                                                        braker_folder, exonerate_wd, args.verbose))
                            t.daemon = True
                            t.start()
                    queue.join()
                    augustus_file, genemark_file = inputEvm.braker_folder_find(braker_folder)
                    augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                    genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                    merged_prot_gff3 = wd + 'exonerate/protein_evidence.gff3'

                else:  # USING PROTEINS AND LONG READS
                    queue = Queue()
                    now = datetime.datetime.now().strftime(fmtdate)
                    sys.stdout.write(('###BRAKER1 (USING LONG READS) AND AAT STARTED AT: \t' + now + '\t###\n'))
                    logistic.check_create_dir(braker_folder)
                    for software in range(2):
                        queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                        for software in range(2):
                            t = Thread(target=handler.braker_exonerate,
                                       args=(queue, ref_rename, long_sorted_bam, args.species, protein_loc,
                                             threads_use, args.fungus, wd, braker_folder, exonerate_wd, args.verbose))
                            t.daemon = True
                            t.start()
                    queue.join()
                    augustus_file, genemark_file = inputEvm.braker_folder_find(braker_folder)
                    augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
                    genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
                    merged_prot_gff3 = wd + 'exonerate/protein_evidence.gff3'
    elif args.species in (err_augustus.decode("utf-8")) or args.species in augustus_species or args.species != "" or args.upgrade != "":
        if os.path.isfile(home + "/.gm_key") and args.proteins != "":
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('###AUGUSTUS, GENEMARK-ES AND AAT STARTED AT:' + now + '\t###\n'))
            queue = Queue()
            for software in range(3):
                queue.put(software)  # QUEUE WITH A ZERO AND A ONE
                for software in range(3):
                    t = Thread(target=handler.august_gmes_exonerate, args=(queue, ref_rename, args.species, protein_loc,
                                                                   threads_use, args.fungus, list_fasta_names, wd, exonerate_wd,
                                                                   args.verbose))
                    t.daemon = True
                    t.start()
            queue.join()
            augustus_file = wd + 'augustus/augustus.gff'
            augustus_gff3 = inputEvm.convert_augustus(augustus_file, wd)
            genemark_file = wd + 'gmes/genemark.gtf'
            genemark_gff3 = inputEvm.convert_genemark(genemark_file, wd)
            merged_prot_gff3 = wd + 'exonerate/protein_evidence.gff3'
    else:
        now = datetime.datetime.now().strftime(fmtdate)
        sys.exit("#####UNRECOGNIZED SPECIES FOR AUGUSTUS AND NO READS\t" + now + "\t#####\n")
    # Prepare EVM input files
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(('###EVM STARTED AT:\t' + now + '\t###\n'))
    # HERE WE CONVERT FILES FOR EVM AND PLACE THEM IN INPUT FOLDER
    round_n = 0
    if args.upgrade == "":
        if not args.short_reads and not long_reads:
            if external_file:
                if external_file.endswith(fasta):
                    external_file_gff3 = mapping.gmap('ext', ref_rename, external_file, threads_use, 'gff3_gene',
                                                      args.min_intron_length, args.max_intron_length, args.end_exon,
                                                      gmap_wd, args.verbose, Fflag=True)
                    external_file_changed = update.external(external_file_gff3, gmap_wd, args.verbose)
                elif external_file.endswith("gff3"):
                    external_file_changed = update.external(external_file, gmap_wd, args.verbose)
                evm_inputs = {'augustus': augustus_gff3, 'genemark': genemark_gff3, 'exonerate': merged_prot_gff3,
                              'external': external_file_changed}
            else:
                evm_inputs = {'augustus': augustus_gff3, 'genemark': genemark_gff3, 'exonerate': merged_prot_gff3}
        elif args.short_reads or long_reads:
            if args.external:
                external_file = args.external
                if external_file.endswith(fasta):
                    external_file_gff3 = mapping.gmap('ext', ref_rename, external_file, threads_use, 'gff3_gene',
                                                      args.min_intron_length, args.max_intron_length, args.end_exon,
                                                      gmap_wd, args.verbose, Fflag=True)
                    external_file_changed = update.external(external_file_gff3, gmap_wd, args.verbose)
                elif external_file.endswith("gff3"):
                    external_file_changed = update.external(external_file, gmap_wd, args.verbose)
                evm_inputs = {'pasa': pasa_gff3, 'augustus': augustus_gff3, 'genemark': genemark_gff3,
                                  'exonerate': merged_prot_gff3, 'gmap': trinity_path,'external': external_file_changed}
            else:
                evm_inputs = {'pasa': pasa_gff3, 'augustus': augustus_gff3, 'genemark': genemark_gff3,
                              'exonerate': merged_prot_gff3, 'gmap': trinity_path}
        # HERE WE RUN EVM; WE PREPARE FILES THAT ARE REQUIRED BY EVM LIKE
        # WEIGTH TABLE

        list_soft, pred_file, transcript_file, protein_file = inputEvm.group_EVM_inputs(evm_inputs_dir, evm_inputs)
        weight_file = inputEvm.evm_weight(evm_inputs_dir, weights_dic, list_soft, pasa_name, gmap_name)
        # EVM PIPELINE


        if args.short_reads or long_reads:  # WE HAVE SHORT READS AND PROTEINS
            evm_gff3 = evmPipeline.evm_pipeline(evm_output_dir, threads_use, ref_rename, weight_file, pred_file,
                                                transcript_file, protein_file, args.segmentSize, args.overlap_size,
                                                args.verbose)
            final_evm = grs.genename_evm(evm_gff3, args.verbose, evm_output_dir, dict_ref_name, args.upgrade)
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(('###UPDATE WITH PASA DATABASE STARTED AT:\t ' + now + '\t###\n'))
            round_n += 1
            final_output = pasa.update_database(threads_use, str(round_n), pasa_dir, pasadb, ref_rename, trinity_out,
                                              final_evm, args.verbose)
            if long_reads == '':
                final_update_all = grs.genename_last(final_output, args.prefix_gene, args.verbose, pasa_dir, dict_ref_name, "pasa")
                final_update_stats = evmPipeline.gff3_stats(final_update_all, pasa_dir)
                final_files.append(final_update_all)
                final_files.append(final_update_stats)
                if "command" not in (iprscan_log.decode("utf-8")) and args.interproscan:
                    annot, bad_models = iprscan.iprscan(masked_ref, final_update_all, interproscan_out_dir, args.threads)
                    final_files.append(annot)
                    final_files.append(bad_models)
                final_output_dir = os.path.join(output_dir, args.species + '_output')
                logistic.check_create_dir(final_output_dir)
                for filename in final_files:
                    if filename != '':
                        logistic.copy_file(filename, final_output_dir)
                cmdstring = "chmod -R 775 %s" % wd
                os.system(cmdstring)
                now = datetime.datetime.now().strftime(fmtdate)
                sys.exit("#####LOREAN FINISHED WITHOUT USING LONG READS\t" + now + "\t. GOOD BYE.#####\n")

            else:
                final_keep = grs.genename_last(final_output, args.prefix_gene, args.verbose, pasa_dir, dict_ref_name, "pasa")
                final_keep_stats = evmPipeline.gff3_stats(final_keep, pasa_dir)
                final_files.append(final_keep)
                final_files.append(final_keep_stats)
        elif not args.short_reads and not long_reads:  # WE HAVE PROTEINS BUT NOT SHORT READS
            transcript_file = ''
            evm_gff3 = evmPipeline.evm_pipeline(evm_output_dir, threads_use, ref_rename, weight_file, pred_file,
                                                transcript_file, protein_file, args.segmentSize, args.overlap_size,
                                                args.verbose)
            final_update_all = grs.genename_last(evm_gff3, args.prefix_gene, args.verbose, pasa_dir, dict_ref_name, "pasa")
            final_update_stats = evmPipeline.gff3_stats(final_update_all, pasa_dir)
            final_files.append(final_update_all)
            final_files.append(final_update_stats)
            now = datetime.datetime.now().strftime(fmtdate)
            if "command" not in (iprscan_log.decode("utf-8")) and args.interproscan:
                annot, bad_models = iprscan.iprscan(masked_ref, final_update_all, interproscan_out_dir, args.threads)
                final_files.append(annot)
                final_files.append(bad_models)
            final_output_dir = os.path.join(output_dir, args.species + '_output')
            logistic.check_create_dir(final_output_dir)
            for filename in final_files:
                if filename != '':
                    logistic.copy_file(filename, final_output_dir)
            cmdstring = "chmod -R 775 %s" % wd
            os.system(cmdstring)
            now = datetime.datetime.now().strftime(fmtdate)
            sys.exit("##### EVM FINISHED AT:\t" + now + "\t#####\n")
    else:
        final_evm = grs.genename_evm(args.upgrade, args.verbose, evm_output_dir, dict_ref_name, args.upgrade)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(('###UPDATE WITH PASA DATABASE STARTED AT:\t ' + now + '\t###\n'))
        round_n += 1
        final_output = pasa.update_database(threads_use, str(round_n), pasa_dir, pasadb, ref_rename, trinity_out,
                                            final_evm, args.verbose)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(('###RUNNING iASSEMBLER\t' + now + '\t###\n'))

    if not long_sorted_bam:
        #print("line 430")
        long_fasta, stranded_value = mseq.filterLongReads(long_reads, args.assembly_overlap_length, args.max_long_read, gmap_wd,
                                          adapter_value, threads_use, args.adapter_match_score, ref_rename,
                                              args.max_intron_length, args.verbose, stranded_value)
        if args.minimap2:
            long_sam = mapping.minimap(ref_rename, long_fasta, threads_use, args.max_intron_length, gmap_wd, args.verbose)
        else:
            long_sam = mapping.gmap('sam', ref_rename, long_fasta, threads_use, 'samse',
                                    args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                    args.verbose, Fflag=False)
        long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, threads_use, wd, args.verbose)
        sam_orig_id = mapping.change_chr(long_sorted_bam, dict_ref_name, gmap_wd, threads_use, args.verbose, "long")
        final_files.append(sam_orig_id)

    # HERE WE MERGE THE GMAP OUTPUT WITH THE EVM OUTPUT TO HAVE ONE            # FILE
    # HERE WE CHECK IF WE HAVE THE PASA UPDATED FILE OR THE EVM
    # ORIGINAL FILE

    mergedmap_gff3 = logistic.catTwoBeds(long_sorted_bam, final_evm, args.verbose, consensus_wd)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("\t###GFFREAD\t" + now + "\t###\n"))

    # HERE WE TRANSFORM THE COODINATES INTO SEQUENCES USING THE
    # REFERENCE
    gffread_fasta_file = consensus.gffread(mergedmap_gff3, ref_rename, consensus_wd, args.verbose)
    # HERE WE STORE THE SEQUENCE IN A DICTIONARY

    gffread_dict = consensus.fasta2Dict(gffread_fasta_file)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("\t#CLUSTERING\t" + now + "\t###\n"))

    # HERE WE CLUSTER THE SEQUENCES BASED ON THE GENOME POSITION
    cluster_list = consensus.cluster_pipeline(mergedmap_gff3, stranded_value, args.verbose)
    now = datetime.datetime.now().strftime(fmtdate)

    sys.stdout.write(("\t#CONSENSUS FOR EACH CLUSTER\t" + now + "\t###\n"))

    # HERE WE MAKE CONSENSUS FOR EACH CLUSTER
    tmp_wd = consensus_wd + 'tmp/'
    logistic.check_create_dir(tmp_wd)
    tmp_assembly_file = tmp_wd + 'assembly.fasta'
    if os.path.isfile(tmp_assembly_file):
        sys.stdout.write('No assembly')
    else:
        consensus.generate_fasta(cluster_list, gffread_dict, args.cluster_min_evidence,
                                 args.cluster_max_evidence, args.assembly_overlap_length, stranded_value, tmp_wd)
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
    sys.stdout.write(("###MAPPING CONSENSUS ASSEMBLIES\t" + now + "\t###\n"))

    # HERE WE MAP ALL THE FASTA FILES TO THE GENOME USING GMAP
    consensus_mapped_gff3 = mapping.gmap('cons', ref_rename, tmp_assembly, threads_use, 'gff3_gene',
                                         args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                         args.verbose, Fflag=True)

    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("###GETTING THE STRAND RIGHT\t" + now + "\t###\n"))
    merged_gff3 = collect.add_EVM(final_output, gmap_wd, consensus_mapped_gff3)
    update2 = grs.exonerate(ref_rename, merged_gff3, threads_use, exonerate_wd, args.verbose)
    update3 = grs.genename_lorean(update2, args.verbose, exonerate_wd)
    # HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN
    # PASA TO CORRECT SMALL ERRORS
    sys.stdout.write(("###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))
    fasta_all = logistic.cat_two_fasta(trinity_out, tmp_assembly_all, long_fasta, pasa_dir)
    round_n += 1
    update5 = pasa.update_database(threads_use, str(round_n), pasa_dir, pasadb,  ref_rename, fasta_all,
                                   update3, args.verbose)
    if args.verbose:
        sys.stdout.write(update5)
    round_n += 1
    update6 = pasa.update_database(threads_use, str(round_n), pasa_dir, pasadb,  ref_rename, fasta_all,
                                   update5, args.verbose)
    if args.verbose:
        sys.stdout.write(update6)
    final_update_update = grs.genename_last(update6, args.prefix_gene, args.verbose, pasa_dir, dict_ref_name, "lorean")
    final_files.append(final_update_update)

    final_update_stats = evmPipeline.gff3_stats(final_update_update, pasa_dir)
    final_files.append(final_update_stats)

    if "command" not in (iprscan_log.decode("utf-8")) and args.interproscan:
        annot, bad_models = iprscan.iprscan(masked_ref, final_update_update, interproscan_out_dir, args.threads)
        final_files.append(annot)
        final_files.append(bad_models)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(('###CREATING OUTPUT DIRECTORY\t' + now + '\t###\n'))

    final_output_dir = os.path.join(output_dir,  args.species + '_output')

    logistic.check_create_dir(final_output_dir)
    now = datetime.datetime.now().strftime(fmtdate)
    sys.stdout.write(("##PLACING OUTPUT FILES IN OUTPUT DIRECTORY\t" + now + "\t###\n"))

    for filename in final_files:
        if filename != '':
            logistic.copy_file(filename, final_output_dir)
            cmdstring = "chmod -R 775 %s" % wd
            os.system(cmdstring)
    if not args.keep_tmp:
        temp_dir.cleanup()
    sys.exit("##### LOREAN FINISHED HERE. GOOD BYE. #####\n")




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
