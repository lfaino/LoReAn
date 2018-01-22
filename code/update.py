#!/usr/bin/env python3

import datetime
import multiprocessing
import os
import subprocess
import sys
import tempfile

# OTHER SCRIPTS
import arguments as arguments
import collectOnly as collect
import consensusIAssembler as consensus
import dirsAndFiles as logistic
import getRightStrand as grs
import manipulateSeq as mseq
import mapping
import pasa as pasa
import reduceUTRs as utrs
import transcriptAssembly as transcripts

#======================================================================================================================

GT_GFF3 = 'gt gff3 -sort -tidy -setsource external %s'

#======================================================================================================================


def external(external, wd, verbose):
    error = tempfile.NamedTemporaryFile(delete=False, mode = "w", prefix="external.", suffix= ".error3")
    tmp_file = tempfile.NamedTemporaryFile(delete=False, mode = "w", prefix="external.", suffix= ".gff3")
    cmd = GT_GFF3 % external
    print(cmd)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % cmd)
        sys.stderr.write('Log file is: %s\n\n' % error.name)
    gt_call = subprocess.Popen(cmd, stdout=tmp_file, stderr=error, cwd=wd, shell=True)
    gt_call.communicate()
    return tmp_file.name


def upgrade():
    '''Core of the program'''

    args = arguments.setting()
    fasta = (".fasta", ".fa", ".fas", ".fsta")
    fastq = (".fastq", ".fq")
    fmtdate = '%H:%M:%S %d-%m'
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
    ref = os.path.join(wd, args.reference)
    if not os.path.exists(ref):
        os.symlink(ref_orig, ref)

    max_threads = multiprocessing.cpu_count()
    if int(args.threads) > max_threads:
        threads_use = str(max_threads)
        sys.stdout.write(('\n### MAX NUMBER OF USED THREADS IS ' + str(max_threads) + ' AND NOT ' + args.threads + ' AS SET ###\n'))
    else:
        threads_use = args.threads

    final_files = []  # STORE THE IMPORTANT OUTPUT FILES

    logistic.check_create_dir(wd)
    logistic.check_file(ref)

    gmap_wd = wd + '/gmap_output/'
    exonerate_wd = wd + '/exonerate/'
    pasa_dir = wd + 'PASA/'
    star_out = wd + '/STAR/'
    trin_dir = wd + '/Trinity/'

    logistic.check_create_dir(trin_dir)
    logistic.check_create_dir(star_out)
    logistic.check_create_dir(pasa_dir)
    logistic.check_create_dir(gmap_wd)
    logistic.check_create_dir(exonerate_wd)
    if args.long_reads:
        consensus_wd = (wd + '/consensus/')
        logistic.check_create_dir(consensus_wd)

    logistic.check_gmap(threads_use, 'samse', args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd, args.verbose)

    if args.repeat_masked:
        genome_gmap = mseq.maskedgenome(gmap_wd, ref, args.repeat_masked)
    else:
        genome_gmap = ref

    if args.short_reads or args.long_reads:
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(('\n###STAR MAPPING  STARTED AT:\t' + now + '\t###\n'))
        if args.short_reads.endswith(fastq):
            if ',' in args.short_reads:
                pairedEndFiles = args.short_reads.split(',')
                short_1 = os.path.abspath(pairedEndFiles[0])
                short_2 = os.path.abspath(pairedEndFiles[1])
                short_reads_file = [short_1, short_2]
            else:
                short_reads_file = os.path.abspath(args.short_reads)
            short_bam = mapping.star(ref, short_reads_file, threads_use, args.max_intron_length, star_out,
                                     args.verbose)
            short_sorted_bam = mapping.samtools_sort(short_bam, threads_use, wd, args.verbose)
            final_files.append(short_sorted_bam)
        # BAM SORTED FILES GET IN HERE
        elif args.short_reads.endswith("bam"):
            logistic.check_create_dir(star_out)
            short_sorted_bam = os.path.abspath(args.short_reads)
            bam_file = short_sorted_bam.split("/")
            short_bam = star_out + "/" + bam_file[-1]
            if not os.path.exists(ref):
                os.symlink(short_sorted_bam, short_bam)
        else:
            short_sorted_bam = False
            sys.stdout.write('No short reads file')
        if args.long_reads.endswith(fastq) or args.long_reads.endswith(fasta):
            now = datetime.datetime.now().strftime(fmtdate)
            sys.stdout.write(("\n###FILTERING OUT LONG READS STARTED AT:\t" + now + "\t###\n"))
            long_fasta = mseq.filterLongReads(args.long_reads, args.assembly_overlap_length,
                                                            args.max_long_read, gmap_wd, args.adapter, threads_use,
                                                            a=True)
            if not short_sorted_bam:
                # If short reads have been mapped dont do it
                now = datetime.datetime.now().strftime(fmtdate)
                sys.stdout.write(('\n###GMAP\t' + now + 't###\n'))
                long_sam = mapping.gmap('sam', genome_gmap, long_fasta, threads_use, 'samse',
                                        args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                        args.verbose, Fflag=False)
                long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, threads_use, wd, args.verbose)
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

    else:
        sys.exit("### NO READS TO USE ###")

    if args.long_reads:
        if not long_sorted_bam:
            long_sam = mapping.gmap('sam', genome_gmap, long_fasta, threads_use, 'samse',
                                    args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                    args.verbose, Fflag=False)
            long_sorted_bam = mapping.sam_to_sorted_bam(long_sam, threads_use, wd, args.verbose)
            final_files.append(long_sorted_bam)
        file_name = consensus_wd + 'mergedGmapEvm.beforeAssembly.gff3'

        mergedmapGFF3 = logistic.catTwoBeds(long_sorted_bam, args.upgrade, file_name, args.verbose)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n\t###GFFREAD\t" + now + "\t###\n"))

        gffread_fasta_file = consensus.gffread(mergedmapGFF3, ref, consensus_wd, args.verbose)
        # HERE WE STORE THE SEQUENCE IN A DICTIONARY
        fake = []
        long_fasta = mseq.filterLongReads(gffread_fasta_file, args.assembly_overlap_length,
                                                        args.max_long_read, consensus_wd, fake, threads_use,
                                                        a=False)

        gffreadDict = consensus.fasta2Dict(gffread_fasta_file)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n\t#CLUSTERING\t" + now + "\t###\n"))

        cluster_list = consensus.cluster_pipeline(mergedmapGFF3, args.assembly_overlap_length, args.stranded, args.verbose)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n\t#CONSENSUS FOR EACH CLUSTER\t" + now + "\t###\n"))
        tmp_wd = consensus_wd + 'tmp/'
        logistic.check_create_dir(tmp_wd)
        tmp_assembly_file = tmp_wd + 'assembly.fasta'
        if os.path.isfile(tmp_assembly_file):
            sys.stdout.write('No assembly')
        else:
            consensus.generate_fasta(cluster_list, gffreadDict, args.cluster_min_evidence,
                                     args.cluster_max_evidence, args.assembly_overlap_length, tmp_wd)
            consensus.assembly(args.assembly_overlap_length, args.assembly_percent_identity, threads_use, tmp_wd,
                               args.verbose)
            utrs.lengthSupport(tmp_wd, threads_use)
        tmp_consensus = os.path.join(consensus_wd , 'tmp/')
        collect.parse_only(args.assembly_read_threshold, tmp_consensus, args.verbose)
        tmp_assembly = collect.cat_assembled(tmp_consensus)
        merged_fasta_filename = consensus_wd + 'assembly.wEVM.fasta'
        collect.add_EVM(gffread_fasta_file, tmp_assembly, merged_fasta_filename)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n###MAPPING CONSENSUS ASSEMBLIES\t" + now + "\t###\n"))
        consensus_mapped_gff3 = mapping.gmap('cons', genome_gmap, merged_fasta_filename, threads_use, 'gff3_gene',
                                             args.min_intron_length, args.max_intron_length, args.end_exon, gmap_wd,
                                             args.verbose, Fflag=True)
        now = datetime.datetime.now().strftime(fmtdate)
        sys.stdout.write(("\n###GETTING THE STRAND RIGHT\t" + now + "\t###\n"))
        strand_mapped_gff3 = grs.strand(args.upgrade, consensus_mapped_gff3, ref, threads_use, gmap_wd, args.verbose)
        gff_pasa = grs.appendID(strand_mapped_gff3)
        no_overl = grs.removeOverlap(gff_pasa, args.verbose)
        no_disc = grs.removeDiscrepancy(no_overl, args.upgrade, args.verbose)
        uniq_gene = grs.newNames(no_disc)

        finalupdate3 = grs.genename(uniq_gene, args.prefix_gene, args.verbose)
        print(("\n###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))
        finalupdate4 = grs.exonerate(ref, finalupdate3, threads_use, exonerate_wd, args.verbose)
        finalupdate5 = grs.genename(finalupdate4, args.prefix_gene, args.verbose)

        # HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN
        # PASA TO CORRECT SMALL ERRORS

        sys.stdout.write(("\n###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))
        round_n = 0
        fasta_all = logistic.cat_two_fasta(trinity_out, merged_fasta_filename, long_fasta, pasa_dir)
        round_n += 1
        pasa.create_pasa_database(pasa_dir, args.pasa_db, args.verbose)
        align_pasa_conf = pasa.pasa_configuration(pasa_dir, args.pasa_db, args.verbose)
        finalupdate = pasa.update_database(threads_use, str(round_n), pasa_dir, args.pasa_db, align_pasa_conf, ref,
                                           fasta_all, finalupdate5, args.verbose)
        round_n += 1
        finalupdate2 = pasa.update_database(threads_use, str(round_n), pasa_dir, args.pasa_db, align_pasa_conf, ref,
                                            fasta_all, finalupdate, args.verbose)
        final_update = grs.genename(finalupdate2, args.prefix_gene, args.verbose)

        final_files.append(final_update)

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