#! /usr/bin/env python3

import argparse


###############
###FUNCTIONS###
###############

def setting():
    '''Parses the arguments from the program invocation'''

    parser = argparse.ArgumentParser(
        prog='lorean',
        usage='%(prog)s [options] reference',
        description='LoReAn - Automated genome annotation pipeline that integrates long reads',
        epilog='Luigi Faino - March 2017')
    parser.add_argument("-pr", "--proteins",
                        nargs="?",
                        help="Path to protein sequences FASTA file []", default="")
    parser.add_argument("reference",
                        help="Path to reference file")
    parser.add_argument("-sp", "--species",
                        nargs="?",
                        help="Species name for AUGUSTUS training. No re-training if species already present in AUGUSTUS config folder", default="")
    parser.add_argument("-d","--stranded",
                        help="Run LoReAn on stranded mode [FALSE]",
                        action='store_true')
    parser.add_argument("-f","--fungus",
                        help="Use this option for fungal species (used in Gene Mark-ES)  [FALSE]",
                        action='store_true')
    parser.add_argument("-k","--keep_tmp",
                        help="Keep temporary files [FALSE]",
                        action='store_true')
    parser.add_argument("-sr","--short_reads",
                        nargs="?",
                        default="",
                        help="Path to short reads FASTQ. If paired end, comma-separated (1-1.fq,1-2.fq). BAM sorted files are allowed; the extension of the file should be filename.sorted.bam []",
                        metavar='FASTQ_file')
    parser.add_argument("-lr","--long_reads",
                        nargs="?", default="",
                        help="Path to long reads FASTQ []",
                        metavar='FASTQ_file')
    parser.add_argument("-a","--adapter",
                        nargs="?",
                        default="",
                        help="FASTA file containing the adapter sequences. Adapter sequences in forward and reverse strain of the same adapter need to be used in the file []",
                        metavar='FASTA_file')
    parser.add_argument("-rp","--repeat_masked",
                        nargs="?",
                        default="",
                        help="GFF or GFF3 or GTF or BED file containing repeats coordinates []",
                        metavar='GFF_file')
    parser.add_argument("-u","--update",
                        nargs="?",
                        default="",
                        help="GFF or GFF3 or GTF or BED file containing repeats coordinates []",
                        metavar='GFF_file')
    parser.add_argument("-m","--max_long_read",
                        nargs="?",
                        default=20000,
                        help="Filter out long reads longer than this value (longer reads may affect mapping and assembling) [20000]",
                        type=int)
    parser.add_argument("-pasa","--pasa_db",
                        nargs="?", default="annotation",
                        help="PASA database name [pipeline_run]")
    parser.add_argument("-n","--prefix_gene",
                        nargs="?",
                        default="species",
                        help="Prefix to add to the final Gff3 gene name [specie]")
    parser.add_argument("-w","--working_dir",
                        nargs="?",
                        default="annotation",
                        help="Working directory (will create if not present) [./]")
    parser.add_argument("-t","--threads",
                        nargs="?", default="3",
                        help="Number of threads [1]",
                        metavar='N')
    parser.add_argument("-cw","--augustus_weigth",
                        nargs="?",
                        default="1",
                        help="Weight assigned to AUGUSTUS evidence for EVM [1]",
                        metavar='N')
    parser.add_argument("-gw","--genemark_weigth",
                        nargs="?",
                        default="1",
                        help="Weight assigned to GENEMARK evidence for EVM [1]",
                        metavar='N')
    parser.add_argument("-tw","--trinity_weigth",
                        nargs="?",
                        default="1",
                        help="Weight assigned to Trinity mapped with GMAP evidence for EVM [1]",
                        metavar='N')
    parser.add_argument("-pw","--pasa_weigth",
                        nargs="?", default="5",
                        help="Weight assigned to PASA evidence for EVM [5]",
                        metavar='N')
    parser.add_argument("-aw","--AAT_weigth",
                        nargs="?",
                        default="1",
                        help="Weight assigned to AAT protein evidence for EVM [1]",
                        metavar='N')
    parser.add_argument("-c","--segmentSize",
                        nargs="?", default="100000",
                        help="Segment size for EVM partitions [100000]",
                        metavar='N')
    parser.add_argument("-e","--overlapSize",
                        nargs="?", default="10000",
                        help="Overlap size for EVM partitions [10000]",
                        metavar='N')
    parser.add_argument("-g","--min_intron_length",
                        nargs="?", default="9",
                        help="Minimal intron length for GMAP [9]",
                        metavar='N')
    parser.add_argument("-q","--max_intron_length",
                        nargs="?",
                        default="1000",
                        help="Maximal intron length for GMAP, STAR and TRINITY [1000]",
                        metavar='N')
    parser.add_argument("-ee", "--end_exon",
                        nargs="?", default="20",
                        help="Minimal length for end exon with GMAP [20]",
                        metavar='N')
    parser.add_argument("-cme","--cluster_min_evidence",
                        nargs="?", default="5",
                        help="Minimal evidence needed to form a cluster [5]",
                        metavar='N')
    parser.add_argument("-cMe","--cluster_max_evidence",
                        nargs="?",
                        default="5000",
                        help="Maximal evidence to form a cluster.Prevents the clustering or rRNA genes i.e. [5000]",
                        metavar='N')
    parser.add_argument("-aol","--assembly_overlapLength",
                        nargs="?",
                        default="200",
                        help="Minimal length (in nt) of overlap for ASSEMBLY [200]",
                        metavar='N')
    parser.add_argument("-api","--assembly_percentIdentity",
                        nargs="?", default="97",
                        help="Minimal identity for the ASSEMBLY (95-100) [97]",
                        metavar='N')
    parser.add_argument("-art","--assembly_readThreshold",
                        nargs="?",
                        default="0.3",
                        help="Fraction of reads supporting an assembled UNITIG to keep on the ASSEMBLY (0.1-1) [0.3]",
                        metavar='F')
    parser.add_argument("-v","--verbose",
                        help="Prints out the commands used in LoReAn[FALSE]",
                        action='store_true')

    args = parser.parse_args()
    return args
