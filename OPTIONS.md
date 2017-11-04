# Minimal options are:

## lorean.py -pr <protein.fasta> -sp \<species> \<genome.fasta>

### The specie name used in LoReAn needs to be in the AUGUSTUS database. 
### However, if RNA-seq is supplied, the tool will calculate the best AUGUSTUS sets using BRAKER1 software

### Other options are:

#### -h, --help
     show help message and exit
#### -pr [PROTEINS], --proteins [PROTEINS]
    Path to protein sequences FASTA file []
#### -sp [SPECIES], --species [SPECIES]
    Species name for AUGUSTUS training. No re-training if species already present in AUGUSTUS config folder
#### -d, --stranded        
    Run LoReAn on stranded mode [FALSE]
#### -f, --fungus
    Use this option for fungal species (used in Gene Mark-ES) [FALSE]
#### -k, --keep_tmp
    Keep temporary files [FALSE]
#### -sr [FASTQ_file], --short_reads [FASTQ_file]
    Path to short reads FASTQ. If paired end, comma- separated (1-1.fq,1-2.fq). BAM sorted files are
    allowed; the extension of the file should be filename.sorted.bam []
#### -lr [FASTQ_file], --long_reads [FASTQ_file]
    Path to long reads FASTQ []
#### -a [FASTA_file], --adapter [FASTA_file]
    FASTA file containing the adapter sequences. Adapter sequences in forward and reverse strain of the same
    adapter need to be used in the file []
#### -rp [GFF_file], --repeat_masked [GFF_file]
    GFF or GFF3 or GTF or BED file containing repeats
    coordinates []
#### -u [GFF_file], --update [GFF_file]
    GFF or GFF3 or GTF or BED file containing repeats
    coordinates []
#### -m [MAX_LONG_READ], --max_long_read [MAX_LONG_READ]
    Filter out long reads longer than this value (longer
    reads may affect mapping and assembling) [20000]
#### -pasa [PASA_DB], --pasa_db [PASA_DB]
    PASA database name [pipeline_run]
#### -n [PREFIX_GENE], --prefix_gene [PREFIX_GENE]
    Prefix to add to the final Gff3 gene name [specie]
#### -w [WORKING_DIR], --working_dir [WORKING_DIR]
    Working directory (will create if not present) [./]
#### -t, --threads
    Number of threads [3]
#### -cw N, --augustus_weigth N
    Weight assigned to AUGUSTUS evidence for EVM [1]
#### -uw N, --update_weigth N
    Weight assigned to GFF3 file to update for EVM [1]
#### -gw, --genemark_weigth
    Weight assigned to GENEMARK evidence for EVM [1]
#### -tw, --trinity_weigth
    Weight assigned to Trinity mapped with GMAP evidence for EVM [1]
#### -pw, --pasa_weigth
    Weight assigned to PASA evidence for EVM [5]
#### -aw, --AAT_weigth
    Weight assigned to AAT protein evidence for EVM [1]
#### -c, --segmentSize
    Segment size for EVM partitions [100000]
#### -e, --overlapSize
    Overlap size for EVM partitions [10000]
#### -g, --min_intron_length
    Minimal intron length for GMAP [9]
#### -q, --max_intron_length
    Maximal intron length for GMAP, STAR and TRINITY [1000]
#### -ee, --end_exon
    Minimal length for end exon with GMAP [20]
#### -cme, --cluster_min_evidence
    Minimal evidence needed to form a cluster [5]
#### -cMe, --cluster_max_evidence
    Maximal evidence to form a cluster.Prevents the
    clustering or rRNA genes i.e. [5000]
#### -aol, --assembly_overlapLength
    Minimal length (in nt) of overlap for ASSEMBLY [200]
#### -api, --assembly_percentIdentity
    Minimal identity for the ASSEMBLY (95-100) [97]
#### -art, --assembly_readThreshold
    Fraction of reads supporting an assembled UNITIG to keep on the ASSEMBLY (0.1-1) [0.3]
#### -v, --verbose
    Prints out the commands used in LoReAn[FALSE]
