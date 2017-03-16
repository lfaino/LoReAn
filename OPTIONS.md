# Minimal options are:
```bash
lorean.py <protein.fasta> <genome.fasta> <species name for augustus>
```

## The specie name used in LoReAn needs to be in the AUGUSTUS database. However, if RNA-seq is supplied,the tool will calculate the best AUGUSTUS sets using Braker1 software

-d or --stranded: 
    
    Run LoReAn on stranded mode [FALSE]
    
-f or --fungus: 
    
    Use this option for fungal species (used in Gene Mark-ES)  [FALSE]
    
-u or --only_unitigs: 
    
    Removes gene models that are not supported by long reads [FALSE]
    
-k or --keep_tmp: 
    
    Keep temporary files [FALSE]
    
-s or --short_reads: 
    Path to short reads FASTQ. If paired end, comma-separated (1-1.fq,1-2.fq). BAM sorted files are allowed; the extension of the file should be filename.sorted.bam []
    
-a or --adapter:    
    FASTA file containing the adapter sequences. Adapter sequences in forward and reverse strain of the same adapter need to be used in the file []
    
-l or --long_reads: 
    Path to long reads FASTQ []
    
-m or --max_long_read: 
    Filter out long reads longer than this value (longer reads may affect mapping and assembling) [20000]
    
-p or --pasa_db: 
    PASA database name [pipeline_run]
    
-n or --prefix_gene: 
    Prefix to add to the final Gff3 gene name [specie]
    
-w or --working_dir: 
    Working directory (will create if not present) [./annotation]
    
-t or --threads: 
    Number of threads [1]
    
-b or --overhang: 
    CAP3 max overhang percent length; this value should be > 3 [20]
    
-cw or --augustus: 
    Weight assigned to AUGUSTUS evidence for EVM [1]
    
-gw or --genemark: 
    Weight assigned to GENEMARK evidence for EVM [1]
    
-tw or --trinity: 
    Weight assigned to Trinity mapped with GMAP evidence for EVM [1]
    
-pw or --pasa: 
    Weight assigned to PASA evidence for EVM [5]
    
-aw or --AAT:
    Weight assigned to AAT protein evidence for EVM [1]
    
-c or --segmentSize:
    Segment size for EVM partitions [100000]
    
-e or --overlapSize:
    Overlap size for EVM partitions [10000]
    
-g or --min_intron_length:
    Minimal intron length for GMAP [9]
    
-q or --max_intron_length:
    Maximal intron length for GMAP, STAR and TRINITY [1000]
    
-ee or --end_exon:
    Minimal length for end exon with GMAP [20]
    
-cme or --cluster_min_evidence: 
    Minimal evidence needed to form a cluster [5]
    
-cMe or --cluster_max_evidence: 
    Maximal evidence to form a cluster. Prevents the clustering or rRNA genes i.e. [5000]
    
-aol or --assembly_overlapLength: 
    Minimal length (in nt) of overlap for ASSEMBLY [200]
    
-api or --assembly_percentIdentity: 
    Minimal identity for the ASSEMBLY (95-100) [97]
    
-art or --assembly_readThreshold: 
    Fraction of reads supporting an assembled UNITIG to keep on the ASSEMBLY (0.1-1) [0.3]
    
-ne or --no_EVM: 
    Run until the preparation of EVM inputs [FALSE]
    
-nc or --no_consensus: 
    Do not run the long reads consensus pipeline [FALSE]
    
-nu or --no_update: 
    Do not run the PASA update[FALSE]
    
-co or --collect_only: 
    Collect only assebmled transcripts [FALSE]

