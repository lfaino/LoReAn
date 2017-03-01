# LoReAn

##HOT TO RUN

LoReAn requires three mandatory files:
* Protein Sequences
* Reference genome 
* Genome name

The software can be run after installatio by:
```bash
lorean.py <protein.fasta> <genome.fasta> <species name for augustus>
```

LoReAn can run BRAKER1 to improve Augustus gene prediction;

To do so, short reads from RNA-seq or long reads RNA-seq need to be provided

##SOFTWARE USED IN THE PIPELINE

- samtools
- bedtools
- bowtie
- bamtools
- AATpackage-r03052011 
- iAssembler-v1.3.2.x64
- gm_et_linux_64 (THIS SOFTWARE IS NOT FREE FOR EVERYONE)
- PASApipeline-2.0.2 
- augustus-3.2.2
- trinityrnaseq-2.4.0
- STAR-2.5.2b
- gmap-gsnap-2015-12-31.v10 
- fasta-36.3.8e
- BRAKER1
- EVidenceModeler1.1.1
- gffread
- genometools-1.5.9

##AUTHORS:
- Luigi Faino
- David Cook 
- Jose A. Espejo
