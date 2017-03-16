# LoReAn (Long-Read Annotation) for automated eukaryotic genome annotation incorporating long-reads

The LoReAn software is an automated annotation pipeline designed for eukaryotic genome annotation. It is built using previously defined annotation rationale and programs, but the key improvement is the incorporation of single-molecule cDNA sequencing data, such as that produced from [Oxford Nanopore](https://nanoporetech.com/) and from [PacBio](http://www.pacb.com/applications/rna-sequencing/). We find this significantly improves automated annotations and reduces the requirments for time-consuming manual annotation. 

We are working to improve LoReAn documentation and getting a paper into bioRxiv. For those familar with the annotation process and with docker, there should be enough infomation to run the program. Otherwise, check back soon ...

## HOW TO RUN

LoReAn requires three mandatory files:
* Protein Sequences
* Reference genome 
* Genome name

To install the software:

Please see the [installation instructions](INSTALL.md) for details. 


The software can be run after installing by:
```bash
lorean.py <protein.fasta> <genome.fasta> <species name for augustus>
```
Additional options can be found at [option instructions](OPTIONS.md)

To see a full list of options type:
```bash
lorean.py --help
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
