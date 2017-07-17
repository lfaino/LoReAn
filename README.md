# LoReAn (Long-Read Annotation) for automated eukaryotic genome annotation incorporating long-reads

The LoReAn software is an automated annotation pipeline designed for eukaryotic genome annotation. It is built using previously defined annotation rationale and programs, but the key improvement is the incorporation of single-molecule cDNA sequencing data, such as that produced from [Oxford Nanopore](https://nanoporetech.com/) and from [PacBio](http://www.pacb.com/applications/rna-sequencing/). We find this significantly improves automated annotations and reduces the requirments for time-consuming manual annotation. 

We are working to improve LoReAn documentation and getting a paper into bioRxiv. For those familar with the annotation process and with docker, there should be enough infomation to run the program. Otherwise, check back soon ...

This is how LoReAn works: [LoReAn schematic view](https://github.com/lfaino/LoReAn/wiki)

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
The full list of options can be found at [option instructions](OPTIONS.md) or by:

```bash
lorean.py --help
```

LoReAn can run BRAKER1 to improve Augustus gene prediction;

To do so, short reads from RNA-seq or long reads RNA-seq need to be provided

## SOFTWARE USED IN THE PIPELINE

- samtools 0.1.19-96b5f2294a
- bedtools v2.25.0
- bowtie  version 1.1.2
- bamtools 2.4.1
- AATpackage r03052011 
- iAssembler v1.3.2.x64
- gm_et_linux_64 (THIS SOFTWARE IS NOT FREE FOR EVERYONE)
- PASApipeline v2.1.0 
- augustus v3.3
- trinityrnaseq v2.4.0
- STAR v2.5.2b
- gmap-gsnap v2017-06-20
- fasta v36.3.8e
- BRAKER v1.11
- EVidenceModeler v1.1.1
- gffread  v0.9.9
- genometools v1.5.9


## AUTHORS:
- Luigi Faino
- David Cook
- Jose Espejo


