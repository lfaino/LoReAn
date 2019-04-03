#!/usr/bin/env python3

import subprocess

def export_path():
    path_list=['export AUGUSTUS_CONFIG_PATH=~/augustus/config/', 'export AUGUSTUS_BIN_PATH=/opt/LoReAn/third_party/software/augustus/bin',
               'export PATH=$PATH:$AUGUSTUS_BIN_PATH:/opt/LoReAn/third_party/software/LoReAn/:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1:'
               '/opt/LoReAn/third_party/software/gm_et_linux_64/gmes_petap:/opt/LoReAn/third_party/software/iAssembler-v1.3.2.x64:'
               '/opt/LoReAn/third_party/software/PASApipeline/scripts:/opt/LoReAn/third_party/software/PASApipeline/:'
               '/opt/LoReAn/third_party/software/Trinity:/opt/LoReAn/third_party/software/AATpackage-r03052011/bin/:'
               '/opt/LoReAn/third_party/software/LoReAn/third_party/:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/:'
               '/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/EvmUtils/misc/:/opt/LoReAn/third_party/scripts/:'
               '/opt/LoReAn/code/:/opt/LoReAn/third_party/software/STAR/bin/Linux_x86_64/:/opt/LoReAn/third_party/software/PASApipeline/bin/:'
               '/opt/LoReAn/third_party/software/genometools-1.5.9/bin/:/opt/LoReAn/third_party/software/BRAKER/scripts/:'
               '/opt/LoReAn/third_party/software/interproscan-5.27-66.0/:/opt/LoReAn/third_party/software/TransDecoder-3.0.1/:'
               '/opt/LoReAn/third_party/software/TransDecoder-3.0.1/util:/opt/LoReAn/third_party/conf_files:/opt/:/usr/local/RepeatMasker:/usr/local/RepeatScout:'
               '/opt/LoReAn/third_party/software/minimap2/:/opt/LoReAn/third_party/software/SE-MEI/:/opt/LoReAn/third_party/software/PASApipeline/misc_utilities/',
               'export GENEMARK_PATH=/opt/LoReAn/third_party/software/gm_et_linux_64/gmes_petap/', 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/SW/src/',
               'export PERL5LIB=$PERL5LIB:/opt/LoReAn/third_party/scripts/:/opt/LoReAn/third_party/software/EVidenceModeler-1.1.1/PerlLib/',
               'export SAMTOOLS_PATH=/usr/local/bin/', 'export BLAST_PATH=/usr/local/bin/']
    for library in path_list:
        com = library
        create_user_call = subprocess.Popen(com, shell=True)
        create_user_call.communicate()