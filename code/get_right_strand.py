#!/usr/bin/env python

###############
###IMPORTS###
###############

from sys import argv
import subprocess
import re
import copy
from BCBio import GFF

def gffread_multiexons(gff3Filename, multiExonFlag = False):
    '''Reports the multiexon genes in a GFF3 file - Output from GMAP'''
    if multiExonFlag:
        args = ['gffread', gff3Filename, '-o-', '-U']
    else:
        args = ['gffread', gff3Filename, '-o-']
    output = subprocess.check_output(args)
    outputList = output.splitlines()    
    return outputList

def gffread_parser(gffreadOutput):
    '''Parses the multiexon output to write a dict with names'''
    multiexon_dict = {}
    for line in gffreadOutput:
        if 'mRNA' in line:
            modline = re.split('\t|;', line)
            for element in modline:
                if element.startswith('ID='):
                    key = element.split('=')[-1]
                    if '_' in key:
                        spl_key = key.split('_')
                        key = spl_key[0]
                        el = spl_key[1]
                    else:
                        el = ''
                    if '.mrna' in el:
                        el = el.split('.mrna')[0]
                    multiexon_dict[key] = el
    return multiexon_dict
                    
def compare_dicts(multiExon, gmapAll, pasaAll):
    '''first compares the pasa dict with the multiexon, then the ones from gmap
    to take all the single exon genes'''
    pasaSingle = {}    
    
    for key in pasaAll.keys():        
        if key not in multiExon.values():            
            pasaSingle[key] = ''
    gmapSingle = {}        
    for key, value in gmapAll.items():
        if not multiExon.has_key(key) and not pasaSingle.has_key(value):
            
            gmapSingle[key] = value
    gmapOut = dict(multiExon.items() + gmapSingle.items())
    return gmapOut, pasaSingle


def combineGff3(gmapDict, pasaDict, gffreadGMAP, gffreadPASA, wd):
    '''the function check for the rigth strand between the GMAP gff3 and the PASA gff3'''
    outputFilename = wd+'finalAnnotation.gff3'
    outputFile = open(outputFilename, 'w')
    for line in gffreadGMAP:   
        line = line.strip()
        if line.startswith('#'):
            continue
        modline = re.split('\t|;', line)  
        if modline[2] == 'mRNA':
            for element in modline:
                if element.startswith('ID='):
                    cand_key = element.split('=')[-1]
                    cand_key = cand_key.split('.mrna')[0] 
                    cand_key = cand_key.split('_')[0]
                    if gmapDict.has_key(cand_key):
                        gene_line = copy.copy(line)
                        gene_line = gene_line.split('\t')
                        gene_line[2] = 'gene'
                        attr = gene_line[-1].split(';')
                        attr_changed = ''
                        for element in attr:
                            if element.startswith('geneID='):
                                attr_changed += 'ID=' + element.split('=')[-1] + ';'
                            elif element.startswith('gene_name='):
                                attr_changed += 'Name=' + element.split('=')[-1]                     
                        gene_line[-1] = attr_changed
                        gene_line = '\t'.join(gene_line)
                        outputFile.write(gene_line+'\n') 
                        mrna_line = copy.copy(line)
                        mrna_line = mrna_line.split('\t')
                        attr = mrna_line[-1].split(';')
                        attr_changed = ''                        
                        for element in attr:
                            if element.startswith('geneID='):
                                attr_changed += 'Parent=' + element.split('=')[-1]
                            elif element.startswith('ID='):
                                attr_changed += element  + ';'                   
                        mrna_line[-1] = attr_changed
                        mrna_line = '\t'.join(mrna_line)
                        outputFile.write(mrna_line+'\n')    
        elif modline[2] == 'CDS' or modline[2] == 'exon':            
            for element in modline:
                if element.startswith('Parent='):
                    cand_key = element.split('=')[-1]
                    cand_key = cand_key.split('.mrna')[0]
                    cand_key = cand_key.split('_')[0]
                    if gmapDict.has_key(cand_key):
                        outputFile.write(line+'\n') 
    for line in gffreadPASA:   
        line = line.strip()
        if line.startswith('#'):
            continue
        modline = re.split('\t|;', line)  
        if modline[2] == 'mRNA':
            
            for element in modline:
                if element.startswith('ID='):
                    
                    cand_key = element.split('=')[-1]
                    cand_key = cand_key.split('.mrna')[0]
                    if pasaDict.has_key(cand_key):  
                        gene_line = copy.copy(line)
                        gene_line = gene_line.split('\t')
                        gene_line[2] = 'gene'
                        attr = gene_line[-1].split(';')
                        attr_changed = ''
                        for element in attr:
                            if element.startswith('geneID='):
                                attr_changed += 'ID=' + element.split('=')[-1] + ';'
                            elif element.startswith('gene_name='):
                                attr_changed += 'Name=' + element.split('=')[-1]                     
                        gene_line[-1] = attr_changed
                        gene_line = '\t'.join(gene_line)
                        outputFile.write(gene_line+'\n') 
                        mrna_line = copy.copy(line)
                        mrna_line = mrna_line.split('\t')
                        attr = mrna_line[-1].split(';')
                        attr_changed = ''                        
                        for element in attr:
                            if element.startswith('geneID='):
                                attr_changed += 'Parent=' + element.split('=')[-1] 
                            elif element.startswith('ID='):
                                attr_changed += element + ';'                   
                        mrna_line[-1] = attr_changed
                        mrna_line = '\t'.join(mrna_line)
                        outputFile.write(mrna_line+'\n')
                               
        elif modline[2] == 'CDS' or modline[2] == 'exon':            
            for element in modline:
                
                if element.startswith('Parent='):
                    cand_key = element.split('=')[-1]
                    cand_key = cand_key.split('.mrna')[0]                    
                    if pasaDict.has_key(cand_key):
                        outputFile.write(line+'\n')            
          
   
    outputFile.close()
    outputFileSort = wd+'finalAnnotation.sort.gff3'
    outputFileLog = wd+'finalAnnotation.sort.gff3.log'
    ###here we use gt software ro reduce the ID length in the gff for each gene. It is important because with long names and ID, PASA will complain
    gt_con = ['gt', 'gff3', '-sort', '-tidy', outputFilename ]
    gt_call = subprocess.Popen(gt_con, stdout = file(outputFileSort , "w"), stderr = file (outputFileLog, "w"))
    gt_call.communicate()
    return outputFileSort

            
            
            
if __name__ == '__main__':
    gff3Pasa = argv[1]
    gff3GMAP = argv[2]
    wd = './'
    outputList_gmap_multi = gffread_multiexons(gff3GMAP, multiExonFlag = True)
    outputList_gmap_all = gffread_multiexons(gff3GMAP)
    outputList_pasa_all = gffread_multiexons(gff3Pasa)
    gmap_dict_multi = gffread_parser(outputList_gmap_multi)
    gmap_dict_all = gffread_parser(outputList_gmap_all)
    pasa_dict_all = gffread_parser(outputList_pasa_all)
        
    gmapOut, pasaOut = compare_dicts(gmap_dict_multi, gmap_dict_all, pasa_dict_all)
        
    combineGff3(gmapOut, pasaOut, outputList_gmap_all, outputList_pasa_all, wd)
        
