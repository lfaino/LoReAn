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


def combineGff3(gmapDict, pasaDict, gffreadGMAP, gffreadPASA, ref, wd):
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
    return outputFilename


def changeName(outputFilename, ref):

    
    outputFileSort = outputFilename+'.changed.Name.sort.gff3'
    outputFileLog = outputFilename+'.changed.Name.sort.gff3.log'
    ###here we use gt software ro reduce the ID length in the gff for each gene. It is important because with long names and ID, PASA will complain
    gt_con = ['gt', 'gff3', '-sort', '-tidy', outputFilename ]
    gt_call = subprocess.Popen(gt_con, stdout = file(outputFileSort , "w"), stderr = file (outputFileLog, "w"))
    gt_call.communicate()
    countGene = ""
    old = ""
    outputSimply = outputFileSort + '.newGeneName' 
    f = open(outputFileSort, 'r')
    o = open(outputSimply, 'w')
    for line in f:
        testLine = line.rstrip()
        fields = testLine.split(";")
        if "merged" in outputFilename:
            newLine = fields[0]

            chrNub = fields[0].split("\t")
            num = chrNub[0]
            if len(chrNub) > 8:
                if "locus" in fields[0]:
                    next
                elif "mRNA" in chrNub[2]:
                    if old == num:
                        #print "inside"
                        countGene += 1
                        chrNub[2] = "gene"
                        chrNub.pop()
                        tab = '\t'
                        name_gene = "Name=" + ref + "_" + old + "g" + str(countGene) 
                        newGeneLine = tab.join(chrNub) + "\tID=" + ref + "_" + old + "g" + str(countGene) + ";" + name_gene + "\n"
                        newName = fields[0] + ";Parent=" + ref + "_" + old + "g" + str(countGene) + ";" + name_gene + "\n"

                        o.write(newGeneLine)
                        o.write(newName)
                    else:
                        countGene = 1
                        old = num
                        chrNub[2] = "gene"
                        chrNub.pop()
                        tab = '\t'
                        name_gene = "Name=" + ref + "_" + old + "g" + str(countGene)
                        newGeneLine = tab.join(chrNub) + "\tID=" + ref + "_" + old + "g" + str(countGene) + ";" + name_gene + "\n"
                        newName = fields[0] + ";Parent=" + ref + "_" + old + "g" + str(countGene) + ";" + name_gene + "\n"
                        o.write(newGeneLine)
                        o.write(newName)
                else: 
                    o.write(newLine + ";" + name_gene + "\n")
                    

        elif "final" in outputFilename:
            if "gene" in fields[0]:
                chrNub = fields[0].split("\t")
                num = chrNub[0]
                if old == num:
                    #print "inside"
                    countGene += 1
                    newName = fields[0] + ";Name=" + ref + "_" + old + "g" + str(countGene) + "\n"
                    o.write(newName)
                else:
                    countGene = 1
                    old = num
                    newName = fields[0] + ";Name=" + ref + "_" + old + "g" + str(countGene) + "\n"
                    o.write(newName)
            else:
                o.write(line)
            
    o.close()
    f.close()
    
    simplifiedGff = cluster_gff(outputSimply)
    return simplifiedGff


def cluster_gff(outputSimply): 
    '''here the cluseter of sequence from the same locus are prepared'''    


    f = open(outputSimply, 'r')
    listGenes = []
    for line in f:
        geneLine = line.rstrip()
        listGenes.append(line)
        
    f.close()


    BTsort1 = ['bedtools', 'sort', '-i', outputSimply]
    BTgrep = ['grep' , 'gene']
    BTmerge1 = ['bedtools', 'merge' , '-o', 'distinct,count' , '-c' , '9,9']

    BTsort1_call = subprocess.Popen(BTsort1, stdout = subprocess.PIPE )
    Grep_call = subprocess.Popen(BTgrep , stdout = subprocess.PIPE, stdin = BTsort1_call.stdout)
    BTmerge1_call = subprocess.Popen(BTmerge1,stdin = Grep_call.stdout , stdout = subprocess.PIPE)
    out = BTmerge1_call.stdout.readlines()

    simplifiedGff = outputSimply + 'simply.gff3'
    o = open(simplifiedGff, 'w')

    for line in out:
        bedLine = line.split('\t')
        count = int(bedLine[4])
        if count < 30:
            fieldAll = re.split(';|,',bedLine[3])
            for element in fieldAll:
                if 'Name=' in element:
                    #print element
                    for lineGene in listGenes:
                        if element in lineGene:
                            bedLine = lineGene.split(';')
                            if "gene" in bedLine[0]:
                                geneLine = lineGene.rstrip()
                                listGenes.remove(lineGene)
                                o.write(lineGene)
                            else:
                                bedLine.pop(-1)
                                newBed = ';'.join(bedLine)
                                listGenes.remove(lineGene)
                                o.write(newBed)
    o.close()
    return simplifiedGff

            
            
            
if __name__ == '__main__':
    change = argv[1]
    ref = "test"
    test = changeName(change, ref)
    
