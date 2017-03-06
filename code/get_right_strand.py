#!/usr/bin/env python

###############
###IMPORTS###
###############

import gffutils
from sys import argv
import re
import subprocess
import gffutils.gffwriter as gffwriter



def genename(gff_filename, prefix):

    gt_call = subprocess.Popen(['gt', 'gff3', '-sort', '-tidy', gff_filename], stdout=subprocess.PIPE, stderr = file(gff_filename + '.gt.log', "w"))
    fields = []
    featureann = []
    chrold = ''
    count = 0
    path = []
    path = gff_filename.split('/')
    path[-1] = prefix + '_LoReAn.annotation.gff3'
    newpath = '/'.join(path)

    o = open(newpath, 'w+')
    o.write ('##version-gff 3\n')
    for ln in gt_call.stdout.readlines():
        lnn = ln.rstrip('\n')
        fields = lnn.split('\t')
        if len(fields) == 9:
            typef = fields[2]
            fields[1] = 'LoReAn'
            if 'gene' in typef:
                featureann = fields[8].split(';')
                for elm in featureann:
                    if 'ID' in elm:
                        idname = elm
                    elif 'Name' in elm:
                        if fields[0] == chrold:
                            count += 10
                            fatureattgene = prefix + '_' + chrold + '_G_' + str(count)
                            fields[8] = idname + ';Name=' + fatureattgene
                        else:
                            chrold = fields[0]
                            count = 0
                            count += 10
                            fatureattgene = prefix + '_' + chrold + '_G_' + str(count)
                            fields[8] = idname + ';Name=' + fatureattgene
                o.write ('\t'.join(fields)+ '\n')
                
            elif 'mRNA' in typef:
                featureann = fields[8].split(';')
                for elm in featureann:
                    if 'ID' in elm:
                        idname = elm
                    elif 'Parent' in elm:
                        parentname = elm
                fields[8] = idname + ';' + parentname
                o.write ('\t'.join(fields)+ '\n')

                    
            else:
                o.write ('\t'.join(fields) + '\n')
    o.close()
    return newpath

def newNames(oldname):

    finalout = oldname + ".sorted.gff"
    errorFile = oldname + ".gt_err.log"
    gt_com = ['gt', 'gff3', '-sort', '-tidy', oldname]
    gt_call = subprocess.Popen(
        gt_com, stdout=file(
            finalout, 'w'), stderr=file(
            errorFile, 'w'))
    gt_call.communicate()

    return finalout

def strand(gff_file1, gff_file2, wd):
    outputFilename = wd + 'finalAnnotation.gff3'
    gff_out = gffwriter.GFFWriter(outputFilename)

    gff_file1_out = gff_file1 + ".intron.tidy.sorted.gff"
    errorFile = gff_file1 + ".gt_err.log"
    gt_com = ['gt', 'gff3', '-sort', '-tidy', '-addintrons', '-retainids', gff_file1]
    gt_call = subprocess.Popen(
        gt_com, stdout=file(
            gff_file1_out, 'w'), stderr=file(
            errorFile, 'w'))
    gt_call.communicate()

    gff_file2_out = gff_file2 + ".intron.tidy.sorted.gff"
    errorFile = gff_file2 + ".gt_err.log"
    gt_com = ['gt', 'gff3', '-sort', '-tidy', '-addintrons', '-retainids', gff_file2]
    gt_call = subprocess.Popen(
        gt_com, stdout=file(
            gff_file2_out, 'w'), stderr=file(
            errorFile, 'w'))
    gt_call.communicate()

    db1 = gffutils.create_db(gff_file1_out, ':memory:',  merge_strategy='create_unique', keep_order=True)
    db2 = gffutils.create_db(gff_file2_out, ':memory:',  merge_strategy='create_unique', keep_order=True)

    listgene1 = []
    listgeneintrons = []
    listgenetotal = []

    for i in db1.features_of_type("intron"):
        g = ' '.join (i.attributes['Parent'])
        listgeneintrons.append(g)
        
    for i in db1.features_of_type("CDS"):
        g = ' '.join (i.attributes['Parent'])
        listgenetotal.append(g)
        
    listgene1 = sorted(set(list(set(listgenetotal)^set(listgeneintrons))))

    listgene2 = []
    listgeneintrons = []
    listgenetotal = []
    newlist = []
    
    for i in db2.features_of_type("intron"):
        g = ' '.join (i.attributes['Parent'])
        listgeneintrons.append(g)
        
    for i in db2.features_of_type("CDS"):
        g = ' '.join (i.attributes['Parent'])
        listgenetotal.append(g)
    
    listgene2 = sorted(set(list(set(listgenetotal)^set(listgeneintrons))))
    
    geneDict = {}
    for a in listgene2:   
        b =  re.split('_|\.', a)
        del b[-1]
        del b[0]
        evm = '.'.join(b)
        newlist.append(evm)
        geneDict[evm]=a

    commonlist = list(set(listgene1).intersection(newlist))
    
    for evm in commonlist:
        if geneDict[evm]:
            #print geneDict[evm]
            del geneDict[evm]

    listuniq = [] 
    listUniqIntrons = []
    listgmap = []

    for a in geneDict:
        listuniq.append(geneDict[a])
    
    listUniqIntrons = list(set(listgeneintrons))
    listgmap = list(set(listuniq))

    for evm in commonlist:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        #for i in db1.parents(evm, featuretype='mRNA', order_by='start'):
            #gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
       
    for evm in listUniqIntrons:
        for i in db2.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i) 
        gff_out.write_rec(db2[evm])
        for i in db2.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        #for i in db2.parents(evm, featuretype='mRNA', order_by='start'):
            #gff_out.write_rec(i)
        for i in db2.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    
    for evm in listgmap:
        for i in db2.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i) 
        gff_out.write_rec(db2[evm])
        for i in db2.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        #for i in db2.parents(evm, featuretype='mRNA', order_by='start'):
            #gff_out.write_rec(i)
        for i in db2.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    return outputFilename

if __name__ == '__main__':
    
    
    gff1 = argv[1]
    gff2 = argv[2]
    wd = argv[3]
    prefix = argv[4]
    out1 = strand(gff1, gff2, wd)
    out2 = newNames(out1)
    genename(out2, prefix)
