#!/usr/bin/env python

###############
###IMPORTS###
###############

import gffutils
from sys import argv
import re
import subprocess
import gffutils.gffwriter as gffwriter
from pybedtools import BedTool

def removeDiscrepancy(gff, evmFile):
    badName = []
    gffVal_call = subprocess.Popen(['pasa_gff3_validator.pl', gff], stdout=subprocess.PIPE)
    for ln in gffVal_call.stdout.readlines():
        name = re.split(' |\.CDS', ln)
        if len(name) > 3 and "ERROR" in name[0]:
            #if "mRNA" in name[2]:
            badName.append(name[2])
    badNameUniq = list(set(badName))
    i = open(gff, 'r')
    listAllName = []
    for line in i:
        fields = line.strip().split('\t')
        if len(fields) > 3:
            if "mRNA" in fields[2]:
                attribute = fields[8].split(';')
                for el in attribute:
                    if "ID" in el:
                        listAllName.append(el.split("=")[1])
    listgene = sorted(set(list(set(listAllName)^set(badNameUniq))))
    outputFilename = gff + '.noProblem.gff3'
    gff_out = gffwriter.GFFWriter(outputFilename)
    db1 = gffutils.create_db(gff, ':memory:',  merge_strategy='create_unique', keep_order=True)
    for evm in listgene:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    bedtools_call = subprocess.Popen(['bedtools', 'intersect', '-v', '-a', evmFile, '-b', gff], stdout=subprocess.PIPE)
    evm_mRNA = []
    for ln in bedtools_call.stdout.readlines():
        if "mRNA" in ln.split('\t')[2]:
            attribute = ln.split('\t')[8].split(';')
            for el in attribute:
                if "ID" in el:
                    mRNA = el.split('=')[1]
                    evm_mRNA.append(mRNA)
    db1 = gffutils.create_db(evmFile, ':memory:',  merge_strategy='create_unique', keep_order=True)
    for evm in evm_mRNA:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
            
    return(outputFilename)

def appendID(gff):
    i = open(gff, 'r')
    outFile = gff + '.reformat.gff'
    o = open (outFile, 'w')
    for line in i:
        chompLine = line.strip()
        listLine = chompLine.split('\t')
        if len(listLine) == 9:
            if "exon" in listLine[2]:
                attribute = listLine[8].split(';')
                for a in attribute:
                    if "Parent" in a:
                        parent = a.split("=")
                        chrName = listLine[0].lower().replace("[a-z]", "")
                        newAtt = a + ";ID=" + parent[1] + ".exon." + str(listLine[3]) + str(listLine[4]) + chrName
                        listLine[8] = newAtt
            elif  "CDS" in listLine[2]:
                attribute = listLine[8].split(';')
                for a in attribute:
                    if "Parent" in a:
                        parent = a.split("=")
                        chrName = listLine[0].lower().replace("[a-z]", "")
                        newAtt = a + ";ID=" + parent[1] + ".CDS." + str(listLine[3]) + str(listLine[4]) + chrName
                        listLine[8] = newAtt
            elif "mRNA" in listLine[2]:
                attribute = listLine[8].split(';')
                attri = []
                for a in attribute:
                    if "Parent" in a or "ID" in a:
                        attri.append(a)
                listLine[8] = ';'.join(attri)
        o.write('\t'.join(listLine) + '\n')
    o.close()
    i.close()
    return(outFile)

def removeOverlap(gff):
    i = open(gff, 'r')
    outFile = gff + '.RNA.gff'
    o = open (outFile, 'w')
    for line in i:
        listLine = line.split('\t')
        if len(listLine) == 9:
            if "mRNA" in listLine[2]:
                o.write(line)
    o.close()
    i.close()
    bedout = outFile + '.bedtools.bed'
    errorFile = outFile + ".bedtools_err.log"
    bedsort = ['bedtools', 'sort', '-i', outFile]
    bedmerge = ['bedtools', 'merge', '-d', '-100', '-s', '-o', 'count,distinct','-c','9,9' ]
    bedsort_call = subprocess.Popen(bedsort, stdout=subprocess.PIPE, stderr=file(errorFile, 'w'))
    bedmerge_call = subprocess.Popen(bedmerge, stdin = bedsort_call.stdout, stdout=file(bedout , 'w'), stderr=file(errorFile, 'w'))
    bedmerge_call.communicate()
    listMultiple = []
    listUniq= []
    count = 0
    dictRNA = {}
    i = open(bedout, 'r')
    for a in i:
        listLine = a.split('\t')
        nameRNA = re.split(',|;', listLine[5])
        count += 1
        locus = "locus" + str(count)
        for elm  in nameRNA:
            if "ID" in elm and int(listLine[4]) > 1:
                mRNAname = elm.split('=')[1]
                listMultiple.append(mRNAname)
                if dictRNA.has_key(locus):
                    dictRNA[locus].append(mRNAname)
                else:
                    dictRNA[locus] = [mRNAname,]
            elif "ID" in elm:
                mRNAname = elm.split('=')[1]
                listUniq.append(mRNAname)
    listMultipleUniq = []
    listMultipleUniq = list(set(listMultiple))
    dictLength= {}
    mRNA = open(gff, 'r')
    for line in mRNA:
        listLine = line.split('\t')
        if len(listLine) == 9:
            if "mRNA" in listLine[2]:
                for key in dictRNA:
                    for el in dictRNA[key]:
                        nameID = "ID=" + el  +  ';'
                        if nameID in line:
                            length = (int(line.split('\t')[4]) - int(line.split('\t')[3]))
                            if dictLength.has_key(key):
                                oldLenght = dictLength[key]
                                if int(oldLenght[1]) < int(length):
                                    dictLength[key] = [el, str(length)]
                            else:
                                dictLength[key] = [el, str(length)]
    for key in dictLength:
        listUniq.append(dictLength[key][0])
    listUniqNew = []
    for mRNA in listUniq:
        mRNAnew = mRNA.strip('\n')
        if mRNAnew not in listUniqNew:
            listUniqNew.append(mRNAnew)
    outputFilename = gff + '.uniq.gff3'
    gff_out = gffwriter.GFFWriter(outputFilename)
    db1 = gffutils.create_db(gff, ':memory:',  merge_strategy='create_unique', keep_order=True)
    for evm in listUniqNew:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)


    return outputFilename

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
    gt_com = ['gt', 'gff3', '-retainids', '-sort', '-tidy', oldname]
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
        #print b
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
    
    
    evmgff = argv[1]
    gmapgff = argv[2]
    wd = argv[3]
    prefix = argv[4]
    out1 = strand(evmgff, gmapgff, wd)
    out2 = appendID(out1)
    out3 = removeOverlap(out2)
    out4 = removeDiscrepancy(out3, evmgff)
    newNames(out4)
    #removeOverlap(out2)
    #genename(out2, prefix)
