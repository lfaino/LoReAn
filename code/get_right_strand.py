#!/usr/bin/env python3

import re
import os
import gffutils
from sys import argv
import subprocess
import gffutils.gffwriter as gffwriter
from pybedtools import BedTool
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from multiprocessing import Pool

def removeDiscrepancy(gff, evmFile):
    badName = []
    gffVal_call = subprocess.Popen(['pasa_gff3_validator.pl', gff], stdout=subprocess.PIPE)
    for ln in gffVal_call.stdout.readlines():
        name = re.split(' |\.CDS', (ln.decode("utf-8")))
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
        lne = ln.decode("utf-8")
        ln = lne
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
            if "CDS" in listLine[2]:
                o.write(line)
    o.close()
    i.close()
    bedout = outFile + '.bedtools.bed'
    bedouffile = open(bedout, "w")
    errorFile = outFile + ".bedtools_err.log"
    errorFilefile = open(errorFile, "w")
    bedsort = ['bedtools', 'sort', '-i', outFile]
    bedmerge = ['bedtools', 'merge', '-d', '-100', '-s', '-o', 'count,distinct','-c','9,9' ]
    bedsort_call = subprocess.Popen(bedsort, stdout=subprocess.PIPE, stderr= errorFilefile)
    bedmerge_call = subprocess.Popen(bedmerge, stdin = bedsort_call.stdout, stdout=bedouffile, stderr=errorFilefile)
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
            if "Parent" in elm and int(listLine[4]) > 1:
                mRNAname = elm.split('=')[1]
                listMultiple.append(mRNAname)
                if locus in dictRNA:
                    dictRNA[locus].append(mRNAname)
                else:
                    dictRNA[locus] = [mRNAname,]
            elif "Parent" in elm:
                mRNAname = elm.split('=')[1]
                listUniq.append(mRNAname)
    listMultipleUniq = []
    listMultipleUniq = list(set(listMultiple))
    dictLength= {}
    mRNA = open(gff, 'r')
    for line in mRNA:
        listLine = line.split('\t')
        if len(listLine) == 9:
            if "CDS" in listLine[2]:
                for key in dictRNA:
                    for el in dictRNA[key]:
                        nameID = "Parent=" + el  +  ';'
                        if nameID in line:
                            length = (int(line.split('\t')[4]) - int(line.split('\t')[3]))
                            if key in dictLength:
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
    errorFile = gff_filename + "error.log"
    errorFilefile = open(errorFile, "w") 
    gt_call = subprocess.Popen(['gt', 'gff3', '-sort', '-tidy', gff_filename], stdout=subprocess.PIPE, stderr = errorFilefile)
    errorFilefile.close()
    
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
        lnn = (ln.decode("utf-8")).rstrip('\n')
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
    finaloutfile = open(finalout, "w")
    errorFilefile = open (errorFile, "w")
    gt_com = ['gt', 'gff3', '-sort', '-tidy', oldname]
    gt_call = subprocess.Popen(gt_com, stdout= finaloutfile, stderr=errorFilefile)
    gt_call.communicate()
    finaloutfile.close()
    errorFilefile.close()
    return finalout

def longest(gff_file, fasta, proc, wd): 
    outputFilename = wd + 'finalAnnotation.strand.gff3'
    outputFilenameLeft = wd + 'finalAnnotation.Left.gff3'
    gff_out = gffwriter.GFFWriter(outputFilenameLeft)
    gff_file_out = gff_file + ".intron.tidy.sorted.gff"
    errorFile = gff_file + ".gt_err.log"
    
    gt_com = ['gt', 'gff3', '-sort', '-tidy', '-addintrons', gff_file]
    gff_file_outfile = open(gff_file_out, "w")
    errorFilefile = open(errorFile, "w")
    gt_call = subprocess.Popen( gt_com, stdout= gff_file_outfile , stderr=errorFilefile)
    gt_call.communicate()
    gff_file_outfile.close()
    errorFilefile.close()
    
    gtf_file_out = gff_file + ".intron.tidy.sorted.gtf"
    errorFile = gff_file + ".gt_err.log"

    gt_com = ['gt', 'gff3_to_gtf', gff_file_out]
    gtf_file_outfile = open(gtf_file_out, "w")
    errorFilefile = open(errorFile, "w")
    gt_call = subprocess.Popen( gt_com, stdout= gtf_file_outfile , stderr=errorFilefile)
    gt_call.communicate()
    gtf_file_outfile.close()
    errorFilefile.close()
    
    db1 = gffutils.create_db(gff_file_out, ':memory:',  merge_strategy='create_unique', keep_order=True)

    fasta_file_out = gff_file + ".intron.tidy.sorted.fasta"
    errorFile = gff_file + ".1.gt_err.log"
    fasta_file_outfile = open(fasta_file_out, "w")
    errorFilefile = open(errorFile, "w")
    com = ['/opt/LoReAn/third_party/software/TransDecoder-3.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl', gtf_file_out, fasta]
    call = subprocess.Popen(com, stdout= fasta_file_outfile , stderr=errorFilefile)
    call.communicate()
    fasta_file_outfile.close()
    errorFilefile.close()
    
    gff_file_out_u = gtf_file_out + ".gff3"
    errorFile = gtf_file_out + ".2.gt_err.log"
    gff_file_outfile = open(gff_file_out_u, "w")
    errorFilefile = open(errorFile, "w")    
    com = ['/opt/LoReAn/third_party/software/TransDecoder-3.0.1/util/cufflinks_gtf_to_alignment_gff3.pl', gtf_file_out]
    call = subprocess.Popen(com, stdout= gff_file_outfile , stderr=errorFilefile)
    call.communicate()
    gff_file_outfile.close()
    errorFilefile.close()
    #/opt/LoReAn/third_party/software/TransDecoder-3.0.1/util/
    
    errorFile = gtf_file_out + ".TrDec_err.2.log"
    gff_file_out = gtf_file_out + ".TrDec_err.stdout"
    gff_file_outfile = open(gff_file_out, "w")
    errorFilefile = open(errorFile, "w")    
    com = ['TransDecoder.LongOrfs', '-m', '10', '-t', fasta_file_out ]
    call = subprocess.Popen(com, stdout = gff_file_outfile, stderr=errorFilefile, cwd = wd)
    call.communicate()
    errorFilefile.close()
    gff_file_outfile.close()
    
    
    errorFile = gtf_file_out + ".TrDec_err.3.log"
    gff_file_out = gtf_file_out + ".TrDec_err.stdout"
    gff_file_outfile = open(gff_file_out, "w")
    errorFilefile = open(errorFile, "w")
    wd_fasta =  fasta_file_out
    com = ['TransDecoder.Predict', '--single_best_orf','--cpu', str(proc), '--retain_long_orfs','10', '-t', wd_fasta]
    call = subprocess.Popen(com, stdout = gff_file_outfile, stderr=errorFilefile, cwd = wd)
    call.communicate()
    errorFilefile.close()
    gff_file_outfile.close()

    errorFile = gtf_file_out + ".TrDec_err.4.log"
    gff_file_out = gtf_file_out + ".TrDec_err.stdout"
    gff_file_outfile = open(outputFilename, "w")
    errorFilefile = open(errorFile, "w")
    wd_fasta = fasta_file_out    
    com = ['/opt/LoReAn/third_party/software/TransDecoder-3.0.1/util/cdna_alignment_orf_to_genome_orf.pl',  wd_fasta + '.transdecoder.gff3',  gff_file_out_u , wd_fasta] 
    call = subprocess.Popen(com, stdout = gff_file_outfile, stderr=errorFilefile, cwd = wd)
    call.communicate()
    errorFilefile.close()
    gff_file_outfile.close()
    
    listErr = []
    err_file = open(errorFile, "r")
    for line in err_file:
        if line.startswith("Warning"):
             listErr.append(("mRNA" + line.split("::")[1]).split(".")[0])
    listErrUniq = list(set(listErr))    
    
    for evm in listErrUniq:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    
    outputFilenameFinal = wd + 'finalAnnotation.Final.gff3'
    outfile = open(outputFilenameFinal, "w")
    com = ['cat', outputFilenameLeft, outputFilename] 
    call = subprocess.Popen(com, stdout = outfile, cwd = wd)
    call.communicate()
    outfile.close()

    return outputFilenameFinal

def strand(gff_file1, gff_file2, fasta, proc, wd): 
    outputFilename = wd + 'finalAnnotation.gff3'
    gff_out = gffwriter.GFFWriter(outputFilename)
    outputFilenameGmap = wd + 'finalAnnotation.gmap.sing.gff3'
    gff_out_s = gffwriter.GFFWriter(outputFilenameGmap)
    
    gff_file1_out = gff_file1 + ".intron.tidy.sorted.gff"
    errorFile = gff_file1 + ".gt_err.log"
    gt_com = ['gt', 'gff3', '-sort', '-tidy', '-addintrons', '-retainids', gff_file1]
    file1 = open(gff_file1_out, 'w')
    err1 = open(errorFile, 'w')
    gt_call = subprocess.Popen( gt_com, stdout=file1, stderr=err1)
    gt_call.communicate()
    
    gff_file2_out = gff_file2 + ".intron.tidy.sorted.gff"
    errorFile = gff_file2 + ".gt_err.log"
    file1 = open(gff_file2_out, 'w')
    err1 = open(errorFile, 'w')
    gt_com = ['gt', 'gff3', '-sort', '-tidy', '-addintrons', '-retainids', gff_file2]
    gt_call = subprocess.Popen( gt_com, stdout=file1, stderr=err1)
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
    #print (len(set(listgenetotal)))
    listgene2 = []
    listgeneintrons = []
    listgenetotal = []

    for i in db2.features_of_type("intron"):
        g = ' '.join (i.attributes['Parent'])
        listgeneintrons.append(g)
    for i in db2.features_of_type("CDS"):
        g = ' '.join (i.attributes['Parent'])
        listgenetotal.append(g)
    listgene2 = sorted(set(list(set(listgenetotal)^set(listgeneintrons))))
    
    newlist = []
    geneDict = {}
    for a in listgene2:
        b =  a.split('_', 1)[1]
        bb = b.split('.')
        del bb[-1]
        evm = '.'.join(bb)
        newlist.append(evm)
        if evm in geneDict:
            z = geneDict[evm]
            geneDict[evm] = z + [a]
        else:
            geneDict[evm] =  [a]
    commonlist = list(set(listgene1).intersection(newlist))
    uniqGmap = sorted(set(list(set(newlist)^set(commonlist))))

    evmList = []
    gmapList = []
    count = 0
    for a in commonlist:
        if geneDict[a] and len(geneDict[a]) < 2:
            evmList.append(a)
        elif geneDict[a] and len(geneDict[a]) > 1:
            gmapList = gmapList + geneDict[a]
    
    for a in uniqGmap:
        if geneDict[a]:
            gmapList = gmapList + geneDict[a]
    listgeneintronsU = []
    listgeneintronsU = (set(listgeneintrons))

    for evm in evmList:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
       
    for evm in gmapList:
        for i in db2.children(evm, featuretype='CDS', order_by='start'):
            gff_out_s.write_rec(i) 
        gff_out_s.write_rec(db2[evm])
        for i in db2.parents(evm, featuretype='gene', order_by='start'):
            gff_out_s.write_rec(i)
        for i in db2.children(evm, featuretype='exon', order_by='start'):
            gff_out_s.write_rec(i)
    
    for evm in listgeneintronsU:
        for i in db2.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i) 
        gff_out.write_rec(db2[evm])
        for i in db2.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db2.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    gff_out.close()
    gff_out_s.close()

    gffOrf = longest(outputFilenameGmap, fasta, proc, wd)
    
    outputFilenameFinal = wd + 'finalAnnotation.Final.Comb.gff3'
    outfile = open(outputFilenameFinal, "w")
    com = ['cat', gffOrf, outputFilename] 
    call = subprocess.Popen(com, stdout = outfile, cwd = wd)
    call.communicate()
    outfile.close()
    return outputFilenameFinal

def exonerate(ref, gff_file, proc, wd):
    exon_file_out = gff_file + ".exons.fasta"
    prot_file_out = gff_file + ".prot.fasta"
    errorFile = gff_file + ".gt_err.log"
    com = ['gffread', '-W','-g', ref, '-w', exon_file_out, '-y', prot_file_out, gff_file]
    fasta_file_outfile = open(exon_file_out, "w")
    errorFilefile = open(errorFile, "w")
    call = subprocess.Popen(com, stdout= fasta_file_outfile , stderr=errorFilefile)
    call.communicate()
    fasta_file_outfile.close()
    errorFilefile.close()
    
    listComplete = []
    listIncomplete = []
    listFasta = []
    dictFastaProt = {}
    for record in SeqIO.parse(prot_file_out, "fasta"):
        if (record.seq).startswith("M") and (record.seq).endswith("."):
                listComplete.append(record.id)
        elif (record.seq).startswith("M") and not (record.seq).endswith("."):
                newseq = (record.seq).split(".")
                if len(newseq) > 1:
                    newseqM = newseq[0] + "."
                    record.seq = newseqM
                    if len(record.seq) > 10:
                        listIncomplete.append(record.id)
                        dictFastaProt[record.id] = record
                        #listFasta.append(record)
        elif not (record.seq).startswith("M") and (record.seq).endswith("."):
                newseq = (record.seq).split("M")
                if len(newseq) > 1:
                    newseqM =  "M" + newseq[1] 
                    record.seq = newseqM
                    if len(record.seq) > 10:
                        listIncomplete.append(record.id)
                        dictFastaProt[record.id] = record
                        #listFasta.append(record)
        elif not (record.seq).startswith("M") and not (record.seq).endswith("."):
                newseq = (record.seq).split("M")
                if len(newseq) > 1:
                    newseqM =  "M" + newseq[1] 
                    record.seq = newseqM
                    newseq = (record.seq).split(".")
                    if len(newseq) > 1:
                        newseqM = newseq[0] + "."
                        record.seq = newseqM
                        if len(record.seq) > 10:
                            listIncomplete.append(record.id)
                            dictFastaProt[record.id] = record
                            #listFasta.append(record)
    #prot_file_out_out = prot_file_out + ".mod.fasta"
    #SeqIO.write(listFasta, prot_file_out_out, "fasta")
    outputFilenameGff = wd + 'Annotation.gff3'
    gff_out = gffwriter.GFFWriter(outputFilenameGff)
    db1 = gffutils.create_db(gff_file, ':memory:',  merge_strategy='create_unique', keep_order=True)
    for evm in listComplete:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    gff_out.close()
    commandList = []
    for record in SeqIO.parse(exon_file_out, "fasta"):
        if record.id in listIncomplete:
            outputFilenameProt = wd + record.id +'.prot.fasta'
            SeqIO.write(dictFastaProt[record.id], outputFilenameProt, "fasta")
            listFields = (record.description).split(' ')
            for elem in listFields:
                outputFilename = wd + record.id +'.genome.fasta'
                bedFile = wd + record.id + '.genome.bed'
                if (elem.startswith('loc') and elem.endswith('+')) or (elem.startswith('loc') and elem.endswith('-')):
                    coordsList = elem.split('|', -2)
                    chrN = coordsList[0].split(':')
                    coord = coordsList[1].split('-')
                    locus = '\t'.join([chrN[1],coord[0],coord[1]])
                    locus = locus + '\n'
                    bedhandler = open(bedFile, 'w')
                    bedhandler.write(locus)
                    bedhandler.close()
                    com = ['bedtools', 'getfasta','-fi', ref, '-bed', bedFile, '-fo', outputFilename]
                    call = subprocess.Popen(com) #, stdout= fasta_file_outfile , stderr=errorFilefile)
                    call.communicate()
                    combList = [outputFilenameProt, outputFilename]
            commandList.append(combList)    

    with Pool(int(proc)) as p:
        p.map(runExonerate, commandList)
    
    listGff3 = []
    for root, dirs, files, in os.walk(wd):
        for fileN in files:
            if fileN.endswith('gff3'):
                listGff3.append(os.path.join(root, fileN))
    
    orintedFIle = wd + '/oriented.gff3'
    dataGff3 = open(orintedFIle, 'w')
    orintedFIleErr = wd + '/oriented.gff3.error'
    dataGff3Err = open(orintedFIleErr, 'w')
    comcat = ['cat'] + listGff3
    gt_com = ['gt', 'gff3', '-sort', '-tidy']
    callcat = subprocess.Popen(comcat, stdout= subprocess.PIPE)
    callgt = subprocess.Popen(gt_com, stdin = callcat.stdout, stdout = dataGff3, stderr=dataGff3Err)
    callgt.communicate()
    dataGff3.close()
    dataGff3Err.close()
    
    return orintedFIle
        
def runExonerate(commandList):
    #print (commandList)
    outputFilenameProt = commandList[0]
    outputFilename = commandList[1]
    protGff = outputFilenameProt + ".exonOut"
    errorFile = outputFilenameProt + ".exonerate_err.log"
    protGff_outfile = open(protGff, "w")
    errorFilefile = open(errorFile, "w")
    com = ['exonerate', '--model', 'protein2genome', '--bestn', '1', '--showtargetgff', 'yes', '--query',  outputFilenameProt, '--target', outputFilename] 
    call = subprocess.Popen(com, stdout= protGff_outfile, stderr=errorFilefile)
    call.communicate()
    fileGff = open(protGff, "r")
    protGff3 = protGff + ".gff3"
    fileFinalGff = open(protGff3, "w")
    gff = fileGff.readlines()
    for line in gff:
        if "exonerate:protein2genome:local" in line:
            splitLine = line.split("\t")
            #print (splitLine)
            if (len(splitLine))> 8:
                chNr = splitLine[0].split(":")
                start = chNr[1].split("-")[0]
                elm = splitLine[8].split(";")
                if "gene" in splitLine[2]:
                    nameGene = elm[1].split(" ")
                    geneList = [chNr[0],"LoReAn" , splitLine[2], str(int(splitLine[3]) + int(start)), str(int(splitLine[4]) + int(start)), '.', splitLine[6], splitLine[7], "ID=" + nameGene[2]]
                    mRNAList = [chNr[0],"LoReAn" , 'mRNA', str(int(splitLine[3]) + int(start)), str(int(splitLine[4]) + int(start)), '.', splitLine[6], splitLine[7], "ID=" + nameGene[2] + ".mRNA;Parent=" + nameGene[2]]
                    fileFinalGff.write(('\t'.join(geneList)) + "\n")
                    fileFinalGff.write(('\t'.join(mRNAList)) + "\n")
                if "cds" in splitLine[2]:
                    cdsList = [chNr[0],"LoReAn" , "CDS", str(int(splitLine[3]) + int(start)), str(int(splitLine[4]) + int(start)), splitLine[5], splitLine[6], splitLine[7], "Parent=" + nameGene[2] + ".mRNA"]
                    fileFinalGff.write(('\t'.join(cdsList)) + "\n")
                if "exon" in splitLine[2]:
                    exonList = [chNr[0],"LoReAn" , splitLine[2], str(int(splitLine[3]) + int(start)), str(int(splitLine[4]) + int(start)), splitLine[5], splitLine[6], splitLine[7], "Parent=" + nameGene[2] + ".mRNA"]
                    fileFinalGff.write(('\t'.join(exonList)) + "\n")
    fileFinalGff.close()




if __name__ == '__main__':
    
    
    evmgff = argv[1]
    #gmap = argv[2]
    fasta = argv[2]
    proc = argv[3]
    #pref = argv[5] 
    wd = argv[4]
    #gffR = strand(evmgff, gmap, fasta, proc, wd)
    #gffPasa = appendID(gffR)
    #noOverl = removeOverlap(gffPasa)
    ##simplified = grs.parseGff(finalOutput)
    #noDisc = removeDiscrepancy(noOverl, evmgff)
    #uniqGene = newNames(noDisc)
    #genename(uniqGene, pref)
    exonerate(fasta, evmgff, proc, wd)
    
