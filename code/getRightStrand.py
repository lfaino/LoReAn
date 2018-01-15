#!/usr/bin/env python3

import re
import shutil
import subprocess
import sys
import tempfile
import time
import warnings
from multiprocessing import Pool, Manager

import gffutils
import gffutils.gffwriter as gffwriter
import progressbar
from Bio import Seq
from Bio import SeqIO

#======================================================================================================================

GT_GFF3 = 'gt gff3 -sort -tidy %s'

BEDTOOLS_GET_FASTA = 'bedtools getfasta -fi %s -bed %s -fo %s'

GFFREAD_W = 'gffread -W -g %s -w %s -y %s %s'

EXONERATE = 'exonerate --model protein2genome --bestn 1 --showtargetgff yes --query %s --target %s'

PASA_VAL = 'pasa_gff3_validator.pl %s'

BEDTOOLS_SORT = 'bedtools sort -i %s'

BEDTOOLS_MERGE = 'bedtools merge -d -100 -s -o count,distinct -c 9,9'

GT_GFF3_INTRON = 'gt gff3 -sort -addintrons -tidy %s'

BEDTOOLS_INTERSECT = 'bedtools intersect -v -a %s -b %s'

GT_RETAINID = 'gt gff3 -sort -tidy -addintrons -retainids %s'

GT_GFF3TOGTF = 'gt gff3_to_gtf %s'

#======================================================================================================================


def removeDiscrepancy(gff, evmFile, verbose):
    badName = []
    comm = PASA_VAL % gff
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % comm)
    gffVal_call = subprocess.Popen(comm, stdout=subprocess.PIPE, shell=True)
    for ln in gffVal_call.stdout.readlines():
        name = re.split(' |\.CDS', (ln.decode("utf-8")))
        if len(name) > 3 and "ERROR" in name[0]:
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
    listgene = sorted(set(list(set(listAllName) ^ set(badNameUniq))))
    outputFilename = gff + '.noProblem.gff3'
    gff_out = gffwriter.GFFWriter(outputFilename)
    db1 = gffutils.create_db(gff, ':memory:', merge_strategy='create_unique', keep_order=True)
    for evm in listgene:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    cmd = BEDTOOLS_INTERSECT % (evmFile, gff)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % cmd)
    bedtools_call = subprocess.Popen(cmd , stdout=subprocess.PIPE, shell=True)
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
    db1 = gffutils.create_db(evmFile, ':memory:', merge_strategy='create_unique', keep_order=True)
    for evm in evm_mRNA:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)

    return outputFilename


def appendID(gff):
    i = open(gff, 'r')
    outFile = gff + '.reformat.gff'
    o = open(outFile, 'w')
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
            elif "CDS" in listLine[2]:
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
    return outFile


def longest(gff_file, fasta, proc, wd, verbose):
    outputFilename = wd + 'finalAnnotation.strand.gff3'
    outputFilenameLeft = tempfile.NamedTemporaryFile(delete=False, dir=wd, prefix="longest.")
    gff_out = gffwriter.GFFWriter(outputFilenameLeft.name)

    gt_com = GT_GFF3_INTRON % gff_file
    gff_file_outfile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".out") #open(gff_file_out, "w")
    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".err") #open(errorFile, "w")
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=gff_file_outfile, stderr=errorFilefile, shell=True)
    gt_call.communicate()

    gt_com = GT_GFF3TOGTF % gff_file_outfile.name
    gtf_file_outfile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".out") #open(gtf_file_out, "w")
    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".err") #open(errorFile, "w")
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=gtf_file_outfile, stderr=errorFilefile, shell=True)
    gt_call.communicate()


    db1 = gffutils.create_db(gff_file_outfile.name, ':memory:', merge_strategy='create_unique', keep_order=True)


    fasta_file_outfile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".out") #open(fasta_file_out, "w")
    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".err") #open(errorFile, "w")
    com = 'cufflinks_gtf_genome_to_cdna_fasta.pl %s %s' % (gtf_file_outfile.name, fasta)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=fasta_file_outfile, stderr=errorFilefile, shell=True)
    call.communicate()

    gff_file_outfile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".out") #open(gff_file_out_u, "w")
    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".err") #open(errorFile, "w")
    com = 'cufflinks_gtf_to_alignment_gff3.pl %s' % gtf_file_outfile.name
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=gff_file_outfile, stderr=errorFilefile, shell=True)
    call.communicate()

    gff_file_outfile_1 = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".out") #open(gff_file_out, "w")
    errorFilefile_1 = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".err") #open(errorFile, "w")
    com = 'TransDecoder.LongOrfs -m 10 -t %s' % fasta_file_outfile.name
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=gff_file_outfile_1, stderr=errorFilefile_1, cwd=wd, shell=True)
    call.communicate()

    gff_file_outfile_2 = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".out") #open(gff_file_out, "w")
    errorFilefile_2 = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.", suffix=".err") #open(errorFile, "w")
    wd_fasta = fasta_file_outfile.name
    com = 'TransDecoder.Predict --single_best_orf --cpu %s --retain_long_orfs 10 -t %s' % (proc, wd_fasta)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=gff_file_outfile_2, stderr=errorFilefile_2, cwd=wd, shell=True)
    call.communicate()

    gff_file_outfile_3 = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.") #open(outputFilename, "w")
    errorFilefile_3 = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=wd, prefix="longest.") #open(errorFile, "w")
    transdecoder = tempfile.NamedTemporaryFile(delete=False)

    com = 'cdna_alignment_orf_to_genome_orf.pl %s %s %s' % (wd_fasta + '.transdecoder.gff3', gff_file_outfile.name, wd_fasta)

    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=gff_file_outfile_3, stderr=errorFilefile_3, cwd=wd, shell=True)
    call.communicate()


    listErr = []
    err_file = open(errorFilefile.name, "r")
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

    gff_files = [outputFilenameLeft.name, gff_file_outfile.name]
    outputFilenameFinal = wd + 'finalAnnotation.Transdecoder.gff3'



    with open(outputFilenameFinal, 'wb') as wfd:
        for f in gff_files:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd, 1024 * 1024 * 10)

    return outputFilenameFinal


def removeOverlap(gff, verbose):
    i = open(gff, 'r')
    #outFile = gff + '.RNA.gff'
    o = tempfile.NamedTemporaryFile(delete=False, mode='w') #open(outFile, 'w')
    for line in i:
        listLine = line.split('\t')
        if len(listLine) == 9:
            if "CDS" in listLine[2]:
                o.write(line)
    i.close()
    bedouffile = tempfile.NamedTemporaryFile()
    #errorFile = outFile + ".bedtools_err.log"
    errorFilefile = tempfile.NamedTemporaryFile() #open(errorFile, "w")
    bedsort = BEDTOOLS_SORT % o.name
    bedmerge = BEDTOOLS_MERGE
    o.close()
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % bedsort)
        sys.stderr.write('Executing: %s\n\n' % bedmerge)
    bedsort_call = subprocess.Popen(bedsort, stdout=subprocess.PIPE, stderr=errorFilefile, shell=True)
    bedmerge_call = subprocess.Popen(bedmerge, stdin=bedsort_call.stdout, stdout=bedouffile, stderr=errorFilefile, shell=True)
    bedmerge_call.communicate()
    errorFilefile.close()
    listMultiple = []
    listUniq = []
    count = 0
    dictRNA = {}
    i = open(bedouffile.name, 'r')
    for a in i:
        listLine = a.split('\t')
        nameRNA = re.split(',|;', listLine[5])
        count += 1
        locus = "locus" + str(count)
        for elm in nameRNA:
            if "Parent" in elm and int(listLine[4]) > 1:
                mRNAname = elm.split('=')[1]
                listMultiple.append(mRNAname)
                if locus in dictRNA:
                    dictRNA[locus].append(mRNAname)
                else:
                    dictRNA[locus] = [mRNAname, ]
            elif "Parent" in elm:
                mRNAname = elm.split('=')[1]
                listUniq.append(mRNAname)
    bedouffile.close()
    listMultipleUniq = []
    listMultipleUniq = list(set(listMultiple))
    dictLength = {}
    mRNA = open(gff, 'r')
    for line in mRNA:
        listLine = line.split('\t')
        if len(listLine) == 9:
            if "CDS" in listLine[2]:
                for key in dictRNA:
                    for el in dictRNA[key]:
                        nameID = "Parent=" + el + ';'
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
    db1 = gffutils.create_db(gff, ':memory:', merge_strategy='create_unique', keep_order=True)
    for evm in listUniqNew:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)

    return outputFilename


def genename(gff_filename, prefix, verbose):


    with tempfile.NamedTemporaryFile(delete=False, mode="w", prefix="rename") as fileout:
        with open(gff_filename, "r") as filein:
            for line in filein:
                if not line.startswith('#'):
                    fileout.write(line)

    errorFilefile = tempfile.NamedTemporaryFile(delete=False, mode = "w", prefix="rename")# open(errorFile, "w")
    tmp_file = tempfile.NamedTemporaryFile(delete=False, mode = "w", prefix="rename")
    cmd = GT_GFF3 % fileout.name
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % cmd)
        sys.stderr.write('Log file is: %s\n\n' % errorFilefile.name)
    gt_call = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=errorFilefile, shell=True)
    chrold = ''
    count = 0
    path = []
    path = gff_filename.split('/')
    path[-1] = prefix + '_LoReAn.annotation.gff3'
    newpath = '/'.join(path)
    o = open(newpath, 'w+')
    o.write('##version-gff 3\n')
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
                o.write('\t'.join(fields) + '\n')

            elif 'mRNA' in typef:
                featureann = fields[8].split(';')
                for elm in featureann:
                    if 'ID' in elm:
                        idname = elm
                    elif 'Parent' in elm:
                        parentname = elm
                fields[8] = idname + ';' + parentname
                o.write('\t'.join(fields) + '\n')


            else:
                o.write('\t'.join(fields) + '\n')
    o.close()
    return newpath


def newNames(oldname):
    finalout = oldname + ".sorted.gff"
    errorFile = oldname + ".gt_err.log"
    finaloutfile = open(finalout, "w")
    errorFilefile = open(errorFile, "w")
    gt_com = GT_GFF3 % oldname
    gt_call = subprocess.Popen(gt_com, stdout=finaloutfile, stderr=errorFilefile, shell=True)
    gt_call.communicate()
    finaloutfile.close()
    errorFilefile.close()
    return finalout


def strand(gff_file1, gff_file2, fasta, proc, gmap_wd, verbose):
    outputFilename = tempfile.NamedTemporaryFile(delete=False, prefix="grs", dir=gmap_wd)
    gff_out = gffwriter.GFFWriter(outputFilename.name)
    outputFilenameGmap = tempfile.NamedTemporaryFile(delete=False, prefix="grs", dir=gmap_wd)
    gff_out_s = gffwriter.GFFWriter(outputFilenameGmap.name)

    gt_com = GT_RETAINID % gff_file1
    file1 = tempfile.NamedTemporaryFile(delete=False, mode="w", prefix="grs", dir=gmap_wd)
    err1 = tempfile.NamedTemporaryFile(delete=False, mode="w", prefix="grs", dir=gmap_wd)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
        sys.stderr.write('Log file is: %s %s\n\n' % (file1.name, err1.name))
    gt_call = subprocess.Popen(gt_com, stdout=file1, stderr=err1, shell=True)
    gt_call.communicate()

    file2 = tempfile.NamedTemporaryFile(delete=False, mode="w")
    err1 = tempfile.NamedTemporaryFile(delete=False, mode="w")
    gt_com = GT_RETAINID % gff_file2
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=file2, stderr=err1, shell=True)
    gt_call.communicate()

    db1 = gffutils.create_db(file1.name, ':memory:', merge_strategy='create_unique', keep_order=True)
    db2 = gffutils.create_db(file2.name, ':memory:', merge_strategy='create_unique', keep_order=True)
    listgeneintrons = []
    listgenetotal = []
    for i in db1.features_of_type("intron"):
        g = ' '.join(i.attributes['Parent'])
        listgeneintrons.append(g)
    for i in db1.features_of_type("CDS"):
        g = ' '.join(i.attributes['Parent'])
        listgenetotal.append(g)
    listgene1 = sorted(set(list(set(listgenetotal) ^ set(listgeneintrons))))
    listgeneintrons = []
    listgenetotal = []

    for i in db2.features_of_type("intron"):
        g = ' '.join(i.attributes['Parent'])
        listgeneintrons.append(g)
    for i in db2.features_of_type("CDS"):
        g = ' '.join(i.attributes['Parent'])
        listgenetotal.append(g)
    listgene2 = sorted(set(list(set(listgenetotal) ^ set(listgeneintrons))))

    newlist = []
    gene_dict = {}
    for a in listgene2:
        b = a.split('_', 1)[1]
        bb = b.split('.')
        del bb[-1]
        evm = '.'.join(bb)
        newlist.append(evm)
        if evm in gene_dict:
            z = gene_dict[evm]
            gene_dict[evm] = z + [a]
        else:
            gene_dict[evm] = [a]
    commonlist = list(set(listgene1).intersection(newlist))
    uniqGmap = sorted(set(list(set(newlist) ^ set(commonlist))))

    evm_list = []
    gmap_list = []
    for a in commonlist:
        if gene_dict[a] and len(gene_dict[a]) < 2:
            evm_list.append(a)
        elif gene_dict[a] and len(gene_dict[a]) > 1:
            gmap_list = gmap_list + gene_dict[a]

    for a in uniqGmap:
        if gene_dict[a]:
            gmap_list = gmap_list + gene_dict[a]
    listgeneintrons_u = (set(listgeneintrons))

    for evm in evm_list:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)

    for evm in gmap_list:
        for i in db2.children(evm, featuretype='CDS', order_by='start'):
            gff_out_s.write_rec(i)
        gff_out_s.write_rec(db2[evm])
        for i in db2.parents(evm, featuretype='gene', order_by='start'):
            gff_out_s.write_rec(i)
        for i in db2.children(evm, featuretype='exon', order_by='start'):
            gff_out_s.write_rec(i)

    for evm in listgeneintrons_u:
        for i in db2.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db2[evm])
        for i in db2.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db2.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    gff_out.close()
    gff_out_s.close()

    gffOrf = longest(outputFilenameGmap.name, fasta, proc, gmap_wd, verbose)

    output_filename_final = gmap_wd + 'finalAnnotation.Final.Comb.gff3'
    gff_files = [gffOrf, outputFilename.name]

    with open(output_filename_final, 'wb') as wfd:
        for f in gff_files:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd, 1024 * 1024 * 10)

    return output_filename_final


def exonerate(ref, gff_file, proc, wd, verbose):

    ##THIS removes the warning. the check of the longest protein was giving a warining. if Biopython change, this could be a problem

    warnings.filterwarnings("ignore")


    exon_file_out = gff_file + ".exons.fasta"
    prot_file_out = gff_file + ".prot.fasta"
    errorFile = gff_file + ".gffread_err.log"
    logFile = exon_file_out + "gffread_log.log"
    com = GFFREAD_W % (ref, exon_file_out, prot_file_out, gff_file)
    fasta_file_outfile = open(logFile, "w")
    errorFilefile = open(errorFile, "w")
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=fasta_file_outfile, cwd = wd, stderr=errorFilefile, shell=True)
    call.communicate()
    fasta_file_outfile.close()
    errorFilefile.close()

    listComplete = []
    dictIncomplete = {}
    dictFastaProt = {}
    longestProt = []
    listTotal = []
    listSingleExons = []
    testDict= {}

    for record in SeqIO.parse(prot_file_out, "fasta"):
        listTotal.append(record.id)
        if record.seq.startswith("M"):# and record.seq.endswith("."):
            listComplete.append(record.id)
        else:
            dictIncomplete[record.id] = record.id
    for record in SeqIO.parse(exon_file_out, "fasta"):
        listFields = record.description.split(' ')
        for elem in listFields:
            if elem.startswith('exons'):
                exonNumber = elem.split(",")
                ## changed here for all genes
                if (len(exonNumber)) > 0:
                    listSingleExons.append(record.id)
                    if record.id in dictIncomplete:
                        newrecord = record.reverse_complement()
                        input_seq = str(record.seq)
                        startP = re.compile('ATG')
                        nuc = input_seq.replace('\n', '')
                        longest = (0,)
                        for m in startP.finditer(nuc):
                            if len(Seq.Seq(nuc)[m.start():].translate(to_stop=True)) > longest[0]:
                                pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
                                longest = [len(pro), m.start(), str(pro), nuc[m.start():m.start() + len(pro) * 3 + 3]]
                                if len(longest) == 4:
                                    record.seq = Seq.Seq(longest[2])
                                    dictFastaProt[record.id] = record
                                else:
                                    dictFastaProt[record.id] = record
                        input_seq = str(newrecord.seq)
                        startP = re.compile('ATG')
                        nuc = input_seq.replace('\n', '')
                        longest = (0,)
                        for m in startP.finditer(nuc):
                            if len(Seq.Seq(nuc)[m.start():].translate(to_stop=True)) > longest[0]:
                                pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
                                longest = [len(pro), m.start(), str(pro), nuc[m.start():m.start() + len(pro) * 3 + 3]]
                                if len(longest) == 4:
                                    if record.id in dictFastaProt:
                                        if (len((dictFastaProt[record.id]).seq)) < (len(longest[2])):
                                            record.seq = Seq.Seq(longest[2])
                                            dictFastaProt[record.id] = record
                                    elif len(longest) == 4:
                                        record.seq = Seq.Seq(longest[2])
                                        dictFastaProt[record.id] = record
                                    else:
                                        dictFastaProt[record.id] = record

    for mod in dictFastaProt:
        longestProt.append(dictFastaProt[mod])
    prot_file_out_mod = prot_file_out + ".mod.fasta"
    SeqIO.write(longestProt, prot_file_out_mod, "fasta")

    manage = Manager()
    queue = manage.Queue()
    pool = Pool(processes=int(proc), maxtasksperchild=10000)

    commandList = []
    listShort = []
    record_dict = SeqIO.to_dict(SeqIO.parse(exon_file_out, "fasta"))
    for key in dictFastaProt:
        if key in record_dict:
            listShort.append(key)
            outputFilenameProt = wd + key + '.prot.fasta'
            SeqIO.write(dictFastaProt[key], outputFilenameProt, "fasta")
            listFields = record_dict[key].description.split(' ')
            for elem in listFields:
                outputFilename = wd + key + '.genome.fasta'
                bedFile = wd + key + '.genome.bed'
                if (elem.startswith('loc') and elem.endswith('+')) or (elem.startswith('loc') and elem.endswith('-')):
                    coordsList = elem.split('|', -2)
                    chrN = coordsList[0].split(':')
                    coord = coordsList[1].split('-')
                    locus = '\t'.join([chrN[1], coord[0], coord[1]])
                    locus = locus + '\n'
                    bedhandler = open(bedFile, 'w')
                    bedhandler.write(locus)
                    bedhandler.close()
                    com = BEDTOOLS_GET_FASTA % (ref, bedFile, outputFilename)
                    if verbose:
                        sys.stderr.write('Executing: %s\n\n' % com)
                    call = subprocess.Popen(com, cwd = wd, shell=True)  # , stdout= fasta_file_outfile , stderr=errorFilefile)
                    call.communicate()
                    combList = [outputFilenameProt, outputFilename, verbose, queue]
            commandList.append(combList)




    results = pool.map_async(runExonerate, commandList)
    with progressbar.ProgressBar(max_value=len(commandList)) as bar:
        while not results.ready():
            size = queue.qsize()
            bar.update(size)
            time.sleep(1)

    outputFilenameGff = wd + 'mRNA_complete_gene_Annotation.gff3'
    exonerate_files = results._value + [outputFilenameGff]

    listInGff = listComplete + listShort
    listAsbent = sorted(set(list(set(listTotal) ^ set(listInGff))))
    listCompleteAll = listAsbent + listComplete


    gff_out = gffwriter.GFFWriter(outputFilenameGff)
    db1 = gffutils.create_db(gff_file, ':memory:', merge_strategy='create_unique', keep_order=True)
    for evm in listCompleteAll:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    gff_out.close()

    orintedFIleN = wd + '/oriented.oldname.gff3'
    with open(orintedFIleN, 'wb') as wfd:
        for f in exonerate_files:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd, 1024 * 1024 * 10)


    return orintedFIleN


def runExonerate(commandList):
    outputFilenameProt = commandList[0]
    outputFilename = commandList[1]
    protGff = outputFilenameProt + ".exonOut"
    errorFile = outputFilenameProt + ".exonerate_err.log"
    protGff_outfile = open(protGff, "w")
    errorFilefile = open(errorFile, "w")
    com = EXONERATE % (outputFilenameProt, outputFilename)
    commandList[3].put(com)
    if commandList[2]:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=protGff_outfile, stderr=errorFilefile, shell=True)
    call.communicate()
    fileGff = open(protGff, "r")
    protGff3 = protGff + ".gff3"
    fileFinalGff = open(protGff3, "w")
    gff = fileGff.readlines()
    for line in gff:
        if "exonerate:protein2genome:local" in line:
            splitLine = line.split("\t")
            # sys.stdout.write (splitLine)
            if (len(splitLine)) > 8:
                chNr = splitLine[0].split(":")
                start = chNr[1].split("-")[0]
                elm = splitLine[8].split(";")
                if "gene" in splitLine[2]:
                    nameGene = elm[1].split(" ")
                    geneList = [chNr[0], "LoReAn", splitLine[2], str(int(splitLine[3]) + int(start)),
                                str(int(splitLine[4]) + int(start)), '.', splitLine[6], splitLine[7],
                                "ID=" + nameGene[2]]
                    mRNAList = [chNr[0], "LoReAn", 'mRNA', str(int(splitLine[3]) + int(start)),
                                str(int(splitLine[4]) + int(start)), '.', splitLine[6], splitLine[7],
                                "ID=" + nameGene[2] + ".mRNA;Parent=" + nameGene[2]]
                    fileFinalGff.write(('\t'.join(geneList)) + "\n")
                    fileFinalGff.write(('\t'.join(mRNAList)) + "\n")
                if "cds" in splitLine[2]:
                    cdsList = [chNr[0], "LoReAn", "CDS", str(int(splitLine[3]) + int(start)),
                               str(int(splitLine[4]) + int(start)), splitLine[5], splitLine[6], splitLine[7],
                               "Parent=" + nameGene[2] + ".mRNA"]
                    fileFinalGff.write(('\t'.join(cdsList)) + "\n")
                if "exon" in splitLine[2]:
                    exonList = [chNr[0], "LoReAn", splitLine[2], str(int(splitLine[3]) + int(start)),
                                str(int(splitLine[4]) + int(start)), splitLine[5], splitLine[6], splitLine[7],
                                "Parent=" + nameGene[2] + ".mRNA"]
                    fileFinalGff.write(('\t'.join(exonList)) + "\n")
    fileFinalGff.close()
    return protGff3
