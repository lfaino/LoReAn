#!/usr/bin/env python3

import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from multiprocessing import Pool

import gffutils
import gffutils.gffwriter as gffwriter
from Bio import Seq
from Bio import SeqIO

#======================================================================================================================

GT_GFF3 = 'gt gff3 -sort -tidy %s'

GT_GFF3_R = 'gt gff3 -retainids -sort -tidy %s'

BEDTOOLS_GET_FASTA = 'bedtools getfasta -fi %s -bed %s -fo %s'

GFFREAD_W = 'gffread -W -g %s -w %s -y %s %s'

GFFREAD_M = 'gffread -F -M -C -K -Q -Z -o %s %s'

EXONERATE = 'exonerate --model coding2genome --bestn 1 --refine region --showtargetgff yes --query %s --target %s'

PASA_VAL = 'pasa_gff3_validator.pl %s'

BEDTOOLS_SORT = 'bedtools sort -i %s'

BEDTOOLS_MERGE = 'bedtools merge -d -100 -s -o count,distinct -c 9,9'

GT_GFF3_INTRON = 'gt gff3 -sort -addintrons -tidy %s'

BEDTOOLS_INTERSECT = 'bedtools intersect -v -a %s -b %s'

GT_RETAINID = 'gt gff3 -sort -tidy -addintrons -retainids %s'

GT_GFF3TOGTF = 'gt gff3_to_gtf %s'

BEDTOOLS_GETFASTA =  'bedtools getfasta -fi %s -bed %s -fo %s -name -split'

CAT = 'cat %s'

#======================================================================================================================

gene_count = 0
exon_cds_count = 0
chrold = ""

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


def genename_last(gff_filename, prefix, verbose, wd, dict_ref_name, step):

    global prefix_name
    prefix_name = prefix
    out = tempfile.NamedTemporaryFile(delete=False, mode="w", dir=wd)
    err = tempfile.NamedTemporaryFile(delete=False, mode="w")
    gt_com = GT_GFF3 % gff_filename
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=out, stderr=err, shell=True)
    gt_call.communicate()

    db1 = gffutils.create_db(out.name, ':memory:', merge_strategy='create_unique', keep_order=True, transform=transform_name)
    gene_count = 0
    list_mrna = [mRNA.attributes["ID"][0] for mRNA in db1.features_of_type('mRNA')]
    out_gff = tempfile.NamedTemporaryFile(delete=False, prefix="gffread", suffix=".gff3", dir=wd)
    gff_out = gffwriter.GFFWriter(out_gff.name)
    gene_name = []
    for evm in list_mrna:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            if i.chrom in dict_ref_name:
                i.chrom = dict_ref_name[i.chrom]
            gff_out.write_rec(i)
        i = db1[evm]
        if i.chrom in dict_ref_name:
            i.chrom = dict_ref_name[i.chrom]
        gff_out.write_rec(i)
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            if i.chrom in dict_ref_name:
                i.chrom = dict_ref_name[i.chrom]
            id_gene = i.attributes['ID'][0]
            if not id_gene in gene_name:
                gff_out.write_rec(i)
                gene_name.append(id_gene)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            if i.chrom in dict_ref_name:
                i.chrom = dict_ref_name[i.chrom]
            gff_out.write_rec(i)
    gff_out.close()
    if "pasa" in step:
        out_name = os.path.join(wd, "Final.evm.update.gff3")
        with open(out_name, "w") as fh:
            gt_com = GT_RETAINID % out_gff.name
            if verbose:
                sys.stderr.write('Executing: %s\n\n' % gt_com)
            gt_call = subprocess.Popen(gt_com, stdout=fh, stderr=err, shell=True)
            gt_call.communicate()
    if "lorean" in step:
        out_name = os.path.join(wd, "Final.LoReAn.update.gff3")
        with open(out_name, "w") as fh:
            gt_com = GT_RETAINID % out_gff.name
            if verbose:
                sys.stderr.write('Executing: %s\n\n' % gt_com)
            gt_call = subprocess.Popen(gt_com, stdout=fh, stderr=err, shell=True)
            gt_call.communicate()
    if verbose:
        print(out_name)
    return out_name


def genename(gff_filename, prefix, verbose, wd):
    global prefix_name
    prefix_name = prefix
    db1 = gffutils.create_db(gff_filename, ':memory:', merge_strategy='create_unique', keep_order=True,
                             transform=transform_name)
    gene_count = 0
    list_gene = [mRNA.attributes["ID"][0] for mRNA in db1.features_of_type('gene')]
    out_gff = tempfile.NamedTemporaryFile(delete=False, prefix="gffread", suffix=".gff3", dir=wd)
    gff_out = gffwriter.GFFWriter(out_gff.name)
    for evm in list_gene:
        for i in db1.children(evm, order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
    gff_out.close()
    out = tempfile.NamedTemporaryFile(delete=False, mode="w",prefix="final.", suffix=".gff3", dir = wd)
    err = tempfile.NamedTemporaryFile(delete=False, mode="w" , prefix="final.", suffix=".err", dir = wd)
    gt_com = GT_GFF3_R % out_gff.name
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=out, stderr=err, shell=True)
    gt_call.communicate()
    if verbose:
        print(out.name)

    return out.name


def transform_name(f):
    global gene_count
    global chrold
    global exon_cds_count
    chrold = ""
    if f.featuretype == 'gene':
        f.source = "LoReAn"
        if chrold == "" or f.chrom in chrold:
            chrold = f.chrom
            gene_count += 10
            fatureattgene = prefix_name + '_' + chrold + '_G_' + str(gene_count)
            f.attributes['Name'] = fatureattgene
        else:
            chrold = f.chrom
            gene_count = 0
            gene_count += 10
            fatureattgene = prefix_name + '_' + chrold + '_G_' + str(gene_count)
            f.attributes['Name'] = fatureattgene
        for field in f.attributes:
            if "ID" in field:
                id = f.attributes["ID"]
            elif "Name" in field:
                name = f.attributes["Name"]
        f.attributes = ""
        gene = {}
        gene["ID"] = id
        gene["Name"] = name
        f.attributes = gene

    elif f.featuretype == 'mRNA':
        f.source = "LoReAn"
        for field in f.attributes:
            if "ID" in field:
                id = f.attributes["ID"]
            elif "Parent" in field:
                parent = f.attributes["Parent"]
        f.attributes = ""
        mRNA = {}
        mRNA["ID"] = id
        mRNA["Parent"] = parent
        f.attributes = mRNA
    else:
        f.source = "LoReAn"
        exon_cds_count += 1
        for field in f.attributes:
            if "Parent" in field:
                parent = f.attributes["Parent"]
        f.attributes = ""
        code = {}
        code["Parent"] = parent
        code["ID"] = "LoReAn-" + str(exon_cds_count)
        f.attributes = code
    return f


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
    call = subprocess.Popen(com, stdout=fasta_file_outfile, cwd=wd, stderr=errorFilefile, shell=True)
    call.communicate()
    fasta_file_outfile.close()
    errorFilefile.close()
    list_complete = []
    dict_incomplete = {}
    dict_fasta_prot = {}
    longest_prot = []
    list_total = []
    list_single_exons = []
    for record in SeqIO.parse(prot_file_out, "fasta"):
        list_total.append(record.id)
        if record.seq.startswith("M"):# and record.seq.endswith("."):
            list_complete.append(record.id)
        else:
            dict_incomplete[record.id] = record.id
    for record in SeqIO.parse(exon_file_out, "fasta"):
        list_fields = record.description.split(' ')
        for elem in list_fields:
            if elem.startswith('exons'):
                exon_number = elem.split(",")
                ## changed here for all genes
                if (len(exon_number)) > 0:
                    list_single_exons.append(record.id)
                    if record.id in dict_incomplete:
                        newrecord = record.reverse_complement()
                        input_seq = str(record.seq)
                        startP = re.compile('ATG')
                        nuc = input_seq.replace('\n', '')
                        longest = (0,)
                        for m in startP.finditer(nuc):
                            if len(Seq.Seq(nuc)[m.start():].translate(to_stop=True)) > longest[0]:
                                pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
                                longest = [len(pro), nuc[m.start():m.start() + len(pro) * 3 + 3]]
                                if len(longest) == 2:
                                    record.seq = Seq.Seq(longest[1])
                                    dict_fasta_prot[record.id] = record
                                else:
                                    dict_fasta_prot[record.id] = record
                        input_seq = str(newrecord.seq)
                        startP = re.compile('ATG')
                        nuc = input_seq.replace('\n', '')
                        longest = (0,)
                        for m in startP.finditer(nuc):
                            if len(Seq.Seq(nuc)[m.start():].translate(to_stop=True)) > longest[0]:
                                pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
                                longest = [len(pro), nuc[m.start():m.start() + len(pro) * 3 + 3]]
                                if len(longest) == 2:
                                    if record.id in dict_fasta_prot:
                                        if (len((dict_fasta_prot[record.id]).seq)) < (len(longest[1])):
                                            record.seq = Seq.Seq(longest[1])
                                            dict_fasta_prot[record.id] = record
                                    elif len(longest) == 2:
                                        record.seq = Seq.Seq(longest[1])
                                        dict_fasta_prot[record.id] = record
                                    else:
                                        dict_fasta_prot[record.id] = record

    for mod in dict_fasta_prot:
        longest_prot.append(dict_fasta_prot[mod])
    prot_file_out_mod = prot_file_out + ".mod.fasta"
    SeqIO.write(longest_prot, prot_file_out_mod, "fasta")

    pool = Pool(processes=int(proc))
    list_get_seq= []
    list_short = []
    record_dict = SeqIO.to_dict(SeqIO.parse(exon_file_out, "fasta"))
    for key in dict_fasta_prot:
        if key in record_dict:
            list_short.append(key)
            output_filename_prot = wd + key + '.prot.fasta'
            SeqIO.write(dict_fasta_prot[key], output_filename_prot, "fasta")
            list_fields = record_dict[key].description.split(' ')
            for elem in list_fields:
                output_filename = wd + key + '.genome.fasta'
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
                    data = [ref, bedFile, output_filename, output_filename_prot, verbose, wd]
                    list_get_seq.append(data)
    results_get = pool.map(get_fasta, list_get_seq, chunksize=1)
    results = pool.map(runExonerate, results_get, chunksize=1)
    output_filename_gff = wd + 'mRNA_complete_gene_Annotation.gff3'
    exonerate_files = results + [output_filename_gff]

    listInGff = list_complete + list_short
    listAsbent = sorted(set(list(set(list_total) ^ set(listInGff))))
    listCompleteAll = listAsbent + list_complete
    gff_out = gffwriter.GFFWriter(output_filename_gff)
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


def collaps_locus(orintedFIleN, wd, verbose):

    log = tempfile.NamedTemporaryFile(delete=False, prefix="gffread", dir=wd)
    err = tempfile.NamedTemporaryFile(delete=False, prefix="gffread", dir=wd)
    out = tempfile.NamedTemporaryFile(delete=False, prefix="gffread", suffix= ".gff3", dir=wd)

    com = GFFREAD_M % (out.name, orintedFIleN)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=log, cwd = wd, stderr=err, shell=True)
    call.communicate()

    db1 = gffutils.create_db(out.name, ':memory:', merge_strategy='create_unique', keep_order=True, transform=transform)
    list_mrna = [mRNA.attributes["ID"][0] for mRNA in db1.features_of_type('mRNA')]
    out = tempfile.NamedTemporaryFile(delete=False, prefix="gffread", suffix= ".gff3", dir=wd)
    gff_out = gffwriter.GFFWriter(out.name)
    gene_name = []
    for evm in list_mrna:
        for i in db1.children(evm, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(db1[evm])
        for i in db1.parents(evm, featuretype='gene', order_by='start'):
            id_gene = i.attributes['ID'][0]
            if not id_gene in gene_name:
                gff_out.write_rec(i)
                gene_name.append(id_gene)
        for i in db1.children(evm, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    gff_out.close()

    return out.name


def transform(f):
    if f.featuretype == 'locus':
        f.featuretype = "gene"
        f.source = "LoReAn"
    elif "mRNA" in f.featuretype:
        f.attributes['Parent'] = f.attributes['locus']
    return f


def get_fasta(list_data):
    err = tempfile.NamedTemporaryFile(delete=True)
    log = tempfile.NamedTemporaryFile(delete=True)

    com = BEDTOOLS_GET_FASTA % (list_data[0], list_data[1], list_data[2])
    if list_data[4]:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, cwd = list_data[5], shell=True, stderr=err, stdout=log)
    call.communicate()
    os.remove(list_data[1])
    combList = [list_data[3], list_data[2], list_data[4], list_data[5]]
    return combList


def runExonerate(commandList):
    outputFilenameProt = commandList[0]
    outputFilename = commandList[1]
    protGff = outputFilenameProt + ".exonOut"
    errorFile = outputFilenameProt + ".exonerate_err.log"
    protGff_outfile = open(protGff, "w")
    errorFilefile = open(errorFile, "w")
    com = EXONERATE % (outputFilenameProt, outputFilename)
    if commandList[2]:
        sys.stderr.write('Executing: %s\n\n' % com)
    call = subprocess.Popen(com, stdout=protGff_outfile, cwd = commandList[3], stderr=errorFilefile, shell=True)
    call.communicate()
    fileGff = open(protGff, "r")
    protGff3 = protGff + ".gff3"
    fileFinalGff = open(protGff3, "w")
    gff = fileGff.readlines()
    for line in gff:
        if "exonerate:" in line:
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
    if os.stat(errorFile).st_size == 0:
        os.remove(errorFile)
        os.remove(protGff)
    os.remove(outputFilenameProt)
    os.remove(outputFilename)

    return protGff3


def genename_evm(gff_filename, verbose, wd):

    gene_evm = "evm.TU."
    mRNA_evm = "evm.model."
    out = tempfile.NamedTemporaryFile(delete=False, mode="w", dir=wd)
    err = tempfile.NamedTemporaryFile(delete=False, mode="w")
    gt_com = GT_GFF3 % gff_filename
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=out, stderr=err, shell=True)
    gt_call.communicate()

    db1 = gffutils.create_db(out.name, ':memory:', merge_strategy='create_unique', keep_order=True)
    list_gene = [gene.attributes["ID"][0] for gene in db1.features_of_type('gene')]
    list_chr = [chr.chrom for chr in db1.features_of_type('gene')]
    chr_count_gene = {}
    chr_count_mRNA = {}

    chrs = list(set(list_chr))
    for elm in chrs:
        chr_count_gene[elm] = 0
        chr_count_mRNA[elm] = 0

    out_gff = tempfile.NamedTemporaryFile(delete=False, prefix="gffread_parse", suffix=".gff3", dir=wd)
    gff_out = gffwriter.GFFWriter(out_gff.name)
    for evm in list_gene:
        exon_count = 0
        cds_count = 0
        gene_chr = db1[evm].chrom
        gene = db1[evm]
        count_gene = chr_count_gene[gene_chr] + 1
        chr_count_gene[gene_chr] = count_gene
        id_new_gene = gene_evm + gene_chr + "." + str(count_gene)
        gene.attributes["ID"][0] = id_new_gene
        gff_out.write_rec(gene)
        for i in db1.children(evm, featuretype='mRNA', order_by='start'):
            mRNA = i
            mRNA_old = i.attributes["ID"][0]
            count_mRNA = chr_count_mRNA[gene_chr] + 1
            chr_count_mRNA[gene_chr] = count_mRNA
            id_new_mRNA = mRNA_evm + gene_chr + "." + str(count_mRNA)
            mRNA.attributes["Parent"][0] = id_new_gene
            mRNA.attributes["ID"][0] = id_new_mRNA
            gff_out.write_rec(mRNA)
            for e in db1.children(mRNA_old, featuretype='exon', order_by='start'):
                exon_count += 1
                exon = e
                exon.attributes["Parent"][0] = mRNA.attributes["ID"][0]
                exon.attributes["ID"] = id_new_mRNA + ".exon" + str(exon_count)
                gff_out.write_rec(exon)
            for c in db1.children(mRNA_old, featuretype='CDS', order_by='start'):
                cds_count += 1
                cds = c
                cds.attributes["Parent"][0] = mRNA.attributes["ID"][0]
                cds.attributes["ID"] = "cds." + str(cds_count) + "." + id_new_mRNA
                gff_out.write_rec(cds)
    gff_out.close()

    out = tempfile.NamedTemporaryFile(delete=False, mode="w", prefix="new_name_update.", suffix= ".gff3", dir=wd)
    err = tempfile.NamedTemporaryFile(delete=False, mode="w", dir=wd)
    gt_com = GT_GFF3_R % out_gff.name
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=out, stderr=err, shell=True)
    gt_call.communicate()
    if verbose:
        print(out.name)
    return out.name


cds_lorean = "LoReAn.CDS."
exons_lorean = "LoReAn.exons."
other_lorean = "LoReAn.other."

def genename_lorean(gff_filename, verbose, wd):

    out = tempfile.NamedTemporaryFile(delete=False, mode="w", dir=wd)
    err = tempfile.NamedTemporaryFile(delete=False, mode="w")
    gt_com = GT_GFF3 % gff_filename
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % gt_com)
    gt_call = subprocess.Popen(gt_com, stdout=out, stderr=err, shell=True)
    gt_call.communicate()

    outfile_gff = tempfile.NamedTemporaryFile(delete=False, prefix="additional.", suffix=".gff3", dir=wd)
    log = tempfile.NamedTemporaryFile(delete=False, prefix="uniq.ID.pasa.", suffix=".log", dir=wd)
    err = tempfile.NamedTemporaryFile(delete=False, prefix="uniq.ID.pasa.", suffix=".err", dir=wd)

    cmd = GFFREAD_M % (outfile_gff.name, out.name)

    if verbose:
        sys.stderr.write('Executing: %s\n\n' % cmd)
    gffread = subprocess.Popen(cmd, cwd=wd, shell=True, stdout=log, stderr=err)
    gffread.communicate()

    db_gffread = gffutils.create_db(out.name, ':memory:', merge_strategy='create_unique', keep_order=True, transform=transform_cds)
    outfile_out = tempfile.NamedTemporaryFile(delete=False, prefix="uniq.ID.pasa.final.", suffix=".gff3", dir=wd)
    gff_out_s = gffwriter.GFFWriter(outfile_out.name)

    for gene in db_gffread.features_of_type("gene"):
        gff_out_s.write_rec(db_gffread[gene])
        for i in db_gffread.children(gene, order_by='start'):
            gff_out_s.write_rec(i)
    return outfile_out.name

cds_count_lorean = 0

def transform_cds(x):
    global cds_count_lorean
    cds_count_lorean += 1
    if 'gene' in x.featuretype:
        x.attributes = {'ID': x.attributes['ID'], 'Name': x.attributes['Name'][0].split("_")[0]}
        x.source = "LoReAn"
        return x
    elif 'mRNA' in x.featuretype:
        x.source = "LoReAn"
        x.attributes = {'ID': x.attributes['ID'], 'Parent': x.attributes['Parent']}
        return x
    elif 'CDS' in x.featuretype:

        x.attributes = {'ID': cds_lorean + str(cds_count_lorean), 'Parent': x.attributes['Parent']}
        x.source = "LoReAn"
        return x
    elif 'exon' in x.featuretype:
        x.attributes = {'ID': exons_lorean + str(cds_count_lorean), 'Parent': x.attributes['Parent']}
        x.source = "LoReAn"
        return x
    else:
        x.attributes = {'ID': other_lorean + str(cds_count_lorean), 'Parent': x.attributes['Parent']}
        x.source = "LoReAn"
        return x



def add_removed_evm(pasa, exon, wd):
    """
    here the clusters of sequence from the same locus are prepared
    """

    db_evm = gffutils.create_db(pasa, ':memory:', merge_strategy='create_unique', keep_order=True)
    ids_evm = [gene.attributes["ID"][0] for gene in db_evm.features_of_type("mRNA")]

    db_gmap = gffutils.create_db(exon, ':memory:', merge_strategy='create_unique', keep_order=True)
    ids_gmap_full = [gene.attributes["ID"][0] for gene in db_gmap.features_of_type("gene")]
    ids_gmap = [gene.attributes["ID"][0].split("_")[0] for gene in db_gmap.features_of_type("mRNA")]


    uniq_evm = [evm for evm in ids_evm if evm not in ids_gmap]
    uniq_gene = [gene.attributes["ID"][0] for mrna in uniq_evm for gene in db_evm.parents(mrna)]
    uniq = list(set(uniq_gene))


    outfile = tempfile.NamedTemporaryFile(delete=False, prefix="additional.", suffix=".gff3", dir=wd)
    gff_out_s = gffwriter.GFFWriter(outfile.name)

    for name in uniq:
        for i in db_evm.children(name, order_by='start'):
            gff_out_s.write_rec(i)
        gff_out_s.write_rec(db_evm[name])
    for name in ids_gmap_full:
        for i in db_gmap.children(name, order_by='start'):
            gff_out_s.write_rec(i)
        gff_out_s.write_rec(db_gmap[name])
    gff_out_s.close()

    return outfile.name


if __name__ == '__main__':
    #strand(*sys.argv[1:])
    #exonerate(fasta, outputFilename, proc, gmap_wd, verbose) genename_evm(gff_filename, verbose, wd)
    genename_lorean(*sys.argv[1:])
