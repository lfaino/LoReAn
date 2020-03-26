#! /usr/bin/env python3
import subprocess as sb

import gffutils
import gffutils.gffwriter as gffwriter
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# final_evm="/home/lfaino/lfainoData/lorean/LoReAn_Example/Crispa/LoReAn_crispa/crispa_output/Final.LoReAn.update.gff3"
# genome = "/home/lfaino/lfainoData/lorean/LoReAn_Example/Crispa/scaffold3.fasta"


def pep_seq(myFasta, final_evm):
    fasta = {}
    db = gffutils.create_db(final_evm, ':memory:', merge_strategy="create_unique", keep_order=True)
    gff_file = final_evm + ".gff3"
    gff_out = gffwriter.GFFWriter(gff_file)
    for t in db.features_of_type('mRNA', order_by='start'):
        position = []
        seq_combined = ''
        j = 0
        for i in db.children(t, featuretype='CDS', order_by='start'):
            j += 1
            if j == 1:
                pphase = i[7]
            seq = i.sequence(myFasta, use_strand=False)
            seq_combined += seq
            position = position + [i.start,i.stop]
        seq_combined = SeqRecord(Seq(seq_combined, generic_dna))
        if t.strand == '-':
            pphase = i[7]
            seq_combined = seq_combined.reverse_complement()
        if pphase == "0" or pphase == ".":
            seq_transl = seq_combined.translate()
        elif pphase == "1":
            seq_transl = seq_combined[1:].translate()
        elif pphase == "2":
            seq_transl = seq_combined[2:].translate()
        seq_transl.description = position
        seq_transl.id = t.id
        position = sorted(position)
        t.start = position[0]
        t.stop = position[-1]
        fasta[str(seq_transl.seq)] = [seq_transl, t]
    for key in fasta:
        count = 0
        for i in db.parents(fasta[key][1], featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
            i.start = fasta[key][1].start
            i.stop = fasta[key][1].stop
        gff_out.write_rec(fasta[key][1])
        for i in db.children(fasta[key][1], featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        for i in db.children(fasta[key][1], featuretype='CDS', order_by='start'):
            count += 1
            i.featuretype="exon"
            i.attributes["ID"][0] = i.attributes["ID"][0] + "-" + str(count)
            gff_out.write_rec(i)
    return (gff_file)

def gff_filter(final_evm, myFasta):
    file_out = final_evm + ".mod.gff3"
    with open(final_evm, "r") as fh, open(file_out, "w") as fhd:
        for line in fh:
            if not line.startswith("#"):
                elm = line.split("\t")
                elm[8] = elm[8].replace("locus", "Parent")
                elm[2] = elm[2].replace("locus", "gene")
                fhd.write("\t".join(elm))
    db = gffutils.create_db(file_out, ':memory:', merge_strategy="create_unique", keep_order=True)

    b = []
    mrna_retain = []
    for t in db.features_of_type('gene', order_by='start'):
        c = 0
        for i in db.children(t, featuretype='mRNA', order_by='start'):
            c += 1
        if c > 1:
            b.append(t)
        else:
            mrna_retain.append(i.attributes["ID"][0])
    mrna_select = []
    for t in b:
        seq_multiple = []
        for a in db.children(t, featuretype='mRNA', order_by='start'):
            seq_combined = ''
            j = 0
            for i in db.children(a, featuretype='CDS', order_by='start'):
                j += 1
                if j == 1:
                    pphase = i[7]
                seq = i.sequence(myFasta, use_strand=False)
                seq_combined += seq
            seq_combined = SeqRecord(Seq(seq_combined, generic_dna))
            if t.strand == '-':
                pphase = i[7]
                seq_combined = seq_combined.reverse_complement()
            if pphase == "0" or pphase == ".":
                seq_transl = seq_combined.translate()
            elif pphase == "1":
                seq_transl = seq_combined[1:].translate()
            elif pphase == "2":
                seq_transl = seq_combined[2:].translate()
            seq_transl.id = a.attributes["ID"][0]
            if seq_multiple:
                seq_multiple_len = len(seq_multiple)
                c = 0
                while seq_multiple_len > c:
                    a1 = str(seq_multiple[c].seq)
                    a2 =      str(seq_transl.seq)
                    if a1.rstrip("*").startswith(a2.rstrip("*")) and len(a1) < len(a2):
                        seq_multiple[c] = seq_transl
                    elif  a2.rstrip("*").startswith(a1.rstrip("*")) and len(a1) < len(a2):
                        seq_multiple[c] = seq_transl
                    else:
                        seq_multiple = seq_multiple + [seq_transl]
                    c+=1
            else:
                seq_multiple = [seq_transl]

        mrna_select = mrna_select + seq_multiple
    mrna_select_name = []
    for seq in mrna_select:
        mrna_select_name.append(seq.id)
    mrna_total = sorted(mrna_retain + mrna_select_name)
    gff_file = final_evm + ".final.gff3"
    gff_out = gffwriter.GFFWriter(gff_file)
    for key in mrna_total:
        for i in db.parents(key, featuretype='gene', order_by='start'):
            gff_out.write_rec(i)
        gff_out.write_rec(key)
        for i in db.children(key, featuretype='CDS', order_by='start'):
            gff_out.write_rec(i)
        for i in db.children(key, featuretype='exon', order_by='start'):
            gff_out.write_rec(i)
    return(gff_file)

def remove_redudant(genome, final_evm):
    fasta = pep_seq(genome, final_evm)
    gff_new = fasta + ".gffread.gff3"
    cm = "gffread -M -H -V -g %s -K %s -F -o %s" % (genome, fasta, gff_new)
    gffread = sb.Popen(cm, shell=True, stdout=sb.PIPE, stderr=sb.PIPE)
    gffread.communicate()
    gff_file = gff_filter(gff_new, genome)
