#! /usr/bin/env python3

import gffutils


def transform(d):
    try:
        d['Parent'] = d['ID']
    except KeyError:
        pass
    return d

def score():
    #args = setting()
    evm_inputs = {}
    evm_inputs["lorean"] = "species_LoReAn.annotation.gff3"
    evm_inputs["augustus"] = "augustus.gff3"
    evm_inputs["genemark"] = "genemark.gtf.gff3"
    evm_inputs["gmap"] = "gmap.trinity.gff3"
    evm_inputs["evm"] = "evm.out.combined.gff3"

    coord_dict = {}

    for gff3 in evm_inputs:
        db = gffutils.create_db(evm_inputs[gff3], dbfn=':memory:', force=True ,keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        list_mRNA = [intron for intron in db.features_of_type('mRNA')]
        for line in list_mRNA:
            start_end = [[cds.chrom, str(cds.start), str(cds.end)] for cds in db.children(line, featuretype='CDS')]
            flat_down = [single for multiple in start_end for single in multiple]
            coords = "_".join(flat_down)
            if coords in coord_dict:
                coord_dict[coords].append(line.source + "-" + line.id)
            else:
                coord_dict[coords] = [line.source + "-" + line.id]

    lorean_names = [mrna.split("-")[1] for key in coord_dict for mrna in coord_dict[key] if "LoReAn" in mrna]

    db = gffutils.create_db(evm_inputs[gff3], dbfn=':memory:', force=True ,keep_order=True, merge_strategy='merge', sort_attribute_values=True)








    return

if __name__ == '__main__':
    score()