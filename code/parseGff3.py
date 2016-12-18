#!/usr/bin/env python

from sys import argv
import subprocess
import re
import sys
import os.path


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
#Chr8    VDAG_Jr2_v4.0.NCBIorder.fasta_GMAPindex gene    3358883 3360550 .       +       .       ID=gene11025;Name=Gene12296_evm.model.Chr8.941
#Chr8    VDAG_Jr2_v4.0.NCBIorder.fasta_GMAPindex mRNA    3358883 3360550 .       +       .       ID=mRNA11134;Parent=gene11025;Name=Gene12296_evm.model.Chr8.941
        
def main():
    '''Main body of the function'''
    gff_filename = os.path.abspath(sys.argv[1])

    prefix = sys.argv[2]
    genename(gff_filename, prefix)
    print '\n\n\n################\n####FINISHED####\n################\n\n'


if __name__ == '__main__':
    main()
