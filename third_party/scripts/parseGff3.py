#!/usr/bin/env python
from gffutils.iterators import DataIterator
import gffutils
from sys import argv
import subprocess
import re
import sys

#def parseGff(input_filename):

input_filename = sys.argv[1]


dictTable = {}
listGene = []


dbName = input_filename + '.db'
db = gffutils.create_db(input_filename, dbfn= dbName, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

p1 = subprocess.Popen(['grep', 'mRNA' , input_filename ,  ], stdout = subprocess.PIPE)
p2 = subprocess.Popen(['awk', r'{print $0"-"$5-$4}'], stdin = p1.stdout, stdout=subprocess.PIPE)
p3 = subprocess.Popen(['bedtools' , 'sort', '-i' ], stdout=subprocess.PIPE, stdin = p2.stdout)
p4 = subprocess.Popen(['bedtools merge -d -200 -s -c  9,9 -o  count,distinct -i - '], stdout = subprocess.PIPE, stdin = p3.stdout, shell = True)
c = 0
for line in p4.stdout:
    listGene = line.split('\t')
    c += 1
    indexName = "seq" + str(c)
    parts = re.split(',|;', listGene[-1])
    for elements in parts:
        if 'Parent' in elements:
            newEl = elements.replace("Parent=", "")
            if dictTable.has_key(indexName):
                cont = dictTable[indexName]
                oldLine = cont.split('-')
                newLine = newEl.split('-')
                if int(oldLine[1]) > (newLine[1]):
                    pass
                else:
                    
                    dictTable[indexName] = newEl.rstrip('\n')
            else:
                
                dictTable[indexName] = newEl.rstrip('\n')

        
for key in dictTable:
    ext = dictTable[key]
    value = ext.split('-')
    db = gffutils.FeatureDB(dbName, keep_order=True)
    gene = db[value[0]]
    
    
    
    print gene
    for i in db.children(gene, featuretype='mRNA', order_by='start'):
        print(i)
    for i in db.children(gene, featuretype='exon', order_by='start'):
        print(i) 
    for i in db.children(gene, featuretype='CDS', order_by='start'):
        print(i) 
    
    
#if __name__ == '__main__':
    
    #input_filename = argv[1]
    #parseGff(input_filename)
