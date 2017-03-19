#!/usr/bin/env python3
''
import os
import subprocess
from sys import argv


def AAT(proteinFastaFile, ref, wd):
    '''Calls genome guided trinity on the BAM file to generate
    assembled transcripts'''
    wd_output = wd + '/protein_evidence.gff3'
    if os.path.isfile(wd_output):
        print(('AAT files exist: ' + wd_output + ' --- skipping\n'))
        return True
    else:
        args = [
            'AAT.pl',
            '-P',
            '-b',
            '-q',
            ref,
            '-s',
            proteinFastaFile,
            r"--dps",
            r"'-f 100 -i 30 -a 200'",
            r"--filter",
            r"'-c 10'",
            r"--nap",
            r"'-x 10'"]
        chromo = ref.split('/')[-1]
        out_name = wd + chromo + '.protein_evidence.gff3'

        log_name = wd + 'AAT.log'
        log = open(log_name, 'w')
        stdout_f = open(wd + 'AAT.stdout', 'w')
        aat_process = subprocess.Popen(
            args, stderr=log, stdout=stdout_f, cwd=wd)
        aat_process.wait()
        log.close()
        stdout_f.close()
        return True


def parseAAT(wd):
    '''all the protein alignemnt files are concatenated in one file. this is becasue the AAT is paralelized'''
    outFilename = wd + '/protein_evidence.out'
    wd_gff = ['cat']
    o_file = open(outFilename, 'w')
    for root, dirs, files in os.walk(wd):
        for name in files:
            if 'btab' in name:
                    wd_gff.append(os.path.join(root, name))
    cat_call = subprocess.Popen(wd_gff, stdout=o_file)
    cat_call.communicate()
    o_file.close()

    outFilenameGff = wd + '/protein_evidence.gff3'
    args_btab = ['AAT_btab_to_gff3.pl', outFilename, 'P', ]
    stdout_file = open(outFilenameGff, 'w')
    subprocess.Popen(args_btab, stdout=stdout_file, cwd=wd)
    stdout_file.close()
    return outFilenameGff


def main():
    protFasta = argv[1]
    ref = argv[2]
    wd = './'
    AAT(protFasta, ref, wd)
    mergedProtGFF3 = parseAAT(wd)
#    print mergedProtGFF3


if __name__ == '__main__':
    main()
