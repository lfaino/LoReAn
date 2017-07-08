#!/bin/python

#
#    Convert SAM to indexed, sorted BAM file with headers
#
#    Copyright (C) 2015  Michael Hamilton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


import os, sys, subprocess
import time
from argparse import ArgumentParser, ArgumentTypeError
FNULL = open(os.devnull, 'w')

def writeStatus(status):
    sys.stderr.write('<%s> %s\n' % (time.asctime(),status))

def convert(samfile, bamfile, verbose, processors, memory):
    start = time.time()
    # Convert to BAM
    if verbose:
        writeStatus('Converting %s to %s' % (samfile,bamfile))
    cmd = 'samtools view -Sbh %s > %s ' %(samfile, bamfile)
    status_convert = subprocess.call(cmd,shell=True,stderr=FNULL)
                                     

    # Sort
    if verbose:
        writeStatus('Sorting BAM file')
    sorted_basename = "%s.sorted" %(bamfile.split(".bam")[0])
    cmd = "samtools sort -@ %d -m %dG %s %s" %(processors, memory,  bamfile, sorted_basename)
    status_sort = subprocess.call(cmd,shell=True)
    
    # mv rename sorted file
    cmd = 'mv %s.bam %s' % ( sorted_basename, bamfile )
    status_rename = subprocess.call(cmd,shell=True)

    # Index
    if verbose:
        writeStatus('Generating BAM index')
    cmd = "samtools index %s" % (bamfile)
    status_index = subprocess.call(cmd,shell=True)

    end = time.time()
    if verbose:
        writeStatus( 'Elapsed time: %.2f minutes.' %((end - start)/60.0))

def main():
    parser = ArgumentParser(description='Generate sorted BAM and index files for given SAM file')
    parser.add_argument('samfile', type=str, help='Samfile to convert')
    parser.add_argument('-o', '--outfile', dest='bamfile', 
                      help='Name of converted BAM file [default=<sambase>.bam]')
    parser.add_argument('-p', '--procs', dest='procs', type=int, default=1,
                      help='Number of processors to use for BAM sorting (default 1)' )
    parser.add_argument('-m', '--memory', dest='memory', type=int, default=2,
                      help='Max memory (in GBs) for each processor used for BAM sorting (default 2)') 


    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default='False',
                      help='Print verbose output')

    args = parser.parse_args()
    path = os.path.split(args.samfile)[0]
    if not args.bamfile:
        args.bamfile = os.path.join(path, os.path.basename(args.samfile).split('.')[0] + '.bam')
    path = os.path.split(args.samfile)[0]
    convert(args.samfile, args.bamfile, args.verbose, args.procs, args.memory)

if __name__ == "__main__":
    main()
    

    
