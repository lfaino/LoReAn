#!/usr/bin/env python3

"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

Porechop makes use of C++ functions which are compiled in cpp_functions.so. This module uses ctypes
to wrap them in similarly named Python functions.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ctypes import CDLL, cast, c_char_p, c_int, c_void_p
from multiprocessing.dummy import Pool as ThreadPool

SO_FILE = 'cpp_functions.so'
SO_FILE_FULL = os.path.join(os.path.dirname(os.path.realpath(__file__)), SO_FILE)
if not os.path.isfile(SO_FILE_FULL):
    sys.exit('could not find ' + SO_FILE + ' - please reinstall')
C_LIB = CDLL(SO_FILE_FULL)

C_LIB.adapterAlignment.argtypes = [c_char_p,  # Read sequence
                                   c_char_p,  # Adapter sequence
                                   c_int,     # Match score
                                   c_int,     # Mismatch score
                                   c_int,     # Gap open score
                                   c_int]     # Gap extension score
C_LIB.adapterAlignment.restype = c_void_p     # String describing alignment


# This function cleans up the heap memory for the C strings returned by the other C functions. It
# must be called after them.
C_LIB.freeCString.argtypes = [c_void_p]
C_LIB.freeCString.restype = None


def adapter_alignment(read_sequence, adapter_sequence, scoring_scheme_vals, alignm_score_out):
    """
    Python wrapper for adapterAlignment C++ function.
    """
    list_adapter = []
    list_run = []
    for adapter in SeqIO.parse(adapter_sequence, "fasta"):
        list_adapter.append(adapter)
        record = SeqRecord(adapter.seq.reverse_complement(), id=adapter.id + "_rev")
        list_adapter.append(record)
    dict_aln = {}
    for sequence in SeqIO.parse(read_sequence, "fastq"):
        dict_aln[sequence.id] = ""
        for adapter in list_adapter:
            match_score = scoring_scheme_vals[0]
            mismatch_score = scoring_scheme_vals[1]
            gap_open_score = scoring_scheme_vals[2]
            gap_extend_score = scoring_scheme_vals[3]
            list_run.append([str(sequence.seq).encode('utf-8'), str(adapter.seq).encode('utf-8'), match_score,
                             mismatch_score, gap_open_score, gap_extend_score, sequence.id, adapter.id])
    print ("start aln")

    with ThreadPool(8) as pool:
        for out in pool.imap(align, list_run):
            out_list = out.split(",")
            if dict_aln[out_list[0]] != "":
                if (float(out.split(",")[9])) > float(dict_aln[out_list[0]].split(",")[9]):
                    dict_aln[out.split(",")[0]] = out
            else:
                dict_aln[out.split(",")[0]] = out

    if alignm_score_out == "":
        alignm_score_mean = np.mean([float(dict_aln[key].split(",")[9]) for key in dict_aln])
        alignm_score_std = np.std([float(dict_aln[key].split(",")[9]) for key in dict_aln])
        alignm_score_out = alignm_score_mean - alignm_score_std

    seq_to_keep = {}
    for key in dict_aln:
        if (float(dict_aln[key].split(",")[9])) > alignm_score_out:
            seq_to_keep[key] = dict_aln[key]
    with open("./example.fasta", "w") as output_handle:
        for sequence in SeqIO.parse(read_sequence, "fastq"):
            if sequence.id in seq_to_keep:
                if seq_to_keep[sequence.id].split(",")[1].endswith("rev"):
                    rev_seq = SeqRecord(sequence.seq.reverse_complement(), id=sequence.id + "_rev")
                    SeqIO.write(rev_seq, output_handle, "fasta")
                else:
                    SeqIO.write(sequence, output_handle, "fasta")
    return


def align(command_in):
    ptr = C_LIB.adapterAlignment(str(command_in[0]).encode('utf-8'), str(command_in[1]).encode('utf-8'),
                                 command_in[2], command_in[3], command_in[4], command_in[5])
    result_string = c_string_to_python_string(ptr)
    single_result_string = result_string.split(",")
    average_score = (float(single_result_string[5]) + float(single_result_string[6])) / 2
    result_string_name = ",".join([command_in[6], command_in[7], result_string, str(average_score)])
    return result_string_name


def c_string_to_python_string(c_string):
    """
    This function casts a C string to a Python string and then calls a function to delete the C
    string from the heap.
    """
    python_string = cast(c_string, c_char_p).value.decode()
    C_LIB.freeCString(c_string)
    return python_string


if __name__ == '__main__':
    scoring = [3, -6, -5, -2]
    alignm_score_out = ""
    adapter_alignment(*sys.argv[1:], scoring, alignm_score_out)
