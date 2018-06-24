#! /usr/bin/env python3


import sys

import collectOnly as collect
import getRightStrand as grs


# OTHER SCRIPTS


def main(ref_rename, final_output, consensus_mapped_gff3, wd):
    threads_use = 2
    verbose = True
    merged_gff3 = collect.add_EVM(final_output, wd, consensus_mapped_gff3)
    update2 = grs.exonerate(ref_rename, merged_gff3, threads_use, wd, verbose)
    update3 = grs.genename_lorean(update2, verbose, wd)

if __name__ == '__main__':
    main(*sys.argv[1:])