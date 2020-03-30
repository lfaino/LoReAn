#! /usr/bin/env python3


###############
###IMPORTS###
###############

# LIBRARIES
import sys
import test as test

# OTHER SCRIPTS
import collectOnly as collect
import dirsAndFiles as logistic
import getRightStrand as grs
import pasa as pasa


def main(ref_rename,):

    merged_gff3 = collect.add_EVM(final_output, gmap_wd, consensus_mapped_gff3)
    #print(merged_gff3)
    update2 = grs.exonerate(ref_rename, merged_gff3, threads_use, exonerate_wd, args.verbose)
    print(ref_rename, update2)
    update3_1 = test.remove_redudant(ref_rename, update2)
    print(update3_1)
    update3 = grs.genename_lorean(update3_1, args.verbose, exonerate_wd)
    print(update3)
    # HERE WE COMBINE TRINITY OUTPUT AND THE ASSEMBLY OUTPUT TO RUN AGAIN
    # PASA TO CORRECT SMALL ERRORS
    sys.stdout.write(("###FIXING GENES NON STARTING WITH MET\t" + now + "\t###\n"))
    fasta_all = logistic.cat_two_fasta(trinity_out, tmp_assembly_all, long_fasta, pasa_dir)
    round_n += 1
    update5 = pasa.update_database(threads_use, str(round_n), pasa_dir, pasadb,  ref_rename, fasta_all,
                                   update3, args.verbose)

if __name__ == '__main__':
    main(*sys.argv[1:])
    # remove_redudant("/home/lfaino/lfainoData/lorean/LoReAn_Example/Crispa/LoReAn_crispa/run/split/scaffold3.fasta.masked.rename.fasta",
    #                "/home/lfaino/lfainoData/lorean/LoReAn_Example/Crispa/LoReAn_crispa/run/exonerate/genename_lorean.4.k9w51gwe.gff3")
