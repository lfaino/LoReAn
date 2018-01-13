#!/usr/bin/env python3

import subprocess
import sys
import tempfile

#======================================================================================================================

GT_GFF3 = 'gt gff3 -sort -tidy -setsource external %s'

#======================================================================================================================


def external(external, wd, verbose):

    error = tempfile.NamedTemporaryFile(delete=False, mode = "w", prefix="external.", suffix= ".error3")
    tmp_file = tempfile.NamedTemporaryFile(delete=False, mode = "w", prefix="external.", suffix= ".gff3")
    cmd = GT_GFF3 % external
    print(cmd)
    if verbose:
        sys.stderr.write('Executing: %s\n\n' % cmd)
        sys.stderr.write('Log file is: %s\n\n' % error.name)
    gt_call = subprocess.Popen(cmd, stdout=tmp_file, stderr=error, cwd=wd, shell=True)
    gt_call.communicate()

    return tmp_file.name


