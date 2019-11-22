#! /usr/bin/env python3

import os
import os.path
import subprocess
import glob

def split_bam(bam='/projects/bsm/alignments/PITT_101/PITT_101_NeuN_pl-22-1Mb.bam'):
    maindir = os.path.dirname(bam)
    bname = os.path.basename(bam).replace('.bam', '')
    subdir = maindir + os.path.sep + bname
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    os.chdir(subdir)
    args = ['samtools', 'split', bam]
    proc = subprocess.run(args)
    splitbams = glob.glob(subdir + os.path.sep + '*.bam')
    return(splitbams)

def replace_rg(bam):
    args = ['samtools', 'view', '-H', bam]
    proc = subprocess.run(args, capture_output=True)
    return(proc)
