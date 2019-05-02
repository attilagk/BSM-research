import sys
import os
import subprocess
import pandas as pd
import numpy as np

def sample_positions(ssize, invcfpath, outvcfpath, seed=19760415):
    args0 = ['bcftools', 'view', '-H', invcfpath]
    args1 = ['cut', '-f1,2']
    proc0 = subprocess.Popen(args0, shell=False, stdout=subprocess.PIPE)
    proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE,
            stdin=proc0.stdout)
    regions = pd.read_csv(proc1.stdout, sep='\t', names=['CHROM', 'POS'],
            dtype={'CHROM': 'category', 'POS': np.int64})
    # downsampling regions
    regions = regions.sample(n=ssize, replace=False, axis=0, random_state=seed)
    regions = regions.sort_index()
    # output: regions file
    outdirpath = os.path.dirname(outvcfpath)
    if not os.path.exists(outdirpath):
        os.makedirs(outdirpath)
    regionspath = outvcfpath + '.regions'
    regions.to_csv(regionspath, sep='\t', header=False, index=False)
    # output: VCF
    args2 = ['bcftools', 'view', '-R', regionspath, '-Oz', '-o', outvcfpath,
            invcfpath]
    subprocess.run(args2)# and os.unlink(regionspath)
    return(regions)
