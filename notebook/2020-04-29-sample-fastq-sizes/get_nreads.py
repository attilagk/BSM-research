#! /usr/bin/env python3
import subprocess
import pandas as pd
import numpy as np
import re
import io
import os.path

testbams = ['/projects/bsm/alignments/MSSM_118/MSSM_118_NeuN_pl.bam',
        '/projects/bsm/alignments/MSSM_118/MSSM_118_muscle.bam']

def nreads_from_bam(bam=testbams[0]):
    args = ['samtools', 'idxstats', bam]
    p = subprocess.run(args, capture_output=True)
    df = pd.read_csv(io.BytesIO(p.stdout), sep='\s+', header=None)
    nreads = np.sum(np.array(df.iloc[:, 1:]).ravel())
    sample = os.path.basename(bam).replace('.bam', '')
    d = {'sample': sample, 'nreads': nreads, 'path': bam}
    df = pd.DataFrame(d, index=[sample])
    return(df)

def nreads_from_bamlist(bamlist=testbams):
    l = [nreads_from_bam(bam) for bam in bamlist]
    df = pd.concat(l, axis=0)
    outdir = '/projects/bsm/attila/results/2020-04-29-sample-fastq-sizes/'
    outfile = outdir + os.path.sep + 'nreads.csv'
    df.to_csv(outfile, sep='\t', index=False, header=True)
    return(df)

if __name__ == "__main__":
        import sys
        bams = sys.argv[1:]
        nreads_from_bamlist(bams)
