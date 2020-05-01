#! /usr/bin/env python3
import subprocess
import pandas as pd
import numpy as np
import re
import io
import os.path
import glob

testbams = ['/projects/bsm/alignments/MSSM_118/MSSM_118_NeuN_pl.bam',
        '/projects/bsm/alignments/MSSM_118/MSSM_118_muscle.bam']

def nreads_from_bam(bam=testbams[0]):
    args = ['samtools', 'idxstats', bam]
    p = subprocess.run(args, capture_output=True)
    idxstats = io.BytesIO(p.stdout)
    df = nreads_from_idxstats(idxstats, bam)
    return(df)

def nreads_from_idxstats(idxstats, bam=None):
    df = pd.read_csv(idxstats, sep='\s+', header=None)
    nreads = np.sum(np.array(df.iloc[:, 2:]).ravel())
    if bam is None:
        ext = '.idxstats'
        fpath = idxstats
    else:
        ext = '.bam'
        fpath = bam
    sample = os.path.basename(fpath).replace(ext, '')
    d = {'sample': sample, 'nreads': nreads, 'path': fpath}
    df = pd.DataFrame(d, index=[sample])
    return(df)

def nreads_from_bamlist(bamlist=testbams):
    l = [nreads_from_bam(bam) for bam in bamlist]
    df = pd.concat(l, axis=0)
    outdir = '/projects/bsm/attila/results/2020-04-29-sample-fastq-sizes/'
    outfile = outdir + os.path.sep + 'nreads.tsv'
    df.to_csv(outfile, sep='\t', index=False, header=True)
    return(df)

def nreads_from_idxstatslist():
    outdir='/home/attila/projects/bsm/results/2020-04-29-sample-fastq-sizes/'
    idxstatslist = glob.glob(outdir + os.path.sep + 'GENEWIZ-idxstats/*.idxstats')
    l = [nreads_from_idxstats(i, bam=None) for i in idxstatslist]
    df = pd.concat(l, axis=0)
    outfile = outdir + os.path.sep + 'nreads-genewiz.tsv'
    df.to_csv(outfile, sep='\t', index=False, header=True)
    return(df)

if __name__ == "__main__":
        import sys
        bams = sys.argv[1:]
        nreads_from_bamlist(bams)
