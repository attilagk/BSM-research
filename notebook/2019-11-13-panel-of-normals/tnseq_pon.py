#! /usr/bin/env python3

import subprocess
import os
import os.path
import subprocess

resdir = '/projects/bsm/attila/results/2019-11-13-panel-of-normals'
NUMBER_THREADS = 16
REFSEQ = '/projects/shared/refgenome/GRCh37/dna/hs37d5.fa'
sentieon_executable = '/opt/sentieon-genomics/current/bin/sentieon'

def start_licsrvr():
    #TODO
    pass

def bam2pon(bam='/projects/bsm/alignments/MSSM_065/MSSM_106_muscle.bam',
        nthreads=NUMBER_THREADS, refseq=REFSEQ, algo='TNhaplotyper'):
    nthreads = str(nthreads)
    vcf = os.path.basename(bam).replace('.bam', '.vcf')
    outdir = resdir + os.path.sep + 'VCFs'
    vcf = outdir + os.path.sep + vcf
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    args = [sentieon_executable, 'driver', '-t', nthreads, '-r', refseq, '-i', bam,
            '--algo', algo, '--detect_pon', '--interval', '22', vcf]
            #'--algo', algo, '--detect_pon', vcf]
    proc = subprocess.run(args)
    return(proc)
    #return(args)

# sentieon driver -t $NUMBER_THREADS -r $REFSEQ -i $BAM --algo TNhaplotyper --detect_pon $VCF
