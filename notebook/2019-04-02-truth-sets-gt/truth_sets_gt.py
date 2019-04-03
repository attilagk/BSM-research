#! /usr/bin/env python3
import pandas as pd
import numpy as np
import subprocess
import os
import re
import shutil

def read_aaf_of_gt(csvpath='/home/attila/projects/bsm/results/2019-03-12-prec-recall-design/aaf_of_gt.csv'):
    names =  ['genotype', 'mix1', 'mix2', 'mix3']
    dtype = {'genotype':np.object}
    aaf_of_gt = pd.read_csv(csvpath, names=names, dtype=dtype, header=0)
    return(aaf_of_gt)

def bcftools_isec_gt(gt,
        indir='/home/attila/projects/bsm/results/2019-03-18-truth-sets/chr22/snp/truthset/individuals',
        outdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets/chr22/snp/truthset/genotypes'):
    def make_flist():
        vcfbns = [('1xxx', '1xxx', '2xxx'),
                ('x1xx', 'x1xx', 'x2xx'),
                ('xx1x', 'xx1x', 'xx2x'),
                ('xxx1', 'xxx1', 'xxx2')]
        res = [y[0][int(y[1])] for y in zip(vcfbns, gt)]
        res = [y for y in res if not y is None]
        res = [indir + os.sep + y + '.vcf.gz' for y in res]
        return(res)
    flist = make_flist()
    nstr = re.sub('2', '1', gt)
    pdirpath = outdir + os.sep + gt
    # arguments for bcftools isec
    bcftool0 = ['bcftools', 'isec']
    p_opt = ['-p', pdirpath]
    c_opt = ['-c', 'none']
    n_opt = ['-n~' + nstr]
    args0 = bcftool0 + p_opt + c_opt + n_opt + flist
    # arguments for copying and gzipping the output VCF to keep
    outfpath = pdirpath + '.vcf.gz'
    bcftool1 = ['bcftools', 'view']
    O_opt = ['-Oz']
    o_opt = ['-o', outfpath]
    vcf2keep = [pdirpath + os.sep + '000' + str(nstr.find('1')) + '.vcf']
    args1 = bcftool1 + O_opt + o_opt + vcf2keep
    # arguments for indexing gzipped output VCF
    args2 = ['bcftools', 'index', '-t', outfpath]
    # execute
    subprocess.run(args0)
    subprocess.run(args1)
    subprocess.run(args2)
    shutil.rmtree(pdirpath)
    return((args0, args1, args2))

if __name__ == '__main__':
    import sys
    bcftools_isec_gt(gt=sys.argv[1], indir=sys.argv[2], outdir=sys.argv[3])
