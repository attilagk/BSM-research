import pandas as pd
import numpy as np
import subprocess
import os
import re

def read_aaf_of_gt(csvpath='/home/attila/projects/bsm/results/2019-03-12-prec-recall-design/aaf_of_gt.csv'):
    names =  ['genotype', 'mix1', 'mix2', 'mix3']
    dtype = {'genotype':np.object}
    aaf_of_gt = pd.read_csv(csvpath, names=names, dtype=dtype, header=0)
    return(aaf_of_gt)

def bcftools_isec_args(gt,
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
    outfpath = outdir + os.sep + gt
    nstr = re.sub('2', '1', gt)
    pdirpath = outdir + os.sep + gt
    # arguments for bcftools isec
    bcftool0 = ['bcftools', 'isec']
    p_opt = ['-p', pdirpath]
    c_opt = ['-c', 'none']
    n_opt = ['-n~' + nstr]
    args0 = bcftool0 + p_opt + c_opt + n_opt + flist
    # arguments for copying VCF to keep
    bcftool1 = ['bcftools', 'view']
    O_opt = ['-Oz']
    o_opt = ['-o', pdirpath + '.vcf.gz']
    vcf2keep = [pdirpath + os.sep + '000' + str(nstr.find('1')) + '.vcf']
    args1 = bcftool1 + O_opt + o_opt + vcf2keep
    return((args0, args1))

# bcftools isec -n~1100 -c all A.vcf.gz B.vcf.gz C.vcf.gz D.vcf.gz
