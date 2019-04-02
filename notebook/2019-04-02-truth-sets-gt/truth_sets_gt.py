import pandas as pd
import numpy as np

def read_aaf_of_gt(csvpath='/home/attila/projects/bsm/results/2019-03-12-prec-recall-design/aaf_of_gt.csv'):
    names =  ['genotype', 'mix1', 'mix2', 'mix3']
    dtype = {'genotype':np.object}
    aaf_of_gt = pd.read_csv(csvpath, names=names, dtype=dtype, header=0)
    return(aaf_of_gt)

def helper(gt):
    vcfbns = [(None, '1xxx', '2xxx'),
            (None, 'x1xx', 'x2xx'),
            (None, 'xx1x', 'xx2x'),
            (None, 'xxx1', 'xxx2')]
    res = [y[0][int(y[1])] for y in zip(vcfbns, gt)]
    res = [y for y in res if not y is None]
    return(res)
