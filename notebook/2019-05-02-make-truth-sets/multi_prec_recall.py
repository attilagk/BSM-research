import pandas as pd
import os.path
import subprocess
import glob
import truth_sets_aaf as tsa

def prec_recall_one_truthset(truthsetVCF, callsetVCFdir='/big/results/bsm/2019-03-22-prec-recall/chr22/snp/'):
    csetVCFs = glob.glob(callsetVCFdir + os.path.sep + '*.vcf.gz')
    #tsetdir = os.path.dirname(truthsetVCF)
    args = ['prec-recall-vcf', '-t', truthsetVCF] + csetVCFs
    proc1 = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    prcsv = pd.read_csv(proc1.stdout)
    return(prcsv)

def prec_recall_all_truthsets(expm, callsetVCFdir='/big/results/bsm/2019-03-22-prec-recall/chr22/snp/'):
    def helper(expm1):
        pathname = tsa.deduce_pathname(expm1)['outdir'] + os.path.sep + 'complete.vcf.gz'
        prcsv = prec_recall_one_truthset(pathname, callsetVCFdir=callsetVCFdir)
        return(prcsv)
    l = tsa.split_up_expm(expm)
    l = [y for y in l if len(y) > 0]
    l = l[0:4] # for testing
    df = [helper(y) for y in l if len(y) > 0]
    return(df)
    pathnames = [tsa.deduce_pathname(y)['outdir'] for y in l if len(y) > 0]
    pathnames = [y + os.path.sep + 'complete.vcf.gz' for y in pathnames]
    y = pathnames[1]
    prcsv = prec_recall_one_truthset(y, callsetVCFdir=callsetVCFdir)
    return(y)
