import pandas as pd
import os.path
import subprocess
import glob
import truth_sets_aaf as tsa

def prepare4prec_recall(region='chr22', vartype='snp'):
    if vartype == 'snp':
        alt_vartype = 'snvs'
    elif vartype == 'indel':
        alt_vartype = 'indels'
    # paths
    outmaindir = '/home/attila/projects/bsm/results/2019-05-02-make-truth-sets/'
    outdir = outmaindir + region + os.path.sep + vartype + os.path.sep
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    csetdir = '/home/attila/projects/bsm/results/calls/benchmark-mix1-mix3/' + alt_vartype
    csetVCFpatt = csetdir + os.path.sep + '*.vcf.gz'
    csetVCFs = glob.glob(csetVCFpatt)
    # helper function
    def helper(csetVCF):
        if region == 'autosomes':
            args1 = ['prepare4prec-recall', '-v', vartype, '-Oz', '-P', csetVCF]
        elif region == 'chr22':
            args1 = ['prepare4prec-recall', '-r22', '-v', vartype, '-Oz', '-P', csetVCF]
        outVCF = outdir + os.path.basename(csetVCF)
        args2 = ['bcftools', 'view', '-Oz', '-o', outVCF]
        proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE)
        proc2 = subprocess.Popen(args2, shell=False, stdout=subprocess.PIPE, stdin=proc1.stdout)
        proc3 = subprocess.run(['bcftools', 'index', '-t', outVCF])
        return(outVCF)
    val = [helper(y) for y in csetVCFs]
    return(val)


def remove_nonvariants_from_callset(region='chr22', vartype='snp', lam='0.04',
        log10s2g='-2', sample='mix1'):
    maindir = '/home/attila/projects/bsm/results/2019-05-02-make-truth-sets/'
    csetVCF = maindir + region + os.path.sep + vartype + os.path.sep + 'strelka2Somatic.vcf.gz'
    truthdir = '/home/attila/projects/bsm/results/2019-03-18-truth-sets/'
    discardedVCF = truthdir + region + os.path.sep + vartype + \
            '/truthset/aaf/exp_model/lambda_' + \
            lam + '/log10s2g_' + log10s2g + os.path.sep + sample + '/discarded-complete.vcf.gz'
    args = ['bcftools', 'isec', '-C', '-Oz', '-p', '/home/attila/Desktop/dir', csetVCF, discardedVCF]
    subprocess.run(args)
    return(args)


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
