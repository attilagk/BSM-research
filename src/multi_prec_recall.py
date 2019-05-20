import pandas as pd
import os.path
import tempfile
import shutil
import subprocess
import glob
import truth_sets_aaf as tsa

__callsetmaindir__ = '/home/attila/projects/bsm/results/calls/benchmark-mix1-mix3/'
__truthsetmaindir__ = '/home/attila/projects/bsm/results/2019-03-18-truth-sets/'
__outmaindir__ = '/home/attila/projects/bsm/results/2019-05-02-make-truth-sets/'
__expmsubdir__ = 'truthset/aaf/exp_model/lambda_'


def getVCFpaths(callsetbn=None, region='chr22', vartype='snp', lam='0.04',
        log10s2g='-2', sample='mix1'):
    '''
    Create pathname for various input, output, intermediate VCF files.

    Parameters:
    callsetbn: None, 'Tnseq.vcf.gz' or ['lofreqSomatic.vcf.gz', 'Tnseq.vcf.gz']
    region: chr22 or autosomes
    vartype: snp or indel
    lam: '0.04' or '0.2'
    log10s2g: '-2', '-3' or '-4'
    sample: 'mix1', 'mix2' or 'mix3'

    Returns:
    
    a dictionary of pathnames

    The keys of the dictionary are as follows:

    callset: the original callset without filtering

    prepared_callset_dir: directory of VCFs filtered for region and vartype
    
    prepared_callset: filtered for region and vartype
    
    reduced_truthset: the truthset according to the exp_model with parameters
    lam, log10s2g, sample
    
    discarded_truthset: the complementer set of reduced_truthset relative to
    the original truthset based on the original mixes

    reduced_callset: created from reduced_truthset by removing the
    "nonvariants" of discarded_truthset

    Details:

    If callsetbn is None then a list of all pathnames are used that match the
    pattern '*.vcf.gz'.  If callsetbn is a list of file basenames then those
    will be extended into pathnames.  If callsetbn is a string of a single
    basename then that will be extended into a single pathname.
    '''
    # some pieces of pathnames
    subdir1 = region + os.path.sep + vartype + os.path.sep
    subdir2 = lam + '/log10s2g_' + log10s2g + os.path.sep + sample + os.path.sep
    if vartype == 'snp':
        alt_vartype = 'snvs'
    elif vartype == 'indel':
        alt_vartype = 'indels'
    # directories
    truthsetdir = __truthsetmaindir__ + subdir1 + __expmsubdir__ + subdir2
    callsetdir = __callsetmaindir__ + alt_vartype + os.path.sep
    prepared_callset_dir = __outmaindir__ + subdir1
    reducedcsetdir = __outmaindir__ + subdir1 + __expmsubdir__ + subdir2
    # filepaths
    reduced_truthset =  truthsetdir + 'complete.vcf.gz'
    discarded_from_truthset =  truthsetdir + 'discarded-complete.vcf.gz'
    if callsetbn is None:
        callset = glob.glob(callsetdir + '*.vcf.gz')
        prepared_callset = [prepared_callset_dir + os.path.basename(y) for y in callset]
        reduced_callset = [reducedcsetdir + os.path.basename(y) for y in callset]
    elif isinstance(callsetbn, list):
        callset = [callsetdir + y for y in callsetbn]
        prepared_callset = [prepared_callset_dir + y for y in callsetbn]
        reduced_callset = [reducedcsetdir + y for y in callsetbn]
    else:
        callset = callsetdir + callsetbn
        prepared_callset = prepared_callset_dir + callsetbn
        reduced_callset = reducedcsetdir + callsetbn
    VCFpaths = {'reduced_truthset': reduced_truthset,
            'discarded_from_truthset': discarded_from_truthset, 'callset':
            callset, 'prepared_callset_dir': prepared_callset_dir,
            'prepared_callset': prepared_callset, 'reduced_callset':
            reduced_callset}
    return(VCFpaths)


def prepare4prec_recall(region='chr22', vartype='snp'):
    '''
    Run prepare4prec-recall shell script on initial callset VCFs for a given
    region and variant type

    Parameter(s):
    region: chr22 or autosomes
    vartype: snp or indel

    Returns: a list of the pathname of output VCFs
    '''
    if vartype == 'snp':
        callers = ['lofreqSomatic', 'somaticSniper', 'strelka2Germline2s', 'strelka2Somatic', 'Tnseq']
    elif vartype == 'indel':
        callers = ['strelka2Germline2s', 'strelka2Somatic', 'Tnseq']
    callsetbn = [c + '.vcf.gz' for c in callers]
    VCFpaths = getVCFpaths(callsetbn=callsetbn, region=region, vartype=vartype)
    outdir = VCFpaths['prepared_callset_dir']
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    csetVCFs = VCFpaths['callset']
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


def reduce_prepared_truthset(region='chr22', vartype='snp', lam='0.04',
        log10s2g='-2', sample='mix1'):
    maindir = '/home/attila/projects/bsm/results/2019-05-02-make-truth-sets/'
    csetVCF = maindir + region + os.path.sep + vartype + os.path.sep + 'strelka2Somatic.vcf.gz'
    discardedVCF = __truthsetmaindir__ + region + os.path.sep + vartype + \
            __expmsubdir__ + \
            lam + '/log10s2g_' + log10s2g + os.path.sep + sample + '/discarded-complete.vcf.gz'
    tempdir = tempfile.TemporaryDirectory()
    tempVCF = tempdir.name + os.path.sep + '0000.vcf.gz'
    tempVCFtbi = tempVCF + '.tbi'
    outdir = '/home/attila/Desktop/'
    outVCF = outdir + os.path.sep + os.path.basename(csetVCF)
    outVCFtbi = outVCF + '.tbi'
    args1 = ['bcftools', 'index', '--tbi', discardedVCF]
    subprocess.run(args1)
    args2 = ['bcftools', 'isec', '-C', '-Oz', '-p', tempdir.name, csetVCF, discardedVCF]
    subprocess.run(args2)
    shutil.move(tempVCF, outVCF)
    shutil.move(tempVCFtbi, outVCFtbi)
    return(args2)


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
