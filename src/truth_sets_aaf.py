import pandas as pd
import numpy as np
import shutil
import subprocess
import os
import re
import joint_gt_ceph as jgc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

def make_ts_aaf(mix='mix1', vartype='snp', region='chr22',
        tsdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets/chr22/snp/truthset'):
    '''
    Create AAF-based truth set partitions given genotype-based partitions and
    count variants.

    Parameters:
    mix: mix1, mix2 or mix3
    vartype: snp or indel
    region: chr22 or autosomes

    Returns:
    a pandas data frame whose rows are truth subsets at various AAFs and whose
    "nvariants" column contains the number of variants
    '''
    def helper(aaf):
        genotypes = gt_of_aaf[mix][aaf]
        invcfs = [indir + os.sep + g + '.vcf.gz' for g in genotypes]
        outvcf = mixdir + os.sep + str(aaf) + '.vcf.gz'
        unsorted_outvcf = mixdir + os.sep + str(aaf) + '-unsorted.vcf.gz'
        args0 = ['bcftools', 'concat', '-o', unsorted_outvcf, '-Oz'] + invcfs
        args1 = ['bcftools', 'sort', '-o', outvcf, '-Oz', unsorted_outvcf]
        args2 = ['bcftools', 'index', '-t', outvcf]
        if not os.path.isfile(outvcf):
            subprocess.run(args0)
            subprocess.run(args1)
            os.remove(unsorted_outvcf)
            subprocess.run(args2)
        # count records
        args3 = ['bcftools', 'view', '-H', outvcf]
        args4 = ['wc', '-l']
        proc1 = subprocess.Popen(args3, shell=False, stdout=subprocess.PIPE)
        proc2 = subprocess.Popen(args4, shell=False, stdout=subprocess.PIPE, stdin=proc1.stdout)
        proc1.stdout.close()
        nrec = proc2.communicate()[0]
        nrec = int(nrec) # turn bytesliteral (e.g. b'5226\n') to integer
        return(nrec)
    gt_of_aaf = jgc.get_gt_of_aaf()
    indir = tsdir + os.sep + 'genotypes'
    aafdir = tsdir + os.sep + 'aaf'
    if not os.path.isdir(aafdir):
        os.mkdir(aafdir)
    unfdir = aafdir + os.sep + 'unfiltered'
    if not os.path.isdir(unfdir):
        os.mkdir(unfdir)
    mixdir = unfdir + os.sep + mix
    if not os.path.isdir(mixdir):
        os.mkdir(mixdir)
    aafs = gt_of_aaf[mix].keys()
    nrec = [helper(f) for f in aafs]
    res = pd.DataFrame({'region': region, 'vartype': vartype, 'sample': mix, 'aaf': list(aafs),
        'nvariants': nrec})
    res = res.astype({'region': 'category', 'vartype': 'category', 'sample': 'category'})
    return(res)

def make_ts_aaf_get_nvariants():
    '''
    Calls make_ts_aaf for all combination of regions, vartypes and samples

    Returns:
    a pandoc data frame with the number of variants at all those combinations and AAFs
    '''
    def tsaaf(mix, vartype, region):
        bdir = '/home/attila/projects/bsm/results/2019-03-18-truth-sets/'
        res = make_ts_aaf(mix=mix, vartype=vartype, region=region, tsdir=bdir + os.sep + region + os.sep + vartype + '/truthset')
        return(res)
    regions = ['chr22', 'autosomes']
    vartypes = ['snp', 'indel']
    samples = ['mix1', 'mix2', 'mix3']
    nvariants = pd.concat([tsaaf(mix=m, vartype=v, region=r) for r in regions
        for v in vartypes for m in samples])
    return(nvariants)

def evaluate_model(Y, xset, model, only_params):
    '''
    Evaluate a model at a set of x values (AAF values)

    Parameters:
    Y: total number of somatic variants
    xset: the set of x values (AAF values) at which somatic variants exist in the truth set
    model: M1, M2,...
    only_params: return only model parameters a and b

    Returns
    if only_params=False: the set of y values (number of modeled somatic variants)
    if only_params=True: the model parameters a (intercept) and b (slope)
    '''
    if model == 'M1':
        x_max = 50
    elif model == 'M2':
        x_max = 100
    elif model == 'M3':
        x_max = 50
    elif model == 'M4':
        x_max = 100
    xsubset = [z for z in xset if z <= x_max]
    n = len(xsubset)
    xsum = sum(xsubset)
    a = Y / (n - xsum / x_max)
    b = - a / x_max
    if model in ['M3', 'M4']:
        a = Y / n
        b = 0
    if only_params:
        return({'a': a, 'b': b})
    def M_function(x):
        y = a + b * x
        return(y)
    yset = [M_function(z) for z in xsubset] + [0] * (len(xset) - n)
    return(yset)


def evalmodel2df(nvariants, sample, vartype, region, model, Y, p_som2germ=None):
    '''
    Calls evaluate_model for a combination of Y (total n.o. variants), sample,
    vartype, region, model

    Parameters:
    nvariants: a data frame returned by make_ts_aaf_get_nvariants
    Y: total number of somatic variants
    sample: mix1, mix2 or mix3
    vartype: snp or indel
    region: chr22 or autosomes
    model: M1, M2,...
    p_som2germ: -log_10 of the ratio of the number of somatic to germline variants

    Returns:
    a pandas data frame with the modeled number of variants
    '''
    xset = nvariants.loc[(nvariants['region'] == region) &
            (nvariants['sample'] == sample) &
            (nvariants['vartype'] == vartype), 'aaf']
    y = evaluate_model(Y=Y, xset=xset, model=model, only_params=False)
    df = pd.DataFrame({'y': y, 'x': xset, 'model': model, 'Y': int(Y),
        'p_som2germ': p_som2germ, 'sample': sample, 'vartype': vartype, 'region': region})
    return(df)


def combine_regions_germ_vars(regions={'autosomes': 2929051733, 'chr22': 51304566},
        germ_vars = {'snp': 4e6, 'indel': 0.8e6}):
    '''
    Scales the number of germline variants according to the length of genomic
    regions

    Parameters:
    regions: a dictionary of lengths in bases for various regions
    germ_vars: the approximate expected number of germline snps and indels

    Returns:
    a pandoc data frame of approx number of germline variants at all combinations of region and variant type
    '''
    gk = germ_vars.keys()
    rk = regions.keys()
    res = {v: [int(germ_vars[v] * regions[r] / regions['autosomes']) for r in rk] for v in gk}
    res = pd.DataFrame(res, index=rk)
    return(res)


def evalmodel2df_all(nvariants, germ_vars=combine_regions_germ_vars()):
    '''
    Calls evalmodel2df for all combination of regions, vartypes samples,
    models and total n.o. variants

    Returns:
    a pandoc data frame with the number of modeled variants at all those combinations and AAFs
    '''
    regions = ['chr22', 'autosomes']
    vartypes = ['snp', 'indel']
    samples = ['mix1', 'mix2', 'mix3']
    models = ['M1', 'M2', 'M3', 'M4']
    p_som2germ = [2, 3, 4]
    l = [evalmodel2df(nvariants, sample=s, vartype=v, region=r, model=m,
        Y=germ_vars.at[r, v] * 10 ** (- s2g), p_som2germ=s2g)
            for s in samples for v in vartypes for r in regions for m in
            models for s2g in p_som2germ]
    return(pd.concat(l))

def get_taejeongs_vaf_sample(sample='S316', scale2pct=True):
    '''
    A single sample version of get_taejeongs_vaf; see details therein
    '''
    csv = '/big/data/bsm/Bae-2018-science/aan8690_TableS1/' + sample + '.csv'
    fr_cx = pd.read_csv(csv)['FR-CX']
    def helper(y):
        if isinstance(y, str) and re.search('somatic', y):
            def splitter(ix):
                return(np.int64(str(y).split(sep=':')[ix]))
            alt = splitter(-1)
            ref = splitter(-2)
            total = alt + ref
            return(alt / total)
        else:
            return(np.nan)
    vaf = [helper(y) for y in fr_cx]
    vaf = [y for y in vaf if not np.isnan(y)]
    if scale2pct:
        vaf = [100 * y for y in vaf]
    df = pd.DataFrame({'VAF': vaf, 'sample': sample})
    df = df.astype({'sample': 'category'})
    return(df)


def get_taejeongs_vaf(samples=['S316', 'S320'], scale2pct=True):
    '''
    Import Taejeong's VAF for all variants in a list of samples from his Science article

    See
    https://science.sciencemag.org/highwire/filestream/703017/field_highwire_adjunct_files/1/aan8690_TableS1.xlsx
    If multiple samples are given, a "pooled" sample is created from them.

    Parameters:
    samples: a list of samples (note that S275 fails to import)
    scale2pct: True if VAF should be in percent

    Returns:
    a pandas DataFrame with a column of VAF values and a sample column
    '''
    l = [get_taejeongs_vaf_sample(s, scale2pct=scale2pct) for s in samples]
    res = pd.concat(l)
    if len(samples) > 1:
        pooled = pd.DataFrame({'VAF': res['VAF'], 'sample': 'pooled'})
        res = pd.concat([res, pooled])
        res = res.astype({'sample': 'category'})
    #res = sum(l, [])
    return(res)


def lambda_hat(vaf):
    '''
    MLE of the rate parameter lambda from a sample of VAF

    Parameters:
    vaf: a list of VAF values

    Returns:
    MLE of lambda
    '''
    lhat = len(vaf) / sum(vaf)
    return(lhat)


def vaf_distplot(vafdf=get_taejeongs_vaf(), fit=stats.expon):
    '''
    Plot the distribution of VAF values 

    If the input vafdf contains multiple samples then a multi-plot grid is
    created.

    Parameters:
    vafdf: the pandas DataFrame output of get_taejeongs_vaf
    fit: None or a scipy.stats distribution

    Returns: a seaborn FacetGrid object
    '''
    sns.set()
    sns.set_context('talk')
    g =  sns.FacetGrid(vafdf, row='sample', aspect=2.5, height=4)
    g.map(sns.distplot, 'VAF', hist=True, rug=True, kde=False,
            fit=fit)
    if fit is None:
        pass
    else:
        pass
    return(g)


def vaf_distplot1(vafdf):
    samples = vafdf['sample'].cat.categories
    nsamples = len(samples)
    fig, ax = plt.subplots(nsamples, 1)
    def plot_sample(sample):
        ix = samples.searchsorted(sample)
        n, bins, patches = ax[ix].hist(vafdf[vafdf['sample'] == sample]['VAF'])
        res = pd.DataFrame({'count': list(n) + [list(n)[-1]], 'bins': bins,
            'sample': sample})
        ax[ix].set_title(sample)
        return(res)
    histo = [plot_sample(s) for s in samples]
    histo = pd.concat(histo)
    g = sns.FacetGrid(histo, row='sample')
    g.map(plt.step, 'bins', 'count', where='post')
    return(histo)
