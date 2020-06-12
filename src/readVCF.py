import pandas as pd
import numpy as np
import subprocess
import io
import re
import synapseclient

def read_annotlist(annotpath='/home/attila/projects/bsm/tables/VCF-HC.annotations', withFORMAT=False):
    with open(annotpath) as f:
        l = f.readlines()
    l = [y.replace('\n', '') for y in l] # remove newline characters
    if not withFORMAT:
        l = [y for y in l if not re.match('^FORMAT', y)]
    return(l)

def make_formatstr(annotlist):
    formatstr = '%' + '\t%'.join(annotlist) + '\n'
    return(formatstr)

def sample_fromVCF(vcfpath):
    cmd = ['bcftools', 'view', '-h', vcfpath]
    p =  subprocess.run(cmd, capture_output=True)
    colnames = p.stdout.splitlines()[-1]
    sample = re.sub('^.*\t([^\t]+)$', '\\1', colnames.decode('utf-8'))
    return(sample)

def convert_sample(bsm_sample):
    s = re.sub('^((MSSM|PITT)_([0-9]+))_(.*)$', '\\1:\\4', bsm_sample)
    l = s.split(':')
    indiv_id = 'CMC_' + l[0]
    tissue = l[1]
    return((indiv_id, tissue))

def readVCF(vcfpath, annotlist=read_annotlist()):
    formatstr = make_formatstr(annotlist)
    colnames = [y.replace('INFO/', '') for y in annotlist]
    cmd = ['bcftools', 'query', '-f', formatstr, vcfpath]
    p =  subprocess.run(cmd, capture_output=True)
    df = pd.read_csv(io.BytesIO(p.stdout), sep='\t', names=colnames)
    return(df)

def add_clinical(vcfpath, clinical):
    vcf = readVCF(vcfpath)
    sample = sample_fromVCF(vcfpath)
    indiv_id = convert_sample(sample)[0]
    df = clinical.loc[[indiv_id]]
    df = df.iloc[np.zeros(vcf.shape[0]), :]
    df.index = vcf.index
    outvcf = pd.concat([vcf, df], axis=1)
    return(outvcf)

def readVCFs(vcflistpath='/big/results/bsm/calls/filtered-vcfs.tsv'):
    # CMC_Human_clinical_metadata.csv
    syn = synapseclient.login()
    wdir='/tmp'
    cmc_clinical_syn = syn.get('syn2279441', downloadLocation=wdir, ifcollision='overwrite.local')
    cmc_clinical = pd.read_csv(cmc_clinical_syn.path, index_col='Individual ID')
    # TODO
    add_clinical('/big/results/bsm/2020-06-10-chromatin-state/MSSM_109_brain.ploidy_50.filtered-annot.vcf',
            clinical)
    return(cmc_clinical)
