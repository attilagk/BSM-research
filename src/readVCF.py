import pandas as pd
import subprocess
import io
import re

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

def readVCF(vcfpath, annotlist=read_annotlist()):
    formatstr = make_formatstr(annotlist)
    colnames = [y.replace('INFO/', '') for y in annotlist]
    cmd = ['bcftools', 'query', '-f', formatstr, vcfpath]
    p =  subprocess.run(cmd, capture_output=True)
    df = pd.read_csv(io.BytesIO(p.stdout), sep='\t', names=colnames)
    return(df)
