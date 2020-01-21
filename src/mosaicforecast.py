#! /usr/bin/env python3

import subprocess
import os
import os.path
import shutil
import pandas as pd
import re

SegDup_and_clustered_bed = '/home/attila/projects/MosaicForecast/resources/SegDup_and_clustered.bed'
mt2_pon_filter_script = '/home/attila/projects/MosaicForecast/MuTect2-PoN_filter.py'

def mt_pon_filter(invcf):
    if re.match('.*\.vcf.gz$', invcf):
        vcf = invcf.replace('.gz', '')
        args = ['bcftools', 'view', '-Ov', '-o', vcf, invcf]
        proc = subprocess.run(args, capture_output=True)
        pass
    elif re.match('.*\.vcf$', invcf):
        vcf = invcf
    else:
        raise Exception('Error: ' + invcf + ' is not a .vcf or .vcf.gz file')
    sample = 'test1'
    tmpbed = sample + '.bed'
    outbed = os.path.dirname(invcf) + os.sep + sample + '.bed'
    args = ['python3', mt2_pon_filter_script, sample, vcf, SegDup_and_clustered_bed]
    proc = subprocess.run(args, capture_output=True)
    shutil.move(tmpbed, outbed)
    #if proc.returncode == 0:
    return(proc)


def bed2regions_file(bedfile='/projects/bsm/attila/results/2020-01-21-mosaicforecast-dev/test.bed'):
    # this depends on the MuTect2-PoN_filter.py script
    colnames = ['chr', 'pos0', 'pos1', 'ref', 'alt', 'sample', 'depth', 'AF']
    bed = pd.read_csv('/projects/bsm/attila/results/2020-01-21-mosaicforecast-dev/test.bed',
            sep='\t', header=None, names=colnames)
    # positions in bed files are 0 based so add 1 to get 1 based positions
    bed['pos'] = bed['pos0'] + 1
    regfname = os.path.dirname(bedfile) + os.sep + os.path.basename(bedfile).replace('.bed', '.regions')
    reg = bed.loc[:, ['chr', 'pos']]
    reg.to_csv(regfname, sep='\t', header=False, index=False)
    return(reg)
