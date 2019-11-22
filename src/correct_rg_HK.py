#! /usr/bin/env python3

import os
import os.path
import subprocess
import glob
import re

diss2subj = {'S_1178': 'MSSM295', 'S_1163': 'MSSM304', 'S_700': 'PITT101', 'S_1307': 'PITT036'}

def split_bam(bam='/projects/bsm/alignments/PITT_101/PITT_101_NeuN_pl-22-1Mb.bam'):
    maindir = os.path.dirname(bam)
    bname = os.path.basename(bam).replace('.bam', '')
    subdir = maindir + os.path.sep + bname
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    mydir = os.getcwd()
    os.chdir(subdir)
    args = ['samtools', 'split', bam]
    proc = subprocess.run(args)
    splitbams = glob.glob(subdir + os.path.sep + '*.bam')
    os.chdir(mydir)
    return(splitbams)


def correct_rg_1bam(bam, oldID='S_700', newID=None):
    '''
    Correct read group for a single BAM by replacing the old sample ID with a new one

    Arguments
    bam: path to the input BAM with incorrect read group
    oldID: the old sample ID
    newID: the new sample ID; if None then get it from diss2subj

    Value
    the path to the corrected BAM
    '''
    if newID is None:
        newID = diss2subj[oldID]
    args = ['samtools', 'view', '-H', bam]
    proc = subprocess.run(args, capture_output=True)
    m = re.search(b'@RG.*\\n', proc.stdout)
    oldrg = m.group(0).decode('utf-8')
    s = re.sub(oldID, newID, oldrg)
    s = re.sub('(@RG)?\\t', ' -r ', s)
    s = re.sub('^ *', '', s)
    s = re.sub('\\n', '', s)
    bambn = os.path.basename(bam)
    bamdir = os.path.dirname(bam)
    corrbam = bamdir + os.path.sep + 'corrected-' + bambn
    args = ['samtools', 'addreplacerg'] + s.split(' ') + ['-o', corrbam, bam]
    if subprocess.run(args, capture_output=True):
        return(corrbam)
    else:
        return(None)

#echo -r "ID:$ID" -r "LB:$LB" -r "SM:$SM" -r "PU:$PU" -r "CN:$CN" -r "PL:$PL"
