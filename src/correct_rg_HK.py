#! /usr/bin/env python3

import os
import os.path
import shutil
import subprocess
import glob
import re

subj2diss = {'MSSM_295': 'S_1178', 'MSSM_304': 'S_1163', 'PITT_101': 'S_700', 'PITT_036': 'S_1307'}

def split_bam(bam='/projects/bsm/alignments/PITT_101/PITT_101_NeuN_pl-22-1Mb.bam'):
    '''
    Splits the BAM with incorrect read group into multiple BAMs in a subdir named after BAM

    Argument
    bam: path to bam

    Value
    the list of pathnames to the split BAMs
    '''
    maindir = os.path.dirname(bam)
    bname = os.path.basename(bam).replace('.bam', '')
    subdir = maindir + os.path.sep + bname
    if os.path.exists(subdir):
        shutil.rmtree(subdir)
    os.mkdir(subdir)
    mydir = os.getcwd()
    os.chdir(subdir)
    args = ['samtools', 'split', bam]
    proc = subprocess.run(args)
    splitbams = glob.glob(subdir + os.path.sep + '*.bam')
    os.chdir(mydir)
    return(splitbams)


def correct_rg_splitbam(bam):
    '''
    Correct read group for a single BAM by replacing the old sample ID with a new one

    Arguments
    bam: path to the input BAM with incorrect read group
    oldID: the old sample ID
    newID: the new sample ID; if None then get it from diss2subj

    Value
    the path to the corrected BAM
    '''
    bambn = os.path.basename(bam)
    bamdir = os.path.dirname(bam)
    indivID = re.search('(MSSM|PITT)_[0-9]+', bambn).group(0)
    newID = indivID.replace('_', '')
    oldID = subj2diss[indivID]
    args = ['samtools', 'view', '-H', bam]
    proc = subprocess.run(args, capture_output=True)
    m = re.search(b'@RG.*\\n', proc.stdout)
    oldrg = m.group(0).decode('utf-8')
    s = re.sub(oldID, newID, oldrg)
    s = re.sub('(@RG)?\\t', ' -r ', s)
    s = re.sub('^ *', '', s)
    s = re.sub('\\n', '', s)
    corrbam = bamdir + os.path.sep + 'corrected-' + bambn
    args = ['samtools', 'addreplacerg'] + s.split(' ') + ['-o', corrbam, bam]
    if subprocess.run(args, capture_output=True):
        return(corrbam)
    else:
        return(None)


def merge_correct_bams(bams, nthreads=4, maxMem='2G'):
    '''
    Merge and sort the corrected BAMs

    Arguments
    bams: a list of paths to the corrected BAMs
    nthreads: the number of total threads to run samtools sort
    maxMem: the max memory per thread for samtools sort

    Value
    path to the merged, sorted BAM

    Details
    An earlier version of this function failed on big BAMs with messages like this:
    > samtools sort: fail to open "PITT_101_NeuN_pl.bam.tmp.1020.bam": Too many open files
    The solution was to increase the number of threads and the maximum memory
    per thread from the default values. See this GitHub issue:
    https://github.com/samtools/samtools/issues/603
    '''
    addthreads = str(nthreads - 1)
    bamdir = os.path.dirname(bams[0])
    mergedbam = bamdir + os.path.sep + 'merged.bam'
    sortedbam = bamdir + os.path.sep + 'sorted.bam'
    args1 = ['samtools', 'merge'] + [mergedbam] + bams
    proc1 = subprocess.run(args1, capture_output=True)
    args2 = ['samtools', 'sort', '-@', addthreads, '-m', maxMem, '-o', sortedbam, mergedbam]
    proc2 = subprocess.run(args2, capture_output=True)
    return(sortedbam)


def correct_rg(bam):
    '''
    Main function that corrects read groups in a BAM

    Argument
    bam: path to the BAM to be corrected

    Value
    a procedure object representing the indexing of the corrected BAM
    '''
    bam = os.path.realpath(bam)
    splitbams = split_bam(bam)
    corrbams = [correct_rg_splitbam(b) for b in splitbams]
    corrbam = merge_correct_bams(corrbams)
    oldbam =  bam + '~'
    os.replace(bam, oldbam)
    os.replace(corrbam, bam)
    args = ['samtools', 'index', '-b', bam]
    proc = subprocess.run(args, capture_output=True)
    bamdir = os.path.dirname(corrbam)
    shutil.rmtree(bamdir) and os.remove(oldbam)
    return(proc)


if __name__ == '__main__':
        import sys
        correct_rg(sys.argv[1])
