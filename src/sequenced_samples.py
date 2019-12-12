#! /usr/bin/env python3

import pandas as pd
import os
import os.path
import subprocess

def sequenced_samples_dissectionID(seqsamples='/projects/bsm/attila/results/2018-09-12-sequenced-individuals/sequenced-samples',
        dnalib='/projects/bsm/data/dnalib/BSM_Project_Chess.csv',
        outcsv='/projects/bsm/attila/results/2018-09-12-sequenced-individuals/sequenced-samples.csv',
        selcols=['Individual ID', 'Dx', 'Tissue', 'Institution Dissection ID']):
    '''
    Get Institution Dissection ID, Individual ID and other info for each
    sequenced sample and write it to a CSV

    Arguments
    seqsamples: the output of sequenced-samples.sh
    dnalib: input CSV: DNA libraries (Mehaa's table)
    outcsv: output CSV: sequenced samples
    selcols: selected columns from dnalib

    Value:
    a table of sequenced samples as a pandas DataFrame

    Details:
    The function uses as input file seqsamples, the ouptut of
    sequenced-samples.sh, which in turn depends on the output of
    sequenced-individuals.sh.  The sequenced_individuals function is a wrapper
    around these two BASH scripts.
    '''
    dnalib = pd.read_csv(dnalib)
    with open(seqsamples) as f:
        ss = f.read().splitlines()
    ss = [s.split('\t') for s in ss]
    ss = [pd.DataFrame({'indIDshort': s[0], 'tissue': s[2:-1]}) for s in ss]
    ss = pd.concat(ss)
    ss.index = range(len(ss.index))
    def helper(ix):
        tissue = ss.iat[ix, 1]
        t2abbrev = {'NeuN_pl': 'np', 'NeuN_mn': 'nn', 'muscle': 'mu'}
        abbrev = t2abbrev[tissue]
        bool0 = dnalib['Individual ID Short'] == ss.iat[ix, 0]
        bool1 = dnalib['Sample'] == abbrev
        val = dnalib.loc[bool0 & bool1, ]
        return(val)
    val = pd.concat([helper(i) for i in ss.index])
    val['Tissue'] = [{'np': 'NeuN_pl', 'nn': 'NeuN_mn', 'mu': 'muscle'}[y] for
            y in val['Sample']]
    val = val.loc[val['Library Replicate'] == 1, ]
    val.index = range(len(val.index))
    def get_fastq_names(indIDshort, tissue):
        fqnames = '/projects/bsm/alignments/' + indIDshort + os.sep + \
                indIDshort + '_' + tissue + '-fastq-names'
        return(fqnames)
    val['FASTQ Names'] = [get_fastq_names(x[0], x[1]) for x in zip(val['Individual ID Short'], val['Tissue'])]
    val['File Exists'] = [os.path.exists(x) for x in val['FASTQ Names']]
    selcols = selcols + ['FASTQ Names', 'File Exists']
    val = val[selcols]
    val.to_csv(outcsv, index=False)
    print('Results written to', outcsv)
    return(val)


def sequenced_individuals():
    args1 = ['/home/attila/projects/bsm/src/sequenced-individuals.sh']
    proc1 = subprocess.run(args1, capture_output=True)
    args2 = ['/home/attila/projects/bsm/src/sequenced-samples.sh']
    proc2 = subprocess.run(args2, capture_output=True)
    return(proc2)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outcsv', help='output CSV: sequenced samples',
            default='/projects/bsm/attila/results/2018-09-12-sequenced-individuals/sequenced-samples.csv')
    parser.add_argument('-d', '--dnalib', help='input CSV: DNA libraries (Mehaa\'s table)',
            default='/projects/bsm/data/dnalib/BSM_Project_Chess.csv')
    args = parser.parse_args()
    sequenced_individuals()
    sequenced_samples_dissectionID(outcsv=args.outcsv, dnalib=args.dnalib)
