import io
import os.path
import subprocess
import pandas as pd

def idxstats(indiv='MSSM_033', celltype='NeuN_pl'):
    '''
    Run samtools idxstats on a BAM and import results

    Arguments
    indiv: the individual ID
    celltype: NeuN_pl | NeuN_mn | muscle

    Value
    a pandas Data Frame

    Details
    The following two columns are added to the 4-column output of idxstats:
    1) % Unmapped (reads)
    2) Sample (= indiv + '_' + celltype)
    '''
    sample = indiv + '_' + celltype
    BAMpath = '/projects/bsm/alignments/' + indiv + os.path.sep + sample + '.bam'
    proc = subprocess.run(['samtools', 'idxstats', BAMpath], capture_output=True)
    string = io.StringIO(proc.stdout.decode('ASCII'))
    df = pd.read_csv(string, sep='\t', header=None, names=['Name', 'Length',
            'Mapped', 'Unmapped'])
    df['% Unmapped'] = df['Unmapped'] / (df['Mapped'] + df['Unmapped']) * 100
    df['Sample'] = sample
    return(df)
