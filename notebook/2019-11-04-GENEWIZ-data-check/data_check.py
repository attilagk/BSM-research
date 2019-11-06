import io
import os.path
import subprocess
import pandas as pd

def idxstats(indiv='MSSM_033', celltype='NeuN_pl'):
    BAMpath = '/projects/bsm/alignments/' + indiv + os.path.sep + indiv + '_' + celltype + '.bam'
    proc = subprocess.run(['samtools', 'idxstats', BAMpath], capture_output=True)
    string = io.StringIO(proc.stdout.decode('ASCII'))
    df = pd.read_csv(string, sep='\t', header=None, names=['Name', 'Length',
            'Mapped', 'Unmapped'])
    df['Unmapped frac'] = df['Unmapped'] / (df['Mapped'] + df['Unmapped'])
    return(df)
