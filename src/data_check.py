import io
import os.path
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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

def selfSMplot(selfSM):
    plt.style.use('seaborn-notebook')
    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True)
    a = sns.stripplot(x='FREEMIX', y='sample', data=selfSM, jitter=False,
            color='blue', ax=ax[0])
    b = sns.stripplot(x='CHIPMIX', y='sample', data=selfSM, jitter=False,
            color='blue', ax=ax[1])
    a.axes.set_xlim([-0.05, 1.05])
    a.axes.set_title('Sample Impurity Evidence')
    a.grid('both')
    b.grid('both')
    b.axes.set_xlim([-0.05, 1.05])
    b.axes.set_title('Sample Swap Evidence')
    b.axes.set_ylabel('')
    fig.suptitle('verifyBamID Results')
    return(fig)
