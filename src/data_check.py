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

def idxstats_contig_plot(idxstats):
    '''
    Plots % of unmapped reads for each contig

    Arguments
    idxstats: the value of idxstats, a pandas DataFrame

    Value
    a matplotlib figure
    '''
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[8, 12], sharey=True)
    sns.stripplot(x='% Unmapped', y='Name', data=idxstats, hue='Sample', ax=ax[0])
    sns.stripplot(x='% Unmapped', y='Name', data=idxstats, hue='Sample', ax=ax[1])
    ax[0].grid('x')
    ax[1].grid('x')
    ax[1].set_xlim([-0.5, 5])
    ax[1].set_ylabel('')
    ax[1].legend().set_visible(False)
    ax[0].set_title('full range')
    ax[1].set_title('zoom')
    return(fig)


def selfSMplot(selfSM):
    '''
    Plot the FREEMIX and CHIPMIX values from VerifyBamID's .selfSM output files

    Arguments
    selfSM: a pandas DataFrame with pd.read_csv('some.selfSM, sep='\t')

    Value
    a matplotlib figure
    '''
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
