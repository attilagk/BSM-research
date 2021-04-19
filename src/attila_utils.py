import os
import os.path
import subprocess

def savefig(fig, basename, dirpath='named-figure/', formats=['png', 'pdf']):
    '''
    Save figure in multiple formats given a base name

    Parameters:
    fig: a figure object whose (base) class is matplotlib.figure.Figure
    basename: basename of the output file
    dirpath: the path to the output directory, typically named-figure/
    formats: list of output formats

    Returns: the list of output filepaths
    '''
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
    bpath = dirpath + basename
    def save_one(fm):
        filepath = bpath + '.' + fm
        fig.savefig(fname=filepath, format=fm, bbox_inches="tight")
        return(filepath)
    res = [save_one(fm) for fm in formats]
    return(res)


def sentieon_licsrvr():
    args0 = ['pgrep', '-f', 'licsrvr']
    proc0 = subprocess.run(args0)
    if proc0.returncode == 0:
        # license server is already running
        return(proc0)
    sentieon_lic = os.environ['SENTIEON_LICENSE']
    sentieon_dir = os.environ['SENTIEON_INSTALL_DIR']
    sentieon_exec = sentieon_dir + os.path.sep + 'bin/sentieon'
    sentieon_log = '/home/attila/tools/sentieon.log'
    args1 = [sentieon_exec, 'licsrvr', '--start', '--log', sentieon_log, sentieon_lic]
    proc1 = subprocess.run(args1)
    return(proc1)
