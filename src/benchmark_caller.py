#! /usr/bin/env python3

import subprocess
import os
import os.path
import shutil

__maindir__ = '/projects/bsm/calls/mixing-experiment'
__BAMdir__ = '/projects/bsm/alignments/ceph-benchmark/'
__init_cfg__ = '/projects/bsm/calls/MSSM_179/NeuN_pl-muscle/JointSNVMix2/MSSM_179-NeuN_pl-muscle.cfg'
__REFSEQ__ = '/projects/shared/refgenome/GRCh37/dna/hs37d5.fa'


def call(caller='somaticSniper', case='Mix1', control='Mix2', nproc=1, overwrite=False):
    '''
    Wrapper around multiCaller that uses the mixing experiment data

    Parameters:
    caller: the name of the caller
    case: case sample e.g. Mix1
    control: control sample e.g. Mix2
    nproc: number of cores
    overwrite: whether to overwrite the caller's subdirectory

    Returns:
    the caller's subdirectory (string)
    '''
    sdir = case + '-' + control
    outdir = __maindir__ + os.path.sep + sdir
    if caller == 'strelka2Germline':
        t_opt = 'Germline'
    else:
        t_opt = 'Somatic'
    callerdir = outdir + os.path.sep + caller
    if os.path.exists(callerdir):
        if overwrite:
            shutil.rmtree(callerdir)
        else:
            print('The directory for the present caller/sample combination already exists.  Exiting...')
            return(None)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=False)
    case_replica = 'A'
    if case == control:
        control_replica = 'B'
    else:
        control_replica = 'A'
    case_sample = case + case_replica 
    control_sample = control + control_replica 
    caseBAM = __BAMdir__ + os.path.sep + case_sample + '.bam'
    controlBAM = __BAMdir__ + os.path.sep + control_sample + '.bam'
    args = ['multiCaller', '-p', nproc, '-r', __REFSEQ__, '-1', caseBAM, '-2', controlBAM,
            '-t', t_opt, '-i', __init_cfg__, '-a', case_sample, '-b', control_sample,
            '-o', outdir, caller]
    subprocess.run(args=args)
    return(sdir)

if __name__ == '__main__':
    import sys
    call(caller=sys.argv[1], case=sys.argv[2], control=sys.argv[3],
            nproc=sys.argv[4], overwrite=False)
