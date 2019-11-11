#!/usr/bin/env python3

import glob
import os
import os.path

readsdir = '/projects/bsm/reads/'
HKdir = readsdir + '2016-10-11-HK/'
DVdir = readsdir + '2016-08-15-DV-X10/'
dissfq_test = HKdir + 'S_1178_NeuN/S_1178_NeuN_DO16090243-D709_HW3C2CCXX_L1_1.fq.gz'
diss2ind = {\
        'S_1178': {'id': 'MSSM295', 'dir': HKdir},
        'S_1163': {'id': 'MSSM304', 'dir': HKdir},
        'S_700': {'id': 'PITT101', 'dir': HKdir},
        'S_1307': {'id': 'PITT036', 'dir': HKdir},
        'S_1164': {'id': 'MSSM193', 'dir': DVdir},
        'S_1181': {'id': 'MSSM056', 'dir': DVdir},
        'S_1148': {'id': 'MSSM331', 'dir': DVdir},
        'S_1542': {'id': 'PITT060', 'dir': DVdir}}

def rename1sample(dissID='S_1178'):
    '''
    Make (renamed) symlinks to FASTQs for dissID

    Arguments
    dissID: dissection ID

    Value: list of pathnames to symlinks named by individual ID
    '''
    dissprefix = dissID + '_NeuN'
    if diss2ind[dissID]['dir'] == HKdir:
        dissprefix1 = dissprefix
    elif diss2ind[dissID]['dir'] == DVdir:
        dissprefix1 = dissID + 'NeuN'
    def rename1fq(dissfq):
        '''
        Make (renamed) symlink to FASTQ

        Arguments
        dissfq: path to FASTQ named by dissection ID

        Value: pathname to symlink named by individual ID
        '''
        indID = diss2ind[dissID]['id']
        indprefix = indID + '_NeuN_pl'
        indbn = os.path.basename(dissfq).replace(dissprefix1, indprefix)
        inddir = readsdir + indprefix + os.path.sep
        indfq = inddir + indbn
        if os.path.exists(inddir) and os.path.isdir(inddir):
            pass
        else:
            os.mkdir(inddir)
        if not os.path.exists(indfq):
            os.symlink(src=dissfq, dst=indfq)
        return(indfq)
    dissdir = diss2ind[dissID]['dir'] + dissprefix + os.path.sep
    it = glob.iglob(dissdir + '*.fq.gz')
    res = [rename1fq(dissfq=x) for x in it]
    return(res)

def rename_all_samples():
    l = [rename1sample(y) for y in diss2ind.keys()]
    return(l)

if __name__ == "__main__":
    rename_all_samples()
