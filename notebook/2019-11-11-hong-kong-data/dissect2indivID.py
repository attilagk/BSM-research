#/usr/bin/env python3

import glob
import os
import os.path

readsdir = '/projects/bsm/reads/'
HKdir = readsdir + '2016-10-11-HK/'
dissfq_test = HKdir + 'S_1178_NeuN/S_1178_NeuN_DO16090243-D709_HW3C2CCXX_L1_1.fq.gz'
diss2ind = {'S_1178': 'MSSM295',
        'S_1163': 'MSSM304',
        'S_700': 'PITT101',
        'S_1307': 'PITT036'}

def rename1fq(dissfq=dissfq_test, dissID='S_1178'):
    indID = diss2ind[dissID]
    dissprefix = dissID + '_NeuN'
    indprefix = indID + '_NeuN_pl'
    indbn = os.path.basename(dissfq).replace(dissprefix, indprefix)
    inddir = readsdir + indprefix + os.path.sep
    indfq = inddir + indbn
    if os.path.exists(inddir) and os.path.isdir(inddir):
        pass
    else:
        os.mkdir(inddir)
    if not os.path.exists(indfq):
        os.symlink(src=dissfq, dst=indfq)
    return(indfq)

def rename1sample(dissID='S_1178'):
    dissprefix = dissID + '_NeuN'
    dissdir = HKdir + dissprefix + os.path.sep
    l = glob.glob(dissdir + '*.fq.gz')
    return(l)
