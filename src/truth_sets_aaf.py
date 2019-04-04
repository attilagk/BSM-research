import shutil
import subprocess
import os
import joint_gt_ceph as jgc

def helper(aaf, mix='mix1',
        tsdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets/chr22/snp/truthset',
        gt_of_aaf=jgc.get_gt_of_aaf()):
    genotypes = gt_of_aaf[mix][aaf]
    indir = tsdir + os.sep + 'genotypes'
    aafdir = tsdir + os.sep + 'aaf'
    if not os.path.isdir(aafdir):
        os.mkdir(aafdir)
    unfdir = aafdir + os.sep + 'unfiltered'
    if not os.path.isdir(unfdir):
        os.mkdir(unfdir)
    mixdir = unfdir + os.sep + mix
    if not os.path.isdir(mixdir):
        os.mkdir(mixdir)
    invcfs = [indir + os.sep + g + '.vcf.gz' for g in genotypes]
    outvcf = mixdir + os.sep + str(aaf) + '.vcf.gz'
    unsorted_outvcf = mixdir + os.sep + str(aaf) + '-unsorted.vcf.gz'
    args0 = ['bcftools', 'concat', '-o', unsorted_outvcf, '-Oz'] + invcfs
    args1 = ['bcftools', 'sort', '-o', outvcf, '-Oz', unsorted_outvcf]
    args2 = ['bcftools', 'index', '-t', outvcf]
    subprocess.run(args0)
    subprocess.run(args1)
    subprocess.run(args2)
    args = (args0, args1, args2)
    os.remove(unsorted_outvcf)
    return(args)
