import pandas as pd
import shutil
import subprocess
import os
import joint_gt_ceph as jgc

def make_ts_aaf(mix='mix1', vartype='snp', region='chr22',
        tsdir='/home/attila/projects/bsm/results/2019-03-18-truth-sets/chr22/snp/truthset'):
        def helper(aaf):
            genotypes = gt_of_aaf[mix][aaf]
            invcfs = [indir + os.sep + g + '.vcf.gz' for g in genotypes]
            outvcf = mixdir + os.sep + str(aaf) + '.vcf.gz'
            unsorted_outvcf = mixdir + os.sep + str(aaf) + '-unsorted.vcf.gz'
            args0 = ['bcftools', 'concat', '-o', unsorted_outvcf, '-Oz'] + invcfs
            args1 = ['bcftools', 'sort', '-o', outvcf, '-Oz', unsorted_outvcf]
            args2 = ['bcftools', 'index', '-t', outvcf]
            if not os.path.isfile(outvcf):
                subprocess.run(args0)
                subprocess.run(args1)
                os.remove(unsorted_outvcf)
                subprocess.run(args2)
            # count records
            args3 = ['bcftools', 'view', '-H', outvcf]
            args4 = ['wc', '-l']
            proc1 = subprocess.Popen(args3, shell=False, stdout=subprocess.PIPE)
            proc2 = subprocess.Popen(args4, shell=False, stdout=subprocess.PIPE, stdin=proc1.stdout)
            proc1.stdout.close()
            nrec = proc2.communicate()[0]
            nrec = int(nrec) # turn bytesliteral (e.g. b'5226\n') to integer
            return(nrec)
        gt_of_aaf = jgc.get_gt_of_aaf()
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
        aafs = gt_of_aaf[mix].keys()
        nrec = [helper(f) for f in aafs]
        res = pd.DataFrame({'region': region, 'vartype': vartype, 'sample': mix, 'aaf': list(aafs),
                'nvariants': nrec})
        res = res.astype({'region': 'category', 'vartype': 'category', 'sample': 'category'})
        return(res)
