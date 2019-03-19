import subprocess

def count_records(filter_expr='', only_biallelic=True):
    '''
    Count records filtered by type

    The Python pipe idiom was lifted from here:
    https://security.openstack.org/guidelines/dg_avoid-shell-true.html#correct
    '''
    biall_filter = []
    if only_biallelic:
        biall_filter = ['-m', '2', '-M', '2']
    #vcf_path = '/big/data/platinum-genomes/ceph-utah-vars/illumina-calls/S1/NA12889_S1.vcf.gz'
    vcf_path = '/home/attila/projects/bsm/results/calls/ceph-utah/illumina/NA12889.vcf.gz'
    args1 = ['bcftools', 'view', '-H', '-i'] + [filter_expr] + biall_filter + [vcf_path]
    args2 = ['wc', '-l']
    proc1 = subprocess.Popen(args1, shell=False, stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(args2, shell=False, stdout=subprocess.PIPE, stdin=proc1.stdout)
    proc1.stdout.close()
    res = proc2.communicate()[0]
    res = int(res) # turn bytesliteral (e.g. b'5226\n') to integer
    return(res)
