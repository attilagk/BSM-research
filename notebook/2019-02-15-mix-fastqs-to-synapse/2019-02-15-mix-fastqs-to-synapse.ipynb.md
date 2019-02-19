Sequencing data from the DNA mixes of CEPH/Utah grandparents are uploaded to the BSMN [Scratch Space](https://www.synapse.org/#!Synapse:syn10964481).

This notebook must be run on server Ada since that's where the FASTQ files are stored.


```python
import synapseclient
import re
syn = synapseclient.login()
```

    Welcome, Attila Gulyás-Kovács!
    


The function below does the following things:

1. Get the pathnames of all FASTQs for a given sample (Mix1 or Mix3); this is based on `*-fastq-names`
1. Given those pathnames upload FASTQs to a Synapse folder


```python
def store_sample(sample, parent, syn=syn):
    '''
    Store FASTQ files for sample in Synapse folder parent.
    '''
    
    def local_files2list(bamdir="/projects/bsm/alignments/ceph-benchmark/"):
        '''
        Create a list of pathnames from the -fastq-names file determined by bamdir and sample.
        '''
        file_list = sample + "-fastq-names"
        with open(bamdir + file_list) as f:
            return([re.sub('\n', '', y) for y in f.readlines()])
        
    def store_files(local_files, folder):
        '''
        Store local files in pathname list local_files in folder in Synapse.
        '''
        stored_files = [ syn.store(synapseclient.File(y, parent=folder)) for y in local_files ]
        return(stored_files)
    
    data_folder = synapseclient.Folder(sample, parent=parent)
    data_folder = syn.store(data_folder)
    return(store_files(local_files2list(), folder = data_folder))
```

Now create and store the target Synapse folder...


```python
main_folder = "syn18233572" # bsmn-pipeline-test/genome_mapping folder on scratch space
parent = syn.store(synapseclient.Folder("ceph_utah_mixes", parent=main_folder))
```

...and finally store files in [the Synapse folder](https://www.synapse.org/#!Synapse:syn18345708)


```python
stored_files = {s: store_sample(sample=s, parent=parent) for s in ("Mix1A", "Mix3A")}
```
