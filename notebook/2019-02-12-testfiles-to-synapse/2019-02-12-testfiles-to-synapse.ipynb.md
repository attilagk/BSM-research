
Small FASTQ files for testing the BSMN mapping workflow are uploaded to various locations in Synapse.


```python
import synapseclient
import glob
syn = synapseclient.login()
```

    Welcome, Attila Gulyás-Kovács!
    


## BSMN Chess Lab

Select all (four) small fastq files sampled from large fastqs for the common sample and store them in the `bsmn-pipeline-test` folder in the BSMN Chess Lab project.

Clean up folder and upload files again.


```python
#folderId = "syn17931318"
for e in syn.getChildren("syn17931318"):
    syn.delete(syn.get(e, downloadFile = False))
    
local_files = glob.glob("/big/data/bsm/MS_02May2016_CommonSample/*.fastq.gz")
def my_upload(local_files = local_files, folderId = "syn17931318"):
    for f in local_files:
        syn.store(synapseclient.File(f, parent = folderId))
        
my_upload()
```

    
    ##################################################
     Uploading file to Synapse storage 
    ##################################################
    
    
    ##################################################
     Uploading file to Synapse storage 
    ##################################################
    
    
    ##################################################
     Uploading file to Synapse storage 
    ##################################################
    
    
    ##################################################
     Uploading file to Synapse storage 
    ##################################################
    


## Scratch Space

The same files are uploaded but the parent folder now is the `small_test_data` on the BSMN Scratch Space.


```python
my_upload(folderId = "syn18233615")
```
