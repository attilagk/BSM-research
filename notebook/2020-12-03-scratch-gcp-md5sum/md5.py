import synapseclient
import synapseutils
import subprocess
import pandas as pd
import numpy as np

syn = synapseclient.login()
converter = '/home/attila/projects/bsm/notebook/2020-12-03-scratch-gcp-md5sum/base128base64.sh'
md5fetcher = '/home/attila/projects/bsm/notebook/2020-12-03-scratch-gcp-md5sum/md5gcp.sh'

def syn_show(synID, syn=syn):
    e = syn.get(synID, downloadFile=False)
    return(e)

def md5_from_entity(e, base64=False):
    md5 = e._file_handle['contentMd5']
    md5 = recode_md5(md5) if base64 else md5
    return(md5)

def recode_md5(md5):
    p = subprocess.run([converter, md5], capture_output=True)
    b = p.stdout # bytes object
    val = b.decode('utf-8').strip()
    return(val)

def gcp_md5_from_entity(e):
    fname = e.name
    p = subprocess.run([md5fetcher, fname], capture_output=True)
    return(p)

def check1file(synID, syn=syn):
    e = syn_show(synID, syn)
    p = gcp_md5_from_entity(e)
    md5dest = np.nan if p.returncode else p.stdout.decode('utf-8').strip()
    md5source = md5_from_entity(e, base64=True)
    match = md5source == md5dest
    d = {'match': match, 'source md5 (base64)': md5source, 'dest md5 (base64)': md5dest, 'synapse ID': synID}
    df = pd.DataFrame(d, index=[e.name])
    return(df)

def check_all_files(synfolderID='syn20735395', syn=syn):
    '''
    Checks all files in source folder Synapse scratch space by comparing md5sum to
    destination GCP store

    Parameters
    synfolderID: source folder
    syn: a synapse object

    Value: a dataframe showing match or mismatch for each file in source

    Details
    Dest md5 is NaN if the file is missing from dest but present in source.
    '''
    l = list(synapseutils.walk(syn, synfolderID))[0][2]
    synIDs = np.array(l)[:, 1]
    dflist = [check1file(y, syn) for y in synIDs]
    df = pd.concat(dflist, axis=0)
    return(df)
