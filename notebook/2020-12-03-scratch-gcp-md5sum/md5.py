import synapseclient
import subprocess

syn = synapseclient.login()
converter = '/home/attila/projects/bsm/notebook/2020-12-03-scratch-gcp-md5sum/base128base64.sh'

def syn_show(synID, syn=syn):
    e = syn.get(synID, downloadFile=False)
    return(e)

def md5_from_entity(e, base64=False):
    md5 = e._file_handle['contentMd5']
    return(md5)

def recode_md5(md5):
    p = subprocess.run([converter, md5], capture_output=True)
    b = p.stdout # bytes object
    val = b.decode('utf-8').strip()
    return(val)
