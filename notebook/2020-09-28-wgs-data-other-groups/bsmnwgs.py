import synapseclient
import pandas as pd
import os.path

syn = synapseclient.login()
wdir = '/home/attila/projects/bsm/results/2020-09-28-wgs-data-other-groups/'

def foo():
    return('Hello World!')

def get_manifest(synid):
    e = syn.get(synid, downloadLocation=wdir + os.path.sep + 'downloads', ifcollision='overwrite.local')
    manifest = pd.read_csv(e.path, skiprows=1)
    #manifest['project'] = project
    return(manifest)

def get_manifests(synids, project, maniftype='genomics_subject02'):
    l = [get_manifest(s) for s in synids]
    manifest = pd.concat(l, axis=0)
    manifest['project'] = project
    fpath = wdir + os.path.sep + project + '-' + maniftype + '.csv'
    manifest.to_csv(fpath)
    return(manifest)
