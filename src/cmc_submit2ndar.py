#! /usr/bin/env python3

import synapseclient
import pandas as pd
import os
import sys

syn = synapseclient.login()

def get_manifest(synapse_id, skiprows=1, download_dir="/tmp/", syn=syn):
    '''
    Download manifest from Synapse and read it into a pandas data frame

    Parameters
    synapse_id: synapse ID of the manifest (string)
    skiprows: number of rows to skip when reading into a data frame
    download_dir: local directory to download manifest file into
    syn: a synapse object returned by synapseclient.login()

    Value
    A pandas data frame and synapse File entity (in a tuple)
    '''
    entity = syn.get(synapse_id, downloadLocation = download_dir)
    df = pd.read_csv(download_dir + entity.properties.name, skiprows=skiprows)
    print("Manifest file path: " + download_dir + entity.properties.name)
    return((df, entity))


def write_manifest(df, template_path, target_path):
    '''
    Creates a manifest file from a data frame adding header from a template 
    Parameters
    df: the manifest data in a pandas data frame with
    column labels
    template_path: path of the template manifest file
    for the first header row
    target_path: path of the manifest file written
    from df with the first header row

    Details
    The first header row is something like this:
    "nichd_btb","02"
    '''
    body_path = target_path + "~"
    df.to_csv(body_path, index=False)
    # get first header row from template manifest
    with open(template_path, "r") as f:
        btb_head = f.readline()
    # get manifest body as string
    with open(body_path, "r") as f:
        btb_body = f.read()
    # write first header row and manifest body into target csv
    with open(target_path, "w") as f:
        f.write(btb_head)
        f.write(btb_body)
    # clean up
    os.remove(body_path)


if False:
    syn = synapseclient.login()
