#! /usr/bin/env python3

import synapseclient
import pandas as pd
import os
import sys

#syn = synapseclient.login()

def get_manifest(synapse_id, syn, skiprows=1, download_dir="/tmp/"):
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
    entity = syn.get(synapse_id, downloadLocation = download_dir, ifcollision="overwrite.local")
    df = pd.read_csv(download_dir + os.sep + entity.properties.name, skiprows=skiprows)
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
    Adds the first header row of the manifest, which is something like this:
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

def extract_subject(template, subject):
    df = template.loc[template["src_subject_id"] == subject, :]
    return(df)

def make_manifests(subject, syn, target_dir="."):
    '''
    Makes all manifest files for a given CMC subject

    Parameters
    subject: a CMC subject_id with or without the "CMC_" prefix
    syn: a synapse object returned by synapseclient.login()
    target_dir: directory for the manifest files created
    '''
    def btb_or_gsubj(template_syn_id):
        manif_temp, manif_syn = get_manifest(template_syn_id, syn, download_dir=download_dir)
        manif = extract_subject(manif_temp, cmc_subject)
        temp_p = download_dir + os.sep + manif_syn.properties.name # template path
        targ_p = target_dir + os.sep + cmc_subject + "-" + os.path.basename(temp_p)
        write_manifest(manif, temp_p, targ_p)
        return(manif)

    subject = subject.replace("CMC_", "") # ensure that subject lacks CMC_ prefix
    cmc_subject = "CMC_" + subject # add CMC_ prefix
    download_dir = target_dir
    return((btb_or_gsubj("syn12154562"), btb_or_gsubj("syn12128754")))

def make_g_sample(gsam_temp, gsubj):
    # obtain shared columns
    is_shared = [y in gsubj.columns for y in gsam_temp.columns]
    shared = gsam_temp.loc[:, is_shared].columns
    # creating genomics sample
    gsam = gsam_temp.reindex(index=list(range(gsubj.shape[0])))
    for col in shared:
        gsam.at[gsam.index[0], col] = gsubj.at[gsubj.index[0], col]
    return(gsam)
