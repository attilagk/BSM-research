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
    btb = btb_or_gsubj("syn12154562")
    gsubj = btb_or_gsubj("syn12128754")
    # genomics samples
    gsam_temp, gsam_syn = get_manifest("syn8464096", syn, download_dir=download_dir)
    gsam = make_g_sample(gsam_temp, gsubj)
    temp_p = download_dir + os.sep + gsam_syn.properties.name
    targ_p = target_dir + os.sep + cmc_subject + "-" + os.path.basename(temp_p)
    write_manifest(gsam, temp_p, targ_p)
    return((btb, gsubj, gsam))

def make_g_sample(gsam_temp, gsubj):
    # obtain shared columns
    is_shared = [y in gsubj.columns for y in gsam_temp.columns]
    shared = gsam_temp.loc[:, is_shared].columns
    # creating genomics sample with missing values
    gsam = gsam_temp.reindex(index=list(range(len(gsubj))))
    # copying values from genomics subject
    for col in shared:
        gsam.at[gsam.index[0], col] = gsubj.at[gsubj.index[0], col]
    # remaining columns
    experiment_id = 1111
    organism = 'human'
    sample_amount = 1 # made up
    sample_unit = 'NA' # made up
    data_file1_type = 'FASTQ'
    data_file1 = '2016-12-15-DV-X10/MSSM106_muscle/MSSM106_muscle_USPD16080279-D702_H7YNMALXX_L6_1.fq.g'
    data_file2_type = 'FASTQ'
    data_file2 = '2016-12-15-DV-X10/MSSM106_muscle/MSSM106_muscle_USPD16080279-D702_H7YNMALXX_L6_2.fq.g'
    storage_protocol = 'NA' # made up
    data_file_location = 'NDAR'
    patient_id_biorepository = gsam.at[gsam.index[0], 'src_subject_id']
    sample_id_biorepository = 'MSSM_DNA_TMPR_69087' # from syn17021773 CMC_Human_WGS_metadata_working.csv
    gsam['experiment_id'] = experiment_id
    gsam['organism'] = organism
    gsam['sample_amount'] = sample_amount
    gsam['sample_unit'] = sample_unit
    gsam['data_file1_type'] = data_file1_type
    gsam['data_file1'] = data_file1
    gsam['data_file2_type'] = data_file2_type
    gsam['data_file2'] = data_file2
    gsam['storage_protocol'] = storage_protocol
    gsam['data_file_location'] = data_file_location
    gsam['patient_id_biorepository'] = patient_id_biorepository
    gsam['sample_id_biorepository'] = sample_id_biorepository
    # path argument for vtcmd -l option
    l = ['/projects/bsm/reads/', '/projects/bsm/alignments/']
    return(gsam)
