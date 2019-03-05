#! /usr/bin/env python3

import synapseclient
import numpy as np
import pandas as pd
import os
import sys
import glob

#syn = synapseclient.login()
experiment_id = 1223

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
    df.to_csv(body_path, index=False, date_format='%m/%d/%Y')
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
        manif_temp, manif_syn = get_manifest(template_syn_id, syn, download_dir=target_dir)
        manif = extract_subject(manif_temp, cmc_subject)
        manif = correct_manifest(manif)
        temp_p = target_dir + os.sep + manif_syn.properties.name # template path
        targ_p = target_dir + os.sep + cmc_subject + "-" + os.path.basename(temp_p)
        write_manifest(manif, temp_p, targ_p)
        return(manif)

    def g_sample(gsubj):
        gsam_temp, gsam_syn = get_manifest("syn8464096", syn, download_dir=target_dir)
        gsam = make_g_sample(gsam_temp, btb, gsubj, syn)
        gsam = correct_manifest(gsam)
        temp_p = target_dir + os.sep + gsam_syn.properties.name
        targ_p = target_dir + os.sep + cmc_subject + "-" + 'genomics_sample03_U01MH106891_Chess.csv'
        write_manifest(gsam, temp_p, targ_p)
        return(gsam)
    subject = subject.replace("CMC_", "") # ensure that subject lacks CMC_ prefix
    cmc_subject = "CMC_" + subject # add CMC_ prefix
    btb = btb_or_gsubj("syn12154562")
    gsubj = btb_or_gsubj("syn12128754")
    gsam = g_sample(gsubj)
    return((btb, gsubj, gsam))

def manifest_type(df):
    '''
    Returns the type of manifest.

    Parameter
    df: the manifest (a pandoc data frame)

    Values
    btb: brain and tissue bank
    gsubj: genomics subjects
    gsam: genomics samples
    None: undetermined
    '''
    if df.columns[1] == 'experiment_id':
        return('gsam') # genomics samples
    else:
        if 'disorder' in df.columns:
            return('btb') # brain and tissue bank
        elif 'phenotype' in df.columns:
            return('gsubj') # genomics subjects

def sample_specifics(sample_id_original):
    '''
    Sample specifics (celltype and br_reg) based on sample_id_original

    Parameter(s)
    sample_id_original: a string like MSSM_106.DLPFC_1399.np1

    Value: a dictionary with keys celltype and br_reg (brain region)
    '''
    sid = sample_id_original
    if '.np1' in sid or '.np2' in sid:
        celltype = 'NeuN+'
        br_reg = 'prefrontal cortex'
    if '.nn1' in sid or '.nn2' in sid:
        celltype = 'NeuN-'
        br_reg = 'prefrontal cortex'
    if '.mu1' in sid or '.mu2' in sid:
        celltype = 'muscle'
        br_reg = np.nan
    val = {'celltype': celltype, 'br_reg': br_reg}
    return(val)

def correct_manifest(df):
    '''
    Corrects manifest of any type

    Corrections are based on the validation results file below
    /projects/bsm/attila/results/2019-02-19-upload-to-ndar/validation_results_20190226T163227.csv
    '''
    def btb_sample_specs(element):
        '''Return a list of sample specifics based on the sample_id_original
        column of df'''
        val = [sample_specifics(y)[element] for y in df['sample_id_original']]
        return(val)

    res = df.copy()
    res['interview_date'] = pd.to_datetime(df['interview_date'])
    if manifest_type(res) in ('btb', 'gsubj'):
        res.loc[res['race'] == 'African American', 'race'] = 'Black or African American'
    if manifest_type(res) == 'btb':
        res['celltype'] = btb_sample_specs('celltype')
        res['br_reg'] = btb_sample_specs('br_reg')
    if manifest_type(res) == 'gsubj':
        res['family_study'] = 'No'
        res['sample_description'] = 'brain'
        res['patient_id_biorepository'] = res['src_subject_id']
        res['sample_id_biorepository'] = res['src_subject_id']
    if manifest_type(res) == 'gsam':
        if any([pd.isna(y) for y in res['sample_amount']]):
            res['sample_amount'] = 'NaN'
    return(res)

def get_sample_id_original(tissue, btb):
    '''
    Returns sample_id_original for a "tissue" and a brain and tissue bank
    manifest

    Parameters
    tissue: one of NeuN_pl, NeuN_mn or muscle
    btb: a pandas data frame, a brain and tissue bank manifest

    Value
    sample_id_original, a string
    '''
    tissue_id = {'NeuN_pl': 'np1', 'NeuN_mn': 'nn1', 'muscle': 'mu1'}
    suffix = tissue_id[tissue]
    ids = list(btb['sample_id_original'])
    res = [s for s in ids if suffix in s]
    return(res[0])

def extract_cmc_wgs(btb, syn):
    '''
    Extract rows from CMC_Human_WGS_metadata_working.csv based on btb

    syn17021773 CMC_Human_WGS_metadata_working.csv 
    '''
    wgs, wgs_syn = get_manifest('syn17021773', syn, skiprows=0)
    ids = list(btb['sample_id_original'])
    wgs = wgs[wgs['Library ID'].isin(ids)]
    return(wgs)

def make_g_sample(gsam_temp, btb, gsubj, syn):
    '''
    Creates a genomics sample manifest based on a genomics sample template and
    two other manifests

    Parameters
    gsam_temp: the genomics sample template, a pandas data frame
    btb: brain and tissue bank manifest, a pandas data frame
    gsubj: genomics subject manifest, a pandas data frame
    syn: a synapse object returned by synapseclient.login()

    Value: genomics sample manifest, a pandas data frame

    Details: The script deduces the sample types ("tissues") based on the BAMs
    found
    '''
    def do_tissue(tissue):
        '''Creates a tissue-specific genomics sample'''
        def do_file(fl):
            '''Creates a file-specific genomics sample'''
            df = gsam.copy(deep=True) # deep copy
            # copying values from genomics subject
            for col in shared:
                df.at[df.index[0], col] = gsubj.at[gsubj.index[0], col]
            df['sample_description'] = sample_description
            df['sample_id_original'] = lib_id
            wgs_lib = wgs[wgs['Library ID'] == lib_id]
            df['sample_id_biorepository'] = wgs_lib.at[wgs_lib.index[0], 'Sample DNA ID']
            df['sample_amount'] = wgs_lib.at[wgs_lib.index[0], 'DNA Amount(ng)']
            df['sample_unit'] = 'ng'
            df['data_file1'] = fl[0]
            if '.fq.gz' in fl[0]:
                df['data_file1_type'] = 'FASTQ'
                df['data_file2_type'] = 'FASTQ'
                df['data_file2'] = fl[1]
            elif '.bam' in fl[0]:
                df['data_file1_type'] = 'BAM'
            return(df)

        fq_names = fastq_names[tissue]
        with open(fq_names) as fqs:
            fastqs = fqs.readlines()
        fastqs_1 = sorted([s.replace('\n', '') for s in fastqs if '_1.fq.gz' in s])
        fastqs_2 = sorted([s.replace('\n', '') for s in fastqs if '_2.fq.gz' in s])
        data_files = list(zip(fastqs_1, fastqs_2)) + [(bams[tissue], )]
        sample_description = 'frontal cortex'
        if tissue == 'muscle':
            sample_description = 'muscle'
        lib_id = get_sample_id_original(tissue, btb)
        val = pd.concat([do_file(f) for f in data_files])
        return(val)

    # from syn17021773 CMC_Human_WGS_metadata_working.csv 
    wgs = extract_cmc_wgs(btb, syn)
    # ensure that gsubj has one and only one row
    if len(gsubj) != 1:
        raise Exception('genomics subjects manifest must have one and only one row')
    # get paths of BAMs and of the fastq-names files
    src_subject_id = gsubj.at[gsubj.index[0], 'src_subject_id']
    simple_id = src_subject_id.replace('CMC_', '')
    aln_p = '/projects/bsm/alignments/' + simple_id + os.sep
    bams = glob.glob(aln_p + simple_id + '*.bam')
    fastq_names = [s.replace('.bam', '-fastq-names') for s in bams]
    tissues = [s.replace(aln_p + simple_id + '_', '').replace('.bam', '') for s in bams]
    # these variables are referred to in inside functions
    bams = dict(zip(tissues, bams))
    fastq_names = dict(zip(tissues, fastq_names))
    # obtain shared columns
    is_shared = [y in gsubj.columns for y in gsam_temp.columns]
    shared = gsam_temp.loc[:, is_shared].columns
    # creating genomics sample with missing values
    gsam = gsam_temp.reindex(index=list(range(len(gsubj))))
    # filling out values shared for all tissues
    gsam['experiment_id'] = experiment_id
    gsam['organism'] = 'human'
    gsam['data_file_location'] = 'NDAR'
    gsam['storage_protocol'] = 'NA' # made up
    gsam['patient_id_biorepository'] = src_subject_id
    val = pd.concat([do_tissue(t) for t in tissues])
    return(val)
