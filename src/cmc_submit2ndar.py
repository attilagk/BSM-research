#! /usr/bin/env python3

import synapseclient
import numpy as np
import pandas as pd
import os
import os.path
import sys
import glob
import subprocess
import io
import re
import datetime

cmc_metadata = {'CMC_Human_clinical_metadata': 'syn2279441',
        'CMC_Human_brainRegion_metadata': 'syn21446693',
        'CMC_Human_isolation_metadata_DNA': 'syn2279444'}

experiment_id = 1223
chess_grant =  'U01MH106891'
manifest_template_synids = {'nichd_btb02': "syn12154562", 'genomics_subject02': "syn12128754", 'genomics_sample03': "syn8464096"}
genewiz_serialn_synid = 'syn21982509' # a.k.a samples-from-Chaggai.csv
chess_s3_bucket = 's3://chesslab-bsmn' 

def read_dfiles(syn, fpath='/home/attila/projects/bsm/results/2021-02-02-submit-to-nda/files_not_uploaded_yet.MSSM.txt'):
    '''
    See 2021-02-02-submit-to-nda: Data files to submit
    '''
    dfiles = pd.read_csv(fpath, sep='\t', names=['synapseID', 'filename'])
    s = dfiles['filename']
    s = s.str.split('_NeuN_pl.').apply(lambda x: ['CMC_' + x[0], x[1]])
    index = pd.MultiIndex.from_tuples(s, names=['indivID', 'filetype'])
    dfiles = pd.DataFrame(dfiles.to_numpy(), index=index, columns=dfiles.columns)
    dfiles = dfiles.sort_index(axis=0, level='filetype')
    def get_dfile_path(synID):
        key = syn.get(entity=synID, downloadFile=False)._file_handle['key']
        l = key.split('/')
        prefix = l.pop(0)
        dfpath = '/'.join(l)
        return(dfpath)
    dfiles['data_file1'] = dfiles['synapseID'].apply(get_dfile_path)
    return(dfiles)


def edit_gsam(gsam, dfiles, gender, dftype='cram'):
    '''
    Edit gsam by replacing data_file1_type with dftyp and data_file1 with synapse S3 path

    Parameters
    gsam: template; the same samples but FASTQ file type
    dfiles: list of data files see read_dfiles and 2021-02-02-submit-to-nda
    gender: a series of gender values (F or M)
    dftype: one of 10 data file types; cram, cram.crai,...

    Value: the edited gsam (copy)
    '''
    columns = gsam.columns
    gsam = gsam.groupby('src_subject_id', as_index=False).first()
    gsam = gsam.set_index('src_subject_id', drop=False)
    gsam = gsam.reindex(columns=columns)
    dfiles = dfiles.xs(dftype, axis=0, level='filetype')
    gsam['data_file1_type'] = dftype
    gsam['data_file1'] = dfiles['data_file1']
    # remove samples with missing file; see CMC_PITT_117
    gsam = gsam.dropna(axis=0, subset=['data_file1'])
    gsam['data_file2_type'] = np.nan
    gsam['data_file2'] = np.nan
    gsam['data_file3_type'] = np.nan
    gsam['data_file3'] = np.nan
    gsam['data_file4_type'] = np.nan
    gsam['data_file4'] = np.nan
    gsam['sex'] = gender
    gsam['sample_unit'] = gsam['sample_unit'].fillna('NA')
    gsam['storage_protocol'] = gsam['storage_protocol'].fillna('NA')
    return(gsam)


def empty_manifest_row(manifest):
    '''
    Returns the NaN-filled last row of a manifest, a pandas DataFrame
    '''
    lastrow = manifest.iloc[-1:, :]
    a = np.array(lastrow)
    a.fill(np.nan)
    new = pd.DataFrame(a, columns=lastrow.columns)
    return(new)

def fillin_gsub_or_btb_row(indiv_id, manif, cmc_clinical, cmc_brainreg, genewiz_serialn):
    '''
    Given a CMC indiv ID fill in a one-rowed gsub or btb manifest

    Parameters
    indiv_id: a CMC individual ID (string)
    manif: the manifest template (pandas DataFrame)
    cmc_clinical: pandas DataFrame holding CMC clinical info
    cmc_brainreg: pandas DataFrame holding CMC brainreg info
    genewiz_serialn: pandas DataFrame with a map GENEWIZ serialnumbers to CMC indiv IDs
    
    Value: a one-rowed gsub or btb manifest
    '''
    manifr = empty_manifest_row(manif)
    cmc = cmc_clinical.loc[indiv_id]
    manifr['src_subject_id'] = indiv_id
    manifr['interview_date'] = datetime.date.today().strftime('%m/%d/%y')
    manifr['interview_age'] = int(cmc['ageOfDeath'] * 12)
    manifr['gender'] = cmc['Reported Gender'][0]
    if cmc['Race'] is np.nan:
        manifr['race'] = 'Unknown or not reported'
    else:
        manifr['race'] = cmc['Race']
    manifr['ethnic_group'] = cmc['Ethnicity']
    instdissectionID = get_instdissectionID(indiv_id, cmc_brainreg, genewiz_serialn)
    manifr['sample_id_original'] = instdissectionID + '.np1'
    return(manifr)

def fillin_gsub_row(indiv_id, gsub, cmc_clinical, cmc_brainreg, genewiz_serialn, syn=None, tissue='NeuN_pl'):
    '''
    Given a CMC indiv ID fill in a one-rowed gsub

    Parameters
    indiv_id: a CMC individual ID (string)
    manif: the manifest template (pandas DataFrame)
    cmc_clinical: pandas DataFrame holding CMC clinical info
    cmc_brainreg: pandas DataFrame holding CMC brainreg info
    genewiz_serialn: pandas DataFrame with a map GENEWIZ serialnumbers to CMC indiv IDs
    syn: a synapse object
    tissue: for the S3 samples it's invariably NeuN_pl
    
    Value: a one-rowed gsubmanifest
    '''
    simple_id = indiv_id.replace('CMC_', '')
    indiv_id = 'CMC_' + simple_id
    gsubr = fillin_gsub_or_btb_row(indiv_id, gsub, cmc_clinical, cmc_brainreg, genewiz_serialn)
    cmc = cmc_clinical.loc[indiv_id]
    gsubr['phenotype'] = cmc['Dx']
    gsubr['phenotype_description'] = 'No'
    gsubr['twins_study'] = 'No'
    gsubr['sibling_study'] = 'No'
    gsubr['sample_taken'] = 'Yes'
    gsubr['sample_description'] = 'PFC'
    if re.match('.*MSSM.*', indiv_id):
        gsubr['biorepository'] = 'MSBB'
    elif re.match('.*PITT.*', indiv_id):
        gsubr['biorepository'] = 'UPittNBB'
    return(gsubr)


def fillin_btb_row(indiv_id, btb, cmc_clinical, cmc_brainreg, genewiz_serialn, syn=None, tissue='NeuN_pl'):
    '''
    Given a CMC indiv ID fill in a one-rowed btb manifest.  For details see
    the doc string of fillin_gsub_row.
    '''
    simple_id = indiv_id.replace('CMC_', '')
    indiv_id = 'CMC_' + simple_id
    btbr = fillin_gsub_or_btb_row(indiv_id, btb, cmc_clinical, cmc_brainreg, genewiz_serialn)
    return(btbr)


def get_instdissectionID(indiv_id, cmc_brainreg, genewiz_serialn):
    '''
    Deduce CMC Institution Dissection ID for a CMC Individual ID

    Details: Because for one individual there can be (and typically are)
    multiple dissections additional info is used here to select a single
    Dissection ID.  These include genewiz_serialn, cmc_brainreg and some
    arbitrary rules.

    Parameters: see fillin_gsub_row docstring for details

    Value: Institution Dissection ID
    '''
    brainr = cmc_brainreg.loc[cmc_brainreg['Individual ID'] == indiv_id, :]
    simple_id = indiv_id.replace('CMC_', '')
    PFCn = genewiz_serialn.loc[simple_id, 'PFC #']
    return(PFCn)
    if re.match('^CMC_.*', PFCn):
        return(PFCn)
    matches = [y for y in brainr['Institution Dissection ID'] if
            re.match('^.*(DRPC|PFC).*' + PFCn + '.*$', y) is not None]
    psych = [y for y in matches if re.match('^.*PsychENCODE.*$', y) is not
            None]
    if len(matches) == 0:
        print('Institution Dissection ID for', indiv_id, 'not found')
        return('TODO')
    if len(psych) >= 1:
        return(psych[0])
    else:
        return(matches[0])

def add_subj_key(manif, pGUIDpath='/home/attila/projects/bsm/results/2020-04-22-upload-to-ndar-from-s3/s3-pseudo-guids'):
    '''
    Add subject keys (pseudo GUIDs) to manifest.

    Details: pseudo GUIDs were created by Andy and sent over to me via email on 2020-04-27
    '''
    wdir = '~/projects/bsm/results/2020-04-22-upload-to-ndar-from-s3/'
    pGUIDs = list(pd.read_csv(wdir + 's3-pseudo-guids', header=None)[0])
    manif['subjectkey'] = pGUIDs[:manif.shape[0]]
    return(manif)

def resources_for_make_manif_s3(wdir):
    '''
    Download and return a bunch of resources to make manifests for data in our S3 bucket
    '''
    syn = synapseclient.login()
    gsubtempl, gsub_syn = get_manifest(manifest_template_synids['genomics_subject02'], syn, download_dir=wdir)
    btbtempl, btb_syn = get_manifest(manifest_template_synids['nichd_btb02'], syn, download_dir=wdir)
    gsamtempl, gsam_syn = get_manifest(manifest_template_synids['genomics_sample03'], syn, download_dir=wdir)
    # CMC_Human_clinical_metadata.csv
    cmc_clinical_syn = syn.get('syn2279441', downloadLocation=wdir, ifcollision='overwrite.local')
    cmc_clinical = pd.read_csv(cmc_clinical_syn.path, index_col='Individual ID')
    # CMC_Human_brainRegion_metadata.csv
    cmc_brainreg_syn = syn.get('syn21446693', downloadLocation=wdir, ifcollision='overwrite.local')
    cmc_brainreg = pd.read_csv(cmc_brainreg_syn.path)
    # originally created by Chaggai but manually edited by Attila with Institution Dissection IDs
    genewiz_serialn_syn = syn.get(genewiz_serialn_synid, downloadLocation=wdir, ifcollision='overwrite.local')
    genewiz_serialn = pd.read_csv(genewiz_serialn_syn.path, index_col='CMC_simple_id')
    return((syn, gsubtempl, gsub_syn, btbtempl, btb_syn, gsamtempl, gsam_syn, cmc_clinical, cmc_brainreg, genewiz_serialn))


def make_manif_s3(wdir = '/home/attila/projects/bsm/results/2020-04-22-upload-to-ndar-from-s3/'):
    '''
    Make all three manifests for the data at s3://chesslab-bsmn/GENEWIZ

    Parameters
    wdir: directory path to download resources and create manifests in

    Value: the three manifests

    Details: The prefix of the manifest files is in the format of YYYY-MM-DD-
    reflecting the date of invocation.
    '''
    syn, gsubtempl, gsub_syn, btbtempl, btb_syn, gsamtempl, gsam_syn, cmc_clinical, cmc_brainreg, genewiz_serialn = resources_for_make_manif_s3(wdir)
    tissue = 'NeuN_pl'
    def helper(fun=fillin_gsub_row, maniftempl=gsubtempl):
        l = [fun(s, maniftempl, cmc_clinical, cmc_brainreg, genewiz_serialn, syn, tissue) for s in genewiz_serialn.index]
        df = pd.concat(l, axis=0)
        df = add_subj_key(df)
        df = correct_manifest(df)
        return(df)
    # make genomics subjects
    gsub = helper(fillin_gsub_row, gsubtempl)
    btb = helper(fillin_btb_row, btbtempl)
    #gsam = helper(fillin_gsam_rows_scratch_space, gsamtempl)
    gsam = helper(fillin_gsam_rows_chess_s3, gsamtempl)
    # retain those subjects/samples only that are present in gsam
    gsub.index = gsub['src_subject_id']
    gsub = gsub.loc[gsam['src_subject_id'], :]
    btb.index = btb['src_subject_id']
    btb = btb.loc[gsam['src_subject_id'], :]
    # write data frames to CSVs
    manifs = {'nichd_btb02': (btb, btb_syn), 'genomics_subject02': (gsub, gsub_syn), 'genomics_sample03': (gsam, gsam_syn)}
    #def writer(manif=gsub, manif_syn=gsub_syn, manifname='gsub'):
    def writer(manifname='nichd_btb02'):
        manif = manifs[manifname][0]
        manif_syn = manifs[manifname][1]
        datestr = datetime.date.today().strftime('%Y-%m-%d')
        csvpath = wdir + os.sep + datestr + '-' + manifname + '.csv'
        csvpath = os.path.normpath(csvpath)
        write_manifest(manif, manif_syn.path, csvpath)
        print(manifname, 'written to', csvpath)
        return(None)
    [writer(k) for k in manifs.keys()]
    return((gsub, btb, gsam))


def fillin_gsam_rows_chess_s3(indiv_id, gsam_temp, cmc_clinical, cmc_brainreg, genewiz_serialn, syn, tissue='NeuN_pl'):
    simple_id = indiv_id.replace('CMC_', '')
    indiv_id = 'CMC_' + simple_id
    fastq_names = get_fastq_names_s3(simple_id, genewiz_serialn)
    gsam = gsam_temp.copy()
    gsam['data_file1'] = [fastq_names[0]]
    gsam['data_file2'] = [fastq_names[1]]
    gsam['data_file1_type'] = ['FASTQ']
    gsam['data_file2_type'] = ['FASTQ']
    gsam = fillin_gsam_rows(indiv_id, gsam, cmc_clinical, cmc_brainreg, genewiz_serialn)
    return(gsam)


def fillin_gsam_rows_scratch_space(indiv_id, gsam_temp, cmc_clinical, cmc_brainreg, genewiz_serialn, syn, tissue='NeuN_pl'):
    '''
    Fill in multiple rows of gsam corresponding to a single CMC individual

    Parameters: see fillin_gsub_row docstring for details

    Value: gsam manifest for the individual
    '''
    simple_id = indiv_id.replace('CMC_', '')
    indiv_id = 'CMC_' + simple_id
    def helper(ftype1='CRAM'):
        if ftype1 == 'CRAM':
            ftype2 = 'CRAM index'
            el1 = get_scratch_space_data_s3(simple_id, syn, ext='.cram', tissue=tissue)
            el2 = get_scratch_space_data_s3(simple_id, syn, ext='.cram.crai', tissue=tissue)
        elif ftype1 == 'VCF':
            ftype2 = 'VCF index'
            el1 = get_scratch_space_data_s3(simple_id, syn, ext='.vcf.gz', tissue=tissue)
            el2 = get_scratch_space_data_s3(simple_id, syn, ext='.vcf.gz.tbi', tissue=tissue)
        gsaml = list()
        for e1, e2 in zip(el1, el2):
            gsam = gsam_temp.copy()
            gsam['data_file1'] = [e1._file_handle['key']]
            gsam['data_file2'] = [e1._file_handle['key']]
            if gsam.shape[0] != 0:
                gsaml.append(gsam)
        if len(gsaml) == 0:
            return(None)
        gsam = pd.concat(gsaml, axis=0)
        gsam['data_file1_type'] = ftype1
        gsam['data_file2_type'] = ftype2
        return(gsam)
    gsaml = [helper(ft) for ft in ['CRAM', 'VCF']]
    if all([y is None for y in gsaml]):
        return(None)
    gsam = pd.concat(gsaml, axis=0)
    gsam = fillin_gsam_rows(indiv_id, gsam, cmc_clinical, cmc_brainreg, genewiz_serialn)
    return(gsam)


def fillin_gsam_rows(indiv_id, gsam, cmc_clinical, cmc_brainreg, genewiz_serialn):
    '''
    Lower level function analogous to fillin_gsub_row et al; called by fillin_gsam_rows_scratch_space
    '''
    cmc = cmc_clinical.loc[indiv_id]
    gsam['experiment_id'] = experiment_id
    gsam['src_subject_id'] = indiv_id
    gsam['interview_age'] = int(cmc['ageOfDeath'] * 12)
    gsam['interview_date'] = datetime.date.today().strftime('%m/%d/%y')
    gsam['sample_description'] = 'frontal cortex'
    instdissectionID = get_instdissectionID(indiv_id, cmc_brainreg, genewiz_serialn)
    gsam['sample_id_original'] = instdissectionID + '.np1'
    gsam['organism'] = 'human'
    gsam['sample_amount'] = 1
    gsam['sample_unit'] = 'NA'
    gsam['storage_protocol'] = 'NA' # made up
    gsam['data_file_location'] = 'NDAR'
    if re.match('.*MSSM.*', indiv_id):
        gsam['biorepository'] = 'MSBB'
    elif re.match('.*PITT.*', indiv_id):
        gsam['biorepository'] = 'UPittNBB'
    gsam['patient_id_biorepository'] = indiv_id
    gsam['sample_id_biorepository'] = gsam['sample_id_original']
    gsam['site'] = chess_grant
    return(gsam)



def get_fastq_names_s3_helper(simple_id, genewiz_serialn, s3prefix='GENEWIZ/30-317737003/'):
    '''
    Get fastq names in S3 that match a CMC subject simple_id
    '''
    prefix = genewiz_serialn.loc[simple_id, 'GENEWIZ_serialn']
    l = ['aws', 's3', 'ls', chess_s3_bucket + '/' + s3prefix]
    p = subprocess.run(l, capture_output=True)
    cnames = ['date', 'time', 'size', 'filename']
    s3ls = pd.read_csv(io.BytesIO(p.stdout), sep='\s+', names=cnames)
    filenames = s3ls['filename']
    def ismatch(fn):
        m = re.match('^' + prefix + '_[RS].*$', fn)
        return(m is not None)
    matches = [fn for fn in filenames if ismatch(fn)]
    # Adding here s3 prefix is not needed because the following syntax takes
    # care of it:
    # $ vtcmd -pre GENEWIZ/30-317737003
    #matches = [s3prefix + m for m in matches]
    return(matches)

def get_fastq_names_s3(simple_id, genewiz_serialn):
    #s3prefixes = ['GENEWIZ/30-317737003/', 'GENEWIZ/30-317737003_12_lanes/']
    s3prefixes = ['GENEWIZ/30-317737003/']
    fastq_names = [get_fastq_names_s3_helper(simple_id, genewiz_serialn, y)
            for y in s3prefixes]
    fastq_names = list(np.array(fastq_names).ravel())
    return(fastq_names)


def get_scratch_space_data_s3(simple_id, syn, ext='.cram', tissue='NeuN_pl'):
    if ext == '.cram' or ext == '.cram.crai' or ext == '.crai':
        syn_folder_id='syn20735395'
    elif ext == '.vcf.gz' or ext == '.vcf.gz.tbi':
        syn_folder_id='syn20812330'
    entity_list = list()
    for item in syn.getChildren(syn_folder_id):
        if re.match('^.*' + simple_id + '_' + tissue + '.*' + ext + '$', item['name']):
            e = syn.get(item['id'], downloadFile=False)
            entity_list.append(e)
    return(entity_list)


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

def get_genewiz_serialn(syn):
    e = syn.get(genewiz_serialn_synid)
    df = pd.read_csv(e['path'], index_col='CMC_simple_id')
    return(df)


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

def make_manifests(subject, syn, target_dir=".", matching_sample_ids=True, tissue=None):
    '''
    Makes all manifest files for a given CMC subject

    Parameters
    subject: a CMC subject_id with or without the "CMC_" prefix
    syn: a synapse object returned by synapseclient.login()
    target_dir: directory for the manifest files created
    matching_sample_ids: whether sample_id_biorepository should match sample_id_original

    Value
    a tuple of the three manifests, each in a pandas DataFrame
    '''
    def btb_or_gsubj(template_syn_id):
        manif_temp, manif_syn = get_manifest(template_syn_id, syn, download_dir=target_dir)
        manif = extract_subject(manif_temp, cmc_subject)
        manif = correct_manifest(manif)
        temp_p = target_dir + os.sep + manif_syn.properties.name # template path
        targ_p = target_dir + os.sep + cmc_subject + "-" + os.path.basename(temp_p)
        if tissue is None:
            write_manifest(manif, temp_p, targ_p)
        return(manif)

    def g_sample(gsubj):
        gsam_temp, gsam_syn = get_manifest(manifest_template_synids['genomics_sample03'], syn, download_dir=target_dir)
        return((gsam_temp, gsam_syn))
        gsam = make_g_sample(gsam_temp, btb, gsubj, syn,
                matching_sample_ids=matching_sample_ids, tissue=tissue)
        gsam = correct_manifest(gsam)
        temp_p = target_dir + os.sep + gsam_syn.properties.name
        targ_p = target_dir + os.sep + cmc_subject + "-" + 'genomics_sample03_U01MH106891_Chess.csv'
        gsam['site'] = chess_grant
        if tissue is None:
            write_manifest(gsam, temp_p, targ_p)
        return(gsam)
    subject = subject.replace("CMC_", "") # ensure that subject lacks CMC_ prefix
    print('processing', subject, tissue)
    cmc_subject = "CMC_" + subject # add CMC_ prefix
    btb = btb_or_gsubj(manifest_template_synids['nichd_btb02'])
    gsubj = btb_or_gsubj(manifest_template_synids['genomics_subject02'])
    # TODO: something is not right below
    gsam_temp, gsam_syn = g_sample(gsubj)
    return((btb, gsubj, gsam_temp, gsam_syn))
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

    def safely_remove_prefix(s, prefix='/projects/bsm/'):
        if pd.isna(s):
            return(s)
        else:
            return(s.replace(prefix, ''))

    res = df.copy()
    res['interview_date'] = pd.to_datetime(df['interview_date'])
    if manifest_type(res) in ('btb', 'gsubj'):
        #res.loc[res['race'] == 'African American', 'race'] = 'Black or African American'
        res.loc[res['ethnic_group'] == 'African-American', 'race'] = 'Black or African American'
        res.loc[res['ethnic_group'] == 'Caucasian', 'race'] = 'White'
        res.loc[res['ethnic_group'] == 'Hispanic', 'race'] = 'White'
        #res.loc[[r is np.nan for r in res['race']], 'race'] = 'Unknown or not reported'
    if manifest_type(res) == 'btb':
        res['celltype'] = btb_sample_specs('celltype')
        res['br_reg'] = btb_sample_specs('br_reg')
    if manifest_type(res) == 'gsubj':
        res['family_study'] = 'No'
        res['sample_description'] = 'brain'
        res['patient_id_biorepository'] = res['src_subject_id']
        res['sample_id_biorepository'] = res['src_subject_id']
    if manifest_type(res) == 'gsam':
        res['data_file1'] = [safely_remove_prefix(s) for s in res['data_file1']]
        res['data_file2'] = [safely_remove_prefix(s) for s in res['data_file2']]
        if any([pd.isna(y) for y in res['sample_amount']]):
            #res['sample_amount'] = 'NaN'
            res['sample_amount'] = 1
            res['sample_unit'] = 'NA'
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

def make_g_sample(gsam_temp, btb, gsubj, syn, matching_sample_ids=True,
        tissue=None, s3prefix=None, genewiz_serialn=None):
    '''
    Creates a genomics sample manifest based on a genomics sample template and
    two other manifests

    Parameters
    gsam_temp: the genomics sample template, a pandas data frame
    btb: brain and tissue bank manifest, a pandas data frame
    gsubj: genomics subject manifest, a pandas data frame
    syn: a synapse object returned by synapseclient.login()
    matching_sample_ids: whether sample_id_biorepository should match sample_id_original
    tissue: a single tissue in {NeuN_pl, NeuN_mn, muscle} can be given

    Value: genomics sample manifest, a pandas data frame

    Details: The script deduces the sample types ("tissues") based on the BAMs
    found when tissue=None or based on tissue when it's one of NeuN_pl,
    NeuN_mn, muscle.
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
            if matching_sample_ids:
                df['sample_id_biorepository'] = lib_id
            else:
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

        # TODO
        # new feature needed: get fastq names from somewhere else than path to
        # BAM files
        if s3prefix is None:
            fq_names = fastq_names[tissue]
            with open(fq_names) as fqs:
                fastqs = fqs.readlines()
            fastqs_1 = sorted([s.replace('\n', '') for s in fastqs if '_1.fq.gz' in s])
            fastqs_2 = sorted([s.replace('\n', '') for s in fastqs if '_2.fq.gz' in s])
            data_files = list(zip(fastqs_1, fastqs_2)) + [(bams[tissue], )]
        else:
            # get fastq names from the S3 bucket using genewiz_serialn
            fq_names = get_fastq_names_s3(simple_id, genewiz_serialn, s3prefix)
            data_files = fq_names
            #fq_names = dict(zip(tissue, fastq_names))
            return(data_files)
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
    src_subject_id = gsubj.at[gsubj.index[0], 'src_subject_id']
    simple_id = src_subject_id.replace('CMC_', '')
    if s3prefix is None:
        # get fastq names from the BAMs and the fastq-names files
        aln_p = '/projects/bsm/alignments/' + simple_id + os.sep
        bams = glob.glob(aln_p + simple_id + '*.bam')
        fastq_names = [s.replace('.bam', '-fastq-names') for s in bams]
        tissues = [s.replace(aln_p + simple_id + '_', '').replace('.bam', '') for s in bams]
        # these variables are referred to in inside functions
        bams = dict(zip(tissues, bams))
        fastq_names = dict(zip(tissues, fastq_names))
    #return(fastq_names)
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
    if tissue is None:
        val = pd.concat([do_tissue(t) for t in tissues])
    else:
        val = do_tissue(tissue)
    return(val)


def make_manifests_main(slistcsv, target_dir=".", prefix='chess-'):
    '''
    Make manifests given a list of samples in a CSV file

    Arguments:
    slistcsv: path to the input CSV file
    target_dir: where the manifest files will be created
    prefix: common prefix to all three manifests

    Value:
    a tuple of the three manifests, each one a pandas DataFrame

    Details:
    The side effect of the function is to write the manifests in the target
    directory with the given prefix.
    '''
    syn = synapseclient.login()
    def do_one_sample(subject, tissue):
        m = make_manifests(subject=subject, tissue=tissue, syn=syn,
                target_dir=target_dir, matching_sample_ids=True)
        return(m)
    samples = pd.read_csv(slistcsv)
    lomanifests = [do_one_sample(subject=samples.iloc[i][0], tissue=samples.iloc[i][1]) for i in samples.index]
    kinds = ['nichd_btb02', 'genomics_subject02', 'genomics_sample03']
    def do_one_manifest(kind):
        m = [row[kind] for row in lomanifests]
        m = pd.concat(m)
        csvpath = target_dir + os.path.sep + prefix + kinds[kind] + '.csv'
        templ_synid = manifest_template_synids[kinds[kind]]
        templ_df, templ_syn = get_manifest(templ_synid, syn, download_dir=target_dir)
        templ_p = target_dir + os.path.sep + templ_syn.properties.name
        write_manifest(m, templ_p, csvpath)
        return(m)
    manifests = tuple(do_one_manifest(k) for k in range(len(kinds)))
    return(manifests)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('slistcsv', help='sample list CSV file')
    parser.add_argument('-d', '--targetdir', help='target directory', default='.')
    parser.add_argument('-p', '--prefix', help='prefix for output files', default='chess-')
    args = parser.parse_args()
    make_manifests_main(slistcsv=args.slistcsv, target_dir=args.targetdir, prefix=args.prefix)
