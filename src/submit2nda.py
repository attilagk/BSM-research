#! /usr/bin/env python3

import subprocess

def submit(btb, gsub, gsam, title, description, build=False, user='andrewchess', pw='Bern1e2017',
        collection='2965', lib='/projects/bsm/', alternateEndpoint='BSMN-S3'):
    '''
    Validate and submit data to NDA

    Arguments
    btb: brain and tissue bank file
    gsub: genomic subjects file
    gsam: genomic samples file
    title: submission title
    description: submission description
    build: whether to build package
    user: NDA user name
    pw: NDA password
    collection: NDA collection
    lib: data file directory
    alternateEndpoint: alternate endpoint

    Value: the process object
    '''
    args = ['/home/attila/.local/bin/vtcmd', btb, gsub, gsam]
    args += ['-u', user, '-p', pw, '-c', collection, '-l', lib, '-a', alternateEndpoint, '-t', title, '-d', description]
    if build:
        args += ['-b']
    proc = subprocess.run(args, capture_output=True)
    return(proc)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('btb', help='brain and tissue bank file')
    parser.add_argument('gsub', help='genomic subjects file')
    parser.add_argument('gsam', help='genomic samples file')
    parser.add_argument('-t', '--title', help='submission title',
            default='Chess lab data')
    parser.add_argument('-d', '--description', help='description',
            default='Illumina reads in FASTQ and BAM')
    parser.add_argument('-u', '--user', help='NDA user name',
            default='andrewchess')
    parser.add_argument('-p', '--password', help='password',
            default='Bern1e2017')
    parser.add_argument('-c', '--collection', help='NDA collection',
            default='2965')
    parser.add_argument('-l', '--lib', help='data file directory',
            default='/projects/bsm/')
    parser.add_argument('-a', '--alternateEndpoint', help='alternate Endpoint',
            default='BSMN-S3')
    parser.add_argument('-b', '--build', help='build package',
            action="store_true")
    args = parser.parse_args()
    proc = submit(btb=args.btb, gsub=args.gsub, gsam=args.gsam,
            title=args.title, description=args.description, build=args.build, user=args.user, pw=args.password,
            collection=args.collection, lib=args.lib, alternateEndpoint=args.alternateEndpoint)
    print(proc.stdout.decode('utf-8'))
