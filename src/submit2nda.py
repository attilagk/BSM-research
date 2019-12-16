#! /usr/bin/env python3

import subprocess

def submit(btb, gsub, gsam, title, build=False, user='andrewchess', pw='Bern1e2017',
        collection='2965', lib='/projects/bsm/'):
    description = 'Unmapped and mapped reads from bulk sequencing; NIH U01MH106891; Chess Lab, Mount Sinai, New York'
    args = ['/home/attila/.local/bin/vtcmd', btb, gsub, gsam]
    args += ['-u', user, '-p', pw, '-c', collection, '-l', lib, '-t', title]
    if build:
        args += ['-b']
    proc = subprocess.run(args, capture_output=True)
    return(proc)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('btb', help='b tissue bank file')
    parser.add_argument('gsub', help='genomic subjects file')
    parser.add_argument('gsam', help='genomic samples file')
    parser.add_argument('-t', '--title', help='submission title',
            default='Chess lab data')
    parser.add_argument('-u', '--user', help='NDA user name',
            default='andrewchess')
    parser.add_argument('-p', '--password', help='password',
            default='Bern1e2017')
    parser.add_argument('-c', '--collection', help='NDA collection',
            default='2965')
    parser.add_argument('-l', '--lib', help='TNseq calling algorithm',
            default='/projects/bsm/')
    parser.add_argument('-b', '--build', help='build package')
    args = parser.parse_args()
    proc = submit(btb=args.btb, gsub=args.gsub, gsam=args.gsam,
            title=args.title, build=args.build, user=args.user, pw=args.password,
            collection=args.collection, lib=args.lib)
    print(proc.stdout.decode('utf-8'))
