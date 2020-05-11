#! /usr/bin/env python3

# Split gzipped fastq file into two halves: one for odd reads and another for
# even reads.

import gzip
import re
import os.path

def main(inpath):
    inpath = os.path.normpath(inpath)
    if re.match('.*\.fastq\.gz', inpath):
        ext = '.fastq.gz'
    elif re.match('.*\.fq\.gz', inpath):
        ext = '.fq.gz'
    else:
        raise ValueError('Extension must be .fastq.gz or .fq.gz')
    bn = inpath.replace(ext, '')
    outpathA = bn + '-A' + ext
    outpathB = bn + '-B' + ext
    with gzip.open(inpath, mode='rb') as infile:
        with gzip.open(outpathA, mode='wb') as outfileA:
            with gzip.open(outpathB, mode='wb') as outfileB:
                for i, line in enumerate(infile):
                    if i % 8 < 4:
                        outfileA.write(line)
                    else:
                        outfileB.write(line)


if __name__ == "__main__":
    import sys
    main(str(sys.argv[1]))
