#! /usr/bin/env bash

usage="usage: `basename $0` [threads_n_1.stderr threads_n_2.stderr ...]\ne.g.\nusage: `basename $0` 1MB/muscle-NeuN_pl-mutect2-?threads/out.stderr"
test $# -eq 0 && { echo -e $usage 2>&1; exit 1; }

echo "nthreads,real.time"
for f; do
    sed 's/.*\([[:digit:]]\+\)threads.*/\1,/' <<< $f | tr -d '\n'
    sed -n '/^.*Total runtime\s*\([0-9.]\+\)\s*secs.*/ {s//\1/; p}' $f
done
