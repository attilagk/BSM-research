#! /usr/bin/env bash

dirprefix=/projects/bsm/attila/ # work on ada
#dirprefix=~/projects/bsm/attila/
wd=$dirprefix/results/2018-12-18-verifyBamID
cd $wd

for id in syn2507223 syn2507168 syn2507166; do
    synapse get $id
done
