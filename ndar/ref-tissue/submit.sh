#! /usr/bin/env bash

d=/home/attila/projects/bsm/ndar/ref-tissue
f1=$d/nichd_btb02_template.csv
f2=$d/genomics_subject02_template.csv
f3=$d/genomics_sample03_template.csv

# Oracle 8
# /usr/lib/jvm/java-8-oracle/bin/java \
# OpenJDK
# java \
/usr/lib/jvm/java-8-oracle/bin/java \
    -jar /opt/validationtool/current/validationtool-console-client-1.10.2.jar \
    --validationAPI https://ndar.nih.gov/api/validationtool/v2 \
    --collectionId C2458 \
    --title mssm-reftissue \
    --description "FASTQ data files for Reference Tissue Project" \
    --username andrewchess \
    --password Bern1e2017 \
    $f1 $f2 $f3
