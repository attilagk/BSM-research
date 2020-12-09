#! /usr/bin/env bash
# recode input checksum from base128 to base64
insum=$1
echo -n $insum | xxd -r -p -l 16 | base64
