#!/bin/bash

GTF=/juno/depot/annotation/M.musculus/ensembl/v83/Mus_musculus.GRCm38.83.gtf

controlFile=$1
targetFile=$2
region=$3

rtag=$(echo $region  | tr ':' '_')

python SparK/SparK.py \
    -pr $region \
    -tf $(cat $targetFile) \
    -cf $(cat $controlFile) \
    -o spark_average_${rtag} \
    -l samples negative_controls \
    -f FF0000 ABB2B9 \
    -gtf $GTF -ps averages -gs yes

