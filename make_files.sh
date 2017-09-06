#!/usr/bin/env bash
# download and extract hg19 reference
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ ./data/hg19/
gunzip ./data/hg19/*
cat ./data/hg19/*.fa > ./data/hg19.fa
rm -rf ./data/hg19/


python ./src/read_files.py