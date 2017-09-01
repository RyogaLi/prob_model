#!/usr/bin/env bash

mRNA="/home/ryogali/data/all_mrna_hg19.txt"
ALEX="/home/ryogali/data/PCAWG_signature_patterns_beta.csv"
HG="/home/ryogali/data/hg19.fa"
TRI="/home/ryogali/data/trinucleotide.txt"
CHROMATIN="/home/ryogali/data/chromatin_data/esophagus-ENCFF430JQX.bed"

python ./src/read_files.py -m $mRNA -c $CHROMATIN -s $ALEX -tri $TRI -hg $HG
