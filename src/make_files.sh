#!/usr/bin/env bash

gunzip ./data/hg19/*
cat ./data/hg19/*.fa > ./data/hg19.fa


python ./src/read_files.py