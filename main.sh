#!/usr/bin/env bash
# v0.0
# output probabilities directory

# train files
#TRAINDIR="/home/q/qmorris/ryogali/pwgs/variants/breastcancer/bc_vcf_unfiltered_features/train/"
# test files
#TESTDIR="/home/q/qmorris/ryogali/pwgs/variants/breastcancer/bc_vcf_unfiltered_features/test/"
# TESTDIR="/home/q/qmorris/ryogali/pwgs/variants/breastcancer/bc_vcf_unfiltered_features/lowsupport/"
# signature exposure files
# SIGDIR="/home/q/qmorris/ryogali/data/mut/"
# mixture files
# MiXTUREDIR="/home/q/qmorris/ryogali/data/mixtures_cnvint_BRCA/"
# MiXTUREDIR="/scratch/q/qmorris/yulia/pwgs/samples/vaf/results_onlyKnownSignatures_vaf/BRCA/"
# PSUB="/home/q/qmorris/yulia/pwgs/ssmvaf/"
#LOWSUPPORTDIR="/home/q/qmorris/ryogali/pwgs/variants/breastcancer/bc_vcf_unfiltered_features/lowsupport/"
#python main.py -o $OUTPUT -train $TRAINDIR -test $TESTDIR -m $MiXTUREDIR -low $LOWSUPPORTDIR

# change the path below
OUTPUT="/home/ryogali/output/test_logging"
INPUT="/home/ryogali/data/consprelim"
MIXTUREDIR="/mnt/raisin/yulia/pwgs/samples/vaf/results_onlyKnownSignatures_vaf/BRCA"
FEATURE="/home/ryogali/data/features"
python ./src/main.py -o $OUTPUT -i $INPUT -m $MIXTUREDIR -f $FEATURE -c "True"
