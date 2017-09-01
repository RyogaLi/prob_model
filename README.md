# Predict probabilities of cancer mutations 

## Dependencies

## Input

1. `VCF file` 
    * File name should follow `tumour_id.vcf`
    * The INFO field should contain `VAF`
2. `tumour_info.csv`

    Example:
    ```
    tumour_id     tumour_type     chromatin_profile     ...
    tumour_1      BRCA            ENCFF001UVV         ...
    ```
    * The first three columns are mandatory. 
    * `tumour_id`: this should be the same as the title of `tumour_id.vcf`.
    * `chromatin_profile`: the profile to use for the input vcf file


## Output





## RUN


#### Step 0: Get all the required feature files 
* `hg19` file can be obtained from [UCSC genome browser](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/)
    * Use the following command to download all the files into ./data/hg19/
      
        ```
        rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ ./data/hg19/
        ```
    * concatenate all the fastq files into one using the following command
        ```
        gunzip ./data/*
        cat ./data/* > hg19.fa
        ```
* `mRNA` annotation file for hg19 can be downloaded from [UCSC annotated databse](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/all_mrna.txt.gz)
* `Chromatin profile` for different tumour types can be downloaded from [ENCODE project - DNase-seq database](https://www.encodeproject.org/matrix/?type=Experiment&status=released&assay_slims=DNA+accessibility&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&award.project=ENCODE)
* `Signature` file here is provided by [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures), other signature files with same format can also be used 
* Save all the feature files into `./data/`

#### Step 1: Specify all the required file path in `conf.py`

#### Step 2: 


## Other functionalities



## Test


