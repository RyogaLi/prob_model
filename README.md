# Predict probabilities of cancer mutations 

## Dependencies

## Input

1. `VCF file` 
    * File name should follow `tumour_id.vcf`
    * The INFO field should contain `VAF`
2. `tumour_info.csv` - more details in **Step 1**-

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
    * concatenate all the fastq files into one using the following commands
        ```
        gunzip ./data/hg19/*.fa.gz
        cat ./data/hg19/*.fa > ./data/hg19.fa
        ```
* `mRNA` annotation file for hg19 can be downloaded from [UCSC annotated databse](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/all_mrna.txt.gz)
* `Signature` file here is provided by [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures), other signature files with same format can also be used 
* Save all the feature files into `./data/`

#### Step 1: Chromatin profile (OPTIONAL)

* It is optional to include the chromatin profile for each tumour.
* Please specify the ID of the `.bed` file for each tumour in `tumour_info.csv`
* if no chromatin profile is provided for a tumour, plase put `N/A` in the `chromatin_profile` column

#### Step 2: Change all the paths and settings in `conf.py`



## Other functionalities
* Simulate variants [See ./simulate_data/]()
* Analysis of output data


## Test
* S

