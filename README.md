# Predict probabilities of cancer mutations 

## Dependencies
- Here we use a universal mixture file to make the prediction. To compute a signature mixture file for each tumour, please goto : [Trackature](https://github.com/YuliaRubanova/Trackature)

## Input

1. `VCF file` 
    * File name should contain `tumour_id`
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


#### Step 0: Get all the required supplementary files 
* `mRNA: all_mrna.txt` annotation file for hg19 can be downloaded from [UCSC annotated databse](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/all_mrna.txt.gz)
* `Signature: alexSignatures.txt` file here is provided by [COSMIC](http://cancer.sanger.ac.uk/cosmic/signatures), other signature files with same format can also be used 
* `Chromatin profile: ID.bed` (**OPTIONAL**) for different tumour types can be downloaded from [ENCODE project - DNase-seq database](https://www.encodeproject.org/matrix/?type=Experiment&status=released&assay_slims=DNA+accessibility&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&award.project=ENCODE)
* Save all the supplementary files into `./data/`
* RUN: 

    ```
    ./make_file.sh
    ```


#### Step 1: Specify all the required file path in `conf.py`

#### Step 2: 


## Other functionalities



## Test


