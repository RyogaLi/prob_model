#Predict probabilities of cancer mutations 

## Input
**Please goto `./src/conf.py` to change the paths to input files**

##### Mandatory input files

1. `VCF file` 
    * File name should follow `tumour_id.vcf`
    * The INFO field should contain `VAF`
2. `tumour_info.csv`

    Example:
    ```
    sample_name     tumour_type     chromatin_profile   ...
    tumour_1        BRCA            ENCFF001UVV         ...
    ```
    * The first three columns are mandatory. 
    * `tumour_id` in the spreadsheet should be the same as the tumour_id in `tumour_id.vcf` file.
    * `chromatin_profile`: the profile to use for the input vcf file
    
3. `./features/`: a directory that contains the following files
    * hg19
    * 