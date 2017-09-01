# Predict probabilities of cancer mutations 

### Dependencies

### Input
---
**Please goto `./src/conf.py` to change the paths to input files**

##### Mandatory input files

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


### Output
---




### RUN
---

#### Step 0: Preprocess the feature files 
* 


#### Step 1: Specify all the required file path in `conf.py`

#### Step 2: 


### Other functionalities
---


### Test
---

