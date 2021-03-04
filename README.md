# COMP 383 MiniProject
## Transcriptome Analysis Pipeline for HCMV

From Wikipedia: 
"Human herpesvirus 5 is also known as Human cytomegalovirus and is typically abbreviated as HCMV. 
Although they may be found throughout the body, HCMV infections are frequently associated with the
salivary glands. HCMV infection is typically unnoticed in healthy people, but can be life-threatening for the immunocompromised, such as HIV-infected persons, organ transplant recipients, or newborn infants.
Congenital cytomegalovirus infection can lead to significant morbidity and even death. After infection, HCMV remains latent within the body throughout life and can be reactivated at any time. Eventually, it may cause mucoepidermoid carcinoma and possibly other malignancies such as prostate cancer."

Cheng et al. 2017 (https://www.ncbi.nlm.nih.gov/pubmed/29158406) produced the transcriptomes of HCMV post
infection. The data produced from this paper will be analyzed in this pipeline.

## Software Dependencies
This script must have the following software installed in the path:
* Python 3
* R v4.0.3
* SRA-toolkit
* Biopython
* kallisto
* sleuth
* Bowtie2
* Spades
* BLAST+

## Options
*  *-t*,*--test:* Use test data to run this pipeline. Will execute in <5 minutes.

### Running the pipeline ###
Run the following in the directory you'd like to clone the repository in:
```shell
git clone https://github.com/MatthewLoffredo/comp-383-miniproject-loffredo.git
```

This will download the repository. To run the pipeline with test data, use the command:

```shell
$ python3 HCMV.py -t
```

To run the pipeline with the full dataset, use the command:

```shell
$ python3 HCMV.py
```

### Output ###
After running the script, all files and subdirectories in the pipeline will be automatically written to a directory called `miniproject_Matt_Loffredo`. The desired pipeline output will be stored in that directory, in a file named `miniProject.log`.

An output of the pipeline using the full dataset is stored in `miniProject.log` in the main repository directory.
