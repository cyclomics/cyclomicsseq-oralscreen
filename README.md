# CycloSeq

Pipelines to process Cyclomics data.

This pipeline uses prior information from the backbone to increase the effectiveness of consensus calling from the circular DNA protocol by Cyclomics.  

## Dependencies

- Nextflow
- Docker
- Acces to the Github repo and a valid PAT token 

### data requirements

- data output by ONT Guppy-hac

## Usage

We assume that you have docker and nextflow installed on your system, if so running the pipeline is easy. You can run the pipeline directly from this repo, or pull it yourself and point nextflow towards it.

If you want to run the pipeline directly from github you need to use a Personal Access Token (PAT) as the password. Click the link to see how to create a [Personal Access Token (PAT)](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). Your PAT should have at least read permissions.
an example whould be:

```bash
nextflow run cyclomics/cycloseq -user <github_username> -r <current_version> -profile docker ...
```

### singularity

If docker is not an option, singularity(or Apptainer as it is called since Q2 2022) is a good alternative that does not require root access and therefor used in shared compute environments.



## Roadmap / Todo:
 ### Functionality:
 1. Add Post_Qc jobs with plots and tables
 
 ### Usability:
 1. Implement warning when the backbone is not present in the data
 1. Improve output folder structure
 1. Remove unused flags
 1. give warning when user sets backbone and its not used
 1. Make the read discovery smarter

### Improvements:
1. Change minimum vaf based on data available


## changelog

## 0.5.2
- Update Cycas version
- set --min-repeat-count default to 5



## 0.5.0
- Make --minimum-repeat-count parameter settable by user
- resolve tag conflict in both x and y tags
- Add csv output to bam accuracy.
- Increase decimal count in vcf to 6.
- Increased sensitivity of variant calling.

## 0.4.6
- Fix annotation file grouping
- Fix region file code
- Add option to select reads based on repeat count (YM bam tag)

### 0.4.5
- Added Y tags to final bam
- Improved QC

### 0.4.4
- Fix issues with bam accuracy plotting

### 0.4.3
- Fix fontsizes

### 0.4.2
- Added many QC steps with plots,
- Removed MionIONQC
- changed default variant calling to validate
- added auto detect of region of interest

### 0.4.1
- changed default read_pattern regex to detect rebasecalled sequencing runs automatically when pointed at the output folder.

### 0.4.0
- Updated Cycas to prevent runtime error with BB41
- Added variant validation optionality.
- Added quick_results flag for a glance of the results in the terminal. 
- Added profiles for conda and singularity support.

### 0.3.1
- Updated cycas version
- disable reporting due to to metadata incompatibility

### 0.3.0
- Added Cycas as the main consensus caller
- Added filtering
- Changed to Varscan variant calling


### Flag descriptions

|flag                           | info  |
|-------------------------------|---|
|--input_read_dir               | Directory where the output of Guppy is located, e.g.: "/data/guppy/exp001/".|
|--read_pattern                 | Regex pattern to look for fastq's in the read directory, defaults to: "fastq_pass/**.{fq,fastq}".|
|--sequencing_quality_summary   | Regex pattern for the summary file, default: "sequencing_summary*.txt".|
|--backbone_fasta               | Path to the fasta file containing backbone sequences.|
|--backbone_name                | Name of the sequence to extract from the backbone fasta excluding the starting ">", e.g.:"BB22". |
|--reference | Path to the reference genome to use, will ingest all index files in the same directory.|
|--output_dir | Directory path where the results, including intermediate files, are stored. |

### Example

```bash
nextflow run cyclomics/cycloseq -user <github_username> -r <current_version> -resume --input_read_dir '/some/path' --read_pattern 'fastq_pass/*.fastq' --output_dir 'testing'
```


### Running on A SLURM cluster such as UMCU HPC

login to the hpc using SSH. there start a sjob with:

```bash
srun --job-name "InteractiveJob" --cpus-per-task 16 --mem=32G --gres=tmpspace:450G --time 24:00:00 --pty bash
```

go to the right project folder

``` bash
cd /hpc/compgen/projects/cyclomics/cycloseq/pipelines/cycloseq/
```


## Developer notes

### Cycas addition to the repo

Cycas was added as a subtree using code from: https://gist.github.com/SKempin/b7857a6ff6bddb05717cc17a44091202.
This was done instead of submodule to make pulling of the repo easier for endusers and to stay compatible with nextflow run <remote> functionallity.


more specifically:
``` bash
git subtree add --prefix Cycas https://github.com/cyclomics/Cycas 0.4.3 --squash
```

To update run
``` bash
git subtree pull --prefix Cycas https://github.com/cyclomics/Cycas <tag> --squash
```
