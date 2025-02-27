# Pipeline_project_COMP383

## Description

## Problem 1: Retrieving Data

The links are from the SRA database under the data access tab

To retrieve Donor 1 (2dpi): "wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030"

To retrieve Donor 1 (6dpi): "wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033"

Then use the command fasterq-dump to convert the SRA files to paired-end fastq files

Donor 1 (2dpi) : 'fasterq-dump SRR5660030'

Donor 1 (6dpi): 'fasterq-dump SRR5660033'

## To run the Script 
clone the repo with 

'git clone https://github.com/evumana/Pipeline_project_COMP383.git'

The repo contains 4 samples of test data and a python script.

Move all 4 samples of test data and the python script to your home.

Then run the command 

"python3 pipeline_proj_final.py -i sampledata_1.fastq sampledata_2.fastq sampledata2_1.fastq sampledata2_2.fastq"

##Dependencies 

-Biopython

-BLAST+

-Bowtie 2

-SPAdes

-NCBI datasets

-SRA toolkit 
