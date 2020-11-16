# farm-tools
A collection of scripts for use on Sanger HPC to manage downloading, relaligning, and pre-processing


## irods_download_remap_pipeline.sh
A script to download from IRODS data using LSF scheduling.
Following download it will extract reads from downloaded cram files using `samtools fastq` then align those
fastq to hg19 reference genome using either `BWA mem` or `star` depending on if provided with argument 'RNA' or 'DNA'.

e.g. to download and process RNA data from lane 1 of run 12345 
```
irods_download_remap_pipeline.sh 12345 1 RNA
```
