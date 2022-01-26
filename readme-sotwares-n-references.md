

# Trouble-shooting

### Hard-coded file paths & folders
1. Hard-coded file paths in `rnaseq-driver.sh`, `rnaesq-main.nf`, and `nextflow.config` are ok.
2. Check file paths in `rnaseq-driver.sh`
	* `scripts_root`: Root directory for script versions (see script version section).  
	* `singularity_container_rscript`: Container used by R-script samplesheet check.
	* `project_root`: Root directory for where nextflow is run, configs and samplesheets are saved etc. This directory is currently synced ("backupped") on the ldb server.
	* `delivery_root`: Root directory for where output directed to and that will be delivered to customer.
	* `ctg_save_root`: Root directory where to save qc data for ctg, i.e. same as `ctg-qc`.
3.  Check file paths in `nextflow.config`. Note that the different profiles have different containers and/or reference files. Make sure that the software versions as installed in .sif are compatible with references defined in nextflow.config, e.g. STAR indexed references.
	* Genome References
	* Singularity containers

\


## demultiplexing: bcl2fastq
[bcl2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf): Converts raw basecalls to fastq, and demultiplex samples based on indexes.  

### notes on bcl2fastq
**demux in ctg-rnaseq pipeline**
In versions ≤2.1.x, demultiplexing is included whithin the main nextflow pipeline.

**the -r -p and -w flags**<br>
The most demanding step is the processing step (-p option). Assign this step the most threads.
u The reading and writing stages are simple and do not need many threads. This consideration is important
for a local hard drive. Too many threads cause too many parallel read-write actions and suboptimal
performance.
**Character limitations**<brr>
- [Data]: The Sample_Project, Sample_ID, and Sample_Name columns accept alphanumeric characters, hyphens (-),
and underscores (_).
- [Header]. Make sure htat no regional characters are used. Other problematic signs are the + sign.


