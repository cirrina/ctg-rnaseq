

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


