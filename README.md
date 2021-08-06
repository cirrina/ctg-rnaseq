# ctg-rnaseq v1.0
Pipeline for demultiplexing, qc, alignment and summarization for RNAseq Illumina sequencing data.
The script can only process samples that are run in different sequencing runs. These have top be processed separately.
Different librarties within a pooled run should be processed separately.


## Installation lsens
Version specific pipeline folders are installed on ls4 at: `/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/`.
**Note** When running `rnaseq-primer`the entire directory with executables and configs are copied to (and run from within) the project folder. Each time the `rnaseq-primer` is run, all script files (including the `nextflow.config` will be overwritten)
**Note** Avoid adding the ctg-rnaseq script directories to $PATH. Instead run the `rnaseq-primer` script usging full path (thus allowing proper vresion control).

## Requirements
 1. Clone and build the **Singularity container** for this pipeline, e.g. https://github.com/cirrina/ctg-singularity/tree/main/rnaseq/v1.1). Add the correct path to the .sif -`singularity_container = ` parameter in `rnaseq-primer` shell sctipt.  
 2. Make SURE that the `scriptsdir`file location, hardcoded in the `rnaseq-primer` shell script, are valid.
 2. Edit your samplesheet to fullill all requirements. See section `SampleSheet` below
 3. Edit the nextflow.config so that file paths are correct. Make sure that the software versions as installed in .sif are compatible with references defined in nextflow.config, e.g. STAR indexed references.


## Run pipeline using Primer & Driver
.. 1-3 from above  
 4. Run the `naseq-primer`from *within the Illumina Sequencing Runfolder*. A project_id -i flag must be supplied as well as a SampleSheet located within the current dir.
 5. Optional: Modify parameters in the `nextflow.params.XXX` file. Individual settings can be  
 6. Run the `rnaseq-driver` from  *within the Project directory generated by the primer*. The rnaseq-driver executable should had been copied into the project folder by the rnaseq-primer script. The -i flag (project_id) must be supplied and match the proiject_id work folder where the script is run.


Example:
```
bash /projects/fs1/shared/ctg-pipelines/ctg-rnaseq/v1.0/rnaseq-primer -i 2021_070 -s 2021_070_SampleSheet.csv
cd /projects/fs1/shared/ctg-projects/rnaseq/2021_070
rnaseq-driver -i 2020_070
```

## Usage Seqonly-Driver
```
usage: prime-ctg-rnaseq [ -i projectid ] [ -s samplesheet ] [ -l force_lane ] [ -l force_index ] [ -h help ]

project_id            -i : Define 'project id'. Typically 202X-XXX. This will define the runfolder, but also define some output folders, such as blc2fastq2 output folders." as suffix
samplesheet           -s : IEM style laboratory samplesheet for this run. Within runfolder. (Default: CTG_SampleSheet.csv)
force_lane            -l : Set to 1 or 2 if to force demux to only run one lane. This if lane divider is used AND if lane is NOT specified in sample sheet. This parameter will overrid the Lane column in sample sheet'
force_index:  -f : Set this flag if to force iem samplesheet Rscript to overwrite index id and sequence informations based on index well id.
help          -h : print help message
```


## Input
- `--project_id`: Required. Typically the project id assigned by CTG LIMS, e.g. 2021_071. Will define name of output files folders that are generated by the pipeline.
- `--samplesheet` (see `SampleSheet` section below)



## Output:
* Pipeeline work folder.
  + `project_root`+`project_id` e.g. `/projects/fs1/shared/ctg-projects/rnaseq/2021_070`.
  + Work folder. Used while pipeline is running. All deliverables and ctg specific logs and qc metrics should be copied from this directory upon pipeline completion and should be safe to delete.
    - `./.nextflow/`. Temporary files used by nextflow pipeline. Will not be archived. Used for debugging faulty runs together with the `/.nextflow.log` and `/nf.log.rnaseq`
    - `./bin/` Executables used by pipeline copied from the `ctg-rnaseq` archive.
    - `./nf-output/` Work foilder in which output from all analyses are (temporarily) written. Files are moved/copied to delivery and ctg archive folders in the last step pipelie execution. Some analysis folders will not be kepot - bull *all* analyses are included when running the `ctg-multiqc`.

* Customer delvery folder
  + primary deliverable. This is the folder that should be delivered to customer. No other file or folder should be needed. In special cases, the outputs such as the
  + `delivery_root`+`project_id` e.g. `/projects/fs1/nas-sync/ctg-delivery/rnaseq/2021_070`.
    - `fastq`
    - `featurecounts`
    - `fastqc`
    - `md5sum.txt`
    - `multiqc` Light version of multiqc for customer delivery (not the same as ctg-multiqc that is more extensive)

* CTG archive folder
  + Here pipeline log files, ececutables, configs and samplesheets are saved together with relevant QC files, such as ctg-multiqc and fastqc.
  + `ctg_save_root`+`project_id`, e.g. `/projects/fs1/shared/ctg-qc/rnaseq/2021_070`
    - `./scripts/nextflow.params.project_id`
    - `iem.rscript.log`
    - `./multiqc-ctg/`. The multiqc-ctg analysis is run on all sequencing and project specific files and folders, including the Illumina runfolder and full bcl2fastq demultiplexing folder. Thus includes more information than the multiqc sent to customer.
  + **note** as of v1.0 not all relevant files are copied/transfered.




## Detailed Pipeline steps:

1. **rnaseq-primer**
  a. Create **work folder** based on the `project_id` flag (-i) and the `project_root` parameter. The pipeline executables and config files are copied and run from whitin the project directory.
  b. Run Rscript `iem-samplesheet-processor.R` to validate and generate modified SampleSheets for downstream analyses.
    + Input: IEM style SampleSheet with some modifications (see SampleSheet below).
    + Output: Tow sample sheets (`SampleSheet-2021_070-ctg.csv` and `SampleSheet-2021_070-demux.csv`) and a logfile (`iem.rscript.log`) that is used by `rnaseq-primer` to generate parameters file `nextflow.params.2021_070` used for the main nextflow sctipt.

    + illegal characters are removed.
    + Column settings are cross checked against database
    + Indexes are cross-checked and replaced if needed
  c. output: nextflow parameters file: `nextflow.params.2021_070`

2. **rnaseq-driver & nextflow-main**
Write run & project specific parameters: `nextflow.params.XXX`
* `Demultiplexing` (bcl2fastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).

* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report (https://multiqc.info/)
* `md5sum`: md5sum of all fastq files



### Samplesheet requirements:
Shpuld be in the foremat of Illumina IEM sample sheet.
The Bold/Italic field below must be correct!
Note that `Species` and `Pooled` are expected and must be added to a sheet generated by the IEM

[Header]
IEMFileVersion,5  
Investigator Name,X  
Experiment Name,X  
Date,YYYY-MM-DD  
Workflow,GenerateFASTQ  
Application,NovaSeq FASTQ Only  
Instrument Type,NovaSeq  
Assay,Nextera XT  
Index Adapters,"Nextera XT v2 Index Kit A"  // Must be specified according
Chemistry,Amplicon
*Species*,Homo sapiens      // curremtly Homo sapiens, Rattus norwegicus or Mus musculus
*Pooled*,true               // IMPORTANT: true or false.

[Reads]  
***26***  
***78***  

[Settings]  
Adapter,***CTGTCTCTTATACACATCT***  

[Data]  
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description  
***S1***,***S1***,,,***N702***,***CGTACTAG***,,,2021_024,  
***S2***,***S2***,,,***N706***,***TAGGCATG***,,,2021_024,  


### Source files


* ./bin
The bin directory is copied into the project directory. The bin contains executables (`iem-samplesheet-processor.R`)
  + `/bin/checklist-iem.csv` :
  + `/bin/checklist-index.csv` :
  * `/bin/iem-samplesheet-processor.R` :



## Important Folders & Variables that are hard-coded (should be moved to .config)
+ `rnaseq-primer` shell sctipt
  + `scriptsdir` : Must match the script file name (version folder)
  + `project_root` : Location where project work folders are generated. This folder is intermittant and used only for analyses. Important files (delivery and ctg save) will be moved from this folder. This foldeer should be safe to delete upon a successful run. (Default location in wich delivery folder will be created default '/projects/fs1/shared/ctg-projects/rnaseq')
  + `delivery_root` : Location where delivery folder will be created '/projects/fs1/nas-sync/ctg-delivery/rnaseq'
  + `ctg_save_root` : Location where files to be archived by ctg will be saved  qc saved for ctg will be created
  + `projectdir` : generate project dir from input -i `project_id` parameter and `project_root`



### Container  
https://github.com/perllb/ctg-seqonly/tree/main/container  
