# ctg-rnaseq v1.2

Pipeline for demultiplexing, qc, alignment and transcript summarization for RNAseq Illumina sequencing data.
The pipeline is designed to handle multiple different RNAseq Assays and Species. Different assays will require differences in read strandness, read trimming etc.     

**Note:** The script can only process samples that are run in one single sequencing run (one Illumina Runfolder). If a project uses multiple sequencing runs, these have top be processed separately.  
**Note** The `project_id` (-i flag) will owerwrite the `Sample_Project` column in sample sheet - again only **one project** is allowed per demux/pipeline run.  
**Note:** Different projects/librarties within a pooled run must be processed separately.  


## Installation lsens
Version specific pipeline folders are installed on ls4 at: `/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/`

**Note:** When running `rnaseq-primer`the entire directory with executables and configs are copied to (and run from within) the project folder. Each time the `rnaseq-primer` is run, all script files (including the `nextflow.config` will be overwritten).  
**Note:** Do **not add** the ctg-rnaseq script directories to **$PATH**. Instead run the `rnaseq-primer` script usging full path - thus allowing proper version control.  

## Running the ctf-rnaseq pipeline
 1. Clone and build the **Singularity container** for this pipeline, e.g. https://github.com/cirrina/ctg-singularity/tree/main/rnaseq/v1.1). Add the correct path to the .sif -`singularity_container = ` parameter in `rnaseq-primer` bash sctipt.
 2. Make sure that the `scriptsdir`file location, hardcoded in the `rnaseq-primer` shell script, are valid and matches the current version.
 2. Edit your samplesheet to fullill all requirements. See section `SampleSheet` below
 3. Run the `naseq-primer`from *within the Illumina Sequencing Runfolder*. A project_id -i flag must be supplied as well as a SampleSheet located within the current dir.
 4. Optional: Modify parameters in the `nextflow.params.XXX` or the `nextflow.config` files. Make sure that the software versions as installed in .sif are compatible with references defined in nextflow.config, e.g. STAR indexed references.
 5. Run the `rnaseq-driver` from  *within the Project directory generated by the primer*. The rnaseq-driver executable should had been copied into the project folder by the rnaseq-primer script. The -i flag (project_id) must be supplied and match the proiject_id work folder in which you execute the the driver script.  


Example:
```
bash /projects/fs1/shared/ctg-pipelines/ctg-rnaseq/v1.0/rnaseq-primer \
  -i 2021_070 \
  -s 2021_070_SampleSheet.csv
cd /projects/fs1/shared/ctg-projects/rnaseq/2021_070
rnaseq-driver -i 2020_070
```

## Usage Seqonly-Driver
```
usage: prime-ctg-rnaseq [ -i projectid ] [ -s samplesheet ] [ -l force_lane ] [ -l force_index ] [ -h help ]

project_id, -i : Define 'project id'. Typically 202X-XXX. This will define the runfolder, but also define some output folders, such as blc2fastq2 output folders." as suffix
samplesheet, -s : IEM style laboratory samplesheet for this run. Within runfolder. (Default: CTG_SampleSheet.csv)
force_lane, -l : Optional. Set to 1 or 2 if to force demux to only run one lane. This if lane divider is used AND if lane is NOT specified in sample sheet. This parameter will overrid the Lane column in sample sheet'
force_index:  -f : Optional. Set this flag if to force iem samplesheet Rscript to overwrite index id and sequence informations based on index well id.
help  -h : print help message
```


## Input
- `--project_id`, `-i` : **Required**. Typically, the project id assigned by CTG LIMS, *e.g.* `2021_071`. Will form the basis for naming nomeclature when generating output files and folders through the pipeline.
- `--samplesheet`, `-s`:  **Required**. (see `SampleSheet` section below)



## Output:
In short, the pipeline will produce (1) customer delivery folder, contatining all relevant deliverables, including a customer-specific multiqc anlalysis (2) ctg archive folder containing an extended multiqc analysis together with logfiles, scripts and samplesheets making it possible to replicate the analysis, and (3) a  workfolder used during pipeline excecution.

1. **Pipeline work folder**  
`project_root`+`project_id` e.g. `/projects/fs1/shared/ctg-projects/rnaseq/2021_070`.  
Temporary work folder - used while pipeline is running. Upon a successful pipeline execution, all deliverables and ctg specific logs and qc metrics should have been copied from this directory, i.e. this folder can be  safely deleted after delivery.

  - `./.nextflow/`. Temporary files used by nextflow pipeline. Will not be archived. Used for debugging faulty runs together with the `/.nextflow.log` and `/nf.log.rnaseq` files.
  - `./bin/` Executables used by pipeline cloned from the `ctg-rnaseq` version specific archive, e.g. `/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/v1.0/bin/`.  
  - `./nf-output/`: Temporary folder to which all analyses are written. Files are moved/copied to delivery and ctg archive folders in the last step of pipelie execution. **All analyses** written to this directory are included when running the `ctg-multiqc`. Not all analyses files are archived though, such as alignmnent files used by fastqscreen analysis.  


2. **Customer delvery folder**  
`delivery_root`+`project_id` e.g. `/projects/fs1/nas-sync/ctg-delivery/rnaseq/2021_070`.  
Primary delivery folder. This is the folder that should be transfered to the delivery server for customer delivery. In a standard pipeline run, no additional files or folders should need to be added.

    - `./fastq/` : bcl2fastq output with fastq files in project id folder (stemming from the project_id column in sample sheet wich is forced to pipeline project id)`./fastq/2021_070` If `pooled, false` (thus a seq run with non pooled data) is specified in sample sheet, then this folder will also incude `./fastq/Undetermined*.fastq` files as well as `./fastq/Reports` and and
    - `./featurecounts/`: transcript summarization matrix, e.g. `2021_070_geneid.featureCounts.txt`
    - `./fastqc/` : fastqc analysis output. One fastqc output per sample and read. e.g. `R09_T1_S22_R1_001_fastqc.html` and `R09_T1_S22_R1_001_fastqc.zip`
    - `./star/` : star alignment output, e.g. `R46_T3_Aligned.sortedByCoord.out.bam` and `R46_T4_Aligned.sortedByCoord.out.bam.csi`
    - `./multiqc/` Lighter version of multiqc for customer delivery (not the same as `ctg-multiqc` that is more extensive). This is basically only run on the contents present in delivery folder (including fastq, fastqc, star and featurecounts, but does not include anlalyses such as fastqscreen, picard, or illumina InterOp)
    - `./md5sum.txt`: md5sum on all files in the customer delivery folder.

3. **CTG archive folder (ctg-qc)**  
`ctg_save_root`+`project_id`, e.g. `/projects/fs1/shared/ctg-qc/rnaseq/2021_070`  
Here pipeline log files, executaables, configs and samplesheets are archived together with relevant QC files, such as ctg-multiqc and fastqc.
  - `./scripts/nextflow.params.project_id`
  - `iem.rscript.log`
  - `./multiqc-ctg/`. The multiqc-ctg analysis is run on all sequencing and project specific files and folders, including the Illumina runfolder and full bcl2fastq demultiplexing folder. Thus includes more information than the multiqc sent to customer.  

**Note:** as of v1.0 all relevant files are not yet copied/transfered.




## Detailed Pipeline steps:

1. **rnaseq-primer**  
  a. Create **work folder**  
  e.g. `project_root`+`project_id` e.g. `/projects/fs1/shared/ctg-projects/rnaseq/2021_070`  
  Based on the `project_id` flag (-i) and the `project_root` parameter.
  he pipeline executables (e.g. `rnaseq-driver`) and config files and scripts (e.g. `nextflow.config` and `rnaseq-main.nf`) are copied and run from whitin the project directory **not** from their primary location on lsens.

  b. Run Rscript `iem-samplesheet-processor.R`  
    To validate and generate modified SampleSheets for downstream analyses.  
    **Input:**  
    Illumina Experiment Manager (IEM) style SampleSheet *modified to fit CTG LIMS* (see SampleSheet below).  
    **Output:**  
    `SampleSheet-2021_070-demux.csv`: used for bcl2fastq  
    `SampleSheet-2021_070-ctg.csv`: used by the main nextflow script. The rsctipt add Assay specific columns, such as, strandness, paired, and the file names for fastq and bam fils etc.  
    `iem.rscript.log`: logfile with paramaters that are imported by nextflow-primer and passed on to `nextflow.params.2021_070` used by the main nextflow script.  
    See `Source` section below for additonal script info.  

  c. rnaseq-primer output:
  nextflow parameters file: `nextflow.params.2021_070`.


2. **rnaseq-driver & nextflow-main**

This driver will initiate nextflow pipeline `nextflow-main` using two config files (`nextflow.config` and `nextflow.params.2021_070`) together with samplesheets (`SampleSheet-2021_070-ctg.csv` and `SampleSheet-2021_070-demux.csv`)  


i. `Demultiplexing`  
bcl2fastq: Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).  

ii. Check files. Chek if all expected fastq files have been generated. The fastq names are hardcoded into the ctg sample sheet (`SampleSheet-2021_070-ctg.csv`) columns `fastq_1`and `fastq_2` by the `iem-samplesheet-processor.R` rscript.  

ii. `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).  

iii. `multiQC`: Compile fastQC and cellranger count metrics in multiqc report (https://multiqc.info/)  


 `md5sum`: md5sum of all fastq files



### Samplesheet requirements:
Should be in the format of Illumina IEM sample sheet.  
The **Bold** fields below must be correctly specified!  
**Note on [Header] `Species` and `Pooled`:** These `[Header]` variables are not included by IEM and **must be added**.  
**Note on [Data] Sample_ID:** : Olny Sample_ID values are required. **Sample_Name** will be forced to the same value as Sample_ID. This to force the file structure output from bcl2fastq demux to be consistant. Sample_ID ≠ Sample_Name will produce additioonal folder structure.    
**Note [Data] Sample_Project:** This column will be force overwritten by the `iem-samplesheet-processor.R` script to the same project id value as defined when executing the pipeline.  

**Note on [Header] Instrument Type:**  Can be NovaSeq or NovaSeq1.0 protocol.  
**Note on [Header] Assay:**  Allowed Assay values are listed in `./bin/checklist-iem.csv`. Assays are specified using same nomenclature as within the IEM software.  
**Note on [Header] Index Adapters:**  Allowed Index Adapters values (and Assay combinations) are listed in `./bin/checklist-iem.csv`. Specified using same nomenclature as within the IEM software.    
**Note on [Settings] Adapter:**  Adapters will be cross-checked using the`./bin/checklist-iem.csv` file. The Adapter value must match the value specified under respective Index Adapter.
**Note on [Reads]:** This section is used to determine if the run is **paired or not**. Note that the actual read lenths are **probably** not used by bcl2fastq, pooling of different Assay types may complicate thisgs,

----
[Header]  
IEMFileVersion,5  
Investigator Name,X  
Experiment Name,X  
Date,YYYY-MM-DD  
Workflow,GenerateFASTQ  
Application,NovaSeq FASTQ Only  
**Instrument Type**,**NovaSeq**
**Assay**,**TruSeq Stranded mRNA**  
**Index Adapters**,**IDT-ILMN TruSeq RNA UD Indexes (96 Indexes)**  
Chemistry,Amplicon  
**Species**,**Homo sapiens**   
**Pooled**,**true**

[Reads]  
**101**    
**101**

[Settings]  
**Adapter**,**CTGTCTCTTATACACATCT**  

[Data]  
**Sample_ID**,Sample_Name,Sample_Plate,Sample_Well,**Index_Plate_Well**,**I7_Index_ID**,**index,I5_Index_ID**,**index2**,Sample_Project,Description  
**R41_C**,,CTGpl_073,A01,**A01**,**UDI0001**,**CCGCGGTT**,**UDI0001**,**CTAGCGCT**,2021_070,CTGpool_0158  
**R41_T1**,,CTGpl_073,B01,**B01**,**UDI0002**,**TTATAACC**,**UDI0002**,**TCGATATC**,2021_070,CTGpool_0158  

----

### Source files

#### ./bin folder
The `./bin/` contains executables used by the pipeline, The bin directory is cloned into the project work directory and files are accessed locally from within there.

1. **`iem-samplesheet-processor.R`**  
This script validate sample sheet, and to generate modified SampleSheets for downstream analyses.

  - **Input:**  
  Illumina Experiment Manager (IEM) style SampleSheet *modified to fit CTG LIMS* (see SampleSheet below).  

  - **Output:**  
  Two sample sheets (`SampleSheet-2021_070-ctg.csv` and `SampleSheet-2021_070-demux.csv`) and a logfile (`iem.rscript.log`).  
  The logfile is used by `rnaseq-primer` to generate parameters file `nextflow.params.2021_070` used for the main nextflow sctipt.

  The Rscript will:
  - check basic integrity of the IEM samplesheet,
  - check that required columns are present and correctly specified. This is done using the `./bin/checklist-iem.csv` file. A given Assay/Index Adapters combionation will have specific requirements, e.g. Read2StartFromCycle or Adapter, all specified in the `checklist-iem.csv`
  - check for illegal sample namimings and duplicated sample names
  - Replace illegal characters.
  - cross check index nmames and sequences against the `./bin/checklist-index.csv` file.
  - add fastq and bam file namings according to defined nomenclature. These are used by `nextflow-main` to chek if all expected files have been generated.  
  Options/flags:
  - Column settings are cross checked against database
  - Indexes are cross-checked and replaced if needed



2. **`/bin/checklist-iem.csv`**  
Used by the `iem-samplesheet-processor` rscript to cross check allowed values and value combinations given in the sample sheet, for example check that the Assay value is allowed, or that Species value is properly specified. Also, conditional rules can here be defined, such as requiring that one paricular Assay must also have Read2StartFromCycle a particular [Settings] defined.

The following columns are expected and are used to define what rules the rscript will use to cross check sample sheet:  
- iem_section
- parameter
- conditional_parameter
- conditional_value
- value
- Comment

Thus, the lines containing `iem_section` *[Header]* and `parameter` *Assay*, such as:  
`[Header],Assay,,,TruSeq Stranded mRNA,TruSeq Stranded mRNA Sample Preparation Guide (Part 15031047 Rev. E) Modification of Enrich DNA Fragment step:from 15 to 13 PCR cycles`
will be used to cross check that the *[Header]* `iem_section` contains one value (row) that specify `Assay` and that the specified value, `value`-column, is mathichg one of the specified, e.g. `TruSeq Stranded mRNA`

Conditions can be added using the `conditional_parameter` and `conditional_value` columns.  
Such as the line:
`[Settings],Read2StartFromCycle,Assay,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,4,`   
will be interpreted as:
IF the `conditional_parameter` *Assay* matches the conditional value `SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian`, the script will require that the `Read2StartFromCycle` parameter is given under `iem_section` *[Settings]* AND that it matches the `value` of *4*.  

3. **`/bin/checklist-index.csv`**  
Used by the `iem-samplesheet-processor` rscript to check that index id and well postions correspond to the sequence given in the samplesheet.  
csv-delimited file
Note that default NovaSeq Instrument type (specified as *Novaseq*) is the NovaSeq1.5 chemistry.  
The following columns are expected:
- Index_Adapters // Adapter kit, e.g. *IDT-Ilmn RNA UD Indexes SetA Ligation*
- product.no // product number, e.g. *20040553*
- Instrument_Type // Novaeq (default) or NovaSeq1.0. The 1.0 chemistry I5 indexes are reverese complement to I7.
- Index_Plate_Well //
- I7_Index_ID
- index
- I5_Index_ID
- index2


## Important Folders & Variables that are hard-coded (should be moved to .config)
+ `rnaseq-primer` shell sctipt
  + `scriptsdir` : Must match the script file name (version folder)
  + `project_root` : Location where project work folders are generated. This folder is intermittant and used only for analyses. Important files (delivery and ctg save) will be moved from this folder. This foldeer should be safe to delete upon a successful run. (Default location in wich delivery folder will be created default '/projects/fs1/shared/ctg-projects/rnaseq')
  + `delivery_root` : Location where delivery folder will be created '/projects/fs1/nas-sync/ctg-delivery/rnaseq'
  + `ctg_save_root` : Location where files to be archived by ctg will be saved  qc saved for ctg will be created
  + `projectdir` : generate project dir from input -i `project_id` parameter and `project_root`



### Container  
https://github.com/perllb/ctg-seqonly/tree/main/container  
