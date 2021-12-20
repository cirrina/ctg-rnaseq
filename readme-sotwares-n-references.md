
# ctg-rnaseq
------------


Primary data processing pipeline for Illumina RNA-seq data produced at CTG. Built to run using Nextflow and (multiple) Singularity containers. The pipeline is designed to start from an Illumina Runfolder to demultiplex data, different qc measures, alignment and transcript summarization.

The pipeline is designed to handle multiple different RNAseq Assays (library preparation methods), and reference Species. Different assays will require differences in read strandness, read trimming etc. These input parameters are either set through the Sample Sheet or profile-specific and specified in the nextflow configurations file (see below).

**Note on Demultiplexing and Fastq:** As of versions 3.x and above, demultiplexing is no longer performed within the pipeline. Instead the pipeline assumes fastq files as starting input. Demultiplexing is performed on a RunFolder basis using the ctg-demux pipeline and may include multiple projects rather than one single project (as of versions ≤2.x. 

**Note on RunFolder(s):** Versions ≤2.x. are designed to process samples in **one single sequencing run** (all run within one Illumina Runfolder). If a project uses multiple sequencing runs, manual adjustments have to be made to fit the pipeline execution requirements (see below).  

**Note on ProjectId/Sample_Project** Only **one project** is allowed per pipeline run. `ProjectId` is supplied through the SampleSheet [Header]. Thus, for a run, all individual `Sample_Project` entries in the SampleSheet [Data] section must be same as `ProjectId` for all **samples**. 

The `project_id` (supplied by `ProjectId` in SampleSheet) will owerwrite the `Sample_Project` column in sample sheet - again  (As of v1.3 the `ProjedId` is supplied through the SampleSheet, **not** using -i flag as for versions 2.x. Thus, different projects (& library pools) within a pooled run must be processed separately.  

**Note on multiple Species:** The pipeline is designed for **one and the same species**. This is due to that transcript summarization i must be performed using the same refence GFF file. If multiple species are used, then run multiple pipeline runs, defined by unique project IDs.


**nextflow profiles {PipelineProfile}**
The pipeline is designed for three different profiles/main configurations.<br>
- rnaseq: Sequencing of a mRNA libraries generated at CTG. 
- rnaseq_total: Sequencing of total RNA libraries generated at CTG. 
- uroscan: The UroScanSeq pipeline. Processing and generation of *bladderreports* as part of the UroSCan project. 
 


# Running the ctg-rnaseq pipeline
The pipeline is intiated by executing a driver script `rnaseq-driver`. The diver will:

* Prime a work directory (copy scripts and configs to this dir)
* Check the SampleSheet (illegal characters, indexes etc). See SampleSheet below.
* Generate a run/project-specific nextflow config. 
* Initiate the `nextflow-main` piepline script.

### Quickstart guide

1. Check file paths in `rnaseq-driver`
	* `scripts_root`: Root directory for script versions (see script version section).  
	* `singularity_container_rscript`: Container used by R-script samplesheet check. 
	* `project_root`: Root directory for where nextflow is run, configs and samplesheets are saved etc. This directory is currently synced ("backupped") on the ldb server. 
	* `delivery_root`: Root directory for where output directed to and that will be delivered to customer.
	* `ctg_save_root`: Root directory where to save qc data for ctg, i.e. same as `ctg-qc`.
2.  Check file paths in `nextflow.config`. Note that the different profiles have different containers and/or reference files. Make sure that the software versions as installed in .sif are compatible with references defined in nextflow.config, e.g. STAR indexed references.
	* Genome References 
	* Singularity containers
3. Edit your samplesheet to fullfill all requirements. See section `SampleSheet`. A SampleSheet must be supplied and present within the pipeline execution dir.  
4. Note that **.fastq** and **.bam** file names must be **defined in the SampleSheet**.
4. Run the `rnaseq-driver` from `Illumina Sequencing Runfolder` (if first time & no project dir has been created)
5. **OR** run the `naseq-driver`from within the `Project work folder`. This requires that path to fastq directory is specified (`-f` flag). Typically this is performed when resuming a failed run or changing a projects paramters using the config files.
6. Optional: A project runfolder can be primed without starting the nextflow pipeline. Use the `-p`, *prime run* flag. This in order to modify parameters in the `nextflow.config.project.XXX` or the `nextflow.config` files. 


Example initiate from within Illumina runfolder (v2.1.x).:

```
/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/2.0.0/rnaseq-driver \
  -s 2021_024_SampleSheet-IEM.csv
```

Re-run a failed run. Done from within the project folder. Filepath to fastq diectory is needed :

```
rnaseq-driver \
  -s 2021_024_SampleSheet-demux.csv \
  -f nf-output/fastq/2021_024
```
Prime a work folder (but do not start a run)

```
/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/2.0.0/rnaseq-driver \
  -s 2021_024_SampleSheet-IEM.csv -p true
```

Resume a failed nextflow run:

```
rnaseq-driver \
  -s 2021_024_SampleSheet-demux.csv \
  -r
```

# Configurations and input parameters
Input parameters and variables are determined by:<br>

*  `SampleSheet`: contains most required input parameters.
*  `nextflow.config`: parameters that define the different nextflow profiles, including genome reference and container files.
*  `nextflow.config.{projectId}`: auto-generated from the rnaseq-driver using parameters from 'SampleSheet' and file paths from the driver file.
*  `rnaseq-driver`: contains file paths etc. 


## Nextflow Profiles
A nextflow config profile must be defined and is  set in `SampleSheet` (below) using `PipelineProfile ` variable. A profile is a set of nextflow configuration attributes that can be activated/chosen when launching a pipeline execution by using the `-profile` option. By convention the `standard` profile is implicitly used when no other profile is specified by the user (not implemeted here - a profile must be defined).

Available profiles:

* rnaseq
* rnaseq_total
* uroscan

**Note:** When using the profiles feature in your config file do **NOT** set attributes in the same scope both inside and outside a profiles context. e.g if using `process` attributes, all of these have to be defined within the different profiles sections. 

## Nextflow Configs
Two different nextflow.configs are used to pass variables to the nextflow pipeline. The project specific `nextflow.config.{projectId}` is auto generated and carries variables from the samoplesheet to the actual nextflow script. `nextflow.config.{projectId}` is initiated by nextflow through the `-c` flag and will override (if reduntant) the `nextflow config`. 


## SampleSheet
The SampleSheet is the main source of input parameters for the pipeline. The SampleSheet is structured in accordance to [Illumnina IEM software](https://support.illumina.com/sequencing/sequencing_software/experiment_manager/downloads.html ) requirements (and can therefore also be used as input for the bcl2fastq software). The ctg samplesheet will include a number of additional placeholders wihtin the [Header] section that are used as holders for metadata and is used by the ctg-pipelines.

- [bcl2fastq](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html) demultixing software compatible (IEM style csv samplesheet)
- Metadata variables added that are used by the ctg-rnaseq pipeline.
- Sample names, including fastq and bam file namings.



### Samplesheet requirements:
Should be in the format of Illumina IEM sample sheet. For a sample sheet with multiple projkects (as outputed from the lab) may include multiple 

**[Header]:** The file must start with a [Header] section. Whithin this section additional metadata can be added.<br>
**PipelineName:** In this case always "ctg-rnaseq".<br>
**PipelineVersion:** Should specify the pipeline version e.g. 2.1.5. This will define what scripts are dowloaded to the work directiory.<br>
**PipelineProfile:** Will determine what nextflow profile to initiate. Nextflow profile parameters are defined in the `nextflow.config`. Avalable profiles are: rnaseq, rnaseq_total or uroscan. See seciton on Nextflow profiles.<br>
**RunFolder:** The Illumina Runfolder Name. <br>
**ProjectId**,2021_145,,,,,,,,,,,,,,,,,,,,,,<br>
**PoolName**,CTGpool_0174,,,,,,,,,,,,,,,,,,,,,,<br>
**Instrument** Type,NovaSeq1.5,,,,,,,,,,,,,,,,,,,,,,<br>
**FlowCell**,"NovaSeq 6000 S1 Reagent Kit, 100 cycles v1.5 (20028319)",,,,,,,,,,,,,,,,,,,,,,<br>
**SharedFlowCell**,false,,,,,,,,,,,,,,,,,,,,,,<br>
**LaneDivider**,false,,,,,,,,,,,,,,,,,,,,,,<br>
**Species**,Other,,,,,,,,,,,,,,,,,,,,,,<br>
**Assay**,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,,,,,,,,,,,,,,,,,,,,,,<br>
**Index** Adapters,SMARTer RNA Unique Dual Index Kit -96U Set A,,,,,,,,,,,,,,,,,,,,,,<br>
**Strandness**,reverse,,,,,,,,,,,,,,,,,,,,,,<br>
**Paired**,false,,,,,,,,,,,,,,,,,,,,,,<br>

The **Bold** fields below must be correctly specified!  
**Note on [Header] `Species` and `Pooled`:** These `[Header]` variables are not included by IEM and **must be added**.  
**Note on [Data] Sample_ID:** : Olny Sample_ID values are required. **Sample_Name** will be forced to the same value as Sample_ID. This to force the file structure output from bcl2fastq demux to be consistant. Sample_ID ≠ Sample_Name will produce additioonal folder structure.    
**Note [Data] Sample_Project:** This column will be force overwritten by the `iem-samplesheet-processor.R` script to the same project id value as defined when executing the pipeline.  

**Note on [Header] Instrument Type:**  Can be NovaSeq or NovaSeq1.0 protocol.  
**Note on [Header] Assay:**  Allowed Assay values are listed in `./bin/checklist-iem.csv`. Assays are specified using same nomenclature as within the IEM software.  
**Note on [Header] Index Adapters:**  Allowed Index Adapters values (and Assay combinations) are listed in `./bin/checklist-iem.csv`. Specified using same nomenclature as within the IEM software.    
**Note on [Settings] Adapter:**  Adapters will be cross-checked using the`./bin/checklist-iem.csv` file. The Adapter value must match the value specified under respective Index Adapter.
**Note on [Reads]:** This section is used to determine if the run is **paired or not**. Note that the actual read lenths are **probably** not used by bcl2fastq, pooling of different Assay types may complicate thisgs,




### IEM samplesheet example:
Parameters in bold are used by the rnaseq pipeline:

**[Header],,,,,,,,,,,,,,,,,,,,,,,**<br>
**Pipeline**,ctg-rnaseq,,,,<br>
**PipelineVersion**,2.1.5,,,,,,,,,,,,,,,,,,,,,,<br>
**PipelineProfile**,rnaseq_total,,,,,,,,,,,,,,,,,,,,,,<br>
**RunFolder**,,,,,,,,,,,,,,,,,,,,,,,<br>
**ProjectId**,2021_145,,,,,,,,,,,,,,,,,,,,,,<br>
**PoolName**,CTGpool_0174,,,,,,,,,,,,,,,,,,,,,,<br>
**Instrument** Type,NovaSeq1.5,,,,,,,,,,,,,,,,,,,,,,<br>
**FlowCell**,"NovaSeq 6000 S1 Reagent Kit, 100 cycles v1.5 (20028319)",,,,,,,,,,,,,,,,,,,,,,<br>
**SharedFlowCell**,false,,,,,,,,,,,,,,,,,,,,,,<br>
**LaneDivider**,false,,,,,,,,,,,,,,,,,,,,,,<br>
**Species**,Other,,,,,,,,,,,,,,,,,,,,,,<br>
**Assay**,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,,,,,,,,,,,,,,,,,,,,,,<br>
**Index** Adapters,SMARTer RNA Unique Dual Index Kit -96U Set A,,,,,,,,,,,,,,,,,,,,,,<br>
**Strandness**,reverse,,,,,,,,,,,,,,,,,,,,,,<br>
**Paired**,false,,,,,,,,,,,,,,,,,,,,,,<br>
IEMFileVersion,5,,,,,,,,,,,,,,,,,,,,,,<br>
Experiment Name,NovaSeq_XP,,,,,,,,,,,,,,,,,,,,,,<br>
Date,,,,,,,,,,,,,,,,,,,,,,,<br>
Workflow,GenerateFASTQ,,,,,,,,,,,,,,,,,,,,,,<br>
Application,NovaSeq FASTQ Only,,,,,,,,,,,,,,,,,,,,,,<br>
Chemistry,Amplicon,,,,,,,,,,,,,,,,,,,,,,<br>
,,,,,,,,,,,,,,,,,,,,,,,<br>
**[Reads],,,,,,,,,,,,,,,,,,,,,,,**<br>
**101**,,,,,,,,,,,,,,,,,,,,,,,<br>
**101**,,,,,,,,,,,,,,,,,,,,,,,<br>
,,,,,,,,,,,,,,,,,,,,,,,<br>
**[Settings],,,,,,,,,,,,,,,,,,,,,,,**<br>
**Adapter**,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,,,,,,,,,,,,,,,<br>
**AdapterRead2**,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,,,,,,,,,,,,,,,<br>
**Read1StartFromCycle**,1,,,,,,,,,,,,,,,,,,,,,,<br>
**Read2StartFromCycle**,4,,,,,,,,,,,,,,,,,,,,,,<br>
,,,,,,,,,,,,,,,,,,,,,,,<br>
**[Data],,,,,,,,,,,,,,,,,,,,,,,**<br>
**Sample_ID**,Sample_Name,**Sample_Project**,Sample_Pool,Species,RIN,concentration,sample_order,pipeline_pofile,Assay,Index_Adapters,Strandness,Read1StartFromCycle,Read2StartFromCycle,Adapter,AdapterRead2,Paired,Read_1_cycles,Read_2_cycles,**Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2**<br>
LUIN_03,LUIN_03,2021_145,CTGpool_0174,Other,NA,,1,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,F11,U086,GAACCGCG,U086,TAAGGTCA<br>
LUIN_06,LUIN_06,2021_145,CTGpool_0174,Other,NA,,2,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,G11,U087,CTCACCAA,U087,TTGCCTAG<br>
LUIN_14,LUIN_14,2021_145,CTGpool_0174,Other,NA,,3,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,H11,U088,TCTGTTGG,U088,CCATTCGA<br>
LUIN_26,LUIN_26,2021_145,CTGpool_0174,Other,NA,,4,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,A12,U089,TATCGCAC,U089,ACACTAAG<br>
LUIN_29,LUIN_29,2021_145,CTGpool_0174,Other,NA,,5,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,B12,U090,CGCTATGT,U090,GTGTCGGA<br>
LUIN_31,LUIN_31,2021_155,CTGpool_0174,Other,NA,,6,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,C12,U091,GTATGTTC,U091,TTCCTGTT<br>
MAAK_03,MAAK_03,2021_155,CTGpool_0174,Other,NA,,7,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,D12,U092,ACGCACCT,U092,CCTTCACC<br>
MAIV_22,MAIV_22,2021_155,CTGpool_0174,Other,NA,,8,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,E12,U093,TACTCATA,U093,GCCACAGG<br>
MAIV_24,MAIV_24,2021_165,CTGpool_0174,Other,NA,,9,rnaseq_total,SMARTer Stranded Total RNA-Seq Kit v2 Pico Input Mammalian,SMARTer RNA Unique Dual Index Kit -96U Set A,reverse,1,4,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,false,101,,F12,U094,CGTCTGCG,U094,ATTGTGAA<br>

----




## Dependencies
Running the ctg-rnaseq pipeline on aurora lsens4 clusters requires:

### N E X T F L O W
[Nextflow](https://www.nextflow.io/) version 19.04.1 build 5072 was used as of ctg-rnaseq 2.1.6.

### Slurm
[The Slurm Workflow Manager](https://slurm.schedmd.com/tutorials.html) version 17.02.1-2 was used as of ctg-rnaseq 2.1.6.

### Singularity
Containers are built using [Singularity](https://sylabs.io/singularity). As of of ctg-rnaseq 2.1.6, singularity version 3.7.0-1.el7 was installed on LS4.



## ctg-rnaseq & ctg-singualrity on Git
Main git repository [https://github.com/cirrina/ctg-rnaseq](https://github.com/cirrina/ctg-rnaseq/tags)


### Git version tags 
Versions are stored as tags in git repository [https://github.com/cirrina/ctg-rnaseq/tags](https://github.com/cirrina/ctg-rnaseq/tags)
Versisons are specified as X.Y.ZZ, e.g. version 2.1.6 and tagged in the git site.

### Download tagged versions from Git

Script for download & extraction og a specified repo (not only ctg-rnaseq):

```
#!/bin/bash
read -p 'Git cirrina repo, e.g. singularity-ctg-rnaseq: ' gitRepo
read -p 'Git cirrina repo tag, e.g. 1.2.0: ' gitTag
cd /Users/david/scripts/git-version-tags
wget https://github.com/cirrina/${gitRepo}/archive/refs/tags/${gitTag}.tar.gz
tar -zxvf ${gitTag}.tar.gz
rm -rf ${gitTag}.tar.gz
```


## LSENS LS4 structure
### ctg-rnasesq versions
The ctg-rnaseq pipeline is deploed onto ls4 in folders:

Development versions: `/projects/fs1/shared/ctg-dev/pipelines/ctg-rnaseq`<br> 
Production ready versions: `/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/`

Versions directories are sprcified only by their respective version.<br>

```
/projects/fs1/shared/ctg-dev/pipelines/ctg-rnaseq/
├── 2.0.0
├── ...
├── 2.1.5
└── 2.1.6
```
### singularity containers


### Genome References 



**Note:** Do **not add** the ctg-rnaseq script directories to **PATH**. Instead run the `rnaseq-primer` script usging full path - thus allowing proper version control, e.g.



# Genome References

The ctg-rnaseq uses multiple references. Some are used as-is whereas others are modified. 
Notes on reference libraries on lsens

# Singularity containers
    - https://github.com/cirrina/singularity-bladderreport
    - https://github.com/cirrina/singularity-ctg-rnaseq
    - https://github.com/cirrina/singularity-uroscan



# Bin


## Uroscan 
The uroscan pipeline uses references prepared by CMD anno-2020.
uroscan. These are non-controlled references used as-is from the Lennart server. All are based on GRCh37. One reference, the RefFlat was compiled as belw.  



# Pipeline Modules



## demultiplexing: bcl2fastq
[bcl2fastq](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf): Converts raw basecalls to fastq, and demultiplex samples based on indexes.  

### SampleSheet
Sample sheet must be given in Illumina IEM csv format. See SampleSheet section below. 

### FASTQ file names
FASTQ files are named by blc2fastq with the following logic: **{samplename}\_{S1}\_{L001}\_{R1}\_001.fastq.gz**

- **{samplename}**: Defined by `Sample_Name` column in [Data] section of SampleSheet. Note that this pipeline will **force** `Sample_Name` to the same as `Sample_ID`. If Sample_Name column is specified but do not match the Sample_ID, the FASTQ files reside in a `Sample_ID` subdirectory where files use the Sample_Name value.
- **{S1}**: The number of the sample based on the order that samples are listed in the sample sheet, starting
with 1. In ctg-rnaseq pipeline the FASTQ names are expected to be supplied in the SampleSheet, under fastq1, fastq2 columns. **NOTE**. Reads that cannot be assigned to any sample are written to a FASTQ file as sample number 0 and
excluded from downstream analysis.
- **{L001}**: The lane number of the flow cell, starting with lane 1, to the number of lanes supported. When using the `--no-lane-splitting` flag, as default in ctg pipelines, lanes will not be split into different files and files always named **L001**. When the Lane column in the Data section is not used, all lanes are converted. Otherwise, only populated lanes are converted.
- **{R1}**: The read. In the example, R1 indicates Read 1. R2 indicates Read 2 of a paired-end run.
- **{_001.fastq.gz}**: The last portion of the file name is always 001.fastq.gz.

When a sample sheet contains multiplexed samples, the software: Places reads without a matching index adapter sequence in the Undetermined_S0 FASTQ file(s).


### bcl2fastsq: ctg-rnaseq 2.1.6 (version & code)
container: '/projects/fs1/shared/ctg-containers/rnaseq/singularity-ctg-rnaseq-1.0.2.sif'
version: bcl2fastq v2.19.0.316

```
mkdir -p ${bcl2fastq_dir}
bcl2fastq -R ${runfolderdir} \\
          --sample-sheet ${samplesheet_demux} \\
          --no-lane-splitting  \\
          -r 1 \\
          -p ${task.cpus}  \\
          -w 1  \\
          --output-dir ${bcl2fastq_dir}

find ${bcl2fastq_dir} -user $USER -exec chmod g+rw {} +
```

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



## RSEM & bowtie2
- [RSEM publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)
- [RSEM at GitHub](https://github.com/deweylab/RSEM)
- [rsem-calculate-expression command description](https://deweylab.github.io/RSEM/rsem-calculate-expression.html)

RSEM is currently (v 2.1.x) only used available for the `uroscan` profile. RSEM can use multiple aligners, e.g. bowtie2, star etc. Uroscan uses bowtie2 as aligner and a hg19 reference genome. 



### RSEM UroScanSeq: ctg-rnaseq 2.1.6 (version, code & references)
- container: 'singularity-uroscan-1.0.1.sif' <br>
- version: RSEM v1.3.0 <br>
- version: bowtie2 version 2.3.3.1
- gtf_hs = `/projects/fs1/shared/references/uroscan/rsem_bowtie2/GRCh37/Homo_sapiens.GRCh37.75.gtf` <br>
- rsem\_bowtie2\_genome_hs  =  `/projects/fs1/shared/references/uroscan/rsem_bowtie2/GRCh37/GRCh37` <br>
- rsem_star_genome_hs =  `/projects/fs1/shared/references/uroscan/star` <br>

Nextflow chunk:

```
mkdir -p ${rsemdir}
rsem-calculate-expression \\
    --num-threads ${task.cpus} \\
    --paired-end \\
    --bowtie2 \\
    --bowtie2-path /opt/software/uroscan_env/bin \\
    --estimate-rspd \\
    --append-names \\
    --no-bam-output \\
    ${rsemfiles} \\
    ${genome} \\
    ${rsemdir}/${sid}.rsem
```

Parameters:<br>

- `--bowtie2` Uses bowtie2 aligner reference genome GRCh37
- `--no-bam-output`: No bam output is produced using the aligneer (bowtie2). Uroscan bam files are produced using STAR. 
- `--paired-end` : **hardcoded at present** - the uroscan pipeline will always be paried end sequencing. 
- `--estimate-rspd` : Set this option if you want to estimate the read start position distribution (RSPD) from data.


Parameters not used:<br>
- `--stranddcness`: For Illumina TruSeq Stranded protocols, please use 'reverse'. (Default: 'none'). **Note** Urocsan uses a stranded protocol (TruSeq). Thus this `--strandness 'reverse'` should be used. This parameter was not used in the original uroscan pipeline set up at CMD and therefore ommmited here.


References were prepared for uoscan pipeline with shell script (from CMD):

```
/data/bnf/sw/RSEM-1.3.0/rsem-prepare-reference \
--gtf Homo_sapiens.GRCh37.75.gtf \
--bowtie2 \
--bowtie2-path /data/bnf/sw/bowtie2/2.3.3/ \
/data/bnf/ref/b37/human_g1k_v37.fasta  \
/data/bnf/ref/rsem/GRCh37
```

### Possible refinements ...
- Prepare bowtie2 rsem reference for hg38 on lsens
- change --paired-end to alternative for non paired end runs. {paired} variable.


<br>
<br>
## STAR
Uroscan and rnaseq is run using different versions of star

### File Names

By default STAR bam files are named


### rnaseq & rnaseq_total (2.1.6)
- container: '/projects/fs1/shared/ctg-containers/rnaseq/singularity-ctg-rnaseq-1.0.2.sif'<br>
- version: STAR version 2.7.6a <br>
- star_genome_hs      =  `/projects/fs1/shared/references/hg38/star/star_2.7.6a/` <br>
- star_genome_rn      =  `/projects/fs1/shared/references/rattus_norvegicus/Rnor_6.0/star/star_2.7.6a` <br>
- star_genome_mm      =  `/projects/fs1/shared/references/mm10/star/star-2.7.6a/` <br>

Genomes are built from: 


```
STAR --genomeDir ${genome} \\
  --readFilesIn ${starfiles} \\
  --runThreadN ${task.cpus}  \\
  --readFilesCommand zcat \\
  --outSAMtype BAM SortedByCoordinate \\
  --limitBAMsortRAM 10000000000 \\
  --outFileNamePrefix ${stardir}/${sid}_
```


### uroscan (2.1.6)
* container: `singularity-uroscan-1.0.1.sif`
* version: STAR_2.5.3a
* star_genome_hs = `/projects/fs1/shared/references/uroscan/star`

```
STAR --genomeDir ${genome} \\
  --readFilesIn ${starfiles} \\
  --runThreadN ${task.cpus}  \\
  --readFilesCommand zcat \\
  --outSAMtype BAM SortedByCoordinate \\
  --limitBAMsortRAM 10000000000 \\
  --outFileNamePrefix ${stardir}/${sid}_
```

### Notes on STAR

**Total RNA & multimapping reads**:<br>
Whemn aligning total RNA seq libraries, multimapping reads are very common resulting in large bam-files (sometimes affecting downstream processes).

Filtering of multimapping reads can be performed with STAR flag `--outFilterMultimapNmax` to output only the best scoring multipmapping read. 

For multi-mappers, all alignments except one are marked with 0x100 (secondary alignment) in the FLAG (column 2 of the SAM). The unmarked alignment is selected from the best ones (i.e. highest scoring). This default behavior can be changed with --outSAMprimaryFlag AllBestScore option, that will output all alignments with the best score as primary alignments (i.e. 0x100 bit in the FLAG unset).

[Picard: bam flags explained](https://broadinstitute.github.io/picard/explain-flags.html)

See [thread](https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools) and [thread](https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools) 

```
singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-samtools-1.9.sif samtools  view -b -F 0x100 15PL19843_01_01_Aligned.sortedByCoord.out.bam > noMultiMap.bam
```
**Use SamTools to filter multimapping reads**<br>
How many primary aligned reads (0x100 = 0) are in the bwa_local.sort.dup.bam file?
samtools view -F 0x104 -c bwa_local.sort.dup.bam





## salmon
[salmon ducumentation](https://salmon.readthedocs.io/en/latest/salmon.html)


### uroscan (2.1.6)
* container: `singularity-uroscan-1.0.1.sif`
* version: 
* salmon_transcripts_hs =  `/projects/fs1/shared/references/uroscan/salmon/GRCh37.transcripts`


```
salmon quant -l A \\
      -i  ${transcripts} \\
      -1  ${fastqdir}/${read1} \\
      -2  ${fastqdir}/${read2} \\
      -p  6 --validateMappings \\
      -o  ${salmondir}/${sid}_0.salmon.salmon \\
      --no-version-check
```
To allow Salmon to automatically infer the library type, simply provide -l A


### Notes on salmon
#### Treads
but there is a point beyond which allocating more threads will not speed up alignment-based quantification. We find that allocating 8 — 12 threads results in the maximum speed, threads allocated above this limit will likely spend most of their time idle / sleeping.



### How to build salmon references

he first is to compute a set of decoy sequences by mapping the annotated transcripts you wish to index against a hard-masked version of the organism’s genome. This can be done with e.g. MashMap2, and we provide some simple scripts to greatly simplify this whole process. Specifically, you can use the generateDecoyTranscriptome.sh script, whose instructions you can find in this README.

[Pre built Salmon references](http://refgenomes.databio.org/)

`> ./bin/salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31`




## featureCounts
featureCounts is a part of the [subread package](http://subread.sourceforge.net/). https://usermanual.wiki/Pdf/SubreadUsersGuide.127634897/html. 
[Rsubread documentation](https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts) for featureCounts

By default multi-mapping reads are not counted included in the mapping *"Due to the mapping ambiguity, it is recommended that multi-mapping
reads should be excluded from read counting (default behavior of featureCounts program) to produce as accurate counts as possible"* (can be turned on e.g. by the `-M` option).

For large bam files (e.g. total RNA libs with high numbers of multi-mapping reads) featureCounts would stall on ls4, taking days. Even though multimappers is not to be counted by featureCounts. Therefore, post 2.1.0 and additional filtering step of bam files is performed before featureCounts. For multi-mappers, all alignments except one are marked with 0x100 (secondary alignment) in the FLAG (column 2 of the SAM). The unmarked alignment is selected from the best ones (i.e. highest scoring). This default behavior can be changed with --outSAMprimaryFlag AllBestScore option, that will output all alignments with the best score as primary alignments (i.e. 0x100 bit in the FLAG unset).

[Picard: bam flags explained](https://broadinstitute.github.io/picard/explain-flags.html)

See [thread](https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools) and [thread](https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools) 

```
singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-samtools-1.9.sif samtools  view -b -F 0x100 15PL19843_01_01_Aligned.sortedByCoord.out.bam > noMultiMap.bam
```
**Use SamTools to filter multimapping reads**<br>
How many primary aligned reads (0x100 = 0) are in the bwa_local.sort.dup.bam file?
samtools view -F 0x104 -c bwa_local.sort.dup.bam


countMultiMappingReads = FALSE (default)
logical indicating if multi-mapping reads/fragments should be counted, FALSE by default. If TRUE, a multi-mapping read will be counted up to N times if it has N reported mapping locations. This function uses the NH tag to find multi-mapping reads.


### uroscan, rnaseq & rnaseq_total (2.1.6)


```
mkdir -p ${stardir_filtered}
cd ${stardir_filtered}
samtools  view -b -F 0x104  ${stardir}/${bam} >  ${stardir_filtered}/${bam}
```
 
```
  mkdir -p ${featurecountsdir}
    # cd ${stardir}
    cd ${stardir_filtered}
    bamstring=\$(echo $bams | sed 's/,/ /g' | sed 's/\\[//g' | sed 's/\\]//g' )
    echo \${bamstring}
    echo "gtf: ${gtf}"
    featureCounts -T ${task.cpus} \\
      -t ${params.fcounts_feature} \\
      --extraAttributes gene_name,gene_type \\
      -a ${gtf} -g gene_id  \\
      -o ${featurecountsdir}/${projectid}_geneid.featureCounts.txt \\
      -p \\
      -s ${strand_numeric} \${bamstring}

    #find ${featurecountsdir} -user $USER -exec chmod g+rw {} +
```




### uroscan (2.1.6)
* container: `singularity-uroscan-1.0.1.sif`
* version: 
* salmon_transcripts_hs =  `/projects/fs1/shared/references/uroscan/salmon/GRCh37.transcripts`
 
  


## RSEQC

```
mkdir -p ${rsemdir}
rsem-calculate-expression \\
    --num-threads ${task.cpus} \\
    --paired-end \\
    --bowtie2 \\
    --bowtie2-path /opt/software/uroscan_env/bin \\
    --estimate-rspd \\
    --append-names \\
    --no-bam-output \\
    ${rsemfiles} \\
    ${genome} \\
    ${rsemdir}/${sid}.rsem
```

rnaseq:  ``


[http://rseqc.sourceforge.net/]()

[https://sourceforge.net/projects/rseqc/]()
[https://sourceforge.net/projects/rseqc/files/BED/]()

inner_distance.py
RSEQC has been problematic to run on lsens. Possibly due to when run frpm conda. As of later versions, rseqc is run in an individual cpontainer and no longer installed through conda.




## qualimap
[http://seqanswers.com/forums/showthread.php?t=66538]()



## fastqscreen



## Known issues

### SAMtools and Conda
[https://github.com/merenlab/anvio/issues/1479]()


### RSEQC & Qualimap
CAN NOT GET RSQEQC to conistently work. Same as qualimap below.
get error messages some times
may have to do with samtools in conda
DO NOT RUN until container with rseqc/samtools/rsamtools is OK


In prevoius versions CAN NOT GET RSQEQC to conistently work. Same as qualimap below.
get error messages some times

For unknown reasons. Possible larger bam files. qualimap, often but not always, crashed. The nodes keep filling up with temporary bam files. May have to do with samtools library problem when initated from a conda environment.

Error messages such as:
WARN: Killing pending tasks (6)
ERROR ~ Error executing process > 'qualimap (20PL05864_01_01)'

Caused by:
  Process `qualimap (20PL05864_01_01)` terminated for an unknown reason -- Likely it has been terminated by the external system

Command executed:

  mkdir -p /projects/fs1/shared/ctg-projects/rnaseq/2021_151/nf-output/qc/qualimap

  ## export JAVA_OPTS="-Djava.io.tmpdir=/data/tmp"
  ## /data/bnf/sw/qualimap_v2.2.1/qualimap --java-mem-size=12G rnaseq -bam /data/bnf/bam/rnaseq/21KF00020.STAR.sort.bam -gtf /data/bnf/ref/rsem/GRCh37/Homo_sapiens.GRCh37.75.gtf -pe -outdir /data/bnf/postmap/rnaseq/21KF00020.STAR.qualimap.folder
  # qualimap --java-mem-size=12G rnaseq -bam /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/star/21KF00082_Aligned.sortedByCoord.out.bam -gtf /projects/fs1/shared/uroscan/references/rsem/GRCh37/Homo_sapiens.GRCh37.75.gtf -pe -outdir /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/qualimap/21KF00082.STAR.qualimap.folder


## BladderReport


## Version 2.2 outline

### ctg-demux script
Demux should now be berformed outside (prior) ctg-rnaseq main script. This by other script `ctg-demux`.
This new scritp shoud

1. Check SampleSheet for required parameters and format - python script. This to replace the R script used in previuos versions. Trigger script using daemon upon completed transfer to ls4?
	* Make option (default) to output .fastq and .bam file names. File namings will be dependent on input order (samplesheet row number), Lane (if used), etc. A second option to set row number after splitting on unique projects. 
2. Output:
	* SampleSheet-original
	* SampleSheet-demux
	* SampleSheet-project - used for nextflow. Split samplesheet based on Sample_Project 
3. Run bcl2fastq demux. Output files to `shared/bcl2fastq-output-fastq`
4. Output file to be recognized by daemon & start project specific pipeline (ctg-rnaseq) 

### ctg-rnaseq changes in structure
Main change in strategy to adapt to Pers pipelines.
* change delivery root to `shared/ctg-delivery`
* The project work folder will be the main source of script backup (as compared to the current ctg-qc folder)
* Output files to be delivered to customer should be put immediately into the delivery foler.
* 

1. fastq files. Default is now to obtain fastq-files from the demux output path. `shared/bcl2fastq-output-fastq`. Check at this location. If there, then prepare to move fastq files to `delivery/projid/fastq`. If no there, use the `-f` flag to define

2. change driver. A new driver strategy should be to move the  `rnaseq-driver` to `\ctg-tools`. This driver shoud be a **light**. The main function is to be a wrapper that identifies the correct script version, create a project dir (if needed) and executes the pipeline. **Q!** Should project scripts be overwritten or NOT?.
	* input SampleSheet. Should contain ALL relevant params (exept for -p prime and -r resume, and -f fastq)
	* identify the version to be run
	* Create the project dir and copy all relevant files to there
	* execute `rnaseq-main.sh` scritp that in much crresponds to the old `rnaseq-driver`
