# ctg-rnaesq

## Requirements

1. Clone and build the Singularity container for this pipeline: (currently. https://github.com/cirrina/ctg-singularity/tree/main/rnaseq/v1.2).
  + Add the path to the .sif in the nextflow.config `container = ` parameter under process{}
2. Edit your samplesheet to fullill the requirements. See section `SampleSheet` below
3. Edit the nextflow.config so that file paths are correct. Make sure that the STAR version installed in .sif matches the verrsion of STAR used to built references in nextflow.config.

```
nohup nextflow run pipe-seqonly-qc.nf > log.pipe-seqonly-qc.txt &
```

## Run pipeline using Primer & Driver
.. 1-3 from above  
4. Start by running the 'rnaseq-driver'
  + 
5.
6. Start driver from runfolder  

```
seqonly-driver
```
- This will use default values (e.g. /path/to/runfolder/CTG_SampleSheet.csv)

## Usage Seqonly-Driver
```
Usage: seqonly-driver [ -i META_ID ] [ -s SAMPLESHEET ] [ -b BCL2FASTQ ARGUMENTS ] [ -r RESUME ] [ -d DEMUX-OFF ] [ -h HELP ]


Optional arguments:
META-ID       -i : Set 'meta-id' for runfolder (e.g. 210330-seqonly). Default: Takes date of runfolder (before first _) and adds '-seqonly' as suffix
SAMPLESHEET   -s : Set samplesheet used for run (Default: CTG_SampleSheet.csv)
BCL2FASTQ arg -b : String with bcl2fastq argument. e.g. '--minimum-trimmed-read-length 20 --mask-short-adapter-reads 20'
RESUME        -r : Set if to resume nf-pipeline
DEMUX-OFF     -d : Set flag to skip bcl2fastq (then fastq must be in FQDIR)
HELP          -h : print help message
```

## Input

- Samplesheet (see `SampleSheet` section below)

### Pipeline steps:

* `Demultiplexing` (bcl2fastq): Converts raw basecalls to fastq, and demultiplex samples based on index (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).
* `FastQC`: FastQC calculates quality metrics on raw sequencing reads (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). MultiQC summarizes FastQC reports into one document (https://multiqc.info/).
* `multiQC`: Compile fastQC and cellranger count metrics in multiqc report (https://multiqc.info/)
* `md5sum`: md5sum of all fastq files


### Output:
* ctg-PROJ_ID-output
    * `qc`: Quality control output.
        * fastqc output (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * multiqc output: Summarizing FastQC output and demultiplexing (https://multiqc.info/)
    * `fastq`: Contains raw fastq files from cellranger mkfastq.


### Samplesheet requirements:

Illumina IEM. The Bold/Italic field below must be correct!

[Header]
IEMFileVersion,5  
Investigator Name,X  
Experiment Name,X  
Date,YYYY-MM-DD  
Workflow,GenerateFASTQ  
Application,NovaSeq FASTQ Only  
Instrument Type,NovaSeq  
Assay,Nextera XT  
Index Adapters,"Nextera XT v2 Index Kit A"  
Chemistry,Amplicon  

[Reads]  
***26***  
***78***  

[Settings]  
Adapter,***CTGTCTCTTATACACATCT***  

[Data]  
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description  
***S1***,***S1***,,,***N702***,***CGTACTAG***,,,2021_024,  
***S2***,***S2***,,,***N706***,***TAGGCATG***,,,2021_024,  


### Container  
https://github.com/perllb/ctg-seqonly/tree/main/container  
