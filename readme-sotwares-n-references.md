
# ctg-rnaseq pipeline
---------------------


Primary data processing pipeline for Illumina RNA-seq data produced at CTG. Built to run using Nextflow and (multiple) Singularity containers. The pipeline is designed to start from an Illumina Runfolder to demultiplex data, different qc measures, alignment and transcript summarization.

The pipeline is designed to handle multiple different RNAseq Assays (library preparation methods), and reference Species. Different assays will require differences in read strandness, read trimming etc. These input parameters are either set through the Sample Sheet or profile-specific and specified in the nextflow configurations file (see below).

**Note on Demultiplexing and Fastq:** As of versions 3.x and above, demultiplexing is no longer performed within the pipeline. Instead the pipeline assumes fastq files as starting input. Demultiplexing is performed on a RunFolder basis using the ctg-demux pipeline and may include multiple projects rather than one single project (as of versions ≤2.x. 

**Note on RunFolder(s):** Versions ≤2.x. are designed to process samples in **one single sequencing run** (all run within one Illumina Runfolder). If a project uses multiple sequencing runs, manual adjustments have to be made to fit the pipeline execution requirements (see below).  

**Note on ProjectId/Sample_Project** Only **one project** is allowed per pipeline run. `ProjectId` is supplied through the SampleSheet [Header]. Thus, for a run, all individual `Sample_Project` entries in the SampleSheet [Data] section must be same as `ProjectId` for all **samples**. 

The `project_id` (supplied by `ProjectId` in SampleSheet) will owerwrite the `Sample_Project` column in sample sheet - again  (As of v1.3 the `ProjedId` is supplied through the SampleSheet, **not** using -i flag as for versions 2.x. Thus, different projects (& library pools) within a pooled run must be processed separately.  

**Note on multiple Species:** The pipeline is designed for **one and the same species**. This is due to that transcript summarization i must be performed using the same refence GFF file. If multiple species are used, then run multiple pipeline runs, defined by unique project IDs.




**nextflow profiles {PipelineProfile}**
The pipeline is designed for three different profiles/main configurations.<br>
- rnaseq. Sequencing of a mRNA libraries. 
- rnaseq_total: 
- uroscan 
 


## Running the ctg-rnaseq pipeline
1. Retrieve (or build) the **Singularity container(s)**. Containers are specified in:
2. Add the correct .sif paths to:
 	- `nextflow.config`
	- `rnaseq-driver`
3. make sure that the `ctg-rnaseq` scriptsdir matches the `scripts_root` in the `rnaseq-driver` script.
4. Make sure that the software versions as installed in .sif are compatible with references defined in nextflow.config, e.g. STAR indexed references.
5. Edit your samplesheet to fullfill all requirements. See section `SampleSheet` below. The SampleSheet (-s) ,must be supplied and present within the pipeline execution dir.  
6. Run the `naseq-driver` from within `Illumina Sequencing Runfolder`
7. *OR* run the `naseq-driver`from within the `Project work folder`. This requires that path to fastq directory is specified. Typically this is performed when resuming a failed run. 
8. Optional: A project runfolder can be primed without starting the pipeline. This in order to modify parameters in the `nextflow.config.project.XXX` or the `nextflow.config` files. 


Example initiate from within Illumina runfolder (v2.x):

```
/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/2.0.0/rnaseq-driver \
  -s 2021_024_SampleSheet-IEM.csv
```

Re-run a failed run, debug etc within a project folder:

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



## Dependencies
Running the ctg-rnaseq pipeline on aurora lsens4 clusters requires:

### N E X T F L O W
[Nextflow](https://www.nextflow.io/) version 19.04.1 build 5072 was used as of ctg-rnaseq 2.1.6.

### Slurm
[The Slurm Workflow Manager](https://slurm.schedmd.com/tutorials.html) version 17.02.1-2 was used as of ctg-rnaseq 2.1.6.

### Singularity
Containers are built using [Singularity](https://sylabs.io/singularity). As of of ctg-rnaseq 2.1.6, singularity version 3.7.0-1.el7 was installed on LSENS.



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


## LSENS structure
### ctg-rnasesq pipeline
The ctg-rnaseq pipeline is deploed onto ls4 in folders:

Development versions: `/projects/fs1/shared/ctg-dev/pipelines/ctg-rnaseq`<br> 
Production ready versions: `/projects/fs1/shared/ctg-pipelines/ctg-rnaseq/`

Versions are sprcified only by their respective version.<br>

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



# PIPELINE 

## Demux FASTQ

- bcl2fastq
- checkfiles fastq
- fastqc

rsem
salmon

star
  + index bam
  + markdups


## bcl2fastq
!! note : Fix option so that Undetermined fastq is NOT outputed if multipke projecs are run!!




## rsem
https://deweylab.github.io/RSEM/rsem-calculate-expression.html

rsem can use multiple aligners, e.g. bowtie2, star etc. Uroscan uses bowtie2 as aligner.

Preparing references for uoscan / shell script from CMD

```
/data/bnf/sw/RSEM-1.3.0/rsem-prepare-reference --gtf Homo_sapiens.GRCh37.75.gtf --bowtie2 --bowtie2-path /data/bnf/sw/bowtie2/2.3.3/ /data/bnf/ref/b37/human_g1k_v37.fasta  /data/bnf/ref/rsem/GRCh37
```

Prepare bowtie2 rsem reference for hg38 on lsens





## UroScanSeq
`--bowtie2` Uses bowtie2 aligner with reference genome GRCh37: `/projects/fs1/shared/uroscan/references/rsem/GRCh37/GRCh37`.

Flags:
`--no-bam-output`: No bam output is produced using the aligneer (bowtie2)
`--paired-end` : hardcoded as is - SHOULD BE CHANGED ACCORDING TO PROTOCOL
`--stranddcness`: For Illumina TruSeq Stranded protocols, please use 'reverse'. (Default: 'none'). NOT USED AS IS
`--estimate-rspd` : Set this option if you want to estimate the read start position distribution (RSPD) from data.



## last working nextflow script. hanging




## STAR
Uroscan and rnaseq is run using different versions of star

#### Versions:

rnaseq pipeline: 
uroscan pipeline:  


```
STAR --genomeDir ${genome} \\
  --readFilesIn ${starfiles} \\
  --runThreadN ${task.cpus}  \\
  --readFilesCommand zcat \\
  --outSAMtype BAM SortedByCoordinate \\
  --limitBAMsortRAM 10000000000 \\
  --outFileNamePrefix ${stardir}/${sid}_
```
NOTES: consiter adding flag to filter multimapping reads. Apparent when aligning total RNA where multimapping are very common.

outFilterMultimapNmax

For multi-mappers, all alignments except one are marked with 0x100 (secondary alignment) in
the FLAG (column 2 of the SAM). The unmarked alignment is selected from the best ones (i.e.
highest scoring). This default behavior can be changed with --outSAMprimaryFlag AllBestScore
option, that will output all alignments with the best score as primary alignments (i.e. 0x100 bit in the FLAG unset).

```
singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-samtools-1.9.sif samtools  view -b -F 0x100 15PL19843_01_01_Aligned.sortedByCoord.out.bam > noMultiMap.bam
```

[https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools]()

How many primary aligned reads (0x100 = 0) are in the bwa_local.sort.dup.bam file?
samtools view -F 0x104 -c bwa_local.sort.dup.bam


[https://broadinstitute.github.io/picard/explain-flags.html]()

## bowtie2 (uroscan only)



## salmon (uroscan only)
[https://salmon.readthedocs.io/en/latest/salmon.html]()

`> ./bin/salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31`


but there is a point beyond which allocating more threads will not speed up alignment-based quantification. We find that allocating 8 — 12 threads results in the maximum speed, threads allocated above this limit will likely spend most of their time idle / sleeping.


he first is to compute a set of decoy sequences by mapping the annotated transcripts you wish to index against a hard-masked version of the organism’s genome. This can be done with e.g. MashMap2, and we provide some simple scripts to greatly simplify this whole process. Specifically, you can use the generateDecoyTranscriptome.sh script, whose instructions you can find in this README.

[http://refgenomes.databio.org/]()



## featurecounts
FeatureCoiunts is run

countMultiMappingReads = FALSE (default)
logical indicating if multi-mapping reads/fragments should be counted, FALSE by default. If TRUE, a multi-mapping read will be counted up to N times if it has N reported mapping locations. This function uses the NH tag to find multi-mapping reads.


## rseqc

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
