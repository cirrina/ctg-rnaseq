# References
The ctg-rnaseq uses multiple references. Some are used as-is whereas others are morified. The uroscan pipelie uses references prepared by CMD pre-2020.

Notes on reference libraries on lsens

uroscan. These are non-controlled references used as-is from the Lennart server. All are based on GRCh37. One reference, the RefFlat was compiled as belw.   



## ORDER OF EVENTS!


## SAMtools
https://github.com/merenlab/anvio/issues/1479

### FASTQ-related

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





### uroscan
`--bowtie2` Uses bowtie2 aligner with reference genome GRCh37: `/projects/fs1/shared/uroscan/references/rsem/GRCh37/GRCh37`.

Flags:
`--no-bam-output`: No bam output is produced using the aligneer (bowtie2)
`--paired-end` : hardcoded as is - SHOULD BE CHANGED ACCORDING TO PROTOCOL
`--stranddcness`: For Illumina TruSeq Stranded protocols, please use 'reverse'. (Default: 'none'). NOT USED AS IS
`--estimate-rspd` : Set this option if you want to estimate the read start position distribution (RSPD) from data.

## last working nextflow script. hanging
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



## STAR
Uroscan and rnaseq is run using different versions of star


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

singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-samtools-1.9.sif samtools  view -b -F 0x100 15PL19843_01_01_Aligned.sortedByCoord.out.bam > noMultiMap.bam

https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+SAMTools

How many primary aligned reads (0x100 = 0) are in the bwa_local.sort.dup.bam file?
samtools view -F 0x104 -c bwa_local.sort.dup.bam


https://broadinstitute.github.io/picard/explain-flags.html

## bowtie2



## salmon
https://salmon.readthedocs.io/en/latest/salmon.html

`> ./bin/salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31`


but there is a point beyond which allocating more threads will not speed up alignment-based quantification. We find that allocating 8 — 12 threads results in the maximum speed, threads allocated above this limit will likely spend most of their time idle / sleeping.


he first is to compute a set of decoy sequences by mapping the annotated transcripts you wish to index against a hard-masked version of the organism’s genome. This can be done with e.g. MashMap2, and we provide some simple scripts to greatly simplify this whole process. Specifically, you can use the generateDecoyTranscriptome.sh script, whose instructions you can find in this README.

http://refgenomes.databio.org/



## featurecounts
countMultiMappingReads = FALSE (default)
logical indicating if multi-mapping reads/fragments should be counted, FALSE by default. If TRUE, a multi-mapping read will be counted up to N times if it has N reported mapping locations. This function uses the NH tag to find multi-mapping reads.


## rseqc
CAN NOT GET RSQEQC to conistently work. Same as qualimap below.
get error messages some times
may have to do with samtools in conda
DO NOT RUN until container with rseqc/samtools/rsamtools is OK

http://rseqc.sourceforge.net/

https://sourceforge.net/projects/rseqc/
https://sourceforge.net/projects/rseqc/files/BED/

inner_distance.py





## qualimap
http://seqanswers.com/forums/showthread.php?t=66538

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



## fastqscreen
