# References
The ctg-rnaseq uses multiple references. Some are used as-is whereas others are morified. The uroscan pipelie uses references prepared by CMD pre-2020.

Notes on reference libraries on lsens

uroscan. These are non-controlled references used as-is from the Lennart server. All are based on GRCh37. One reference, the RefFlat was compiled as belw.   


## rsem
https://deweylab.github.io/RSEM/rsem-calculate-expression.html

rsem can use multiple aligners, e.g. bowtie2, star etc. Uroscan uses bowtie2 as aligner.

### uroscan
`--bowtie2` Uses bowtie2 aligner with reference genome GRCh37: `/projects/fs1/shared/uroscan/references/rsem/GRCh37/GRCh37`.

Flags:
`--no-bam-output`: No bam output is produced using the aligneer (bowtie2)
`--paired-end` : hardcoded as is - SHOULD BE CHANGED ACCORDING TO PROTOCOL
`--strandness`: For Illumina TruSeq Stranded protocols, please use 'reverse'. (Default: 'none'). NOT USED AS IS
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


## bowtie2



## salmon


## featurecounts


## rseqc



## qualimap


## fastqscreen
