

## CTG save


Scripts
configs
sample sheets &
logs (log.rscript, log.nextflow x2. .nextflow)

## ctg-qc

multiqc ctg
fastqc +++

other qc
(interop)


## delivery
fastqc +++




#markdups errorReport

singularity exec --bind /projects/fs1/  /projects/fs1/shared/ctg-containers/rnaseq/singularity-ctg-rnaseq-1.0.2.sif picard MarkDuplicates \
     INPUT=/projects/fs1/shared/ctg-projects/rnaseq/2021_095/nf-output/delivery/star/7162_13_Aligned.sortedByCoord.out.bam \
     OUTPUT=/projects/fs1/shared/ctg-projects/rnaseq/2021_095/nf-output/markdups_bam_tmp/7162_13_Aligned.sortedByCoord.out.bam \
     METRICS_FILE=/projects/fs1/shared/ctg-projects/rnaseq/2021_095/nf-output/markdups/7162_13_bam.MarkDuplicates.metrics.txt \
     TAGGING_POLICY=All \
     REMOVE_DUPLICATES=false \
     ASSUME_SORTED=true \
     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000 \
     QUIET=true \
     VERBOSITY=WARNING





multiqc ctg

THEN

stage_delivery
copy - scripts - logs
-------
setup_deliverytemp -
move_fastq


deliverysamplesheets
deliveryscripts
deliveryconfigs
deliverylogs


COPY THESE TO CTG SAVE


THEN
setup_ctg_save



THEN FINALLY
multiqc_delivery
md5sum delivery


Final
finalize_delivery
copy to nas sync
chmods etc
