
Bug Error when starting run to resume with flags, sample sheet is specified but forgot to specify fastq FILES
i.e. runfolder mode ... but f is not supplied ...
wanrs about

... ... expecting fastq_custom dir through -f flag ...
RunFolder is not properly supplied in samplesheet.
Must be supplied as 'RunFolder' within the [Header] section of sample sheet
If set to 'NA', then fastq_path (-f) must be supplied.


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


ERROR
Sed cannot use a path that include /
sed "s/fastqdir .*/fastqdir            =  ${fastqdir}/g" $nf_conf > tmp.txt ; mv tmp.txt $nf_conf
sed "s|fastqdir .*|fastqdir|g"



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

deliveryqc

COPY THESE TO CTG SAVE


THEN
setup_ctg_save

ctg_save_samplesheets = ctg_save_dir+'/samplesheets'
ctg_save_scripts = ctg_save_dir+'/scripts'
ctg_save_configs = ctg_save_dir+'/configs'
ctg_save_logs =  ctg_save_dir+'/logs'



THEN FINALLY
multiqc_delivery
md5sum delivery


Final
finalize_delivery
copy to nas sync
chmods etc



RANDOM error
processing file: bladderreport-ctg-1.1.1.Rmd
output file: bladderreport-ctg-1.1.1.knit.md


Output created: /projects/fs1/shared/ctg-projects/rnaseq/2021_024/nf-output/delivery/bladderreport/21KF00103.STAR.bladderreport_anonymous.html
[0929/123841.331119:ERROR:bus.cc(393)] Failed to connect to the bus: Failed to connect to socket /var/run/dbus/system_bus_socket: No such file or directory
[0929/123841.360704:WARNING:headless_browser_main_parts.cc(83)] Cannot create Pref Service with no user data dir.
[0929/123841.368472:ERROR:gpu_init.cc(426)] Passthrough is not supported, GL is disabled
[0929/123922.183096:INFO:headless_shell.cc(620)] Written to file 21KF00103.STAR.bladderreport.pdf.
chmod: cannot access '/projects/fs1/shared/ctg-projects/rnaseq/2021_024/nf-output/delivery/bladderreport/tmp_21KF00081': No such file or directory
