# bladderreport

## test on ls4  server using singularity container
exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-bladderreport-1.0.1.sif \

Rscript -e "library(rmarkdown); rmarkdown::render('/projects/fs1/shared/ctg-projects/uroscan/2021_024/bin/bladderreport/bladderreport-ctg-1.1.0..Rmd', params=list( sampleid = '21KF00082',  rsem_in = '/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/rsem/21KF00082.rsem.genes.results', star_qc = '/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/star/21KF00082_Log.final.out', clarity_id = '21KF00082', RIN = '9.5', concentration = '130.4', date= '2021-09-14'), output_file='/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/bladderreport/21KF00082.bladderreport.html' )"


singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-bladderreport-1.1.3.sif Rscript -e "library('rmarkdown'); \
  rmarkdown::render( \
    '/projects/fs1/shared/ctg-projects/uroscan/2021_024/bin/bladderreport/bladderreport-ctg-1.1.0.Rmd',  \
    params = list(   \
      sampleid='21KF00106', \
      rsem_in='/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/rsem/21KF00106.rsem.genes.results', \
      star_qc='/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/star/21KF00106_Log.final.out', \
      RIN='7.8', \
      concentration='35.4', \
      clarity_id='21KF00106'),  \
      output_file='/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/bladderreport/21KF00106.STAR.bladderreport_anonymous.html')"


singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-bladderreport-1.1.2.sif chromium-browser --headless --disable-gpu --print-to-pdf=/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/bladderreport/21KF00082.bladderreport.pdf /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/bladderreport/21KF00082.bladderreport.html

singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-bladderreport-1.1.2.sif /projects/fs1/shared/ctg-projects/uroscan/2021_024/bin/bladderreport/bladder_noreport2txt.pl /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/bladderreport/21KF00082.bladderreport.html > /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/bladderreport/21KF00082.bladderreport.txt

#    chmod 664 ${sample_id}.STAR.bladderreport.pdf ${sample_id}.STAR.bladderreport.txt


cho === `date +"%T, %Y-%m-%d"` '# Qualimap 21KF00020'
export JAVA_OPTS="-Djava.io.tmpdir=/data/tmp"
/data/bnf/sw/qualimap_v2.2.1/qualimap --java-mem-size=12G rnaseq -bam /data/bnf/bam/rnaseq/21KF00020.STAR.sort.bam -gtf /data/bnf/ref/rsem/GRCh37/Homo_sapiens.GRCh37.75.gtf -pe -outdir /data/bnf/postmap/rnaseq/21KF00020.STAR.qualimap.folder

singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-ctg-rnaseq-1.0.2.sif qualimap --java-mem-size=12G rnaseq -bam /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/star/21KF00082_Aligned.sortedByCoord.out.bam -gtf /projects/fs1/shared/uroscan/references/rsem/GRCh37/Homo_sapiens.GRCh37.75.gtf -pe -outdir /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/qualimap/21KF00082.STAR.qualimap.folder
