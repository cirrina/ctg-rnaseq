# bladderreport

## test on ls4  server using singularity container
exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/rnaseq/singularity-bladderreport-1.0.1.sif \

Rscript -e "library(rmarkdown); rmarkdown::render('/projects/fs1/shared/ctg-projects/uroscan/2021_024/bin/bladderreport/bladderreport-ctg-1.1.0..Rmd', params=list( sampleid = '21KF00082',  rsem_in = '/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/rsem/21KF00082.rsem.genes.results', star_qc = '/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/star/21KF00082_Log.final.out', clarity_id = '21KF00082', RIN = '9.5', concentration = '130.4', date= '2021-09-14'), output_file='/projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/bladderreport/21KF00082.bladderreport.html' )"
