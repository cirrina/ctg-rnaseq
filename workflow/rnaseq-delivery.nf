/* ===============================================================
  *      -- RNASEQ DELIVERY --
  =============================================================== */


//  project specific config
// -----------------------------

// root directories
project_root        =  params.project_root
delivery_root       =  params.delivery_root
qc_root             =  params.qc_root
//log_root            =  params.log_root

//  project  and run folders
projectid           =  params.projectid
projectdir          =  params.projectdir

//  samplesheets
samplesheet         =  params.samplesheet
samplesheet_demux   =  params.samplesheet_demux

//  demux specific
pooled              =  params.pooled
runfolderdir        =  params.runfolderdir
fastqdir            =  params.fastqdir
deliver_raw         =  params.deliver_raw
deliver_fastq       =  params.deliver_fastq


//  Other Fixed Directories (must be identical to dirnames in rnaseq.nf file)
// -----------------------------
outputdir =  projectdir+'/nf-output' // main ooutput directory for files genetated with the Pipeline

featurecountsdir = outputdir+'/featurecounts'
stardir = outputdir+'/star'
markdupsdir = outputdir+'/markdups_bam_tmp'
fastqcdir = outputdir+'/fastqc'
markdupsqcdir = outputdir+'/markdups'
rnaseqmetricsdir = outputdir+'/rnaseqmetrics'
multiqcctgdir = outputdir+'/multiqc_ctg'
fastqscreendir = outputdir+'/fastqscreen'

ctg_qc_dir      = qc_root+projectid
deliverydir     = delivery_root + '/' + projectid



// set up the ctg qc folder
// multiqc and fastqc files
process setup_ctg_qc {
  cpus 4
  tag "$id"
  memory '32 GB'
  time '3h'

  script:
  file(ctg_qc_dir).mkdir()
  file(ctg_qc_dir+'multiqc').mkdir()
  file(ctg_qc_dir+'fastqc').mkdir()

  """
  cp -i ${multiqcctgdir}/* ${ctg_qc_dir}/multiqc
  cp -i ${fastqcdir}/* ${ctg_qc_dir}/fastqc
  """
}



process setup_delivery {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'

  script:
  file(deliverydir).mkdir()

  """
    cp -i samplesheet ${deliverydir}/

  """


}

process rsync_fastq {

  when:
  deliver_fastq

  script:
  """
  """
}

process md5sum_delivery {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'
}



process multiqc_delivery {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'
}




process cleanup_projectdir {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'

}
