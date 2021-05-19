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

ctg_qc_dir      = qc_root+'/rnaseq'
deliverydir = delivery_root + '/' + projectid



// set up the ctg qc folder
//
process setup_ctg_qc {

  file().mkdir()

  script:

  """

  """


}



process setup_delivery {



}

process md5sum_delivery {

}



process multiqc_delivery {

}




process cleanup_projectdir {


}
