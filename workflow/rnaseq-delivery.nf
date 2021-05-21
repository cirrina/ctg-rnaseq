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

ctg_qc_dir      = qc_root + '/' + projectid
deliverydir     = delivery_root + '/' + projectid

// Illumina runfolder stats
interopdir_ilm = runfolderdir + '/InterOp'
interopdir_ctg = runfolderdir + '/ctg-interop'

// create project specific delivery dir and ctg qc dir
// -----------------------------
file(deliverydir).mkdir()
file(ctg_qc_dir).mkdir()


// set up the ctg-qc folder ( multiqc and fastqc files )
process setup_ctg_qc {
  cpus 4
  tag "$id"
  memory '32 GB'
  time '3h'

  output:
  val "x" into ctg_qc_complete_ch

  script:
  //#mqdir = ctg_qc_dir+'multiqc'
  //#fqdir  = ctg_qc_dir+'fastqc'

  """
  mv ${multiqcctgdir} ${ctg_qc_dir}
  cp -r ${fastqcdir} ${ctg_qc_dir}
  """


}



process setup_delivery {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'

  input:
  val x from ctg_qc_complet_che.collect()

  output:
  val "x" into setup_delivery_complete_ch

  script:
  """
    cp ${samplesheet} ${deliverydir}
    cp ${samplesheet_demux} ${deliverydir}

  """
  // move fastqc, star and featurecounts dirs
  """
    if [ -d ${stardir} ]; then; mv ${stardir} ${deliverydir}; fi
    if [ -d ${fastqcdir} ]; then; mv ${fastqcdir} ${deliverydir}; fi
    if [ -d ${featurecountsdir} ]; then; mv ${featurecountsdir} ${deliverydir}; fi
  """

// fastq file dilivery. if not pooled data, deliver the complete bcl2fastq directory including stats and undetermined fastq
  if ( params.deliver_fastq ){
    if ( params.pooled )
      """
        mv ${fastqdir} ${deliverydir}
      """
    else
      """
      if [ -d ${fastqdir_bcl2fastq} ]; then
        mv ${fastqdir_bcl2fastq} ${deliverydir}
      else
        mv ${fastqdir} ${deliverydir}
      fi
      """
  }
}



process multiqc_delivery {

  tag "$id"
  cpus 6
  memory '32 GB'
  time '3h'
  echo debug_mode

  input:
  val x from setup_delivery_complete_ch.collect()

  output:
  val "x" into multiqc_complete_ch

  script:
  // if a pooled run do not add interop stats to multiqc
  if ( !params.pooled)
    """
      if [ -d ${interopdir_ilm} ]; then
        # cp -r ${interopdir_ilm} ${ctg_qc_dir}/InterOp
        cp -r ${runfolderdir}/RunInfo.xml ${ctg_qc_dir}
        cp -r ${runfolderdir}/RunParameters.xml ${ctg_qc_dir}
      fi
      if [ -d ${interopdir_ctg} ]; then
        cp -r ${interopdir_ctg} ${ctg_qc_dir}/InterOp
      fi
    """
    else
    """
      cd ${deliverydir}
      multiqc -n ${projectid}_multiqc_report \\
        --interactive \\
        -o ${multiqcctgdir} .
    """
}



// ADD TO INBOX




// # if correct runfolder is specified then copy xml files
// if [ -d ${illumina_interopdir} ]; then
//   cp -r ${illumina_interopdir} ${ctg_qc_dir}/InterOp
//   cp -r ${runfolderdir}/RunInfo.xml ${ctg_qc_dir}
//   cp -r ${runfolderdir}/RunParameters.xml ${ctg_qc_dir}
// fi




// process md5sum_delivery {
//   cpus 8
//   tag "$id"
//   memory '64 GB'
//   time '3h'
// }
//
//
//



//
//
// process cleanup_projectdir {
//   cpus 8
//   tag "$id"
//   memory '64 GB'
//   time '3h'
//
// }
