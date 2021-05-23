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
multiqcdeliverydir   = deliverydir+'/multiqc'

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
  cp ${fastqcdir} ${ctg_qc_dir}
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
  // move stardir if present
  """
    if [ -d ${stardir} ]; then
      mv ${stardir} ${deliverydir}
    fi
  """
  // move featurecounts dir
  """
    if [ -d ${featurecountsdir} ]; then; mv ${featurecountsdir} ${deliverydir}; fi
  """
  // fastqc dir
  """
    if [ -d ${fastqcdir} ]; then; mv ${fastqcdir} ${deliverydir}; fi
  """

  // fastq files.
  // if not pooled deliver the complete bcl2fastq directory including stats and undetermined fastq
  if ( params.deliver_fastq ){
    if ( params.pooled )
      """
        mkdir ${deliverydir}/fastq
        mv ${fastqdir} ${deliverydir}/fastq
      """
    else
      """
      if [ -d ${fastqdir_bcl2fastq} ]; then
        mv ${fastqdir_bcl2fastq} ${deliverydir}
      else
        mkdir ${deliverydir}/fastq
        mv ${fastqdir} ${deliverydir}/fastq
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
  val "x" into multiqc_complete_1_ch
  val "x" into multiqc_complete_2_ch

  script:
    """
      cd ${deliverydir}
      multiqc -n ${projectid}_multiqc_report \\
        --interactive \\
        -o ${multiqcdeliverydir} .
    """
}



// ADD TO OUTBOX
process add_outbox {
  tag "$id"
  cpus 6
  memory '32 GB'
  time '3h'
  echo debug_mode

  input:
  val x from multiqc_complete_1_ch.collect()
  when: deliver_outbox

  """
  mkdir ~/outbox/${projectid}
  cp  ${multiqcdeliverydir} ~/outbox/${projectid}
  cp  ${deliverydir}/fastqc ~/outbox/${projectid}
  cp  ${ctg_qc_dir}/multiqc_ctg ~/outbox/${projectid}

  """

}

process md5sum_delivery {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'

  input:
  val x from multiqc_complete_2_ch.collect()

  script:
  """
  if [ -d ${fastqdir_bcl2fastq} ]; then
    cd ${deliverydir}
    find -type f -exec md5sum '{}' \; > md5sum.txt ; echo
  fi
  """
}






// # if correct runfolder is specified then copy xml files
// if [ -d ${illumina_interopdir} ]; then
//   cp -r ${illumina_interopdir} ${ctg_qc_dir}/InterOp
//   cp -r ${runfolderdir}/RunInfo.xml ${ctg_qc_dir}
//   cp -r ${runfolderdir}/RunParameters.xml ${ctg_qc_dir}
// fi







//
//
// process cleanup_projectdir {
//   cpus 8
//   tag "$id"
//   memory '64 GB'
//   time '3h'
//
// }
