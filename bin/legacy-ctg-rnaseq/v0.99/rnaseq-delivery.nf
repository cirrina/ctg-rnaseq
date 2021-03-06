/* ===============================================================
  *      -- RNASEQ DELIVERY --
  =============================================================== */


// get and define all directory variablres from  project specific config
// ----------------------------------------------------------------------

// root directories
project_root        =  params.project_root
delivery_root       =  params.delivery_root
qc_root             =  params.qc_root
completed_root      =  params.completed_root
//log_root            =  params.log_root

//  project  and run folders
projectid           =  params.projectid
projectdir          =  params.projectdir

//  samplesheets
samplesheet         =  params.samplesheet
samplesheet_demux   =  params.samplesheet_demux
samplesheet_original =  params.samplesheet_original

//  demux specific
runfolderdir        =  params.runfolderdir
fastqdir            =  params.fastqdir
fastqdir_bcl2fastq  =  params.fastqdir_bcl2fastq
deliver_raw         =  params.deliver_raw
deliver_fastq       =  params.deliver_fastq

// extras
outboxdir           =  params.outboxdir


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

ctg_qc_dir          =  qc_root + '/' + projectid
deliverydir         =  delivery_root + '/' + projectid
multiqcdeliverydir  =  deliverydir+'/multiqc'
completeddir        =  completed_root + '/' + projectid

// Illumina runfolder stats
interopdir_ilm = runfolderdir + '/InterOp'
interopdir_ctg = runfolderdir + '/ctg-interop'

// create project specific delivery dir and ctg qc dir
// -----------------------------
file(deliverydir).mkdir()
file(ctg_qc_dir).mkdir()

// readme deliverydir
readme = deliverydir +'/README_ctg_delivery_' + projectid



//
debug_mode = false // will turn echo to true if applicaple





/* ===============================================================
  *      -- START PROCESSES --
  =============================================================== */
//
//
// // set up the ctg-qc folder ( multiqc and fastqc files )
// process setup_ctg_qc {
//   cpus 4
//   tag "$id"
//   memory '32 GB'
//   time '3h'
//
//   output:
//   val "x" into ctg_qc_complete_ch
//
//   script:
//   //#mqdir = ctg_qc_dir+'multiqc'
//   //#fqdir  = ctg_qc_dir+'fastqc'
//
//   """
//   if [ -d ${multiqcctgdir} ]
//     then
//     mv ${multiqcctgdir} ${ctg_qc_dir}
//   fi
//   if [ -d ${fastqcdir} ]
//     then
//     cp -r ${fastqcdir} ${ctg_qc_dir}
//   fi
//   """
// }
//
//
//
// process setup_delivery {
//   cpus 8
//   tag "$id"
//   memory '64 GB'
//   time '3h'
//
//   input:
//   val x from ctg_qc_complete_ch.collect()
//
//   output:
//   val "x" into setup_delivery_complete_ch
//
//   script:
//   """
//     cp ${samplesheet} ${deliverydir}
//     if [ -f ${samplesheet_demux} ]; then
//      cp ${samplesheet_demux} ${deliverydir}
//     fi
//
//     if [ -d ${stardir} ]; then
//       mv ${stardir} ${deliverydir}
//     fi
//
//     if [ -d ${featurecountsdir} ]; then
//       mv ${featurecountsdir} ${deliverydir}
//     fi
//
//     if [ -d ${fastqcdir} ]; then
//      mv ${fastqcdir} ${deliverydir}
//     fi
//   """
// }
//
// process move_fastq {
//
//   cpus 8
//   tag "$id"
//   memory '64 GB'
//   time '3h'
//   echo debug_mode
//
//   input:
//   val x from setup_delivery_complete_ch.collect()
//
//   output:
//   val "x" into move_fastq_complete_ch
//   // fastq files.
//   // if not pooled deliver the complete bcl2fastq directory including stats and undetermined fastq
//
//   script:
//   if ( params.deliver_fastq ){
//     if ( params.pooled )
//       """
//         mkdir -p ${deliverydir}/fastq
//         if [ -d ${fastqdir} ]; then
//           echo "pooled data. moving fastq foldler only."
//           mv ${fastqdir} ${deliverydir}/fastq
//         fi
//       """
//     else
//       """
//       if [ -d ${fastqdir_bcl2fastq} ]; then
//         echo "non pooled data. moving comlplete bcl2fastq output foldler."
//         mv ${fastqdir_bcl2fastq} ${deliverydir}
//       elif [ -d ${fastqdir} ]; then
//         echo "non pooled data. moving fastq foldler only."
//         mkdir -p ${deliverydir}/fastq
//         mv ${fastqdir} ${deliverydir}/fastq
//       fi
//       """
//   }
//   else
//   """
//   """
// }
//
//
//
// process multiqc_delivery {
//
//   tag "$id"
//   cpus 6
//   memory '32 GB'
//   time '3h'
//   echo debug_mode
//
//   input:
//   val x from move_fastq_complete_ch.collect()
//
//   output:
//   val "x" into multiqc_complete_1_ch
//   val "x" into multiqc_complete_2_ch
//
//   script:
//   mqcreport = ${projectid} + '_multiqc_report'
//
//   if (! new File( mqcreport+'.html' ).exists() )
//     """
//       cd ${deliverydir}
//       multiqc -n ${mqcreport} \\
//         --interactive \\
//         -o ${multiqcdeliverydir} .
//     """
//   else
//   """
//     echo "${mqcreport} already exists - skipping"
//   """
// }
//
//
//
//
// process md5sum_delivery {
//   cpus 8
//   tag "$id"
//   memory '64 GB'
//   time '3h'
//
//   input:
//   val x from multiqc_complete_2_ch.collect()
//
//   when:
//   params.run_md5sum
//
//   script:
//   md5sumfile= ${deliverydir} + '/md5sum.txt'
//
//   if (! new File( md5sumfile ).exists() )
//   """
//     find ${deliverydir} -type f -exec md5sum '{}' \\; > ${md5sumfile} ; echo
//   """
//   else
//   """
//     echo "${md5sumfile} already exists. skipping."
//   """
// }
//
//
//
//
//
// /// provess add README with dir size
// process genereate_readme {
//   cpus 2
//   tag "$id"
//   memory '16 GB'
//   time '3h'
//
//   input:
//   val x from multiqc_complete_1_ch.collect()
//
//   output:
//   val "x" into genereate_readme_complete_ch
//
//   script:
//
//   """
//   cd ${deliverydir}
//   echo "CTG Delivery"                        > $readme
//   echo "Project:   ${projectid}"             >> $readme
//   du -ch -d 0 . | grep 'total' > $readme
//   """
//
// }


// // ADD TO OUTBOX
// process add_outbox {
//   tag "$id"
//   cpus 6
//   memory '32 GB'
//   time '3h'
//   echo debug_mode
//
//   input:
//   val x from genereate_readme_complete_ch.collect()
//
//   output:
//   val "x" into add_outbox_complete_ch
//
//
//   script:
//
//   if ( params.copy_to_outbox ){
//     """
//     # userid=\$(whoami)
//     userid="percebe"
//     boxdir="/box/outbox/percebe/${projectid}"
//
//     mkdir -p \${boxdir}
//     cp  -r ${multiqcdeliverydir} \${boxdir}
//
//     """}
//   else{
//     """
//     """}
//
// }

/* ===============================================================
  *      ADD TO OUTBOX FOR CONVENIANT DOWNLOAD
  =============================================================== */

// I CAN NOT GET OUTBOX TO MOUNT PROPERLY TO SINGULARITY CONTAINER.
// IF Copy stuff to outbox - then add this to shell script instead.
// process add_outbox {
//   tag "$id"
//   cpus 6
//   memory '32 GB'
//   time '3h'
//   echo debug_mode
//
//   input:
//   val x from genereate_readme_complete_ch.collect()
//
//   output:
//   val "x" into add_outbox_complete_ch
//
//
//   script:
//
//   if ( params.copy_to_outbox ){
//     """
//     mkdir -p ${outboxdir}
//     cp -r ${multiqcdeliverydir} ${outboxdir}
//
//     """}
//   else{
//     """
//     """}
//
// }


// SETUP COMPLETED DIR
// CLEANUP OF PROJECT DIR
// everythiing must now have been copied or moved from the project dir
// should include sample sheets, nextflow scripts,
// rscirpt log from runFolder
// .log.completed

// process setup_completed {
//   tag "$id"
//   cpus 2
//   memory '16 GB'
//   time '1h'
//   echo debug_mode
//
//   input:
//   val x from genereate_readme_complete_ch.collect()
//
//   output:
//   val "x" into cleanup_complete_ch
//
//
//   script:
//   if ( params.run_cleanup )
//     """
//     mkdir -p ${completeddir}
//     cd ${projectdir}
//
//     if [[ -f "${runfolderdir}/iem.rscript.log" ]]; then
//       cp ${samplesheet_original} ${completeddir}
//     fi
//     if [ -f ${} ]; then
//       cp ${samplesheet_original} ${completeddir}
//     fi
//
//     if [ -f ${samplesheet_demux} ]; then
//       mv ${samplesheet_demux} ${completeddir}
//     fi
//     if [ -f ${samplesheet} ]; then
//       mv ${samplesheet} ${completeddir}
//     fi
//     if [ -f ${samplesheet} ]; then
//       mv ${samplesheet} ${completeddir}
//     fi
//     mv ./rnaseq-delivery.nf ${completeddir}
//     mv ./rnaseq-main.nf ${completeddir}
//     mv ./nextflow.config ${completeddir}
//     mv ./${projectid}.log.complete ${completeddir}
//     mv ./bin ${completeddir}
//     mv ./nextflow.params.${projectid} ${completeddir}
//     """
//   else
//     """
//     """
//
// }
//




/* ===============================================================
  *     ++++ CLEANUP AND FINALIZATION ++++
  =============================================================== */

// move/copy relevant files and folders tpo completed dir
// ctg completed will contain scripts and sample sheeets efc for all runs

process setup_completed {
  tag "$id"
  cpus 2
  memory '16 GB'
  time '3h'
  echo debug_mode

  input:
  val x from genereate_readme_complete_ch.collect()

  output:
  val "x" into cleanup_complete_ch


  script:
  if ( params.run_cleanup ){
    if (! new File( runfolderdir + '/iem.rscript.log' ).exists() ){
        myFile.copyTo( completeddir + '/iem.rscript.log')
    }



    """

    cd ${projectdir}

    ## sample sheets
    if [[ -f "${samplesheet_original}" ]]; then
      cp ${samplesheet_original} ${completeddir}
    fi
    if [[ -f ${samplesheet_demux} ]]; then
      cp ${samplesheet_demux} ${completeddir}
    fi
    if [ -f ${samplesheet} ]; then
      cp ${samplesheet} ${completeddir}
    fi

    ## nextflow sripts and configs
    if [ -f ${} ]; then
      cp ${} ${completeddir}
    fi

    ## logs
    if [[ -f "${runfolderdir}/iem.rscript.log" ]]; then
      cp ${samplesheet_original} ${completeddir}
    fi

    if [ -f ${samplesheet} ]; then
      cp ${samplesheet} ${completeddir}
    fi
    cp ./rnaseq-delivery.nf ${completeddir}
    cp ./rnaseq-main.nf ${completeddir}
    cp ./nextflow.config ${completeddir}
    cp ./${projectid}.log.complete ${completeddir}
    cp ./bin ${completeddir}
    cp ./nextflow.params.${projectid} ${completeddir}

    """
  }else{
    """
    """
  }
}


// CLEANUP OF PROJECT FOLDER
//  move all keeper-fiels into one single


// // final chmods
// process final_chmods {
//   tag "$id"
//   cpus 2
//   memory '16 GB'
//   time '3h'
//   echo debug_mode
//
//   input:
//   val x from genereate_readme_complete_2_ch.collect()
//
//   output:
//   val "x" into final_chmods_ch
//
//
//   script:
//
//   if ( params.copy_to_outbox )
//   """
//   chmod 770 -R ${deliverydir}
//   chmod 770 -R ${ctg-qc-dir}
//   """
//   else
//   """
//   """
//
// }





//
//
// process cleanup_projectdir {
//   cpus 8
//   tag "$id"
//   memory '64 GB'
//   time '3h'
//
// }
