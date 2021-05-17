
// Basics:
// Primer
// driver
// nf main script


// Configuration implicit variables. These are implicitly defined in the Nextflow configuration file
  // baseDir    - projectDir in v20.04.x: The directory where the main workflow script is located
  // workDir    - the work folder for nextflow temporary work-files
  // launchDir  - the directory where the workflow is run (requires version 20.04.0 or later).

  // [nas-sync/upload]: ls4 Illumina runfolder upload sync dir, shared/ctg-projects/nas-sync/upload/
  //    └──– [runfolder] = Illumina runfolder.
  //      |--- Data ....
  //      |
  //      └──-  outputdir: shared/ctg-projects/ctg-rnaseq/nf-output

  // [basedir]: ls4 ctg projets base directory, e.g. shared/ctg-projects/ctg-rnaseq
  //   └──–– [projectdir] = workdir = nf execution dir = baseDir. e.g. shared/ctg-projects/ctg-rnaseq/<projectid>
  //      |--- [fastq] (if demux)
  //      |      |- <project_id>
  //      |      |      └── fastq-files.fastq.gz
  //      |      |- Reports
  //      |      |- Stats
  //      |      └─ "Undetermined ... fastq.gz ". Remember to NOT COPY these if pooled sample
  //      |---  [bindir]: /nf-output: shared/ctg-projects/ctg-rnaseq/bin
  //      |---  "nextflow.config"
  //      |---  "ctg-rnaesq.nf"
  //      |---  "sample sheet original IEM"
  //      |---  [nfworkdir] = workDir: shared/ctg-projects/ctg-rnaseq/work; used by Nextflow
  //      └──-  [outputdir]: shared/ctg-projects/ctg-rnaseq/nf-output






//  project specific config
// -----------------------------

// root directories
project_root        =  params.project_root
delivery_root       =  params.delivery_root
qc_root             =  params.qc_root
log_root            =  params.log_root

//  project  and run folders
projectid           =  params.projectid
// n_samples           =  params.n_samples
projectdir          =  params.projectdir
//execdir             =  params.execdir
bindir              =  params.bindir


//  samplesheets
samplesheet         =  params.samplesheet
samplesheet_demux   =  params.samplesheet_demux


//  assay specific
// assay               =  params.assay
// instrument_type     =  params.instrument_type
// index_adapters      =  params.index_adapters
//paired              =  params.paired
// strandness          =  params.strandness

//  demux specific
runfolderdir        =  params.runfolderdir
runfolder           =  params.runfolder
fastqdir            =  params.fastqdir
//pooled              =  params.pooled


//  module specific
run_demux             =  params.run_demux
run_align             =  params.run_align
run_fastqc            =  params.run_fastqc
run_multiqc           =  params.run_multiqc
run_multiqc_ctg       =  params.run_multiqc_ctg
run_fastq_screen      =  params.run_fastq_screen
run_bam_indexing      =  params.bam_indexing
run_markdups          =  params.run_markdups
run_rnaseqmetrics     =  params.run_rnaseqmetrics
run_checkfiles        =  params.run_checkfiles


//  log files
logdir               =  params.logdir
logfile              =  file( logdir + '/' + projectid + '.log.complete' )




//  Other Parameters
// -----------------------------
outputdir =  projectdir+'/nf-output' // main ooutput directory for files genetated with the Pipeline

featurecountsdir = outputdir+'/featurecounts'
stardir = outputdir+'/star'
markdupsdir = outputdir+'/markdups_bam_tmp'
// qcdir = outputdir+'/qc'
fastqcdir = outputdir+'/fastqc'
markdupsqcdir = outputdir+'/markdups'
rnaseqmetricsdir = outputdir+'/rnaseqmetrics'



// ctg qc storage  dir
ctg_qc_dir      = qc_root+'/ctg-rnaseq'

// delivery directory
deliverydir = delivery_root + '/' + projectid


//  create output and logdirs
// -----------------------------
file(outputdir).mkdir()
file(fastqdir).mkdir()
//file(qcdir).mkdir()
file(fastqcdir).mkdir()
if( run_align ) file(stardir).mkdir()
if( run_align ) file(markdupsdir).mkdir()
if( run_align ) file(markdupsqcdir).mkdir()
if( run_align ) file(rnaseqmetricsdir).mkdir()
if( run_align && params.run_featurecounts ) file(featurecountsdir).mkdir()


// featurecounts
fcounts_feature     =  params.fcounts_feature

// --------------------------------------------------------------------------------------





// // Debug & test params
// // -----------------------------
// debug_mode = true
// if( debug_mode == true){
//   run_demux             =  false
//   run_fastqc            =  false
//   run_multiqc           =  false
//   run_multiqc_ctg       =  false
//   run_fastq_screen      =  false
//   bam_indexing          =  false
//   run_markdups          =  false
//   run_rnaseqmetrics     =  false
// }


//  Check paramters
// -----------------------------
if (projectid      == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet    == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}


// if not pooled - copy runfolder io-stats to project folder

//file(logdir).mkdir()


// Check if these exist
checkPathParamList = [
  project_root, delivery_root, qc_root,
  log_root, projectdir, outputdir, bindir, samplesheet
]
for (param in checkPathParamList) {
    if (param) {
	file(param, checkIfExists: true)
    }
}

// Demux specific (bcl2fastq2 )
// -----------------------------
// Check if runfolder is defined. If not set demux to false and assume that a custom fastq dir is supplied
if ( run_demux == true ){
  file(runfolderdir, checkIfExists: true)
  file(samplesheet_demux, checkIfExists: true)
}


// Define messages to print and for logfiles
def msg_startup = """\

    Workflow execution parameters
    ---------------------------------
    project id              :  ${projectid}
    project work dir        :  ${projectdir}
    nextflow execution dir  :  ${baseDir}
    nextflow output dir     :  ${outputdir}
    nextflow work dir       :  ${workDir}
    nextflow log dir        :  ${logdir}
    delivery dir            :  ${deliverydir}
    ctg qc dir              :  ${ctg_qc_dir}
    sample sheet ctg        :  ${samplesheet}
    sample sheet demux      :  ${samplesheet_demux}
   """
       .stripIndent()
println( msg_startup )


//     nf baseDir   : ${workflow.baseDir}


workflow.onComplete {

  def msg_completed = """\

  	Pipeline execution summary
  	---------------------------
  	Completed at : ${workflow.complete}
  	Duration     : ${workflow.duration}
  	Success      : ${workflow.success}
  	scriptFile   : ${workflow.scriptFile}
    exit status  : ${workflow.exitStatus}
  	errorMessage : ${workflow.errorMessage}
  	errorReport  :
  	"""
  	.stripIndent()
  def error = """\
		${workflow.errorReport}
	   """


  // base = csv.getBaseName()
  logfile.text = msg_startup.stripIndent()
  logfile.append( msg_completed.stripIndent() )
  logfile.append( error )

  println( msg_completed )
}


def msg_modules = """\

    Run modules
    ---------------------------------
    demux       :  ${run_demux}
    fastqc      :  ${run_fastqc}

   """
   .stripIndent()

println( msg_modules )




// =====================
//  Create input channels
// =====================

// all samplesheet info
// if ( params.paired == true ) {
Channel
  .fromPath(samplesheet)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.fastq_1, row.fastq_2, row.Species ) }
  .tap{ infoall }
  .set { fastq_ch }

Channel
  .fromPath(samplesheet)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.bam, row.Strandness, row.Species) }
  .tap { infobam }
  .into { bam_checkbam_ch; bam_indexbam_ch; bam_rnaseqmetrics_ch; bam_markdups_ch; bam_featurecounts_ch }


    // .set { fastq_ch }
  // println("running params.paired")
// }
// else {
//   Channel
//       .fromPath(samplesheet)
//       .splitCsv(header:true)
//       .map { row -> tuple( row.Sample_ID, row.fastq_1) }
//       .tap{ infoall }
//       .set { fastqc_ch }
//   }

/// Channels input files


println " > Samples to process: "
println "[Sample_ID,fastq1,fastq2]"
infoall.subscribe { println "Info: $it" }
println " > Projects to process : "
//println "[agg,Sample_Project]"
//infoProject.subscribe { println "Info Projects: $it" }





// Run bcl2fastq if run_demux
process bcl2fastq {
  // -w must be lower than number of samples
  publishDir "${fastqdir}", mode: 'copy', overwrite: 'true'
  cpus 4
  tag "$id"
  memory '8 GB'
  time '3h'

  when:
  run_demux

  input:
  val samplesheet_demux

  output:
  val "x" into fastq_check_ch

  script:
  """
   bcl2fastq -R ${runfolderdir} \\
            --sample-sheet ${samplesheet_demux} \\
            --no-lane-splitting  \\
            -r 1 \\
            -p ${task.cpus}  \\
            -w 1  \\
            --output-dir ${fastqdir}
   """
}

// Channel to start count if demux == 'n'
// Projects
if ( run_demux == false ) {
   Channel
	 .from("x")
   .set{ fastq_check_ch }
}





process checkfiles_fastq {
  // Run fastqc. Also check if all expected files, defined in the ctg samplesheet, are present in fastqdir
  tag "$id"
  cpus 1
  memory '5 GB'
  time '3h'
  echo true

  input:
  val x from fastq_check_ch.collect()
  set sid, read1, read2, species from fastq_ch

  output:
  set
  val "x" into run_star_ch
  set sid, read1, read2, species into fastqc_ch
  set sid, read1, read2, species into star_ch

  //input:
  //  set file(samples_csv) from sheet_ctg_ch
  script:

  if( params.paired && run_checkfiles )
    """
      echo "running fastqc in paired reads "
      if [ ! -f ${fastqdir}/${read1} ]; then
        echo "Warning: Cannot locate fastq_1 file ${fastqdir}/${read1}"
        exit 2
      fi

      if [ ! -f ${fastqdir}/${read1} ]; then
        echo "Warning: Cannot locate fastq_2 file ${fastqdir}/${read2}"
        exit 2
      fi
    """
  else if ( params.paired == false && run_checkfiles)
    """
      echo "running fastqc in non paired mode "

      if [ ! -f ${fastqdir}/${read1} ]; then
        "Cannot locate fastq_1 file ${read2}"
        exit 2
      fi
    """
  else
    """
      echo "file check overridden"
    """

}




/* ===============================================================
  *      -- ALIGMENT SECTION --
  =============================================================== */

  // Run STAR: ls4 ctg projets base directory, e.g. shared/ctg-projects/ctg-rnaseq
  //   └──–– check_bam: uses sample sheet to check if expected bams are generated
  //      |-  index_bam : optional
  //      |-  markdups  : optional
  //      └── rnaseqmetrics : optional
  //   The three latter are always run to generate flag - but may be run with no script if set ti false
  //    When all three are run, process




// Set align chanel complete is !run_align
if ( run_align == false ) {
   Channel
	 .from("x")
   .set{ align_complete_ch }
}



// Run STAR
process star  {
  tag "$id"
  cpus 1
  memory '5 GB'
  time '3h'
  echo true
  publishDir "${stardir}", mode: 'copy', overwrite: true

  input:
  val x from run_star_ch.collect()
  set sid, read1, read2, species from star_ch // from checkfiles_fastq

  output:
  val "x" into checkbam_ch
  // file "${sid}_Aligned.sortedByCoord.out.bam" into bam_featurecounts_ch // channel defined start instead

  when:
  run_align

  script:
  if ( species == "Homo sapiens" ){
    genome=params.star_genome_hs }
  else if ( species == "Mus musculus" ){
    genome=params.star_genome_mm }
  else{
    genome = ""
    println( "Warning: Species not recognized." )}

  """
#  STAR --genomeDir ${genome} \\
#    --readFilesIn ${fastqdir}/${read1} ${fastqdir}/${read2} \\
#    --runThreadN ${task.cpus}  \\
#    --readFilesCommand zcat \\
#    --outSAMtype BAM SortedByCoordinate \\
#    --limitBAMsortRAM 10000000000 \\
#    --outFileNamePrefix ${sid}_
  """

}




// Check STAR bam files against names in sample sheet
// checkfiles_star
process check_bam {
  // Run fastqc. Also check if all expected files, defined in the ctg samplesheet, are present in fastqdir
  tag "$id"
  cpus 1
  memory '1 GB'
  time '3h'
  echo true

  input:
  val x from checkbam_ch.collect() // checkbam_ch - when star is completed
  set sid, bam, strand, species from bam_checkbam_ch

  output:
  //set sid, bam into bam_markdups_ch
  //set sid, bam into bam_index_ch
  val "x" into indexbam_ch
  val "x" into rnaseqmetrics_ch
  val "x" into markdups_ch
  val "x" into featurecounts_ch

  //input:
  //  set file(samples_csv) from sheet_ctg_ch
  script:
  if( run_checkfiles )
    """
      if [ ! -f ${stardir}/${bam} ]; then
        echo "Warning: Cannot locate fastq_1 file ${stardir}/${bam}"
        exit 2
      fi
    """
  else
    """
    echo "file check overridden"
    """
}


// samtools index bamfile
// ml Java; ml nextflow/19.04.1
// ml Singularity
// ml GCC/7.3.0-2.30
// ml SAMtools/1.9
// samtools index bamfile

process index_bam {
  tag "$id"
  cpus 1
  memory '1 GB'
  time '3h'
  echo true
  publishDir "${stardir}", mode: 'copy', overwrite: 'true'

  input:
  val x from indexbam_ch.collect()
  set sid, bam, strand, species from bam_indexbam_ch


  output:
  val "x" into indexbam_complete_ch

  when:
  run_align

  script:
  if ( run_bam_indexing )
    """
    echo "${stardir}/${bam}"
    # samtools index -bc ${stardir}/${bam}
    """
  else
    """
    """
}


// picard mark duplicates

process markdups {
  tag "$id"
  cpus 1
  memory '1 GB'
  time '3h'
  echo true

  input:
  val x from markdups_ch.collect()
  set sid, bam, strand, species from bam_markdups_ch

  output:
  val "x" into markdups_complete_ch
  // val "x" into move

  when:
  run_align

  script:
  if (run_markdups)
  """
  echo "bam: ${bam}"
  echo "markdupsdir: ${markdupsdir}/${bam}"
    #java -jar /usr/local/bin/picard.jar MarkDuplicates \\
    #  INPUT=${stardir}/${bam} \\
    #  OUTPUT=${markdupsdir}/${bam} \\
    #  METRICS_FILE=${markdupsqcdir}/${sid}_bam.MarkDuplicates.metrics.txt \\
    #  TAGGING_POLICY=All \\
    #  REMOVE_DUPLICATES=false \\
    #  ASSUME_SORTED=true \\
    #  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000 \\
    #  QUIET=true \\
    #  VERBOSITY=WARNING

    # mv -f ${markdupsdir}/${bam} ${stardir}/${bam}
  """
  else
  """
  """
}


process rnaseqmetrics {
  tag "$id"
  cpus 1
  memory '1 GB'
  time '3h'
  echo true

  input:
  val x from rnaseqmetrics_ch.collect()
  set sid, bam, strand, species from bam_rnaseqmetrics_ch

  output:
  val "x" into rnaseqmetrics_complete_ch

  when:
  run_align

  script:
  // NONE, FIRST_READ_TRANSCRIPTION_STRAND, and SECOND_READ_TRANSCRIPTION_STRAND.
  if ( strand == "forward" )
    strand="FIRST_READ_TRANSCRIPTION_STRAND"
  else if ( strand == "reverse" )
    strand="SECOND_READ_TRANSCRIPTION_STRAND"
  else
    strand="NONE"


  if ( species == "Homo sapiens" ){
    refflat = params.picard_refflat_hs
    rrna = params.picard_rrna_hs}
  else if ( species == "Mus musculus" ){
    refflat = params.picard_refflat_mm
    rrna = params.picard_rrna_mm}
  else{
    refflat = ""
    rrna = ""
  }



  if ( run_rnaseqmetrics )
  """
  #  if [ "${species}" == "Homo sapiens" ]; then
  #    refflat=${params.picard_refflat_hs}
  #    rrna=${params.picard_rrna_hs}
  #  elif [ "${species}" == "Mus musculus" ]; then
  #    refflat=${params.picard_refflat_mm}
  #    rrna=${params.picard_rrna_mm}
  #  else
  #    echo "Warning: rnaseqmetrics, Species not recognized"
  #    refflat="APAs"
  #    rrna="nbajsss"
  #  fi

    echo "strand: ${strand}"

    echo "rrna file: ${rrna}"
    echo "refflat file: ${refflat}"

      ## java -jar /usr/local/bin/picard.jar CollectRnaSeqMetrics \\
      ##  INPUT=${stardir}/${bam} \\
      ##  OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
      ##  REF_FLAT=\${refflat} \\
      ##  STRAND=${strand} \\
      ##  RIBOSOMAL_INTERVALS=\${rrna}
    """
  else
    """
    """

}



// featureCounts on bam-fililes

// NOTE: featurecounts will use the SAME Starndness for all samples, ie.e paaram.strandness not individual sample strandness
// NOTE p flag: assumes paied (not applicabe if not paired data)

process featurecounts {
  tag "$id"
  cpus 1
  memory '1 GB'
  time '3h'
  echo true

	input:
  val x from featurecounts_ch.collect()
	file bams from bam_featurecounts_ch.collect()

	output:
	val "x" into featurecounts_complete_ch

	when:
	run_align

  script:
  // Global settings - for ALL Samples
  if( params.strandness == "forward" )
    strand_numeric = 1
  else if ( params.strandness == "reverse" )
    strand_numeric = 2
  else
    strand_numeric = 0

  // gtf used for featurecounts
  if ( params.species_global == "Homo sapiens" ){
    gtf = params.gtf_hs}
  else if  ( params.species_global == "Mus musculus" ){
    gtf = params.gtf_mm}
  else{
    gtf=""}



  if( params.run_featurecounts )
    """
      echo "gtf: ${gtf}"
      ## featureCounts -T ${task.cpus} \\
      ##  -t ${params.fcounts_feature} \\
      ##  --extraAttributes gene_name,gene_type \\
      ##  -a ${gtf} -g gene_id  \\
      ##  -o ${featurecountsdir}/${projectid}_geneid.featureCounts.txt \\
      ##  -p \\
      ##  -s ${strand_numeric} ${bams}

    """
  else
    """
    """


}


// Collect processes and prepare for MultiQC
process collect_align {
  tag "$id"
  cpus 1
  memory '5 GB'
  time '3h'
  // echo true

  input:
  val x from indexbam_complete_ch.collect()
  val x from markdups_complete_ch.collect()
  val x from rnaseqmetrics_complete_ch.collect()
  val x from featurecounts_complete_ch.collect()

  output:
  val "x" into align_complete_ch

  when:
  run_align

  script:
  """
  """

}





/* ===============================================================
  *      FASTQC
  =============================================================== */

// Customer multi QC - not same as CTG multiQC
// Customer multi QC should be run on the delivery dir?

// Pooled sequencing run - will not provide sequen


process fastqc {
  // Run fastqc. Also check if all expected files, defined in the ctg samplesheet, are present in fastqdir
  tag "$id"
  cpus 1
  memory '5 GB'
  time '3h'
  // echo true

  input:
  set sid, read1, read2, species from fastqc_ch  // from check fastq

  output:
  val "x" into multiqc_ctg_ch

  when:
  run_fastqc

  //input:
  //  set file(samples_csv) from sheet_ctg_ch
  script:
  if(params.paired)
    """
      echo "running fastqc in paired reads mode"
      ## fastqc ${fastqdir}/${read1} ${fastqdir}/${read2}  --outdir ${fastqcdir}
  """
  else
    """
      echo "running fastqc in non paired reads mode "
      ## fastqc ${fastqdir}/${read1}  --outdir ${fastqcdir}
    """

}



/* ===============================================================
  *      MULTIQC CTG - ON NF DIRECTORY
  =============================================================== */
// This is not carried on to the customer.



process multiqc_ctg {
  //publishDir "${multiqcctgdir}", mode: 'copy', overwrite: 'true'
  tag "$id"
  cpus 1
  memory '5 GB'
  time '3h'
  echo true

  input:
  val x from multiqc_ctg_ch.collect()

  output:
  val "x" into multiqc_ctg_complete_ch

  when:
  run_multiqc_ctg

  script:
  """
    ##  cd ${outputdir}
    ## multiqc -n ${projectid}_multiqc_report \\
    ## --interactive \\
    ## -o ${multiqcctgdir} .
  """

}


// process move_files {
//
//   tag "$id"
//   cpus 1
//   memory '5 GB'
//   time '3h'
//   echo true
//
//   input:
//   val x from multiqc_ctg_complete_ch.collect()
//
//   output:
//   val "x" into deliverydir_comple
//
//
//   script:
//
//
//
//   """
//
//
//   """
//
// }


 // process multiqc_delivery {
 //
 //
 //
 // }
