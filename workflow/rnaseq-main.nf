
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
  //      └──-  outputdir: shared/ctg-projects/rnaseq/nf-output

  // [basedir]: ls4 ctg projets base directory, e.g. shared/ctg-projects/rnaseq
  //   └──–– [projectdir] = workdir = nf execution dir = baseDir. e.g. shared/ctg-projects/rnaseq/<projectid>
  //      |--- [fastq] (if demux)
  //      |      |- <project_id>
  //      |      |      └── fastq-files.fastq.gz
  //      |      |- Reports
  //      |      |- Stats
  //      |      └─ "Undetermined ... fastq.gz ". Remember to NOT COPY these if pooled sample
  //      |---  [bindir]: /nf-output: shared/ctg-projects/rnaseq/bin
  //      |---  "nextflow.config"
  //      |---  "ctg-rnaesq.nf"
  //      |---  "sample sheet original IEM"
  //      |---  [nfworkdir] = workDir: shared/ctg-projects/rnaseq/work; used by Nextflow
  //      └──-  [outputdir]: shared/ctg-projects/rnaseq/nf-output




//  project specific config
// -----------------------------

// root directories
project_root        =  params.project_root
delivery_root       =  params.delivery_root
qc_root             =  params.qc_root
completed_root      =  params.completed_root
//log_root            =  params.log_root

//  project  and run folders
projectid           =  params.projectid
projectdir          =  params.projectdir
bindir              =  params.bindir
// n_samples           =  params.n_samples
//execdir             =  params.execdir


//  samplesheets
samplesheet         =  params.samplesheet
samplesheet_demux   =  params.samplesheet_demux

//  demux specific
runfolderdir        =  params.runfolderdir
runfolder           =  params.runfolder
fastqdir            =  params.fastqdir
fastqdir_bcl2fastq  =  params.fastqdir_bcl2fastq
//pooled              =  params.pooled

outboxdir           = params.outboxdir

//  module specific
run_demux             =  params.run_demux
run_align             =  params.run_align
run_fastqc            =  params.run_fastqc
run_multiqc           =  params.run_multiqc
run_multiqc_ctg       =  params.run_multiqc_ctg
run_fastqscreen       =  params.run_fastqscreen
run_bam_indexing      =  params.run_bam_indexing
run_markdups          =  params.run_markdups
run_rnaseqmetrics     =  params.run_rnaseqmetrics
run_checkfiles        =  params.run_checkfiles
run_featurecounts     =  params.run_featurecounts

//  log files
logdir               =  params.logdir
logfile              =  file( logdir + '/' + projectid + '.log.complete' )





//  Other Fixed Directories
// -----------------------------
outputdir =  projectdir+'/nf-output' // main ooutput directory for files genetated with the Pipeline

featurecountsdir = outputdir+'/featurecounts'
stardir = outputdir+'/star'
markdupsdir = outputdir+'/markdups_bam_tmp'
fastqcdir = outputdir+'/fastqc'
markdupsqcdir = outputdir+'/markdups'
rnaseqmetricsdir = outputdir+'/rnaseqmetrics'
multiqcctgdir = outputdir+'/multiqc-ctg'
fastqscreendir = outputdir+'/fastqscreen'


//  create output and logdirs
// -----------------------------
file(outputdir).mkdir()

if ( run_demux ) file(fastqdir_bcl2fastq).mkdir()
if ( run_demux ) file(fastqdir).mkdir()

//file(qcdir).mkdir()
file(fastqcdir).mkdir()
file(multiqcctgdir).mkdir()

if( run_align ) file(stardir).mkdir()
if( run_align ) file(markdupsdir).mkdir()
if( run_align ) file(markdupsqcdir).mkdir()
if( run_align ) file(rnaseqmetricsdir).mkdir()
if( run_align && run_featurecounts ) file(featurecountsdir).mkdir()
if( run_fastqscreen ) file(fastqscreendir).mkdir()
if( copy_to_outbox ) file(outboxdir).mkdir()

// featurecounts
fcounts_feature     =  params.fcounts_feature

// --------------------------------------------------------------------------------------



// Check if files and directories exist
checkPathParamList = [
  project_root, delivery_root, qc_root,
  projectdir, bindir,
  fastqdir, logdir, outputdir,
  samplesheet
]
for (param in checkPathParamList) {
    if (param) {
	file(param, checkIfExists: true)
    }
}
if ( run_demux ) {
  file(fastqdir_bcl2fastq, checkIfExists: true)
  file(samplesheet_demux, checkIfExists: true)
}


// // Debug & test params
// // -----------------------------
debug_mode = false // will turn echo to true


//  Check paramters
// -----------------------------
if (projectid      == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet    == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}


// if not pooled - copy runfolder io-stats to project folder

//file(logdir).mkdir()



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
    sample sheet ctg        :  ${samplesheet}
    sample sheet demux      :  ${samplesheet_demux}
   """
       .stripIndent()
println( msg_startup )



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
  .into { bam_checkbam_ch; bam_indexbam_ch; bam_rnaseqmetrics_ch; bam_markdups_ch }

Channel
    .fromPath(samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.bam) }
    .tap{ infoallfcounts }
    .set { bam_featurecounts_ch }

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
  // publishDir "${fastqdir}", mode: 'copy', overwrite: 'true'
  cpus 4
  tag "$id"
  memory '110 GB'
  time '24h'

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
            --output-dir ${fastqdir_bcl2fastq}
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
  echo debug_mode

  input:
  val x from fastq_check_ch.collect()
  set sid, read1, read2, species from fastq_ch

  output:
  val "x" into run_star_ch
  set sid, read1, read2, species into fastqc_ch
  set sid, read1, read2, species into star_ch
  set sid, read1, read2, species into fastqscreen_ch

  //input:
  //  set file(samples_csv) from sheet_ctg_ch
  script:

  if( params.paired && run_checkfiles )
    """
      echo "running fastqc with paired reads "
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

  // Run STAR: ls4 ctg projets base directory, e.g. shared/ctg-projects/rnaseq
  //   └──–– check_bam: uses sample sheet to check if expected bams are generated
  //      |-  index_bam : optional
  //      |-  markdups  : optional
  //      └── rnaseqmetrics : optional
  //   The three latter are always run to generate flag - but may be run with no script if set ti false
  //    When all three are run, process




// Set align chanel complete is !run_align. It is now possible to run processes  downstream star...
//if ( run_align == false ) {
//   Channel
//	 .from("x")
//   .set{ align_complete_ch }
//}



// Run STAR
process star  {
  tag "$id"
  cpus 16
  memory '100 GB'
  time '36h'
  echo debug_mode
  //publishDir "${stardir}", mode: 'copy', overwrite: true

  input:
  val x from run_star_ch.collect()
  set sid, read1, read2, species from star_ch // from checkfiles_fastq

  output:
  val "x" into checkbam_ch
  // file "${sid}_Aligned.sortedByCoord.out.bam" into bam_featurecounts_ch // channel defined start instead

  //when:
  //run_align

  script:
  if ( species == "Homo sapiens" ){
    genome=params.star_genome_hs }
  else if ( species == "Mus musculus" ){
    genome=params.star_genome_mm }
  else{
    genome = ""
    println( "Warning: Species not recognized." )}

  if ( params.paired ){
      starfiles = "${fastqdir}/${read1} ${fastqdir}/${read2}" }
  else{
      starfiles = "${fastqdir}/${read1}" }


  if ( run_align )
  """
  STAR --genomeDir ${genome} \\
    --readFilesIn ${starfiles} \\
    --runThreadN ${task.cpus}  \\
    --readFilesCommand zcat \\
    --outSAMtype BAM SortedByCoordinate \\
    --limitBAMsortRAM 10000000000 \\
    --outFileNamePrefix ${stardir}/${sid}_
  """
  else
  """
  """

}




// Check STAR bam files against names in sample sheet
// checkfiles_star
process checkfiles_bam {
  // Run fastqc. Also check if all expected files, defined in the ctg samplesheet, are present in fastqdir
  tag "$id"
  cpus 1
  memory '1 GB'
  time '1h'
  echo debug_mode

  input:
  val x from checkbam_ch.collect() // checkbam_ch - when star is completed
  set sid, bam, strand, species from bam_checkbam_ch

  output:
  //set sid, bam into bam_markdups_ch
  //set sid, bam into bam_index_ch
  //val "x" into indexbam_ch
  val "x" into rnaseqmetrics_ch
  //val "x" into markdups_ch
  //val "x" into featurecounts_ch

  //input:
  //  set file(samples_csv) from sheet_ctg_ch
  script:
  if( run_checkfiles )
    """
      if [ ! -f ${stardir}/${bam} ]; then
        echo "Warning: Cannot locate bam file ${stardir}/${bam}"
        exit 2
      fi
    """
  else
    """
    echo "file check overridden"
    """
}




process rnaseqmetrics {
  tag "$id"
  cpus 8
  memory '48 GB'
  time '24h'
  echo debug_mode

  input:
  val x from rnaseqmetrics_ch.collect()
  set sid, bam, strand, species from bam_rnaseqmetrics_ch

  output:
  val "x" into rnaseqmetrics_complete_ch

  //when:
  //run_align

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
    echo "strand: ${strand}"
    echo "rrna file: ${rrna}"
    echo "refflat file: ${refflat}"

    java -jar /usr/local/bin/picard.jar CollectRnaSeqMetrics \\
      INPUT=${stardir}/${bam} \\
      OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
      REF_FLAT=${refflat} \\
      STRAND=${strand} \\
      RIBOSOMAL_INTERVALS=${rrna}
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
  cpus 20
  memory '100 GB'
  time '24h'
  echo debug_mode

	input:
  //val x from featurecounts_ch.collect()
  val x from rnaseqmetrics_complete_ch.collect()
	val bams from bam_featurecounts_ch.collect()

  output:
	val "x" into featurecounts_complete_ch

	//when:
	//run_align

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



  if( run_featurecounts )
    """
      cd ${stardir}
      bamstring=\$(echo $bams | sed 's/,/ /g' | sed 's/\\[//g' | sed 's/\\]//g' )
      echo \$bamstring

      echo "gtf: ${gtf}"
      featureCounts -T ${task.cpus} \\
        -t ${params.fcounts_feature} \\
        --extraAttributes gene_name,gene_type \\
        -a ${gtf} -g gene_id  \\
        -o ${featurecountsdir}/${projectid}_geneid.featureCounts.txt \\
        -p \\
        -s ${strand_numeric} \${bamstring}

    """
  else
    """
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
  cpus 4
  memory '32 GB'
  time '3h'
  echo debug_mode
  //publishDir "${stardir}", mode: 'copy', overwrite: 'true'

  input:
  //val x from align_complete_ch()
  val x from featurecounts_complete_ch.collect()
  set sid, bam, strand, species from bam_indexbam_ch


  output:
  val "x" into indexbam_complete_ch

  // when:
  // run_align

  script:
  if ( run_bam_indexing )
    """
    cd ${stardir}
    echo "${stardir}/${bam}"
    samtools index -bc ${stardir}/${bam}
    """
  else
    """
    """
}



// picard mark duplicates

process markdups {
  tag "$id"
  cpus 4
  memory '32 GB'
  time '24h'
  echo debug_mode

  input:
  //val x from markdups_ch.collect()
  val x from indexbam_complete_ch.collect()
  set sid, bam, strand, species from bam_markdups_ch

  output:
  val "x" into markdups_complete_ch
  // val "x" into move

  //when:
  //run_align

  script:
  if ( run_markdups )
  """
  echo "bam: ${bam}"
  echo "markdupsdir: ${markdupsdir}/${bam}"
    java -jar /usr/local/bin/picard.jar MarkDuplicates \\
      INPUT=${stardir}/${bam} \\
      OUTPUT=${markdupsdir}/${bam} \\
      METRICS_FILE=${markdupsqcdir}/${sid}_bam.MarkDuplicates.metrics.txt \\
      TAGGING_POLICY=All \\
      REMOVE_DUPLICATES=false \\
      ASSUME_SORTED=true \\
      MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000 \\
      QUIET=true \\
      VERBOSITY=WARNING

    mv -f ${markdupsdir}/${bam} ${stardir}/${bam}
  """
  else
  """
  """
}



// Collect processes and prepare for MultiQC
process collect_align {
  tag "$id"
  cpus 1
  memory '8 GB'
  time '3h'
  // echo true

  input:
  val x from markdups_complete_ch.collect()

  output:
  val "x" into align_complete_ch


  script:
  if ( run_align)
  """
  """
  else
  """
  """

}

/* ===============================================================
  *      FASTQSCREEN
  =============================================================== */

// fastq_screen
process fastqScreen {
    tag "$id"
    cpus 16
    memory '32 GB'
    time '24h'

    input:
    set sid, read1, read2, species from fastqscreen_ch //

    output:
    val "x" into fastqscreen_complete_ch

    script:
    if ( params.paired ){
        fqsfiles = "${fastqdir}/${read1} ${fastqdir}/${read2}" }
    else{
        fqsfiles = "${fastqdir}/${read1}" }

    if ( run_fastqscreen)
    """
      /usr/local/bin/FastQ-Screen-0.14.1/fastq_screen \\
        --conf ${params.fastqscreen_config} \\
        --subset 500000 \\
        --outdir ${fastqscreendir} \\
        ${fqsfiles}
    """
    else
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
  cpus 6
  memory '32 GB'
  time '3h'
  // echo true

  input:
  set sid, read1, read2, species from fastqc_ch  // from check fastq

  output:
  val "x" into fastqc_complete_ch

  //input:
  //  set file(samples_csv) from sheet_ctg_ch
  script:
  if ( params.paired && run_fastqc)
    """
      echo "running fastqc in paired reads mode"
      fastqc ${fastqdir}/${read1} ${fastqdir}/${read2}  --outdir ${fastqcdir}
  """
  else if ( !params.paired && run_fastqc)
    """
      echo "running fastqc in non paired reads mode "
      fastqc ${fastqdir}/${read1}  --outdir ${fastqcdir}
    """
  else
    """
    """

}




/* ===============================================================
  *      CTG MULTIQC - MULTIQC ON ALL ANALYSES - runfolder included
  =============================================================== */
// This multiQC is for CTG infouse and not to the customer.
// The customer will obtain a lighter multiQC carried out below

// run the inhouse multiqc analysis on entire nf-output dir as well as the runfolderdir (to get interop stats)
process multiqc_ctg {
  //publishDir "${multiqcctgdir}", mode: 'copy', overwrite: 'true'
  tag "$id"
  cpus 6
  memory '32 GB'
  time '3h'
  echo debug_mode

  input:
  val x from fastqc_complete_ch.collect()
  val x from align_complete_ch.collect()
  val x from fastqscreen_complete_ch.collect()

  output:
  val "x" into multiqc_ctg_complete_ch

  when:
  run_multiqc_ctg

  script:
  """
    cd ${outputdir}
    multiqc -n ${projectid}_multiqc_report \\
      --interactive \\
      -o ${multiqcctgdir} . ${runfolderdir}
  """

}







// /* ===============================================================
//   *      ADD TO OUTBOX FOR CONVENIANT DOWNLOAD
//   =============================================================== */
// process add_outbox {
//   tag "$id"
//   cpus 6
//   memory '32 GB'
//   time '3h'
//   echo debug_mode
//
//   input:
//   val x from multiqc_ctg_complete_ch.collect()
//
//   output:
//   val "x" into add_outbox_complete_ch
//
//
//   script:
//
//   if ( params.copy_to_outbox ){
//     """
//     cp -r ${multiqcctgdir} ${outboxdir}
//     cp -r ${fastqcdir} ${outboxdir}
//
//     """}
//   else{
//     """
//     """}
//
// }
