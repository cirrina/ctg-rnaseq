
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
n_samples           =  params.n_samples
projectdir          =  params.projectdir
//execdir             =  params.execdir
bindir              =  params.bindir


//  samplesheets
samplesheet         =  params.samplesheet
samplesheet_demux   =  params.samplesheet_demux

//  assay specific
assay               =  params.assay
instrument_type     =  params.instrument_type
index_adapters      =  params.index_adapters
paired              =  params.paired
strandness          =  params.strandness

//  demux specific
runfolderdir        =  params.runfolderdir
runfolder           =  params.runfolder
fastqdir            =  params.fastqdir
pooled              =  params.pooled

//  module specific
run_demux             =  params.run_demux
run_fastqc            =  params.run_fastqc
run_multiqc           =  params.run_multiqc
run_multiqc_ctg       =  params.run_multiqc_ctg
run_fastq_screen      =  params.run_fastq_screen
bam_indexing          =  params.bam_indexing
run_markdups      =  params.run_markdups
run_rnaseqmetrics  =  params.run_rnaseqmetrics

//  log files
logdir               =  '/Users/david/tasks/rnaseq_test/ctg-projects/2021_030/'

//  Other Parameters
// -----------------------------
outputdir =  projectdir+'/nf-output' // main ooutput directory for files genetated with the Pipeline
deliverydir = delivery_root + '/' + projectid

ctg_qc_dir      = qc_root+'/ctg-rnaseq'

// featurecounts
fcounts_feature     =  params.fcounts_feature


//  Check paramters
// -----------------------------
if (projectid      == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet    == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}




// create output and logdirs
file(outputdir).mkdir()
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

    delivery dir            :  ${deliverydir}
    ctg qc dir              :  ${ctg_qc_dir}

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
  	nf workDir   : ${workflow.workDir}
    nf baseDir   : ${workflow.baseDir}
  	exit status  : ${workflow.exitStatus}
  	errorMessage : ${workflow.errorMessage}
  	errorReport  :

  	"""
  	.stripIndent()

  def error = """\
		${workflow.errorReport}
	   """


  // base = csv.getBaseName()
  logFile = file( logdir + '/' + projectid + '.log.complete' )
  logFile.text = msg_startup.stripIndent()
  logFile.append( msg_completed.stripIndent() )
  logFile.append( error )

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
if ( paired == true ) {
  Channel
    .fromPath(samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.fastq_1, row.fastq_2) }
    .tap{ infoall }
    .set { fastq_ch }
  println("running paired")
  }
else {
  Channel
      .fromPath(samplesheet)
      .splitCsv(header:true)
      .map { row -> tuple( row.Sample_ID, row.fastq_1) }
      .tap{ infoall }
      .set { fastq_ch }
  }

/// Channels input files


println " > Samples to process: "
println "[Sample_ID,fastq1,fastq2]"
infoall.subscribe { println "Info: $it" }
println " > Projects to process : "
//println "[agg,Sample_Project]"
//infoProject.subscribe { println "Info Projects: $it" }





// Run bcl2fastq
process bcl2fastq {
  // -w must be lower than number of samples
  cpus 4
  tag "$id"
  memory '8 GB'
  time '3h'



  when:
    run_demux
  input:
    val samplesheet_demux
  output:
    val "y" into demux_complete_ch
    val "x" into moveFastq

  script:

    """

    bcl2fastq -R $runfolderdir \\
              --sample-sheet $samplesheet_demux \\
              --no-lane-splitting  \\
              -r 1 \\
              -p $task.cpus  \\
              -w 1  \\
              --output-dir $fastqdir

     """
}


// process moveFastq {
//
//     input:
//     val x from moveFastq
//     val projid from mvFastq_ch
//
//     output:
//     val projid into fqc_ch, star_go
//
//     """
//     ## mkdir -p ${OUTDIR}/${projid}
//     ## mkdir -p ${OUTDIR}/${projid}/Fastq_Raw
//     ## mv ${FQDIR}/${projid}/* ${OUTDIR}/${projid}/Fastq_Raw/
//     """
//
// }



// Here insert process to check the output from bcl2fastq2
// process post_demux {
//     input: val x from demux_complete_ch
//
//     """
//     """
//
// }


// Channel to start count if demux == 'n'
// Projects
if ( run_demux == false ) {

   Channel
	 .from("y")
    	 .set{ demux_complete_ch }

}



// sample info from the ctg sample sheet
// samplesheet_ctg = outputdir+'/'+samplesheet_ctg
// Channel
//     .fromPath(samplesheet_ctg)
//     .splitCsv(header:true)
//     .map { row -> tuple( row.Sample_ID ) }
//     .unique()
//     .tap{ infoSamples }
//     .into{ fastqc_ch; star_ch }
//
// infoSamples.subscribe{ println "Samples: $it" }

// sheet_ctg_ch
//   .println { "File: ${it.name} => ${it.text}" }
//   .splitCsv


process fastqc {
  publishDir "${fastqcDir}", mode: 'copy', overwrite: 'true'
  tag "$id"
  cpus 1
  memory '5 GB'
  time '3h'

  input:
    val "y" from demux_complete_ch.collect()

  when:
    run_fastqc
  //input:
  //  set file(samples_csv) from sheet_ctg_ch
  script:
  println("Running fastqc")

    """

    """

}
