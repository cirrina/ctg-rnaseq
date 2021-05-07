
// Configuration implicit variables. These are implicitly defined in the Nextflow configuration file
  // baseDir    - projectDir in v20.04.x: The directory where the main workflow script is located
  // workDir    - the work folder for nextflow temporary work-files
  // launchDir  - the directory where the workflow is run (requires version 20.04.0 or later).


//  project specific config
// -----------------------------

// root directories
project_root        =  params.project_root
delivery_root       =  params.delivery_root
qc_root             =  params.qc_root
log_root            =  params.log_root

//  project  and run folders
projectid           =  params.projectid
projectdir          =  params.projectdir
execdir             =  params.execdir
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
runfolder           =  params.runFolder
fastqdir            =  params.fastqdir
pooled              =  params.pooled

//  module specific
run_demux             =  params.run_demux
run_fastqc            =  params.run_fastqc
run_multiqc           =  params.run_multiqc
run_multiqc_ctg       =  params.run_multiqc_ctg
run_fastq_screen      =  params.run_fastq_screen
bam_indexing          =  params.bam_index
piccard_markdups      =  params.run_markdups
picard_rnaseqmetrics  =  params.run_rnaseqmetrics

//  log files
logdir               =  '/Users/david/tasks/rnaseq_test/ctg-projects/2021_030'




//  Other Parameters
// -----------------------------



//  Check paramters
// -----------------------------
if (projectid      == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet    == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}


// Auto genaerate dir-parameters. Note that baseDir and workDir are reserved by nextflow and not used here.
// -----------------------------
//basedir            =  '/Users/david/tasks/rnaseq/nftest' // determines base ctg workfolder. previusly basedir
projectdir         =  params.projectdir// prev. set by basedir+'/'+projectid  // defined by project id set in config
nfworkdir          =  workDir // mind the difference between workdir and workDir that is reserved for nextflow
outputdir          =  projectdir+'/nf-output' // main ooutput directory for files genetated with the Pipeline
logdir             = projectdir+'/logfiles'

// create output and logdirs
file(outputdir).mkdir()
file(logdir).mkdir()


// Check if these exist
checkPathParamList = [
    projectdir,
    outputdir,
    bindir,
    samplesheet
]
for (param in checkPathParamList) {
    if (param) {
	file(param, checkIfExists: true)
    }
}


// Demux specific (bcl2fastq2 )
// -----------------------------
run_demux          =  params.demux        // Flag true or false if to run demux.
runfolderdir          =  params.runfolderdir    // Illumina runfolder. full path.
fastqdir           =  params.fastqdir     // Output directory for the bcl2fastq2 process. Defaults in driver to projectdir+fastq
samplesheet_demux     =  params.samplesheet_demux

// Check if runfolder is defined. If not set demux to false and assume that a custom fastq dir is supplied
if ( run_demux == 'true' ){
  checkPathParamList = [
      runfolderdir
      // samplesheet_demux
  ]
  for (param in [checkPathParamList]) {
    if (param) {
	     file(param, checkIfExists: true)
	    }
    }
}

// NOT YET IMPLEMENTED
// deliverydir        =  params.delivery_basedir+'/'+projectid // project specific sub-dir within the base dir above.
// ctgqcdir           =  params.ctg_qc_basedir+'/'+projectid // project specifc dir wihthin qc basedir above


// set names for sample sheets
//samplesheet_ctg      =  "SampleSheet-${projectid}-ctg.csv"
// samplesheet_demux    =  "SampleSheet-${projectid}-demux.csv"


// featurecounts
fcounts_feature     =  params.fcounts_feature



// Define parameters from the profiles section
// --------------------------------------------
species             =  params.species
gtf                 =  params.gtf
star_genome         =  params.star_genome
run_fastqc          =  params.fastqc
run_multiqc         =  params.multiqc
run_fastqc_ctg      =  params.multiqc_ctg

// ADDONS
run_fastqscreen         =  params.fastq_screen
run_bamindex            =  params.bam_index
run_markdups            =  params.piccard_markdups
run_rnaseqmetrics       =  params.picard_rnaseqmetrics



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





// Check if defined project work dirextory (workdir) from nextflow implicit baseDir (exection dir). If so exit on error
//    i.e. project must be executed in the folder used as project specific workfolder.
//
// if ( ! workdir  == baseDir ) {
//     exit 1, "You must exexute the workflow from wiithin a folder with the same name as the project_id parameter in nextflow.comfig, i.e. workdir (basedir+projectid) must be same as nextflow implicit baseDir"
// }


//



// Check if custom Fastq directory exists. If so run pipe using fastqs in this directory and set demux mode to false
// fastqcustom  =  params.fastq_custom_dir
// for (param in [fastqcustom]) {
//     if (param) {
// 	file(param, checkIfExists: true)
// 	run_demux = false
// 	println("")
// 	println("custom fastq dir supplied, setting demux to false")
// 	fastqdir =  fastqcustom  // set fastqdir to fastcustom
//     }
//     if (param == '' && runfolder == '') {
// 	exit 1, "You have to specify either illumina runfolder or fastq cutom dir in the nextflow.config"
//     }
// }
// if (run_demux) fastqdir     =  workdir+'/fastq'  // FastQ directory will be set into working directory
// //println("Fastq directory: "+params.fastqDirCustom)




// SET MODULES :: THIS IS TEMP FOR TESTING
// run_fastqc = true
// run_demux = true
// run_sheet_check = true


// Define messages to print and for logfiles
def msg_startup = """\

    Workflow execution parameters
    ---------------------------------
    project id              :  ${projectid}
    project work dir        :  ${projectdir}
    species                 :  ${species}

    nextflow execution dir  :  ${baseDir}
    nextflow output dir     :  ${outputdir}
    nextflow work dir       :  ${nfworkdir}

    delivery dir            :
    ctg qc dir              :

    sample sheet ctg        :  ${samplesheet}
    sample sheet demux      :  ${samplesheet_demux}
    custom_samplesheet     :  ${custom_samplesheet}

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
  	NF workDir   : ${workflow.workDir}
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
if ( paired == 'false' ) {
  Channel
    .fromPath(samplesheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.Sample_ID, row.fastq_1_path, row.fastq_2_path) }
    .tap{infoall}
    .set { fastqc_ch }
  }
else {
  Channel
      .fromPath(samplesheet)
      .splitCsv(header:true)
      .map { row -> tuple( row.Sample_ID, row.fastq_1, row.fastq_1_path) }
      .tap{infoall}
      .set { fastqc_ch }
  }

/// Channels input files


// =====================
//  CREATE DIRECTORIES
// =====================
// process init_runfolders {
//   """
//   mkdir -p ${outputdir}
//   """
//
// }

// myDir = file(outputdir)
// myDir.mkdir()


// =====================
//  CHECK SAMPLE SHEET
// =====================
// R script. Check sample sheet. Generate two new sample sheets used for demux and downstream analyses, respectively
// process check_sheet {
//   cpus 1
//   publishDir "${outputdir}", mode: 'copy', overwrite: 'true', pattern: '*.csv'
//   tag "$id"
//   memory '5 GB'
//   time '3h'
//
//   when:
//     run_sheet_check
//   output:
//     file "${samplesheet_ctg}" into sheet_ctg_ch
//     file "${samplesheet_demux}" into sheet_demux_ch
//
//     """
//     ${projectdir}/bin/ctg_sampleSheetCheck.R \\
//         --sample_sheet ${samplesheet} \\
//         --output_ctg_sheet ${samplesheet_ctg} \\
//         --output_demux_sheet ${samplesheet_demux}  \\
//         --paired ${params.paired} \\
//         --strandness ${params.strandness}
//     """
// }
//




process demux_bcl2fastq2 {
  cpus 1
  tag "$id"
  memory '5 GB'
  time '3h'

  when:
    run_demux
  input:
    val samplesheet_demux
  output:
    val "y" into demux_complete_ch

  script:
    println("Running demux")

  """

  """
}




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
if ( run_demux == 'false' ) {

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
    val y from demux_complete_ch.collect()

  when:
    run_fastqc
  //input:
  //  set file(samples_csv) from sheet_ctg_ch


    """

    """

}
