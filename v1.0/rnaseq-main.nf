
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





/* ===============================================================
  *      PARAMS FROM CONFIGS
  =============================================================== */

//  project specific config
// -----------------------------

// root directories
project_root        =  params.project_root
delivery_root       =  params.delivery_root
ctg_save_root       =  params.ctg_save_root

//  project and run folders
projectid           =  params.projectid
projectdir          =  params.projectdir
bindir              =  params.bindir
// n_samples           =  params.n_samples


//  samplesheets
samplesheet_ctg       =  params.samplesheet           // simple sample sheet used for pipeline. Must includes file paths to fastq and bamsm, as well as species etc.
samplesheet_demux     =  params.samplesheet_demux     // IEM style samplesheet used for bcl2fastq. often generated by iem-samplesheet-processor.R
samplesheet_original  =  params.samplesheet_original // The original (non corrected) sample sheet obtained from lab


//  demux specific
runfolderdir        =  params.runfolderdir        // illumina raw data runfolder (full path)
runfolder           =  params.runfolder           // illumina raw data runfolder (folder name only)
fastqdir_bcl2fastq  =  params.fastqdir_bcl2fastq  // base directry where blc2fastq will write its output to (including undetermined fastq and stats folder).
fastqdir            =  params.fastqdir            // subdirectory where blc2fastq will write fastq files to. fastq-files will be read according to sample sheet. Defaults to <fastqdir_bcl2fastq>/<projectid>


// Delivery params
outboxsyncdir             = params.outboxsyncdir     // folder in ls4 private folder (e.g. /projects/fs1/percebe/outboxsync) used to sync with /box/outbox/percebe (as /box/.. is not mounted on lfs)
//deliver_raw         =  params.deliver_raw    // if to transfer raw data to delivery folder. Detaults to FALSE for RNAseq pipe
//deliver_fastq       =  params.deliver_fastq  // if to transfer fastq to output. Defaults to TRUE




/* ===============================================================
  *      DEFINE DIRECTORIES FROM PARAMS
  =============================================================== */

outputdir =  projectdir+'/nf-output' // main ooutput directory for files genetated with the Pipeline
file(outputdir).mkdir() // main nexttlow work dir for output of analyses. Subdirectory of the project foilder. Files and folders will be moved and copiued from this folder upon pipeline  completion.

// the deliverytemp will be used to save analyses that are bound for delivery t ocustomer
deliverytemp        =  outputdir+'/delivery' // this temp deliverydir is used within the nf workfolder/outputdir to store files that are comitted for delivery. A customer multiqc will be run only on this dir. Upon completion of all analyses this will be moved to delivery dir

featurecountsdir = outputdir+'/featurecounts'

stardir = deliverytemp+'/star'
fastqcdir = deliverytemp+'/fastqc'

markdupsdir = outputdir+'/markdups_bam_tmp'
markdupsqcdir = outputdir+'/markdups'
rnaseqmetricsdir = outputdir+'/rnaseqmetrics'
multiqcctgdir = outputdir+'/multiqc-ctg'
fastqscreendir = outputdir+'/fastqscreen'


// delivery
ctg_save_dir        =  ctg_save_root + '/' + projectid
multiqcdeliverydir  =  deliverytemp+'/multiqc'
deliverydir         =  delivery_root + '/' + projectid  // final delivery dir (on ... /nas-sync/. Note that delivery is prepared in "deliverytemp" is used in projectfolder)


// Illumina runfolder stats
interopdir_ilm = runfolderdir + '/InterOp'
interopdir_ctg = runfolderdir + '/ctg-interop'



/* ===============================================================
  *       create output and logdirs
  =============================================================== */




// if( params.run_align ) file(stardir).mkdir()
/// if( params.run_align ) file(markdupsdir).mkdir()
// if( params.run_align ) file(rnaseqmetricsdir).mkdir()
//if( params.run_align && params.run_featurecounts ) file(featurecountsdir).mkdir()
// if( params.run_fastqscreen ) file(fastqscreendir).mkdir()
// if( params.sync_outbox ) file(outboxsyncdir).mkdir()


// create project specific delivery dir and ctg qc dir
// -----------------------------
//file(deliverydir).mkdir()
//file(ctg_save_dir).mkdir()
// readme deliverydir
readme = deliverydir +'/README_ctg_delivery_' + projectid

// log file for nextflow .onComplete
logfile              =  file( projectdir + '/' + projectid + '.nextflow.log.complete' )
logfile_sav          =  file( ctg_save_dir + '/' + projectid + '.nextflow.log.complete' )


/* ===============================================================
  *       CHECKS FILES AND PARAMS
  =============================================================== */


//  Check paramters
// -----------------------------
if (projectid         == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet_ctg   == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}


// Check if files and directories exist
checkPathParamList = [
  project_root, delivery_root, ctg_save_root,
  projectdir, bindir,
  outputdir,
  samplesheet_ctg
]
for (param in checkPathParamList) {
    if (param) {
	file(param, checkIfExists: true)
    }
}

// Demux specific (bcl2fastq2 )
// -----------------------------
// Check if runfolder is defined. If not set demux to false and assume that a custom fastq dir is supplied
if ( params.run_demux == true ) {
  file(runfolderdir, checkIfExists: true)
  file(samplesheet_demux, checkIfExists: true)
}


// // Debug & test params
// // -----------------------------
debug_mode = false // will turn echo to true



/* ===============================================================
  *       MESSAGES
  =============================================================== */

// Define messages to print and for logfiles
def msg_startup = """\

    Workflow execution parameters
    ---------------------------------
    project id              :  ${projectid}
    project work dir        :  ${projectdir}
    nextflow execution dir  :  ${baseDir}
    nextflow output dir     :  ${outputdir}
    nextflow work dir       :  ${workDir}
    sample sheet ctg        :  ${samplesheet_ctg}
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

  if ( new File( logfile ).exists() && ! new File( logfile_sav ).exists())  { new File( logfile_sav ) << new File( logfile ).text }

  println( msg_completed )
}


def msg_modules = """\

    Run modules
    ---------------------------------
    demux       :  ${params.run_demux}
    fastqc      :  ${params.run_fastqc}

   """
   .stripIndent()

println( msg_modules )





// all samplesheet info
// if ( params.paired == true ) {
Channel
  .fromPath(samplesheet_ctg)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.fastq_1, row.fastq_2, row.Species ) }
  .tap{ infoall }
  .set { fastq_ch }

Channel
  .fromPath(samplesheet_ctg)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.bam, row.Strandness, row.Species) }
  .tap { infobam }
  .into { bam_checkbam_ch; bam_indexbam_ch; bam_rnaseqmetrics_ch; bam_markdups_ch }

Channel
    .fromPath(samplesheet_ctg)
    .splitCsv(header:true)
    .map { row -> tuple( row.bam) }
    .tap{ infoallfcounts }
    .set { bam_featurecounts_ch }


println " > Samples to process: "
println "[Sample_ID,fastq1,fastq2]"
infoall.subscribe { println "Info: $it" }
println " > Projects to process : "




/* ===============================================================
  *    ++++ START PROCESSES ++++
  =============================================================== */



/* ===============================================================
  *    --- Demux and fastq files section ---
  =============================================================== */

// Run bcl2fastq if run_demux
process bcl2fastq {
  // -w must be lower than number of samples
  // publishDir "${fastqdir}", mode: 'copy', overwrite: 'true'
  cpus 4
  tag "$id"
  memory '110 GB'
  time '24h'

  when:
  params.run_demux

  input:
  val samplesheet_demux

  output:
  val "x" into fastq_check_ch

  script:
  """
  mkdir -p ${fastqdir_bcl2fastq}

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
if ( params.run_demux == false ) {
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

  if( params.paired && params.run_checkfiles )
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
  else if ( params.paired == false && params.run_checkfiles)
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


  if ( params.run_align )
  """
  mkdir -p ${stardir}
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
  if( params.run_checkfiles )
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



  if ( params.run_rnaseqmetrics )
  """
    echo "strand: ${strand}"
    echo "rrna file: ${rrna}"
    echo "refflat file: ${refflat}"
    mkdir -p ${rnaseqmetricsdir}

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



  if( params.run_featurecounts )
    """
      mkdir -p ${featurecountsdir}
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
  if ( params.run_bam_indexing )
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
  if ( params.run_markdups )
  """
  echo "bam: ${bam}"
  mkdir -p ${markdupsdir}
  mkdir -p ${markdupsqcdir}

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
  if ( params.run_align)
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

    if ( params.run_fastqscreen)
    """
    mkdir -p ${fastqscreendir}

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
  if ( params.paired && params.run_fastqc)
    """
      mkdir -p ${fastqcdir}
      echo "running fastqc in paired reads mode"
      fastqc ${fastqdir}/${read1} ${fastqdir}/${read2}  --outdir ${fastqcdir}
  """
  else if ( !params.paired && params.run_fastqc)
    """
      mkdir -p ${fastqcdir}
      echo "running fastqc in non paired reads mode "
      fastqc ${fastqdir}/${read1}  --outdir ${fastqcdir}
    """
  else
    """
    """

}



/* ===============================================================
  *     ++++ POST RUN CRUNCHING SECTION ++++
  =============================================================== */



/* ===============================================================
  *     ctg multiqc - on all analyses - runfolder included
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
  val "x" into multiqc_ctg_complete_2_ch

  when:
  params.run_multiqc_ctg

  script:
  """
    mkdir -p ${multiqcctgdir}
    cd ${outputdir}
    multiqc -n ${projectid}_multiqc_report \\
      --interactive \\
      -o ${multiqcctgdir} . ${runfolderdir}
  """

}




/* ===============================================================
  *     Genaerate Delivery folder (temp folder within project dir)
  =============================================================== */
// generate a delivery folder and collect all files n folders to deliver
// run additional multiqc and md5 summ
// This temp delivery folder can be moved to delivery site on ls4 (nas sync) after nextflow sctipt completion.
// this delivery in additional shell script.
// move fastq files to delivery folder
// if not pooled, deliver the complete bcl2fastq directory including stats and undetermined fastq

// setup a deliveryfolder in nextflow dir.
// This will later be moved to /nas-sync/ctg-delivery

process setup_deliverytemp {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'

  input:
  val x from multiqc_ctg_complete_2_ch.collect()

  output:
  val "x" into setup_deliverytemp_complete_ch

  script:
  """
    mkdir -p ${deliverytemp}

    ## copy sample sheets to delivery
    cp ${samplesheet_ctg} ${deliverytemp}
    if [ -f ${samplesheet_demux} ]; then
     cp ${samplesheet_demux} ${deliverytemp}
    fi
    ## deliver star folder if exists
    if [ -d ${stardir} ]; then
      mv ${stardir} ${deliverytemp}
    fi
    ## deliver featurecounts folder if exists
    if [ -d ${featurecountsdir} ]; then
      mv ${featurecountsdir} ${deliverytemp}
    fi
    ## deliver fastqc folder if exists
    if [ -d ${fastqcdir} ]; then
     mv ${fastqcdir} ${deliverytemp}
    fi
  """
}

process move_fastq {

  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'
  echo debug_mode

  input:
  val x from setup_deliverytemp_complete_ch.collect()

  output:
  val "x" into move_fastq_complete_ch
  val "x" into move_fastq_complete_2_ch

  script:
    if ( params.pooled_run &&  params.deliver_fastq)
      """
        mkdir -p ${deliverytemp}/fastq
        if [ -d ${fastqdir} ]; then
          echo "pooled run. moving fastq foldler only."
          mv ${fastqdir} ${deliverytemp}/fastq
        fi
      """
    else if ( !params.pooled_run &&  params.deliver_fastq )
      """
      if [ -d ${fastqdir_bcl2fastq} ]; then
        echo "non pooled data. moving comlplete bcl2fastq output foldler."
        mv ${fastqdir_bcl2fastq} ${deliverytemp}
      elif [ -d ${fastqdir} ]; then
        echo "non pooled run but cannot locate fastqdir_bcl2fastq. moving fastq foldler only."
        mkdir -p ${deliverytemp}/fastq
        mv ${fastqdir} ${deliverytemp}/fastq
      fi
      """
    else
     """
     """
}



// Run customer multiQC (on delivery temp folder only).
// Not to be confused with ctg-multiqc that is run on all analyses and runfolder
process multiqc_delivery {

  tag "$id"
  cpus 6
  memory '32 GB'
  time '3h'
  echo debug_mode

  input:
  val x from move_fastq_complete_ch.collect()

  output:
  val "x" into multiqc_complete_ch

  script:
  mqcreport = multiqcdeliverydir + '/' + projectid + '_multiqc_report'

  if (! new File( mqcreport+'.html' ).exists() )
    """
      mkdir -p ${multiqcdeliverydir}
      cd ${deliverytemp}
      multiqc -n ${mqcreport} \\
        --interactive \\
        -o ${multiqcdeliverydir} .
    """
  else
  """
    echo "${mqcreport} already exists - skipping"
  """
}


// generate md5 sum on all files included in delivery temp folder
process md5sum_delivery {
  cpus 8
  tag "$id"
  memory '64 GB'
  time '3h'

  input:
  val x from multiqc_complete_ch.collect()

  output:
  val "x" into md5sum_complete_ch

  when:
  params.run_md5sum

  script:
  md5sumfile = deliverytemp + '/md5sum.txt'

  if (! new File( md5sumfile ).exists() )
    """
     cd ${deliverytemp}
     find . -type f -exec md5sum {} \\; > ${md5sumfile} ; echo
    """
  else
    """
      echo "${md5sumfile} already exists. skipping."
    """
}




/// provess add README with dir size
process finalize_delivery {
  cpus 2
  tag "$id"
  memory '16 GB'
  time '3h'

  input:
  val x from md5sum_complete_ch.collect()

  output:
  val "x" into finalize_delivery_ch

  script:

  """
  mv ${deliverytemp} ${deliverydir}
  cd ${deliverydir}
  echo "ctg delivery complete"               > $readme
  echo "Project:   ${projectid}"             >> $readme
  du -ch -d 0 . | grep 'total'               >> $readme
  """
}







/* ===============================================================
  *     FINALIZE CTG SAVE - save qc files and scripts
  =============================================================== */
// CTG should store multiQC (ctg-multiqc) and fastqc for all samples
// as of this version files are copied to ctg-qc dir. could change to the folder that also keep scripts and configs and sample sheets
// save

//  - multiQC (ctg_multiqc)
//  - fastQC
//
//  - Sample Sheets in ./samplesheets/
//  - Nextflow scripts, nextflow.params., rnaseq-main, nextflow.config, drivers, ./bin files etc (./scripts)

//  - logs,  (the final nextflow genereated onComplee is copied in that segion)

process setup_ctg_save {
  cpus 4
  tag "$id"
  memory '32 GB'
  time '3h'

  input:
  val x from multiqc_ctg_complete_ch.collect()
  val x from move_fastq_complete_2_ch.collect()


  output:
  val "x" into setup_ctg_save_complete_ch

  script:
  """
  mkdir -p ${ctg_save_dir}/samplesheets
  mkdir -p ${ctg_save_dir}/scripts

  cd ${projectdir}

  ## ctg specific multiQC and fastQCs
  if [ -d ${multiqcctgdir} ]; then
    cp -r ${multiqcctgdir} ${ctg_save_dir}
  fi
  if [ -d ${fastqcdir} ]; then
    cp -r ${fastqcdir} ${ctg_save_dir}
  fi

  ## all sample sheets
  if [ -f ${samplesheet_ctg} ]; then
    cp ${samplesheet_ctg} ${ctg_save_dir}/samplesheets
  fi
  if [[ -f ${samplesheet_demux} ]]; then
    cp ${samplesheet_demux} ${ctg_save_dir}/samplesheets
  fi
  if [ -f ${samplesheet_original} ]; then
    cp ${samplesheet} ${ctg_save_dir}/samplesheets
  fi

  ## rhe samplesheet check rscript output
  if [[ -f "${runfolderdir}/iem.rscript.log" ]]; then
    cp ${runfolderdir}/iem.rscript.log ${ctg_save_dir}
  fi

  ## project specic parameters file
  if [[ -f "${projectdir}/nextflow.params.${projectid}" ]]; then
    cp ${projectdir}/nextflow.params.${projectid} ${ctg_save_dir}/scripts
  fi

  ## copy the entire scripts dir into the ctg save dir
  if [[ -d "${params.scriptsdir}/rnaseq-main.nf" ]]; then
    cp -r ${params.scriptsdir} ${ctg_save_dir}/scripts
  fi

  """
}
