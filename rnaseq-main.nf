
]
/* ===============================================================
  *      PARAMS FROM CONFIGS
  =============================================================== */

//  re-assign some params from nextflow.configs
// ----------------------------------------------------
//  project and run folders
projectid           =  params.projectid    // ctg project id. e.g. 2021_024
project_dir         =  params.project_dir   // .../shared/ctg-rnaseq/uroscan/2021_024 // NOT to be confused with project_dir
delivery_dir        =  params.delivery_dir
samplesheet         =  params.samplesheet           // name of simple sample sheet used for pipeline. Must includes file paths to fastq and bamsm, as well as species etc.
fastq_dir           =  params.fastq_dir            // subdirectory fastq files are located. For a default run the output location from blc2fastq. fastq-files will be read according to sample sheet. Defaults to <bcl2fastq_dir>/<projectid>
ctg_qc_dir          =  params.ctg_qc_dir

/* ===============================================================
  *      DEFINE DIRECTORIES FROM PARAMS
  =============================================================== */
output_dir =  project_dir+'/nf-output' // nextflow (temp) output directory for files genetated with the Pipeline that are NOT to be delivered
file(output_dir).mkdir() // main nexttlow work dir for output of analyses. Subdirectory of the project foilder. Files and folders will be moved and copiued from this folder upon pipeline  completion.

deliverysamplesheets = delivery_dir+'/samplesheet'
deliveryscripts = delivery_dir+'/scripts'
deliveryconfigs = delivery_dir+'/configs'
deliverylogs = delivery_dir+'/logs'




// the deliverytemp will be used to save analyses that are bound for delivery t ocustomer

// output dirs to delivery
deliverytemp  =  output_dir+'/delivery' // this temp delivery_dir is used within the nf workfolder/output_dir to store files that are comitted for delivery. A customer multiqc will be run only on this dir. Upon completion of all analyses this will be moved to delivery dir

stardir = deliverytemp+'/star'
stardir_filtered = deliverytemp+'/star_filtered'
salmondir = deliverytemp+'/salmon'
rsemdir = deliverytemp+'/rsem'
bladderreportdir = deliverytemp+'/bladderreport'
featurecountsdir = deliverytemp+'/featurecounts'

deliveryqc = deliverytemp+'/qc'
fastqcdir = deliverytemp+'/qc/fastqc'
multiqcdelivery_dir  =  deliverytemp+'/qc/multiqc'
mqcreport = deliverytemp+'/qc/multiqc' + '/' + projectid + '_multiqc_report'
readme = delivery_dir +'/README_ctg_delivery_' + projectid


// output for qc - temp workdirs
qcdir = output_dir+'/qc'
qualimapdir = qcdir+'/qualimap'
rseqcdir = qcdir+'/rseqc'


markdupstempdir = qcdir+'/markdups_bam_tmp'
markdupsqcdir = qcdir+'/markdups'
rnaseqmetricsdir = qcdir+'/rnaseqmetrics'
multiqcctgdir = qcdir+'/multiqc-ctg'
fastqscreendir = qcdir+'/fastqscreen'


/// ctg sav dirs
ctg_save_samplesheets = ctg_save_dir+'/samplesheets'
ctg_save_scripts = ctg_save_dir+'/scripts'
ctg_save_configs = ctg_save_dir+'/configs'
ctg_save_logs =  ctg_save_dir+'/logs'


// Illumina runfolder stats
interopdir_ilm = runfolderdir + '/InterOp'
interopdir_ctg = runfolderdir + '/ctg-interop'






/* ===============================================================
  *       create output and logdirs
  =============================================================== */

// log file for nextflow .onComplete
logfile   =  file( project_dir + '/' + 'log.nextflow.complete' )
// logfile_sav          =  file( ctg_save_dir + '/' + 'log.nextflow.complete' )






/* ===============================================================
  *       CHECKS FILES AND PARAMS
  =============================================================== */


//  Check paramters
// -----------------------------
if (projectid         == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet   == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}


// Check if files and directories exist
checkPathParamList = [
  params.project_root, params.delivery_root, ctg_qc_root,
  project_dir,
  output_dir,
  samplesheet
]
for (param in checkPathParamList) {
    if (param) {
	file(param, checkIfExists: true)
    }
}



// // Debug & test params
// // -----------------------------
debug_mode = false // will turn echo to true



/* ===============================================================
  *       MESSAGES
  =============================================================== */



def msg_deliverymail = """\


 """


// Define messages to print and for logfiles
def msg_startup = """\

    Workflow execution parameters
    ---------------------------------
    project id              :  ${projectid}
    project work dir        :  ${project_dir}
    nextflow execution dir  :  ${baseDir}
    nextflow output dir     :  ${output_dir}
    nextflow work dir       :  ${workDir}
    samplesheet             :  ${samplesheet}

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

  // if ( new File( logfile ).exists() && ! new File( logfile_sav ).exists())  { new File( logfile_sav ) << new File( logfile ).text }

  println( msg_completed )
}


def msg_modules = """\

    Run modules
    ---------------------------------
    run_blcl2fastq    :  ${params.run_blcl2fastq}
    fastqc        :  ${params.run_fastqc}

   """
   .stripIndent()

println( msg_modules )



/* ===============================================================
  *    PROCESS SAMPLE SHEET & DEFINE CHANNELS
  =============================================================== */

// Process ctg IEM style sample sheet & save
// samplesheet_nextflow="${project_dir}/SampleSheet-nexflow.csv"

// Read and process sample sheet. Save to SampleSheet-nexflow.csv
// samplesheet to be parsed as input channel (take everything below [Data] section).
sheet = file(params.samplesheet)
all_lines = sheet.readLines()
write_row = false // if next lines has sample info
sheet_nf = file("${project_dir}/SampleSheet-nexflow.csv")
sheet_nf.text=""

for ( line in all_lines ) {
  if ( write_row ) {
    sheet_nf.append(line + "\n")
  }
  if (line.contains("[Data]")) {
    write_row = true
  }
}

// Define Channels based from SampleSheet
Channel
  .fromPath(chsheet)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.fastq_1, row.fastq_2, row.Species ) }
  .tap{ infoall }
  .set { fastq_ch }

Channel
  .fromPath(chsheet)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.bam, row.Strandness, row.Species, row.RIN, row.concentration ) }
  .tap { infobam }
  .into { bam_checkbam_ch; bam_indexbam_ch; bam_rnaseqmetrics_ch; bam_markdups_ch; bam_filter_multimap_ch; bam_qualimap_ch; bam_rseqc_ch;  bam_bladderreport_ch }

Channel
    .fromPath(chsheet)
    .splitCsv(header:true)
    .map { row -> tuple( row.bam ) }
    .tap{ infoallfcounts }
    .set { bam_featurecounts_ch }

println " > Samples to process: "
println "[Sample_ID,fastq1,fastq2,species]"
infoall.subscribe { println "Info: $it" }




/* ===============================================================
  *    --- CHECK FASTQ FILES ---
  =============================================================== */

process checkfiles_fastq {
  // Run fastqc. Also check if all expected files, defined in the ctg samplesheet, are present in fastq_dir

  tag  "$sid"
  cpus params.cpu_min
  memory params.mem_min

  input:
  set sid, read1, read2, species from fastq_ch

  output:
  val "x" into checkfiles_fastq_complete_ch
  set sid, read1, read2, species into fastqc_ch

  script:
  if( params.paired )
    """
      echo "checking fastq files - paired reads "
      if [ ! -f ${fastq_dir}/${read1} ]; then
        echo "Warning: Cannot locate fastq_1 file ${fastq_dir}/${read1}"
        exit 2
      fi

      if [ ! -f ${fastq_dir}/${read1} ]; then
        echo "Warning: Cannot locate fastq_2 file ${fastq_dir}/${read2}"
        exit 2
      fi
    """

  else
    """
      echo "checking fastq files - non paired  "

      if [ ! -f ${fastq_dir}/${read1} ]; then
        "Cannot locate fastq_1 file ${read2}"
        exit 2
      fi
    """

}


/* ===============================================================
  *      FASTQC
  =============================================================== */

process fastqc {
  // Run fastqc. Also check if all expected files, defined in the ctg samplesheet, are present in fastq_dir
  tag  { params.run_fastqc  ? "$sid" : "blank_run"  }
  cpus { params.run_fastqc  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_fastqc  ? params.mem_standard : params.mem_min  }


  input:
  val x from checkfiles_fastq_complete_ch.collect()
  set sid, read1, read2, species from fastqc_ch  // from check fastq

  output:
  val "x" into fastqc_complete_ch
  val "x" into run_star_ch
  val "x" into run_salmon_ch
  val "x" into run_rsem_ch
  set sid, read1, read2, species into star_ch
  set sid, read1, read2, species into salmon_ch
  set sid, read1, read2, species into rsem_ch
  set sid, read1, read2, species into fastqscreen_ch


  // when: params.run_fastqc

  script:
  if ( params.paired && params.run_fastqc)
    """
      mkdir -p ${fastqcdir}
      echo "running fastqc in paired reads mode"
      fastqc ${fastq_dir}/${read1} ${fastq_dir}/${read2}  --outdir ${fastqcdir}

      ## find ${fastqcdir} -user $USER -exec chmod g+rw {} +

  """
  else if ( !params.paired && params.run_fastqc)
    """
      mkdir -p ${fastqcdir}
      echo "running fastqc in non paired reads mode "
      fastqc ${fastq_dir}/${read1}  --outdir ${fastqcdir}

      ## find ${fastqcdir} -user $USER -exec chmod g+rw {} +
    """
  else
    """
    echo "run_fastqc skipped"
    """
}



/* ===============================================================
  *      FASTQSCREEN
  =============================================================== */

process fastqscreen {

    tag  { params.run_fastqscreen  ? "$sid" : "blank_run"  }
    cpus { params.run_fastqscreen  ? params.cpu_standard : params.cpu_min  }
    memory { params.run_fastqscreen  ?  params.mem_standard : params.mem_min  }


    input:
    set sid, read1, read2, species from fastqscreen_ch //

    output:
    val "x" into fastqscreen_complete_ch

    script:
    if ( params.paired ){
        fqsfiles = "${fastq_dir}/${read1} ${fastq_dir}/${read2}" }
    else{
        fqsfiles = "${fastq_dir}/${read1}" }

    if ( params.run_fastqscreen)
      """
      mkdir -p ${fastqscreendir}

      fastq_screen \\
          --conf ${params.fastqscreen_config} \\
          --subset 500000 \\
          --outdir ${fastqscreendir} \\
          ${fqsfiles}


      """
    else
      """
      echo "run_fastqscreen skipped"
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


//
// if ( run_star == false ) {
//   Channel
// 	 .from("x")
//   .into{ align_complete_ch ; align_complete_report_ch}
// }



// Run salmon
// if ( params.run_salmon == false ) {
//    Channel
// 	 .from("x")
//    .set{ salmon_complete_ch }
// }
process salmon  {
  tag  { params.run_salmon  ? "$sid" : "blank_run"  }
  cpus { params.run_salmon  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_salmon  ?  params.mem_standard : params.mem_min  }

  // cpus 6
  // memory '48 GB'
  // time '36h'
  // echo debug_mode
  // //publishDir "${stardir}", mode: 'copy', overwrite: true

  input:
  val x from run_salmon_ch.collect()
  set sid, read1, read2, species from salmon_ch // from checkfiles_fastq

  output:
  val "x" into salmon_complete_ch

  // when: params.run_salmon

  script:
  if ( species == "Homo sapiens" ){
    transcripts=params.salmon_transcripts_hs }
  else if ( species == "Mus musculus" ){
    transcripts=params.salmon_transcripts_mm  }
  else if ( species == "Rattus norvegicus" ){
    transcripts=params.salmon_transcripts_rn  }
  else{
    genome = ""
    println( "Warning: Species not recognized." )}

  // if ( params.paired )


  if ( params.paired && params.run_salmon )
    """
    salmon quant -l A \\
      -i  ${transcripts} \\
      -1  ${fastq_dir}/${read1} \\
      -2  ${fastq_dir}/${read2} \\
      -p  6 --validateMappings \\
      -o  ${salmondir}/${sid}_0.salmon.salmon \\
      --no-version-check

    ## find ${salmondir} -user $USER -exec chmod g+rw {} +

    """
  else if ( !params.paired && params.run_salmon )
    """
    salmon quant -l A \\
      -i  ${transcripts} \\
      -1  ${fastq_dir}/${read1} \\
      -p  6 --validateMappings \\
      -o  ${salmondir}/${sid}_0.salmon.salmon \\
      --no-version-check

    #find ${salmondir} -user $USER -exec chmod g+rw {} +
    """
  else
    """
    echo "skipping salmon"
    """

}


/* ===============================================================
  *      -RSEM SECTION
  =============================================================== */
process rsem {
  tag  { params.run_rsem  ? "$sid" : "blank_run"  }
  cpus { params.run_rsem  ? params.cpu_high : params.cpu_min  }
  memory { params.run_rsem  ?  params.mem_high : params.mem_min  }

  // cpus 20
  // memory '100 GB'
  // time '36h'
  // //publishDir "${stardir}", mode: 'copy', overwrite: true

  input:
  val x from run_rsem_ch.collect()
  set sid, read1, read2, species from rsem_ch // from checkfiles_fastq

  output:
  val "x" into rsem_complete_ch
  val "x" into rsem_complete_report_ch
  // file "${sid}_Aligned.sortedByCoord.out.bam" into bam_featurecounts_ch // channel defined start instead

  script:

  // species and references (bowtie2 refs)
  if ( species == "Homo sapiens" ){
    genome=params.rsem_bowtie2_genome_hs }
  else if ( species == "Mus musculus" ){
    genome=params.star_genome_mm }
  else if ( species == "Rattus norvegicus" ){
      genome=params.star_genome_rn }
  else{
    genome = ""
    println( "Warning: Species not recognized." )
  }

  // paired end
  if ( params.paired ){
    rsemfiles = "${fastq_dir}/${read1} ${fastq_dir}/${read2}"
    paired='--paired-end'}
  else{
    rsemfiles = "${fastq_dir}/${read1}"
    paired=''}

  // strand
  if( params.strandness == "forward" )
    strand = 'forward'
  else if ( params.strandness == "reverse" )
    strand = 'reverse'
  else
    strand = 'none'

    //the uroscan pipe is run without strandness flag.
  if ( params.run_rsem && params.pipelineProfile == "uroscan" )
    """
    mkdir -p ${rsemdir}
    rsem-calculate-expression \\
        --num-threads ${task.cpus} \\
        --paired-end \\
        --bowtie2 \\
        --bowtie2-path /opt/software/uroscan_env/bin \\
        --estimate-rspd \\
        --append-names \\
        --no-bam-output \\
        ${rsemfiles} \\
        ${genome} \\
        ${rsemdir}/${sid}.rsem
    """
  else
    """
    echo "rsem not run"
    """
}




/* ===============================================================
  *      -STAR AND BAM SECTION
  =============================================================== */

process star  {
  tag  { params.run_star  ? "$sid" : "blank_run"  }
  cpus { params.run_star  ? params.cpu_high : params.cpu_min  }
  memory { params.run_star  ?  params.mem_high : params.mem_min  }

  // cpus 20
  // memory '100 GB'
  // time '36h'
  // echo debug_mode
  //publishDir "${stardir}", mode: 'copy', overwrite: true

  input:
  val x from run_star_ch.collect()
  set sid, read1, read2, species from star_ch // from checkfiles_fastq

  output:
  val "x" into star_complete_ch

  // when: params.run_star

  script:
  if ( species == "Homo sapiens" ){
    genome=params.star_genome_hs }
  else if ( species == "Mus musculus" ){
    genome=params.star_genome_mm }
  else if ( species == "Rattus norvegicus" ){
      genome=params.star_genome_rn }
  else{
    genome = ""
    println( "Warning: Species not recognized." )}

  if ( params.paired ){
      starfiles = "${fastq_dir}/${read1} ${fastq_dir}/${read2}" }
  else{
      starfiles = "${fastq_dir}/${read1}" }


  if ( params.run_star )
  """
  mkdir -p ${stardir}

  ### added genomeLoad remove - star crashes if not for version 2.5x uroscan pipeline
  ## STAR  --genomeDir ${genome} --genomeLoad Remove

  STAR --genomeDir ${genome} \\
    --readFilesIn ${starfiles} \\
    --runThreadN ${task.cpus}  \\
    --readFilesCommand zcat \\
    --outSAMtype BAM SortedByCoordinate \\
    --limitBAMsortRAM 10000000000 \\
    --outFileNamePrefix ${stardir}/${sid}_

  #find ${stardir} -user $USER -exec chmod g+rw {} +
  """
  else
  """
  echo "skipping star"
  """

}



process checkfiles_bam {
  // Run fastqc. Also check if all expected files, defined in the ctg samplesheet, are present in fastq_dir
  tag  { params.run_checkfiles_bam  ? "$sid" : "blank_run"  }
  cpus params.cpu_min
  memory params.mem_min

  input:
  val x from star_complete_ch.collect() // checkbam_ch - when star is completed
  set sid, bam, strand, species, RIN, concentration from bam_checkbam_ch

  output:
  val "x" into checkfiles_bam_complete_ch

  // when: params.run_checkfiles_bam

  script:
  if( params.run_checkfiles_bam )
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



/// INDEX BAMs


// samtools index bamfile
// ml Java; ml nextflow/19.04.1
// ml Singularity
// ml GCC/7.3.0-2.30
// ml SAMtools/1.9
// samtools index bamfile



process index_bam {
  tag  { params.run_index_bam  ? "$sid" : "blank_run"  }
  cpus { params.run_index_bam  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_index_bam  ?  params.mem_standard : params.mem_min  }



  input:
  val x from checkfiles_bam_complete_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_indexbam_ch

  output:
  val "x" into indexbam_complete_ch

  // when: params.run_index_bam

  script:
  if ( params.run_index_bam )
    """
    cd ${stardir}
    echo "${stardir}/${bam}"
    samtools index -bc ${stardir}/${bam}
    """
   // else if ( params.run_index_bam  && params.pipelineProfile == "uroscan"  )
   //  """
   //  cd ${stardir}
   //  echo "${stardir}/${bam}"
   //  sambamba index ${stardir}/${bam}
   //  """
  else
    """
    echo "skipped indexing"
    """
}



process markdups {
  tag  { params.run_markdups  ? "$sid" : "blank_run"  }
  cpus { params.run_markdups  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_markdups  ?  params.mem_standard : params.mem_min  }

  input:
  val x from indexbam_complete_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_markdups_ch

  output:
  val "x" into markdups_complete_ch
  val "x" into markdups_complete_fcounts_ch
  val "x" into markdups_complete_report_ch
  // val "x" into move

  // when: params.run_markdups

  script:
  if ( params.run_markdups )
    """
    echo "bam: ${bam}"
    mkdir -p ${markdupstempdir}
    mkdir -p ${markdupsqcdir}

    echo "markdupstempdir: ${markdupstempdir}/${bam}"
    # java -jar picard MarkDuplicates \\
    picard MarkDuplicates \\
        INPUT=${stardir}/${bam} \\
        OUTPUT=${markdupstempdir}/${bam} \\
        METRICS_FILE=${markdupsqcdir}/${sid}_bam.MarkDuplicates.metrics.txt \\
        TAGGING_POLICY=All \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=2000 \\
        QUIET=true \\
        VERBOSITY=WARNING

    mv -f ${markdupstempdir}/${bam} ${stardir}/${bam}

    ## find ${stardir} -user $USER -exec chmod g+rw {} +
    ## find ${markdupstempdir} -user $USER -exec chmod g+rw {} +
    """
  else
    """
    echo "run markdups skipped"
    """
}






/* ===============================================================
  *     FEATURE COUNTS SECTION
  =============================================================== */
// Filter bams on only primary mapped sequence using 0x104 flag.
// Featurecounts will only use primary mapped read anlyway - but tends to crash if multi mapped sequences are included

process filter_multimap {
  tag  { params.run_featurecounts  ? "$sid" : "blank_run"  }
  cpus { params.run_featurecounts  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_featurecounts  ?  params.mem_standard : params.mem_min  }

  input:
  val x from markdups_complete_fcounts_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_filter_multimap_ch

  output:
  val "x" into filter_multimap_complete_ch

  // when: params.run_markdups

  script:
  if ( params.run_featurecounts )
    """
    echo "bam: ${bam}"
    mkdir -p ${stardir_filtered}

    cd ${stardir_filtered}
    samtools  view -b -F 0x104  ${stardir}/${bam} >  ${stardir_filtered}/${bam}

    """
  else
    """
    echo "run filter_multimap skipped"
    """
}


process featurecounts {
  tag  { params.run_featurecounts  ? "$projectid" : "blank_run"  }
  cpus { params.run_featurecounts  ? params.cpu_max : params.cpu_min  }
  memory { params.run_featurecounts  ?  params.mem_max : params.mem_min  }

  // cpus 20
  // memory '350 GB'
  // time '96h'

	input:
  val x from filter_multimap_complete_ch.collect()
	val bams from bam_featurecounts_ch.collect()


  output:
	val "x" into featurecounts_complete_ch


  script:
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
  else if  ( params.species_global == "Rattus norvegicus" ){
      gtf = params.gtf_rn}
  else{
    gtf=""}



  if( params.run_featurecounts && params.paired)
    """
    mkdir -p ${featurecountsdir}
    # cd ${stardir}
    cd ${stardir_filtered}
    bamstring=\$(echo $bams | sed 's/,/ /g' | sed 's/\\[//g' | sed 's/\\]//g' )
    echo \${bamstring}
    echo "gtf: ${gtf}"
    featureCounts -T ${task.cpus} \\
      -t ${params.fcounts_feature} \\
      --extraAttributes gene_name,gene_type \\
      -a ${gtf} -g gene_id  \\
      -o ${featurecountsdir}/${projectid}_geneid.featureCounts.txt \\
      -p \\
      -s ${strand_numeric} \${bamstring}

    #find ${featurecountsdir} -user $USER -exec chmod g+rw {} +
    """
  else if( params.run_featurecounts && !params.paired)
    """
    mkdir -p ${featurecountsdir}
    cd ${stardir_filtered}
    bamstring=\$(echo $bams | sed 's/,/ /g' | sed 's/\\[//g' | sed 's/\\]//g' )
    echo \${bamstring}
    echo "gtf: ${gtf}"
    featureCounts -T ${task.cpus} \\
      -t ${params.fcounts_feature} \\
      --extraAttributes gene_name,gene_type \\
      -a ${gtf} -g gene_id  \\
      -o ${featurecountsdir}/${projectid}_geneid.featureCounts.txt \\
      -s ${strand_numeric} \${bamstring}
    """
  else
    """
    echo "featurecounts skipped"
    """
}

// Remove temooraryfiltered bams
// instead of changing directories that goes into multiqc-ctg_
// should not be in this multiqc sice file names are the same as non filtered bams

process rm_filter_mmap {
  tag  { params.run_featurecounts  ? "$sid" : "blank_run"  }
  cpus { params.run_featurecounts  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_featurecounts  ?  params.mem_standard : params.mem_min  }

  input:
  val x from featurecounts_complete_ch.collect()

  output:
  val "x" into rm_mmapfiltered_complete_ch

  // when: params.run_markdups

  script:
  if ( params.run_featurecounts )
    """
    if [ -d ${stardir_filtered} ]; then
      rm -rf ${stardir_filtered}
    fi
    """
  else
    """
    echo "rm_mmapfiltered"
    """
}




/* ===============================================================
  *     OTHER QC APPS
  =============================================================== */

process rnaseqmetrics {
  tag  { params.run_rnaseqmetrics  ? "$sid" : "blank_run"  }
  cpus { params.run_rnaseqmetrics  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_rnaseqmetrics  ?  params.mem_standard : params.mem_min  }


  input:
  val x from markdups_complete_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_rnaseqmetrics_ch

  output:
  val "x" into rnaseqmetrics_complete_ch

  // when: params.run_rnaseqmetrics

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
    rrna = params.picard_rrna_mm }
  else if ( species == "Rattus norvegicus" ){
      refflat = params.picard_refflat_rn
      rrna = params.picard_rrna_rn }
  else{
    refflat = ""
    rrna = ""
  }

  // if ( params.run_rnaseqmetrics && species == "Rattus norvegicus" )
  // else if ( pecies == "Rattus norvegicus" )
  // changed to NOT use the rrna file. not working anyway?
  if ( params.run_rnaseqmetrics && params.pipelineProfile == "uroscan")
    """
    echo "strand: ${strand}"
    echo "refflat file: ${refflat}"
    mkdir -p ${rnaseqmetricsdir}

    ## java -jar picard.jar CollectRnaSeqMetrics \\ ## old line
    picard CollectRnaSeqMetrics \\
        INPUT=${stardir}/${bam} \\
        OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
        REF_FLAT=${refflat} \\
        STRAND=${strand}

    #find ${rnaseqmetricsdir} -user $USER -exec chmod g+rw {} +

    """
  else if ( params.run_rnaseqmetrics && species == "Rattus norvegicus")
    """
    echo "strand: ${strand}"
    echo "refflat file: ${refflat}"
    mkdir -p ${rnaseqmetricsdir}

    ## java -jar picard.jar CollectRnaSeqMetrics \\ ## old line
    picard CollectRnaSeqMetrics \\
        INPUT=${stardir}/${bam} \\
        OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
        REF_FLAT=${refflat} \\
        STRAND=${strand}

    #find ${rnaseqmetricsdir} -user $USER -exec chmod g+rw {} +
    """
  else if ( params.run_rnaseqmetrics && params.pipelineProfile == "rnaseq")
    """
    echo "strand: ${strand}"
    echo "rrna file: ${rrna}"
    echo "refflat file: ${refflat}"
    mkdir -p ${rnaseqmetricsdir}

    picard CollectRnaSeqMetrics \\
      INPUT=${stardir}/${bam} \\
      OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
      REF_FLAT=${refflat} \\
      STRAND=${strand} \\
      RIBOSOMAL_INTERVALS=${rrna}

    #find ${rnaseqmetricsdir} -user $USER -exec chmod g+rw {} +
    """
  // temp workaround - ribosomal intervals file for Rat is not workling in v1.0

  else
    """
    echo "picard rnaseqmetrics skipped"
    """

}





// qualimap us quite redundant given the other apps Used
// If to be used a container that does uses Qualimap inastalled through pip3 NOT conda
process qualimap {
  tag  { params.run_qualimap  ? "$sid" : "blank_run"  }
  cpus { params.run_qualimap  ? params.cpu_high : params.cpu_min  }
  memory { params.run_qualimap  ?  params.mem_high : params.mem_min  }


  input:
  val x from rnaseqmetrics_complete_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_qualimap_ch

  output:
  val "x" into qualimap_complete_ch

  // when: params.run_qualimap

  script:
  // gtf used for featurecounts
  if ( species == "Homo sapiens" ){
    gtf = params.gtf_hs}
  else if  ( species == "Mus musculus" ){
    gtf = params.gtf_mm}
  else if  ( species == "Rattus norvegicus" ){
      gtf = params.gtf_rn}
  else{
    gtf=""}

  if ( params.run_qualimap )
    """
    mkdir -p ${qualimapdir}

    ## export JAVA_OPTS="-Djava.io.tmpdir=/data/tmp"
    ## /data/bnf/sw/qualimap_v2.2.1/qualimap --java-mem-size=12G rnaseq -bam /data/bnf/bam/rnaseq/21KF00020.STAR.sort.bam -gtf /data/bnf/ref/rsem/GRCh37/Homo_sapiens.GRCh37.75.gtf -pe -outdir /data/bnf/postmap/rnaseq/21KF00020.STAR.qualimap.folder
    # qualimap --java-mem-size=12G rnaseq -bam /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/star/21KF00082_Aligned.sortedByCoord.out.bam -gtf /projects/fs1/shared/uroscan/references/rsem/GRCh37/Homo_sapiens.GRCh37.75.gtf -pe -outdir /projects/fs1/shared/ctg-projects/uroscan/2021_024/nf-output/delivery/qualimap/21KF00082.STAR.qualimap.folder

    qualimap --java-mem-size=90G rnaseq -bam ${stardir}/${bam} -gtf ${gtf} -pe -outdir ${qualimapdir}/${sid}.STAR.qualimap.folder

    """
  else
    """
    echo "qualimap skipped"
    """
}



process rseqc {
  tag  { params.run_rseqc  ? "$sid" : "blank_run"  }
  cpus { params.run_rseqc  ? params.cpu_high : params.cpu_min  }
  memory { params.run_rseqc  ?  params.mem_high : params.mem_min  }


  input:
  val x from qualimap_complete_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_rseqc_ch

  output:
  val "x" into rseqc_complete_ch
  val "x" into rseqc_complete_report_ch
  // when: params.run_rseqc

  script:
  // gtf used for featurecounts
  if ( species == "Homo sapiens" ){
    rcqc_bed = params.rcqc_bed
    rcqc_housekeeping = params.rcqc_housekeeping
  }
  else if  ( species == "Mus musculus" ){
    rcqc_bed = params.rcqc_bed_mm
    rcqc_housekeeping = params.rcqc_housekeeping_mm
  }
  else{
    rcqc_bed=""}

  if ( params.run_rseqc && ( species == "Mus musculus" || species == "Homo sapiens") )
    """
    mkdir -p ${rseqcdir}

    geneBody_coverage.py \\
      -i ${stardir}/${bam} \\
      -r ${rcqc_housekeeping}\\
      -o ${rseqcdir}/${sid}.genebodycov

    inner_distance.py \\
      -i ${stardir}/${bam} \\
      -r ${rcqc_bed} \\
      -o ${rseqcdir}/${sid}.innerdistance

    """
  else
    """
    echo "skipped rseeqc"
    """
}













/* ===============================================================
  *      BLADDER REPORT
  =============================================================== */

// /usr/bin/Rscript -e "library(rmarkdown, lib='/home/petter/R/x86_64-pc-linux-gnu-library/3.4'); rmarkdown::render('/data/bnf/scripts/ort/bladderreport_ctg_anonymous.Rmd',params=list(sampleid='21KF00020', rsem_in='/data/bnf/premap/rnaseq/21KF00020_0.rsem',star_qc='/data/bnf/tmp/rnaseq/21KF00020_0.sort.bam.folder/Log.final.out', clarity_id='ALL557A36'),output_file='/data/bnf/postmap/rnaseq/21KF00020.STAR.bladderreport_anonymous.html')"


//
// if ( params.run_bladderreport == false ) {
//    Channel
// 	 .from("x")
//    .set{ bladderreport_complete_ch }
// }
process bladderreport {

  tag  { params.run_bladderreport  ? "$sid" : "blank_run"  }
  cpus { params.run_bladderreport  ? params.cpu_mid : params.cpu_min  }
  memory { params.run_bladderreport  ?  params.mem_mid : params.mem_min  }



  input:
  val x from rseqc_complete_report_ch.collect()
  val x from rsem_complete_report_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_bladderreport_ch

  output:
  val "x" into bladderreport_complete_ch

  // when: params.run_bladderreport

  script:
  bladderreport_scriptsdir = project_dir+'/bin/bladderreport'
  bladderreport_scriptname= params.bladderreport_scriptname


  if ( params.run_bladderreport )
  """
    mkdir -p ${bladderreportdir}/tmp_${sid}
    cp -r ${bladderreport_scriptsdir} ${bladderreportdir}/tmp_${sid}/
    cd ${bladderreportdir}/tmp_${sid}/bladderreport

    Rscript -e "library('rmarkdown'); \\
      rmarkdown::render( \\
        '${bladderreportdir}/tmp_${sid}/bladderreport/${bladderreport_scriptname}',  \\
        params = list(   \\
          sampleid='${sid}', \\
          rsem_in='${rsemdir}/${sid}.rsem.genes.results', \\
          star_qc='${stardir}/${sid}_Log.final.out', \\
          RIN='${RIN}', \\
          koncentration='${concentration}'),  \\
        output_file='${bladderreportdir}/${sid}.STAR.bladderreport_anonymous.html')"

    cd ${bladderreportdir}
    chromium --headless --disable-gpu --no-sandbox --print-to-pdf=${sid}.STAR.bladderreport.pdf ${bladderreportdir}/${sid}.STAR.bladderreport_anonymous.html

    mv -f ${bladderreportdir}/tmp_${sid}/bladderreport/${sid}.LundClassifier.rds ${bladderreportdir}/${sid}.LundClassifier.rds

    ## find ${bladderreportdir} -user $USER -exec chmod g+rw {} +

  """
  else
    """
    echo "run_bladderreport skipped"
    """
}










/* ===============================================================
* ===============================================================
  *     ----------- POST ANALYSIS SECTION -------
  ===============================================================
  =============================================================== */





/* ===============================================================
  *     ctg multiqc - on all analyses - runfolder included
  =============================================================== */
// This multiQC is for CTG infouse and not to the customer.
// The customer will obtain a lighter multiQC carried out below


process multiqc_ctg {
  //publishDir "${multiqcctgdir}", mode: 'copy', overwrite: 'true'
  tag  { params.run_multiqc_ctg  ? "$projectid" : "blank_run"  }
  // cpus 8
  // memory '64 GB'
  // time '3h'
  // echo debug_mode
  cpus { params.run_multiqc_ctg  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_multiqc_ctg  ?  params.mem_standard : params.mem_min  }



  input:
  val x from rm_mmapfiltered_complete_ch.collect()
  val x from fastqc_complete_ch.collect()
  //val x from markdups_complete_ch.collect()
  val x from rseqc_complete_ch.collect()
  val x from fastqscreen_complete_ch.collect()
  val x from rsem_complete_ch.collect()
  val x from bladderreport_complete_ch.collect()
  val x from salmon_complete_ch.collect()

  output:
  val "x" into multiqc_ctg_complete_ch
  val "x" into multiqc_ctg_complete_2_ch

  // when: params.run_multiqc_ctg

  script:
  if ( params.run_multiqc_ctg )
    """
      mkdir -p ${multiqcctgdir}
      cd ${output_dir}
      multiqc -n ${projectid}_multiqc_report \\
        --interactive \\
        -o ${multiqcctgdir} . ${runfolderdir}

      find ${multiqcctgdir} -user $USER -exec chmod g+rw {} +
    """
  else
    """
    echo "run_multiqc_ctg skipped"
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
// if ( params.run_setup_deliverytemp == false ) {
//    Channel
// 	 .from("x")
//    .set{ setup_deliverytemp_complete_ch }
// }




process stage_delivery {
  // cpus 4
  tag  { params.run_stage_delivery  ? "$projectid" : "blank_run"  }
  // memory '64 GB'
  // time '3h'
  cpus { params.run_stage_delivery  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_stage_delivery  ?  params.mem_standard : params.mem_min  }



  input:
  val x from multiqc_ctg_complete_ch.collect()

  output:
  val "x" into stage_delivery_complete_ch

  // when: params.run_setup_deliverytemp

  // add scripts
  // add sample sheets
  // move fastq

  script:
  if ( params.run_stage_delivery )
    """
    ## additional cleanups (star and bladderreport)
    ## -------------------------------------------
    if [ -d ${bladderreportdir} ]; then
      cd ${bladderreportdir}
      find . -type d -name "tmp_*" -exec rm -r {} +
    fi

    if [ -d ${stardir} ]; then
      cd ${stardir}
      find . -type d -name "*__STARtmp" -exec rm -r {} +
    fi

    if [ -d ${markdupstempdir} ]; then
      rm -rf ${markdupstempdir}
    fi



    ##  copy sample sheet to delivery
    ## ------------------------------
    mkdir -p ${deliverysamplesheets}
    if [ -f ${samplesheet} ]; then
      cp ${samplesheet} ${deliverysamplesheets}/
    fi

    ##  logs
    ## -----------------
    ## mkdir -p ${deliverylogs}


    ##  scripts dir (executables bins etc, version specific) and configs (project specific)
    ## --------------------------------------------------------------
    mkdir -p ${deliveryscripts}

    if [[ -d "${params.scriptsdir}" ]]; then
      cp -r ${params.scriptsdir} ${deliveryscripts}
    fi


    ## configs (project specific) as well as the rscript log config
    ##   --------------------------------------------------------------
    mkdir -p ${deliveryconfigs}
    if [[ -f "${project_dir}/nextflow.config.project.${projectid}" ]]; then
      cp ${project_dir}/nextflow.config.project.${projectid} ${deliveryconfigs}
    fi

    if [[ -f "${project_dir}/nextflow.config" ]]; then
      cp ${project_dir}/nextflow.config ${deliveryconfigs}
    fi
    if [[ -f "${runfolderdir}/log.rscript.samplesheet" ]]; then
      cp ${runfolderdir}/log.rscript.samplesheet ${deliveryconfigs}
    fi





    ## chmods
    ## --------
    find ${deliverytemp} -user $USER -exec chmod g+rw {} +


    """
  else
    """
    echo "run_setup_deliverytemp skipped"
    """


}

// move fastq files to delivery temp
// -----------------------------
process move_fastq {

  //cpus 6
  tag  { params.run_move_fastq  ? "$projectid" : "blank_run"  }
  // memory '64 GB'
  // time '3h'
  // echo debug_mode
  cpus { params.run_move_fastq  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_move_fastq  ?  params.mem_standard : params.mem_min  }


  input:
  val x from stage_delivery_complete_ch.collect()

  output:
  val "x" into move_fastq_complete_ch
  val "x" into move_fastq_complete_2_ch

  script:
    if ( params.sharedflowcell &&  params.run_move_fastq)
      """
        mkdir -p ${deliverytemp}/fastq
        if [ -d ${fastq_dir} ]; then
          echo "pooled run. moving fastq foldler only."
          mv -f ${fastq_dir} ${deliverytemp}/fastq
        fi
      """
    else if ( !params.sharedflowcell &&  params.run_move_fastq )
      """
      if [ -d ${bcl2fastq_dir} ]; then
        echo "non pooled data. moving comlplete bcl2fastq output foldler."
        mv -f ${bcl2fastq_dir} ${deliverytemp}
      elif [ -d ${fastq_dir} ]; then
        echo "non pooled run but cannot locate bcl2fastq_dir. moving fastq foldler only."
        mkdir -p ${deliverytemp}/fastq
        mv -f ${fastq_dir} ${deliverytemp}/fastq
      fi
      """
    else
     """
     echo "run_move_fastq skipped"
     """
}



// Run customer multiQC (on delivery temp folder only).
// -----------------------------
// Not to be confused with ctg-multiqc that is run on all analyses and runfolder

// if ( params.run_multiqc_delivery == false ) {
//    Channel
// 	 .from("x")
//    .set{ multiqc_complete_ch }
// }
process multiqc_delivery {

  tag  { params.run_multiqc_delivery  ? "${projectid}" : "blank_run"  }
  // cpus 6
  // memory '32 GB'
  // time '3h'
  // echo debug_mode
  cpus { params.run_multiqc_delivery  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_multiqc_delivery  ?  params.mem_standard : params.mem_min  }



  input:
  val x from move_fastq_complete_ch.collect()

  output:
  val "x" into multiqc_complete_ch

  // when: params.run_multiqc_delivery

  script:
  // if (! new File( mqcreport+'.html' ).exists() && params.run_multiqc_delivery)
 if ( params.run_multiqc_delivery  && params.pipelineProfile == "rnaseq"  )
  """
    ## remove if multiqc is already present from failed run. Will not overwrite ...
    rm -rf ${multiqcdelivery_dir}
    mkdir -p ${multiqcdelivery_dir}

    cd ${deliverytemp}
    multiqc -n ${mqcreport} \\
      --interactive \\
      -o ${multiqcdelivery_dir} .
  """
  else if ( params.run_multiqc_delivery  && params.pipelineProfile == "uroscan"  )
  """
    mkdir -p ${multiqcdelivery_dir}

    if [[ -d "${multiqcctgdir}" ]]; then
        cp -r ${multiqcctgdir}/* ${multiqcdelivery_dir}
    elif [[ -d "${ctg_save_dir}/qc/multiqc-ctg" ]]; then
        cp -r ${ctg_save_dir}/qc/multiqc-ctg/* ${multiqcdelivery_dir}
    fi
  """
  else
  """
    echo "skipping run_multiqc_delivery"
  """

}


//    md5 sum on delivery_dir
// -----------------------------
// generate md5 sum on all files included in delivery temp folder
// if ( params.run_md5sum_delivery == false ) {
//    Channel
// 	 .from("x")
//    .set{ md5sum_complete_ch }
// }
process md5sum_delivery {
  //cpus 8
  tag  { params.run_md5sum_delivery  ? "$projectid" : "blank_run"  }
  // memory '64 GB'
  // time '3h'
  cpus { params.run_md5sum_delivery  ? params.cpu_high : params.cpu_min  }
  memory { params.run_md5sum_delivery  ?  params.mem_high : params.mem_min  }


  input:
  val x from multiqc_complete_ch.collect()

  output:
  val "x" into md5sum_complete_ch

  // when: params.run_md5sum_delivery

  script:
  md5sumfile = delivery_dir + '/md5sum.txt'

  // if (! new File( md5sumfile ).exists() && params.run_md5sum_delivery)
  if ( params.run_md5sum_delivery )
  """
   cd ${deliverytemp}
   find . -type f -exec md5sum {} \\; > ${md5sumfile} ; echo
  """
  else
  """
  echo "skipping run_md5sum_delivery"
  """
  // else
  //   """
  //     echo "${md5sumfile} already exists. skipping."
  //   """
}





/* ===============================================================
  *     FINALIZE CTG SAVE & MOVE DELIVERY - save qc files and scripts
  =============================================================== */
// CTG should store multiQC (ctg-multiqc) and fastqc for all samples
// as of this version files are copied to ctg-qc dir. could change to the folder that also keep scripts and configs and sample sheets
// save

//  - multiQC (ctg_multiqc)
//  - fastQC
//
//  - Sample Sheets in ./samplesheets/
//  - Nextflow scripts, nextflow.config.project., rnaseq-main, nextflow.config, drivers, ./bin files etc (./scripts)

//  - logs,  (the final nextflow genereated onComplee is copied in that segion)
// if ( params.run_stage_ctg_save == false ) {
//    Channel
// 	 .from("x")
//    .set{ stage_ctg_save_complete_ch }
// }
process stage_ctg_save {
  //cpus 4
  tag "${projectid}"
  // memory '32 GB'
  // time '3h'
  cpus { params.run_stage_ctg_save  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_stage_ctg_save  ?  params.mem_standard : params.mem_min  }


  input:
  val x from multiqc_ctg_complete_2_ch.collect()
  val x from move_fastq_complete_2_ch.collect()

  output:
  val "x" into stage_ctg_save_complete_ch

  // when: params.run_stage_ctg_save


  script:
  if (params.run_stage_ctg_save)
  """

  ##  sample sheets
  ## -----------------
  cd ${project_dir}

  mkdir -p ${ctg_save_samplesheets}

  if [ -f ${samplesheet} ]; then
    cp ${samplesheet} ${ctg_save_samplesheets}
  fi

  ##  logs
  ##   -----------------
  mkdir -p ${ctg_save_logs}


  ## scripts dir (executables bins etc, version specific) and configs (project specific)
  ##   --------------------------------------------------------------
  mkdir -p ${ctg_save_scripts}

  if [[ -d "${params.scriptsdir}" ]]; then
    cp -r ${params.scriptsdir} ${ctg_save_scripts}
  fi


  ## configs (project specific) as well as the rscript log config
  ##   --------------------------------------------------------------
  mkdir -p ${ctg_save_configs}
  if [[ -f "${project_dir}/nextflow.config.project.${projectid}" ]]; then
    cp ${project_dir}/nextflow.config.project.${projectid} ${ctg_save_configs}
  fi

  if [[ -f "${project_dir}/nextflow.config" ]]; then
    cp ${project_dir}/nextflow.config ${ctg_save_configs}
  fi
  if [[ -f "${runfolderdir}/log.rscript.samplesheet" ]]; then
    cp ${runfolderdir}/log.rscript.samplesheet ${ctg_save_configs}
  fi


  ##  Move ctg qc dir from project folder to delivery location
  ##   --------------------------------------------------------------
  mv -f ${qcdir} ${ctg_save_dir}

  ##  duplicate the fastqc analyses from delivery dir
  ## -----------------------------------------
  cp -r ${fastqcdir} ${ctg_save_dir}/qc/


  ## CHOD
  find ${ctg_save_dir} -user $USER -exec chmod g+rw {} +

  """
  else
  """
    echo "skipping run_stage_ctg_save"
  """
}




// Finalize delivery_dir
// -----------------------------
// provess add README with dir size


process finalize_pipeline {

  tag  { params.run_finalize_pipeline  ? "${projectid}" : "blank_run"  }
  memory params.mem_min
  cpus params.cpu_min

  input:
  val x from md5sum_complete_ch.collect()
  val x from stage_ctg_save_complete_ch.collect()

  //val x from checkfiles_rscript_predelivery_ch.collect()

  output:
  val "x" into finalize_pipeline_ch



  script:
  if (params.run_finalize_pipeline)
  """


    ##  Move delivery temp dir from project folder to delivery location
    ##   --------------------------------------------------------------
    mkdir -p ${delivery_dir}
    mv -f ${deliverytemp}/* ${delivery_dir}


    cd ${delivery_dir}

    ## echo "ctg delivery complete"               > $readme
    ## echo "Project:   ${projectid}"             >> $readme
    ## du -ch -d 0 . | grep 'total'               >> $readme

    find ${delivery_dir} -user $USER -exec chmod g+rw {} +
    find ${project_dir} -user $USER -exec chmod g+rw {} +
    find ${ctg_save_dir} -user $USER -exec chmod g+rw {} +
  """
  else
  """
    echo "skipping run_finalize_pipeline"
  """
}
