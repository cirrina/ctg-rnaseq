

/* ===============================================================
  *      PARAMS FROM CONFIGS
  =============================================================== */

//  re-assign some params from nextflow.configs
// ----------------------------------------------------
//  project and run folders
projectid           =  params.projectid    // ctg project id. e.g. 2021_024
project_dir         =  params.project_dir   // .../shared/ctg-rnaseq/uroscan/2021_024 // NOT to be confused with project_dir
samplesheet         =  params.samplesheet           // name of simple sample sheet used for pipeline. Must includes file paths to fastq and bamsm, as well as species etc.
fastq_input_dir     =  params.fastq_input_dir            // subdirectory fastq files are located. For a default run the output location from blc2fastq. fastq-files will be read according to sample sheet. Defaults to <bcl2fastq_dir>/<projectid>
delivery_dir        =  params.delivery_dir
file(delivery_dir).mkdir() // main nexttlow work dir for output of analyses. Subdirectory of the project foilder. Files and folders will be moved and copiued from this folder upon pipeline  completion.
ctg_qc_dir          =  params.ctg_qc_dir
file(ctg_qc_dir).mkdir() // main nexttlow work dir for output of analyses. Subdirectory of the project foilder. Files and folders will be moved and copiued from this folder upon pipeline  completion.


/* ===============================================================
  *      DEFINE PATHS FROM INPUT PARAMS
  =============================================================== */
// Note that the primary output dir is the delivery_dir (to store all files to be delivered )
// the nf_tmp_dir is used to store qc-files and temp files that are not included in the delivery
nf_tmp_dir =  project_dir+'/nf-output' // nextflow (temp) output directory for files genetated with the Pipeline that are NOT to be delivered
file(nf_tmp_dir).mkdir() // main nexttlow work dir for output of analyses. Subdirectory of the project foilder. Files and folders will be moved and copiued from this folder upon pipeline  completion.


// Process & Module specific paths (Delivered to customer)
fastq_dir = delivery_dir+'/fastq'+'/'+projectid // To where fastqfiles are moved (delivery path). Follows bcl2fastq output naming. This to allow for transfer of a complete blc2fastq output (incuding stats folders etc) into the delivery dir.
stardir = delivery_dir+'/star'
salmondir = delivery_dir+'/salmon'
rsemdir = delivery_dir+'/rsem'
bladderreportdir = delivery_dir+'/bladderreport'
featurecountsdir = delivery_dir+'/featurecounts'
samplesheetsdir = delivery_dir+'/samplesheets'
deliveryscripts = delivery_dir+'/scripts'

// QC dirs
qcdir = delivery_dir+'/qc'

fastqcdir = qcdir+'/qc/fastqc'
multiqcdir  =  qcdir+'/qc/multiqc'
qualimapdir = qcdir+'/qualimap'
rseqcdir = qcdir+'/rseqc'
markdupsqcdir = qcdir+'/markdups'
rnaseqmetricsdir = qcdir+'/rnaseqmetrics'
fastqscreendir = qcdir+'/fastqscreen'

// tmp output dirs
markdupstempdir = nf_tmp_dir+'/markdups_bam_tmp'
stardir_filtered = nf_tmp_dir+'/star_filtered' // temporary folder used for filtered bam-files (no multi-map reads). These are used for featureCounts only due to problems with handling od large bam files.



/* ===============================================================
  *       create output and logdirs
  =============================================================== */
// log file for nextflow .onComplete
logfile   =  file( project_dir + '/' + 'log.nextflow.complete' )



/* ===============================================================
  *       CHECKS FILES AND PARAMS
  =============================================================== */


//  Check paramters
// -----------------------------
if (projectid         == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet   == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}


// Check if files and directories exist
checkPathParamList = [
  project_dir,
  nf_tmp_dir,
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
    nextflow tmp output dir :  ${nf_tmp_dir}
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
  logfile.text = msg_startup.stripIndent()
  logfile.append( msg_completed.stripIndent() )
  logfile.append( error )

  println( msg_completed )
}


// def msg_modules = """\
//     Run modules
//     ---------------------------------
//     fastqc            :    ${params.run_fastqc}
//     """
//     .stripIndent()
//
// println( msg_modules )



/* ===============================================================
  *    PROCESS SAMPLE SHEET & DEFINE CHANNELS
  =============================================================== */
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


/* ===============================================================
  *     Define Channels based from SampleSheet
  =============================================================== */
Channel
  .fromPath(sheet_nf)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.fastq_1, row.fastq_2, row.Species ) }
  .tap{ infoall }
  .set { fastq_ch }

Channel
  .fromPath(sheet_nf)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.bam, row.Strandness, row.Species, row.RIN, row.concentration ) }
  .tap { infobam }
  .into { bam_checkbam_ch; bam_qualimap_ch; bam_rseqc_ch; bam_bladderreport_ch; bam_rnaseqmetrics_ch }

Channel
    .fromPath(sheet_nf)
    .splitCsv(header:true)
    .map { row -> tuple( row.bam ) }
    .tap{ infoallfcounts }
    .set { bam_featurecounts_ch }

println " > Samples to process: "
println "[Sample_ID,fastq1,fastq2,species]"
infoall.subscribe { println "Info: $it" }




/* ===============================================================
  *    --- CHECK & MOVE FASTQ FILES ---
  =============================================================== */
// check if all expected files provided by samplesheet (columns fastq_1 and fastq_2) in fastq_input_dir
// if paired run, check bit
process checkfiles_fastq {

  tag  "$sid"
  cpus params.cpu_min
  memory params.mem_min

  input:
  set sid, read1, read2, species from fastq_ch

  output:
  // val "x" into checkfiles_fastq_complete_ch
  set sid, read1, read2, species into move_fastqc_ch

  script:
  if( params.paired_global )
    """
    file1=\$(find ${fastq_input_dir} -type f -name ${read1})
    file2=\$(find ${fastq_input_dir} -type f -name ${read2})
    if [[ -z \${file1} ]]; then
      echo "Warning: Cannot locate fastq_1 file supplied dir: ${fastq_input_dir}/${read1}"
      exit 2
    fi
    if [[ -z \${file2} ]]; then
    echo "Warning: Cannot locate fastq_1 file supplied dir: ${fastq_input_dir}/${read2}"
      exit 2
    fi
    """
  else
    """
    file1=\$(find ${fastq_input_dir} -type f -name ${read1})
    if [[ -z \${file1} ]]; then
      echo "Warning: Cannot locate fastq_1 file supplied dir: ${fastq_input_dir}/${read1}"
      exit 2
    fi
    """
}


// move fastq files to delivery dir {delivery_dir}/fastq/{project_id}
// --------------------------------
// Follow bcl2fastq file path convention, i.e. /fastq/${projectid}
// do NOT move if they already reside in a (subdirectory) of the delivery dir. IF so change the fastq_dir (used by nextflow) to same as fastq_input_dir
// Check if present in fastq_dir rather than fastq_input_dir - if file is already in expected folder and do not move
process move_fastq {

  tag  { params.run_move_fastq  ? "$projectid" : "blank_run"  }
  cpus { params.run_move_fastq  ? params.cpu_min : params.cpu_min  }
  memory { params.run_move_fastq  ?  params.mem_min : params.mem_min  }

  input:
  set sid, read1, read2, species from move_fastqc_ch  // from check fastq

  output:
  set sid, read1, read2, species into fastqc_ch
  //set sid, read1, read2, species into fastqc_ch

  script:
  if( params.paired_global && params.run_move_fastq)
    """
    mkdir -p ${fastq_dir}
    file1=\$(find ${fastq_dir} -type f -name ${read1})
    file2=\$(find ${fastq_dir} -type f -name ${read2})
    if [[ -z \${file1} ]]; then
      mv ${fastq_input_dir}/${read1} ${fastq_dir}
    fi
    if [[ -z \${file2} ]]; then
      mv ${fastq_input_dir}/${read2} ${fastq_dir}
    fi
    """
  else if( !params.paired_global && params.run_move_fastq)
    """
    mkdir -p ${fastq_dir}
    file1=\$(find ${fastq_dir} -type f -name ${read1})
    if [[ -z \${file1} ]]; then
      mv ${fastq_input_dir}/${read1} ${fastq_dir}
    fi
    """
  else
    """
    """
}



/* ===============================================================
  *      FASTQC
  =============================================================== */

process fastqc {
  tag  { params.run_fastqc  ? "$sid" : "blank_run"  }
  cpus { params.run_fastqc  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_fastqc  ? params.mem_standard : params.mem_min  }

  input:
  set sid, read1, read2, species from fastqc_ch  // from move_fastqc_ch

  output:
  val "x" into fastqc_complete_ch_1
  val "x" into fastqc_complete_ch_2
  val "x" into fastqc_complete_ch_3
  val "x" into fastqc_complete_ch_4
  set sid, read1, read2, species into star_ch
  set sid, read1, read2, species into salmon_ch
  set sid, read1, read2, species into rsem_ch
  set sid, read1, read2, species into fastqscreen_ch

  script:
  if ( params.paired_global && params.run_fastqc)
    """
    mkdir -p ${fastqcdir}
    echo "running fastqc in paired reads mode"
    fastqc ${fastq_dir}/${read1} ${fastq_dir}/${read2}  --outdir ${fastqcdir}
    """
  else if ( !params.paired_global && params.run_fastqc)
    """
    mkdir -p ${fastqcdir}
    echo "running fastqc in non paired reads mode "
    fastqc ${fastq_dir}/${read1}  --outdir ${fastqcdir}
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
    if ( params.paired_global ){
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
  *      -- SALMON  --
  =============================================================== */

process salmon  {
  tag  { params.run_salmon  ? "$sid" : "blank_run"  }
  cpus { params.run_salmon  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_salmon  ?  params.mem_standard : params.mem_min  }

  input:
  val x from fastqc_complete_ch_1.collect()
  set sid, read1, read2, species from salmon_ch // from checkfiles_fastq

  output:
  val "x" into salmon_complete_ch

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

  if ( params.paired_global && params.run_salmon )
    """
    salmon quant -l A \\
      -i  ${transcripts} \\
      -1  ${fastq_dir}/${read1} \\
      -2  ${fastq_dir}/${read2} \\
      -p  6 --validateMappings \\
      -o  ${salmondir}/${sid}_0.salmon.salmon \\
      --no-version-check
    """
  else if ( !params.paired_global && params.run_salmon )
    """
    salmon quant -l A \\
      -i  ${transcripts} \\
      -1  ${fastq_dir}/${read1} \\
      -p  6 --validateMappings \\
      -o  ${salmondir}/${sid}_0.salmon.salmon \\
      --no-version-check
    """
  else
    """
    echo "skipping salmon"
    """
}


/* ===============================================================
  *      RSEM ALIGNMENT
  =============================================================== */
process rsem {
  tag  { params.run_rsem  ? "$sid" : "blank_run"  }
  cpus { params.run_rsem  ? params.cpu_high : params.cpu_min  }
  memory { params.run_rsem  ?  params.mem_high : params.mem_min  }

  input:
  val x from fastqc_complete_ch_2.collect()
  set sid, read1, read2, species from rsem_ch // from checkfiles_fastq

  output:
  val "x" into rsem_complete_ch_1
  val "x" into rsem_complete_ch_2

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
  if ( params.paired_global ){
    rsemfiles = "${fastq_dir}/${read1} ${fastq_dir}/${read2}"
    paired='--paired-end'}
  else{
    rsemfiles = "${fastq_dir}/${read1}"
    paired=''}

  // strand - NOTE the uroscan pipe is run without strandness flag.
  if( params.strandness_global == "forward" )
    strand = 'forward'
  else if ( params.strandness_global == "reverse" )
    strand = 'reverse'
  else
    strand = 'none'

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
  *      STAR AND BAM SECTION
  =============================================================== */

process star  {
  tag  { params.run_star  ? "$sid" : "blank_run"  }
  cpus { params.run_star  ? params.cpu_high : params.cpu_min  }
  memory { params.run_star  ?  params.mem_high : params.mem_min  }

  input:
  val x from fastqc_complete_ch_3.collect()
  set sid, read1, read2, species from star_ch // from checkfiles_fastq

  output:
  val "x" into star_complete_ch

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

  if ( params.paired_global ){
      starfiles = "${fastq_dir}/${read1} ${fastq_dir}/${read2}" }
  else{
      starfiles = "${fastq_dir}/${read1}" }


  if ( params.run_star )
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
    echo "skipping star"
    """

}


process checkfiles_bam {
  tag  { params.run_checkfiles_bam  ? "$sid" : "blank_run"  }
  cpus params.cpu_min
  memory params.mem_min

  input:
  val x from star_complete_ch.collect() // checkbam_ch - when star is completed
  set sid, bam, strand, species, RIN, concentration from bam_checkbam_ch

  output:
  val "x" into checkfiles_bam_complete_ch
  set sid, bam, strand, species, RIN, concentration into bam_indexbam_ch

  script:
  if( params.run_checkfiles_bam )
    """
    if [[ ! -f ${stardir}/${bam} ]]; then
      echo "Warning: Cannot locate bam file ${stardir}/${bam}"
      exit 2
    fi
    """
  else
    """
    echo "file check overridden"
    """
}


process index_bam {
  tag  { params.run_index_bam  ? "$sid" : "blank_run"  }
  cpus { params.run_index_bam  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_index_bam  ?  params.mem_standard : params.mem_min  }

  input:
  val x from checkfiles_bam_complete_ch.collect()
  set sid, bam, strand, species, RIN, concentration from bam_indexbam_ch

  output:
  val "x" into indexbam_complete_ch
  set sid, bam, strand, species, RIN, concentration into bam_markdups_ch

  script:
  if ( params.run_index_bam )
    """
    cd ${stardir}
    echo "${stardir}/${bam}"
    samtools index -bc ${stardir}/${bam}
    """
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
  val "x" into markdups_complete_ch_1
  val "x" into markdups_complete_ch_2
  set sid, bam, strand, species, RIN, concentration into bam_filter_multimap_ch

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
    """
  else
    """
    echo "run markdups skipped"
    """
}



/* ===============================================================
  *     FEATURECOUNTS SECTION
  =============================================================== */
// Filter bams on only primary mapped sequence using 0x104 flag.
// Featurecounts will only use primary mapped read anlyway - but tends to crash if multi mapped sequences are included

process filter_multimap {
  tag  { params.run_featurecounts  ? "$sid" : "blank_run"  }
  cpus { params.run_featurecounts  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_featurecounts  ?  params.mem_standard : params.mem_min  }

  input:
  val x from markdups_complete_ch_1.collect()
  set sid, bam, strand, species, RIN, concentration from bam_filter_multimap_ch

  output:
  val "x" into filter_multimap_complete_ch

  script:
  if ( params.run_featurecounts )
    """
    echo "bam: ${bam}"
    mkdir -p ${stardir_filtered}
    cd ${stardir_filtered}
    samtools view -b -F 0x104 ${stardir}/${bam} >  ${stardir_filtered}/${bam}
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

	input:
  val x from filter_multimap_complete_ch.collect()
	val bams from bam_featurecounts_ch.collect()

  output:
	val "x" into featurecounts_complete_ch

  script:
  fcounts_feature   =  'exon'

  if( params.strandness_global == "forward" ) {
    strand_numeric = 1 }
  else if ( params.strandness_global == "reverse" ) {
    strand_numeric = 2 }
  else {
    strand_numeric = 0 }

  // gtf used for featurecounts
  if ( params.species_global == "Homo sapiens" ) {
    gtf = params.gtf_hs }
  else if  ( params.species_global == "Mus musculus" ) {
    gtf = params.gtf_mm }
  else if  ( params.species_global == "Rattus norvegicus" ) {
    gtf = params.gtf_rn }
  else {
    gtf="" }

  if( params.run_featurecounts && params.paired_global )
    """
    mkdir -p ${featurecountsdir}
    cd ${stardir_filtered}
    bamstring=\$(echo ${bams} | sed 's/,/ /g' | sed 's/\\[//g' | sed 's/\\]//g' )
    echo \${bamstring}
    echo "gtf: ${gtf}"
    featureCounts -T ${task.cpus} \\
      -t ${fcounts_feature} \\
      --extraAttributes gene_name,gene_type \\
      -a ${gtf} -g gene_id  \\
      -o ${featurecountsdir}/${projectid}_geneid.featureCounts.txt \\
      -p \\
      -s ${strand_numeric} \${bamstring}
    """
  else if( params.run_featurecounts && !params.paired_global )
    """
    mkdir -p ${featurecountsdir}
    cd ${stardir_filtered}
    bamstring=\$(echo $bams | sed 's/,/ /g' | sed 's/\\[//g' | sed 's/\\]//g' )
    echo \${bamstring}
    echo "gtf: ${gtf}"
    featureCounts -T ${task.cpus} \\
      -t ${fcounts_feature} \\
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



/* ===============================================================
  *     OTHER QC APPS
  =============================================================== */

process rnaseqmetrics {
  tag  { params.run_rnaseqmetrics  ? "$sid" : "blank_run"  }
  cpus { params.run_rnaseqmetrics  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_rnaseqmetrics  ?  params.mem_standard : params.mem_min  }

  input:
  val x from markdups_complete_ch_2.collect()
  set sid, bam, strand, species, RIN, concentration from bam_rnaseqmetrics_ch

  output:
  val "x" into rnaseqmetrics_complete_ch

  script:
  if ( strand == "forward" ) {
    strand_input="FIRST_READ_TRANSCRIPTION_STRAND" }
  else if ( params.strandness_global == "reverse" ) {
    strand_input="SECOND_READ_TRANSCRIPTION_STRAND" }
  else {
    strand_input="NONE" }

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

  if ( params.run_rnaseqmetrics && params.pipelineProfile == "uroscan")
    """
    echo "strand: ${strand_input}"
    echo "refflat file: ${refflat}"
    mkdir -p ${rnaseqmetricsdir}

    picard CollectRnaSeqMetrics \\
        INPUT=${stardir}/${bam} \\
        OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
        REF_FLAT=${refflat} \\
        STRAND=${strand_input}
    """
  else if ( params.run_rnaseqmetrics && species == "Rattus norvegicus")
    """
    echo "strand: ${strand_input}"
    echo "refflat file: ${refflat}"
    mkdir -p ${rnaseqmetricsdir}

    picard CollectRnaSeqMetrics \\
        INPUT=${stardir}/${bam} \\
        OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
        REF_FLAT=${refflat} \\
        STRAND=${strand_input}
    """
  else if ( params.run_rnaseqmetrics && (params.pipelineProfile == "rnaseq" || params.pipelineProfile == "rnaseq_total"))
    """
    echo "strand: ${strand_input}"
    echo "rrna file: ${rrna}"
    echo "refflat file: ${refflat}"
    mkdir -p ${rnaseqmetricsdir}

    # note that rrna ribosomal intervals file seem not to work in present state.
    picard CollectRnaSeqMetrics \\
      INPUT=${stardir}/${bam} \\
      OUTPUT=${rnaseqmetricsdir}/${sid}_bam.collectRNAseq.metrics.txt \\
      REF_FLAT=${refflat} \\
      STRAND=${strand_input}
      ## RIBOSOMAL_INTERVALS=${rrna}
    """
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

  script:
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
  *      UROSCAN - BLADDER REPORT
  =============================================================== */
// Bladderreport need a temporary folder for each analysis since temp filles with non unique names are generated

process bladderreport {
  tag  { params.run_bladderreport  ? "$sid" : "blank_run"  }
  cpus { params.run_bladderreport  ? params.cpu_mid : params.cpu_min  }
  memory { params.run_bladderreport  ?  params.mem_mid : params.mem_min  }

  input:
  val x from rsem_complete_ch_1.collect()
  set sid, bam, strand, species, RIN, concentration from bam_bladderreport_ch

  output:
  val "x" into bladderreport_complete_ch

  script:
  bladderreport_scriptsdir = project_dir+'/bin/bladderreport'
  bladderreport_scriptname = params.bladderreport_scriptname

  if ( params.run_bladderreport )
    """
    mkdir -p ${bladderreportdir}/tmp_${sid}
    cp -r ${bladderreport_scriptsdir} ${bladderreportdir}/tmp_${sid}/
    cd ${bladderreportdir}/tmp_${sid}/bladderreport

    ## Run Rscript to generate markdown pdf report
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

    ## Run Chromium to generate html from pdf
    cd ${bladderreportdir}
    chromium --headless --disable-gpu --no-sandbox --print-to-pdf=${sid}.STAR.bladderreport.pdf ${bladderreportdir}/${sid}.STAR.bladderreport_anonymous.html
    mv -f ${bladderreportdir}/tmp_${sid}/bladderreport/${sid}.LundClassifier.rds ${bladderreportdir}/${sid}.LundClassifier.rds
    """
  else
    """
    echo "run_bladderreport skipped"
    """
}



/* ===============================================================
* ===============================================================
  *     ----------- POST ANALYSIS - MULTIQC  -------
  ===============================================================
  =============================================================== */



process multiqc {
  tag  { params.run_multiqc  ? "$projectid" : "blank_run"  }
  cpus { params.run_multiqc  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_multiqc  ?  params.mem_standard : params.mem_min  }

  input:
  val x from featurecounts_complete_ch.collect()
  val x from fastqc_complete_ch_4.collect()
  val x from rseqc_complete_ch.collect()
  val x from fastqscreen_complete_ch.collect()
  val x from rsem_complete_ch_2.collect()
  val x from salmon_complete_ch.collect()

  output:
  val "x" into multiqc_complete_ch

  script:
  if ( params.run_multiqc )
    """
    ## use -f flag to overwrite if multiqc is already present from failed run.
    cd ${delivery_dir}
    multiqc -n ${mqcreport} \\
      --interactive \\
      -f \\
      -o ${multiqc_dir} .

    """
  else
    """
    echo "run_multiqc skipped"
    """
}


/* ===============================================================
* ===============================================================
  *     ----------- POST ANALYSIS - nd5sum  -------
  ===============================================================
  =============================================================== */
// perform md5sums on FASTQ and BAM directories

process md5sum {
  tag  { params.run_md5sum_delivery  ? "$projectid" : "blank_run"  }
  cpus { params.run_md5sum_delivery  ? params.cpu_high : params.cpu_min  }
  memory { params.run_md5sum_delivery  ?  params.mem_high : params.mem_min  }

  input:
  val x from multiqc_complete_ch.collect()

  output:
  val "x" into md5sum_complete_ch

  script:
  if ( params.run_md5sum )
    """
    if [[ -d "${fastq_dir}" ]]; then
      find ${fastq_dir} -type f -exec md5sum {} \\; > ${fastq_dir}/md5sum_fastq.txt ; echo
    fi
    if [[ -d "${stardir}" ]]; then
      find ${stardir} -type f -exec md5sum {} \\; > ${stardir}/md5sum_star.txt ; echo
    fi
    if [[ -d "${featurecountsdir}" ]]; then
      find ${featurecountsdir} -type f -exec md5sum {} \\; > ${featurecountsdir}/md5sum_featurecounts.txt ; echo
    fi
    """
  else
    """
    echo "skipping run_md5sum_delivery"
    """
}



/* ===============================================================
* ===============================================================
  *     ----------- POST ANALYSIS COPY N CLEANUP  -------
  ===============================================================
  =============================================================== */


/* ===============================================================
  *     Stage the delivery folder - Scripts SampleSheets & logs
  =============================================================== */
// remove tmp files
// copy SampleSheets & Scripts

process stage_delivery {
  tag  { params.run_stage_delivery  ? "$projectid" : "blank_run"  }
  cpus { params.run_stage_delivery  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_stage_delivery  ?  params.mem_standard : params.mem_min  }

  input:
  val x from multiqc_complete_ch.collect()
  val x from bladderreport_complete_ch.collect()

  output:
  val "x" into stage_delivery_complete_ch

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

    if [[ -d ${stardir_filtered} ]]; then
      rm -rf ${stardir_filtered}
    fi

    ##  copy sample sheet to delivery
    ## ------------------------------
    mkdir -p ${samplesheetsdir}
    cp -r ${samplesheet} ${samplesheetsdir}
    cp -r ${sheet_nf} ${samplesheetsdir}


    ##  scripts & configs  (executables bins etc, version specific) and configs (project specific)
    ## --------------------------------------------------------------
    mkdir -p ${deliveryscripts}
    cp -r ${project_dir}/rnaseq-driver ${deliveryscripts}
    cp -r ${project_dir}/rnaseq-main.nf ${deliveryscripts}
    cp -r ${project_dir}/bin ${deliveryscripts}
    cp -r ${project_dir}/nextflow.config.* ${deliveryscripts}



    """
  else
    """
    echo "run_setup_deliverytemp skipped"
    """

}





/* ===============================================================
  *     FINALIZE CTG SAVE & MOVE DELIVERY - save qc files and scripts
  =============================================================== */

process finalize_pipeline {

  tag  { params.run_finalize_pipeline  ? "${projectid}" : "blank_run"  }
  memory params.mem_min
  cpus params.cpu_min

  input:
  val x from  stage_delivery_complete_ch.collect()

  output:
  val "x" into finalize_pipeline_complete_ch

  script:
  if (params.run_finalize_pipeline)
    """

    ## Copy QC files to ctg-qc
    ## -----------------------
    cp -r ${multiqcdir} ${ctg_qc_dir}
    cp -r ${fastqcdir} ${ctg_qc_dir}
    cd ${delivery_dir}

    ## Chmod all dirs
    ## -----------------------
    find ${delivery_dir} -user $USER -exec chmod g+rw {} +
    find ${project_dir} -user $USER -exec chmod g+rw {} +
    find ${ctg_qc_dir} -user $USER -exec chmod g+rw {} +
    """
  else
    """
    echo "skipping run_finalize_pipeline"
    """
}
