## v1.0 to-do
+ Add 'Rattus norvegicus' references to ls4
+ Add 'Rattus norvegicus' species file paths in nextflow.config



file(ctg_save_dir).mkdir()





## Known issues

- Should add checks for multiple prpcesses that they do not run if expected files are present. This featuer should not clash with the resume of nextflow. i.e. if resume is true, folders/.files may exist but not comlpete. this check must therefore be done on complete expected output. May be lesss of a problem with running e.g. multiqc etc
Possible sollution is resume mode that would turn all initial flags to true.
One other solutione is to add "start_from" flags. to use when -resume is needed  but not working.
- nf run flags. as of now fiddeling with true & false parameters can be done. nf should log all paramaters somehow.


- 1. -resume flag.
- 2. -if no resume set to FALSE if dir/file does not exist
- 3. optional flags. module process flags. Could be set by modes. (mode full, mode algn, mode mininal)

Delivery mode:
alignment or NOT aligned data.


deliverydir

/data
/qc
perform only md5sum on the data folder




OTHER stuff
´´delivery on lfs. deliver per PROJECT rather than customer/user


# ctg-rnaseq

This pipeline will process RNA seq sequencing data from Illumina raw data runfolders.

A typical run is executed through two bash scripts in sequential order.

1) 'rnaseq-primer':
his script will prepare an analysis for the 'rnaseq-driver' script. It will further retrieve the nextflow sctipts directory that include the main nf script 'ctg-rnaseq' as well as the 'nextflow.con fig' and /bin directory that contain script files such as the 'iem-sampleseet-processor.R'. This Rscript cross-checks the inputed IEM style samplesheet and generates two sample sheets (demux and plain) that are in correct format. Finally, the primer will create a project specific nextflow configuration file that carry project-specific paramaters of the project ('nextflow.params.<project_id>'). This will be the primary file for changing nextflow parameters.

This script shoud ideally be initiated within the Illumina runfolder, bu can also be initiated in a folder with the same name as set by the -i parameter (e.g. if the script is re-run with modified samplesheet or if no demux shall be performed)


2) 'rnaesq-driver'
This script will run the ctg-rnaseq nextflow workflow. It uses samplesheets generated by the rnaseq-primer. The script also uses the nextflow parameter file 'nextflow.params.<project_id>' generated by the primer.

This script must be initiated in the project directory (same project directory as set by the -i flag).




## Specification of individual scripts





Input for rnaseq-primer sctipt is.
a) project id. (-i)
a) Sample sheet (-s)


### rnaseq-primer

- Script location:
  /projects/fs1/shared/ctg-drivers

- Initiation:
  Preferrably withinIllumina runfolder but can also be initiated in the exact project directory that should be run. e.g. '/projects/fs1/shared/ctg-projects/ctg-rnaseq/2021_044'

- Usage:
  -i, projectid, e.g. '2021_044'
  -s, samplesheet. Should be in IEM format and preferably generated by Illumina IEM software.

  -l, force_lane (optional): Set to 1 or 2 if to force demux to only run one lane. This if lane divider is used AND if lane is NOT specified in sample sheet. This parameter will override the Lane column in sample sheet"
  -f, force_replace_index (optional): Set this flag if to force iem samplesheet Rscript to overwrite index id and sequence informations based on index well id. The script uses '/ctg-rnaseq/workflow/bin/checklist-index.csv' to cross-check indexes.

- Script details.

This script will prepare rnaseq analysis pipeline by:

  a) setting up a project folder and copy nextflow ctg-rnaseq sctripts from 'ctg-pipelines/ctg-rnaseq/workflow/' folder.
  b) run 'ien-sampleseet-processor.R' Rscript and generate two sample sheets.
  c) generate project specific parameters file, 'nextflow.params.<project_id>'.



### iem-sampleseet-processor.R


  i)   check the format of the input IEM style sample sheet
  ii)  replace illegal characters (this section could be refined/expanded)
  iii) check indexes against a database of valid indexes. Sample sheet [Header] must specify 'Index Adapters' and 'Instrument Type', e.g. NovaSeq or NovaSeq1.0 to get the correct index pairs.
  iv)  check adapters and additional trimming. Sample sheet [Header] must specify 'Index Adapters' and 'Assay'
  v)   generate a bcl2fastq demux sample sheet (default name: 'SampleSheet-<project_id>-demux.csv').
  vi)  generate plain sample sheet (default name: 'SampleSheet-<project_id>-ctg.csv').
    This sample sheet carries information of the project, in paticular fastq and bam file names that are checked used in the downstream nextflow script.
  iv)



### rnaseq-driver

- Script location:
  /projects/fs1/shared/ctg-drivers

- Initiation:
    Within the exact directory as project id that should be run. e.g. '/projects/fs1/shared/ctg-projects/ctg=rnaseq/2021_044/
    This directory is usially created when running the 'rnaseq-primer' script

- Usage:
  -i projectid, e.g. '2021_044'
  -r resume

- Required files:
  - 'SampleSheet-<project_id>-ctg.csv', e.g. SampleSheet-2021_044_test-ctg.cs. Sample sheet generated by the primer script. Plain sample sheet used by e.g. fastqc, star etc.
  - 'SampleSheet-<project_id>-demux.csv', e.g. SampleSheet-2021_044_test-demux.cs. Sample sheet generated by the primer script. Used for bcl2fastq demux.
  - 'nextflow.params.<project_id>'. Nextflow parameters file. Used when initating the main nextflow sctipt.






### iem-sampleseet-processor.R
i)   check the format of the input IEM style sample sheet
ii)  replace illegal characters (this section could be refined/expanded)
iii) check indexes against a database of valid indexes. Sample sheet [Header] must specify 'Index Adapters' and 'Instrument Type', e.g. NovaSeq or NovaSeq1.0 to get the correct index pairs.
iv)  check adapters and additional trimming. Sample sheet [Header] must specify 'Index Adapters' and 'Assay'
iv)  generate a bcl2fastq demux sample sheet (default name: 'SampleSheet-<project_id>-demux.csv').
iv)



- Input files
    -

The ctg-rnaseq pipeline is




## TESTS
Test Lanes demux. How does no-lane splitting work (). e.g. 2021_030.



## Questions
+ Chmods of generated files and fodlers !! should be full access by all within the lsens4 group. Add chmods to script? What about executables...
+ Log files should first be generated in project dir. Then on completeion moved? Or always genarated elsewhere?


## Random Comments to script


## Log files



## Processes for bcl2fastq

When threading is supported,the software uses the follow defaults to manage threads for processing: u Fourth reads for reading thedata.
u Fourth reads for writing thedata.
u Twenty percent of threads for demultiplexing data.
u One hundred percent of threads for processing demultiplexed data.

The file io threads aretypicallyinactiveandconsumeminimalprocessingtime.Processingdemultiplexeddataisallocatedonethreadpercentralprocessingunit(CPU)topreventidleCPUs,resultinginmorethreadsthanCPUsbydefault.ConsiderationsforMultipleThreadsWhenusingprocessingoptionstoassignmultiplethreads,considerthefollowinginformation:uThemostdemandingstepistheprocessingstep(-poption).Assignthisstepthemostthreads.uThereadingandwritingstagesaresimpleanddonotneedmanythreads.Thisconsiderationisimportantforalocalharddrive.Toomanythreadscausetoomanyparallelread
-writeactionsandsuboptimalperformance.uUseonethreadperCPUcoreplussomeextra.ThismethodpreventsCPUsfrombeingidleduetoathreadbeingblockedwhilewaitingforanotherthread.uThenumberofthreadsdependsonthedata.Ifyouspecifymorewritingthreadsthansamples,theextrathreadsdonoworkbutcosttimeduetocontextswitching




bam_rnaseqmetrics_ch
rnaseqmetrics_ch ... x


bam_indexbam_ch
indexbam_ch

bam_checkbam_ch
checkbam_ch

bam_markdups_ch
markdups_ch



bam_index_complete_ch
markdups_complete_ch
rnaseqmetrics_complete_ch

bam_featurecounts_ch



if [ "${species}" == "Homo sapiens" ]
then
  gtf=${params.gtf_hs}
elif [ "${species}" == "Mus musculus" ]
then
  gtf=${params.gtf_mm}
else
  echo "Warning: featurecounts, Species not recognized. ${species}"
fi



featurecountsdir = outputdir+'/featurecounts'
stardir = outputdir+'/star'
fastqcdir = outputdir+'/fastqc'
markdupsqcdir = outputdir+'/markdups'
rnaseqmetricsdir = outputdir+'/rnaseqmetrics'

fastqdir

pooled
runfolderdir

markdupsdir = outputdir+'/markdups_bam_tmp'




### OLD ST

# Set Species

if [ "${species}" == "Homo sapiens" ]
then
  echo $species
  genome=${params.star_genome_hs}
elif [ "${species}" == "Mus musculus" ]
then
  genome=${params.star_genome_mm}
else
  echo "Warning: Species not recognized. ${species}"
fi



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

## log_root='/Users/david/tasks/rnaseq_test/nftest/logs' # logfolder='/Users/david/tasks/rnaseq/rnaseq_test/logs'


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


// Define SLURM specs
// process {
//   executor='slurm'
//   container = 'singularity/rnaSeqTools/rnaseqtools.dl.0.1.sif'
//   time='3h'
//   cpus='4'
// //
// //   withName:mkfastq {
// //     time='24h'
// //     cpus='16'
// //     memory='110 GB'
// //   }
// //   withName:count {
// //     time='2d'
// //     cpus='20'
// //     memory='120 GB'
// //   }
// //   withName:count_nuc {
// //     time='2d'
// //     cpus='20'
// //     memory='120 GB'
// //   }
// //   withName:aggregate {
// //     time='2d'
// //     cpus='16'
// //     memory='120 GB'
// //
// //   }
//  }


Caused by:
  Process `bcl2fastq (null)` terminated with an error exit status (127)

Command executed:

  bcl2fastq -R /projects/fs1/nas-sync/upload/210510_A00681_0363_BH3MN5DRXY \
           --sample-sheet /projects/fs1/shared/ctg-projects/ctg-rnaseq/2021_044_test/SampleSheet-2021_044_test-demux.csv \
           --no-lane-splitting  \
           -r 1 \
           -p 4  \
           -w 1  \
           --output-dir /projects/fs1/shared/ctg-projects/ctg-rnaseq/2021_044_test/fastq/2021_044_test


           cd $run
     singularity exec --bind /projects/fs1/ /projects/fs1/shared/ctg-containers/bulkRNA/bulkRNA_STARv2.7.6a.sif bcl2fastq \
         -R $run \
         --sample-sheet $ssheet \
         --no-lane-splitting  \
         -r 1 \
         -p 16 \
         -w 1 \
         -o $run/Fastq_Raw/



singularity exec -bind /projects/fs1/ /projects/fs1/shared/ctg-containers/bulkRNA/bulkRNA_STARv2.7.6a.sif bcl2fastq \
    -R /projects/fs1/nas-sync/upload/210510_A00681_0363_BH3MN5DRXY \
    --sample-sheet /projects/fs1/shared/ctg-projects/ctg-rnaseq/2021_044_test/SampleSheet-2021_044_test-demux.csv \
    --no-lane-splitting  \
    -r 1 \
    -p 4  \
    -w 1  \
    --output-dir /projects/fs1/shared/ctg-projects/ctg-rnaseq/2021_044_test/fastq/2021_044_test



    ## Check force_lane !! NOT WORKING (Rscript will test this. ..)
    # allowed_lane=('0' '1' '2')
    # if [[ " "${allowed_lane[@]}" " != *" "$force_lane" "* ]] ;then
    #     #echo "$assay: ok"
    # # else
    #     echo ""
    #     echo "Failed check for '-l' flag. "
    #     echo "$force_lane: not recognized. Valid names are:"
    #     echo "${allowed_lane[@]/%/,}"; echo ""; echo ""
    #     exit 1
    # fi
