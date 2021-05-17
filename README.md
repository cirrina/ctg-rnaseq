# ctg-rnaseq

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
