

# Samplesheet names - outputed from samplesheet rscript
# samplesheet_nextflow="${project_dir}/SampleSheet-${projectid}-demux.csv"
# samplesheet_ctg="${project_dir}/SampleSheet-${projectid}-ctg.csv"
# samplesheet_="${project_dir}/SampleSheet-${projectid}-original.csv"




# project specific config file


pipeline="${pipelineName} ${pipelineVersion}"



## Read PoolName from samplesheet
## ----------------------------------------------
poolname=$(awk -F, '$1 == "PoolName"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... PoolName: $poolname"

## Read SharedFlowCell from samplesheet
## ----------------------------------------------






# #############################################################################
#  1. IF RESUME FLAG - Try initiate Nextflow Pipeline using RESUME
# #############################################################################
if [[ $resume == "true" ]]; then
  if [[ $runfolder_mode == "true" ]]; then
    echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
    echo " Nextflow resume flag -r must not be used while executing sctipt in Illumina runfolder."
    echo " Move to project directory and initiate resume"
    echo""; echo ""
    exit_abnormal
  fi
  echo " RESUME mode"
  echo " ... trying to resume nextflow"
  echo " ... Hail mary !!! "
  nohup nextflow run $nf_script -c $nf_config_project -profile ${pipelineProfile} --resume > log.nextflow.rnaseq &
  echo " ... :  nextflow run $nf_script -c $nf_config_project --resume"
fi








###################################################################
# 1. IF `projectfolder_setup_mode` - Prime Project Work folder & chek sample sheet
###################################################################
#
if [[ $resume == "false" ]]; then
  if [[ $projectfolder_setup_mode == "true" ]]; then




    # parse and generate sample sheets using Rscript in singularity container
    #########################################################################
    echo ""
    echo ""
    echo " Running Rscript: iem-SampleSheet-processor "
    # echo " -------------------------------------------  "
    echo " ... force_lane            : $force_lane";
    echo " ... force_index           : $force_index";



  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHANGE WHEN DEPOLY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #${project_dir}/bin/iem-samplesheet-processor.R \
    $singcmd_rscript ${project_dir}/bin/iem-samplesheet-processor.R \
          --project_id ${projectid} \
          --sample_sheet ${samplesheet} \
          --output_demux_sheet ${samplesheet_demux}  \
          --output_ctg_sheet ${samplesheet_ctg} \
          --bin_dir "${project_dir}/bin/" \
          --force_lane ${force_lane} \
          --force_replace_index ${force_index}
    #echo " -------------------------------------------  "


    if [ ! -f "${execdir}/log.rscript.samplesheet" ]; then
      echo ""; echo "Error:"
      echo " R script 'bin/iem-samplesheet-processor.R' did not produce log.rscript.samplesheet"; echo ""; echo ""
      exit_abnormal
      # echo "'${project_dir}' already exists. Overwriting this folder."
    fi

    ## copy logdfile from runfolder to project_dir
    chmod 770 ${execdir}/log.rscript.samplesheet
    cp ${execdir}/log.rscript.samplesheet ${project_dir}


    echo ""
    echo ""
    echo " Parsing rscript logfile "
    echo ""

    logfile_iem="${execdir}/log.rscript.samplesheet" ## path and file suffix 'log.rscript.samplesheet' is dedfined within 'iem-samplesheet-processor.R'.
    assay=$(cat ${logfile_iem} | grep 'Assay' |  cut -f2 -d",")
    echo "  ... Assay:  ${assay}"
    instrument_type=$(cat ${logfile_iem} | grep 'Instrument_type' |  cut -f2 -d",")
    echo "  ... Instrument_type:  ${instrument_type}"
    index_adapters=$(cat ${logfile_iem} | grep 'Index_Adapters' |  cut -f2 -d",")
    echo "  ... Index_Adapters:  ${index_adapters}"
    strandness=$(cat ${logfile_iem} | grep 'Strandness' |  cut -f2 -d",")
    echo "  ... strandness:  ${strandness}"
    paired=$(cat ${logfile_iem} | grep 'Paired' |  cut -f2 -d",")
    echo "  ... paired:  ${paired}"
    n_samples=$(cat ${logfile_iem} | grep 'number_samples' |  cut -f2 -d",")
    echo "  ... n_samples:  ${n_samples}"
    echo ""


    #  End with Warning if samplesheet was not accepted
    #######################################################################
    echo " Checking validity of samplesheet"
    samplesheet_accepted=$(cat ${logfile_iem} | grep 'samplesheet_accepted' |  cut -f2 -d",")

    if [ $samplesheet_accepted == "FALSE" ]; then
      echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
      echo " ... IEM Sample Sheet was NOT accepted"
      echo " ... Please check logfile for details: "
      echo " ... ... ${logfile_iem}"
      echo ""; echo ""
      cat ${logfile_iem}
      echo ""
      #echo "- Please specify correct samplesheet, or create a CTG_SampleSheet.csv in current runfolder"
      exit_abnormal
    fi
    echo "  ...  OK!"



    # Prime nextflow configuration file -- nextflow.params_${projectid} --
    #######################################################################
    echo ""
    echo " Writing nextflow parameters to project config"

    ## Write nextflow params to file
    echo ""  > $nf_config_project
    echo "//  nextflow configuration file"                              >> $nf_config_project
    echo "//  Project:  ${projectid}"                                   >> $nf_config_project
    echo ""                                                             >> $nf_config_project
    echo "//  project specific parameters"                              >> $nf_config_project
    echo "//  will override params in 'nextflow.config' "               >> $nf_config_project
    echo "//"                                                           >> $nf_config_project
    echo ""                                                             >> $nf_config_project
    echo " params {"                                                      >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  // Pipeline & container                                     " >> $nf_config_project
    echo "  pipelineProfile   =  '${pipelineProfile}'                   " >> $nf_config_project
    echo "  scripts_dir        =  '${scripts_dir}'                        " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  // Project & assay specififc                                " >> $nf_config_project
    echo "  projectid           =  '${projectid}'                       " >> $nf_config_project
    echo "  project_dir          =  '${project_dir}'                      " >> $nf_config_project
    echo "  runfolder           =  '${runfolder}'                       " >> $nf_config_project
    echo "  fastqdir            =  '${fastqdir}'                        " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  sharedflowcell      =   ${sharedflowcell}                   " >> $nf_config_project
    echo "  species_all_samples      =  '${species_all_samples}'                  " >> $nf_config_project
    echo "  n_samples           =  '${n_samples}'                       " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  assay               =  '${assay}'                           " >> $nf_config_project
    echo "  instrument_type     =  '${instrument_type}'                 " >> $nf_config_project
    echo "  index_adapters      =  '${index_adapters}'                  " >> $nf_config_project
    echo "  paired              =   ${paired}                           " >> $nf_config_project
    echo "  strandness          =  '${strandness}'                      " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  //  samplesheets                                            " >> $nf_config_project
    echo "  samplesheet           =  '${samplesheet_ctg}'               " >> $nf_config_project
    echo "  samplesheet_demux     =  '${samplesheet_demux}'             " >> $nf_config_project
    echo "  samplesheet_original  =  '${samplesheet_original}'          " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  // bcl2fastq specific                                       " >> $nf_config_project
    echo "  runfolderdir        =  '${runfolderdir}'                    " >> $nf_config_project
    echo "  bcl2fastq_dir       =  '${bcl2fastq_dir}'              " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  // root directories                                         " >> $nf_config_project
    echo "  project_root        =  '${project_root}'                    " >> $nf_config_project
    echo "  delivery_root       =  '${delivery_root}'                   " >> $nf_config_project
    echo "  ctg_save_root       =  '${ctg_save_root}'                   " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  // Flags for individual nexflow processes                     " >> $nf_config_project
    echo "  //   NOTE! debugging only                                     " >> $nf_config_project
    echo "  //   'fastq_custom' flag will determine the run bcl2fastq flag   " >> $nf_config_project
    echo "  //   All other modules are 'true' by default                  " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo "  run_blcl2fastq          =  ${run_blcl2fastq}                 " >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    echo " }"                                                             >> $nf_config_project
    echo ""                                                               >> $nf_config_project
    ## Singularity container depending on rnaseq or rnaseq
    # echo " // Define SLURM specs                                        " >> $nf_config_project
    # echo " process {                                                    " >> $nf_config_project
    # echo "  executor = 'slurm'                                          " >> $nf_config_project
    # echo "  container = '${singularity_container}'                      " >> $nf_config_project
    # echo "  time = '48h'                                                " >> $nf_config_project
    # echo "  cpus = '16'                                                 " >> $nf_config_project
    # echo "  memory = '100 GB'                                           " >> $nf_config_project
    # echo " } " >> $nf_config_project


    # If skip demox flag is FALSE- override the config file
    #################################################
    # If skip demux flag is given in execution.
  echo " ... profile set to: ${pipelineProfile}'"
    # if [[ ${pipelineProfile} == "rnaseq" ]]; then
      # echo " ... ... forcing 'run_rsem' to 'false'"
      # sed "s/run_rsem.*/run_rsem            =  false/g" $nf_config > tmp.txt ; mv tmp.txt $nf_config
      # echo " ... ... forcing 'run_bladderreport' to 'false'"
      #sed "s/run_bladderreport.*/run_bladderreport            =  false/g" $nf_config > tmp.txt ; mv tmp.txt $nf_config
    # fi

    #if [[ ${pipelineProfile} == "uroscan" ]]; then
  #    echo " ... ... forcing 'run_featurecounts' to 'false'"
  #    sed "s/run_featurecounts.*/run_featurecounts            =  false/g" $nf_config > tmp.txt ; mv tmp.txt $nf_config
#    fi




  fi ## end priming section -  if projectfolder_setup_mode is true



  ###################################################################
  # 2. `projectfolder_setup_mode` FALSE
  ###################################################################
  ## NOTE: THIS IS NOTE SAME AS -RESUME. resume mode is only to try re-starting a faulty nextflow run.
  # if script is not initiated in Illumina runfolder - assume project is already proimed and ready to start.
  # This will include:
  # i. samplesheets
  # ii. nextflow sctipts, configs and bins
  # iii. fastq file path supplied
  # This

  if [[ $projectfolder_setup_mode == "false" ]]; then
    echo ""
    echo " projectfolder_setup_mode is FALSE "
    echo " ... executiondir is same as ProjectId in SampleSheet"
    echo " ... assuming all configuration and sctipt files are present "
    echo " ... assuming valid ctg-style sample sheet is supplied in project config : $nf_config_project "

    # If skip demox flag is FALSE- override the config file
    #################################################
    # If skip demux flag is given in execution.
    if [[ ${run_blcl2fastq} == "false" ]]; then
      echo " ... forcing 'run_blcl2fastq' to 'false'"
      sed "s/run_blcl2fastq.*/run_blcl2fastq            =  false/g" $nf_config_project > tmp.txt ; mv tmp.txt $nf_config_project
    fi
    if [[ ${fastq_custom} == "true" ]]; then
      echo " ... custom fastq dir. setting this fastqdir in config file ${nf_config_project}"
      sed "s|fastqdir.*|fastqdir            =  '${fastqdir}'|g" $nf_config_project > tmp.txt ; mv tmp.txt $nf_config_project
    fi
  fi
