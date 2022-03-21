#!/bin/bash -l


#######################################
# driver for  ctg-rnaseq pipeline
#######################################

## INPUT: CTG IEM Sample Sheet
## INPUT: FASTQ files (default fastq parth is output from bcl2fastq analysis)


# scripts_root='/Users/david/scripts' #

# script execution dir. If generate workfolder & copy scripts, this WILL be checked & warn against the script version rovided in samplesheet
# if initiating script in project workfolder, the script here is used & version will not be checked
script_exec_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # Where is this script located...
script_exec_dir=$(cd ${script_exec_dir} && pwd -P) # needed in case script exec dir is a symlink
# echo ${script_exec_dir}

## Root directories (based on ourr set folder naming conventions on lsens)

scripts_root="/projects/fs1/shared/ctg-pipelines"
project_root='/projects/fs1/shared/ctg-projects'
delivery_root='/projects/fs1/shared/ctg-delivery' ## should be added pipelineProfile/ProjectID
ctg_qc_root='/projects/fs1/shared/ctg-qc/' ## should be added pipelineProfile/ProjectID

## Script & config names
nf_script="rnaseq-main.nf"


#######################
#  == Initiation ==
#######################
# fastq_custom=true ## set to false if no fastq_input_dir is supplied (will defalt to {project_})
prime_projectfolder_mode='false' ## -p true,  if to generate project folder and configs, but not initiate the nextgflow
resume='false' ## used to resume nextflow pipeline, Need to be executed within project workfolder.
exec_dir=$(pwd)

# usage message
usage() {
    echo ""
    echo " rnaseq-driver [ -s samplesheet ] [ -p prime_projectfolder_only ] [ -r resume ] [ -f fastq_input_dir ] \
[ -h help ] "  1>&2

    echo "------------------- "
    echo " samplesheet                -s : IEM style laboratory SampleSheet. Project and pipeline specific parameters must be added "
    echo " fastq_input_dir            -f : Set to a full path where all fastq files are located. individual filenames should be specified in ctg sample sheet. This flag will set 'run_blcl2fastq' to 'false' and all flags related to demux to false "
    echo " prime_projectfolder_mode:  -p : Set to 'true' if you do not start nextflow. A dry run that will generate project folder, config files etc but will not initiate nextflow."
    echo " resume                     -r : Nextflow-specific. If to resume nextflow run. Can only be used when executing in Project work dir, NOT from Illumina runfolder"
    echo " help                       -h : print help message"
    echo " (see https://github.com/cirrina/rnaseq/blob/main/README.md for details)"
    echo "------------------- "
    echo ""
}

exit_abnormal() {
  usage
  exit 1
}


################################################
# == 1 ==   Check input arguments
################################################

while getopts "s:f:rph" opt; do
    case $opt in
      s) samplesheet=$OPTARG
	          ;;
      f) fastq_input_dir=$OPTARG
            ;;
      r) resume='true'
       	    ;;
      p) prime_projectfolder_mode='true'
            ;;
      h) exit_abnormal
        ;;
      \?) echo echo ""; echo " Error ctg-rnaseq : ";"Invalid option -$OPTARG" >&2
        exit_abnormal ;;
      :) echo ""; echo "Error:"; echo " -${OPTARG} requires an argument!"echo ""; echo ""
	     exit 1 ;;
    esac
done



shift "$(( OPTIND -1 ))"

## Check Sample Sheet. if file is present in work directory.
if [[ -z ${samplesheet} ]]; then
  echo ""; echo ""; echo "Error:"
  echo "You must specify sample sheet (in current dir): '-s' flag. "; echo ""
  exit 1
fi
if [ ! -f ${samplesheet} ]; then
  echo ""; echo ""; echo "Error:"
  echo "Sample Sheet does not exist (in current dir)"
  #echo "- Please specify correct samplesheet, or create a CTG_SampleSheet.csv in current runfolder"
  exit 1
fi

## If -f is defined (fastq_input_dir), set fastq_custom to true
if [[ -z ${fastq_input_dir} ]] && [[ ${prime_projectfolder_mode} != "true"  ]]; then
  echo ""; echo ""; echo "Error:"
  echo "You must specify path where fastq files are located: '-f' flag. "; echo ""
  exit 1
fi
## check if supplied fastq_input_dir is valid
fastq_input_dir=$(realpath $fastq_input_dir)
if [[ ! -d ${fastq_input_dir} ]] && [[ ${prime_projectfolder_mode} != "true"  ]]; then
  echo ""; echo ""; echo "Error:"
  echo "fastq_input_dir (-f) does not exist: ${fastq_input_dir} "; echo ""
  exit 1
fi



################################################
##   == 2 == Read paramters from SampleSheet
################################################

# check that samplesheet is a proper CTG IEM style sheet
if [[ $(cat ${samplesheet} | grep "\[Header\]" | cut -f1 -d ",") != "[Header]" ]]
  then
    echo " ERROR: ${samplesheet} does not seem to be a valid SampleSheet - lacking first row [Header]."
    exit 0
fi
if [[ $(cat ${samplesheet} | grep "\[Data\]" | cut -f1 -d ",") != "[Data]" ]]
  then
    echo " ERROR: ${samplesheet} does not seem to be a valid SampleSheet - lacking [Data] section."
    exit 0
fi
echo " Reading SampleSheet params:"
echo "  ... samplesheet: $samplesheet"


## Read PROJECT ID from samplesheet
## ----------------------------------------------
projectid=$(awk -F, '$1 == "ProjectID"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... project id: $projectid"

if [ -z "$projectid" ]; then
  echo ""; echo "Error: "
  echo " ProjectId is not properly supplied in samplesheet."
  echo " CTG Project id must be given as 'ProjectID' within the [Header] section of sample sheet (mind the semantics)"; echo""; echo ""
  exit 1
fi
if [[ $projectid =~ ' ' ]]; then
  echo ""; echo ""; echo "Error: "
  echo " Warning: project id must not inlude white space ' '"; echo""; echo ""
  exit 1
fi
if [[ $projectid =~ \\. ]]; then
  echo ""; echo ""; echo "Error: "
  echo " Warning: project id must not inlude a dot '.'"; echo""; echo ""
  exit 1
fi


## Read CTG PIPELINE parameters from samplesheet.
## ----------------------------------------------
pipelineName=$(awk -F, '$1 == "PipelineName"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... pipelineName: $pipelineName"
pipelineVersion=$(awk -F, '$1 == "PipelineVersion"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... pipelineVersion: $pipelineVersion"
pipelineProfile=$(awk -F, '$1 == "PipelineProfile"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... pipelineProfile: $pipelineProfile"
runFolder=$(awk -F, '$1 == "RunFolder"' ${samplesheet} | awk -F, '{print $2}')

# PipelineName
if [ -z "$pipelineName" ]; then
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
  echo "Warning: 'PiepelineName' is not properly supplied in samplesheet. \
 Must be given as 'PipelineName, rnaseq' within the [Header] section of sample sheet"; echo""; echo ""
  exit 1
fi

# PipelineVersion
## --------------------------
if [ -z "$pipelineVersion" ]; then
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR : "; echo " ---------------------------  ";
  echo "Warning: 'PiepelineVersion' is not properly supplied in samplesheet. \
 Must be given as 'PipelineVersion', 'current' or '2.1.5' within the [Header] section of sample sheet"; echo""; echo ""
  exit 1
fi

## Profile uroscan or rnaseq or rnaseq_total
if [ -z "$pipelineProfile" ]; then
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR : "; echo " ---------------------------  ";
  echo "Warning: 'pipelineProfile' is not properly supplied in samplesheet.";
  echo " Should be 'uroscan, rnaseq or rnaseq_total' within the [Header] section of sample sheet"; echo""; echo ""
  exit 1
fi
if [[ $pipelineProfile == "rnaseq_mrna" ]]; then
  echo "  ... pipelineProfile: $pipelineProfile"
elif [[ $pipelineProfile == "rnaseq_total"  ]]; then
  echo "  ... pipelineProfile: $pipelineProfile"
elif [[ $pipelineProfile == "uroscan"  ]]; then
  echo "  ... pipelineProfile: $pipelineProfile"
else
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR : "; echo " ---------------------------  ";
  echo "Warning: 'pipelineProfile' is not properly supplied in samplesheet.";
  echo " Should be 'uroscan, rnaseq_mrna or rnaseq_total' within the [Header] section of sample sheet"; echo""; echo ""
  exit 1
fi


## Read params from samplesheet [Header].
## -------------------------------

## species = ReferenceGenome
# indicate shared species, i.e. for all samples. Will only be used for modules that input ALL samples. Default is to use the species option supplied as column header in [Data] section.
species_global=$(awk -F, '$1 == "Species"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... species_global, from samplesheet [Header]: ${species_global}"
if [[ ${species_global} != "Homo sapiens" ]] && [[ ${species_global} != "Mus musculus"  ]] && [[ ${species_global} != "Rattus norwegicus"  ]]; then
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
  echo " species_global is not properly supplied in samplesheet.";
  echo " Must be supplied as 'Species' within the [Header] section of sample sheet";
  echo " Accepted values are Homo sapiens, Mus musculus and Rattus norwegicus";
  echo""; echo ""
  exit 1
fi

paired_global=$(awk -F, '$1 == "Paired"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... paired_global, from samplesheet [Header]: ${paired_global}"
if [[ ${paired_global} != "true" ]] && [[ ${paired_global} != "false"  ]]; then
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
  echo " paired_global is not properly supplied in samplesheet."; echo " Must be supplied as 'Paired' 'true' OR 'false' within the [Header] section of sample sheet"; echo""; echo ""
  exit 1
fi

strandness_global=$(awk -F, '$1 == "Strandness"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... strandness_global, from samplesheet [Header]: ${strandness_global}"
if [[ ${strandness_global} != "forward" ]] && [[ ${strandness_global} != "reverse"  ]]; then
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
  echo " strandness_global is not properly supplied in samplesheet."; echo " Must be supplied as 'Strandness' 'forward' OR 'reverse' within the [Header] section of sample sheet"; echo""; echo ""
  exit 1
fi

sharedflowcell=$(awk -F, '$1 == "SharedFlowCell"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... SharedFlowCell: ${sharedflowcell}"
if [[ ${sharedflowcell} != "true" ]] && [[ ${sharedflowcell} != "false"  ]]; then
  echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
  echo " SharedFlowCell is not properly supplied in samplesheet. \
  Must be supplied as 'SharedFlowCell' 'true' OR 'false' within the [Header] section of sample sheet"; echo" This to indicate wether multiple customer projects are present within the runfolder data"; echo ""
  exit 1
fi



######################################################
##  == 3 == Set paths and filenames based on input params
######################################################

# set project work dir
project_dir="${project_root}/${pipelineName}/${pipelineProfile}/${projectid}"
delivery_dir="${delivery_root}/${pipelineName}/${pipelineProfile}/${projectid}"
ctg_qc_dir="${ctg_qc_root}/${pipelineName}/${pipelineProfile}/${projectid}"

## scripts - where and what version of pipelie sctips to copy to project work dir

## As from version 2.2.10, scripts are primarily run using 'current'. The script dir path not used as from samplesheet. the path is saved in the project specific params file
# e.g. scripts_dir='/projects/fs1/shared/ctg-dev/pipelines/rnaseq/current'
# scripts_dir="${scripts_root}/${pipelineName}/${pipelineVersion}"
script_exec_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # Where is this script located... may be a symlink 'current'
script_exec_dir=$(cd ${script_exec_dir} && pwd -P) # needed in case script exec dir is a symlink


## if -p true, prime_projectfolder_mode. must be executed outside execution dir.
if [[ ${prime_projectfolder_mode} == "true" && "${exec_dir}" == "${project_dir}" ]]; then
  echo ""; echo " Can not prime workfolder with the same path as the one where script is executed"; echo ""
  echo " You must specify a different ProjectID wihtin your SampleSheet"
  exit 1
fi

nf_config_project="${project_dir}/nextflow.config.params.${projectid}"


# #############################################################################
#  == 4.a == IF RESUME FLAG - Try Resume Nextflow Pipeline using RESUME
# #############################################################################
if [[ $resume == "true" ]]; then
  if [[ "${script_exec_dir}" != "${project_dir}" ]]; then
    echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
    echo " Nextflow resume flag -r must used when executing script from the project folder defined by ProjectID in samplesheet"
    echo " Move to project directory and initiate resume"
    echo""; echo ""
    exit 1
  fi
  echo " RESUME mode"
  echo " ... trying to resume nextflow"
  echo " ... Hail mary !!! "
  nohup nextflow run ${nf_script} -c ${nf_config_project} -profile ${pipelineProfile} --resume > log.nextflow.rnaseq &
  echo " ... :  nextflow run ${nf_script} -c ${nf_config_project} --resume"
  exit 0
fi




################################################
##  == 4.b ==  projectfolder_setup_mode = true
################################################
# if executed outside the project dir, then initiate in projectfolder_setup_mode
# This will generate a project folder and copy the scrips (of specfied version) to work folder
# IF workfolder and scripts exists these will be overwritten.

if [[ "${exec_dir}" != "${project_dir}" ]]; then

  echo ""; echo ""; echo " Executing script outside defined project workfolder."
  projectfolder_setup_mode=true
  echo "  ... ... projectfolder_setup_mode: ${projectfolder_setup_mode}"

  ##  Check if script_exec_dir exist
  if [[ ! -d ${script_exec_dir} ]]; then
    echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
    echo " script_exec_dir does not exist: ${script_exec_dir} "; echo ""
    echo " Make sure PipelineName and PipelineVersion are correctly supplied in SampleSheet"; echo " AND that they match a directory whtin the scripts_root folder: ${scripts_root}."
    exit 1
  fi

  ## verify that script execution dir matches the path in which this script was executed
  # if [[ "${script_exec_dir}" != "${scripts_dir}" ]]; then
  #   echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
  #   echo " PipelineVersion in SampleSheet does not match the scripts_dir where this script is initiated : "
  #   echo " ... scripts version specific dir (given in SampleSheet) : ${scripts_dir} ";
  #   echo " ... script execution dir : ${script_exec_dir} "; echo ""
  #   exit 1
  # fi

  ## Warnings & Prompt.
  echo ""
  echo "  ... ... scripts_dir:  ${script_exec_dir}"
  echo "  ... ... project workfolder:  ${project_dir} "
  echo ""

  # Prompt user to approve running in current directory and input
  # read -p "
  # proceed? (y/n)
  #  ... " prompt
  #
  # if [[ $prompt != "y" ]]; then
  #     echo ""; echo " Exiting!! "
  #     exit 0
  # fi
  # echo ""

  if [[ -d ${project_dir} ]]; then
    echo ""
    echo "... ... WARNING! project work folder already exists!"
    echo "... ... remove this to continue."
    echo ""; echo "Exiting!! "
    exit 0
    #
    #   read -p "  ... Warning!!
    # ... ... project work folder already exists!
    # ... ... scripts & configs will be overwritten!
    #
    #
    # proceed? (y/n)
    #  ... " prompt
    #   if [[ $prompt != "y" ]]; then
    #       echo ""; echo "Exiting!! "
    #       exit 0
    #   fi
  fi
  echo ""
  echo ""

  ## Create project directory
  ## ------------------------
  mkdir -p ${project_dir}
  ## Copy scripts from version specific pipeline scripts dir
  cp -r ${script_exec_dir}/* ${project_dir}/ # copy all scripts to workfolder. Will overwrite netflow.config
  ## Copy samplesheet to project workfolder
  cp ${samplesheet} ${project_dir}
  cd ${project_dir}


  # Create nextflow configuration file -- nextflow.params_${ProjectID} --
  #######################################################################


  echo ""
  echo " ... Writing nextflow parameters to project-specific config: ${nf_config_project}"

  ## Write nextflow params to file
  echo ""  > ${nf_config_project}
  echo "//  nextflow configuration file"                              >> ${nf_config_project}
  echo "//  Project:  ${projectid}"                                   >> ${nf_config_project}
  echo ""                                                             >> ${nf_config_project}
  echo "//  project specific parameters"                              >> ${nf_config_project}
  echo "//  will override params in 'nextflow.config' "               >> ${nf_config_project}
  echo "//"                                                           >> ${nf_config_project}
  echo ""                                                             >> ${nf_config_project}
  echo " params {"                                                      >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  // Pipeline                                                 " >> ${nf_config_project}
  echo "  pipelineName       =  '${pipelineName}'                      " >> ${nf_config_project}
  echo "  pipelineProfile    =  '${pipelineProfile}'                   " >> ${nf_config_project}
  echo "  pipelineVersion    =  '${pipelineVersion}'                   " >> ${nf_config_project}
  echo "  script_execution_dir  =  '${script_exec_dir}'               " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  // Project and run directories                               " >> ${nf_config_project}
  echo "  projectid          =  '${projectid}'                       " >> ${nf_config_project}
  echo "  project_dir        =  '${project_dir}'                      " >> ${nf_config_project}
  echo "  delivery_dir       =  '${delivery_dir}'                      " >> ${nf_config_project}
  echo "  ctg_qc_dir         =  '${ctg_qc_dir}'                   " >> ${nf_config_project}
  echo "  fastq_input_dir    =  '${fastq_input_dir}'                        " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  runFolder          =   '${runFolder}'                        " >> ${nf_config_project}
  echo "  sharedflowcell     =   ${sharedflowcell}                   " >> ${nf_config_project}
  echo "  species_global     =  '${species_global}'                     " >> ${nf_config_project}
  echo "  paired_global      =   ${paired_global}                     " >> ${nf_config_project}
  echo "  strandness_global  =  '${strandness_global}'              " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  //  samplesheets                                            " >> ${nf_config_project}
  echo "  samplesheet        =  '${samplesheet}'               " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  // root directories                                         " >> ${nf_config_project}
  echo "  project_root       =  '${project_root}'                    " >> ${nf_config_project}
  echo "  delivery_root      =  '${delivery_root}'                   " >> ${nf_config_project}
  echo "  ctg_qc_root      =  '${ctg_qc_root}'                   " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo " }"                                                             >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}

  #  Priming of project folder complete
  ## --------------------------------
  echo "";echo ""
  echo "  Project primed"
  echo "  ... Project dir        :  ${project_dir}"
  echo "  ... Nextflow config    :  ${nf_config_project}"
  echo "  ... samplesheet        :  ${samplesheet}"
  echo "  ... Scripts dir        :  ${script_exec_dir}"
  echo ""

  chmod -R 775 ${project_dir}
  cd ${project_dir}

  ## Set profile defaults in nextflow.config
  ## also, remove the bladderrepoort folder for profiles rnaseq_total and rnaseq_mrna. Only used for uroscan
  if [[ $pipelineProfile == "rnaseq_mrna" ]]  ; then
    ## update the fastq_input_dir (-f argument) in config file
    rm -rf ${project_dir}/bin/bladderreport
    sed "s|run_rsem.*|run_rsem            =  false|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
    sed "s|run_salmon.*|run_salmon            =  false|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
    sed "s|run_bladderreport.*|run_bladderreport            =  false|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
  elif [[ $pipelineProfile == "rnaseq_total"  ]]; then
    rm -rf ${project_dir}/bin/bladderreport
    sed "s|run_rsem.*|run_rsem            =  false|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
    sed "s|run_salmon.*|run_salmon           =  false|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
    sed "s|run_bladderreport.*|run_bladderreport       =  false|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
  elif [[ $pipelineProfile == "uroscan" ]]; then
    sed "s|run_rsem.*|run_rsem            =  true|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
    sed "s|run_salmon.*|run_salmon            =  true|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
    sed "s|run_bladderreport.*|run_bladderreport            =  true|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
    sed "s|run_featurecounts.*|run_featurecounts            =  false|g" nextflow.config > tmp.txt ; mv tmp.txt nextflow.config
  fi


  if [[ $prime_projectfolder_mode == "true" ]]; then
    echo ""
    echo " prime_projectfolder_mode is TRUE "
    echo " ... exiting without submitting nextflow run"
    echo " ... go to Project folder:  ${project_dir}"
    echo " ... inspect & edit your configs, then start run using rnaseq-drive, "
    exit 0
  fi
fi # end if projectfolder_setup_mode


#
################################################
##  == 4.c ==  Script executed within project folder (projectfolder_setup_mode = false)
################################################
if [[ "${exec_dir}" == "${project_dir}" ]]; then
  ## add expected files - nextflow.main config.project ...
  echo " Script executed within the defined project workfolder."
  echo " ... will NOT overwrite scripts and configs. Assume all scripts are present."
  echo " ... -s SampleSheet must be same as defined in  ${nf_config_project}"

  ## make sure that specified SampleSheet matches the one specified in nextflow config file
    ## For completion these must be the same.
    ## Note that parameters in the config files are the only that matters.
    ## If changing samplesheet params - make sure that params in config files have been changed as\
    ## If a new samplesheet is applied, consider priming a new project folder with -p glag

    samplesheet_cfg="$(grep ".*samplesheet.*=.*" ${nf_config_project} | cut -d \' -f2)"
    echo "${samplesheet_cfg}"

    if [ "${samplesheet_cfg}" != "${samplesheet}" ]; then
      echo " Warning: "
      echo " SampleSheet argument (-s): ${samplesheet}"
      echo " ${samplesheet}"
      echo " does not match samplesheet supplied in ${nf_config_project}"
      echo " ${samplesheet_cfg}"
      echo ""
      echo " make sure that these are the same"
      echo " Note that parameters in the config files are the only that matters when (re-)initiating pipeline whithin project folder "
      echo " If changing samplesheet params - make sure that params in config files also have been changed"
      echo " If changing multiple samplesheet params - you may want to delete and re-prime the project folder using -p flag"
      exit 1
    fi

  ## update the fastq_input_dir (-f argument) in config file
  sed "s|fastq_input_dir.*|fastq_input_dir            =  \'$fastq_input_dir\'|g" ${nf_config_project} > tmp.txt ; mv tmp.txt ${nf_config_project}
fi




###################################################################
# == 5a == Final Check expected f\les needed for nextflow initiation
###################################################################
## Check if nextflow config files are present
cd ${project_dir}
if [ ! -f ${nf_config_project} ]; then # | ! -f "${exec_dir}/nextflow.config" | ! -f "${exec_dir}/rnaseq.nf"
  echo ""; echo "Error:"
  echo "'${nf_config_project}' does not exist in current dir. Generate this in projectfolder_setup_mode or crete one from template."; echo ""; echo ""
  exit 1
fi
if [ ! -f "${project_dir}/nextflow.config" ]; then # | ! -f "${exec_dir}/nextflow.config" | ! -f "${exec_dir}/rnaseq.nf"
  echo ""; echo "Error:"
  echo "'nextflow.config' file does not exist in project work directory"; echo ""; echo ""
  exit 1
fi

###################################################################
# == 5b == Update run_modules deoending on PipelineProfile
###################################################################



################################################
##  == 6 ==  Execute the NextFlow rnaseq-main.nf sctipt
################################################
echo ""; echo "";
echo " Initiating nextflow pipeline"
echo " ... pipelineName      : ${pipelineName}";
echo " ... pipelineVersion   : ${pipelineVersion}";
echo " ... pipelineProfile   : ${pipelineProfile}";
echo " ... project id    : ${projectid}";
echo " ... project dir   : ${project_dir}";
echo " ... delivery dir  : ${delivery_dir}";
echo ""; echo "";

# Prompt user to approve running in current directory and input
# echo "";
# read -p "
#   Initiating nextflow pipeline - proceed? (y/n)
#   ... " prompt
# if [[ $prompt != "y" ]]; then
#   echo ""; echo "Exiting!! "
#   exit 0
# fi

## intiate the nextflow command. include project specific config & profile -p
cd ${project_dir}
module load Java
module load nextflow/19.04.1
module load Singularity
nextflow run ${nf_script} -c ${nf_config_project} -profile ${pipelineProfile} > log.nextflow.progress &

echo ""; echo ""
echo "  Running :   nextflow run ${nf_pipe} -c ${nf_config_project} -profile ${pipelineProfile}"
echo "";
echo "  ########################## "
echo "      S U B M I T T E D "
echo "  ########################## "
echo "  Logfile :  ${project_dir}/log.nextflow.progress "
echo ""
