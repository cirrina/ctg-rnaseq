#!/bin/bash

#####################
# sc-rna-10x driver #
#####################

### This script will
### * Run the sc-rna-10x pipeline on data in current runfolder
### * Modify standard nextflow.config to project specific
### * Generate project folder in shared/ctg-delivery/sc-rna-10x
###  -> Here it will store nextflow.config, nf-pipeline, samplesheet in ctg-log
###  -> Write pipeline output
### * nf-Pipeline writes qc to shared/ctg-qc/sc-rna-10x

# Initialize variables
runfolder=$(pwd)
run=$(basename $runfolder)
demux="ON"
resume='n'

# usage message
usage() {

    echo ""
    echo "Usage: sc-rna-10x [ -i META_ID ] [ -s SAMPLESHEET ] [ -r RESUME ] [ -c CUSTOM-GENOME ] [ -d DEMUX-OFF ] [ -h HELP ] "  1>&2
    echo ""
    echo ""
    echo "Optional arguments: "
    echo "META-ID    -i : Set 'meta-id' for runfolder (e.g. 210330-10x). Default: Takes date of runfolder (before first _) and adds sc-rna-10x as suffix "
    echo "SAMPLESHEET   -s : Set samplesheet used for run (Default: CTG_SampleSheet.csv) "
    echo "RESUME        -r : Set if to resume nf-pipeline"
    echo "CUSTOM-GENOME -c : Path to custom reference genome if needed. Skip if human/mouse defined in samplesheet "
    echo "DEMUX-OFF     -d : Set flag to skip mkfastq (then fastq must be in FQDIR) "
    echo "HELP          -h : print help message"

}

exit_abnormal() {
    usage
    exit 1
}

# Read and control input arguments
while getopts i:s:c:dh opt; do
    case $opt in
	i) id=$OPTARG
	    ;;
	s) sheet=$OPTARG
	    ;;
	r) resume="y"
	    ;;
	c) custom_genome=$OPTARG
	    ;;
	d) demux="OFF"
	    ;;
	h) exit_abnormal
	    ;;
	\?) echo "> Error: Invalid option -$OPTARG" >&2
	    exit_abnormal ;;
	:) echo "> Error: -${OPTARG} requires an argument: -i needs project-id and -s need samplesheet name! "
	    exit_abnormal ;;
    esac
done

## Check arguments
shift "$(( OPTIND -1 ))"

if [ -z $id ]; then
    echo ""; echo "> WARNING! No meta-ID specified"
    metaid=$(echo $run | cut -f1 -d"_")
    id="${metaid}-sc-rna-10x"
    echo "- Using: '${id}'"
    if [ -f /projects/fs1/shared/ctg-projects/sc-rna-10x/$id ]; then
	echo "> Error: $id has been used before. "
	echo "/projects/fs1/shared/ctg-projects/sc-rna-10x/$id already exists. Please chose another meta-id."
	exit_abnormal
    fi
fi
if [ -z $sheet ]; then
    echo ""; echo "> WARNING! No samplesheet specified"
    sheet="CTG_SampleSheet.csv"
    echo "- Using 'CTG_SampleSheet.csv'"
    if [ ! -f $sheet ]; then
	echo "> Error: CTG_SampleSheet.csv does not exist (in current dir)"
	echo "- Please specify correct samplesheet, or create a CTG_SampleSheet.csv in current runfolder"
	exit_abnormal
    fi
fi

##############
# Print info #
##############
echo ""
echo "> The following arguments are entered:"
echo "ID               : $id"
echo "Sheet            : $sheet";
if [ -z $custom_genome ]; then
    echo "Custom Genome    : NONE "
else
    echo "Custom Genome    : $custom_genome "
fi
if [ $demux == "ON" ]; then
    echo "Demux            : YES "
else
    echo "Demux            : NO "
fi
if [ $resume == "y" ]; then
    echo "Resume           : YES "
else
    echo "Resume           : NO "
fi

echo ""
echo "> Running on runfolder: "
echo "Runfolder       : $runfolder "
echo ""

echo ""
echo "> Project log folder: "
echo "Logfolder       : /projects/fs1/shared/ctg-projects/sc-rna-10x/$id"
echo ""

echo ""
echo "> Output folder : "
echo "Output          : /projects/fs1/nas-sync/ctg-deliver/sc-rna-10x/$id"
echo ""

# Prompt user to approve running in current directory and input
read -p "> WARNING: Can only be run from within runfolder! Are you in runfolder in which you want run? (And are input described above correct) ?  (y/n)  ... " prompt
if [[ $prompt != "y" ]]
then
    echo "> Exiting: Go to runfolder!"
    exit 0
fi

################
# Set up files #
################

# Creating project dir for logging pipeline
projdir="/projects/fs1/shared/ctg-projects/sc-rna-10x/$id/"
mkdir -p $projdir

# Copy nextflow script and config to project folder
nf_pipe="/projects/fs1/shared/ctg-pipelines/ctg-sc-rna-10x/pipe-sc-rna-10x-multiproj/pipe-sc-rna-10x-multiproj.v1.nf"
nf_config="/projects/fs1/shared/ctg-pipelines/ctg-sc-rna-10x/pipe-sc-rna-10x-multiproj/nextflow.config"
cp $nf_pipe $projdir
cp $nf_config $projdir

# Edit config file
proj_conf=$projdir/nextflow.config
sed "s/xmetaidx/$id/g" $proj_conf > tmp.txt; mv tmp.txt $proj_conf
sed "s/xrunfolderx/$run/g" $proj_conf > tmp.txt; mv tmp.txt $proj_conf
if [ ! -z $custom_genome ]; then
    sed "s/xcustomgenomex/$custom_genome/g" $proj_conf > tmp.txt; mv tmp.txt $proj_conf
fi
if [ $demux == "ON" ];then
    sed "s/xdemuxx/y/g" $proj_conf > tmp.txt; mv tmp.txt $proj_conf
else
    sed "s/xdemuxx/n/g" $proj_conf > tmp.txt; mv tmp.txt $proj_conf
fi

# Copy edited config to runfolder
cp $proj_conf $runfolder
# Copy samplesheet to project folder (ctg-projects..)
cp $sheet $projdir

#####################
# Start nf-pipeline #
#####################
if [ $resume == "y" ]; then
    nohup nextflow run $nf_pipe -resume > nf.log.sc-rna-10x &
else
    nohup nextflow run $nf_pipe > nf.log.sc-rna-10x &
fi
echo ""; echo ""
echo "#################################"
echo "# sc-rna-10x pipeline submitted #"
echo "#################################"
echo ""; echo "";

##################
# ctg-interop-qc #
##################
# Will do interop analysis on runfolder, and compile multiqc report.
# This is written to
# 1. runfolder/ctg-interop
# 2. ctg-qc/interop
ctg-interop-qc
