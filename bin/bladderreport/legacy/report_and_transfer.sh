#!/bin/bash
# Run bladderreport_ctg_anonymous.Rmd and tranfer files with lfsscp.
# Created by PP 20191104.

if [[ $# != 1 ]] ; then
    echo 'How to use this script:'
    echo 'report_and_transfer.sh <ctg sample identifer>'
    exit 0
fi

sampleid=$1

echo "Fetching paths for CTG identifier: ${sampleid}"

rsemin_path=`locate ${sampleid} | grep "/data/bnf/premap/rnaseq/" | grep ".rsem$"`

if [ -z "${rsemin_path}" ]
then
    echo "${rsemin_path} is not valid, check if ${sampleid} exists!"
    exit 0
fi

starqc_path=`locate ${sampleid} | grep "/data/bnf/tmp/rnaseq/" | grep "Log.final.out$"`

if [ -z "${starqc_path}" ]
then
    echo "${starqc_path} is not valid, check if ${sampleid} exists!"
    exit 0
fi


clarity_ID=`locate 19KFU0140-CTG | \
grep "/data/bnf/batch/" | \
grep ".sh$" | \
perl -na -e 'chomp; $sh_extract=\`grep clarity_id $_\`; $quote=chr(39); $sh_extract=~m/clarity_id=$quote(.*)$quote\),/; $clarity_id=$1; print "$clarity_id\n";'`

if [ -z "${clarity_ID}" ]
then
    echo "${clarity_ID} is not valid, check if ${sampleid} exists!"
    exit 0
fi

echo "Path rsem_in: ${rsemin_path}"
echo "Path star_qc: ${starqc_path}"
echo "Name Clarity: ${clarity_ID}"

echo "####################################"
echo "Creating anonymous report..."

# function to scp results to file server lfs603.srv.lu.se
lfsscp ()
{
    if [[ ( -z $1 ) || ( $1 == "-h" ) ]]; then
        echo 'Syntax: lfsscp [-r] <file> [<destination>]';
        echo 'Destination can be a directory or a file name depending on';
        echo 'the prescence of a trailing slash or not';
        echo '   -r copies directories recursively.';
        echo;
        echo '(The primary large disk is located at /srv/data/)';
    else
        if [[ $1 == "-r" ]]; then
            if [[ $# -gt 3 ]]; then
                echo 'The maximum number of args are 2. If you like to use a';
                echo 'wild card then quote the expression:';
                echo "lfsscp -r 'foo*.json' /some/destination/";
                exit;
            fi;
            recursive=$1;
            copyfrom=$2;
            copyto=$3;
        else
            if [[ $# -gt 2 ]]; then
                echo 'The maximum number of args are 2. If you like to use a';
                echo 'wild card then quote the expression:';
                echo "lfsscp 'foo*.json' /some/destination/";
                exit;
            fi;
            recursive="";
            copyfrom=$1;
            copyto=$2;
        fi;
        scp -o ProxyCommand="ssh -W %h:%p rs-fs1.lunarc.lu.se" -P 22022 $recursive $copyfrom lfs603.srv.lu.se:$copyto;
    fi
}

# run R markdown to generate report
/usr/bin/Rscript -e "library(rmarkdown,\
lib='/home/petter/R/x86_64-pc-linux-gnu-library/3.4');\
rmarkdown::render('/data/bnf/scripts/bladderreport/bladderreport_ctg_anonymous.Rmd',\
params=list(\
sampleid='${sampleid}',\
rsem_in='${rsemin_path}',\
star_qc='${starqc_path}',\
clarity_id='${clarity_ID}'\
),\
output_file='/data/bnf/postmap/rnaseq/${sampleid}.STAR.bladderreport_anonymous.html')"

echo "Done"

echo "Transfer report to fileserver..."

lfsscp /data/bnf/postmap/rnaseq/${sampleid}.STAR.bladderreport_anonymous.html /srv/data/pontus_e/anonymous_bladder_reports

echo "Done"
