#!/usr/bin/env bash

usage() { echo "Usage: $0 [-r rsem_file] [-s sampleid] [-q starqc]
example: 
 bash $0 -r example_data/19KFU0006_0.rsem -s 19KFU0006_0.rsem -q example_data/Log.final.out
" 1>&2; exit 1; }

while getopts ":r:s:q:" o; do
    case "${o}" in
	r)
	    r=${OPTARG}
	    ;;
	s)
	    s=${OPTARG}
	    ;;
	q)
	    q=${OPTARG}
	    ;;
	
	*)
	    usage
	    ;;
    esac
done

shift $((OPTIND-1))

if [ -z "${r}" ] || [ -z "${s}" ]; then
    usage
fi

if [ "$(uname)" == "Darwin" ]; then
    export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc
fi

echo "Rscript -e 'library(rmarkdown); rmarkdown::render(\"./bladderreport_ctg.Rmd\", params=list(sampleid=\"$s\", rsem_in=\"$r\", star_qc=\"$q\"), \"html_document\")'" | bash


