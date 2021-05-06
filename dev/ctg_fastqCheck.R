#!/usr/bin/env Rscript

##  R Script that checks whether the output from demux is ok or not.
##  Uses fi;e names supplied in the proveded sample sheet


## ==========================================
##  LOAD PACKAGES
## ==========================================
library(optparse)
library(dplyr)



## ==========================================
##  PARSE COMMAND-LINE PARAMETERS
## ==========================================

option_list <- list(
  make_option(c("-i", "--sample_sheet"  ), type="character" , default=NULL   , metavar="path"   , help="Sample Sheet, colData, where samples are rows and columns are sample annotations."    ),
  make_option(c("-p", "--paired"        ), type="logical"   , default=TRUE   , metavar="boolean", help="If pariend end or not."                                            ),
)


opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt        <- parse_args(opt_parser)

if (is.null(opt$sample_sheet)){
  print_help(opt_parser)
  stop("Please provide a --sample_sheet file.", call.=FALSE)
}

# check strandness
paired.types <- c(TRUE, FALSE)
if(!opt$paired %in% paired.types) stop("--paired must be one of: ", paste(paired.types, collapse=", "))


## ========================================== ##
## PROCESS COUNTS FILE  - IEM STYLE
## ========================================== ##
## comment table is read and hashes are commented, hashes used when featureCounts produce output

# read data sheet as data frame (from below Data section header)
data_table <- read.delim(file=opt$sample_sheet, header = TRUE, sep = ",")
