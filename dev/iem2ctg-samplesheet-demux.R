#!/usr/bin/env Rscript

## ========================================== ##
##    Sctipt usage
## ========================================== ##

# Input:    An Illumina IEM style .csv sample sheet as input. "-s", "--sample_sheet" AND a project ID ("-i", "--projectid")

# Output:   Illiumina IEM style csv. 'ctg-demux' sample sheet. "-o", "--output_ctg_sheet" AND a

# Output:   CTG style sample sheet for delivery and nextflow processes. Will include file names for e.g. fastq and bam.
#    Columns:
#       - paired: TRUE/FALSE
#       - fastq_1


# 1: Check and substitute for illegal characters
# 2: Forces Sample_Name to Sample_ID
# 3: Forces Sample_Project same as input flag '-i --projectid'.

# Optional:
# 1: Check if adapter sequence is present in Settings. This should ideally be updated and checked by other means
# 2: Add Lana column (1 or 2). This e.g. if lane divider is used.

## Develiopment possibilities / known issues
## 1: cross check index sequences
## 2: cross check start read X from. e.g. Read2StartFrom, 4


# ========================================== ##
## LOAD PACKAGES                             ##
## ========================================== ##
library(optparse)
library(dplyr)


## ========================================== ##
## PARSE COMMAND-LINE PARAMETERS              ##
## ========================================== ##

option_list <- list(
  ## global params
  make_option(c("-i", "--projectid"     ), type="character" , default=NULL   , metavar="character"   , help="Project id. Required. NOTE: This will overwrite Sample_Project column in sample sheet"    ),
  make_option(c("-s", "--sample_sheet"  ), type="character" , default=NULL   , metavar="path"   , help="Sample Sheet, colData, where samples are rows and columns are sample annotations."    ),
  ## demux specific params
  make_option(c("-d", "--output_demux_sheet"  ), type="character"   , default=NULL  , metavar="path", help="Name of demux sample sheet"                          ),
  make_option(c("-l", "--force_lane"          ), type="numeric"     , default=0  , metavar="numeric", help="If to force lane, this is mostly only used if lane divider is used. Default is 0 which means no action. If '1' or '2' an exta column is added in the demux sample sheet."),
  make_option(c("-a", "--require_adapter"), type="logical"  , default=TRUE  , metavar="boolean", help="If to check for Illumina adapter. Will only check if there is a character vector here or not"),
  ## ctg sample sheet specific params
  make_option(c("-f", "--fastq_path"  ), type="character" , default=NULL   , metavar="path"   , help="Full path to fastq files."    ),
  make_option(c("-b", "--bam_path"  ), type="character" , default=NULL   , metavar="path"   , help="Optional. If specified, used to define full path to bam files."    ),
  make_option(c("-o", "--output_ctg_sheet"    ), type="character"   , default=NULL  , metavar="path", help="Name of output sample sheet"                          ),
  make_option(c("-p", "--paired"        ), type="logical"   , default=TRUE   , metavar="boolean", help="If pariend end or not. Used only for ctg style sample sheet, not for demux."                                            )

)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt        <- parse_args(opt_parser)

if (is.null(opt$projectid)){
  print_help(opt_parser)
  stop("Please provide p --projectid.", call.=FALSE)
}
if (is.null(opt$sample_sheet)){
  print_help(opt_parser)
  stop("Please provide a --sample_sheet file.", call.=FALSE)
}
if (is.null(opt$output_demux_sheet)){
  print_help(opt_parser)
  stop("Please provide --output_demux_sheet file.", call.=FALSE)
}
if (is.null(opt$output_ctg_sheet)){
  print_help(opt_parser)
  stop("Please provide --output_ctg_sheet file.", call.=FALSE)
}
if (is.null(opt$fastq_path)){
  print_help(opt_parser)
  stop("Please provide --fastq_path ", call.=FALSE)
}
# if (is.null(opt$bam_path)){
#   print_help(opt_parser)
#   stop("Please provide a --bam_path file.", call.=FALSE)
# }

if(!opt$force_lane %in% c(0,1,2)){
  print_help(opt_parser)
  stop("--force_lane must be set to '0', '1' or '2'. Default is '0'.", call.=FALSE)
}


cat("\n ... Rscript: Parsing sample sheet. \n\n")

# ========================================== ##
## PROCESS - IEM STYLE SAMPLE SHEET
## ========================================== ##
all.lines <- scan(file = opt$sample_sheet, what = "character", sep="\n", nmax = 250)

## Gsub all regibon specific CHARACTERS

## Start replacing non ascii characters if present
all.lines=iconv(all.lines, "latin1", "ASCII", sub="")

## check nlines of document, max is set to 500
if(length(all.lines)>=500) stop("Number of rows in sample sheet exceeds what is read by the scan function ('nmax' set to 500)")

## Start replacing all foreign caharacters if present


## Check that IEM sheet contains correct section
iem.headers <- list(Header="[[]Header[]]", Reads="[[]Reads[]]",Settings="[[]Settings[]]",Data="[[]Data[]]")
iem.index <- sapply(iem.headers, function(x) grep(x, all.lines))
if(!any(unlist(lapply(iem.index, length)==1))){
  stop("Sample Sheet Not IEM format - check Illumina IEM sections, must include - [Header], [Reads], [Settings], [Data] ...\n")
}
iem.index <- sapply(iem.headers, function(x) grep(x, all.lines))

## Check Adapter if adapter sequence slot is NA or not
# Note, will not check if different Adapter reads for forward or reverse
if(opt$require_adapter){
  settings.section <- read.delim(file=opt$sample_sheet, header = F, sep = ",",
                                 skip = iem.index["Settings"], nrows = iem.index["Data"]-iem.index["Settings"]-1)
  u <- grep("Adapter", settings.section[,1])
  if(!length(u)) stop("Cant find 'Adapter' row in [Settings] section ")
  if(is.na(settings.section[u[1],1])) stop("No Adapter sequence detected. Define adapter or set --require_adapter FALSE")
}


# read data sheet as data frame (from below Data section header)
data_table <- read.delim(file=opt$sample_sheet, header = TRUE, sep = ",",
                           skip = iem.index["Data"])

## Start replacing non ascii characters if present
# if(!identical(data_table, iconv(data_table, "latin1", "ASCII", sub=""))) stop("non ascii characters in data section!")


# Check Data section.
## 1. Columns must include Sample_ID,Sample_Project columns
### NOTE THAT Samlple_Name is not needed. Sample_Name will be forced to Sample_ID
### NOTE that Sample_Project is not longer needed/used. This will be forced from opt projectid
data_table$Sample_Project <- opt$projectid
data_table$Sample_Name <- data_table$Sample_ID

required.columns <- c("Sample_ID","Sample_Project", "index")
if(!all (required.columns %in% colnames(data_table))){
  stop("not all required [Data] section header names are present")
}

## Allow only Sample_ID - Force Sample_Name to Sample_ID
## Force Sample_Project to input project id



# 2a. Check for non permitted special characters
## The Description column, if present, is treated a bit more gentle. Only type_a specials are disallowed here amd spaces, dashes, parenthese, dots are allowed
specials <- list(
  type_a=c("[*]","[/]","[+]","[']","[`]","[?]","[=]"),
  type_b=c("[:]","[-]","[(]","[)]","[.]"," ")
  )

# matrix(data_table)
special.flag <- F
data_table.out <- data_table

## Check and replace the specific columnns in data section for stricter chacacter rules.
cat("\n\n ... Checking data columns and replacing special characters. ")
special.flag <- F # used to keep track if special characters are present in any column

test_columns = unique(c(required.columns, colnames(data_table.out)[grepl("index|Index",colnames(data_table.out))]))

for(i in 1:ncol(data_table.out)){
  # test_specials <- unlist(specials)
  # if(colnames(data_table.out)[i]=="Description") test_specials <- specials$type_a
  ## the test_columns will be tested for the most strict character list, i.e type_a
  if(colnames(data_table.out)[i] %in% test_columns){
    test_specials <-  unlist(specials)
  }else{
    test_specials <- specials$type_a
    }

  for(ii in 1:length(test_specials)){
    u<-grepl(test_specials[ii], data_table.out[,i])
    if(any(u)){
      special.flag <- T
      cat("\n ... ... Warning: column  '", colnames(data_table.out)[i], "'  contains: ", ifelse(test_specials[ii]==" ", "white space. Replacing with ''", test_specials[ii]))
      ## replace these with underscore in out-file
      data_table.out[,i] <- gsub(test_specials[ii], "", data_table.out[,i])
      }
  }
}


## Check if Sample_ID contain duplicates
u <- duplicated(data_table$Sample_ID)
if(any(u)) stop("Duplicated Sample_ID present: ", data_table$Sample_ID[u])
cat("\n\n ... OK, no duplicate Sample_IDs found.")


## add force_lane if specified
if(opt$force_lane != 0) data_table.out$Lane <- opt$force_lane
if(opt$force_lane != 0) data_table.out <- dplyr::select(.data = data_table.out, c(Lane, everything())) # strandness

## save output iem file for demux
out.file.name <- opt$output_demux_sheet

# if(grepl("[.]csv", opt$sample_sheet)) out.file.name <- gsub(pattern = "[.]csv", replacement = "_iem_corrected.csv", x = opt$sample_sheet)
# if(file.exists(out.file.name)) stop("Out file already exists: ",out.file.name)
cat("\n\n ... writing checked ctg-demux sample sheet file to output: \n ",out.file.name)

cat(all.lines[1:iem.index["Data"]], file = out.file.name, sep = "\n", append = F)
cat(colnames(data_table.out), file = out.file.name, sep = ",", append = T)
cat("\n", file = out.file.name, sep = ",", append = T)
write.table(x = data_table.out, file = out.file.name, sep = ",", append = T, quote = F, row.names = F, col.names = F, na="", )



cat("\n ... #   Generate ctg style sample sheet "            )
cat("\n ... " )

# nf core style sample sheets must have following columns
## Strandneess and paired should preferrably be obtained from IEM, but as for now supplied as input variables
## group replicate fastq_1 fastq_2 strandness
## group: try to get group from sample sheet, if not available, then use Sample_Project
## replicate: use if presernt in sample sheet, else create one using Sample_Project
## For each Sample_Name processed, a suffiux is added as:
## _nS_R1_001.fastq.qz. The S-counter  will just be numeric order 1:nrow
## The L00X wil be differnent if multiple libraries, as for now 001 is the only one supported


data_table.nfcore <- data_table.out

## create group column
if(!("group" %in% colnames(data_table.nfcore))){
  cat("\n ... creating 'group' column from ¨Sample_Project¨")
  data_table.nfcore$group <- data_table.nfcore$Sample_Project
}
# generate replicate column
if(!("replicate" %in% colnames(data_table.nfcore))){
  cat("\n ... creating 'replicate' column from 'group'")
  u <- sapply(unique(data_table.nfcore$group), function(x) length(which(data_table.nfcore$group==x)))
  rep.vec <- data_table.nfcore$group
  for(i in 1:length(u)){
    rep.vec[which(rep.vec==names(u)[i])] <- 1:u[i]
  }
  data_table.nfcore$replicate <- as.integer(rep.vec)
}
# generate batch column
if(!("batch" %in% colnames(data_table.nfcore))){
  cat("\n ... creating 'batch' column from ¨Sample_Project¨")
  data_table.nfcore$batch <- data_table.nfcore$Sample_Project
}


# ========================================== ##
## PROCESS - GENERATE CTG STYLE SAMPLE SHEET
## ========================================== ##
## Create sample-sheet-ctg-<project_id>.csv
## For the rna seq pipeline, the sample sheet should include -
## bcl2Fastq will name fastq files: "FASTQ files are named with i) the sample name 2) and number, 3) the flow cell lane, 4) and read. 5) After that # 001.
## The file extension is *.fastq.gz. Forexample:<samplename>_S1_L001_R1_001.fastq.gz"
## bcl2fastq path:
##     path is supplied to bcl2fastq with the -o flag, in our case something like Fastq_Raw, defined in
## Note: if Sample_ID != Sample_Name each fastq

    # All input (and checked) columns
    # fastq_1: Name of fastq R1 file. This according to default bcl2fastq settings.

fastq.sufifx <- "_001.fastq.qz"
bam.suffix <- "_Aligned.sortedByCoord.out.bam"

data_table.nfcore$fastq_1 <- paste(data_table.nfcore$Sample_Name, "_S",1:nrow(data_table.nfcore), "_R1", fastq.sufifx, sep="")
## Note that when path is defined it is assumed that Sample_ID == Sample_Name (blc2fastq2 generate flat output). If name != ID folder is created for Name
data_table.nfcore$fastq_1_path <- file.path(opt$fastq_path,data_table.nfcore$fastq_1)

if(opt$paired) data_table.nfcore$fastq_2 <- paste(data_table.nfcore$Sample_Name, "_S",1:nrow(data_table.nfcore), "_R2", fastq.sufifx, sep="")
if(opt$paired) data_table.nfcore$fastq_2_path <- file.path(opt$fastq_path,data_table.nfcore$fastq_2)


# if(opt$paired) data_table.nfcore$fastq_2_id_path <- paste(data_table.nfcore$Sample_Name, "_S",1:nrow(data_table.nfcore), "_R2", fastq.sufifx, sep="")

## if bam - supply bam path
if(!is.null(opt$bam_path)) data_table.nfcore$bam <- paste0(data_table.nfcore$Sample_Name, opt$bam_suffix)
if(!is.null(opt$bam_path)) data_table.nfcore$bam_path <- file.path(opt$bam_path, data_table.nfcore$bam)

data_table.nfcore$counts_matrix_id <- paste0(data_table.nfcore$Sample_Name, opt$count_mat_siuffix)
data_table.nfcore$paired <- opt$paired
#data_table.nfcore$strandness <- opt$strandness
data_table.nfcore$counts_matrix_id <- paste0(data_table.nfcore$Sample_Name, opt$count_mat_siuffix)

## re-order sheet
if(opt$paired) data_table.nfcore <- dplyr::select(.data = data_table.nfcore, c(Sample_ID, Sample_Name, group, replicate, paired, fastq_1, fastq_2, everything())) # strandness
if(!opt$paired) data_table.nfcore <- dplyr::select(.data = data_table.nfcore, c(Sample_ID, Sample_Name, group, replicate, paired, fastq_1, everything())) # strandness


# If to save updated sample sheet to use in next run/itteratve process
#if(opt$output_ctg){

  # out.file.name <- paste(opt$sample_sheet, "_ctg.csv", sep="")
 out.file.name <- opt$output_ctg_sheet


  # if(grepl("[.]csv", opt$sample_sheet)) out.file.name <- gsub(pattern = "[.]csv", replacement = "_ctg.csv", x = opt$sample_sheet)
  cat("\n\n ... writing ctg style sample sheet with expected filenames: \n ", out.file.name, "\n\n")
  write.table(x = data_table.nfcore, file = out.file.name, sep = ",", append = F, quote = F, row.names = F, col.names = T, na="", )
