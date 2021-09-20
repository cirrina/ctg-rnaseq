#!/usr/bin/env Rscript
#source('http://bioconductor.org/biocLite.R') 

list.of.packages <- c("ggplot2", "cowplot", "reshape", "knitr","kableExtra","dplyr","ranger","jsonlite")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]


if(length(new.packages)>1){
  cat("Installing", new.packages)
  install.packages(new.packages)
}else{
  cat("all packages already installed\n")
}