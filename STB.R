# Check proximity of DEGs to QTLs for a particular stress

library(grid)
library(gridExtra)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

all_markers_positions_filtered <- read.csv("~/Documents/InnoVar/QTLs/QTL_data/all_markers_positions_filtered.csv", row.names=1)

disease="Fg" # change here to say which disease you want to analyse

# read in expression data for that disease
files<-list.files("~/Documents/InnoVar/QTLs/Expression_data/", pattern=disease)

for(file in files){
  study<-substr(file, 5, 13)
  df<-read.csv(paste("~/Documents/InnoVar/QTLs/Expression_data/", file, sep=""))
  assign(paste(study), df)
}

setwd(paste("~/Documents/InnoVar/QTLs/QTL_data/), disease)