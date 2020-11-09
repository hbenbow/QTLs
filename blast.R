# see what the genes are from arab, brachy, rice
library(grid)
library(gridExtra)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)

setwd("~/Documents/InnoVar/QTLs/")
# read in blast results 
arab<-read.table("wheat_arab.txt")
brachy<-read.table("wheat_brachy.txt")
rice<-read.table("wheat_rice.txt")


colnames(brachy)<-c(
  "qseqid",
  "sseqid",
  "pident",
  "length",
  "mismatch",	
  "gapopen",
  "qstart",
  "qend",
  "sstart",	 
  "send",
  "evalue",
  "bitscore")	
arab_genes <- read.delim("~/Documents/InnoVar/QTLs/arab_genes.txt", header=FALSE)
brachy_genes <- read.delim("~/Documents/InnoVar/QTLs/brachy_genes.txt", header=FALSE)
rice_genes <- read.delim("~/Documents/InnoVar/QTLs/rice_genes.txt", header=FALSE)

colnames(rice_genes)<-c("Gene", "Function")

arab2<-merge(arab, arab_genes, by.x="sseqid", by.y="Gene", all.x = T)
brachy2<-merge(brachy, brachy_genes, by.x="sseqid", by.y="Gene", all.x = T)
rice2<-merge(rice, rice_genes, by.x="sseqid", by.y="Gene", all.x = T)

GOIs<-c("TraesCS1A02G157900.1",
        "TraesCS1A02G403300.1",
        "TraesCS1B02G265600.1",
        "TraesCS2D02G082700.1",
        "TraesCS4A02G264400.1",
        "TraesCS4B02G043100.1",
        "TraesCS4B02G285000.1",
        "TraesCS5B02G396600.1",
        "TraesCS5D02G407600.1"
)

goisA<-subset(arab2, arab2$qseqid %in% GOIs) 
goisB<-subset(brachy2, brachy2$qseqid %in% GOIs) 
goisR<-subset(rice2, rice2$qseqid %in% GOIs) 
write.csv(goisA, file="GOIsA.csv")
write.csv(goisB, file="GOIsB.csv")
write.csv(goisR, file="GOIsR.csv")



