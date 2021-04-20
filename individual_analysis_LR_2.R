# de testing for individual datasets that were filtered before

library(WGCNA)
library(tximport)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(gplots)
library(tidyr)
library(Hmisc)
library(corrplot)
library(VennDiagram)

setwd("~/Documents/InnoVar/QTLs/Leaf_rust/SRR12998/")

samples<-as.data.frame(dir())
colnames(samples)<-"samples"
files <- file.path(samples$samples, "abundance.h5")
names(files) <- paste0(samples$samples)
all(file.exists(files))
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
colData<-read.delim("metadata.txt")
colData$Factor<-paste(colData$Genotype, colData$Isolate)
colData<-colData[order(colData$Run),]

# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Isolate + Genotype)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("Isolate", "Genotype"), returnData=TRUE)

# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(colour=Isolate, shape=Genotype), size=4, alpha=0.7) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black"))
# Import kallisto files with txi 
# ==================================================================================


dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Factor)

summary(dds)
dds<-DESeq(dds)

RES1<-as.data.frame(results(dds, contrast=c("Factor", "ThatcherLr14a Lr_96209",  "ThatcherLr14a Lr_mock"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES2<-as.data.frame(results(dds, contrast=c("Factor", "Thatcher Lr_96209", "Thatcher Lr_mock"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES3<-as.data.frame(results(dds, contrast=c("Factor", "ThatcherLr14a Lr_95037", "ThatcherLr14a Lr_mock"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES4<-as.data.frame(results(dds, contrast=c("Factor", "Thatcher Lr_95037", "Thatcher Lr_mock"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))

RES1$Factor<-"ThatcherLr14a Lr_96209"
RES2$Factor<-"Thatcher Lr_96209"
RES3$Factor<-"ThatcherLr14a Lr_95037"
RES4$Factor<-"Thatcher Lr_95037"
RES<-rbind(RES1,
           RES2,
           RES3,
           RES4)
sig<-RES[(RES$padj<=0.05),]
sig$abs<-abs(sig$log2FoldChange)
sig<-sig[(sig$abs >=1),]
sig<-na.omit(sig)
table(sig$Factor)

sig<-sig%>% 
  mutate(Genotype=case_when(
Factor == "ThatcherLr14a Lr_96209" ~ "Thatcher Lr14a",
Factor == "Thatcher Lr_96209" ~ "Thatcher",
Factor == "ThatcherLr14a Lr_95037" ~ "Thatcher Lr14a",
Factor == "Thatcher Lr_95037" ~ "Thatcher"
  ))

sig<-sig%>% 
  mutate(Isolate=case_when(
    Factor == "ThatcherLr14a Lr_96209" ~ "96209",
    Factor == "Thatcher Lr_96209" ~ "96209",
    Factor == "ThatcherLr14a Lr_95037" ~ "95037",
    Factor == "Thatcher Lr_95037" ~ "95037"
  ))

sig$conditon<-1
sig2<-sig[c(1, 8, 12)]
sig3<-spread(sig2, key="Factor", value="conditon", fill=0)
row.names(sig3)<-sig3$row
sig3$row<-NULL

likes <- function(animals) {
  ppl <- sig3
  names(ppl) <- colnames(sig3)
  for (i in 1:length(animals)) {
    ppl <- subset(ppl, ppl[animals[i]] == T)
  }
  nrow(ppl)
}

# How many people like dogs?
likes(colnames(sig3))
plotAnimals <- function(a, ...) {
  grid.newpage()
  if (length(a) == 1) {
    out <- draw.single.venn(likes(a), ...)
  }
  if (length(a) == 2) {
    out <- draw.pairwise.venn(likes(a[1]), likes(a[2]), likes(a[1:2]), ...)
  }
  if (length(a) == 3) {
    out <- draw.triple.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[1:2]), 
                            likes(a[2:3]), likes(a[c(1, 3)]), likes(a), ...)
  }
  if (length(a) == 4) {
    out <- draw.quad.venn(likes(a[1]), likes(a[2]), likes(a[3]), likes(a[4]), 
                          likes(a[1:2]), likes(a[c(1, 3)]), likes(a[c(1, 4)]), likes(a[2:3]), 
                          likes(a[c(2, 4)]), likes(a[3:4]), likes(a[1:3]), likes(a[c(1, 2, 
                                                                                     4)]), likes(a[c(1, 3, 4)]), likes(a[2:4]), likes(a), ...)
  }
  if (!exists("out")) 
    out <- "Oops"
  return(out)
}


plotAnimals(colnames(sig3), category = c("Susceptible + Virulent", " Susceptible + Avirluent", 
                                         "Resistant + Virulent", "Resistant + Avirulent"), 
          fill = c("skyblue", "pink1", "mediumorchid", "beige"),
            scaled=F)



setwd("~/Documents/InnoVar/QTLs/Leaf_rust/")

write.csv(sig, "~/Documents/InnoVar/QTLs/Leaf_rust/SRR12998/all_sig.csv")


degs<-as.character(unique(all_filtered_sig$row))
degs_all_data<-subset(all_filtered, all_filtered$row %in% degs)

degs_all_data_M<-degs_all_data[c(1, 3, 8)]
degs_all_data_M<-spread(degs_all_data_M, key="comparison", value="log2FoldChange", fill=0)

write.csv(all_filtered, "~/Documents/S_L_DH/DEtesting/data/DE_tests/separate_filtered.csv")
write.csv(degs_all_data_M, file="DEGs_all_data_matrix.csv")
write.csv(degs_all_data, file="DEGs_all_data.csv")



