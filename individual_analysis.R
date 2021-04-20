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

setwd("~/Documents/InnoVar/QTLs/Leaf_rust/SRR3929/")

samples<-as.data.frame(dir())
colnames(samples)<-"samples"
files <- file.path(samples$samples, "abundance.h5")
names(files) <- paste0(samples$samples)
all(file.exists(files))
txi.kallisto.tsv <- tximport(files, type = "kallisto", countsFromAbundance = "scaledTPM", ignoreAfterBar = TRUE, txIn=TRUE, txOut=TRUE)
colData<-read.delim("metadata.txt")
colData$Factor<-paste(colData$Genotype, colData$HPI)
colData<-colData[order(colData$Run),]

# colData is metadata with factor column. Make dds object with deseq
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ HPI + Genotype)
# transform using variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

# generate PC1 and PC2 for visualisation
pcaData <- plotPCA(vsd, intgroup=c("HPI", "Genotype", "replicate"), returnData=TRUE)

# plot PC1 vs PC2
ggplot(pcaData, aes(x=PC1, y=PC2)) + geom_point(aes(colour=HPI, shape=Genotype), size=4, alpha=0.7) +
  theme_classic() +
  theme(text = element_text(size=20, colour="black"),
        axis.text.x = element_text(colour="black"))
# Import kallisto files with txi 
# ==================================================================================


dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, colData, ~ Factor)

summary(dds)
dds<-DESeq(dds)

RES1<-as.data.frame(results(dds, contrast=c("Factor","WL711 12hpi", "WL711 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES2<-as.data.frame(results(dds, contrast=c("Factor","WL711 24hpi",	"WL711 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES3<-as.data.frame(results(dds, contrast=c("Factor","WL711 48hpi",	"WL711 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES4<-as.data.frame(results(dds, contrast=c("Factor","WL711 72hpi",	"WL711 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES5<-as.data.frame(results(dds, contrast=c("Factor","WL711_Lr57 12hpi",	"WL711_Lr57 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES6<-as.data.frame(results(dds, contrast=c("Factor","WL711_Lr57 24hpi",	"WL711_Lr57 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES7<-as.data.frame(results(dds, contrast=c("Factor","WL711_Lr57 48hpi",	"WL711_Lr57 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))
RES8<-as.data.frame(results(dds, contrast=c("Factor","WL711_Lr57 72hpi",	"WL711_Lr57 0hpi"), pAdjustMethod='BH', alpha=0.05, format='DataFrame', tidy=T))

RES1$Factor<-"WL711_12hpi"
RES2$Factor<-"WL711_24hpi"
RES3$Factor<-"WL711_48hpi"
RES4$Factor<-"WL711_72hpi"
RES5$Factor<-"WL711_Lr57_12hpi"
RES6$Factor<-"WL711_Lr57_24hpi"
RES7$Factor<-"WL711_Lr57_48hpi"
RES8$Factor<-"WL711_Lr57_72hpi"
RES<-rbind(RES1,
           RES2,
           RES3,
           RES4,
           RES5,
           RES6,
           RES7,
           RES8)
sig<-RES[(RES$padj<=0.05),]
sig<-na.omit(sig)
table(sig$Factor)

sig<-sig%>% 
  mutate(Genotype=case_when(
    Factor == "WL711_12hpi" ~ "WL711",
    Factor == "WL711_24hpi" ~ "WL711",
    Factor == "WL711_48hpi"~ "WL711",
    Factor == "WL711_72hpi" ~ "WL711",
    Factor =="WL711_Lr57_12hpi" ~ "WL711_Lr57",
    Factor =="WL711_Lr57_24hpi" ~ "WL711_Lr57",
    Factor =="WL711_Lr57_48hpi" ~ "WL711_Lr57",
    Factor =="WL711_Lr57_72hpi" ~ "WL711_Lr57"
  ))

WL711<-sig[(sig$Genotype=="WL711"),]
WL711_Lr57<-sig[(sig$Genotype=="WL711_Lr57"),]
WL711_genes<-as.data.frame(unique(WL711$row))
WL711_Lr57_genes<-as.data.frame(unique(WL711_Lr57$row))

WL711_genes$Condition<-1
WL711_Lr57_genes$Condition<-1
colnames(WL711_genes)<-c("Gene", "Condition")
colnames(WL711_Lr57_genes)<-c("Gene", "Condition")

match<-full_join(WL711_genes, WL711_Lr57_genes, by="Gene")
match[is.na(match)] <- 0
colnames(match)<-c("Gene", "WL711", "WL711_Lr57")
match$cross<-match$WL711 + match$WL711_Lr57

match<-match%>% 
  mutate(Both=case_when(
    cross == 2 ~ 1,
    cross == 1 ~ 0,
  ))

colSums(match)
sum(match$WL711)
sum(match$WL711_Lr57)
sum(match$Both)

draw.pairwise.venn(area1 = sum(match$WL711),                        # Create pairwise venn diagram
                   area2 = sum(match$WL711_Lr57),
                   cross.area = sum(match$Both),
                   fill = c("pink", "orange"),
                   category = c("WL711", "WL711 + Lr57"),
                   scaled=F, cat.pos = 0)

setwd("~/Documents/InnoVar/QTLs/Leaf_rust/")
write.csv(sig, file="~/Documents/InnoVar/QTLs/Leaf_rust/SRR3929/all_sig.csv")



degs<-as.character(unique(all_filtered_sig$row))
degs_all_data<-subset(all_filtered, all_filtered$row %in% degs)

degs_all_data_M<-degs_all_data[c(1, 3, 8)]
degs_all_data_M<-spread(degs_all_data_M, key="comparison", value="log2FoldChange", fill=0)

write.csv(all_filtered, "~/Documents/S_L_DH/DEtesting/data/DE_tests/separate_filtered.csv")
write.csv(degs_all_data_M, file="DEGs_all_data_matrix.csv")
write.csv(degs_all_data, file="DEGs_all_data.csv")



