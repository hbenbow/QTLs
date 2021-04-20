# Check proximity of DEGs to QTLs for a particular stress

library(grid)
library(gridExtra)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(knitr)

# Read in all files
setwd("~/Documents/InnoVar/QTLs/QTL_data/")
all_markers_positions <- read.csv("~/Documents/InnoVar/QTLs/marker_positions/all_markers_positions.csv", row.names=1)
QTL_database <- read.csv("~/Documents/InnoVar/QTLs/QTL_data/LR/qtl_database.csv")
QTL_database$Chromosome<-paste("chr", QTL_database$Chromoome, sep="")

qtl_indb<-subset(all_markers_positions, all_markers_positions$Feature %in% QTL_database$Linked.markers)
write.csv(qtl_indb, file="~/Documents/InnoVar/QTLs/QTL_data/LR/available_markers.csv")

# ========================================================================================================================================
# this chunk will identify markers that have multiple hits on the same chromosome
list<-list()
qtl_indb<-subset(all_markers_positions, all_markers_positions$Feature %in% QTL_database$Linked.markers)
for(marker in qtl_indb$Feature){
  data<-qtl_indb[(qtl_indb$Feature==marker),]
  table<-spread(as.data.frame(table(data[,c(1,2)])), key=Start, value=Freq)
  row.names(table)<-table$Chromosome
  table$Chromosome<-NULL
  table<-as.matrix(table)
  table[table>0]<-1
  table<-as.data.frame(table)
  table$sum<-rowSums(table)
  table$marker<-paste(marker)
  test<-table[(table$sum>1),]
  list[[length(list)+1]]<-test$marker}

mismatches<-do.call(rbind.data.frame, list)
mismatches<-as.data.frame(unique(mismatches[,1]))
mismatches<-subset(all_markers_positions, all_markers_positions$Feature %in% mismatches[,1])
write.csv(mismatches, file="~/Documents/InnoVar/QTLs/QTL_data/LR/mismatches.csv")

# ========================================================================================================================================
# this section makes a new file that has no mismatches (i.e. markers that have multiple hits on the same chromosome)

all_markers_positions$ID<-paste(all_markers_positions$Chromosome, all_markers_positions$Start, all_markers_positions$Feature)
no_dups<-all_markers_positions[!(duplicated(all_markers_positions$ID)),]
no_dups$ID<-paste(no_dups$Chromosome, no_dups$Feature)
no_dups<-no_dups[!(duplicated(no_dups$ID)),]
all_markers_positions_filtered<-no_dups
write.csv(all_markers_positions_filtered, file="~/Documents/InnoVar/QTLs/QTL_data/LR/all_markers_positions_filtered.csv")

qtl_indb_filtered<-subset(all_markers_positions_filtered, all_markers_positions_filtered$Feature %in% QTL_database$Linked.markers)



exp_data<-read.csv("~/Documents/InnoVar/QTLs/Leaf_rust/all_sig_lr.csv")
wheat_bed<-read.csv("~/Documents/InnoVar/QTLs/QTL_data/wheat_all.csv")
colnames(wheat_bed)<-c( "Chromosome", "start", "end", "GeneID","Score", "strand")

exp_data2<-merge(exp_data, wheat_bed, by.x="Gene", by.y="GeneID")

qtls_clusters<-list()
all_qtls<-list()
qtls<-QTL_database$QTL
for(qtl in qtls){
  data2<-QTL_database[(QTL_database$QTL==qtl),]
  marker<-data2$Linked.markers
  for(i in marker){
    data<-data2[(data2$Linked.markers==i),]
    chr<-data$Chromosome
    position<-subset(all_markers_positions_filtered, all_markers_positions_filtered$Chromosome %in% chr)
    position<-subset(position, position$Feature %in% i)
    clusters<-subset(exp_data2, exp_data2$Chromosome %in% chr)
    len<-nrow(clusters)
    lenp<-nrow(position)
    if(lenp>0 && len>0){
      clusters$qtl<-paste(qtl)
      clusters$marker<-paste(i)
      start<-as.numeric(position[(position$Feature==i),2])
      clusters$marker_position<-paste(start)
      qtls_clusters[[length(qtls_clusters)+1]]<-clusters
    }
  }
  # clusters$marker_position<-as.numeric(paste(start))
  # clusters$difference<-abs(clusters$start-clusters$marker_position)/1000000
  # all_qtls[[length(all_qtls)+1]]<-clusters
}

qtls_clusters<-do.call(rbind.data.frame, qtls_clusters)
qtls_clusters$start<-as.numeric(qtls_clusters$start)
qtls_clusters$end<-as.numeric(qtls_clusters$end)
qtls_clusters$marker_position<-as.numeric(qtls_clusters$marker_position)
qtls_clusters<- qtls_clusters %>%
  mutate(within = case_when(marker_position < end & marker_position > start ~ "Within", marker_position > end & marker_position < start ~ "Within"))
qtls_clusters$Difference<-abs(as.numeric(qtls_clusters$start) - as.numeric(qtls_clusters$marker_position)) /1000000
write.csv(qtls_clusters, file="~/Documents/InnoVar/QTLs/QTL_data/LR/Lr_DEGs_QTLs.csv")

# ====================================================================================================================
# A histogram of distance between QTLs and FRGCs
q<-quantile(qtls_clusters$Difference, probs = c(0.01,0.05, 0.02, 0.05))
q
summary(qtls_clusters$Difference)


ggplot(qtls_clusters, aes(x=Difference)) + 
  geom_histogram(fill="grey60", colour="black") +
  theme_classic()+
  theme(text = element_text(size=20, colour="black")) +
  xlab("Distance between DEG and QTL markers (Mbp)") +
  geom_vline(aes(xintercept=q[2]), colour="firebrick", size=1)


# lets see those where the qtl marker is within a DEG
withins<-qtls_clusters[(qtls_clusters$within=="Within"),]
withins<-na.omit(withins)
# lets examine qtls that are <5 mbp from a FRGC

threshold<-q[2]
data<-qtls_clusters[(qtls_clusters$Difference<=threshold),]
kable(table(data$Gene))

ggplot(data, aes(x=Difference)) + 
  geom_histogram(binwidth=20, fill="grey60", colour="black") +
  theme_classic()+
  theme(text = element_text(size=20, colour="black")) +
  xlab("Distance between DEG and QTL markers")

withins$filter<-paste(withins$Gene, withins$qtl)
filter<-withins[!(duplicated(withins$filter)),]

write.csv(data, file="~/Documents/InnoVar/QTLs/QTL_data/LR/LR_QTLs_less_that_1%kb.csv")
save.image("~/Documents/InnoVar/QTLs/QTL_data/LR/LR_data.RData")
