library(grid)
library(gridExtra)
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)

setwd("~/Documents/Hotspots/Paper_version_4/QTLs/")
Fg_final <- read.csv("~/Documents/Hotspots/Paper_version_4/Fg_final.csv")
all_markers_positions <- read.csv("~/Documents/Hotspots/Paper_version_4/all_markers_positions.csv", row.names=1)
final_hotspots <- read.csv("~/Documents/Hotspots/Paper_version_4/final_hotspots.csv")
QTL_database <- read.csv("~/Documents/Hotspots/Paper_version_4/QTLs/QTL_database.csv")
QTL_database$Chromosomes<-paste("chr", QTL_database$Chromosomes, sep="")


# check how many qtl markers are present in the marker database file, and check location of them
qtl_indb<-subset(all_markers_positions, all_markers_positions$Feature %in% QTL_database$Linked.markers)
write.csv(qtl_indb, file="~/Documents/Hotspots/Paper_version_4/QTLs/available_markers.csv")

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
write.csv(mismatches, file="~/Documents/Hotspots/Paper_version_4/QTLs/mismatches.csv")

# ========================================================================================================================================
# this section makes a new file that has no mismatches (i.e. markers that have multipe hits on the same chromosome)

all_markers_positions$ID<-paste(all_markers_positions$Chromosome, all_markers_positions$Start, all_markers_positions$Feature)
no_dups<-all_markers_positions[!(duplicated(all_markers_positions$ID)),]
no_dups$ID<-paste(no_dups$Chromosome, no_dups$Feature)
no_dups<-no_dups[!(duplicated(no_dups$ID)),]
all_markers_positions_filtered<-no_dups
write.csv(all_markers_positions_filtered, file="~/Documents/Hotspots/Paper_version_4/all_markers_positions_filtered.csv")

qtl_indb_filtered<-subset(all_markers_positions_filtered, all_markers_positions_filtered$Feature %in% QTL_database$Linked.markers)

# ========================================================================================================================================
qtls_clusters<-list()
all_qtls<-list()
qtls<-QTL_database$General.number
for(qtl in qtls){
  data2<-QTL_database[(QTL_database$General.number==qtl),]
  marker<-data2$Linked.markers
  for(i in marker){
    data<-data2[(data2$Linked.markers==i),]
    chr<-data$Chromosomes
    position<-subset(all_markers_positions_filtered, all_markers_positions_filtered$Chromosome %in% chr)
    position<-subset(position, position$Feature %in% i)
    clusters<-subset(final_hotspots, final_hotspots$Chromosome %in% chr)
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
qtls_clusters$Difference<-abs(as.numeric(qtls_clusters$start) - as.numeric(qtls_clusters$marker_position))/1000000
write.csv(qtls_clusters, file="~/Documents/Hotspots/Paper_version_4/QTLs/qtls_and_clusters.csv")

# ====================================================================================================================
# A histogram of distance between QTLs and FRGCs

ggplot(qtls_clusters, aes(x=Difference)) + 
  geom_histogram(binwidth=20, fill="grey60", colour="black") +
  theme_classic()+
  theme(text = element_text(size=20, colour="black"))

# identify percentiles
quantile(qtls_clusters$Difference, probs = c(0.01, 0.02, 0.05))

# lets examine qtls that are <5 mbp from a FRGC

threshold<-5
data<-qtls_clusters[(qtls_clusters$Difference<=threshold),]
write.csv(data, file="~/Documents/Hotspots/Paper_version_4/QTLs/qtls_lessthat_5mbp.csv")

for(marker in data$marker){
  newdata<-data[(data$marker==marker),]
  min<-min(newdata$marker_position) - 5000000
  max<-max(newdata$marker_position) + 5000000
  chr<-as.character(unique(newdata$Chromosome))
  
  df<-Fg_final[(Fg_final$Chromosome==chr),]
  df<-df[(df$start>=min),]
  df<-df[(df$end<=max),]
  markers<-qtls_clusters[(qtls_clusters$Chromosome==chr),]
  markers<-markers[(markers$marker_position>=min),]
  markers<-markers[(markers$marker_position<=max),]
  
  ggplot(df, aes(x=end, y=density)) +
    geom_jitter(size=3, aes(colour=Colour), alpha=0.6) +
    xlab("Position (bp)") + 
    ylab("Gene Density") + theme_bw() +
    theme(text = element_text(size=16, colour="black")) +
    scale_color_manual( values=c("grey60", "orangered2")) +
    # coord_cartesian(ylim=c(0,1), xlim=c(min, max))
    geom_vline(data =markers, aes(xintercept=marker_position), colour="green")+
    ggtitle(paste(chr))
}

data = qtls_clusters[(qtls_clusters$Chromosome=="chr2D"),]
data<-data[,c(14,15)]
data<-data[!duplicated(data$marker),]
ggplot(Fg_final[(Fg_final$Chromosome=="chr2D"),], aes(x=end, y=density)) +
  geom_jitter(size=3, aes(colour=Colour), alpha=0.6) +
  xlab("Position (bp)") + 
  ylab("Gene Density") + theme_bw() +
  theme(text = element_text(size=16, colour="black")) +
  # scale_color_manual( values=c("grey60", "orangered2", )) +
  # coord_cartesian(ylim=c(0,1), xlim=c(0, max(Fg_final[(Fg_final$Chromosome=="chr2D"),])))
  geom_vline(data = data, aes(xintercept=marker_position, color = marker))+
  ggtitle("Chromosome 2D")
}
# ====================================================================================================================