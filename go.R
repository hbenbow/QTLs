# Go graphs

stb<-read.delim("~/Documents/InnoVar/QTLs/QTL_data/STB/direct_go_count_bp_biomart_data_taestivum_eg_gene_v1_1_070319_corrected_extraction_20210211_154100.txt")
stb$BP<-row.names(stb)

ggplot(stb[1:30,], aes(x=reorder(BP, GO), y=GO)) +
  geom_col() +
  coord_flip() +
  theme_classic()+
  theme(text = element_text(size=20, colour='black'), axis.text.x = element_text(colour="black")) +
  ylab("Number of Genes") +
  xlab("High level biological process")
ggsave("~/Documents/InnoVar/QTLs/raw_figs/stb.pdf")


fhb<-read.delim("~/Documents/InnoVar/QTLs/QTL_data/Fg/direct_go_count_bp_biomart_data_taestivum_eg_gene_v1_1_070319_corrected_extraction_20210211_151714.txt")
fhb$BP<-row.names(fhb)

ggplot(fhb[1:30,], aes(x=reorder(BP, GO), y=GO)) +
  geom_col() +
  coord_flip() +
  theme_classic()+
  theme(text = element_text(size=20, colour='black'), axis.text.x = element_text(colour="black")) +
  ylab("Number of Genes") +
  xlab("High level biological process")
ggsave("~/Documents/InnoVar/QTLs/raw_figs/fhb.pdf")
