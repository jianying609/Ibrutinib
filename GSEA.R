library(Seurat)
library(dplyr)
library(future)
library(sctransform)
library(ggplot2)
library(Matrix)
library(cowplot)
library(patchwork)
library(here)
library(qs)
library(Polychrome)
library(RColorBrewer)
library(tidyverse)
library(enrichR)
library(EnhancedVolcano)
library(ggrepel)
library(Rcpp)
library(fgsea)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi) 
library(org.Mm.eg.db)

setwd("/home/OSUMC.EDU/li176/work/Jianying/Wenlab")
dir.create("./GSEA")
setwd("./GSEA")


# compare cluster ko and cluster wt



markers <- read.csv("mac_DEG_results_filted_wt_vs_ko.csv")

## Take the ave_logFC out MDSC cells

stats_mac_LFC <- markers$logfoldchanges
head(stats_mac_LFC)## Take the variable gene names out

names(stats_mac_LFC) <- markers$names
stats_mac_LFC <- stats_mac_LFC[!is.na(names(stats_mac_LFC))]
head(stats_mac_LFC)
## here we are looking into the database for human, where gene names are saved into ENTREZID, 
# so by the follwoing step we are converting genenames in our sample saved in Symbol format into ENTREZID format using mapIds command
names(stats_mac_LFC) <- mapIds(org.Mm.eg.db, keys = names(stats_mac_LFC), keytype = "SYMBOL", column = "ENTREZID")
head(stats_mac_LFC)

#How to look at how many duplicated genes
#Here we are removing duplicated genes 
stats_mac_LFC <- stats_mac_LFC[!duplicated(names(stats_mac_LFC))]
#check how many NA I have in gene.diff_gsea
sum(is.na(names(stats_mac_LFC)))
#[1] 1
length(stats_mac_LFC)
#[1] 4309

#Here we are removing rows that fial to switch SYMBOL to ENTREZID
gene.diff_Mac_gsea <- stats_mac_LFC[!is.na(names(stats_mac_LFC))]#300
sum(is.infinite(gene.diff_Mac_gsea))

# 0
#remove rows that have Inf
gene.diff_Mac_gsea <- stats_mac_LFC[!is.infinite(stats_mac_LFC)]

#download Genesets
#H hallmark gene sets (rdata file)
#C2 curated gene sets (rdata file)
#C3 motif gene sets (rdata file)
#C4 computational gene sets (rdata file)
#C5 GO gene sets (rdata file)
#C6 oncogenic signatures (rdata file)
#C7 immunologic signatures (rdata file)

load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c7_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c3_v5p2.rdata"))
load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c5_v5p2.rdata"))

## run the fgsea anaylses using pathways
## run the fgsea anaylses using pathways
fgsea <- fgsea(Mm.H, stats = gene.diff_Mac_gsea,  eps = 0.0, minSize=15, maxSize = 500)
fgsea_c7 <- fgsea(Mm.c7, stats = gene.diff_Mac_gsea,  eps = 0.0, minSize=15, maxSize = 500)
fgsea_c2<- fgsea(Mm.c2, stats = gene.diff_Mac_gsea,  eps = 0.0, minSize=15, maxSize = 500)
fgsea_c3 <- fgsea(Mm.c3, stats = gene.diff_Mac_gsea,  eps = 0.0, minSize=15, maxSize = 500)
fgsea_c5 <- fgsea(Mm.c5, stats = gene.diff_Mac_gsea,  eps = 0.0, minSize=15, maxSize = 500)

#select the sig pathway
fgsea_c7_sig <- fgsea_c7[which(fgsea_c7$pval <= 0.05), ]
fgsea_c2_sig <- fgsea_c2[which(fgsea_c2$pval <= 0.05), ]
fgsea_c5_sig <- fgsea_c5[which(fgsea_c5$pval <= 0.05), ]
fgsea_c3_sig <- fgsea_c3[which(fgsea_c3$pval <= 0.05), ]
fgsea_sig <- fgsea[which(fgsea$pval <= 0.05), ]
#convert entrezid to gene symbol mac
#convert entrezid to gene symbol myloid
for (i in 1:length(fgsea$leadingEdge)) {
  fgsea$leadingEdge[[i]] <- mapIds(org.Mm.eg.db, keys = fgsea_sig$leadingEdge[[i]], column = "SYMBOL", keytype = "ENTREZID")
}

for (i in 1:length(fgsea_c7_sig$leadingEdge)) {
  fgsea_c7_sig$leadingEdge[[i]] <- mapIds(org.Mm.eg.db, keys = fgsea_c7_sig$leadingEdge[[i]], column = "SYMBOL", keytype = "ENTREZID")
}

for (i in 1:length(fgsea_c2_sig$leadingEdge)) {
  fgsea_c2_sig$leadingEdge[[i]] <- mapIds(org.Mm.eg.db, keys = fgsea_c2_sig$leadingEdge[[i]], column = "SYMBOL", keytype = "ENTREZID")
}

for (i in 1:length(fgsea_c5_sig$leadingEdge)) {
  fgsea_c5_sig$leadingEdge[[i]] <- mapIds(org.Mm.eg.db, keys = fgsea_c5_sig$leadingEdge[[i]], column = "SYMBOL", keytype = "ENTREZID")
}

for (i in 1:length(fgsea_c3_sig$leadingEdge)) {
  fgsea_c3_sig$leadingEdge[[i]] <- mapIds(org.Mm.eg.db, keys = fgsea_c3_sig$leadingEdge[[i]], column = "SYMBOL", keytype = "ENTREZID")
}



#for (i in 1:length(fgsea_c3$leadingEdge)) {
#  fgsea_c3$leadingEdge[[i]] <- mapIds(org.Mm.eg.db, keys = fgsea_c3$leadingEdge[[i]], column = "SYMBOL", keytype = "ENTREZID")
#}
fgsea$leadingEdge <- as.character(fgsea$leadingEdge)
fgsea_c7$leadingEdge <- as.character(fgsea_c7$leadingEdge)
fgsea_c2$leadingEdge <- as.character(fgsea_c2$leadingEdge)
fgsea_c5$leadingEdge <- as.character(fgsea_c5$leadingEdge)
fgsea_c3$leadingEdge <- as.character(fgsea_c3$leadingEdge)
write.csv(fgsea, file = "mac.fgsea.csv")
write.csv(fgsea_c7, file = "mac.fgsea_c7.csv")
write.csv(fgsea_c2, file = "mac.fgsea_c2.csv")
write.csv(fgsea_c3, file = "mac.fgsea_c3.csv")
write.csv(fgsea_c5, file = "mac.fgsea_c5.csv")
mac.pathway <- read.csv("/fs/ess/PCON0022/Jianyingli/project_experiment/data/macpathway.csv")
head(mac.pathway)


#delete the duplicate
# MDSC.pathway <- MDSC.pathway[-10,]
MDSC.pathway$leadingEdge <- as.character(MDSC.pathway$leadingEdge)
fgesaResTidy <- MDSC.pathway%>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) 
conflict_prefer("filter", "dplyr")
fgesaResTidy$pathway <- factor(fgesaResTidy$pathway, 
                               levels = fgesaResTidy$pathway[order(fgesaResTidy$NES)]) 



# order by normalized enrichment score (NES)
area.color <- c("DMSO","DMSO","IBR","IBR","IBR") # to specify pos and neg; need to adjust according to the figure
p2 <- ggplot(data=fgesaResTidy, aes(x=pathway, y=NES, fill=area.color)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("red" ,"blue")) +
  theme_classic() + 
  coord_flip()+
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="MDSC pathways Enrichment Score from GSEA")

png('MDSC pathways Enrichment Score from GSEA.png',width=1480,height=1000)
print(p2)
dev.off()





