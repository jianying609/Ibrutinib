getwd()
setwd("/home/OSUMC.EDU/li176/work/Jianying/Carson/result")
#dir.create("./UMAP")
#setwd("./UMAP")
## choose the library
library(Seurat)
library(dplyr)
library(future)
library(sctransform)
library(ggplot2)
library(Matrix)
library(Seurat)
library(cowplot)
library(ggrepel)
library(patchwork)
library(here)
library(qs)
library(Polychrome)
library(RColorBrewer)
library(tidyverse)
library(enrichR)
library(EnhancedVolcano)
library(ggrepel)
library(msigdbr) 
library(fgsea)
library(tibble)
library(gghighlight)
library(pacman)
library(DESeq2)
#### Reintergrate the data
# Check the treatment 
# Way1: Doing it using Seurat function
dat.3_rename <- PercentageFeatureSet(dat.3_rename, "^MT-", col.name = "percent_mito")

# Way2: Doing it manually
total_counts_per_cell <- colSums(  dat.3_rename@assays$RNA@counts  )
mito_genes <- rownames(dat.3_rename)[grep("^MT-",rownames(dat.3_rename))]
dat.3_rename$percent_mito <- colSums(  dat.3_rename@assays$RNA@counts[mito_genes,]  ) / total_counts_per_cell

head(mito_genes,10)
# Way1: Doing it using Seurat function
dat.3_rename <- PercentageFeatureSet(dat.3_rename, "^RP[SL]", col.name = "percent_ribo")

# Way2: Doing it manually
ribo_genes <- rownames(dat.3_rename)[grep("^RP[SL]",rownames(dat.3_rename))]
head(ribo_genes,10)
dat.3_rename$percent_ribo <- colSums(  dat.3_rename@assays$RNA@counts[ribo_genes,]  ) / total_counts_per_cell
selected_mito <- WhichCells(dat.3_rename, expression = percent_mito < 0.25)
selected_ribo <- WhichCells(dat.3_rename, expression = percent_ribo > 0.05)

# and subset the object to only keep those cells
dat.3_rename <- subset(dat.3_rename, cells = selected_mito)
dat.3_rename <- subset(dat.3_rename, cells = selected_ribo)
dim(dat.3_rename)
head(dat.3_rename)







################### doublet finder
library(DoubletFinder)
library(scater)
pk.list <- paramSweep(dat.3_rename, PCs = 1:20, sct = T)
pk.stats <- summarizeSweep(pk.list, GT = FALSE)
bcmvn <- find.pK(pk.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

DoubletRate = 0.039                     # 5000细胞对应的doublets rate是3.9%
homotypic.prop <- modelHomotypic(dat.3_rename$Cluster)   # 最好提供celltype
nExp_poi <- round(DoubletRate*ncol(dat.3_rename)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## 使用确定好的参数鉴定doublets
dat.3_rename <- doubletFinder(dat.3_rename, PCs = 20, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F, sct = T)

## 结果展示，分类结果在dat.3_rename@meta.data中
DimPlot(dat.3_rename, reduction = "umap", group.by = "DF.classifications_0.25_0.3_171")
## 使用确定好的参数鉴定doublets
dat.3_rename <- doubletFinder_v3(dat.3_rename, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F, sct = T)















table(OSU_20210821$patient)
dat1 <- OSU_20210821
# get the count information from the R object
counts <- GetAssayData(dat1, assay = "RNA",slot = "counts")
#  check the row names of the count
rownames(counts)
# split the object
tissue.list <- SplitObject(dat1, split.by = "patient")
# normalize and identify variable features for each dataset independently
tissue.list<- lapply(X = tissue.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across dat asets for integration
features <- SelectIntegrationFeatures(object.list = tissue.list)
##########################
immune.anchors <- FindIntegrationAnchors(object.list = tissue.list, anchor.features = features)
all.genes <- rownames(experiment.aggregate)
#filter genes of experiment.aggregat


DefaultAssay(dat.3_rename) <- "RNA"
dat.3_rename[["percent.mt"]]<- PercentageFeatureSet(dat.3_rename, pattern = "^MT-")
dat.3_rename <- subset(dat.3_rename, subset= nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <10)



DimPlot(experiment.aggregate, reduction = "umap", group.by = "condition",pt.size=0.2)
dim_plot <- DimPlot(experiment.aggregate, reduction = "umap", split.by = "condition",pt.size=0.2)
DimPlot(exp_rename, reduction = "umap",group.by = "seurat_clusters",label = T, pt.size = 0.2)

# this command creates an 'integrated' data assay
dat.3<- IntegrateData(anchorset = immune.anchors)
#dat.3[["percent.mt"]]<- PercentageFeatureSet(dat.3, pattern = "^MT-")
#dat.3 <- subset(dat.3, subset= nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <10)
dat.3 <- NormalizeData(dat.3)
dat.3 <- FindVariableFeatures(dat.3, selection.method = "vst", nFeature = 2000)
# Run PCA
all.genes <- rownames(dat.3)
all.genes
DefaultAssay(dat.3) <- "integrated"
dat.3 <- ScaleData(dat.3, verbose = FALSE)
dat.3 <- RunPCA(dat.3, npcs = 20, verbose = FALSE)
# PCA
ElbowPlot(dat.3)
dat.3_rename <- RunUMAP(dat.3_rename, reduction ="pca", dims =1:20)
dat.3 <- FindNeighbors(dat.3, reduction = "pca",dims =1:20)
dat.3 <- FindClusters(dat.3, resolution = 0.5)

dat.3_rename[["percent.mt"]]<- PercentageFeatureSet(dat.3_rename, pattern = "^MT-")
dat.3 <- subset(dat.3_rename, subset= nFeature_RNA >200 & nFeature_RNA < 4000 & percent.mt <10)
#UMAP plot
DimPlot(dat.3, reduction = "umap", label = T, pt.size = 0.2)

# Find markers
markers <- FindAllMarkers(dat.3, only.pos = TRUE, logfc.threshold = 0.25 )
write.csv(x = markers, file = "postive.csv")
head(markers)
table(markers$cluster)
# Rename Cluster
# 9 FCGR3A VMO1 # 10
dat.3_rename <- subset(dat.3, idents = c("0","1","2","3","4","5","6","7","8", "9", "10","12","13","14","15","16","17","18","19","20","21"))
tmp_ident <- as.factor(droplevels(dat.3_rename@meta.data$seurat_clusters))
#"Plac8", "Irf7","Klrc2","Cxcl10","Cd8b1","Il7r","Cd69","Fn1","H2-Ab1","Igkc","Cd4","Foxp3","Fscn1","S100a9","Stmn1"
levels(tmp_ident) <- c("MDSC", "CD4 EM","CD4 naive", "NK CD56 dim", "CD14+ Monocyte","CD8 effector","MDSC", "Exhausted CD8","Mature B cells","CD16+ Monocyte","Megakaryocyte",
                       "NK CD56 Bright","Treg","cDC","Naive CD8","MDSC","pDC","Plasma cells","HSPC","Platelet","Erythroid")
dat.3_rename <- AddMetaData(dat.3_rename, tmp_ident, "Cluster")

table(dat.3_rename$Cluster)
dat.3_rename$integrated_snn_res.0.5
Idents(dat.3_rename)<-"integrated_snn_res.0.5"
DimPlot(dat.3_rename, reduction = "umap", label = T, pt.size = 0.2)
Idents(dat.3)<-"seurat_clusters"
DimPlot(dat.3, reduction = "umap", label = T, pt.size = 0.2)
   #
FeaturePlot(dat.3_rename, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                                                  "CD8A"))
#
markers.to.plot <- c("S100A8", "S100A9", "", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57","FOXP3")
DotPlot(dat.3_rename, features = markers.to.plot, dot.scale = 8) +
  RotatedAxis()
dat.3_rename$response
table(dat.3_rename$patient)
Idents <- dat.3_rename$Cluster
DimPlot(dat.3_rename, reduction = "umap", label = F, pt.size = 0.2, split.by  = "response")

############################################################
##############Benefit vs progressive#################
########################################################

Idents(dat.3_rename) <- "patient"
dat.3_re_1 <- subset(dat.3_rename, idents = c("1","2","4","5","6","10","11","15","16",
                                            "17", "20", 
                                            "21", "22", "23"))
Idents(dat.3_re_1) <- "patient"
dat.3_re <- subset(dat.3_re_1,  idents = c("1","2","4","5","6","10","11","15","16",
                                           "17", "20", 
                                           "21", "22", "23"))
tmp_ident <- as.factor(droplevels(dat.3_re@meta.data$patient))
#"Plac8", "Irf7","Klrc2","Cxcl10","Cd8b1","Il7r","Cd69","Fn1","H2-Ab1","Igkc","Cd4","Foxp3","Fscn1","S100a9","Stmn1"
levels(tmp_ident1) <- c("B", "A","B","B","A",
                       "B","A","B",
                       "A","A","B",
                       "B","B","A")
dat.3_rename <- AddMetaData(dat.3_re, tmp_ident, "Benefit")
table(dat.3_rename$Benefit)






















#An object of class Seurat 
#13714 features across 2638 samples within 1 assay 
#Active assay: RNA (13714 features, 2000 variable features)
# 3 dimensional reductions calculated: pca, umap, tsne
table(dat.3_rename$treatment)
dat6 <- dat.3_rename

Idents(dat6) <- "patient"
dat.5 <- subset(dat6, idents = c("1","4","5","10","15","20","21","22"))
#dat.5 <- subset(dat6, idents = c("2","6","11","16","17","23"))
dat.3_rename <- dat.5

umap = dat.3_rename@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(patient = dat.3_rename@meta.data$patient)%>%
  cbind(disease = dat.3_rename@meta.data$disease)%>%# 注释后的label信息 ，改为res.0.5 
  cbind(timepoint = dat.3_rename@meta.data$timepoint)%>%
  cbind(response = dat.3_rename@meta.data$response)%>%
  cbind(res.0.5 = dat.3_rename@meta.data$integrated_snn_res.0.5)%>% 
  cbind(cell_type = dat.3_rename@meta.data$Cluster)%>% 
  cbind(treatment = dat.3_rename@meta.data$treatment)%>% 
  cbind(cloneType = dat.3_rename@meta.data$cloneType)%>% 
  cbind(tcr_clonetype = dat.3_rename@meta.data$tcr_clonetype)
  
dat.3_rename$cdr3s_aa
dat.3_rename$unique_tcr
table(dat.3_rename$unique_tcr)
table(dat.3_rename$CTstrict)
table(dat.3_rename$tcr_clonetype,dat.3_rename$CTstrict )



allcolour=c("#DC143C","#0000FF","#20B2AA","#F08080","#FFA500","#228B22","#1E90FF","#32CD32","#E9967A",
            "#808000","#FF00FF","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#40E0D0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#E9967A","#4682B4","#FFFFE0","#EE82EE"
            ,"#6A5ACD","#9932CC","#8B008B","#8B4513")
p <- ggplot(umap,aes(x= umap_1 , y = umap_2 ,color = cell_type)) +
  geom_point(size = 0.2 , alpha =1 )  +
  scale_color_manual(values = allcolour)
p

p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
p2
p3 <- p2 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 
p3

p4 <- p3 + 
  geom_segment(aes(x = min(umap$umap_1) , y = min(umap$umap_2) ,
                   xend = min(umap$umap_1) +3, yend = min(umap$umap_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$umap_1)  , y = min(umap$umap_2)  ,
                   xend = min(umap$umap_1) , yend = min(umap$umap_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$umap_1) +1.5, y = min(umap$umap_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$umap_1) -1, y = min(umap$umap_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p4
#3.4 调整umap图 - repel - labels
#1）计算每个cluster的median 坐标位置
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
library(ggrepel)
p4 +
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines"))
p4 +
  geom_label_repel(aes(label=cell_type), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.2, "lines")) +
  theme(legend.position = "none")
p4

p3 + facet_wrap(~response , ncol = 3)

p3 + facet_wrap(~patient , ncol = 3)

p3 + facet_wrap(vars(patient,timepoint), labeller = "label_both")
p3 + facet_wrap(vars(disease,timepoint), labeller = "label_both", ncol=3)


p3 + facet_wrap(~timepoint, ncol = 3)
p3 + facet_wrap(~disease, ncol = 3)
p3 + facet_wrap(~treatment, ncol = 3)


as.data.frame(table(dat.3_rename$Cluster))
as.data.frame(table(dat.3_rename$timepoint))
dat.3_rename$patient


tab.1=table(dat.3_rename$timepoint,dat.3_rename$cloneType)

################### 
Cellratio <- prop.table(tab.1, margin = 2)
tab.1
Cellratio <- as.data.frame(Cellratio)
Cellratio
allcolour=c("#DC143C","#0000FF","#20B2AA","#F08080","#FFA500","#228B22","#1E90FF","#32CD32","#E9967A",
            "#808000","#FF00FF","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#40E0D0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#E9967A","#4682B4","#FFFFE0","#EE82EE"
            ,"#6A5ACD","#9932CC","#8B008B","#8B4513")

ggplot(Cellratio) + 
  geom_bar(aes(x =Var1, y= Freq, fill = Var2),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='timepoint',y = 'Relative abundance')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.2, linetype="solid"))     

##################### find the markers in the time point##################
#########Find the marekrs of the different timepoint############
######################
Idents(dat.3) <- "timepoint"
#Idents <- dat.3_rename@meta.data$timepoint
DefaultAssay(dat.3_rename) <- "RNA"
table(dat.3_rename$Clusters)
marker1 <- FindMarkers(object = dat.3_rename, ident.1 = "C1D+1", ident.2 = "C1D-7") 

##### find the mareks from the different group.###############
Idents(dat.3_rename) <- "Cluster"
table(dat.3_rename$Cluster)
dat.4 <- subset(dat.3_rename, idents = c("Mature B cells"))
                                         #"CD4 EM","CD4 naive", "NK CD56 dim", "CD14+ Monocyte","CD8 effector","Exhausted CD8","Mature B cells","CD16+ Monocyte",
                                        #"NK CD56 Bright","Treg","cDC","Naive CD8","pDC"))
Idents(dat.4) <- "patient"

marker2 <-  FindMarkers(object = dat.4, ident.1 = c("2","6","11","16","17","23"),ident.2 = c("1","4","5","15","10","20","21","22") )

write.csv(marker2,"Group A vs Group B_Mature B cells.csv")

################ check the gene list 
DefaultAssay(dat.3_rename) <- "RNA"

Idents(dat.3_rename) <- "timepoint"
dat.3_rename_CD <- subset(dat.3_rename, idents = c("C1D-7"))
Idents(dat.3_rename_CD) <- "Cluster"


dat.3_rename_new <- subset(dat.3_rename_CD, idents = c("MDSC","CD4 EM","CD4 naive", "NK CD56 dim", "CD14+ Monocyte","CD8 effector","Exhausted CD8","Mature B cells","CD16+ Monocyte",
                                        "NK CD56 Bright","Treg","cDC","Naive CD8","pDC","HSPC"))

getwd()
DefaultAssay(dat.3_rename_CD) <- "RNA"
png("DotPlot.png", width = 5, height = 10, units = "in", res = 300)
DotPlot(dat.3_rename_CD, 
        #features = c("BID", "BBC3", "PMAIP1", "BAX", "BCL2L11",  "BCL2L1", "BCL2L2", "XIAP", "MCL1", "BCL2A1","CFLIP","BCL2","BRD4"), 
        features = ("BRD4"),
        group.by = "Cluster", 
        scale = T,
        dot.scale = 12, 
    )

dev.off()

# you can plot raw counts as well
Idents(dat.3_rename_CD) <- "Cluster"
png("VlnPlot_BRD4.png", width = 10, height = 5, units = "in", res = 300)
VlnPlot(dat.3_rename_CD, features = c("BRD4"))
dev.off()

# you can plot raw counts as well

dat.3_1 <- ScaleData(dat.3_rename, features=rownames(dat.3_rename))

marker <- c("BRD4")
Idents(dat.3_1) <- "Cluster"
png("heatmap_BRD4.png", width = 10, height = 5, units = "in", res = 300)
# DoHeatmap(dat.3_1, features = c("BRD4"))

marker = as.character(marker$gene)

DoHeatmap(dat.3_1, features = c("BRD4"), label = TRUE,slot = "scale.data",disp.min = -2.5,group.by = "Cluster" ) + scale_fill_gradient(low = "white", high = "red") +
  theme( axis.text.y = element_text(face = "bold", size = 20, angle = 0),legend.key.height= unit(1, 'cm'), legend.title = element_text(size=30),legend.text = element_text(face="bold",size = 20))
dev.off()

RidgePlot(dat.3_rename_CD, features = c("BRD4"))
DoHeatmap(dat.3_1, features = ("BRD4"))


