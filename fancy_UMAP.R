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
library(scRepertoire)

umap = dat.3_rename@reductions$umap@cell.embeddings %>%  #坐标信息
    as.data.frame() %>% 
    cbind(patient = dat.3_rename@meta.data$patient)%>%
    cbind(disease = dat.3_rename@meta.data$disease)%>%# 注释后的label信息 ，改为res.0.5 
    cbind(timepoint = dat.3_rename@meta.data$timepoint)%>%
    cbind(response = dat.3_rename@meta.data$response)%>%
    cbind(res.0.5 = dat.3_rename@meta.data$integrated_snn_res.0.5)%>% 
    cbind(cell_type = dat.3_rename@meta.data$Cluster)%>% 
    cbind(treatment = dat.3_rename@meta.data$treatment)%>%
    cbind(tcr_clonetype=dat.3_rename@meta.data$trc_clonetype)
    
head(dat.3_rename)

table(dat.3_rename$tcr_clonetype)

head(umap)


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


tab.1=table(dat.3_rename$response,dat.3_rename$patient)
balloonplot(tab.1)
################### 
Cellratio <- prop.table(tab.1, margin = 2)
Cellratio <- as.data.frame(Cellratio)
allcolour=c("#DC143C","#0000FF","#20B2AA","#F08080","#FFA500","#228B22","#1E90FF","#32CD32","#E9967A",
            "#808000","#FF00FF","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#40E0D0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#E9967A","#4682B4","#FFFFE0","#EE82EE"
            ,"#6A5ACD","#9932CC","#8B008B","#8B4513")

ggplot(Cellratio) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
    theme_classic() +
    labs(x='timepoint',y = 'Ratio')+
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
