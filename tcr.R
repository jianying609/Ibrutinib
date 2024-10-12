library(Seurat)
library(scRepertoire)
library(ggsci)
library(vegan)
quantContig(dat.T, cloneCall = "CTgene+CTnn",scale = T)
getwd()
setwd("/home/OSUMC.EDU/li176/work/Jianying/Carson/result")
dir.create("./TCR")
setwd("./TCR")

p2 + p3
Idents(dat.3_rename) <- "Cluster"
table(dat.3_rename$Cluster)
dat.T <- subset(dat.3_rename,idents=c("CD4 EM","CD4 naive","CD8 effector","Treg","Naive CD8","Exhausted CD8"))
occupiedscRepertoire(dat.T, x.axis = "timepoint")
occupiedscRepertoire(dat.3_rename, x.axis = "patient")
table(dat.T$patient)

alluvialClonotypes(dat.T, cloneCall = "gene", 
                   y.axes = c("treatment","timepoint",  "response"), 
                   color = "ident")
head(dat.3_rename)
?clonalDiversity()

diversity_matrics <- clonalDiversity(dat.T, "treatment")
clonalDiversity(dat.T, 
                split.by = "treatment", 
                group.by ="patient",
                #x.axis = "tretment",
                cloneCall = "nt")
Idents(dat.T) <- "patient"
table(dat.T$patient,dat.T$cloneType)

table(dat.T$patient, dat.T$cloneType)

datA.T <- subset(dat.T, idents=c("2","6","11","16","17","23"))
datB.T <- subset(dat.T, idents=c("1","4","5","15","10","20","21","22"))

png("Proportion_groupA.png", width = 5, height = 10, units = "in", res = 300)

combined2 <- expression2List(datA.T, 
                             split.by = "timepoint")
table(datA.T$timepoint,datA.T$cloneType)
length(combined2)
#now listed by cluster
ow listed by cluster
clonalHomeostasis(combined2, 
                  cloneCall = "nt",
                  split.by = "timepoint")

dev.off()

png("Proportion_groupB.png", width = 5, height = 10, units = "in", res = 300)


combined2 <- expression2List(datB.T, 
                             split.by = "timepoint")

clonalHomeostasis(combined2, 
                  cloneCall = "nt",
                  split.by = "timepoint")
dev.off()

Idents(dat.T) <- "Cluster"
clonalDiversity(dat.T, 
                split.by = "treatment", 
                group.by ="Cluster",
                cloneCall = "nt")
##############3

# 基于 meta 属性中的 ident 列
png("Cellresponse_groupA.png", width = 5, height = 10, units = "in", res = 300)
alluvialClonotypes(datA.T, cloneCall = "gene", 
                   y.axes = c("timepoint", "Cluster"), 
                   color = "timepoint")
dev.off()

png("Cellresponse_groupB.png", width = 5, height = 10, units = "in", res = 300)
alluvialClonotypes(datB.T, cloneCall = "gene", 
                   y.axes = c("timepoint", "Cluster"), 
                   color = "timepoint")
dev.off()
#############3SUBCLONE
library(circlize)
library(scales)

circles <- getCirclize(dat.T, 
                       group.by = "ident")

#Just assigning the normal colors to each cluster
grid.cols <- scales::hue_pal()(length(unique(dat.T@active.ident)))
names(grid.cols) <- levels(dat.T@active.ident)

#Graphing the chord diagram
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)

alluvialClonotypes(OSU_20210821, cloneCall = "gene", 
                   y.axes = c("patient", "ident", "timepoint"), 
                   color = "ident")

clonotypeBias(datA.T, 
              cloneCall = "aa", 
              split.by = "patient", 
              group.by = "ident",
              n.boots = 20, 
              min.expand = 10)
clonotypeBias(datB.T, 
              cloneCall = "aa", 
              split.by = "patient", 
              group.by = "ident",
              n.boots = 20, 
              min.expand = 10)
clonotypeBias(dat.T, 
              cloneCall = "aa", 
              split.by = "patient", 
              group.by = "ident",
              n.boots = 20, 
              min.expand = 10)
table(dat.3_rename$cloneType, dat.3_rename$)
