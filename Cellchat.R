library(conflicted)
library(CellChat)
library(patchwork)
library(tidyr)
library(AnnotationDbi)
library(base)
library(stats)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomicRanges)
library(ggplot2)
#library(monocle3)
library(dplyr)
library(tibble)
library(purrr)
library(future)
library(igraph)
library(SummarizedExperiment)
library(Seurat)
library(SeuratObject)
setwd("/home/OSUMC.EDU/li176/work/Jianying/Carson/result/cellchat")
dir.create("./A")
setwd("./B")

table(dat.3_rename$timepoint)

Idents(dat.3_rename) <- "Cluster"

table(dat.3_rename$Cluster)

dat.4 <- subset(dat.3_rename, idents=c("MDSC","Exhausted CD8","Treg",  "CD8 effector", "CD4 EM", "CD4 naive", "Naive CD8","NK CD56 Bright","NK CD56 dim"))

Idents(dat.4) <- "patient"
dat.4A <- subset(dat.4, idents = c("2","6","11","16","17","23"))
dat.4B <- subset(dat.4, idents = c("1","4","5","15","10","20","21","22"))

cellchat <- createCellChat(object = dat.4A,group.by = "Cluster",assay = "RNA")
# group.by = "timepoint",
head(cellchat)
#Step 1: Access the ligand-receptor interaction information in CellChatDB
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
interaction_input <- cellchat@DB$interaction
complex_input <- cellchat@DB$complex
cofactor_input <- cellchat@DB$cofactor
geneInfo <- cellchat@DB$geneInfo
write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
write.csv(complex_input, file = "complex_input_CellChatDB.csv")
write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")
write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")
#Step 2: Update the required files by adding users’ curated ligand-receptor pairs
#Step 3: Update CellChatDB


showDatabaseCategory(cellchat@DB)


dplyr::glimpse(cellchat@DB$interaction)
#use a subset of CellChatDB for cell-cell communication analysis
cellchat@DB.use <- subsetDB(cellchat@DB, search = "Secreted Signaling") # use Secreted Signaling
#subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat,features = NULL) # This step is necessary even if using the whole database

cellchat<- identifyOverExpressedGenes(cellchat)

cellchat<- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)
length(cellchat@idents)
#[1] 14480
dim(cellchat@data.signaling)
#[1]   606 14480
unique(cellchat@idents)
conflicts_prefer(base::setdiff)
#The function droplevels is used to drop unused levels from a factor or, more commonly, from factors in a data frame.
cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))
#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Extract the inferred cellular communication network as a data frame
## We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,

df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
##Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
table(dat.3_rename$Cluster)


df.net <- subsetCommunication(cellchat, sources.use = c("MDSC"), targets.use = c("Exhausted CD8","Treg",  "CD8 effector", "CD4 EM", "CD4 naive", "Naive CD8","NK CD56 Bright","NK CD56 dim")) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
cellchat <- computeCommunProbPathway(cellchat)
#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
ggsave('netVisual.Number of interactions.pdf',width = 6, height = 4, units = "in")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
ggsave('Interaction weights.pdf', width = 6, height = 4, units = "in")
#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
#Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
#par("mar")
#par(mar=c(1,1,1,1))
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}



netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(2:15), lab.cex = 0.5,legend.pos.y = 30) 



#pathway
pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Chord diagram
group.cellType <- c(rep("MDSC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

