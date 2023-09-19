
rm(list=ls())

library(Seurat)
library(ggplot2)
library(dplyr)

setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/")

integrated <- readRDS("logNormalize_combined_mTmG_SalineAngII.rds")

DimPlot(integrated, label = T, raster=T)
integrated$Sample <- factor(integrated$Sample, levels=c("Saline.GFP", "Saline.RFP", "nonIMH.GFP", "nonIMH.RFP", "IMH.GFP", "IMH.RFP"))
DimPlot(integrated, split.by = "Sample", raster = T)

integrated$Lineage <- ifelse(integrated$Sample %in% c("Saline.GFP","nonIMH.GFP", "IMH.GFP"), "GFP", "RFP")

meta <- integrated@meta.data
meta <- meta %>% mutate(Group=case_when(
  Sample %in% c("Saline.GFP", "Saline.RFP") ~ "Control",
  Sample %in% c("nonIMH.GFP", "nonIMH.RFP") ~ "nonIMH",
  Sample %in% c("IMH.GFP", "IMH.RFP") ~ "IMH"
))
integrated$Group <- meta$Group
integrated$Group <- factor(integrated$Group, levels=c("Control","nonIMH","IMH"))

integrated$SMC <- ifelse(integrated$seurat_clusters %in% c(1,3,4,6,7,8,9), "SMC", "Others")
DimPlot(integrated, split.by = "orig.ident", group.by = "SMC", raster = T, cols = c("grey", "darkblue"))+NoLegend()

integrated$SMC2 <- ifelse(integrated$seurat_clusters %in% c(1,3,4,6,7,8,9), 
                          as.character(integrated$seurat_clusters), NA)
integrated$SMC2 <- factor(integrated$SMC2, levels=c(1,3,4,6,7,8,9,NA))
DimPlot(integrated, group.by = "SMC2", raster=T, label=T)
DimPlot(integrated, group.by = "SMC2", raster=T, label=T, split.by = "Sample")

FeaturePlot(integrated, features = c("Col1a2","Dcn", "S100a4","Pdgfra"), raster = T) #split.by = "Sample", 
ggsave("Fig4_featurePlot_TFs_Sample.pdf", width = 16, height = 18)

FeaturePlot(integrated, features = c("Nfia"), raster = T, split.by = "Sample") #, 
ggsave("Fig4_featurePlot_TFs_Sample.pdf", width = 16, height = 18)

############################# separate analysis of each dataset ###############################
setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/SeparateAnalysis/")

Saline.RFP <- Read10X_h5("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/CellRangerOut/Basal_RFP/filtered_feature_bc_matrix.h5")
Saline.RFP <- CreateSeuratObject(counts = Saline.RFP, project = "Saline.RFP", min.cells = 3, min.features = 200)

########### basic QC
Saline.RFP[["percent.mt"]] <- PercentageFeatureSet(Saline.RFP, pattern = "^mt-")
VlnPlot(Saline.RFP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("Saline.RFP_Fig1_QC1.pdf", width=8, height=5)

plot1 <- FeatureScatter(Saline.RFP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Saline.RFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
ggsave("Saline.RFP_Fig1_QC2.pdf", width=8, height=4)

Saline.RFP <- subset(Saline.RFP, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 12)

########### normalize, identify variable, scale 
Saline.RFP <- NormalizeData(Saline.RFP)
Saline.RFP <- FindVariableFeatures(Saline.RFP, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Saline.RFP), 10)
plot1 <- VariableFeaturePlot(Saline.RFP)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
ggsave("Saline.RFP_Fig2_variable genes.pdf", width=10, height=4)

all.genes <- rownames(Saline.RFP)
Saline.RFP <- ScaleData(Saline.RFP, features = all.genes)


########### lineal dimentional reduction, cluster cells, non-lineal dimentional reduction 
Saline.RFP <- RunPCA(Saline.RFP, features = VariableFeatures(object = Saline.RFP))
ElbowPlot(Saline.RFP)
ggsave("Saline.RFP_Fig3_ElbowPlot.pdf", width=8, height=4)

Saline.RFP <- FindNeighbors(Saline.RFP, dims = 1:18)
Saline.RFP <- FindClusters(Saline.RFP, resolution = 0.6)
Saline.RFP <- RunUMAP(Saline.RFP, dims = 1:18)
DimPlot(Saline.RFP, reduction = "umap", label=TRUE)+ggtitle("Saline.RFP")
ggsave("Saline.RFP_Fig4_UMAP.pdf", width=5.5, height=5)

########## find differentially expressed genes
markers <- FindAllMarkers(Saline.RFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "Saline.RFP_markersAll.csv", quote = F)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Saline.RFP, features = top10$gene)
ggsave("Saline.RFP_Fig5_top10Heatmap.pdf", width=20, height=20)

######### Feature plot
FeaturePlot(Saline.RFP, features = c("Myh11", "Acta2", "Dcn", "Lum", "Fbln1","Ly6a ", 
                                "Pecam1", "Lyz2", "Cd3d"))
ggsave("Saline.RFP_Fig6_FeaturePlot1.pdf", width=12, height=10)

saveRDS(Saline.RFP, file = "Saline.RFP_dim18res0.6.rds")


##################################### DoubletFinder ###########################################
library(DoubletFinder)
## pK Identification (no ground-truth) 
sweep.res.list_DMT <- paramSweep_v3(DMT, PCs = 1:18, sct = FALSE)
sweep.stats_DMT <- summarizeSweep(sweep.res.list_DMT, GT = FALSE)
bcmvn_DMT <- find.pK(sweep.stats_DMT)

## pK Identification (ground-truth) 
sweep.res.list_DMT <- paramSweep_v3(DMT, PCs = 1:18, sct = FALSE)
gt.calls <- DMT@meta.data[rownames(sweep.res.list_DMT[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_DMT <- summarizeSweep(sweep.res.list_DMT, GT = TRUE, GT.calls = gt.calls)
bcmvn_DMT <- find.pK(sweep.stats_DMT)

## Homotypic Doublet Proportion Estimate 
annotations <- DMT@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- DMT@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(DMT@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies 
DMT <- doubletFinder_v3(DMT, PCs = 1:18, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DMT <- doubletFinder_v3(DMT, PCs = 1:18, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_1001", sct = FALSE)

saveRDS(DMT, "DMT_dim18res0.6_rmDoublet.rds")



################################### integrate 6 datasets on Linux #######################################
#move to Linux
#srun --tasks=1 --cpus-per-task=2 --mem=210G --time=72:00:00 --pty /bin/bash
#module load R/4.2.0
#R
setwd("/project/lemaire/6_mouse_scRNAseq/")
library(Seurat)
library(dplyr)

files <- list.files(path=".", pattern = "*_dim18res0.6.rds$")
files <- files[c(1:4,8:9)]
#[1] "CMG_dim18res0.6.rds"        "CMT_dim18res0.6.rds"
#[3] "DMG_dim18res0.6.rds"        "DMT_dim18res0.6.rds"
#[5] "Saline.GFP_dim18res0.6.rds" "Saline.RFP_dim18res0.6.rds"

aorta.list <- lapply(files, readRDS)
names <- c("nonIMH.GFP", "nonIMH.RFP", "IMH.GFP", "IMH.RFP", "Saline.GFP", "Saline.RFP")
for (i in 1:length(aorta.list)) {
  aorta.list[[i]]$Sample <- names[i]
  aorta.list[[i]] <- NormalizeData(aorta.list[[i]], verbose = FALSE)
  aorta.list[[i]] <- FindVariableFeatures(aorta.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = aorta.list)

anchors <- FindIntegrationAnchors(object.list = aorta.list, dims = 1:30)

integrated <- IntegrateData(anchorset = anchors)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

DefaultAssay(integrated) <- "RNA"
saveRDS(integrated, "logNormalize_combined_mTmG_SalineAngII.rds")

############################# define the clusters #########################
markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "AllMarker_cluster_res0.5_mTmG_SalineAngII.csv", quote = F)

markers <- read.csv("AllMarker_cluster_res0.5_mTmG_SalineAngII.csv", header=T, row.names = "X", stringsAsFactors = F)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DotPlot(integrated, features = unique(top10$gene), cols=c("lightgrey", "darkblue"))+
  RotatedAxis()+ggtitle("Top10 integrate RNA")
ggsave("Fig2_dotplot_top10RNA.pdf", width=36, height=6)

VlnPlot(integrated, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0, ncol=1)+NoLegend()
VlnPlot(integrated, features = c("Bax", "Bad","Casp3"), pt.size = 0, ncol=1)+NoLegend()

FeaturePlot(integrated, features = c("Myh11", "Acta2", "Igfbp2","Tnnt2", 
                                     "Dcn", "Lum", "Ly6a", "Col8a1",  
                                     "Lyz2", "Cd74","Cdh5",
                                     "Top2a", "Mki67","Thbs1","Itk","Plp1"),
            raster=T)
FeaturePlot(integrated, features = c("nCount_RNA"),raster=T, split.by = "orig.ident",
            cols=c("darkblue", "gold"))

################################ majorCelltype of the combined dataset ############################
meta <- integrated@meta.data
meta <- meta %>% mutate(majorCelltype=case_when(
  seurat_clusters %in% c(1,3,4,6,7,8,9) ~ "SMC",
  seurat_clusters %in% c(0,5,10,11,16) ~ "FB",
  seurat_clusters %in% c(12) ~ "EC",
  seurat_clusters %in% c(2,14,15) ~ "Macrophage",
  seurat_clusters %in% c(13) ~ "Cycling",
  seurat_clusters %in% c(17) ~ "Gliocyte",
  seurat_clusters %in% c(18) ~ "T lymphocyte",
  seurat_clusters %in% c(19) ~ "19",
  seurat_clusters %in% c(20) ~ "20"
))
integrated$majorCelltype <- meta$majorCelltype
integrated$majorCelltype <- factor(integrated$majorCelltype, levels = c("SMC", "FB", "EC", "Macrophage","Gliocyte","T lymphocyte","Cycling","19","20"))
Idents(integrated) <- "majorCelltype"

DimPlot(integrated, label = T, raster = T, group.by = "majorCelltype")
DimPlot(integrated, split.by = "Sample", raster = T)
DimPlot(integrated, split.by = "Group", raster = T)


################################# basic visulize ####################################
DotPlot(integrated, features = c("Tph1", "Tph2","Ddc","Htr3a","Htr1a"))+coord_flip()+
  RotatedAxis() 


VlnPlot(integrated, features = "nCount_RNA", pt.size = 0)+NoLegend()
FeaturePlot(integrated, features = "Cd86", split.by = "Group", raster = T)

FeaturePlot(integrated, features = c("Spp1", "Cd44", "Itga8","Itgb1"), split.by = "Sample", raster = T)

FeaturePlot(integrated, features = c("Myh11", "Acta2"), raster = T)

VlnPlot(integrated, features = "nCount_RNA", group.by = "majorCelltype", split.by ="Lineage" ,pt.size = 0)+
  stat_summary(fun.y=median, geom="point", size=2, color="black", position = "jitter")+
  scale_fill_manual(values=c("LightSeaGreen","Salmon"))

integrated$majorCelltype_sample <- paste0(integrated$majorCelltype, "_", integrated$Sample)
DotPlot(integrated, features = c("Nfia"), group.by = "majorCelltype_sample")


######################################## cluster propotion  #############################################
library(ggbreak)
meta <- integrated@meta.data
meta1 <- meta[meta$Sample %in% c("Saline.GFP","nonIMH.GFP","IMH.GFP"),]
meta.pct <- meta1 %>% group_by(Group,majorCelltype) %>% 
  summarise(Percentage=n()) %>% 
  group_by(Group) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)
meta.pct$Group <- factor(meta.pct$Group, levels=c("Control","nonIMH","IMH"))
p1 <- ggplot(meta.pct, aes(x=majorCelltype, y=Percentage, fill=Group))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+xlab("Cluster")+ggtitle("GFP cluster proportion")+
  scale_y_break(c(5,70))+
  scale_fill_manual(values=c("darkgrey","darkgoldenrod1","firebrick3"))

meta2 <- meta[meta$Sample %in% c("Saline.RFP", "nonIMH.RFP", "IMH.RFP"),]
meta.pct <- meta2 %>% group_by(Group,majorCelltype) %>% 
  summarise(Percentage=n()) %>% 
  group_by(Group) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)
meta.pct$Group <- factor(meta.pct$Group, levels=c("Control","nonIMH","IMH"))
p2 <- ggplot(meta.pct, aes(x=majorCelltype, y=Percentage, fill=Group))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+xlab("Cluster")+ggtitle("RFP cluster proportion")+
  #scale_y_break(c(5,70))+
  scale_fill_manual(values=c("darkgrey","darkgoldenrod1","firebrick3"))

p1+p2

########################## differential analysis ##############################################
setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/Differential/") 

Idents(integrated) <- paste0(integrated$majorCelltype, "_", integrated$Group)
for(i in unique(integrated@meta.data$majorCelltype)){
  each <- FindMarkers(integrated, ident.1 = paste0(i,"_IMH"), ident.2 = paste0(i,"_nonIMH"))
  write.csv(each, paste0("wilcox_",i,"_IMHvsnonIMH.csv"), quote = F)
}
for(i in unique(integrated@meta.data$majorCelltype)){
  each <- FindMarkers(integrated, ident.1 = paste0(i,"_nonIMH"), ident.2 = paste0(i,"_Control"))
  write.csv(each, paste0("wilcox_",i,"_nonIMHvscontrol.csv"), quote = F)
}

# avg
avg <- AverageExpression(integrated, assays = "RNA")
avg <- avg$RNA
write.csv(avg, "avg_majorCelltype_Group.csv", quote=F)

#################################### SMC DEGs, classified as categories ###########################################
files <- list.files(path=".", pattern = "wilcox_SMC_*")
dfs <- lapply(files, function(x) read.csv(x, header=T, stringsAsFactors=F, row.names = "X"))
dfs <- lapply(dfs, function(x) x[which(x$p_val_adj<0.05 & abs(x$avg_log2FC)>0),])
DEGs <- unique(unlist(lapply(dfs, row.names)))

setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/SMC/")
d <- avg[row.names(avg) %in% DEGs, c(19,1,10)]
p <- pheatmap(d, scale="row",show_rownames = F, cluster_cols = F, cutree_rows = 5,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main=paste0("SMC DEGs"))

plot(p$tree_row)
abline(h=1.8, col="red", lty=2, lwd=1)
categories <- sort(cutree(p$tree_row, h=1.8))
categories <- as.data.frame(cbind(names(categories), categories))
write.csv(categories,"DEG_category.csv", quote = F)

for(i in unique(categories$categories)){
  pheatmap(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]), 
           scale = "row", show_rownames=F, cluster_cols = F, border_color=NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main=paste0("DEGs in Category", i),
           filename = paste0("Fig6_heatmap_DEG_C",i, ".pdf"), width=4, height=4)
  print(
    dim(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]))
  )
}

## sub category 1
avg <- read.csv("../avg_majorCelltype_Group.csv", header=T, row.names = "X", stringsAsFactors = F)
categories <- read.csv("DEG_category.csv",header=T, row.names = "X", stringsAsFactors = F)

d <- avg[row.names(avg) %in% categories[categories$categories==1,1], c(19,1,10)]
p <- pheatmap(d, scale="row",show_rownames = F, cluster_cols = F, cutree_rows = 3,
              color = colorRampPalette(c("blue", "white", "red"))(50),
              main=paste0("SMC DEGs category1"))

plot(p$tree_row)
abline(h=1.0, col="red", lty=2, lwd=1)
categories <- sort(cutree(p$tree_row, h=1))
categories <- as.data.frame(cbind(names(categories), categories))
write.csv(categories,"DEG_category1_subcategory.csv", quote = F)

for(i in unique(categories$categories)){
  pheatmap(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]), 
           scale = "row", show_rownames=F, cluster_cols = F, border_color=NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main=paste0("DEGs in Category1_", i),
           filename = paste0("Fig6_heatmap_DEG_C1_",i, ".pdf"), width=4, height=4)
  print(
    dim(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]))
  )
}


## GO analysis
library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

data_background <- toTable(org.Mm.egSYMBOL)

test1 = bitr(categories[categories$categories==1,1], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, paste0("GO_DEG_C1.csv"), quote=F)

test1 = bitr(categories[categories$categories==2,1], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, paste0("GO_DEG_C2.csv"), quote=F)


###################################### extract SMC #####################################
SMC <- subset(integrated, subset=majorCelltype=="SMC")
saveRDS(SMC, "SMC_from_integrated.rds")

SMC <- readRDS("SMC_from_integrated.rds") #../

meta <- SMC@meta.data
meta <- meta %>% mutate(Group=case_when(
  Sample %in% c("Saline.GFP", "Saline.RFP") ~ "Control",
  Sample %in% c("nonIMH.GFP", "nonIMH.RFP") ~ "nonIMH",
  Sample %in% c("IMH.GFP", "IMH.RFP") ~ "IMH"
))
SMC$Group <- meta$Group
SMC$Group <- factor(SMC$Group, levels=c("Control","nonIMH","IMH"))

SMC$Sample <- factor(SMC$Sample, levels=c("Saline.GFP", "Saline.RFP", "nonIMH.GFP", "nonIMH.RFP", "IMH.GFP", "IMH.RFP"))

varGene <- FindVariableFeatures(SMC)
varGenes <- varGene@assays$RNA@var.features
write.csv(varGenes, "SMC_variableGenes.csv",quote = F)

#SMC.GFP <- subset(SMC, subset=Sample %in% c("Saline.GFP", "nonIMH.GFP", "IMH.GFP"))
DotPlot(SMC, features = c("Tmem173", "Irf3","Ezh2"), group.by = "Sample")+coord_flip()+RotatedAxis()

FeaturePlot(SMC, features = c("Tnfrsf11b","Lum","Dcn","Col3a1"), raster=T)

#highlight cluster3
SMC$SMCcluster <- ifelse(SMC$seurat_clusters==3, "3", NA)
DimPlot(SMC, group.by = "SMCcluster", raster=T, label = T)

################################ define SMC clusters ######################
Idents(SMC) <- "seurat_clusters"
markers <- FindAllMarkers(SMC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "AllMarker_SMCcluster_mTmG_SalineAngII.csv", quote = F)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DotPlot(SMC, features = unique(top10$gene), cols=c("lightgrey", "darkblue"))+
  RotatedAxis()+ggtitle("Top10 integrate RNA")

SMC$seurat_clusters <- factor(SMC$seurat_clusters, levels = c(1,7,4,6,3,8,9))
Idents(SMC) <- "seurat_clusters"
DotPlot(SMC, features = c("Myh11","Rock1","Tnnt2","Des","Col4a6","Acta2","Myl9","Tpm2","Cnn1",
                         "Cmss1","S100a4","Lgals1","Mustn1", "Spp1","Dcn","Lum","Serpinf1", "Tnfrsf11b", 
                         "Gpc6","Neat1", "Nfkb1",
                         "Atf3","Fos","Junb","Egr1","Klf2"),
        cols=c("lightgrey", "mediumblue"))+coord_flip()
  

VlnPlot(SMC, features = "nCount_RNA", pt.size = 0)+NoLegend()

## GO of marker genes
data_background <- toTable(org.Mm.egSYMBOL)

markers <- read.csv("AllMarker_SMCcluster_mTmG_SalineAngII.csv", header=T, row.names = "X", stringsAsFactors = F)
for(i in unique(markers$cluster)){
  each <- markers[markers$cluster==i,]
  test1 = bitr(each$gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                         qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
  each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
  write.csv(each_up_GO, paste0("GO_allMarkers_SMCcluster_",i,".csv"), quote=F)
}

## heatmap of selected GO
files <- list.files(path=".", pattern = "^GO_allMarkers_SMCcluster.*\\.csv$")
dfs <- lapply(files, function(x) read.csv(x, header=T, stringsAsFactors=F))

top10 <- lapply(dfs, function(x) x[1:10,4])
top10 <- unique(unlist(top10))

data <- c()
for(i in 1:length(dfs)){
  each <- dfs[[i]]
  each <- each[each$Description %in% top10, c(4,9)]
  each$Cluster <- str_split_fixed(files[i], "_", 3)[,3]
  each$Cluster <- substr(each$Cluster, 1, nchar(each$Cluster)-4)
  data <- rbind(data, each)
}
data$qvalue_log <- -log10(data$qvalue)
data <- data[,c(1,3,4)]
data <- reshape::cast(data, Description~Cluster)
row.names(data) <- data$Description
data <- data[,-1]
data[is.na(data)] = 0

pheatmap(data, scale = "row", color = colorRampPalette(c("blue", "white", "red"))(50))

#
selected <- read.table("clipboard", header=F, sep="\t", stringsAsFactors = F)
data <- c()
for(i in 1:length(dfs)){
  each <- dfs[[i]]
  each <- each[each$Description %in% selected$V2, c(4,9)]
  each$Cluster <- str_split_fixed(files[i], "_", 3)[,3]
  each$Cluster <- substr(each$Cluster, 1, nchar(each$Cluster)-4)
  data <- rbind(data, each)
}
data$qvalue_log <- -log10(data$qvalue)
data <- data[,c(1,3,4)]
data <- reshape::cast(data, Description~Cluster)
row.names(data) <- data$Description
data <- data[,-1]
data[is.na(data)] = 0
data <- data[match(selected$V2, row.names(data)),c(1,5,3,4,2,6,7)]

pheatmap(data, scale = "row", cluster_cols = F, cluster_rows = F,
         color = colorRampPalette(c("blue", "white", "red"))(50))

#################################### compare RFP vs. GFP on each SMC clusters ####################################
SMC$Lineage <- ifelse(SMC$Sample %in% c("Saline.GFP", "nonIMH.GFP","IMH.GFP"), "GFP", "RFP")
Idents(SMC) <- paste0(SMC$seurat_clusters, "_", SMC$Lineage)
for (i in unique(SMC$seurat_clusters)){
  each <- FindMarkers(SMC, ident.1 = paste0(i,"_RFP"), ident.2 = paste0(i,"_GFP"))
  write.csv(each, paste0("wilcox_",i,"_RFPvsGFP.csv"), quote = F)
}


################################ SMC cluster proportion #######################
meta <- SMC@meta.data
meta1 <- meta[meta$Sample %in% c("Saline.GFP","nonIMH.GFP","IMH.GFP"),]
meta.pct <- meta1 %>% group_by(Group,seurat_clusters) %>% 
  summarise(Percentage=n()) %>% 
  group_by(Group) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)
meta.pct$Group <- factor(meta.pct$Group, levels=c("Control","nonIMH","IMH"))
p1 <- ggplot(meta.pct, aes(x=seurat_clusters, y=Percentage, fill=Group))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+xlab("Cluster")+ ggtitle("GFP cluster proportion")+
  scale_fill_manual(values=c("darkgrey","darkgoldenrod1","firebrick3"))

meta2 <- meta[meta$Sample %in% c("Saline.RFP", "nonIMH.RFP", "IMH.RFP"),] 
meta.pct <- meta2 %>% group_by(Group,seurat_clusters) %>% 
  summarise(Percentage=n()) %>% 
  group_by(Group) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)
meta.pct$Group <- factor(meta.pct$Group, levels=c("Control","nonIMH","IMH"))
p2 <- ggplot(meta.pct, aes(x=seurat_clusters, y=Percentage, fill=Group))+
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()+xlab("Cluster")+ ggtitle("RFP cluster proportion")+
  scale_fill_manual(values=c("darkgrey","darkgoldenrod1","firebrick3"))

p1+p2

## cluster proportion chi-square test of goodness-of-fit
meta.count <- meta1 %>% group_by(Sample,seurat_clusters) %>% 
  summarise(Count=n())
meta.count <- reshape::cast(meta.count, seurat_clusters~Sample) 

# nonIMH vs. control
for(i in 1:dim(meta.count)[1]){
  each1 <- meta.count[i,2:4]
  each2 <- colSums(meta.count)-meta.count[i,2:4]
  each <- rbind(each1, each2)
  test <- chisq.test(each[,2], 
                     p = c(each[1,3]/colSums(meta.count)[3],each[2,3]/colSums(meta.count)[3]))
  print(meta.count[i,1])
  print(test$p.value*7) # Bonferroni correction of p Value
}

# IMH vs nonIMH
for(i in 1:dim(meta.count)[1]){
  each1 <- meta.count[i,2:4]
  each2 <- colSums(meta.count)-meta.count[i,2:4]
  each <- rbind(each1, each2)
  test <- chisq.test(each[,1], 
                     p = c(each[1,2]/colSums(meta.count)[2],each[2,2]/colSums(meta.count)[2]))
  print(meta.count[i,1])
  print(test$p.value*7) # Bonferroni correction of p Value
}

## two-proportions z-test
# nonIMH vs. control
for(i in 1:dim(meta.count)[1]){
  each1 <- meta.count[i,2:4]
  each2 <- colSums(meta.count)-meta.count[i,2:4]
  each <- rbind(each1, each2)
  test <- prop.test(x = each[,3], n = each[,2])
  print(meta.count[i,1])
  print(test$p.value*7) # Bonferroni correction of p Value
}


################################################ visulize SMC gene expression ################################
avg <- AverageExpression(SMC, assays = "RNA", group.by = "Sample")
avg <- avg$RNA
write.csv(avg, "avg_SMC_Sample.csv", quote=F)
avg <- read.csv("avg_SMC_Sample.csv",header=T, row.names = "X", stringsAsFactors = F)

## violin plot
genes <- c("Myh11","Myl9","Mylk", "Tagln","Acta2","Cnn1","Tpm2","Lmod1",
           "Il1b","Ccl2","Nfkb1","Spp1","Cxcl1","Cxcl12","Il6","Cxcl10")
genes <- c("Batf3","Cic", "Esr2", "Lyl1",
           "Nfia", "Rest","Tal1") #"Ascl1","Myod1","Myog","Myf5",

i=1
p <- VlnPlot(SMC, features =genes[i], group.by = "Sample", pt.size = 0, 
             )+ #cols = c("darkgrey","darkgoldenrod1","firebrick3")
  ggthemes::theme_few()+
  stat_summary(fun.y=median, geom="point", size=2, color="black")
for(i in 2:length(genes)){
  p1 <- VlnPlot(SMC, features = genes[i], group.by = "Sample", pt.size = 0, 
                )+ #cols = c("darkgrey","darkgoldenrod1","firebrick3")
    ggthemes::theme_few()+
    stat_summary(fun.y=median, geom="point", size=2, color="black")
  
  p <- p+p1
  p + patchwork::plot_layout(ncol = 6)
}
p + patchwork::plot_layout(ncol = 6)
ggsave("Fig3_Vln_TFs_Sample.pdf", width = 30, height = 6)

## dotplot
gene1 <- c("Myh11","Myl9","Mylk", "Tagln","Acta2","Tpm2","Lmod1")
gene2 <- c("Nfkb1","Spp1","Cxcl1","Cxcl12","Il1b","Ccl2","Il6","Cxcl10")
DotPlot(SMC, features = gene1, group.by = "Group")+coord_flip()+RotatedAxis()
DotPlot(SMC, features = gene2, group.by = "Group")+coord_flip()+RotatedAxis()

## heatmap
d <- avg[row.names(avg) %in% genes, c(5,3,1)] #,6,4,2
pheatmap(d, scale="row",show_rownames = T, cluster_cols = F, cutree_rows = 2,gaps_col = 3,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main=paste0("SMC contractile inflammatory"))

## feature plot
FeaturePlot(SMC, features = c("Spp1", "Cd44", "Itga8","Itgb1"), split.by = "Sample", raster = T)
FeaturePlot(SMC, features = genes, split.by = "Sample", raster = T)
ggsave("Fig4_featurePlot_TFs_Sample.pdf", width = 16, height = 18)


############################################## SMC_GFP only ###########################
SMC.GFP <- subset(SMC, subset=Lineage=="GFP")

## violin plot
genes <- c("Myh11","Myl9","Mylk", "Tagln","Acta2","Cnn1","Tpm2","Lmod1",
           "Il1b","Ccl2","Nfkb1","Spp1","Cxcl1","Cxcl12","Il6","Cxcl10")

i=1
p <- VlnPlot(SMC.GFP, features =genes[i], group.by = "Group", pt.size = 0, 
             cols = c("darkgrey","darkgoldenrod1","firebrick3"))+
  ggthemes::theme_few()+
  stat_summary(fun.y=median, geom="point", size=2, color="black")
for(i in 2:length(genes)){
  p1 <- VlnPlot(SMC.GFP, features = genes[i], group.by = "Group", pt.size = 0, 
                cols = c("darkgrey","darkgoldenrod1","firebrick3"))+
    ggthemes::theme_few()+
    stat_summary(fun.y=median, geom="point", size=2, color="black")
  
  p <- p+p1
  p + patchwork::plot_layout(ncol = 6)
}
p + patchwork::plot_layout(ncol = 6)
ggsave("Fig3_Vln_contractile_inflammatory_GFP.pdf", width = 20, height = 8)

## differntial
Idents(SMC.GFP) <- paste0(SMC.GFP$majorCelltype, "_", SMC.GFP$Group)
for(i in unique(SMC.GFP@meta.data$majorCelltype)){
  each <- FindMarkers(SMC.GFP, ident.1 = paste0(i,"_IMH"), ident.2 = paste0(i,"_nonIMH"))
  write.csv(each, paste0("wilcox_GFP_",i,"_IMHvsnonIMH.csv"), quote = F)
}
for(i in unique(SMC.GFP@meta.data$majorCelltype)){
  each <- FindMarkers(SMC.GFP, ident.1 = paste0(i,"_nonIMH"), ident.2 = paste0(i,"_Control"))
  write.csv(each, paste0("wilcox_GFP_",i,"_nonIMHvscontrol.csv"), quote = F)
}

## elatic fiber genes
genes <- read.table("Z:/CT Surgery Lab/Yanming Li/AHA_scRNAseq/human/6.4_DissectionAneurysmControl/SMC/marker_RNAmodification.txt", 
header=F, sep="\t", stringsAsFactors = F)
genes <- read.table("marker_ElasticFiber.txt", header=F, sep="\t", stringsAsFactors = F)
genes <- genes[,1]
library(nichenetr)
genes = genes %>% convert_human_to_mouse_symbols()
genes <-unname(genes)

SMC.GFP$Sample <- factor(SMC.GFP$Sample, levels=c("Saline.GFP", "nonIMH.GFP", "IMH.GFP"))
DotPlot(SMC.GFP, features = rev(genes), cols=c("lightgrey", "mediumblue"), group.by = "Sample")+
  coord_flip()+RotatedAxis()


#################################### SMC_GFP DEGs, classified as categories, GO functions ###########################################
files <- list.files(path=".", pattern = "wilcox_GFP_SMC_*")
dfs <- lapply(files, function(x) read.csv(x, header=T, stringsAsFactors=F, row.names = "X"))
dfs <- lapply(dfs, function(x) x[which(x$p_val_adj<0.05 & abs(x$avg_log2FC)>0),])
DEGs <- unique(unlist(lapply(dfs, row.names)))
avg <- read.csv("avg_SMC_Sample.csv", header=T, row.names = "X", stringsAsFactors = F)

setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/SMC/")
d <- avg[row.names(avg) %in% DEGs, c(5,3,1)]
p <- pheatmap(d, scale="row",show_rownames = F, cluster_cols = F, cutree_rows = 8,
              color = colorRampPalette(c("blue", "white", "red"))(50),
              main=paste0("SMC_GFP DEGs"))

plot(p$tree_row)
abline(h=1.35, col="red", lty=2, lwd=1)
categories <- sort(cutree(p$tree_row, h=1.35))
categories <- as.data.frame(cbind(names(categories), categories))
write.csv(categories,"DEG_GFP_category.csv", quote = F)
categories <- read.csv("DEG_GFP_category.csv", header = T, stringsAsFactors = F, row.names = "X")

for(i in unique(categories$categories)){
  pheatmap(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]), 
           scale = "row", show_rownames=F, cluster_cols = F, border_color=NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main=paste0("GFP DEGs in Category", i),
           filename = paste0("Fig6_heatmap_GFP_DEG_C",i, ".pdf"), width=4, height=4)
  print(
    dim(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]))
  )
}

## re-organize each categories
categories$categories <- factor(categories$categories, level=c(3,2,5,8,1,4,6,7))
categ2 <- categories %>% arrange(categories)

d2 <- d[match(categ2$V1, row.names(d)),]
pheatmap(d2, scale="row", cluster_rows = F, cluster_cols = F, show_rownames = F,
         gaps_row = c(dim(categ2[categ2$categories %in% c(3,2),])[1],
                      dim(categ2[categ2$categories %in% c(3,2,5,8),])[1],
                      dim(categ2[categ2$categories %in% c(3,2,5,8,1),])[1],
                      dim(categ2[categ2$categories %in% c(3,2,5,8,1,4),])[1]),
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main="SMC GFP DEGs re-organized")

## GO analysis
data_background <- toTable(org.Mm.egSYMBOL)

test1 = bitr(categories[categories$categories %in% c(6,7),1], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, paste0("GO_GFP_DEG_C67.csv"), quote=F)

test1 = bitr(categories[categories$categories==1,1], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, paste0("GO_GFP_DEG_C1.csv"), quote=F)

test1 = bitr(categories[categories$categories==4,1], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, paste0("GO_GFP_DEG_C4.csv"), quote=F)

## barplot of selected GO
GOs <- read.table("clipboard", header=T, sep="\t", stringsAsFactors = F)
GOs$Description <- factor(GOs$Description, levels=rev(GOs$Description))

ggplot(GOs, aes(x=Description, y=qvalue_converted, fill=Category))+geom_bar(stat="identity", width=0.75)+
  coord_flip()+theme_classic()+ylab("-log10(q value)")+ggtitle("GFP DEG category")+
  scale_fill_manual(values=c("orange","darkgreen","firebrick3","firebrick3"))


########################################## SMC_RFP DEGs, as categories #########################################
## differential
Idents(SMC) <- paste0(SMC$majorCelltype, "_", SMC$Sample)

each <- FindMarkers(SMC, ident.1 = "SMC_IMH.RFP", ident.2 = "SMC_nonIMH.RFP")
write.csv(each, paste0("wilcox_RFP_SMC_IMHvsnonIMH.csv"), quote = F)

each <- FindMarkers(SMC, ident.1 = "SMC_nonIMH.RFP", ident.2 = "SMC_Saline.RFP")
write.csv(each, paste0("wilcox_RFP_SMC_nonIMHvsControl.csv"), quote = F)

## set as categories
files <- list.files(path=".", pattern = "wilcox_RFP_SMC_*")
dfs <- lapply(files, function(x) read.csv(x, header=T, stringsAsFactors=F, row.names = "X"))
dfs <- lapply(dfs, function(x) x[which(x$p_val_adj<0.05 & abs(x$avg_log2FC)>0),])
DEGs <- unique(unlist(lapply(dfs, row.names)))
avg <- read.csv("avg_SMC_Sample.csv", header=T, row.names = "X", stringsAsFactors = F)

setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/SMC/")
d <- avg[row.names(avg) %in% DEGs, c(6,4,2)]
p <- pheatmap(d, scale="row",show_rownames = F, cluster_cols = F, cutree_rows = 4,
              color = colorRampPalette(c("blue", "white", "red"))(50),
              main=paste0("SMC_RFP DEGs"))

plot(p$tree_row)
abline(h=2.5, col="red", lty=2, lwd=1)
categories <- sort(cutree(p$tree_row, h=2.5))
categories <- as.data.frame(cbind(names(categories), categories))
write.csv(categories,"DEG_RFP_category.csv", quote = F)

for(i in unique(categories$categories)){
  pheatmap(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]), 
           scale = "row", show_rownames=F, cluster_cols = F, border_color=NA,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main=paste0("RFP DEGs in Category", i),
           filename = paste0("Fig6_heatmap_RFP_DEG_C",i, ".pdf"), width=4, height=4)
  print(
    dim(as.matrix(d[row.names(d)%in% categories[categories$categories==i,1],]))
  )
}

## GO analysis
data_background <- toTable(org.Mm.egSYMBOL)
for (i in unique(categories$categories)){
  test1 = bitr(categories[categories$categories==i,1], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
  each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                         qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
  each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
  write.csv(each_up_GO, paste0("GO_RFP_DEG_C",i,".csv"), quote=F)
}






######################################### trajectory of GFP+ cells #############################################
setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_Challenge/WithoutAngII3d/GFPtrajectory/")
GFP <- subset(integrated, subset = orig.ident %in% c("CMG", "DMG"))
DimPlot(GFP, group.by = "seurat_clusters", label=T, raster=T)


library(monocle3)

cds <- readRDS("SMC_monocle3_cds.rds") 
branch1 <- readRDS("SMC_monocle3_branch1_cds.rds")
branch2 <- readRDS("SMC_monocle3_branch2_cds.rds")
branch3 <- readRDS("SMC_monocle3_branch3_cds.rds")

# analysis
expression_matrix<- GFP@assays$RNA@counts
cell_metadata <- GFP@meta.data
gene_annotation <- data.frame("gene_short_name"=row.names(expression_matrix),
                              "row.names"=row.names(expression_matrix))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

## Step 2: Remove batch effects with cell alignment
#cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by="seurat_clusters")
plot_cells(cds, color_cells_by="Group")
plot_cells(cds, color_cells_by="orig.ident")
# no need for batch correction
#cds = align_cds(cds, num_dim = 100, alignment_group = "orig.ident")
#cds = reduce_dimension(cds)
p1 <- plot_cells(cds, color_cells_by="seurat_clusters")+ggtitle("clustered by Seurat")
p3 <- plot_cells(cds, color_cells_by="orig.ident")+ggtitle("Colored by orig.ident")
plot_cells(cds, color_cells_by="orig.ident")

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
plot_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds, reduction_method = "UMAP")
cds@colData@listData$Group <- factor(cds@colData@listData$orig.ident, levels=c("CMG","DMG"))

saveRDS(cds, "GFP_monocle3_cds.rds")


p2 <- plot_cells(cds, color_cells_by = "pseudotime")+ggtitle("Pseudotime by Monocle3")
p4 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = F)+
  ggtitle("Pseudotime by Monocle3")
p1+p2+p3+p4

p1 <- plot_cells(cds[,colData(cds)$Group=="Control"], color_cells_by = "pseudotime", show_trajectory_graph = F)+ggtitle("Control")+ylim(-5,7.5)
p2 <- plot_cells(cds[,colData(cds)$Group=="AngII_nonDissection"], color_cells_by = "pseudotime", show_trajectory_graph = F)+ggtitle("AngII_nonDissection")
p3 <- plot_cells(cds[,colData(cds)$Group=="AngII_Dissection"], color_cells_by = "pseudotime", show_trajectory_graph = F)+ggtitle("AngII_Dissection")
p1+p2+p3

## define branches
branch1 <- choose_graph_segments(cds)
saveRDS(branch1, "SMC_monocle3_branch1_cds.rds")
branch2 <- choose_graph_segments(cds)
saveRDS(branch2, "SMC_monocle3_branch2_cds.rds")
branch3 <- choose_graph_segments(cds)
saveRDS(branch3, "SMC_monocle3_branch3_cds.rds")


colData(cds)$Branch1 <- ifelse(colnames(cds) %in% colnames(branch1), cds@principal_graph_aux$UMAP$pseudotime, Inf)
colData(cds)$Branch2 <- ifelse(colnames(cds) %in% colnames(branch2), cds@principal_graph_aux$UMAP$pseudotime, Inf)
colData(cds)$Branch3 <- ifelse(colnames(cds) %in% colnames(branch3), cds@principal_graph_aux$UMAP$pseudotime, Inf)
p1 <- plot_cells(cds, color_cells_by = "Branch1", show_trajectory_graph = F, label_cell_groups=FALSE)+
  ggtitle("Branch1")+viridis::scale_color_viridis(option="viridis")
p2 <- plot_cells(cds, color_cells_by = "Branch2", show_trajectory_graph = F, label_cell_groups=FALSE)+
  ggtitle("Branch2")+viridis::scale_color_viridis(option="viridis")
p3 <- plot_cells(cds, color_cells_by = "Branch3", show_trajectory_graph = F, label_cell_groups=FALSE)+
  ggtitle("Branch3")+viridis::scale_color_viridis(option="viridis")
p1+p2+p3

p1 <- plot_cells(cds[,colData(cds)$Group=="Control"], color_cells_by = "Branch3", show_trajectory_graph = F)+ggtitle("Control")+
  viridis::scale_color_viridis(option="viridis")
p2 <- plot_cells(cds[,colData(cds)$Group=="AngII_nonDissection"], color_cells_by = "Branch3", show_trajectory_graph = F)+ggtitle("AngII_nonDissection")+
  viridis::scale_color_viridis(option="viridis")
p3 <- plot_cells(cds[,colData(cds)$Group=="AngII_Dissection"], color_cells_by = "Branch3", show_trajectory_graph = F)+ggtitle("AngII_Dissection")+
  viridis::scale_color_viridis(option="viridis")
p1+p2+p3

summary(cds@principal_graph_aux$UMAP$pseudotime)

## gene expression on a branch in UMAP
colData(cds)$ACTA2 <- cds@assays@data$counts[row.names(cds@assays@data$counts)=="ACTA2",]
colData(cds)$ACTA2 <- ifelse(colnames(cds) %in% colnames(branch1), colData(cds)$ACTA2, Inf)
p1 <- plot_cells(cds, color_cells_by = "ACTA2", show_trajectory_graph = F, label_cell_groups=FALSE)+
  ggtitle("ACTA2")+viridis::scale_color_viridis(option="turbo")

colData(cds)$IL6 <- cds@assays@data$counts[row.names(cds@assays@data$counts)=="IL6",]
colData(cds)$IL6 <- ifelse(colnames(cds) %in% colnames(branch1), colData(cds)$IL6, Inf)
p2 <- plot_cells(cds, color_cells_by = "IL6", show_trajectory_graph = F, label_cell_groups=FALSE)+
  ggtitle("IL6")+viridis::scale_color_viridis(option="turbo")

p1+p2

## cell density along pseudotime by branch 
out2 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
colnames(out2) <- "Pseudotime"
out2$Group <- cds@colData$Group
out2 <- out2 %>% mutate(Branches=case_when(
  row.names(out2) %in% colnames(branch2) ~ "Branch2",
  row.names(out2) %in% colnames(branch1) ~ "Branch1",
  TRUE ~ "Others"
))

ggplot(out2, aes(x=Pseudotime, color=Group))+geom_density(lwd=1)+
  theme_classic()+
  scale_color_manual(values=c("darkgrey","darkgoldenrod1","firebrick3","magenta4"))+
  facet_grid(.~Branches)

## gene expression by time
data <- SMC[,colnames(SMC) %in% colnames(cds)]
data <- NormalizeData(data) #one time thing

out2 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
colnames(out2) <- "Pseudotime"
out2 <- out2[row.names(out2) %in% colnames(branch3),,drop=F]

genes <- c("Acta2","Tagln","Myh11",
           "Eln","Fn1","Lum","Col8a1","Col1a2",  
           "Atf3","Hif1a","Myc","Jun",
           "Ccl2","Il6","Spp1","S100a4") 
gene_activity <- t(data@assays$RNA@data[row.names(data@assays$RNA@data) %in% genes,])
gene_activity <- gene_activity[row.names(gene_activity) %in% row.names(out2),]

out2 <- cbind(out2, gene_activity)
out2$Cluster <- data@meta.data[row.names(data@meta.data)%in%row.names(out2), ]$seurat_clusters
out2 <- reshape2::melt(out2, id=c("Pseudotime", "Cluster"))
colnames(out2)[3:4] <- c("Gene", "Expression")
out2$Gene <- factor(out2$Gene, levels=genes)

ggplot(out2, aes(x=Pseudotime, y=Expression, colour=Pseudotime))+ggrastr::geom_point_rast(size=0.5)+
  facet_wrap(.~Gene, ncol=8)+viridis::scale_color_viridis(option="viridis")+ #"C"
  geom_smooth(se=F, colour="black")+theme_classic()+ggtitle("Branch3")



ggplot(out2, aes(x=Pseudotime, y=Expression, colour=Cluster))+geom_point(size=0.5)+
  facet_wrap(.~Gene, ncol=8)+
  geom_smooth(se=F, colour="black")+theme_classic()+ggtitle("Branch1")


##################################### Macrophage in the background of all cells ###########################
integrated$Maph2 <- ifelse(integrated$seurat_clusters %in% c(2,14,15), 
                          as.character(integrated$seurat_clusters), NA)
integrated$Maph2 <- factor(integrated$Maph2, levels=c(2,14,15,NA))
DimPlot(integrated, group.by = "Maph2", raster=T, label=T)
DimPlot(integrated, group.by = "Maph2", raster=T, label=T, split.by = "Sample")
DimPlot(integrated, group.by = "Maph2", raster=T, label=T, split.by = "Lineage")

# Macrophage proportion in GFP cells
meta <- integrated@meta.data
meta1 <- meta[meta$Sample %in% c("Saline.GFP","nonIMH.GFP","IMH.GFP"),]
meta.pct <- meta1 %>% group_by(Group,majorCelltype) %>% 
  summarise(Percentage=n()) %>% 
  group_by(Group) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)
meta.pct$Group <- factor(meta.pct$Group, levels=c("Control","nonIMH","IMH"))
ggplot(meta.pct[meta.pct$majorCelltype=="Macrophage",], aes(x=majorCelltype, y=Percentage, fill=Group))+
  geom_bar(stat="identity", position=position_dodge(), width = 0.8)+
  theme_classic()+xlab("Cluster")+ ggtitle("Macrophage in GFP cells")+
  scale_fill_manual(values=c("darkgrey","darkgoldenrod1","firebrick3"))


################################# extract macrophage #######################################
Maph <- subset(integrated, subset=majorCelltype=="Macrophage")
saveRDS(Maph, "Maph_from_integrated.rds")

Maph <- readRDS("../Maph_from_integrated.rds")

DimPlot(Maph, group.by = "seurat_clusters", split.by = "Sample")

############################ compare GFP Maph Vs. RFP Maph ##################################
setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/Macrophage/")
meta <- Maph@meta.data
meta <- meta %>% mutate(Lineage=case_when(
  Sample %in% c("Saline.GFP", "nonIMH.GFP","IMH.GFP") ~ "GFP",
  Sample %in% c("Saline.RFP", "nonIMH.RFP", "IMH.RFP") ~ "RFP"
))
Maph$Lineage <- meta$Lineage

Idents(Maph) <- paste0(Maph$majorCelltype, "_", Maph$Lineage)
each <- FindMarkers(Maph, ident.1 = paste0("Macrophage_GFP"), ident.2 = paste0("Macrophage_RFP"))
write.csv(each, paste0("wilcox_Maph_GFPvsRFP.csv"), quote = F)

## GO and KEGG 
library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

data_background <- toTable(org.Mm.egSYMBOL)

each_up <- each[which(each$avg_log2FC>0 & each$p_val_adj < 0.05),]
test1 = bitr(row.names(each_up), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                         qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, paste0("GO_up_Maph_GFPvsRFP.csv"), quote=F)
  #each_up_KEGG <- enrichKEGG(test1$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
  #                           qvalueCutoff = 0.2)
  #each_up_KEGG@result$Description <- str_remove(each_up_KEGG@result$Description, ",")
  #write.csv(each_up_KEGG, paste0("KEGG_allDEGs_", i, ".csv"),  quote=F)
each_down <- each[which(each$avg_log2FC<0 & each$p_val_adj < 0.05),]
test1 = bitr(row.names(each_down), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_down_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_down_GO@result$Description <- str_remove(each_down_GO@result$Description, ",")
write.csv(each_down_GO, paste0("GO_down_Maph_GFPvsRFP.csv"), quote=F)

## visulize GO
GOs <- read.table("clipboard", header=T, sep="\t", stringsAsFactors = F)
GOs$Description <- factor(GOs$Description, levels=rev(GOs$Description))
GOs$qvalue_converted_2 <- ifelse(GOs$Trend=="GFP_Down", GOs$qvalue_converted*-1, GOs$qvalue_converted)

ggplot(GOs, aes(x=Description, y=qvalue_converted_2, fill=Trend))+geom_bar(stat="identity", width=0.85)+
  coord_flip()+theme_classic()+ylab("-log10(q value)")+ggtitle("GFP.Marophage vs. RFP.Macrophage")+
  scale_fill_manual(values=c("LightSeaGreen","Salmon"))



############################## Fibroblast in the background of all cells ########################
integrated$FB2 <- ifelse(integrated$seurat_clusters %in% c(0,5,10,11,16), 
                           as.character(integrated$seurat_clusters), NA)
integrated$FB2 <- factor(integrated$FB2, levels=c(0,5,10,11,16,NA))
DimPlot(integrated, group.by = "FB2", raster=T, label=T)
DimPlot(integrated, group.by = "FB2", raster=T, label=T, split.by = "Sample")
DimPlot(integrated, group.by = "FB2", raster=T, label=T, split.by = "Lineage")

# FB proportion in GFP cells
meta <- integrated@meta.data
meta1 <- meta[meta$Sample %in% c("Saline.GFP","nonIMH.GFP","IMH.GFP"),]
meta.pct <- meta1 %>% group_by(Group,majorCelltype) %>% 
  summarise(Percentage=n()) %>% 
  group_by(Group) %>% 
  mutate(Percentage=Percentage/sum(Percentage)*100)
meta.pct$Group <- factor(meta.pct$Group, levels=c("Control","nonIMH","IMH"))
ggplot(meta.pct[meta.pct$majorCelltype %in% c("Macrophage", "FB"),], aes(x=majorCelltype, y=Percentage, fill=Group))+
  geom_bar(stat="identity", position=position_dodge(), width = 0.8)+
  theme_classic()+xlab("Cluster")+ ggtitle("Proportion in GFP cells")+
  scale_fill_manual(values=c("darkgrey","darkgoldenrod1","firebrick3"))

## cluster proportion Fisher's test
meta.count <- meta1 %>% group_by(Sample,majorCelltype) %>% 
  summarise(Count=n())
meta.count <- reshape::cast(meta.count, majorCelltype~Sample) 
meta.count <- meta.count[1:6,]

# nonIMH vs. control
for(i in 1:dim(meta.count)[1]){
  each1 <- meta.count[i,2:3]
  each2 <- colSums(meta.count[,2:4])[1:2]-meta.count[i,2:3]
  each <- rbind(each1, each2)
  test <- fisher.test(each)
  print(meta.count[i,1])
  print(test$p.value*6)
}

# IMH vs nonIMH
for(i in 1:dim(meta.count)[1]){
  each1 <- meta.count[i,3:4]
  each2 <- colSums(meta.count[,2:4])[2:3]-meta.count[i,3:4]
  each <- rbind(each1, each2)
  test <- fisher.test(each)
  print(meta.count[i,1])
  print(test$p.value*6)
}

# IMH vs control
for(i in 1:dim(meta.count)[1]){
  each1 <- meta.count[i,c(2,4)]
  each2 <- colSums(meta.count[,2:4])[c(1,3)]-meta.count[i,c(2,4)]
  each <- rbind(each1, each2)
  test <- fisher.test(each)
  print(meta.count[i,1])
  print(test$p.value*6)
}


############################ compare GFP FB Vs. RFP FB ##################################
setwd("Z:/CT Surgery Lab/Yanming Li/mouse_scRNAseq/mTmG_basal_angII/Fibroblast/")
FB <- subset(integrated, subset=majorCelltype=="FB")

Idents(FB) <- paste0(FB$majorCelltype, "_", FB$Lineage)
each <- FindMarkers(FB, ident.1 = paste0("FB_GFP"), ident.2 = paste0("FB_RFP"))
write.csv(each, paste0("wilcox_FB_GFPvsRFP.csv"), quote = F)

## GO and KEGG 
library(DOSE)
library(org.Mm.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

data_background <- toTable(org.Mm.egSYMBOL)

each_up <- each[which(each$avg_log2FC>0 & each$p_val_adj < 0.05),]
test1 = bitr(row.names(each_up), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_up_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                       qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_up_GO@result$Description <- str_remove(each_up_GO@result$Description, ",")
write.csv(each_up_GO, paste0("GO_up_FB_GFPvsRFP.csv"), quote=F)
#each_up_KEGG <- enrichKEGG(test1$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
#                           qvalueCutoff = 0.2)
#each_up_KEGG@result$Description <- str_remove(each_up_KEGG@result$Description, ",")
#write.csv(each_up_KEGG, paste0("KEGG_allDEGs_", i, ".csv"),  quote=F)
each_down <- each[which(each$avg_log2FC<0 & each$p_val_adj < 0.05),]
test1 = bitr(row.names(each_down), fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")
each_down_GO <- enrichGO(test1$ENTREZID,OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "ALL", 
                         pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = data_background$gene_id, 
                         qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE)
each_down_GO@result$Description <- str_remove(each_down_GO@result$Description, ",")
write.csv(each_down_GO, paste0("GO_down_FB_GFPvsRFP.csv"), quote=F)

## visulize GO
GOs <- read.table("clipboard", header=T, sep="\t", stringsAsFactors = F)
GOs$Description <- factor(GOs$Description, levels=rev(GOs$Description))
GOs$qvalue_converted_2 <- ifelse(GOs$Trend=="GFP_Down", GOs$qvalue_converted*-1, GOs$qvalue_converted)

ggplot(GOs, aes(x=Description, y=qvalue_converted_2, fill=Trend))+geom_bar(stat="identity", width=0.85)+
  coord_flip()+theme_classic()+ylab("-log10(q value)")+ggtitle("GFP.FB vs. RFP.FB")+
  scale_fill_manual(values=c("LightSeaGreen","Salmon"))



################## convert seurat object to loom ###################
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

SaveH5Seurat(integrated, filename = "logNormalize_combined_mTmG_SalineAngII.h5Seurat")
Convert("logNormalize_combined_mTmG_SalineAngII.h5Seurat", dest = "h5ad")


