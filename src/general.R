##=======================================================================##
##====================Chapter01：创建Seurat对象并质控====================##
##=======================================================================##



setwd("/home/zhangyiwei/cd45/7_12_52week")
library(ggplot2)
library(Seurat)
library(harmony)
# library(SeuratWrappers)
library(tidyverse)
library(patchwork)
library(ggpubr)
#rm(list = ls())
## 批量读取数据
### 设置数据路径与样本名称

#control <- readRDS("/home/yiweizhang/test/肾脏/DoubletFinder/control.rds")
#controlorig <- 
#IgAN13_1 <- readRDS("/home/yiweizhang/test/肾脏/DoubletFinder/IgAN13_1.rds")


assays <- dir("/home/yiweizhang/test/肾脏/CD45_filtered_doublet/")
dir <- paste0("/home/yiweizhang/test/肾脏/CD45_filtered_doublet/", assays)

# 按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
samples_name = c("control","IgAN14" ,"IgAN21","IgAN52")

scRNAlist <- list()
for(i in 1:length(dir)){
  #不设置min.cells过滤基因会导致CellCyclescRNAring报错：
  #Insufficient data values to produce 24 bins.  
  scRNAlist[[i]] <- readRDS( dir[i])
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])  
  #计算线粒体基因比例
  if(T){    
    scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^mt-") 
  }
  #计算核糖体基因比例
  if(T){
    scRNAlist[[i]][["percent.rb"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^Rp[sl]")
  }
  # #计算红细胞基因比例
  if(T){
    HB.genes <- c("Hbegf","Hbb-bs","Hba-a2","Hba-a1","Hbb-bt","Hbq1b","Hbq1a")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]]))
    scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes)
  }
}


names(scRNAlist) <- samples_name
#system.time(save(scRNAlist, file = "Integrate/scRNAlist0.Rdata")) 
#system.time(saveRDS(scRNAlist, file = "scRNAlist0.rds"))

## 合并
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])

#查看过滤前细胞数目
table(scRNA$orig.ident)
saveRDS(scRNA,"beforeqc.rds")


### 设置质控标准
minGene=200
maxGene=5000
maxUMI=30000 
pctMT=50
pctHB=1

### 看看细胞数过滤的如何
scRNAafter <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                       nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)
##
table(scRNAafter$orig.ident)

### 换成scRNA往后走了
scRNA <- subset(scRNA, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)


scRNA$together <- recode(scRNAafter$orig.ident, 
                         "control" = "control", 
                         "IgAN14" = "IgAN14",
                         "IgAN21" = "IgAN21",
                         "IgAN52" = "IgAN52")








## 查看批次效应
### 降维聚类
# library(future)
# options(future.globals.maxSize = 80 * 1024^3) #将全局变量上限调至80G
# plan(multiprocess, workers = 30)
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),verbose = T)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:40

#scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
#scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters()

###Harmony整合
system.time({scRNA <- RunHarmony(scRNA, group.by.vars = "together")})
#scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:15)%>% 
#FindClusters(resolution = 0.5)
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims=pc.num) %>% FindClusters()
scRNA<- RunUMAP(scRNA, reduction = "harmony", dims = pc.num)








#FeaturePlot(scRNA,features = 'Cd3e')
#DimPlot(scRNA, label = T)

p <- DimPlot(scRNA, group.by = "orig.ident",raster=FALSE)
p
ggsave("UMAP_Samples.pdf", p, width = 8, height = 6)

p <- DimPlot(scRNA, group.by = "orig.ident", split.by = "together", ncol = 3,raster=FALSE)
p
ggsave("UMAP_Samples_Split.pdf", p, width = 18, height = 12)

saveRDS(scRNA, "scRNA.rds")







#看不同resolution对分群的影响
#scRNA <- FindClusters(object = scRNA, resolution = c(seq(.1,1.6,.2)))
#library(clustree)
#clustree(scRNA@meta.data, prefix = "RNA_snn_res.")
#dev.off()
#colnames(scRNA@meta.data)
#pdf(file="clustree.pdf",width=10,height=14)  #可视化不同resolution对分群的影响



#选取0.9
scRNA <- FindClusters(object = scRNA, resolution = 0.9)

marker <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#11.8  运行到此处
write.csv(marker,file = 'cluster_marker.csv')

saveRDS(scRNA, "scRNAresolution0.9noharmony.rds")

##################################tsne & umap分群########################
p <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = T, label.size = 3, raster=FALSE)
p
ggsave("UMAP_Clusters.pdf",p, width = 8, height = 6)
p1 = DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", label = T, label.size = 2)

ggsave("Clusters.pdf", p2, width = 12, height = 5)



###########################点图#################

#免疫细胞
Neutrophils <- c('Ly6g','S100a9','S100a8','Retnlg','Ifitm1')
Dendritic_cell <- c('Flt3','Itgax','H2-Ab1')
Macrophage <- c('C1qa','C1qb','Cd74','Adgre1')
Monocyte <- c('Csf1r','Itgam')
B_cell <- c('Cd19','Cd79a','Cd79b')
T_cell <- c('Cd3d','Cd3e','Cd3g')
NK_cell <- c('Klrb1c','Klra8','Nkg7')
resident <- c('Cd69','Runx3','Itga1','Itgae','Cxcr3','Il7r','Cx3cr1')
Fibroblast<- c('Cd209a','Clec10a','F13a1','Plac8')
plasmacell <- c('Igha', 'Iglc1','Ighg2b','Ighm')
mastcell <- c('Kit','Fcer1a')

#dotplot1
genelist1 <- list(T_cell, resident, NK_cell, Macrophage,
                  B_cell,
                  Neutrophils,plasmacell,mastcell )
names(genelist1) <- c('T cell','resident','NK cell', 'Macrophage','B cell', 
                      'Neutrophils','Plasma cell','Mast cell')
#开始画1！

DotPlot(scRNA, features =genelist1,dot.scale = 8,col.min = 0,group.by = "seurat_clusters")+theme_bw()+#去除背景，旋转图片  
  theme(panel.grid = element_blank()) +
  #  axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
  scale_color_gradient2( low = "darkblue", mid = "white", high = "darkred")+#颜色渐变设置  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

ggsave("Dotplotmerge1.pdf",width = 30, height = 6)

p <- DimPlot(scRNA, group.by = "seurat_clusters", split.by = "orig.ident", ncol = 3, label = T,raster=FALSE)
p

p <- DimPlot(scRNA, group.by = "celltype.main", split.by = "orig.ident", ncol = 3, label = T,raster=FALSE)
p

p <- DimPlot(scRNA, group.by = "seurat_clusters", ncol = 1, label = T,raster=FALSE)
p



scRNA$celltype.main <- recode(scRNA$seurat_clusters,
                              "2" = "T cell",
                              "21" = "T cell",
                              "5" = "T cell",
                              "10" = "T cell",
                              "14" = "T cell",
                              "24" = "NK cell",
                              "18" = "NK cell",
                              "26" = "NK cell",
                              "0" = "Macrophage",
                              "3" = "Macrophage",
                              "4" = "Macrophage",
                              "0" = "Macrophage",
                              "1" = "Neutrophil",
                              "31" = "Neutrophil",
                              "32" = "Neutrophil",
                              "17" = "Neutrophil",
                              "25" = "Neutrophil",
                              "7" = "B cell",
                              "22" = "B cell",
                              "13" = "Plasma cell",
                              "33" = "Mast cell"
)
table(scRNA$celltype.main, scRNA$together)
table(scRNA$together)
p <- DimPlot(scRNA, group.by = "celltype.main",  ncol = 1, label = T,raster=FALSE)
p

table(scRNA$orig.ident)

Idents(scRNA) <- scRNA$celltype.main 
Cellratio <- prop.table(table(Idents(scRNA), scRNA$together), margin = 2)


