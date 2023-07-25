
setwd("/home/yiweizhang/test/肾脏/CD45_filtered_doublet")
library(ggplot2)
library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)

scRNA <- Read10X(data.dir ="/home/yiweizhang/test/肾脏/CD45_raw_data/control/") 
scRNA <- CreateSeuratObject(scRNA, project = "control", min.cells = 3, min.features = 200)  


### 降维聚类
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),verbose = T)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:40

scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters()

#选取0.5
scRNA <- FindClusters(object = scRNA, resolution = 0.3)#七个的时候0.5

##过滤双细胞
# DoubletFinder查看·过滤双细胞
library(DoubletFinder)

#寻找最优pK值
#这是一个测试最佳参数的过程，运行速度慢
sweep.res.list <- paramSweep_v3(scRNA, PCs = pc.num, sct = F)
#使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值

##确定nExp_poi.adj
DoubletRate = ncol(scRNA)*8*1e-6 #更通用
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
homotypic.prop <- modelHomotypic(scRNA$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi <- round(DoubletRate*ncol(scRNA)) 
# 使用同源双细胞比例对计算的双细胞比例进行校正 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##进行操作
scRNA <- doubletFinder_v3(scRNA, PCs =pc.num, pN = 0.25, pK = pK_bcmvn, 
                          nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
DimPlot(scRNA, group.by = names(scRNA@meta.data)[8])

dim(scRNA)


##过滤doublet
control <- scRNA[, scRNA@meta.data[[names(scRNA@meta.data)[8]]] == "Singlet"]
dim(control) 

saveRDS(control , "control.rds")




scRNA <- Read10X(data.dir ="/home/yiweizhang/test/肾脏/CD45_raw_data/IgAN14/") 
scRNA <- CreateSeuratObject(scRNA, project = "IgAN14", min.cells = 3, min.features = 200)  


### 降维聚类
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),verbose = T)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:40

scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters()

#选取0.5
scRNA <- FindClusters(object = scRNA, resolution = 0.3)#七个的时候0.5

##过滤双细胞
# DoubletFinder查看·过滤双细胞
library(DoubletFinder)

#寻找最优pK值
#这是一个测试最佳参数的过程，运行速度慢
sweep.res.list <- paramSweep_v3(scRNA, PCs = pc.num, sct = F)
#使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值

##确定nExp_poi.adj
DoubletRate = ncol(scRNA)*8*1e-6 #更通用
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
homotypic.prop <- modelHomotypic(scRNA$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi <- round(DoubletRate*ncol(scRNA)) 
# 使用同源双细胞比例对计算的双细胞比例进行校正 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##进行操作
scRNA <- doubletFinder_v3(scRNA, PCs =pc.num, pN = 0.25, pK = pK_bcmvn, 
                          nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
DimPlot(scRNA, group.by = names(scRNA@meta.data)[8])

dim(scRNA)


##过滤doublet
IgAN14 <- scRNA[, scRNA@meta.data[[names(scRNA@meta.data)[8]]] == "Singlet"]
dim(IgAN14) 

saveRDS(IgAN14 , "IgAN14.rds")



scRNA <- Read10X(data.dir ="/home/yiweizhang/test/肾脏/CD45_raw_data/IgAN21/") 
scRNA <- CreateSeuratObject(scRNA, project = "IgAN21", min.cells = 3, min.features = 200)  


### 降维聚类
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),verbose = T)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:40

scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters()

#选取0.5
scRNA <- FindClusters(object = scRNA, resolution = 0.3)#七个的时候0.5

##过滤双细胞
# DoubletFinder查看·过滤双细胞
library(DoubletFinder)

#寻找最优pK值
#这是一个测试最佳参数的过程，运行速度慢
sweep.res.list <- paramSweep_v3(scRNA, PCs = pc.num, sct = F)
#使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值

##确定nExp_poi.adj
DoubletRate = ncol(scRNA)*8*1e-6 #更通用
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
homotypic.prop <- modelHomotypic(scRNA$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi <- round(DoubletRate*ncol(scRNA)) 
# 使用同源双细胞比例对计算的双细胞比例进行校正 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##进行操作
scRNA <- doubletFinder_v3(scRNA, PCs =pc.num, pN = 0.25, pK = pK_bcmvn, 
                          nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
DimPlot(scRNA, group.by = names(scRNA@meta.data)[8])

dim(scRNA)


##过滤doublet
IgAN21 <- scRNA[, scRNA@meta.data[[names(scRNA@meta.data)[8]]] == "Singlet"]
dim(IgAN21) 

saveRDS(IgAN21 , "IgAN21.rds")



scRNA <- Read10X(data.dir ="/home/yiweizhang/test/肾脏/CD45_raw_data/IgAN52/") 
scRNA <- CreateSeuratObject(scRNA, project = "IgAN52", min.cells = 3, min.features = 200)  


### 降维聚类
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData()
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),verbose = T)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:40

scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
scRNA <- FindNeighbors(scRNA, dims=pc.num) %>% FindClusters()

#选取0.5
scRNA <- FindClusters(object = scRNA, resolution = 0.3)#七个的时候0.5

##过滤双细胞
# DoubletFinder查看·过滤双细胞
library(DoubletFinder)

#寻找最优pK值
#这是一个测试最佳参数的过程，运行速度慢
sweep.res.list <- paramSweep_v3(scRNA, PCs = pc.num, sct = F)
#使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值

##确定nExp_poi.adj
DoubletRate = ncol(scRNA)*8*1e-6 #更通用
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
homotypic.prop <- modelHomotypic(scRNA$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi <- round(DoubletRate*ncol(scRNA)) 
# 使用同源双细胞比例对计算的双细胞比例进行校正 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

##进行操作
scRNA <- doubletFinder_v3(scRNA, PCs =pc.num, pN = 0.25, pK = pK_bcmvn, 
                          nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
DimPlot(scRNA, group.by = names(scRNA@meta.data)[8])

dim(scRNA)


##过滤doublet
IgAN52 <- scRNA[, scRNA@meta.data[[names(scRNA@meta.data)[8]]] == "Singlet"]
dim(IgAN52) 

saveRDS(IgAN52 , "IgAN52.rds")

