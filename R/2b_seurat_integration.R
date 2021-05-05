library(tidyverse)
library(viridis)
library(Seurat)
library(Matrix)
library(clustree)

date = Sys.Date()
load(file.path("data", "prdata.RData"))

all <- prdata %>%
  CreateSeuratObject(min.cells = 10, min.features = 500)

all$Region <- case_when(
  grepl("Ret", colnames(all)) ~ "eye",
  T ~ "brain"
)

#mitochondrial gene number
all[["percent.mt"]] <-PercentageFeatureSet(all,pattern="^mt-")

#split object
all_list <- SplitObject(all, split.by = "Region")

#normalize dataset 
for (i in 1:length(all_list)) {
  all_list[[i]] <- SCTransform(all_list[[i]], verbose = FALSE)
}

all_features <- SelectIntegrationFeatures(object.list = all_list, nfeatures = 10000)
all_list <- PrepSCTIntegration(object.list = all_list, anchor.features = all_features, 
                                    verbose = FALSE)


all_anchors <- FindIntegrationAnchors(object.list = all_list, normalization.method = "SCT", 
                                           anchor.features = all_features, verbose = FALSE)
all_integrated <- IntegrateData(anchorset = all_anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

all_integrated <- RunPCA(all_integrated, verbose = FALSE)
all_integrated <- RunUMAP(all_integrated, dims = 1:30)
DimPlot(all_integrated, group.by = c("Region"), combine = FALSE)

all <- all_integrated

all<- all %>% 
  FindNeighbors(dims=1:30) %>%
  FindClusters(resolution=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(all)

ggsave("plots/overview-cluster-resolutions.pdf", useDingbats=F)

all<-FindClusters(all,resolution=.6)


DimPlot(all, label = TRUE) + NoLegend()

save(all, file = file.path("data", "all.RData"))

#find cluster markers
all.markers<-FindAllMarkers(all,only.pos=T,min.pct=.25,logfc.threshold=.25,return.thresh = 0.05)

save(all.markers, file = file.path("data","diffgenes.RData"))
write_csv(all.markers, file.path("data", "diffgenes.csv"))
