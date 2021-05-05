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

#normalize dataset 
all <- all %>%
  SCTransform(vars.to.regress = c("percent.mt", "Region"),
              variable.features.n = 10000) %>%
  RunPCA() 

#
ElbowPlot(all)

all<- all %>% 
  RunUMAP(dims=1:20) %>%
  FindNeighbors(dims=1:20) %>%
  FindClusters(resolution=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(all)

ggsave("plots/overview-cluster-resolutions.pdf", useDingbats=F)

all<-FindClusters(all,resolution=.3)


DimPlot(all, label = TRUE) + NoLegend()
DimPlot(all, label = TRUE, group.by = "Region") + NoLegend()

save(all, file = file.path("data", "all.RData"))

#find cluster markers
all.markers<-FindAllMarkers(all,only.pos=T,min.pct=.25,logfc.threshold=.25,return.thresh = 0.05)

save(all.markers, file = file.path("data","diffgenes.RData"))
write_csv(all.markers, file.path("data", "diffgenes.csv"))
