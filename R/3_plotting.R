library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)
library(Seurat)
library(ggrepel)
library(ggExtra)

date = Sys.Date()
load(file.path("data", "all.RData"))
source("R/functions.R")

#export sessionInfo
writeLines(capture.output(sessionInfo()), "data/sessionInfo.txt")

if (!file.exists("data/metadata.RData")) {
  #build data frame
  metadata <- data.frame(Cluster = all@active.ident, Region = all$Region, all@reductions$umap@cell.embeddings) 
  metadata$ID <- rownames(metadata)
  
  if (!file.exists("data/retain_cl.RData")){
    retain_cl <- levels(all)
    save(retain_cl, file = "data/retain_cl.RData")
  } else {
    load("data/retain_cl.RData")
  }
  
  if (!file.exists('data/ord_clust.RData')) {
    ord_clust <- retain_cl
    save(ord_clust, file = 'data/ord_clust.RData')
  } else {
    load('data/ord_clust.RData')
  }
  
  #metadata <- metadata[metadata$Cluster %in% retain_cl,]
  
  save(metadata, file="data/metadata.RData")
} else {
  load("data/metadata.RData")
  load("data/retain_cl.RData")
  load('data/ord_clust.RData')
}

#Region umap no legends
umap <- DimPlot(all, group.by = "Region", pt.size = 5) +
  scale_color_brewer(palette = "Dark2") +
  NoLegend() +
  theme_void()

umap

ggsave(paste0('plots/region-umap-plot.pdf'), width = 8.57, height = 5.79)  

#marimekko plot - patients
mosaicGG2(metadata, "Cluster", "Region", c(colors_pat, colors_many), rect_col = 'black', line_width = 0.1) +
  scale_fill_brewer(palette = "Accent")
ggsave(paste0('plots/Region-marimekko-cluster-plot.pdf'))

#heatmap of marker genes
load(file.path("data","diffgenes.RData"))

top20 <- all.markers %>% 
  filter(!grepl("(Htra|Lin|EEF|CTC-|MIR|CTD-|AC0|RP|Fos|Jun|MTRNR|MT-|Xist|Dusp|Zfp36|Rgs|Pmaip1|Hsp|Neat1|Hist|Malat1|Gm|mt-|Rik$)", gene)) %>%
  group_by(cluster) %>% 
  top_n(n = 15, wt = avg_log2FC) %>%
  filter(cluster %in% c("0","1","2"))


#define genes

genes <- str_to_sentence(c("Cx3cr1", "Tmem119", "Hexb", "P2ry12", "SLC2A5", "OLFML3", "CSF1R","Rho", "Mid1", "Scoc","Luc7l3", "Prdx1","FAXC", "CALM2","RBFOX3", "PDGFRA", "MOG", "MAG", "CADPS2", "SYN3", "GAD2", "AQP4", "FGFR3", "MFGE8", "PRDX6", "SOX9", "SLC1A3", "CTPS", "BMP4", "NEU4", "OPALIN","PLP1"))
all2 <- subset(all, cells = colnames(all)[all$seurat_clusters %in% c("0","1","2")])
Idents(all2) <- all2$Region
all.markers2 <- FindAllMarkers(all2,only.pos=T,min.pct=.25,logfc.threshold=.25,return.thresh = 0.05)

heat <- DoHeatmap(all2,features = genes, group.colors = rev(colors_pat)) 
heat + 
  scale_fill_viridis(option = "B")

ggsave(paste0("plots/heatmaps/",date,"-early-micr-top15-gene-heatmap.pdf"), width = 30, height = 20)

#violin plots
violin <- all2[["SCT"]]@data[genes[genes %in% rownames(all2)],] %>%
  as.matrix() %>%
  t %>%
  as.data.frame() %>%
  mutate(Cluster=all2$Region,
         ID=colnames(all2))

violin <- violin %>%
  pivot_longer(Cx3cr1:Plp1 ,"Gene", values_to="Expression")

#line plots
p <- gene_line_plot(violin, "Ctps") +
  scale_color_brewer(palette = "Dark2") +
  theme(text = element_text(size=20)) +
  geom_point(size=0)
ggMarginal(p, data = violin[violin$Gene=="Hexb",], type = "violin")
gene_line_plot(violin, "Rho") +
  scale_color_brewer(palette = "Dark2") +
  theme(text = element_text(size=20))
gene_line_plot(violin, "Mid1") +
  scale_color_brewer(palette = "Dark2") +
  theme(text = element_text(size=20))

walk(genes, function(x) {
  tryCatch({
  p <- gene_line_plot(violin, x) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme(text = element_text(size=10))
  print(p)
  ggsave(file.path("plots", paste0(x,"_gene_line_plot.pdf")), useDingbats=F, width = 20, height = 2)
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

walk(genes, function(x) {
  tryCatch({
  p1 <- ggplot(violin[violin$Gene == x,], aes(Cluster, Expression, fill=Cluster)) +
    geom_violin(scale = "width", draw_quantiles = .5) +
    scale_fill_brewer(palette = "Dark2") +
    theme_void() +
    theme(text = element_text(size=10))
  print(p1)
  ggsave(file.path("plots", paste0(x,"_gene_violin_plot.pdf")), useDingbats=F, width = 10, height = 10)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})


if (!file.exists(file.path("data","micr_eye_vs_brain_diffgenes.RData"))) {
  micr <- FindAllMarkers(all2,
                       logfc.threshold = 0.01,
                       min.pct = 0.01) %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 2.225074e-308, p_val_adj))

  save(micr, file = file.path("data","micr_eye_vs_brain_diffgenes.RData"))
} else {
    load("data/micr_eye_vs_brain_diffgenes.RData")
  }


top5_wt <- micr %>% 
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", gene)) %>%
  top_n(30, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes)

top_5_both <- micr %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", gene)) %>%
  top_n(-30, avg_log2FC) %>%
  mutate(show_genes = gene) %>%
  select(gene, show_genes) %>%
  bind_rows(top5_wt) 

micr <- micr %>%
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", gene)) %>%
  left_join(top_5_both) %>%
  filter(avg_log2FC > 0) %>%
  mutate(
    genes_sig = ifelse(p_val_adj < .05 & avg_log2FC > .2, "sig.", "not sig."),
    #show_genes = ifelse(genes_sig == "sig.", gene, NA),
    avg_log2FC = case_when(
      cluster == "brain"  ~ -1 * avg_log2FC,
      TRUE ~ avg_log2FC
    )
  ) %>%
  distinct()

ggplot(micr, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  geom_point(size=5) + #aes(size=avg_logFC)
  geom_text_repel(size=7, box.padding=1.15, max.overlaps = 150) +
  expand_limits(x=c(-2.25, 2.25)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. logFC", y="-log10 transf. adj. p-value")



