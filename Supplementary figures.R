####SUPPLEMENTARY FIGURES####

####Required packages###
library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggrepel)
library(EnhancedVolcano)
library(presto)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(DALI)

####Load data####

#CITE-seq
ATRIP <- readRDS("ATRIP_18062024.rds")
metadata <- ATRIP@meta.data

#TCR-seq
LBA001_LBA002_T_reclustered <- readRDS("LBA001_LBA002_T_reclustered.rds")
ATRIP_T_reclustered <- readRDS("ATRIP_T_reclustered.rds")
HC1_T_reclustered <- readRDS("HC1_T_reclustered.rds")
HC2_T_reclustered <- readRDS("HC2_T_reclustered.rds")
HC3_T_reclustered <- readRDS("HC3_T_reclustered.rds")

CD8_TEM_all_patients <- readRDS("CD8_TEM_all_patients.rds")
CD8_TEM_ATRIP <- readRDS("CD8_TEM_ATRIP.rds")
CD8_TEM_HC1 <- readRDS("CD8_TEM_HC1.rds")
CD8_TEM_HC2 <- readRDS("CD8_TEM_HC2.rds")
CD8_TEM_HC3 <- readRDS("CD8_TEM_HC3.rds")

#BCR-seq
LBA003_LBA004_BCR_merge_Azimuth_subset <- readRDS(
  "LBA003_LBA004_BCR_merge_Azimuth_subset.rds")
ATRIP_B <- readRDS("ATRIP_B.rds")
HC1_B <- readRDS("HC1_B.rds")
HC2_B <- readRDS("HC2_B.rds")
HC3_B <- readRDS("HC3_B.rds")

####Set colours####
annotation_cols <-
  c("#238B45FF",
    "#00441BFF",
    "#9ECAE1FF",
    "#2171B5FF",
    "#8C6BB1FF",
    "#88419DFF",
    "#810F7CFF",
    "#FC8D59FF",
    "#EF6548FF",
    "#D7301FFF",
    "#7F0000FF",
    "#878787FF",
    "#FA9FB5FF",
    "#F768A1FF",
    "#DD3497FF",
    "#80CDC1FF",
    "#8C510AFF",
    "#BF812DFF",
    "#DFC27DFF",
    "#FD8D3CFF",
    "#FC4E2AFF",
    "#E31A1CFF",
    "#BD0026FF",
    "#477B95FF",
    "#315B88FF")

HC_color <- "#606060"
F1Pt_color <- "#FF8000"

####Supplementary figure 2d: dotplot signature genes used for annotation####

##DE testing an Azimuth cluster: top 3 RNA
Idents(ATRIP) <- 'predicted.celltype.l2'
DefaultAssay(ATRIP) <- "RNA"
Azimuth_l2_RNA <- FindAllMarkers(ATRIP, only.pos = TRUE, min.pct = 0.1)
Azimuth_l2_RNA %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3_Azimuth_l2_RNA
write.csv2(top3_Azimuth_l2_RNA, file = "top3_Azimuth_l2_RNA.csv")

##Dotplot
Idents(ATRIP) <- "predicted.celltype.l2"
DefaultAssay(ATRIP) <- "RNA"
unique_genes <- unique(top3_Azimuth_l2_RNA$gene)
cell_type_order <- c(
  "B naive",
  "B intermediate",
  "B memory",
  "CD14 Mono",
  "CD16 Mono",
  "CD4 CTL",
  "CD4 Naive",
  "CD4 TCM",
  "CD4 TEM",
  "CD8 Naive",
  "CD8 TCM",
  "CD8 TEM",
  "cDC1",
  "cDC2",
  "dnT",
  "gdT",
  "ILC",
  "MAIT",
  "NK",
  "NK Proliferating",
  "NK_CD56bright",
  "pDC",
  "Plasmablast",
  "Platelet",
  "Treg"
)

dotplot_rna <- DotPlot(
  ATRIP,
  features = unique_genes,
  group.by = "predicted.celltype.l2",
  dot.scale = 4,
  cols = "RdYlBu"
) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 9,
      family = "Arial"
    ),
    axis.text.y = element_text(family = "Arial"),
    axis.title = element_text(family = "Arial"),
    plot.title = element_text(family = "Arial")
  ) +
  scale_y_discrete(limits = rev(cell_type_order))
dotplot_rna
ggsave("dotplot_rna.svg",
       plot = dotplot_rna,
       width = 14,
       height = 8)

###OPTIONAL: Create the DotPlot with improved aesthetics and light gridlines
dotplot_rna_gridlines <- DotPlot(
  ATRIP,
  features = unique_genes,
  group.by = "predicted.celltype.l2",
  dot.scale = 4,
  cols = "RdYlBu"
) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10,
      family = "Arial"
    ),
    axis.text.y = element_text(size = 10, family = "Arial"),
    axis.title.x = element_text(size = 12, family = "Arial"),
    axis.title.y = element_text(size = 12, family = "Arial"),
    plot.title = element_text(
      size = 14,
      family = "Arial",
      face = "bold"
    ),
    legend.title = element_text(size = 12, family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    panel.background = element_blank()
  ) +
  labs(x = "Predicted Cell Types",
       y = "Genes",
       title = "Dot Plot of Gene Expression Across Predicted Cell Types",
       fill = "Expression Level") +
  scale_y_discrete(limits = rev(cell_type_order))
ggsave(
  "dotplot_rna_gridlines.svg",
  plot = dotplot_rna_gridlines,
  width = 14,
  height = 8
)

####Supplementary figure 2e: circosplots showing TRBV & TRAV pairing####

###See also:  Verstaen K, Lammens I, Roels J, Saeys Y, Lambrecht BN, Vandamme N,
###Vanhee S. DALI (Diversity AnaLysis Interface): a novel tool for the
###integrated analysis of multimodal single cell RNAseq data and immune receptor
###profiling. bioRxiv 2021.12.07.471549; doi:
###https://doi.org/10.1101/2021.12.07.471549

svg("circos_ATRIP_allT_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(ATRIP_T_reclustered, group.by = "predicted.celltype.l2", 
                seed = 123)
dev.off()

svg("circos_HC1_allT_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(HC1_T_reclustered, group.by = "predicted.celltype.l2", 
                seed = 123)
dev.off()

svg("circos_HC2_allT_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(HC2_T_reclustered, group.by = "predicted.celltype.l2", 
                seed = 123)
dev.off()

svg("circos_HC3_allT_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(HC3_T_reclustered, group.by = "predicted.celltype.l2", 
                seed = 123)
dev.off()

####Supplementary figure 2f: frequency of unique T cell clones####
clonotypefreq_all_T_splitbysample <- ClonotypeFrequency(
  LBA001_LBA002_T_reclustered,
  chain = c("VDJ", "VJ"),
  group.by = "clinical.sample",
  subset = NULL,
  use.sequence = F,
  sequence.type = c("AA", "NT"),
  clonotype.column = NULL,
  bulk = F,
  show.missing = F,
  plot.type = "bar"
)
svg("clonotypefreq_all_T_splitbysample.svg", width = 10, height = 8)

####Supplementary figure 2g: distribution of CDR3 region lengths for T cells####
svg("CDR3_length_distribution_all_T_ridge.svg", width = 12, height = 8)
CDR3Plot(
  LBA001_LBA002_T_reclustered,
  group.by = "clinical.sample",
  subset = NULL,
  plot.type = "ridge",
  sequence.type = c("AA"),
  color.theme = ColorThemes(),
  colors = c("#FF8000", "#E9ECEF","#CED4DA","#6C757D")
)
dev.off()

####Supplementary figure 2h: circosplots showing IGH, IGK and IGL pairing####
svg("circos_ATRIP_allB_FTG.svg", width = 10, height = 8)
CircosPlotGenes(ATRIP_B, group.by = "predicted.celltype.l2",seed = 123) 
dev.off()

svg("circos_HC1_allB_FTG.svg", width = 10, height = 8)
CircosPlotGenes(HC1_B, group.by = "predicted.celltype.l2",seed = 123) 
dev.off()

svg("circos_HC2_allB_FTG.svg", width = 10, height = 8)
CircosPlotGenes(HC2_B, group.by = "predicted.celltype.l2",seed = 123) 
dev.off()

svg("circos_HC3_allB_FTG.svg", width = 10, height = 8)
CircosPlotGenes(HC3_B, group.by = "predicted.celltype.l2",seed = 123) 
dev.off()

####Supplementary figure 2i: frequency of unique B cell clones####
svg("clonotypefreq_all_B_splitbysample.svg", width = 10, height = 8)
barplot_sample <- ClonotypeFrequency(
  LBA003_LBA004_BCR_merge_Azimuth_subset,
  chain = c("VDJ", "VJ"),
  group.by = "clinical.sample",
  subset = NULL,
  use.sequence = F,
  sequence.type = c("AA", "NT"),
  clonotype.column = NULL,
  bulk = F,
  show.missing = F,
  plot.type = "bar"
)
dev.off()

####Supplementary figure 2j: distribution of CDR3 region length for B cells####
svg("CDR3_length_distribution_allB_ridge.svg", width = 12, height = 8)
CDR3Plot(
  LBA003_LBA004_BCR_merge_Azimuth_subset,
  group.by = "clinical.sample",
  subset = NULL,
  plot.type = c("ridge"), #or line
  sequence.type = c("AA"),
  color.theme = ColorThemes(),
  colors = c("#FF8000", "#E9ECEF","#CED4DA","#6C757D")
)
dev.off()

####Supplementary figure 3a: rank plot####

##Functional analysis: GSEA
DEStats <- wilcoxauc(ATRIP,group_by='condition')
head(DEStats,2)
DEStats[DEStats$feature=="STAT1",]
DES <- DEStats[DEStats$group=="ATRIP",] 
head(DES)

##Hallmark (H) gene set: get functional annotations
msigdbr_species()
print(msigdbr_collections(),n=23)
mdf <- msigdbr(species = "Homo sapiens", category = "H") 
head(mdf)
fgsea_sets <- split(mdf$gene_symbol,f=mdf$gs_name) 

##Prepare GSEA
DES_ord <- DES[order(DES$auc,decreasing=TRUE),] 
head(DES_ord)
ranks <- DES_ord[,"auc"] 
names(ranks) <- DES_ord[,"feature"]
head(ranks)
tail(ranks)

##Perform GSEA
fgseaRes <- fgsea(fgsea_sets,stats=ranks, eps=0, scoreType="pos")

##Order results according to normalized enrichment score (NES)
fgseaRes_ord <- fgseaRes[order(fgseaRes$NES,decreasing=TRUE),]
head(fgseaRes_ord)
length(unique(fgseaRes_ord$padj)) 

##Plot the results
fgseaRes_ord$pathway <- factor(fgseaRes_ord$pathway, 
                               levels = rev(fgseaRes_ord$pathway))
GSEA_top10_hallmark_ATRIPvsHC_all <- ggplot(fgseaRes_ord[1:10, ], 
                                            aes(NES, pathway)) +
  stat_summary(
    geom = "point",
    fun = identity,
    color = "#FF8000",
    size = 4
  ) + xlab("Normalized Enrichment Score (NES)") + theme_bw()
GSEA_top10_hallmark_ATRIPvsHC_all
ggsave(
  "GSEA_top10_hallmark_ATRIPvsHC_all.svg",
  plot = GSEA_top10_hallmark_ATRIPvsHC_all,
  width = 8,
  height = 6
)

####Supplementary figure 3b: enrichment plots####
enrichment_plot_DNA_repair <- plotEnrichment(
  fgsea_sets$HALLMARK_DNA_REPAIR, ranks) + 
  labs(title = "HALLMARK_DNA_REPAIR")
enrichment_plot_DNA_repair
ggsave(
  "enrichment_plot_DNA_repair.svg",
  plot = enrichment_plot_DNA_repair,
  width = 8,
  height = 8
)
enrichment_plot_IFNG_response <- plotEnrichment(
  fgsea_sets$HALLMARK_INTERFERON_GAMMA_RESPONSE, ranks) + 
  labs(title = "HALLMARK_INFγ_RESPONSE")
enrichment_plot_IFNG_response
ggsave(
  "enrichment_plot_IFNG_response.svg",
  plot = enrichment_plot_IFNG_response,
  width = 8,
  height = 8
)

####Supplementary figure 3c: heat map T effector cells####

##Subset T_effector
table(ATRIP$annotation_DE_analysis_effector_T_pooled)
Idents(ATRIP) <- 'annotation_DE_analysis_effector_T_pooled'
T_effector_subset <- subset(ATRIP, idents = "Effector T")
plot <- DimPlot(T_effector_subset, reduction = "wnn.umap")
plot

##Perform wilcoxauc DE test
DEStats <- wilcoxauc(T_effector_subset,group_by='condition', 
                     seurat_assay = "RNA")
head(DEStats,2)
DEStats[DEStats$feature=="STAT1",]
DES <- DEStats[DEStats$group=="ATRIP",] 
head(DES)

##MSigDB hallmark gene set: get functional annotations
mdf <- msigdbr(species = "Homo sapiens", category = "H") 
fgsea_sets <- split(mdf$gene_symbol,f=mdf$gs_name) 

##Prepare GSEA
DES_ord <- DES[order(DES$auc,decreasing=TRUE),] 
ranks <- DES_ord[,"auc"] 
names(ranks) <- DES_ord[,"feature"]

##Perform GSEA
fgseaRes <- fgsea(fgsea_sets,stats=ranks, scoreType="pos")

##Order results according to normalized enrichment score (NES)
fgseaRes_ord <- fgseaRes[order(fgseaRes$NES,decreasing=TRUE),]
head(fgseaRes_ord)
length(unique(fgseaRes_ord$padj)) 

##Heat maps
AE <- AggregateExpression(T_effector_subset,group.by="clinical.sample")
head(AE$RNA)
ncells <- table(T_effector_subset$clinical.sample) 
AE$RNA <- sweep(AE$RNA,2,ncells,"/")
AE$RNA <- AE$RNA*5000

##Create data frame with sample annotation for the heatmap
samples <- T_effector_subset$clinical.sample
condition <- T_effector_subset$condition
df <- unique(data.frame(samples,condition))
row.names(df) <- df$samples
df$samples <- NULL
df

##Create heatmap 
ann_colors <- c("ATRIP" = "#FFB771", "HC" = "#A0A0A4")

TopPathways <- as.character(fgseaRes_ord$pathway[1:20])
L <- lapply(fgsea_sets[TopPathways],FUN=function(x)
  AE$RNA[rownames(AE$RNA) %in% x,])
PathSum <- t(sapply(L,FUN=colMeans))

heatmap_hallmark_T_effector <- pheatmap(
  PathSum,
  annotation_col = df,
  scale = "row",
  fontsize = 8,
  annotation_colors = list(condition = ann_colors),
  border_color = "black",
  cluster_cols = FALSE,
  main = "T effector cells"
)
ggsave("heatmap_hallmark_T_effector.svg", 
       plot = heatmap_hallmark_T_effector, width = 8, height = 8)
heatmap_hallmark_T_effector
dev.off()

####Supplementary figure 3d: heatmap NK cells####

##GSEA: subset NK
Idents(ATRIP) <- 'annotation_DE_analysis_effector_T_pooled_NK_pooled'
NK_subset <- subset(ATRIP, idents = "NK")
plot <- DimPlot(NK_subset, reduction = "wnn.umap")
plot

##Perform wilcoxauc
DEStats <- wilcoxauc(NK_subset,group_by='condition', seurat_assay = "RNA")
head(DEStats,2)
DEStats[DEStats$feature=="STAT1",]
DES <- DEStats[DEStats$group=="ATRIP",]
head(DES)

##MSigDB hallmark gene set: get functional annotations
mdf <- msigdbr(species = "Homo sapiens", category = "H") 
fgsea_sets <- split(mdf$gene_symbol,f=mdf$gs_name) 

##Prepare GSEA
DES_ord <- DES[order(DES$auc,decreasing=TRUE),] 
ranks <- DES_ord[,"auc"] 
names(ranks) <- DES_ord[,"feature"]

##Perform GSEA
fgseaRes <- fgsea(fgsea_sets,stats=ranks, scoreType="pos")

###See also: An algorithm for fast preranked gene set enrichment analysis using
###cumulative statistic calculation Alexey A. Sergushichev, bioRxiv 060012; doi:
###https://doi.org/10.1101/060012

##Order results according to normalized enrichment score (NES)
fgseaRes_ord <- fgseaRes[order(fgseaRes$NES,decreasing=TRUE),]
head(fgseaRes_ord)
length(unique(fgseaRes_ord$padj)) 

##Heat maps
AE <- AggregateExpression(NK_subset,group.by="clinical.sample")
head(AE$RNA)
ncells <- table(NK_subset$clinical.sample)
AE$RNA <- sweep(AE$RNA,2,ncells,"/")
AE$RNA <- AE$RNA*5000

##Create data frame with sample annotation for the heatmap
samples <- NK_subset$clinical.sample
condition <- NK_subset$condition
df <- unique(data.frame(samples,condition))
row.names(df) <- df$samples
df$samples <- NULL
df

##Create heatmap 
ann_colors <- c("ATRIP" = "#FFB771", "HC" = "#A0A0A4")

TopPathways <- as.character(fgseaRes_ord$pathway[1:20])
L <- lapply(fgsea_sets[TopPathways],FUN=function(x)
  AE$RNA[rownames(AE$RNA) %in% x,])
PathSum <- t(sapply(L,FUN=colMeans))

heatmap_hallmark_NK <- pheatmap(
  PathSum,
  annotation_col = df,
  scale = "row",
  fontsize = 8,
  annotation_colors = list(condition = ann_colors),
  border_color = "black",
  cluster_cols = FALSE,
  main = "NK cells"
)
ggsave("heatmap_hallmark_NK.svg", plot = heatmap_hallmark_NK, 
       width = 8, height = 8)
heatmap_hallmark_NK
dev.off()
