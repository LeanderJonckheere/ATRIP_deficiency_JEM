####MAIN FIGURES####

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

####Figure 3a: visualization of annotated CITE-seq data####

##Umap with Azimuth l2 annotation, no grid lines, no legend, high quality
Azimuth_l2_no_grid <- DimPlot(
  ATRIP,
  reduction = "wnn.umap",
  label = FALSE,
  repel = TRUE,
  shuffle = TRUE,
  group.by = "predicted.celltype.l2",
  alpha = 0.6,
  pt.size = 0.5, 
  cols = annotation_cols,
  label.size = 4
) #+ theme_nothing()  #additional argument to remove legend, axes and grid
Azimuth_l2_no_grid
ggsave(
  "Azimuth_l2_merge_all_samples_no_grid_no_legend.svg",
  plot = Azimuth_l2_no_grid,
  width = 8,
  height = 6
)

##Weighted-nearest neighbour UMAP with Azimuth level 2 automated annotation
##using reference-based mapping, split between condition (patient vs HC).

##See also: Hao Y, Hao S, Andersen-Nissen E, Mauck WM 3rd, Zheng S, Butler A,
##Lee MJ, Wilk AJ, Darby C, Zager M, Hoffman P, Stoeckius M, Papalexi E, Mimitou
##EP, Jain J, Srivastava A, Stuart T, Fleming LM, Yeung B, Rogers AJ, McElrath
##JM, Blish CA, Gottardo R, Smibert P, Satija R. Integrated analysis of
##multimodal single-cell data. Cell. 2021 Jun 24;184(13):3573-3587.e29. doi:
##10.1016/j.cell.2021.04.048. Epub 2021 May 31. PMID: 34062119; PMCID:
##PMC8238499.

umap_HC <- DimPlot(
  ATRIP,
  reduction = "wnn.umap",
  label = FALSE,
  repel = TRUE,
  shuffle = TRUE,
  group.by = "condition",
  cells = WhichCells(ATRIP, expression = condition == "HC"),
  cols = HC_color,
  pt.size = 0.5,
  alpha = 0.6
) + NoLegend() + ggtitle("HC")
umap_HC
ggsave("umap_HC.svg", plot=umap_HC, width=8, height=6)

umap_F1pt <- DimPlot(
  ATRIP,
  reduction = "wnn.umap",
  label = FALSE,
  repel = TRUE,
  shuffle = TRUE,
  group.by = "condition",
  cells = WhichCells(ATRIP, expression = condition == "ATRIP"),
  cols = F1Pt_color,
  pt.size = 0.5,
  alpha = 0.6
) + NoLegend() + ggtitle("ATRIP")
umap_F1pt
ggsave("umap_F1pt.svg", plot=umap_F1pt, width=8, height=6)

##Barplot with celltype ratios per patient sample
cell_types <- c(
  "B intermediate",
  "B memory",
  "B naive",
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
names(annotation_cols) <- cell_types
ATRIP$predicted.celltype.l2 <- factor(ATRIP$predicted.celltype.l2, 
                                      levels = cell_types)

####Prepare the data for the plot
metadata <- ATRIP@meta.data
cell_type_ratios <- metadata %>%
  group_by(clinical.sample, predicted.celltype.l2) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  group_by(clinical.sample) %>%
  mutate(ratio = count / sum(count)) %>%
  ungroup()

###Calculate the total ratios for each cell type for the ATRIP patient
atrip_ratios <- cell_type_ratios %>%
  filter(grepl("ATRIP", clinical.sample)) %>%
  group_by(predicted.celltype.l2) %>%
  summarise(total_ratio = sum(ratio), .groups = 'drop') %>%
  arrange(desc(total_ratio))

###Ensure that all cell types are included in the final ordering
final_order <- union(atrip_ratios$predicted.celltype.l2,
                     cell_type_ratios$predicted.celltype.l2)

###Order the factor levels of 'predicted.celltype.l2' according to most
###prevalent celltype in ATRIP
cell_type_ratios$predicted.celltype.l2 <- factor(
  cell_type_ratios$predicted.celltype.l2, levels = final_order)

###Create the plot

####NoLegend
barplot_celltype_ratio_no_legend <- ggplot(cell_type_ratios,
                                           aes(x = clinical.sample, 
                                               y = ratio, 
                                               fill = predicted.celltype.l2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = NULL, y = "Ratio", fill = "Cell Type") +
  scale_fill_manual(values = annotation_cols) + theme_cowplot() + NoLegend()
barplot_celltype_ratio_no_legend
ggsave(
  "barplot_celltype_ratio_no_legend.svg",
  plot = barplot_celltype_ratio_no_legend,
  width = 4,
  height = 6
)

####Legend
barplot_celltype_ratio_legend <- ggplot(cell_type_ratios,
                                        aes(x = clinical.sample, 
                                            y = ratio, 
                                            fill = predicted.celltype.l2)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = NULL, y = "Ratio", fill = "Cell Type") +
  scale_fill_manual(values = annotation_cols) + theme_cowplot()
barplot_celltype_ratio_legend
ggsave(
  "barplot_celltype_ratio_legend.svg",
  plot = barplot_celltype_ratio_legend,
  width = 6,
  height = 6
)

####Figure 3b: circosplots showing TRBV and TRAV pairing in CD8 TEM####

###See also:  Verstaen K, Lammens I, Roels J, Saeys Y, Lambrecht BN, Vandamme N,
###Vanhee S. DALI (Diversity AnaLysis Interface): a novel tool for the
###integrated analysis of multimodal single cell RNAseq data and immune receptor
###profiling. bioRxiv 2021.12.07.471549; doi:
###https://doi.org/10.1101/2021.12.07.471549

##CircosPlotGenes for family to gene distribution
svg("circos_ATRIP_CD8TEM_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(
  ATRIP_T_reclustered,
  group.by = "predicted.celltype.l2",
  subset = "CD8 TEM",
  seed = 123
)
dev.off()

svg("circos_HC1_CD8TEM_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(
  HC1_T_reclustered,
  group.by = "predicted.celltype.l2",
  subset = "CD8 TEM",
  seed = 123
)
dev.off()

svg("circos_HC2_CD8TEM_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(
  HC2_T_reclustered,
  group.by = "predicted.celltype.l2",
  subset = "CD8 TEM",
  seed = 123
)
dev.off()

svg("circos_HC3_CD8TEM_FTG.svg",
    width = 10,
    height = 8)
CircosPlotGenes(
  HC3_T_reclustered,
  group.by = "predicted.celltype.l2",
  subset = "CD8 TEM",
  seed = 123
)
dev.off()

####Figure 3c: frequency plot of unique cell clones for CD8 TEM per sample####
ClonotypeFrequency(
  CD8_TEM_all_patients,
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
svg("clonotypefreq_CD8_TEM_splitbysample.svg", width = 10, height = 8)

####Figure 3d: distribution of CDR3 region length in CD8 TEM per sample####
svg("CDR3_length_distribution_CD8_TEM_ridge.svg", width = 12, height = 8)
CDR3Plot(
  CD8_TEM_all_patients,
  group.by = "clinical.sample",
  subset = NULL,
  plot.type = c("ridge"), #or ridge
  sequence.type = c("AA"),
  color.theme = ColorThemes(),
  colors = c("#FF8000", "#E9ECEF","#CED4DA","#6C757D")
)
dev.off()

####Figure 4a: heatmap for top10 enriched hallmark gene sets in F1Pt vs HCs####

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

###See also: An algorithm for fast preranked gene set enrichment analysis using
###cumulative statistic calculation Alexey A. Sergushichev, bioRxiv 060012; doi:
###https://doi.org/10.1101/060012

##Order results according to normalized enrichment score (NES)
fgseaRes_ord <- fgseaRes[order(fgseaRes$NES,decreasing=TRUE),]
head(fgseaRes_ord)
length(unique(fgseaRes_ord$padj))

##Heat maps
ATRIP$clinical.sample <- factor(ATRIP$clinical.sample, 
                                      levels = c('ATRIP', 'HC1', 'HC2', 'HC3'))
levels(ATRIP$clinical.sample)
AE <- AggregateExpression(ATRIP,group.by="clinical.sample")
head(AE$RNA)
ncells <- table(ATRIP$clinical.sample)
AE$RNA <- sweep(AE$RNA,2,ncells,"/")
AE$RNA <- AE$RNA*5000

##Create data frame with sample annotation for the heatmap
samples <- ATRIP$clinical.sample
condition <- ATRIP$condition
df <- unique(data.frame(samples,condition))
row.names(df) <- df$samples
df$samples <- NULL
df

##Look at expression patterns by scaling the genes
ann_colors <- c("ATRIP" = "#FFB771", "HC" = "#A0A0A4") 
desired_order <- c('ATRIP', 'HC1', 'HC2', 'HC3')
AEPath <- AEPath[, desired_order]

##Check for NAs
any(is.nan(AEPath))

##Check for rows with zero variance
zero_variance_genes <- apply(AEPath, 1, var) == 0
zero_variance_gene_names <- rownames(AEPath)[zero_variance_genes]
AEPath_filtered <- AEPath[!zero_variance_genes, ]

##Create a heatmap for all gene sets combined (top 10)
TopPathways <- as.character(fgseaRes_ord$pathway[1:10])
L <- lapply(fgsea_sets[TopPathways],FUN=function(x)
  AE$RNA[rownames(AE$RNA) %in% x,])
head(L[[1]])
PathSum <- t(sapply(L,FUN=colMeans))
head(PathSum)
pheatmap(log10(PathSum+1),annotation_col=df,fontsize=6) 
heatmap_msigdb_hallmark_all_PBMCs_split_by_sample <- pheatmap(  
  PathSum,
  annotation_col = df,
  scale = "row",
  fontsize = 8,
  annotation_colors = list(condition = ann_colors),
  main = "Gene set expression heatmap of top 10 enriched MSigDB 
  hallmark gene sets in PBMCs\n\n"
)
ggsave(
  "heatmap_msigdb_hallmark_all_PBMCs_split_by_sample.svg",
  plot = heatmap_msigdb_hallmark_all_PBMCs_split_by_sample,
  width = 8,
  height = 8
)

####Figure 4b: volcano plot showing DEGs in PBMCs of F1Pt compard to HCs####

##Calculate DEGs for all PBMCs, F1pt vs all healthy controls (pooled)
DefaultAssay(ATRIP) <- "RNA"
all_markers_ATRIPvsHC_RNA <- FindMarkers(
  ATRIP,
  ident.1 = "ATRIP",
  ident.2 = "HC",
  group.by = "condition",
  min.pct = 0.1
)

##Generate volcanoplot
volcano_allmarkers_atripvsHC <- EnhancedVolcano(all_markers_ATRIPvsHC_RNA,
                                    lab = rownames(all_markers_ATRIPvsHC_RNA),
                                                x = 'avg_log2FC',
                                                y = 'p_val',
                                                subtitle = NULL,
                                                pCutoff = 10e-32,
                                                FCcutoff = 1,
                                                labSize = 3,
                                                pointSize = 2.5,
                                                drawConnectors = TRUE,
                            col = c("#CED4DA", "#212529","#6C757D", "#FF8000"),
                                                colAlpha = 0.7,
                                                xlim = c(-3,3)
)  +   ggplot2::coord_cartesian(xlim=c(-3, 3)) +
  ggplot2::scale_x_continuous(
    breaks=seq(-3,3, 1))
volcano_allmarkers_atripvsHC
ggsave("volcano_allmarkers_atripvsHC.svg", plot = volcano_allmarkers_atripvsHC,
       width = 8, height = 8)
