# motif_analysis.R
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(motifmatchr)
library(GenomicRanges)
library(ggseqlogo)
library(patchwork)
library(ggplot2)

#####  Inputs / Outputs via Snakemake #####
object_file    <- snakemake@input$rdata_norm
da_peaks_file  <- snakemake@input$da_peaks
motif_plot_file <- snakemake@output$motif_plot

object <- readRDS(object_file)
da_peaks <- read.csv(da_peaks_file, header = TRUE, stringsAsFactors = FALSE)

#####Cluster UMAP plot #####
p <- DimPlot(object, label = TRUE, pt.size = 0.1) + NoLegend()
ggsave(file.path(dirname(motif_plot_file), "UMAP_clusters.pdf"), plot = p, width = 6, height = 5)

##### Get JASPAR motifs #####
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

##### Add motif information (if not already added) #####
if (!"motifs" %in% names(object@assays$peaks)) {
  object <- AddMotifs(
    object = object,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )
}

##### Select top DA peaks #####
top.da.peak <- da_peaks$peak[da_peaks$p_val_adj < 0.05 & da_peaks$pct.1 > 0.2]

if (length(top.da.peak) == 0) {
  stop("No DA peaks passed threshold for motif analysis.")
}

#####  Find enriched motifs #####
enriched.motifs <- FindMotifs(
  object = object,
  features = top.da.peak
)

##### Save results #####
write.csv(enriched.motifs, file = file.path(dirname(motif_plot_file), "EnrichedMotifs.csv"))

#####  Plot top motifs #####
enriched.motifs <- enriched.motifs[order(enriched.motifs$pvalue), ]
top_motifs <- head(rownames(enriched.motifs), 5)
p_motif <- MotifPlot(object, motifs = top_motifs)
ggsave(filename = motif_plot_file, plot = p_motif, width = 6, height = 4)

cat("Motif analysis completed. Top motifs plotted.\n")
