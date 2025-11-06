#########################################
# Generate multiple QC plots 
#########################################

# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
library(future)
library(BSgenome.Hsapiens.UCSC.hg38)

library(Signac)

data("blacklist_hg38_unified", package = "Signac")

# --- Inputs/Outputs from Snakemake ---

seurat_file <- snakemake@input[['rdata']]
output_fragmenthist <- snakemake@output[["FragmentHistogram"]]
output_violin <- snakemake@output[["ViolinPlot"]]
output_density <- snakemake@output[["DensityScatter"]]

object <- readRDS(seurat_file)

#compute nucleosome signal score per cell
object <- NucleosomeSignal(object = object)

print("compute TSS enrichment score per cell")

object <- TSSEnrichment(object = object)

###### add fraction of reads in peaks ###### 
object$pct_reads_in_peaks <- object$peak_region_fragments / object$passed_filters * 100

# add blacklist ratio
object$blacklist_ratio <- FractionCountsInRegion(
  object = object, 
  assay = 'peaks',
  regions = blacklist_hg38_unified
)

######  Generate QC Plots ###### 


p1 <- DensityScatter(object, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

pdf(output_density, width=7, height=5)

print(p1)


dev.off()

###### Fragments per Cell###### 

object$nucleosome_group <- ifelse(object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p2 <- FragmentHistogram(object = object, group.by = 'nucleosome_group')

# Save Plot 2
pdf(output_fragmenthist, width=7, height=5)
print(p2)
dev.off()


p3 <- VlnPlot(
  object = object,
  features = c('nCount_peaks', 'TSS.enrichment', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

# Save Plot 3
pdf(output_violin , width=10, height=5)
print(p3)
dev.off()


