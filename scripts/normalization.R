# norm_signac.R
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)

# --- Inputs/Outputs from Snakemake ---

rdata_file <- snakemake@input$rdata
vireo_done <- snakemake@input$vireo_done
output_rds <- snakemake@output$rdata_norm

sample <- basename(dirname(output_rds))

output_dir <- dirname(output_rds)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

##### LOAD SEURAT OBJECT ##### 
object <- readRDS(rdata_file)


##### READ VIREO OUTPUT ##### 
vireo_file <- file.path(dirname(vireo_done), "donor_assignment.tsv")
vireo_assign <- read.table(vireo_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Ensure barcodes match (Cell Ranger ATAC appends -1)
vireo_assign$barcode <- paste0(vireo_assign$barcode, "-1")

##### ADD DONOR METADATA ##### 
object$donor <- vireo_assign$donor[match(colnames(object), vireo_assign$barcode)]
object$doublet <- vireo_assign$doublet_prob[match(colnames(object), vireo_assign$barcode)]

##### TF-IDF NORMALIZATION ##### 
object <- RunTFIDF(object)
object <- FindTopFeatures(object, min.cutoff = 'q0')
object <- RunSVD(object)

##### CLUSTERING AND UMAP ##### 
object <- FindNeighbors(object, reduction = 'lsi', dims = 2:30)
object <- FindClusters(object, resolution = 0.5)

set.seed(1234)
object <- RunUMAP(object, reduction = 'lsi', dims = 2:30)

# --- SAVE PLOTS ---
#####  UMAP colored by donor ##### 
p1 <- DimPlot(object, group.by = 'donor', label = TRUE) + ggtitle("UMAP by Donor")
ggsave(file.path(output_dir, paste0(sample, "_UMAP_by_donor.pdf")), plot = p1, width = 6, height = 5)

##### UMAP colored by cluster ##### 
p2 <- DimPlot(object, label = TRUE) + ggtitle("UMAP by Cluster")
ggsave(file.path(output_dir, paste0(sample, "_UMAP_by_cluster.pdf")), plot = p2, width = 6, height = 5)

#####  Violin plot of fragments per cluster #####  
p3 <- VlnPlot(object, features = "nCount_peaks", group.by = "seurat_clusters") +
      ggtitle("Fragments per Cluster")
ggsave(file.path(output_dir, paste0(sample, "_Fragments_per_Cluster.pdf")), plot = p3, width = 6, height = 5)

##### SAVE NORMALIZED OBJECT ##### 
saveRDS(object, output_rds)
