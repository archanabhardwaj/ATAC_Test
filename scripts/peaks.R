# peak.r
library(Signac)
library(Seurat)
library(dplyr)
library(ggplot2)
library(presto)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


#####  Snakemake inputs/outputs #####
rdata_norm <- snakemake@input$rdata_norm
output_csv <- snakemake@output$da_csv
plots_dir <- snakemake@output$plots_dir
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

sample <- basename(dirname(rdata_norm))

object <- readRDS(rdata_norm)

#####  Differential accessibility (DA)#####  
DefaultAssay(object) <- "peaks"
da_results <- list()
clusters_label <- sort(unique(Idents(object)))

for (clust in clusters_label ) {
    da_peaks <- FindMarkers(object,
                            ident.1 = as.character(clust),
                            test.use = 'wilcox',
                            min.pct = 0.1)
    da_peaks$cluster <- clust
    da_results[[paste0("cluster_", clust)]] <- da_peaks
}

da_peaks_all <- do.call(rbind, da_results)
write.csv(da_peaks_all, file = output_csv, row.names = TRUE)

#####   Top peak plots #####  
if(nrow(da_peaks_all) > 0){
	top_peak <- rownames(da_peaks_all)[which.min(da_peaks_all$p_val_adj)]
    plot1 <- VlnPlot(object, features = top_peak, pt.size = 0.1)
    plot2 <- FeaturePlot(object, features = top_peak, pt.size = 0.1)
    ggsave(file.path(plots_dir, paste0("VlnPlot_top_peak.pdf")), plot1, width = 6, height = 4)
    ggsave(file.path(plots_dir, paste0("FeaturePlot_top_peak.pdf")), plot2, width = 6, height = 4)
}

cat("DA analysis done for", sample, "\n")
