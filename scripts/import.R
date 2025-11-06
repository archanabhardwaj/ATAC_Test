library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(future)
library(BSgenome.Hsapiens.UCSC.hg38)


out_dir <- dirname(snakemake@output[['rdata']])
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Inputs/Outputs from Snakemake ---
counts_file <- snakemake@input[['h5']]
metadata_file <- snakemake@input[['metadata']]
fragments_file <- snakemake@input[['fragments']]
output_rds <- snakemake@output[['rdata']]
sample_name <- snakemake@wildcards[['sample']]

# --- Load data ---
counts <- Read10X_h5(filename = counts_file)
metadata <- read.csv(
  file = metadata_file,
  header = TRUE,
  row.names = 1
)

head(metadata)

# --- Create ChromatinAssay ---
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = fragments_file,
  #min.cells = 10,
  #min.features = 200
)

chrom_assay

# --- Create Seurat Object ---
object <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata,
  project = sample_name
)

object


# --- Filter peaks to standard chromosomes ---
peaks.keep <- seqnames(granges(object)) %in% standardChromosomes(granges(object))
object <- object[as.vector(peaks.keep), ]
print("import finished")


## Annotation phase ### 
# extract gene annotations from EnsDb

print("starting annotation download")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
print("finished annotation download")


# change to UCSC style since the data was mapped to hg38
#seqlevelsStyle(annotations) <- 'UCSC'

seqlevels(annotations) <- paste0("chr", seqlevels(annotations))
print("seq annotations done ")
genome(annotations) <- "hg38"
print("genome annotation done ")


# add the gene information to the object
Annotation(object) <- annotations

object 

saveRDS(object, file = output_rds)


