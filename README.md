ATAC_Test
Configuration
All parameters are stored in config.yaml.

Update this file according to your input data and analysis preferences.

Create the Conda environment
conda env create -f env.yaml

conda activate signac_pipeline

Dependencies
All dependencies are included in env.yaml

About the Pipeline
This snakemake workflow provides a complete, reproducible framework for processing ATAC-seq datasets and performing downstream analyses.

-> Cell Ranger ATAC : Alignment, peak calling, and generation of fragment/barcode matrices.

-> Data import and quality control : Loading Cell Ranger outputs into R/Signac and generating QC metrics and visualizations.

-> Demultiplexing with cellSNP-lite and Vireo : Assigning cells to donors based on genotype information.

-> Signac normalization : Normalizing ATAC data.

-> Differential accessibility analysis : Identifying differentially accessible peaks across clusters.

-> Motif enrichment analysis : examining transcription factor motif activity in differentially accessible regions.

