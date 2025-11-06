import glob
import os
import re

# --- CONFIG ---
configfile: "config.yaml"
POOLED = config["pooled"]
NUM_DONORS = config.get("num_donors", 3)
CORES = config["cores"]
MEM = config["mem"]
FASTQ_DIR = config["fastq_dir"]
REF = config["ref"]
OUTPUT_DIR = "cellranger_output"
OUTPUT_DIR_QC = "cellranger_output/QC"
DONOR_SNP_FILE = config["donor_snp_file"]

#### SAMPLES ####
fastq_files = glob.glob(f"{FASTQ_DIR}/*_R1_*.fastq.gz")
SAMPLES = sorted(list({re.match(r"(.*)_S\d+_", os.path.basename(f)).group(1) for f in fastq_files}))

####  RULE ALL ####
rule all:
    input:
        # Cell Ranger outputs
        expand(os.path.join(OUTPUT_DIR, "{sample}", "outs", "summary.csv"), sample=SAMPLES),
        # R import/QC outputs
        expand(os.path.join(OUTPUT_DIR_QC, "{sample}_pbmc.rds"), sample=SAMPLES),
        expand(f"{OUTPUT_DIR_QC}/{{sample}}_FragmentHistogram.pdf", sample=SAMPLES),
        expand(f"{OUTPUT_DIR_QC}/{{sample}}_ViolinPlot.pdf", sample=SAMPLES),
        expand(f"{OUTPUT_DIR_QC}/{{sample}}_DensityScatter.pdf", sample=SAMPLES),
        # cellSNP + Vireo
        expand(f"{OUTPUT_DIR_QC}/cellSNP_noGT/{{sample}}/cellSNP.done", sample=SAMPLES),
        expand(f"{OUTPUT_DIR_QC}/vireo/{{sample}}/vireo.done", sample=SAMPLES),
        # Signac normalized
        expand(f"{OUTPUT_DIR_QC}/{{sample}}_Signac_normPlots/{{sample}}_Signac_norm.rds", sample=SAMPLES),
        expand(f"{OUTPUT_DIR_QC}/{{sample}}_SignacNorm_DA/DA_peaks_clusters.csv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR_QC}/{{sample}}_SignacNorm_DA/motif_top.pdf", sample=SAMPLES)



####  Cell Ranger ATAC ####
rule cellranger_atac:
    input:
        fastqs=lambda wc: glob.glob(f"{FASTQ_DIR}/{wc.sample}*_R*_*.fastq.gz")
    output:
        summary=os.path.join(OUTPUT_DIR, "{sample}", "outs", "summary.csv")
    params:
        outid=lambda wc: wc.sample,
        outdir=lambda wc: os.path.join(OUTPUT_DIR, wc.sample),
        ref=REF
    threads: CORES
    resources:
        mem_mb=MEM
    shell:
        """
        mkdir -p {params.outdir}
        cellranger-atac count \
            --id={params.outid} \
            --fastqs={FASTQ_DIR} \
            --sample={wildcards.sample} \
            --reference={params.ref} \
            --localcores={threads} \
            --localmem={resources.mem_mb}
        mv {params.outid}/* {params.outdir}/
        """

####  Import Data ####
rule import_data:
    input:
        h5=lambda wc: f"{OUTPUT_DIR}/{wc.sample}/outs/filtered_peak_bc_matrix.h5",
        metadata=lambda wc: f"{OUTPUT_DIR}/{wc.sample}/outs/singlecell.csv",
        fragments=lambda wc: f"{OUTPUT_DIR}/{wc.sample}/outs/fragments.tsv.gz"
    output:
        rdata=os.path.join(OUTPUT_DIR_QC, "{sample}_pbmc.rds")
    script:
        "scripts/import.R"

# --- QC ---
rule qc:
    input:
        rdata=os.path.join(OUTPUT_DIR_QC, "{sample}_pbmc.rds")
    output:
        FragmentHistogram=os.path.join(OUTPUT_DIR_QC, "{sample}_FragmentHistogram.pdf"),
        ViolinPlot=os.path.join(OUTPUT_DIR_QC, "{sample}_ViolinPlot.pdf"),
        DensityScatter=os.path.join(OUTPUT_DIR_QC, "{sample}_DensityScatter.pdf")
    script:
        "scripts/QC.R"

#### cellSNP-lite ####
rule cellsnp_snps:
    input:
        bam=lambda wc: f"{OUTPUT_DIR}/{wc.sample}/outs/possorted_bam.bam",
        barcodes=lambda wc: f"{OUTPUT_DIR}/{wc.sample}/outs/filtered_peak_bc_matrix/barcodes.tsv",
        snp=DONOR_SNP_FILE
    output:
        done=f"{OUTPUT_DIR_QC}/cellSNP_noGT/{{sample}}/cellSNP.done"
    threads: CORES
    shell:
        """
        mkdir -p {os.path.dirname(output.done)}
        cellsnp-lite -s {input.bam} -b {input.barcodes} -O {os.path.dirname(output.done)} -p {threads} --genotype -R {input.snp}
        touch {output.done}
        """

#### Vireo ####
rule vireo:
    input:
        cellsnp_done=f"{OUTPUT_DIR_QC}/cellSNP_noGT/{{sample}}/cellSNP.done"
    output:
        done=f"{OUTPUT_DIR_QC}/vireo/{{sample}}/vireo.done"
    params:
        N=NUM_DONORS,
        randSeed=2
    threads: CORES
    shell:
        """
        mkdir -p {os.path.dirname(output.done)}
        vireo -c {os.path.dirname(input.cellsnp_done)} -N {params.N} -o {os.path.dirname(output.done)} --randSeed {params.randSeed}
        touch {output.done}
        """

#### Signac normalization ####
rule norm_signac:
    input:
        rdata=os.path.join(OUTPUT_DIR_QC, "{sample}_pbmc.rds"),
        vireo_done=f"{OUTPUT_DIR_QC}/vireo/{{sample}}/vireo.done"
    output:
        rdata_norm=os.path.join(OUTPUT_DIR_QC, "{sample}_Signac_normPlots", "{sample}_Signac_norm.rds")
    script:
        "scripts/normalization.R"


#####  Differential accessibility (peak analysis) #####
rule peak:
    input:
       rdata_norm=os.path.join(OUTPUT_DIR_QC, "{sample}_Signac_normPlots", "{sample}_Signac_norm.rds")
    output:
        da_csv=os.path.join(OUTPUT_DIR_QC, "{sample}_SignacNorm_DA", "DA_peaks_clusters.csv"),
        plots_dir=os.path.join(OUTPUT_DIR_QC, "{sample}_SignacNorm_DA", "plots")
    script:
        "scripts/peaks.R"

#####  Motif  analysis #####

rule motif_analysis:
    input:
        rdata_norm=os.path.join(OUTPUT_DIR_QC, "{sample}_Signac_normPlots", "{sample}_Signac_norm.rds"),
        da_peaks=os.path.join(OUTPUT_DIR_QC, "{sample}_SignacNorm_DA", "DA_peaks_clusters.csv")
    output:
        motif_plot=os.path.join(OUTPUT_DIR_QC, "{sample}_SignacNorm_DA", "motif_top.pdf")
    script:
        "scripts/motif.R"

