#!/bin/bash

# ==============================================================================
# Script Name:   vg_unified_v1.1.sh
# Version:       1.1 (Pure VG Workflow)
# Logic:         v0.8 (Augmentation/Creation) + v1.0 (Standardized VCF Calling)
# Description:   Full pangenome pipeline using 'vg call' for all variant types.
# ==============================================================================

set -e
set -o pipefail

# --- 1. GRANULAR TOGGLES ---
DO_ALIGN=true      # vg giraffe
DO_AUGMENT=true    # v0.8: Creation of novel SV/CNVs
DO_PACK=true       # Coverage calculation
DO_CALL=true       # v1.0: VCF generation (SNP/SV/CNV)
DO_BAM=true        # Surjection for visualization (IGV)

# --- 2. ENVIRONMENT & DIRECTORIES ---
get_environment() {
    echo "--- Personal Genomics Pipeline v1.1 (Pure VG) ---"
    read -p "Enter Sample Name: " SAMPLE
    
    # Pathing logic from v0.8/v1.0 feedback
    BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    INPUT_DIR="${BASE_DIR}/input"
    OUTPUT_DIR="${BASE_DIR}/output"
    
    # Centralized Reference Volume
    REF_DIR="/data/pangenomes/human"
    
    echo "Select Reference: 1) GRCh38 2) CHM13"
    read -p "Selection: " REF_CHOICE
    case $REF_CHOICE in
        1) REF="GRCh38"; GBZ="${REF_DIR}/GRCh38.gbz" ;;
        2) REF="CHM13";  GBZ="${REF_DIR}/CHM13.gbz" ;;
        *) echo "Invalid selection."; exit 1 ;;
    esac

    THREADS=$(nproc)
    mkdir -p "${OUTPUT_DIR}/vcf" "${OUTPUT_DIR}/bam" "${OUTPUT_DIR}/logs"
}

# --- 3. PROCESSING MODULES ---

run_alignment() {
    local gam="${OUTPUT_DIR}/${SAMPLE}.gam"
    if [[ "$DO_ALIGN" == "true" && ! -f "$gam" ]]; then
        echo "[Step 1] Aligning reads with vg giraffe..."
        vg giraffe -Z "$GBZ" \
            -d "${GBZ%.gbz}.dist" -m "${GBZ%.gbz}.min" \
            -f "${INPUT_DIR}/${SAMPLE}_R1.fastq.gz" \
            -f "${INPUT_DIR}/${SAMPLE}_R2.fastq.gz" \
            -t "$THREADS" > "$gam"
    fi
}

run_v08_creation() {
    local aug_graph="${OUTPUT_DIR}/${SAMPLE}_augmented.pg"
    local aug_gam="${OUTPUT_DIR}/${SAMPLE}_augmented.gam"
    if [[ "$DO_AUGMENT" == "true" && ! -f "$aug_graph" ]]; then
        echo "[Step 2/v0.8] Augmenting graph (Creating novel SV/CNVs)..."
        # This modifies the graph structure based on sample-specific reads
        vg augment "$GBZ" "${OUTPUT_DIR}/${SAMPLE}.gam" \
            -t "$THREADS" -A "$aug_gam" > "$aug_graph"
    fi
}

run_v10_calling() {
    local aug_graph="${OUTPUT_DIR}/${SAMPLE}_augmented.pg"
    local aug_gam="${OUTPUT_DIR}/${SAMPLE}_augmented.gam"
    local pack="${OUTPUT_DIR}/${SAMPLE}.pack"
    local final_vcf="${OUTPUT_DIR}/vcf/${SAMPLE}_unified.vcf.gz"

    if [[ "$DO_CALL" == "true" && ! -f "$final_vcf" ]]; then
        echo "[Step 3] Calculating coverage (vg pack)..."
        vg pack -x "$aug_graph" -g "$aug_gam" -o "$pack" -t "$THREADS"

        echo "[Step 4/v1.0] Calling SVs/CNVs/SNPs to VCF..."
        # Using -a to ensure augmented (created) variants are included
        vg call "$aug_graph" -k "$pack" -r "$REF" -t "$THREADS" -a > "tmp.vcf"
        
        echo "[Step 5] Standardizing VCF format..."
        bcftools sort "tmp.vcf" -O z -o "$final_vcf"
        bcftools index -t "$final_vcf"
        rm "tmp.vcf" "$pack"
    fi
}

run_surjection() {
    local final_bam="${OUTPUT_DIR}/bam/${SAMPLE}_${REF}_ucsc.bam"
    if [[ "$DO_BAM" == "true" && ! -f "$final_bam" ]]; then
        echo "[Step 6] Surjecting to BAM for visualization..."
        vg surject -x "$GBZ" "${OUTPUT_DIR}/${SAMPLE}.gam" -b -p "$REF" -t "$THREADS" > "tmp_raw.bam"
        
        # Preserving v0.8 UCSC naming convention (chr1, chr2...)
        samtools view -H "tmp_raw.bam" | sed -e 's/SN:\([0-9XYM]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' > header.sam
        samtools reheader header.sam "tmp_raw.bam" > "$final_bam"
        samtools index "$final_bam"
        rm "tmp_raw.bam" header.sam
    fi
}

# --- 4. EXECUTION ---
main() {
    get_environment
    run_alignment
    run_v08_creation
    run_v10_calling
    run_surjection
    echo "Pipeline Finished Successfully."
    echo "VCF Location: ${OUTPUT_DIR}/vcf/${SAMPLE}_unified.vcf.gz"
}

main "$@"
