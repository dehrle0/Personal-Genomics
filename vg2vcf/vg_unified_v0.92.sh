#!/usr/bin/env bash
#
# vg_unified_v0.92.sh
# Status: TESTING / MEMORY-OPTIMIZED
#
# Logic: Bypasses HashGraph (.hg) and Augmentation to stay under 40GB RAM.
# Uses .d9.gbz (Downsampled) + .snarls for variant calling.
# Includes full BAM header repair (Stripping prefixes + matching LN to .fai).
#
# v0.92 Fixes:
# - Changed 'vg index -m' to 'vg minimizer' for .min generation.
# - Changed 'vg index -d' to 'vg index -j' for .dist generation.

set -euo pipefail

# --- 1. GRANULAR TOGGLES ---
DO_DOWNLOAD=false
DO_INDEX=true       # Rebuilds .dist, .min, and .snarls
DO_QC=true
DO_ALIGN=true       # vg giraffe
DO_SURJECT=true     # vg surject + Header/LN Fix
DO_VCF_DV=true      # DeepVariant
DO_VCF_CNV=true     # CNVkit
DO_VCF_SV=true      # vg call (Direct from GBZ)
DO_MULTIQC=true

# --- 2. CONFIGURATION ---
SAMPLE="TEST"
SAMPLE_SEX="female" # Options: male, female
THREADS="12"        # IO and Mapping threads
MEM_THREADS="4"     # Throttled threads for RAM-heavy indexing/calling
REF_MODE="GRCh38"   # GRCh38, CHM13, BOTH

# Host Paths
REF_DIR_HOST="/mnt/data/work/references"
HG38_LOCAL="${REF_DIR_HOST}/hg38.fa"
HS1_LOCAL="${REF_DIR_HOST}/chm13v2.0_maskedY_rCRS.fa"

GRAPH_DIR_HOST="${REF_DIR_HOST}/hrpc"
# Reverted to .d9.gbz (Downsampled graph, lower memory)
GBZ_HOST="${GRAPH_DIR_HOST}/hprc-v1.1-mc-chm13.d9.gbz"

READS_DIR_HOST="/mnt/c/Genomes/${SAMPLE}/Data/Source"
FQ1_HOST="${READS_DIR_HOST}/ME_chr20_R1.fastq.gz"
FQ2_HOST="${READS_DIR_HOST}/ME_chr20_R2.fastq.gz"

# Working Directories
WORK_DIR="/mnt/data/work/vg_pangenome"
LOG_DIR="${WORK_DIR}/logs"
QC_DIR="${WORK_DIR}/qc"
CNV_DIR="${WORK_DIR}/cnv"
VAR_DIR="${WORK_DIR}/variants"
COMPLETION_FILE="${WORK_DIR}/completed_steps_${SAMPLE}.txt"

# Docker Images
VG_IMAGE="quay.io/vgteam/vg:v1.70.0"
DV_IMAGE="google/deepvariant:1.10.0-beta"
CNVKIT_IMAGE="etal/cnvkit:latest"
BCFTOOLS_IMAGE="quay.io/biocontainers/bcftools:1.19--h8b25389_0"
SAMTOOLS_IMAGE="staphb/samtools:1.19"
FASTQC_IMAGE="staphb/fastqc:0.12.1"
MULTIQC_IMAGE="ewels/multiqc:v1.19"

# --- 3. UTILITIES ---
mkdir -p "${WORK_DIR}" "${LOG_DIR}" "${QC_DIR}/fastqc" "${CNV_DIR}" "${VAR_DIR}"
touch "${COMPLETION_FILE}"

if [[ "$REF_MODE" == "BOTH" ]]; then REFS=("GRCh38" "CHM13"); elif [[ "$REF_MODE" == "CHM13" ]]; then REFS=("CHM13"); else REFS=("GRCh38"); fi

log() { printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$1" | tee -a "${LOG_DIR}/pipeline_${SAMPLE}.log"; }
step_done() { grep -qx "$1" "${COMPLETION_FILE}" 2>/dev/null || return 1; }
mark_done() { echo "$1" >> "${COMPLETION_FILE}"; }

# --- 4. PRE-FLIGHT & INDEXING ---

download_refs() {
    if [[ "$DO_DOWNLOAD" == "false" ]]; then return; fi
    log "Step 00: Checking references..."
    mkdir -p "${GRAPH_DIR_HOST}" "${REF_DIR_HOST}"

    [[ ! -f "${GBZ_HOST}" ]] && log "Downloading GBZ..." && wget -c -P "${GRAPH_DIR_HOST}" https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/release/v1.1/mc_chm13/hprc-v1.1-mc-chm13.d9.gbz
    
    if [[ ! -f "${HG38_LOCAL}" ]]; then
        log "Downloading hg38..."
        wget -c -O "${HG38_LOCAL}.gz" http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
        gunzip -f "${HG38_LOCAL}.gz"
    fi
    
    [[ ! -f "${HG38_LOCAL}.fai" ]] && docker run --rm -v "${REF_DIR_HOST}":/refs "${SAMTOOLS_IMAGE}" samtools faidx "/refs/$(basename "${HG38_LOCAL}")"
}

rebuild_indexes() {
    if [[ "$DO_INDEX" == "false" ]]; then return; fi
    local BASE_NAME=$(basename "${GBZ_HOST}")
    local SNARLS_HOST="${GBZ_HOST%.gbz}.snarls"

    # 1. Dist Index (Using -j for Snarl Distance Index)
    if [[ ! -f "${GBZ_HOST%.gbz}.dist" ]]; then
        log "Step 01a: Rebuilding .dist index..."
        docker run --rm -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg index -t "${THREADS}" -j "/graph/${BASE_NAME%.gbz}.dist" "/graph/${BASE_NAME}"
    fi

    # 2. Minimizer Index (Using vg minimizer, not vg index)
    if [[ ! -f "${GBZ_HOST%.gbz}.min" ]]; then
        log "Step 01b: Rebuilding .min index..."
        docker run --rm -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg minimizer -t "${THREADS}" \
            -d "/graph/${BASE_NAME%.gbz}.dist" \
            -o "/graph/${BASE_NAME%.gbz}.min" "/graph/${BASE_NAME}"
    fi

    # 3. Snarls (Critical for calling)
    if [[ ! -f "${SNARLS_HOST}" ]]; then
        log "Step 01c: Generating Snarls (Using throttled threads for RAM safety)..."
        docker run --rm -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg snarls -t "${MEM_THREADS}" "/graph/${BASE_NAME}" > "${SNARLS_HOST}"
    fi
}

# --- 5. CORE PIPELINE ---

run_fastqc() {
    if [[ "$DO_QC" == "false" ]] || step_done "05_FastQC"; then return; fi
    log "Step 05: FastQC..."
    docker run --rm -v "${READS_DIR_HOST}":/reads:ro -v "${QC_DIR}/fastqc":/output "${FASTQC_IMAGE}" \
        fastqc -t "${THREADS}" -o /output "/reads/$(basename "${FQ1_HOST}")" "/reads/$(basename "${FQ2_HOST}")"
    mark_done "05_FastQC"
}

run_alignment() {
    if [[ "$DO_ALIGN" == "false" ]] || step_done "10_Align"; then return; fi
    log "Step 10: vg giraffe alignment (GAM output)..."
    
    # Removed -H (haplotype) flag as requested for d9 usage
    docker run --rm -v "${GRAPH_DIR_HOST}":/graph -v "${READS_DIR_HOST}":/reads:ro -v "${WORK_DIR}":/work "${VG_IMAGE}" \
        vg giraffe \
        -Z "/graph/$(basename "${GBZ_HOST}")" \
        -m "/graph/$(basename "${GBZ_HOST%.gbz}.min")" \
        -d "/graph/$(basename "${GBZ_HOST%.gbz}.dist")" \
        -f "/reads/$(basename "${FQ1_HOST}")" -f "/reads/$(basename "${FQ2_HOST}")" \
        -t "${THREADS}" -o GAM > "${WORK_DIR}/${SAMPLE}.gam"
    mark_done "10_Align"
}

run_surjection_and_dv() {
    local ref=$1
    local d_ref_path=$2
    if [[ "$DO_SURJECT" == "false" ]]; then return; fi
    
    if ! step_done "30_Surject_${ref}"; then
        log "Step 30: Surjecting and Repairing BAM for ${ref}..."
        
        local LOCAL_FAI_NAME="$(basename "${HG38_LOCAL}.fai")"
        [[ "$ref" == "CHM13" ]] && LOCAL_FAI_NAME="$(basename "${HS1_LOCAL}.fai")"

        # Surject
        docker run --rm -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg surject -x "/graph/$(basename "${GBZ_HOST}")" -b -t "${THREADS}" --into-ref "${ref}" \
            "/work/${SAMPLE}.gam" > "${WORK_DIR}/tmp.${ref}.bam"

        # Complex Header repair (Strip Prefix + Fix LN)
        docker run --rm -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/refs:ro "${SAMTOOLS_IMAGE}" /bin/bash -c "
            samtools view -H /work/tmp.${ref}.bam > /work/old.sam
            awk -v ref=\"${ref}\" '
                BEGIN { while ((getline < \"/refs/${LOCAL_FAI_NAME}\") > 0) { split(\$0, a, \"\t\"); len[a[1]] = a[2]; } }
                {
                    if (\$1 == \"@SQ\") {
                        for (i=2; i<=NF; i++) {
                            if (\$i ~ /^SN:/) {
                                n=substr(\$i,4); gsub(ref\"#0#\",\"\",n); gsub(ref\"#\",\"\",n);
                                if (n in len) { printf \"@SQ\tSN:%s\tLN:%s\n\", n, len[n]; }
                                else { sub(substr(\$i,4),n); print \$0; }
                                next;
                            }
                        }
                    }
                    print \$0;
                }
            ' /work/old.sam > /work/new.sam
            samtools reheader /work/new.sam /work/tmp.${ref}.bam | samtools sort -@ ${THREADS} -o /work/${SAMPLE}.${ref}.bam
            samtools index /work/${SAMPLE}.${ref}.bam
            rm /work/tmp.${ref}.bam /work/old.sam /work/new.sam
        "
        mark_done "30_Surject_${ref}"
    fi

    if [[ "$DO_VCF_DV" == "true" ]] && ! step_done "40_DV_${ref}"; then
        log "Step 40: DeepVariant (${ref})..."
        docker run --rm -v "${WORK_DIR}":/input -v "${VAR_DIR}":/output -v "${REF_DIR_HOST}":/refs:ro "${DV_IMAGE}" \
            run_deepvariant --model_type=WGS --ref="${d_ref_path}" \
            --reads="/input/${SAMPLE}.${ref}.bam" \
            --output_vcf="/output/${SAMPLE}.${ref}.dv.vcf.gz" \
            --num_shards="${THREADS}"
        mark_done "40_DV_${ref}"
    fi
}

run_sv_calling() {
    local ref=$1
    if [[ "$DO_VCF_SV" == "false" ]]; then return; fi
    local SNARLS_NAME="$(basename "${GBZ_HOST%.gbz}.snarls")"

    if ! step_done "50_Pack_${ref}"; then
        log "Step 50: Packing Graph (RAM Throttled)..."
        docker run --rm -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg pack -x "/graph/$(basename "${GBZ_HOST}")" -g "/work/${SAMPLE}.gam" \
            -t "${MEM_THREADS}" -o "/work/${SAMPLE}.pack"
        mark_done "50_Pack_${ref}"
    fi

    if ! step_done "55_SV_${ref}"; then
        log "Step 55: vg call (Direct from GBZ)..."
        docker run --rm -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph -v "${VAR_DIR}":/variants "${VG_IMAGE}" \
            /bin/sh -c "vg call /graph/$(basename "${GBZ_HOST}") -k /work/${SAMPLE}.pack -r /graph/${SNARLS_NAME} -t ${MEM_THREADS} | bgzip -c > /variants/${SAMPLE}.${ref}.sv.vcf.gz"
        mark_done "55_SV_${ref}"
    fi
}

run_cnvkit() {
    local ref=$1
    if [[ "$DO_VCF_CNV" == "false" ]] || step_done "60_CNVkit_${ref}"; then return; fi
    log "Step 60: CNVkit..."
    local D_REF="/refs/hg38.fa"; [[ "$ref" == "CHM13" ]] && D_REF="/refs/hs1.fa"
    local SEX_FLAG=""; [[ "${SAMPLE_SEX}" == "male" ]] && SEX_FLAG="--male-reference"

    docker run --rm -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/refs:ro -v "${CNV_DIR}":/output "${CNVKIT_IMAGE}" \
        cnvkit.py batch "/work/${SAMPLE}.${ref}.bam" --fasta "${D_REF}" --method wgs -d /output/ ${SEX_FLAG}
    mark_done "60_CNVkit_${ref}"
}

run_multiqc() {
    if [[ "$DO_MULTIQC" == "false" ]]; then return; fi
    log "Step 90: MultiQC..."
    docker run --rm -v "${WORK_DIR}":/data "${MULTIQC_IMAGE}" multiqc /data -o /data/qc
}

# --- 6. MAIN ---
main() {
    log "Starting vg_unified v0.92..."
    download_refs
    rebuild_indexes
    run_fastqc
    run_alignment
    for R in "${REFS[@]}"; do
        log ">>> Ref: $R <<<"
        D_PATH="/refs/hg38.fa"; [[ "$R" == "CHM13" ]] && D_PATH="/refs/hs1.fa"
        run_surjection_and_dv "$R" "$D_PATH"
        run_sv_calling "$R"
        run_cnvkit "$R"
    done
    run_multiqc
    log "All steps complete."
}

main "$@"