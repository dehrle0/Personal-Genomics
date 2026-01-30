#!/usr/bin/env bash
# vg_unified_v0.94.sh
# Separated Logic & Config | Memory Optimized | Auto-Cleanup

set -euo pipefail

# --- 1. LOAD CONFIGURATION ---
CONFIG_FILE="./pipeline.env"
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
    echo "Loaded configuration from $CONFIG_FILE"
else
    echo "Error: $CONFIG_FILE not found. Please create it based on your settings."
    exit 1
fi

# --- 2. DERIVED VARIABLES & PATHS ---
GBZ_HOST="${GRAPH_DIR_HOST}/${GBZ_NAME}"
HAPL_HOST="${GRAPH_DIR_HOST}/${HAPL_NAME}"
HG38_LOCAL="${REF_DIR_HOST}/hg38.fa"
HS1_LOCAL="${REF_DIR_HOST}/chm13v2.0_maskedY_rCRS.fa"

# Detect FastQs - Assuming standard naming if not explicitly set in env
FQ1_HOST="${READS_DIR_HOST}/${FQ1_NAME:-MelindaEhrle-SQXB2888-30x-WGS-Sequencing_com-2024-08-18.1.fq.gz}"
FQ2_HOST="${READS_DIR_HOST}/${FQ2_NAME:-MelindaEhrle-SQXB2888-30x-WGS-Sequencing_com-2024-08-18.2.fq.gz}"

WORK_DIR="/mnt/data/work/vg_pangenome"
LOG_DIR="${WORK_DIR}/logs"
QC_DIR="${WORK_DIR}/qc"
CNV_DIR="${WORK_DIR}/cnv"
VAR_DIR="${WORK_DIR}/variants"
COMPLETION_FILE="${WORK_DIR}/completed_steps_${SAMPLE}.txt"

# Docker Images
VG_IMAGE="quay.io/vgteam/vg:v1.70.0"
DV_IMAGE="google/deepvariant:1.10.0-beta"
CNVKIT_IMAGE="etal/cnvkit:0.9.10"
SAMTOOLS_IMAGE="staphb/samtools:1.19"
FASTQC_IMAGE="staphb/fastqc:0.12.1"
MULTIQC_IMAGE="ewels/multiqc:v1.19"

# Determine Reference List
if [[ "$REF_MODE" == "BOTH" ]]; then REFS=("GRCh38" "CHM13"); elif [[ "$REF_MODE" == "CHM13" ]]; then REFS=("CHM13"); else REFS=("GRCh38"); fi

# --- 3. UTILITIES ---
mkdir -p "${WORK_DIR}" "${LOG_DIR}" "${QC_DIR}/fastqc" "${CNV_DIR}" "${VAR_DIR}"
touch "${COMPLETION_FILE}"

log() { printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$1" | tee -a "${LOG_DIR}/pipeline_${SAMPLE}.log"; }
step_done() { grep -qx "$1" "${COMPLETION_FILE}" 2>/dev/null || return 1; }
mark_done() { echo "$1" >> "${COMPLETION_FILE}"; }

# --- 4. INDEXING & PRE-FLIGHT ---

rebuild_indexes() {
    if [[ "$DO_INDEX" == "false" ]]; then return; fi
    local BASE=$(basename "${GBZ_HOST}")
    local DIST="${GBZ_HOST%.gbz}.dist"
    local MIN="${GBZ_HOST%.gbz}.min"
    local SNARLS="${GBZ_HOST%.gbz}.snarls"

    if [[ ! -f "$DIST" ]]; then
        log "Index: Building .dist (vg index -j)..."
        docker run --rm -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" vg index -t "${THREADS}" -j "/graph/$(basename "$DIST")" "/graph/$BASE"
    fi

    if [[ ! -f "$MIN" ]]; then
        log "Index: Building .min (vg minimizer -d)..."
        docker run --rm -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" vg minimizer -t "${THREADS}" -d "/graph/$(basename "$DIST")" -o "/graph/$(basename "$MIN")" "/graph/$BASE"
    fi

    if [[ ! -f "$SNARLS" ]]; then
        log "Index: Building .snarls (RAM Throttled to ${MEM_THREADS} threads)..."
        docker run --rm -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" vg snarls -t "${MEM_THREADS}" "/graph/$BASE" > "$SNARLS"
    fi
}

# --- 5. CORE PIPELINE ---

run_alignment() {
    if [[ "$DO_ALIGN" == "false" ]] || step_done "10_Align"; then return; fi
    log "Step 10: vg giraffe alignment..."

    local giraffe_cmd=(vg giraffe
        -Z "/graph/$(basename "${GBZ_HOST}")"
        -m "/graph/$(basename "${GBZ_HOST%.gbz}.min")"
        -d "/graph/$(basename "${GBZ_HOST%.gbz}.dist")"
        -f "/reads/$(basename "${FQ1_HOST}")" 
        -f "/reads/$(basename "${FQ2_HOST}")"
        -t "${THREADS}" -o GAM)

    if [[ "$USE_HAPL" == "true" ]] && [[ -n "$HAPL_NAME" ]]; then
        log "Haplotype mode enabled: Using $HAPL_NAME"
        giraffe_cmd+=(-H "/graph/${HAPL_NAME}")
    fi

    docker run --rm -v "${GRAPH_DIR_HOST}":/graph -v "${READS_DIR_HOST}":/reads:ro -v "${WORK_DIR}":/work "${VG_IMAGE}" \
        "${giraffe_cmd[@]}" > "${WORK_DIR}/${SAMPLE}.gam"
    mark_done "10_Align"
}

run_surjection_and_dv() {
    local ref=$1
    local d_ref_path=$2
    if [[ "$DO_SURJECT" == "false" ]]; then return; fi
    
    if ! step_done "30_Surject_${ref}"; then
        log "Step 30: Surjecting and Repairing BAM for ${ref}..."
        local LOCAL_FAI_NAME="$(basename "${HG38_LOCAL}.fai")"; [[ "$ref" == "CHM13" ]] && LOCAL_FAI_NAME="$(basename "${HS1_LOCAL}.fai")"

        docker run --rm -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg surject -x "/graph/$(basename "${GBZ_HOST}")" -b -t "${THREADS}" --into-ref "${ref}" "/work/${SAMPLE}.gam" > "${WORK_DIR}/tmp.${ref}.bam"

        docker run --rm -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/refs:ro "${SAMTOOLS_IMAGE}" /bin/bash -c "
            samtools view -H /work/tmp.${ref}.bam > /work/old.sam
            awk -v ref=\"${ref}\" 'BEGIN { while ((getline < \"/refs/${LOCAL_FAI_NAME}\") > 0) { split(\$0, a, \"\t\"); len[a[1]] = a[2]; } } { if (\$1 == \"@SQ\") { for (i=2; i<=NF; i++) { if (\$i ~ /^SN:/) { n=substr(\$i,4); gsub(ref\"#0#\",\"\",n); gsub(ref\"#\",\"\",n); if (n in len) { printf \"@SQ\tSN:%s\tLN:%s\n\", n, len[n]; } else { sub(substr(\$i,4),n); print \$0; } next; } } } print \$0; }' /work/old.sam > /work/new.sam
            samtools reheader /work/new.sam /work/tmp.${ref}.bam | samtools sort -@ ${THREADS} -o /work/${SAMPLE}.${ref}.bam
            samtools index /work/${SAMPLE}.${ref}.bam && rm /work/tmp.${ref}.bam /work/old.sam /work/new.sam"
        mark_done "30_Surject_${ref}"
    fi

    if [[ "$DO_VCF_DV" == "true" ]] && ! step_done "40_DV_${ref}"; then
        log "Step 40: DeepVariant (${ref})..."
        docker run --rm -v "${WORK_DIR}":/input -v "${VAR_DIR}":/output -v "${REF_DIR_HOST}":/refs:ro "${DV_IMAGE}" run_deepvariant --model_type=WGS --ref="${d_ref_path}" --reads="/input/${SAMPLE}.${ref}.bam" --output_vcf="/output/${SAMPLE}.${ref}.dv.vcf.gz" --num_shards="${THREADS}"
        mark_done "40_DV_${ref}"
    fi
}

run_sv_calling() {
    local ref=$1
    if [[ "$DO_VCF_SV" == "false" ]]; then return; fi
    local SNARLS_NAME="$(basename "${GBZ_HOST%.gbz}.snarls")"

    if ! step_done "50_Pack_${ref}"; then
        log "Step 50: Packing Graph (RAM Throttled: ${MEM_THREADS})..."
        docker run --rm -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" vg pack -x "/graph/$(basename "${GBZ_HOST}")" -g "/work/${SAMPLE}.gam" -t "${MEM_THREADS}" -o "/work/${SAMPLE}.pack"
        mark_done "50_Pack_${ref}"
    fi

    if ! step_done "55_SV_${ref}"; then
        log "Step 55: vg call (${ref})..."
        docker run --rm -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph -v "${VAR_DIR}":/variants "${VG_IMAGE}" /bin/sh -c "vg call /graph/$(basename "${GBZ_HOST}") -k /work/${SAMPLE}.pack -r /graph/${SNARLS_NAME} -S ${ref} -z -t ${MEM_THREADS} | bgzip -c > /variants/${SAMPLE}.${ref}.sv.vcf.gz"
        rm -f "${WORK_DIR}/${SAMPLE}.pack"
        mark_done "55_SV_${ref}"
    fi
}

run_cnvkit() {
    local ref=$1
    if [[ "$DO_VCF_CNV" == "false" ]] || step_done "60_CNVkit_${ref}"; then return; fi
    log "Step 60: CNVkit (${ref})..."
    
    local D_REF="/refs/hg38.fa"; [[ "$ref" == "CHM13" ]] && D_REF="/refs/hs1.fa"
    local SEX_FLAG=""; [[ "${SAMPLE_SEX}" == "male" ]] && SEX_FLAG="--male-reference"
    local REF_CNN_BASE="${SAMPLE}.${ref}.reference.cnn"

    docker run --rm -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/refs:ro -v "${CNV_DIR}":/output "${CNVKIT_IMAGE}" \
        cnvkit.py batch "/work/${SAMPLE}.${ref}.bam" --fasta "${D_REF}" --method wgs --output-reference "/output/${REF_CNN_BASE}" -d /output/ -n ${SEX_FLAG}
    
    rm -f "${CNV_DIR}/${REF_CNN_BASE}"
    mark_done "60_CNVkit_${ref}"
}

# --- 6. MAIN ---
main() {
    log "Starting vg_unified v0.94..."
    rebuild_indexes
    run_alignment
    for R in "${REFS[@]}"; do
        log ">>> Processing Reference: $R <<<"
        D_PATH="/refs/hg38.fa"; [[ "$R" == "CHM13" ]] && D_PATH="/refs/hs1.fa"
        run_surjection_and_dv "$R" "$D_PATH"
        run_sv_calling "$R"
        run_cnvkit "$R"
    done
    [[ "$DO_MULTIQC" == "true" ]] && log "Step 90: MultiQC..." && docker run --rm -v "${WORK_DIR}":/data "${MULTIQC_IMAGE}" multiqc /data -o /data/qc
    log "All steps complete for $SAMPLE."
}

main "$@"