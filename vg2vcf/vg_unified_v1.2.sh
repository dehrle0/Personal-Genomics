#!/usr/bin/env bash
#
# vg_unified_v1.2.sh
# Logic: v0.8 Augment + v1.0 AWK Header + vg 1.70 Index/Version handling
# Structure: Modular functions with non-root Docker execution

set -euo pipefail

# --- 1. GRANULAR TOGGLES ---
DO_DOWNLOAD=false   # wget HPRC references if missing
DO_INDEX=false      # Set to false if you already have .dist and .min
DO_ALIGN=false      # Set to false if you already have the .gam
DO_AUGMENT=true     # v0.8 logic: Create novel variants in graph
DO_BAM=true         # vg surject + Surgical AWK Fix
DO_VCF_SNP=true     # DeepVariant 1.10-beta
DO_VCF_SV=true      # vg call -a (Structural Variants)
DO_VCF_CNV=true     # vg call (Copy Number Variants)
DO_STATS=true       # Generate .stats for MultiQC

# --- 2. CONFIGURATION ---
SAMPLE="ME"
THREADS="12"
REF_MODE="GRCh38"  # Options: GRCh38, CHM13, BOTH

# Host Paths
REF_DIR_HOST="/mnt/data/work/references"
HG38_LOCAL="${REF_DIR_HOST}/hg38.fa"
HS1_LOCAL="${REF_DIR_HOST}/chm13v2.0_maskedY_rCRS.fa"
GRAPH_DIR_HOST="${REF_DIR_HOST}/hrpc"
GBZ_HOST="${GRAPH_DIR_HOST}/hprc-v1.1-mc-chm13.d9.gbz"
READS_DIR_HOST="/mnt/c/Genomes/${SAMPLE}/Data/Source"
FQ1_HOST="${READS_DIR_HOST}/MelindaEhrle-SQXB2888-30x-WGS-Sequencing_com-2024-08-18.1.fq.gz"
FQ2_HOST="${READS_DIR_HOST}/MelindaEhrle-SQXB2888-30x-WGS-Sequencing_com-2024-08-18.2.fq.gz"

# Working Directories
WORK_DIR="/mnt/data/work/vg_pangenome"
LOG_DIR="${WORK_DIR}/logs"
QC_DIR="${WORK_DIR}/qc"
COMPLETION_FILE="${WORK_DIR}/completed_steps_${SAMPLE}.txt"

# Docker Images
VG_IMAGE="quay.io/vgteam/vg:v1.70.0"
DV_IMAGE="google/deepvariant:1.10.0-beta"
BCFTOOLS_IMAGE="quay.io/biocontainers/bcftools:1.19--h8b25389_0"
SAMTOOLS_IMAGE="staphb/samtools:1.19"

# --- 3. UTILITIES & SETUP ---
mkdir -p "${WORK_DIR}" "${LOG_DIR}" "${QC_DIR}"/stats
touch "${COMPLETION_FILE}"

if [[ "$REF_MODE" == "BOTH" ]]; then REFS=("GRCh38" "CHM13"); elif [[ "$REF_MODE" == "CHM13" ]]; then REFS=("CHM13"); else REFS=("GRCh38"); fi

log() { printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$1" | tee -a "${LOG_DIR}/pipeline_${SAMPLE}.log"; }
step_done() { grep -qx "$1" "${COMPLETION_FILE}" 2>/dev/null || return 1; }
mark_done() { echo "$1" >> "${COMPLETION_FILE}"; }

# --- 4. PRE-FLIGHT FUNCTIONS ---

download_refs() {
    if [[ "$DO_DOWNLOAD" == "false" ]]; then return; fi
    log "Checking for HPRC references..."
    mkdir -p "${GRAPH_DIR_HOST}"
    if [ ! -f "${GBZ_HOST}" ]; then
        wget -c -P "${GRAPH_DIR_HOST}" https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/release/v1.1/mc_chm13/hprc-v1.1-mc-chm13.d9.gbz
    fi
}

rebuild_indexes() {
    if [[ "$DO_INDEX" == "false" ]]; then return; fi
    local BASE="${GBZ_HOST%.gbz}"
    log "Rebuilding .dist and .min indexes for vg 1.70 compatibility..."
    docker run --rm --user $(id -u):$(id -g) -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        vg index -t "${THREADS}" -d "/graph/$(basename "${BASE}").dist" "/graph/$(basename "${GBZ_HOST}")"
    docker run --rm --user $(id -u):$(id -g) -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        vg index -t "${THREADS}" -m "/graph/$(basename "${BASE}").min" "/graph/$(basename "${GBZ_HOST}")"
    mark_done "01_Index"
}

# --- 5. CORE PIPELINE FUNCTIONS ---

run_alignment() {
    if step_done "10_Align" || [[ "$DO_ALIGN" == "false" ]]; then return; fi
    log "Step 10: Alignment (vg giraffe)..."
    docker run --rm --user $(id -u):$(id -g) -v "${GRAPH_DIR_HOST}":/graph -v "${READS_DIR_HOST}":/reads:ro -v "${WORK_DIR}":/work "${VG_IMAGE}" \
        /bin/sh -c "vg giraffe -x /graph/$(basename "${GBZ_HOST}") -f /reads/$(basename "${FQ1_HOST}") -f /reads/$(basename "${FQ2_HOST}") -t ${THREADS} -o GAM > /work/${SAMPLE}.gam"
    mark_done "10_Align"
}

run_augmentation() {
    if step_done "20_Augment" || [[ "$DO_AUGMENT" == "false" ]]; then return; fi
    
    local GBZ_NAME=$(basename "${GBZ_HOST}")
    local HG_REF_NAME="${GBZ_NAME%.gbz}.hg" # Changed from .pg to .hg

    log "Step 20: Preparing HashGraph (Mutable) and Augmenting..."
    
    # 1. Convert GBZ to HG (HashGraph is the required 'MutablePathMutableHandleGraph' type)
    if [ ! -f "${GRAPH_DIR_HOST}/${HG_REF_NAME}" ]; then
        log "Converting reference GBZ to HashGraph format for augmentation..."
        docker run --rm --user $(id -u):$(id -g) -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg convert -H "/graph/${GBZ_NAME}" > "${GRAPH_DIR_HOST}/${HG_REF_NAME}"
    fi

    # 2. Perform Augmentation using the .hg reference
    # Note: Output is still saved as .aug.pg (PackedGraph) for efficiency in calling
    docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        /bin/sh -c "vg augment /graph/${HG_REF_NAME} /work/${SAMPLE}.gam -t ${THREADS} -A /work/${SAMPLE}.aug.gam > /work/${SAMPLE}.aug.pg"
    
    mark_done "20_Augment"
}

run_surjection() {
    local ref=$1
    if step_done "30_Surject_${ref}" || [[ "$DO_BAM" == "false" ]]; then return; fi
    log "Step 30: Surjecting to $ref with AWK Header Fix..."

    local LOCAL_FAI="${REF_DIR_HOST}/hg38.fa.fai"
    [[ "$ref" == "CHM13" ]] && LOCAL_FAI="${REF_DIR_HOST}/chm13v2.0_maskedY_rCRS.fa.fai"

    docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        vg surject -x "/graph/$(basename "${GBZ_HOST}")" "/work/${SAMPLE}.gam" -b -p "${ref}" -t "${THREADS}" > "${WORK_DIR}/tmp.${ref}.bam"

    docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro "${SAMTOOLS_IMAGE}" /bin/bash -c "
        samtools view -H /work/tmp.${ref}.bam > /work/old.sam
        awk -v ref=\"${ref}\" '
            BEGIN { while ((getline < \"/references/$(basename "${LOCAL_FAI}")\") > 0) { split(\$0, a, \"\t\"); len[a[1]] = a[2]; } }
            {
                if (\$1 == \"@SQ\") {
                    for (i=2; i<=NF; i++) {
                        if (\$i ~ /^SN:/) {
                            n = substr(\$i, 4); gsub(ref \"#0#\", \"\", n); gsub(ref \"#\", \"\", n);
                            if (n in len) { printf \"@SQ\tSN:%s\tLN:%s\n\", n, len[n]; }
                            else { sub(substr(\$i, 4), n); print \$0; }
                            next;
                        }
                    }
                }
                print \$0;
            }
        ' /work/old.sam > /work/new.sam
        samtools reheader /work/new.sam /work/tmp.${ref}.bam | samtools sort -@ ${THREADS} -o /work/${SAMPLE}.${ref}.bam
        samtools index /work/${SAMPLE}.${ref}.bam && rm /work/tmp.${ref}.bam /work/old.sam /work/new.sam
    "
    mark_done "30_Surject_${ref}"
}

run_vcf_calls() {
    local ref=$1
    local d_ref_path=$2
    
    if [[ "$DO_VCF_SNP" == "true" ]] && ! step_done "40_SNP_${ref}"; then
        log "Step 40: DeepVariant SNPs ($ref)..."
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/input -v "${REF_DIR_HOST}":/references:ro "${DV_IMAGE}" \
            run_deepvariant --model_type=WGS --ref="${d_ref_path}" --reads="/input/${SAMPLE}.${ref}.bam" --output_vcf="/input/${SAMPLE}.${ref}.snps.vcf.gz" --num_shards="${THREADS}"
        mark_done "40_SNP_${ref}"
    fi

    if ! step_done "50_Pack_${ref}"; then
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${VG_IMAGE}" \
            vg pack -x "/work/${SAMPLE}.aug.pg" -g "/work/${SAMPLE}.aug.gam" -t "${THREADS}" -o "/work/${SAMPLE}.${ref}.pack"
        mark_done "50_Pack_${ref}"
    fi

    if [[ "$DO_VCF_SV" == "true" ]] && ! step_done "60_SV_${ref}"; then
        log "Step 60: SV Call ($ref)..."
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${VG_IMAGE}" \
            /bin/sh -c "vg call /work/${SAMPLE}.aug.pg -k /work/${SAMPLE}.${ref}.pack -S ${ref} -a -t ${THREADS} | bgzip -c > /work/${SAMPLE}.${ref}.sv.vcf.gz"
        mark_done "60_SV_${ref}"
    fi

    if [[ "$DO_VCF_CNV" == "true" ]] && ! step_done "70_CNV_${ref}"; then
        log "Step 70: CNV Call ($ref)..."
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${VG_IMAGE}" \
            /bin/sh -c "vg call /work/${SAMPLE}.aug.pg -k /work/${SAMPLE}.${ref}.pack -S ${ref} -t ${THREADS} | awk '\$0 ~ /^#/ || \$0 ~ /TYPE=CNV/ || \$0 ~ /SVTYPE=DUP/ || \$0 ~ /SVTYPE=DEL/' | bgzip -c > /work/${SAMPLE}.${ref}.cnv.vcf.gz"
        mark_done "70_CNV_${ref}"
    fi
}

# --- 6. MAIN ---
main() {
    log "Starting vg_unified v1.2 Flow..."
    download_refs
    rebuild_indexes
    run_alignment
    run_augmentation
    
    for R in "${REFS[@]}"; do
        D_PATH="/references/hg38.fa"; [[ "$R" == "CHM13" ]] && D_PATH="/references/hs1.fa"
        run_surjection "$R"
        run_vcf_calls "$R" "$D_PATH"
    done
    log "Pipeline Complete."
}

main "$@"