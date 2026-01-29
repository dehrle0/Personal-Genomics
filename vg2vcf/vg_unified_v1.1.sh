#!/usr/bin/env bash
#
# vg_unified_v1.1.sh
# Logic: v0.8 Creation + v1.0 reporting + DV 1.10-beta
# Updates: No sudo/root requirements, Modular functions, Granular Toggles

set -euo pipefail

# --- 1. GRANULAR TOGGLES ---
DO_QC_FQ=false      # FastQC (Read quality)
DO_ALIGN=false      # vg giraffe (Mapping)
DO_AUGMENT=true    # v0.8 logic: Create novel variants in graph
DO_BAM=true        # vg surject + UCSC fix (Linear visualization)
DO_VCF_SNP=true    # DeepVariant 1.10-beta (SNPs/Indels)
DO_VCF_SV=true     # vg call -a (Structural Variants)
DO_VCF_CNV=true    # vg call (Copy Number Variants)
DO_STATS=true      # Generate .stats for MultiQC

# --- 2. CONFIGURATION ---
SAMPLE="ME"
SAMPLE_SEX="female"  # male needed for CNVkit on chrY - not used in this script
THREADS="12"
REF_MODE="GRCh38"

# Host Paths (v1.0 Standards)
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
MULTIQC_IMAGE="ewels/multiqc:v1.19"
BCFTOOLS_IMAGE="quay.io/biocontainers/bcftools:1.19--h8b25389_0"
FASTQC_IMAGE="staphb/fastqc:0.12.1"
SAMTOOLS_IMAGE="staphb/samtools:1.19"

# --- 3. UTILITIES & SETUP ---
mkdir -p "${WORK_DIR}" "${LOG_DIR}" "${QC_DIR}"/stats "${QC_DIR}"/fastqc
touch "${COMPLETION_FILE}"

if [[ "$REF_MODE" == "BOTH" ]]; then REFS=("GRCh38" "CHM13"); elif [[ "$REF_MODE" == "CHM13" ]]; then REFS=("CHM13"); else REFS=("GRCh38"); fi

log() { printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$1" | tee -a "${LOG_DIR}/pipeline_${SAMPLE}.log"; }
step_done() { grep -qx "$1" "${COMPLETION_FILE}" 2>/dev/null || return 1; }
mark_done() { echo "$1" >> "${COMPLETION_FILE}"; }

# --- 4. MODULAR FUNCTIONS ---

run_fastqc() {
    if step_done "05_FastQC" || [[ "$DO_QC_FQ" == "false" ]]; then return; fi
    log "Step 05: Running FastQC..."
    docker run --rm --user $(id -u):$(id -g) -v "${READS_DIR_HOST}":/reads:ro -v "${QC_DIR}/fastqc":/output "${FASTQC_IMAGE}" \
        fastqc -t "${THREADS}" -o /output "/reads/$(basename "${FQ1_HOST}")" "/reads/$(basename "${FQ2_HOST}")"
    mark_done "05_FastQC"
}

run_alignment() {
    if step_done "10_Align" || [[ "$DO_ALIGN" == "false" ]]; then return; fi
    log "Step 10: vg giraffe alignment..."
    # We use /work/ inside the container to ensure the file is written directly to the mapped volume
    docker run --rm --user $(id -u):$(id -g) \
        -v "${GRAPH_DIR_HOST}":/graph \
        -v "${READS_DIR_HOST}":/reads:ro \
        -v "${WORK_DIR}":/work "${VG_IMAGE}" \
        /bin/sh -c "vg giraffe -x /graph/$(basename "${GBZ_HOST}") \
            -f /reads/$(basename "${FQ1_HOST}") \
            -f /reads/$(basename "${FQ2_HOST}") \
            -t ${THREADS} -o GAM > /work/${SAMPLE}.gam"
    mark_done "10_Align"
}

run_augmentation() {
    if step_done "20_Augment" || [[ "$DO_AUGMENT" == "false" ]]; then return; fi
    
    local GBZ_NAME=$(basename "${GBZ_HOST}")
    local PG_REF_NAME="${GBZ_NAME%.gbz}.pg"
    
    log "Step 20 (v0.8): Preparing Mutable Graph and Augmenting..."

    # 1. Convert GBZ to PG (Modern vg convert logic)
    if [ ! -f "${GRAPH_DIR_HOST}/${PG_REF_NAME}" ]; then
        log "Converting reference GBZ to PG format for augmentation..."
        # We use vg convert -p to explicitly target PackedGraph format
        docker run --rm --user $(id -u):$(id -g) \
            -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            vg convert -p "/graph/${GBZ_NAME}" > "${GRAPH_DIR_HOST}/${PG_REF_NAME}"
    fi

    # 2. Run Augmentation using the Mutable PG graph
    log "Augmenting graph with sample reads..."
    docker run --rm --user $(id -u):$(id -g) \
        -v "${WORK_DIR}":/work \
        -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        /bin/sh -c "vg augment /graph/${PG_REF_NAME} /work/${SAMPLE}.gam -t ${THREADS} -A /work/${SAMPLE}.aug.gam > /work/${SAMPLE}.aug.pg"
    
    mark_done "20_Augment"
}

run_surjection() {
    local ref=$1
    local step_id="30_Surject_${ref}"
    if step_done "$step_id" || [[ "$DO_BAM" == "false" ]]; then return; fi
    log "Step 30: Surjecting to $ref (Surgical Header Fix)..."

    # Define the correct FAI path for this reference
    local LOCAL_FAI="${HG38_LOCAL}.fai"
    [[ "$ref" == "CHM13" ]] && LOCAL_FAI="${HS1_LOCAL}.fai"

    # 1. Surject to temp BAM
    docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        vg surject -x "/graph/$(basename "${GBZ_HOST}")" "/work/${SAMPLE}.gam" -b -p "${ref}" -t "${THREADS}" > "${WORK_DIR}/tmp.${ref}.bam"

    # 2. Run your proven AWK-based Header Manipulation
    docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro "${SAMTOOLS_IMAGE}" /bin/bash -c "
        samtools view -H /work/tmp.${ref}.bam > /work/old_header.sam
        
        awk -v ref=\"${ref}\" '
            BEGIN {
                while ((getline < \"/references/$(basename "${LOCAL_FAI}")\") > 0) {
                    split(\$0, a, \"\t\");
                    len[a[1]] = a[2];
                }
            }
            {
                if (\$1 == \"@SQ\") {
                    for (i=2; i<=NF; i++) {
                        if (\$i ~ /^SN:/) {
                            orig_name = substr(\$i, 4);
                            clean_name = orig_name;
                            gsub(ref \"#0#\", \"\", clean_name);
                            gsub(ref \"#\", \"\", clean_name);
                            
                            if (clean_name in len) {
                                printf \"@SQ\tSN:%s\tLN:%s\n\", clean_name, len[clean_name];
                            } else {
                                # Fallback if not in FAI, still clean the name
                                sub(orig_name, clean_name);
                                print \$0;
                            }
                            next;
                        }
                    }
                }
                print \$0;
            }
        ' /work/old_header.sam > /work/new_header.sam && \
        
        samtools reheader /work/new_header.sam /work/tmp.${ref}.bam | \
        samtools sort -@ ${THREADS} -o /work/${SAMPLE}.${ref}.bam && \
        samtools index /work/${SAMPLE}.${ref}.bam && \
        rm /work/tmp.${ref}.bam /work/old_header.sam /work/new_header.sam
    "
    
    if [[ "$DO_STATS" == "true" ]]; then
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${SAMTOOLS_IMAGE}" \
            samtools stats -@ "${THREADS}" "/work/${SAMPLE}.${ref}.bam" > "${QC_DIR}/stats/${SAMPLE}.${ref}.bam.stats"
    fi
    mark_done "$step_id"
}

run_vcf_calls() {
    local ref=$1
    local d_ref_path=$2
    
    # SNP (DeepVariant)
    if [[ "$DO_VCF_SNP" == "true" ]] && ! step_done "40_SNP_${ref}"; then
        log "Step 40: DeepVariant SNPs for $ref..."
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/input -v "${REF_DIR_HOST}":/references:ro "${DV_IMAGE}" \
            run_deepvariant --model_type=WGS --ref="${d_ref_path}" --reads="/input/${SAMPLE}.${ref}.bam" \
            --output_vcf="/input/${SAMPLE}.${ref}.snps.vcf.gz" --num_shards="${THREADS}"
        mark_done "40_SNP_${ref}"
    fi

    # Pack for SV/CNV
    if ! step_done "50_Pack_${ref}"; then
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${VG_IMAGE}" \
            vg pack -x "/work/${SAMPLE}.aug.pg" -g "/work/${SAMPLE}.aug.gam" -t "${THREADS}" -o "/work/${SAMPLE}.${ref}.pack"
        mark_done "50_Pack_${ref}"
    fi

    # SV Calling
    if [[ "$DO_VCF_SV" == "true" ]] && ! step_done "60_SV_${ref}"; then
        log "Step 60: vg call SVs for $ref..."
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${VG_IMAGE}" \
            /bin/sh -c "vg call /work/${SAMPLE}.aug.pg -k /work/${SAMPLE}.${ref}.pack -S ${ref} -a -t ${THREADS} | bgzip -c > /work/${SAMPLE}.${ref}.sv.vcf.gz"
        mark_done "60_SV_${ref}"
    fi

    # CNV Calling
    if [[ "$DO_VCF_CNV" == "true" ]] && ! step_done "70_CNV_${ref}"; then
        log "Step 70: vg call CNVs for $ref..."
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${VG_IMAGE}" \
            /bin/sh -c "vg call /work/${SAMPLE}.aug.pg -k /work/${SAMPLE}.${ref}.pack -S ${ref} -t ${THREADS} | awk '\$0 ~ /^#/ || \$0 ~ /TYPE=CNV/ || \$0 ~ /SVTYPE=DUP/ || \$0 ~ /SVTYPE=DEL/' | bgzip -c > /work/${SAMPLE}.${ref}.cnv.vcf.gz"
        mark_done "70_CNV_${ref}"
    fi
    
    # VCF Stats
    if [[ "$DO_STATS" == "true" ]]; then
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "bcftools stats /work/${SAMPLE}.${ref}.snps.vcf.gz > /qc/stats/${SAMPLE}.${ref}.snps.stats && \
                        bcftools stats /work/${SAMPLE}.${ref}.sv.vcf.gz > /qc/stats/${SAMPLE}.${ref}.sv.stats && \
                        bcftools stats /work/${SAMPLE}.${ref}.cnv.vcf.gz > /qc/stats/${SAMPLE}.${ref}.cnv.stats"
    fi
}

# --- 5. MAIN SECTION ---
main() {
    log "Starting Main Pipeline Flow..."
    
    
    run_fastqc
    run_alignment
    run_augmentation # v0.8 logic
    
    for R in "${REFS[@]}"; do
        log "Processing Reference: $R"
        D_PATH="/references/hg38.fa"; [[ "$R" == "CHM13" ]] && D_PATH="/references/hs1.fa"
        
        run_surjection "$R"
        run_vcf_calls "$R" "$D_PATH"
    done
    
    log "All steps complete. Results in ${WORK_DIR} and stats in ${QC_DIR}/stats."
}

main "$@"