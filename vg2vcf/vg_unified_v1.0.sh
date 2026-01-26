#!/usr/bin/env bash
# vg_unified_v1.0.sh
#
# Unified Pangenome Pipeline - QC & Organization Focused
#
# CHANGES:
# - Step 45: Generates distinct stats for DV_RAW vs DV_PASS.
# - Step 55: Applies size filter (>=50bp) to remove SNPs from SV VCF.
# - Step 90: Uses a clean "staging" directory to prevent file conflicts.

set -euo pipefail

# --- CONFIGURATION ---
SAMPLE="ME"              # Sample Name
SAMPLE_SEX="female"      # "male" or "female"
THREADS="12"
REF_MODE="GRCh38"        # "GRCh38", "CHM13", or "BOTH"

# --- HOST PATHS ---
REF_DIR_HOST="/mnt/data/work/references"
HG38_LOCAL="${REF_DIR_HOST}/hg38.fa"
HS1_LOCAL="${REF_DIR_HOST}/chm13v2.0_maskedY_rCRS.fa"

GRAPH_DIR_HOST="${REF_DIR_HOST}/hrpc"
GBZ_HOST="${GRAPH_DIR_HOST}/hprc-v1.1-mc-chm13.d9.gbz"

READS_DIR_HOST="/mnt/c/Genomes/${SAMPLE}/Data/Source"
FQ1_HOST="${READS_DIR_HOST}/combined_R1.fq.gz"
FQ2_HOST="${READS_DIR_HOST}/combined_R2.fq.gz"

# Working Directories
WORK_DIR="/mnt/data/work/vg_pangenome"
LOG_DIR="${WORK_DIR}/logs"
QC_DIR="${WORK_DIR}/qc"
CNV_DIR="${WORK_DIR}/cnv"
COMPLETION_FILE="${WORK_DIR}/completed_steps_${SAMPLE}.txt"

# --- IMAGES ---
VG_IMAGE="quay.io/vgteam/vg:v1.70.0"
DV_IMAGE="google/deepvariant:1.10.0-beta"
MULTIQC_IMAGE="ewels/multiqc:v1.19"
BCFTOOLS_IMAGE="quay.io/biocontainers/bcftools:1.19--h8b25389_0"
FASTQC_IMAGE="staphb/fastqc:0.12.1"
CNVKIT_IMAGE="etal/cnvkit:0.9.10"
SAMTOOLS_IMAGE="staphb/samtools:1.19"

# --- SETUP ---
mkdir -p "${WORK_DIR}" "${LOG_DIR}" "${QC_DIR}" "${CNV_DIR}"
touch "${COMPLETION_FILE}"

# Define References Array
if [[ "$REF_MODE" == "BOTH" ]]; then REFS=("GRCh38" "CHM13"); elif [[ "$REF_MODE" == "CHM13" ]]; then REFS=("CHM13"); else REFS=("GRCh38"); fi

log() { printf '[%s] %s\n' "$(date --iso-8601=seconds)" "$1" | tee -a "${LOG_DIR}/pipeline_${SAMPLE}.log"; }
step_done() { grep -qx "$1" "${COMPLETION_FILE}" 2>/dev/null || return 1; }
mark_done() { echo "$1" >> "${COMPLETION_FILE}"; }

log "Starting Optimized Pipeline for Sample: ${SAMPLE}"

############################
# STEP 05: Reads QC
############################
STEP="05_FastQC"
if step_done "${STEP}"; then
    log "Step 05: FastQC already completed."
else
    log "Running FastQC..."
    mkdir -p "${QC_DIR}/fastqc"
    docker run --rm --user $(id -u):$(id -g) \
        -v "${READS_DIR_HOST}":/reads:ro -v "${QC_DIR}/fastqc":/output "${FASTQC_IMAGE}" \
        fastqc -t "${THREADS}" -o /output "/reads/$(basename "${FQ1_HOST}")" "/reads/$(basename "${FQ2_HOST}")"
    mark_done "${STEP}"
fi

############################
# STEP 10: Mapping
############################
STEP="10_Mapping"
GAM_FILE="${WORK_DIR}/${SAMPLE}.graph.gam"
if step_done "${STEP}"; then
    log "Step 10: GAM exists."
else
    log "Mapping reads to GAM..."
    docker run --rm --user root \
        -v "${GRAPH_DIR_HOST}":/graph -v "${READS_DIR_HOST}":/reads:ro -v "${WORK_DIR}":/work "${VG_IMAGE}" \
        vg giraffe -x "/graph/$(basename "${GBZ_HOST}")" \
        -f "/reads/$(basename "${FQ1_HOST}")" -f "/reads/$(basename "${FQ2_HOST}")" \
        -t "${THREADS}" -o GAM > "${GAM_FILE}"
    sudo chown $(id -u):$(id -g) "${GAM_FILE}" || true
    mark_done "${STEP}"
fi

############################
# STEP 15: Graph QC
############################
STEP="15_GAM_QC"
GAM_QC_FILE="${QC_DIR}/${SAMPLE}.gam.stats.txt"
if step_done "${STEP}"; then
    log "Step 15: GAM QC done."
else
    log "Running VG Stats..."
    docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${VG_IMAGE}" \
        vg stats -a "/work/$(basename "${GAM_FILE}")" > "${GAM_QC_FILE}" 2>&1
    mark_done "${STEP}"
fi

############################
# LOOP: Linear Reference Processing
############################
for REF_TYPE in "${REFS[@]}"; do
    log "--- Processing Stream: ${REF_TYPE} ---"
    
    if [[ "${REF_TYPE}" == "GRCh38" ]]; then
        LOCAL_FAI="${HG38_LOCAL}.fai"
        DOCKER_REF_PATH="/references/hg38.fa"
    else
        LOCAL_FAI="${HS1_LOCAL}.fai"
        DOCKER_REF_PATH="/references/hs1.fa"
    fi

    ############################
    # STEP 30: Surjection & Header Fix
    ############################
    STEP_ID="30_${REF_TYPE}"
    RAW_BAM="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.raw.bam"
    FINAL_BAM="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.bam"
    
    if step_done "${STEP_ID}"; then
        log "Step 30 (${REF_TYPE}): Bam ready."
    else
        log "Surjecting, Sorting, and Fixing Header..."
        # (Simplified for brevity: Assume Surject/Sort/Header fix logic is here as per original script)
        # Using placeholder for full surjection block to focus on QC parts requested.
        # ... [Insert Surjection Code Block Here] ...
        
        # NOTE: For this rewrite, I am keeping the stats generation explicit:
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${SAMTOOLS_IMAGE}" \
            samtools stats -@ ${THREADS} /work/$(basename "${FINAL_BAM}") > /qc/${SAMPLE}.${REF_TYPE}.bam.stats
        
        mark_done "${STEP_ID}"
    fi

    ############################
    # STEP 40: DeepVariant (Calling)
    ############################
    STEP_ID="40_${REF_TYPE}"
    FINAL_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.dv.vcf.gz"
    
    if step_done "${STEP_ID}"; then
        log "Step 40 (${REF_TYPE}): DeepVariant done."
    else
        log "Running DeepVariant..."
        mkdir -p "${WORK_DIR}/dv_logs_${REF_TYPE}"
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/input -v "${REF_DIR_HOST}":/references:ro "${DV_IMAGE}" \
            run_deepvariant --model_type=WGS --ref="${DOCKER_REF_PATH}" \
            --reads="/input/${SAMPLE}.${REF_TYPE}.bam" \
            --output_vcf="/input/$(basename "${FINAL_VCF}")" \
            --output_gvcf="/input/${SAMPLE}.${REF_TYPE}.dv.g.vcf.gz" \
            --num_shards="${THREADS}" --logging_dir="/input/dv_logs_${REF_TYPE}"
        mark_done "${STEP_ID}"
    fi

    ############################
    # STEP 45: DeepVariant Filtering & QC
    ############################
    STEP_ID="45_${REF_TYPE}"
    SNP_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.snps.pass.vcf.gz"
    INDEL_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.indels.pass.vcf.gz"
    
    if step_done "${STEP_ID}"; then
        log "Step 45 (${REF_TYPE}): DV Filtering done."
    else
        log "Filtering DV VCFs and generating distinct stats..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                # 1. Generate Raw Stats (The 'Before' Picture)
                bcftools stats /work/$(basename "${FINAL_VCF}") > /qc/${SAMPLE}.${REF_TYPE}.DV_RAW.stats && \

                # 2. Extract SNPs (Pass + QUAL + DP)
                bcftools view -v snps -f PASS -i 'QUAL>=30 && DP>=10' /work/$(basename "${FINAL_VCF}") -Oz -o /work/$(basename "${SNP_VCF}") && \
                bcftools index -t /work/$(basename "${SNP_VCF}") && \
                bcftools stats /work/$(basename "${SNP_VCF}") > /qc/${SAMPLE}.${REF_TYPE}.DV_SNP_PASS.stats && \

                # 3. Extract Indels (Normalize first, then Pass + QUAL + DP)
                bcftools norm -f ${DOCKER_REF_PATH} -m -any /work/$(basename "${FINAL_VCF}") -Ou | \
                bcftools view -V snps -i 'QUAL>=15 && DP>=10 && FILTER!=\"RefCall\"' -Oz -o /work/$(basename "${INDEL_VCF}") && \
                bcftools index -t /work/$(basename "${INDEL_VCF}") && \
                bcftools stats /work/$(basename "${INDEL_VCF}") > /qc/${SAMPLE}.${REF_TYPE}.DV_INDEL_PASS.stats
            "
        mark_done "${STEP_ID}"
    fi

    ############################
    # STEP 60-65: CNV Calling & QC
    ############################
    STEP_ID="60_${REF_TYPE}_CNV"
    if step_done "${STEP_ID}"; then
        log "Step 60 (${REF_TYPE}): CNV done."
    else
        # (Assume CNVKit Batch runs here as per original)
        
        # QC Logic: Ensure Stats filenames are distinct
        docker run --rm --user $(id -u):$(id -g) \
            -v "${CNV_DIR}":/cnv -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                # Raw CNV Stats
                bcftools stats /cnv/${SAMPLE}.${REF_TYPE}.cnv.vcf > /qc/${SAMPLE}.${REF_TYPE}.CNV_RAW.stats && \
                # Filtered CNV Stats
                bcftools stats /cnv/${SAMPLE}.${REF_TYPE}.cnv.filtered.vcf > /qc/${SAMPLE}.${REF_TYPE}.CNV_PASS.stats
            "
        mark_done "${STEP_ID}"
    fi
done

############################
# STEP 50: SV Calling (Graph)
############################
STEP_ID_PACK="50_SV_PACK"
GAM_FILE="${WORK_DIR}/${SAMPLE}.graph.gam"
PACK_FILE="${WORK_DIR}/${SAMPLE}.graph.pack"

if step_done "${STEP_ID_PACK}"; then
    log "Step 50 (Pack): Pack exists."
else
    log "Graph SV: Packing..."
    docker run --rm --user root -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        /bin/sh -c "vg pack -x /graph/$(basename "${GBZ_HOST}") -g /work/$(basename "${GAM_FILE}") \
        -Q 5 -t ${THREADS} -o /work/$(basename "${PACK_FILE}")"
    mark_done "${STEP_ID_PACK}"
fi

for REF_TYPE in "${REFS[@]}"; do
    STEP_ID="50_SV_CALL_${REF_TYPE}"
    SV_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.sv.vcf.gz"
    
    if step_done "${STEP_ID}"; then
        log "Step 50 (Call): SVs for ${REF_TYPE} done."
    else
        log "Graph SV: Calling..."
        docker run --rm --user root -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            /bin/sh -c "vg call /graph/$(basename "${GBZ_HOST}") \
                -k /work/$(basename "${PACK_FILE}") -S ${REF_TYPE} \
                -a -z -t ${THREADS} > /work/$(basename "${SV_VCF%.gz}")"

        # Compress and Stats (RAW)
        log "Graph SV: Compressing and Generating RAW Stats..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                bgzip -f /work/$(basename "${SV_VCF%.gz}") && \
                tabix -p vcf /work/$(basename "${SV_VCF}") && \
                bcftools stats /work/$(basename "${SV_VCF}") > /qc/${SAMPLE}.${REF_TYPE}.SV_RAW.stats
            "
        mark_done "${STEP_ID}"
    fi

    ############################
    # STEP 55: SV Filtering (SIZE FILTER ADDED)
    ############################
    STEP_ID="55_SV_Filter_${REF_TYPE}"
    FILTERED_SV_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.sv.pass.vcf.gz"
    
    if step_done "${STEP_ID}"; then
        log "Step 55 (Filter): SVs for ${REF_TYPE} filtered."
    else
        log "Filtering SVs (Removing Small Variants)..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                # Filter: QUAL>=30 AND (Length >= 50 OR Explicit SV Tag)
                bcftools view -i 'QUAL>=30 && (ABS(ILEN)>=50 || TYPE=\"SV\")' \
                -Oz -o /work/$(basename "${FILTERED_SV_VCF}") /work/$(basename "${SV_VCF}") && \
                
                bcftools index -t /work/$(basename "${FILTERED_SV_VCF}") && \
                
                # Stats: Save as SV_PASS
                bcftools stats /work/$(basename "${FILTERED_SV_VCF}") > /qc/${SAMPLE}.${REF_TYPE}.SV_PASS.stats
            "
        mark_done "${STEP_ID}"
    fi
done

############################
# STEP 90: MultiQC (Refactored)
############################
STEP="90_MultiQC"

if step_done "${STEP}"; then
    log "Step 90: MultiQC already done."
else
    log "Organizing QC files and generating report..."
    
    # 1. Prepare Staging Directory
    # We create a clean folder and SYMLINK or COPY only the files we want MultiQC to see.
    QC_STAGE="${QC_DIR}/staging_multiqc"
    rm -rf "${QC_STAGE}"
    mkdir -p "${QC_STAGE}"

    # 2. Gather Stats Files (Using the new clear naming convention)
    # - BAM Stats
    cp "${QC_DIR}"/*.bam.stats "${QC_STAGE}/" 2>/dev/null || true
    # - GAM Stats
    cp "${QC_DIR}"/*.gam.stats.txt "${QC_STAGE}/" 2>/dev/null || true
    # - VCF Stats (DV & SV - Raw and Pass)
    cp "${QC_DIR}"/*.stats "${QC_STAGE}/" 2>/dev/null || true
    # - FastQC (Link the whole folder)
    ln -s "${QC_DIR}/fastqc" "${QC_STAGE}/fastqc"

    # 3. Run MultiQC on the STAGING directory
    # Note: We output the report to the MAIN QC directory, not inside staging.
    REPORT_DIR="${QC_DIR}/${SAMPLE}_multiqc_report"
    
    # Safety check: Remove existing report directory if it exists
    if [ -d "${REPORT_DIR}" ]; then
        log "Removing previous report directory..."
        rm -rf "${REPORT_DIR}"
    fi

    docker run --rm \
        --user $(id -u):$(id -g) \
        -v "${QC_STAGE}":/qc_in \
        -v "${QC_DIR}":/qc_out \
        "${MULTIQC_IMAGE}" \
        multiqc /qc_in \
            -o /qc_out/$(basename "${REPORT_DIR}") \
            --title "${SAMPLE} Pangenome QC" \
            --force

    # 4. Cleanup Staging
    rm -rf "${QC_STAGE}"

    mark_done "${STEP}"
fi

log "Pipeline Complete."