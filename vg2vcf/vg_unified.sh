#!/usr/bin/env bash
# vg_unified.sh
#
# Unified Pangenome Pipeline
# Updates:
# - Step 50: Added BCFTools stats generation for SV VCFs.
# - Step 55: Added SV filtering based on QUAL >= 30.
# - Step 60: Added VCF export and stats generation for CNVkit results.
# - Step 65: Added CNV filtering based on FILTER=PASS.
# - Step 90: Updated MultiQC file gathering to include filtered stats.

set -euo pipefail

# --- CONFIGURATION ---
SAMPLE="ME"              # Change to "ME" for next sample
SAMPLE_SEX="female"      # "male" or "female" (Critical for CNVkit)
THREADS="12"
REF_MODE="GRCh38"        # "GRCh38", "CHM13", or "BOTH"

# --- HOST PATHS ---
REF_DIR_HOST="/mnt/data/work/references"
HG38_LOCAL="${REF_DIR_HOST}/hg38.fa"
HS1_LOCAL="${REF_DIR_HOST}/chm13v2.0_maskedY_rCRS.fa"

GRAPH_DIR_HOST="${REF_DIR_HOST}/hrpc"
GBZ_HOST="${GRAPH_DIR_HOST}/hprc-v1.1-mc-chm13.d9.gbz"

# Reads (Assumes /mnt/c/Genomes/SAMPLE/Data/Source structure)
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

log "Starting v2.2 for Sample: ${SAMPLE} (${SAMPLE_SEX})"

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
        log "Surjecting to ${REF_TYPE}..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${GRAPH_DIR_HOST}":/graph -v "${WORK_DIR}":/work "${VG_IMAGE}" \
            vg surject -x "/graph/$(basename "${GBZ_HOST}")" -b -t "${THREADS}" --into-ref "${REF_TYPE}" \
            "/work/$(basename "${GAM_FILE}")" > "${RAW_BAM}"

        log "Sorting and Fixing Header..."
        TEMP_SORTED="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.sorted.bam"
        
        # Dockerized Samtools Sort
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${SAMTOOLS_IMAGE}" \
            samtools sort -@ "${THREADS}" -o "/work/$(basename "${TEMP_SORTED}")" "/work/$(basename "${RAW_BAM}")"

        # AWK-based Header Manipulation
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro "${SAMTOOLS_IMAGE}" /bin/bash -c "
            samtools view -H /work/$(basename "${TEMP_SORTED}") > /work/old_header.sam
            
            awk -v ref=\"${REF_TYPE}\" '
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
                                    gsub(orig_name, clean_name);
                                    print \$0;
                                }
                                next;
                            }
                        }
                    }
                    print \$0;
                }
            ' /work/old_header.sam > /work/new_header.sam
        "

        # Dockerized Reheader, Index, and Stats
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${SAMTOOLS_IMAGE}" /bin/bash -c "
            samtools reheader /work/new_header.sam /work/$(basename "${TEMP_SORTED}") > /work/$(basename "${FINAL_BAM}") && \
            samtools index /work/$(basename "${FINAL_BAM}") && \
            samtools stats -@ ${THREADS} /work/$(basename "${FINAL_BAM}") > /qc/${SAMPLE}.${REF_TYPE}.bam.stats"

        # Cleanup
        rm -f "${WORK_DIR}/new_header.sam" "${WORK_DIR}/old_header.sam" "${RAW_BAM}" "${TEMP_SORTED}"
        mark_done "${STEP_ID}"
    fi
 
    ############################
    # STEP 40: DeepVariant
    ############################
    STEP_ID="40_${REF_TYPE}"
    FINAL_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.dv.vcf.gz"
    if step_done "${STEP_ID}"; then
        log "Step 40 (${REF_TYPE}): DeepVariant done."
    else
        log "Running DeepVariant for ${REF_TYPE}..."
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

    # STEP 45: Specialized Filtering
    STEP_ID="45_${REF_TYPE}"
    SNP_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.snps.vcf.gz"
    INDEL_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.indels_complex.vcf.gz"
    if step_done "${STEP_ID}"; then
        log "Step 45 (${REF_TYPE}): Filtering done."
    else
        log "Splitting VCFs for ${REF_TYPE}..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                bcftools view -v snps -f PASS -i 'QUAL>=30 && DP>=10' /work/$(basename "${FINAL_VCF}") -Oz -o /work/$(basename "${SNP_VCF}") && \
                bcftools index -t /work/$(basename "${SNP_VCF}") && \
                bcftools norm -f ${DOCKER_REF_PATH} -m -any /work/$(basename "${FINAL_VCF}") -Ou | \
                bcftools view -V snps -i 'QUAL>=15 && DP>=10 && FILTER!=\"RefCall\"' -Oz -o /work/$(basename "${INDEL_VCF}") && \
                bcftools index -t /work/$(basename "${INDEL_VCF}") && \
                bcftools stats /work/$(basename "${SNP_VCF}") > /qc/$(basename "${SNP_VCF}").stats && \
                bcftools stats /work/$(basename "${INDEL_VCF}") > /qc/$(basename "${INDEL_VCF}").stats
            "
        mark_done "${STEP_ID}"
    fi

    ############################
    # STEP 60: CNV Calling (UPDATED)
    ############################
    STEP_ID="60_${REF_TYPE}"
    if step_done "${STEP_ID}"; then
        log "Step 60 (${REF_TYPE}): CNV done."
    else
        log "Running CNVkit for ${REF_TYPE}..."
        
        SEX_FLAG=""
        if [[ "${SAMPLE_SEX}" == "male" ]]; then SEX_FLAG="--male"; fi
        
        # 1. CNVkit Batch
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro -v "${CNV_DIR}":/output "${CNVKIT_IMAGE}" \
            cnvkit.py batch "/work/$(basename "${FINAL_BAM}")" \
            --normal --fasta "${DOCKER_REF_PATH}" --method wgs \
            --output-dir /output/ ${SEX_FLAG} --drop-low-coverage --scatter --diagram

        # 2. Export to VCF for Stats
        log "Exporting CNV calls to VCF..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${CNV_DIR}":/output "${CNVKIT_IMAGE}" \
            cnvkit.py export vcf "/output/${SAMPLE}.${REF_TYPE}.cns" \
            -i "${SAMPLE}" -o "/output/${SAMPLE}.${REF_TYPE}.cnv.vcf"

        # 3. Generate Stats
        log "Generating CNV stats..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${CNV_DIR}":/cnv -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "bcftools stats /cnv/${SAMPLE}.${REF_TYPE}.cnv.vcf > /cnv/${SAMPLE}.${REF_TYPE}.cnv.stats"
        
        mark_done "${STEP_ID}"
    fi

    ############################
    # STEP 65: CNV Filtering
    ############################
    STEP_ID="65_CNV_Filter_${REF_TYPE}"
    FILTERED_CNV_VCF="${CNV_DIR}/${SAMPLE}.${REF_TYPE}.cnv.filtered.vcf"
    if step_done "${STEP_ID}"; then
        log "Step 65 (Filter): CNVs for ${REF_TYPE} already filtered."
    else
        log "Filtering CNVs for ${REF_TYPE}..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${CNV_DIR}":/cnv "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                bcftools view -f PASS -Ov -o /cnv/$(basename "${FILTERED_CNV_VCF}") /cnv/${SAMPLE}.${REF_TYPE}.cnv.vcf && \
                bcftools stats /cnv/$(basename "${FILTERED_CNV_VCF}") > /cnv/${SAMPLE}.${REF_TYPE}.cnv.filtered.stats
            "
        mark_done "${STEP_ID}"
    fi
done

############################
# STEP 50: SV Calling (Multi-Ref: GRCh38 & CHM13)
############################
STEP_ID_PACK="50_SV_PACK"
GAM_FILE="${WORK_DIR}/${SAMPLE}.graph.gam"
PACK_FILE="${WORK_DIR}/${SAMPLE}.graph.pack"

# 1. Generate Pack (Once)
if step_done "${STEP_ID_PACK}"; then
    log "Step 50 (Pack): Pack file exists."
else
    log "Graph SV: Packing GAM..."
    docker run --rm --user root -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
        /bin/sh -c "vg pack -x /graph/$(basename "${GBZ_HOST}") -g /work/$(basename "${GAM_FILE}") \
        -Q 5 -t ${THREADS} -o /work/$(basename "${PACK_FILE}")"
    mark_done "${STEP_ID_PACK}"
fi

# 2. Call Loop (GRCh38 & CHM13)
for REF_TYPE in "${REFS[@]}"; do
    STEP_ID="50_SV_CALL_${REF_TYPE}"
    SV_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.sv.vcf.gz"
    
    if step_done "${STEP_ID}"; then
        log "Step 50 (Call): SVs for ${REF_TYPE} already done."
    else
        log "Graph SV: Calling variants for ${REF_TYPE}..."
        docker run --rm --user root -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            /bin/sh -c "vg call /graph/$(basename "${GBZ_HOST}") \
                -k /work/$(basename "${PACK_FILE}") \
                -S ${REF_TYPE} \
                -a -z -t ${THREADS} > /work/$(basename "${SV_VCF%.gz}")"

        # Compress/Index AND Generate Stats (UPDATED)
        log "Graph SV: Compressing and generating stats for ${REF_TYPE}..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                bgzip -f /work/$(basename "${SV_VCF%.gz}") && \
                tabix -p vcf /work/$(basename "${SV_VCF}") && \
                bcftools stats /work/$(basename "${SV_VCF}") > /qc/${SAMPLE}.${REF_TYPE}.sv.vcf.stats
            "
            
        mark_done "${STEP_ID}"
    fi

    ############################
    # STEP 55: SV Filtering
    ############################
    STEP_ID="55_SV_Filter_${REF_TYPE}"
    FILTERED_SV_VCF="${WORK_DIR}/${SAMPLE}.${REF_TYPE}.sv.filtered.vcf.gz"
    if step_done "${STEP_ID}"; then
        log "Step 55 (Filter): SVs for ${REF_TYPE} already filtered."
    else
        log "Filtering SVs for ${REF_TYPE}..."
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${QC_DIR}":/qc "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                bcftools view -i 'QUAL>=30' -Oz -o /work/$(basename "${FILTERED_SV_VCF}") /work/$(basename "${SV_VCF}") && \
                bcftools index -t /work/$(basename "${FILTERED_SV_VCF}") && \
                bcftools stats /work/$(basename "${FILTERED_SV_VCF}") > /qc/${SAMPLE}.${REF_TYPE}.sv.filtered.vcf.stats
            "
        mark_done "${STEP_ID}"
    fi
done

############################
# STEP 90: MultiQC
############################

STEP="90_MultiQC"

if step_done "${STEP}"; then
    log "Step 90: MultiQC already done."
else
    log "Generating MultiQC..."

    # Fresh workspace
    mkdir -p "${QC_DIR}/current_run"
    rm -f "${QC_DIR}/current_run/"*

    ###############################################
    # Copy QC artifacts explicitly (no directories)
    ###############################################

    # 1. Main QC Stats (BAM, VCF, GAM, SV)
    # Catches .stats, .txt, and .vchk files in the QC directory
    find "${QC_DIR}" -maxdepth 1 -type f \
        \( -name "${SAMPLE}*.stats" -o -name "${SAMPLE}*.txt" -o -name "${SAMPLE}*.vchk" \) \
        -exec cp {} "${QC_DIR}/current_run/" \;

    # 2. CNV Stats (located in CNV_DIR)
    # Catches the specific CNV VCF stats we generated in Step 60
    if [ -d "${CNV_DIR}" ]; then
        find "${CNV_DIR}" -maxdepth 1 -type f \
            -name "${SAMPLE}*.stats" \
            -exec cp {} "${QC_DIR}/current_run/" \;
    fi

    # 3. FastQC HTML + ZIP files
    if [ -d "${QC_DIR}/fastqc" ]; then
        find "${QC_DIR}/fastqc" -maxdepth 1 -type f \
            \( -name "*.html" -o -name "*.zip" \) \
            -exec cp {} "${QC_DIR}/current_run/" \;
    fi

    ###############################################
    # Run MultiQC
    ###############################################

    docker run --rm \
        --user $(id -u):$(id -g) \
        -v "${QC_DIR}/current_run":/qc \
        "${MULTIQC_IMAGE}" \
        multiqc /qc \
            -o /qc/multiqc_report \
            --title "${SAMPLE} Pangenome QC" \
            --force

    ###############################################
    # Move final report into stable location
    ###############################################

    mv "${QC_DIR}/current_run/multiqc_report" \
       "${QC_DIR}/${SAMPLE}_multiqc_report"

    ###############################################
    # Cleanup current_run
    ###############################################

    rm -rf "${QC_DIR}/current_run"

    mark_done "${STEP}"
fi

log "Pipeline Complete."