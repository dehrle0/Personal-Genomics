#!/usr/bin/env bash
# vg_unified_v0.8.sh
#
# Unified Pangenome Pipeline v2.1
# Updates:
# - Step 50 now loops to call SVs for BOTH GRCh38 and CHM13.
# - Tracks progress in 'completed_steps_${SAMPLE}.txt'.
# - Fixes MultiQC to isolate current sample stats.

set -euo pipefail

# --- CONFIGURATION ---
SAMPLE="ME"              # Change to "ME" for next sample
SAMPLE_SEX="female"        # "male" or "female" (Critical for CNVkit)
THREADS="12"
REF_MODE="BOTH"          # "GRCh38", "CHM13", or "BOTH"

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
README_FILE="${WORK_DIR}/README_${SAMPLE}_v2.1.md"

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

log "Starting v2.1 for Sample: ${SAMPLE} (${SAMPLE_SEX})"

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

        # AWK-based Header Manipulation (No Python required)
        # We mount the FAI file to read lengths
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro "${SAMTOOLS_IMAGE}" /bin/bash -c "
            # 1. Dump header
            samtools view -H /work/$(basename "${TEMP_SORTED}") > /work/old_header.sam
            
            # 2. Process with AWK
            # We pass the reference type (GRCh38 or CHM13) and the FAI file path
            awk -v ref=\"${REF_TYPE}\" '
                BEGIN {
                    # Load FAI file into associative array: length[name] = len
                    while ((getline < \"/references/$(basename "${LOCAL_FAI}")\") > 0) {
                        split(\$0, a, \"\t\");
                        len[a[1]] = a[2];
                    }
                }
                {
                    if (\$1 == \"@SQ\") {
                        # Extract SN:name
                        for (i=2; i<=NF; i++) {
                            if (\$i ~ /^SN:/) {
                                orig_name = substr(\$i, 4);
                                # Clean the name
                                clean_name = orig_name;
                                gsub(ref \"#0#\", \"\", clean_name);
                                gsub(ref \"#\", \"\", clean_name);
                                
                                # Check if we have a length for this clean name
                                if (clean_name in len) {
                                    printf \"@SQ\tSN:%s\tLN:%s\n\", clean_name, len[clean_name];
                                } else {
                                    # Fallback: just print the line with the name cleaned but old length (rare)
                                    gsub(orig_name, clean_name);
                                    print \$0;
                                }
                                next;
                            }
                        }
                    }
                    # Print all other lines (HD, PG, CO) as is
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
            -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "
                bcftools view -v snps -f PASS -i 'QUAL>=30 && DP>=10' /work/$(basename "${FINAL_VCF}") -Oz -o /work/$(basename "${SNP_VCF}") && \
                bcftools index -t /work/$(basename "${SNP_VCF}") && \
                bcftools norm -f ${DOCKER_REF_PATH} -m -any /work/$(basename "${FINAL_VCF}") -Ou | \
                bcftools view -V snps -i 'QUAL>=15 && DP>=10 && FILTER!=\"RefCall\"' -Oz -o /work/$(basename "${INDEL_VCF}") && \
                bcftools index -t /work/$(basename "${INDEL_VCF}") && \
                bcftools stats /work/$(basename "${SNP_VCF}") > /work/$(basename "${SNP_VCF}").stats && \
                bcftools stats /work/$(basename "${INDEL_VCF}") > /work/$(basename "${INDEL_VCF}").stats
            "
        mv "${WORK_DIR}"/*.stats "${QC_DIR}/"
        mark_done "${STEP_ID}"
    fi

    # STEP 60: CNV Calling
    STEP_ID="60_${REF_TYPE}"
    if step_done "${STEP_ID}"; then
        log "Step 60 (${REF_TYPE}): CNV done."
    else
        log "Running CNVkit for ${REF_TYPE}..."
        
        # CORRECTED SEX LOGIC:
        # CNVkit batch assumes female by default. 
        # We only add '--male' if the sample is male.
        SEX_FLAG=""
        if [[ "${SAMPLE_SEX}" == "male" ]]; then SEX_FLAG="--male"; fi
        
        docker run --rm --user $(id -u):$(id -g) \
            -v "${WORK_DIR}":/work -v "${REF_DIR_HOST}":/references:ro -v "${CNV_DIR}":/output "${CNVKIT_IMAGE}" \
            cnvkit.py batch "/work/$(basename "${FINAL_BAM}")" \
            --normal --fasta "${DOCKER_REF_PATH}" --method wgs \
            --output-dir /output/ ${SEX_FLAG} --drop-low-coverage --scatter --diagram
        
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
        # We use -S ${REF_TYPE} because HPRC graphs contain paths like GRCh38#... and CHM13#...
        docker run --rm --user root -v "${WORK_DIR}":/work -v "${GRAPH_DIR_HOST}":/graph "${VG_IMAGE}" \
            /bin/sh -c "vg call /graph/$(basename "${GBZ_HOST}") \
                -k /work/$(basename "${PACK_FILE}") \
                -S ${REF_TYPE} \
                -a -z -t ${THREADS} > /work/$(basename "${SV_VCF%.gz}")"

        # Compress/Index
        docker run --rm --user $(id -u):$(id -g) -v "${WORK_DIR}":/work "${BCFTOOLS_IMAGE}" \
            /bin/sh -c "bgzip -f /work/$(basename "${SV_VCF%.gz}") && tabix -p vcf /work/$(basename "${SV_VCF}")"
            
        mark_done "${STEP_ID}"
    fi
done
# rm -f "${PACK_FILE}" # Optional cleanup

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

    # 1. Stats files (BAM, VCF, GAM)
    find "${QC_DIR}" -maxdepth 1 -type f \
        \( -name "${SAMPLE}*.stats" -o -name "${SAMPLE}*.txt" \) \
        -exec cp {} "${QC_DIR}/current_run/" \;

    # 2. FastQC HTML + ZIP files
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
    # Cleanup unless debugging
    ###############################################

    if [[ -z "${DEBUG_QC}" ]]; then
        rm -rf "${QC_DIR}/current_run"
    else
        log "DEBUG_QC enabled — preserving current_run/"
    fi

    mark_done "${STEP}"
fi


log "Pipeline Complete."
