# Pangenome WGS Unified Pipeline
**Results only spot checked - validation needed**

A robust, "batteries-included" Bash pipeline for Whole Genome Sequencing (WGS) analysis using the **HPRC Pangenome Graph** (v1.1).

This pipeline performs read mapping against a pangenome graph (VG Giraffe) and projects the data onto **two linear references simultaneously** (GRCh38 and CHM13). It outputs fully compliant BAMs, VCFs (SNPs/Indels), Structural Variants (SVs), and Copy Number Variants (CNVs).

## 🚀 Key Features
* **Graph-Based Mapping:** Uses `vg giraffe` to map against the complex HPRC graph.
* **Dual-Reference Output:** Automatically surjects reads to both **GRCh38** and **CHM13** (T2T) coordinate systems.
* **State-of-the-Art Callers:**
    * **Small Variants:** Google DeepVariant (WGS model).
    * **Structural Variants:** `vg call` (Graph-based SV detection).
    * **CNVs:** CNVkit (Moving average approach).
* **Resumable:** Tracks progress in a `.txt` file; automatically skips completed steps.
* **Dockerized:** All tools run in isolated containers (VG, DeepVariant, BCFTools, CNVkit, MultiQC).

## 🛠️ Pipeline Steps

### **Phase 1: Pre-processing**
* **Step 05: Reads QC**
    * Runs **FastQC** on raw `.fq.gz` inputs to assess sequencing quality.
* **Step 10: Mapping (Pangenome)**
    * Maps reads to the `.gbz` graph (HPRC v1.1) using `vg giraffe`.
    * **Output:** `.gam` (Graph Alignment Map).
* **Step 15: Graph QC**
    * Generates alignment statistics (coverage, mapping quality) from the GAM file.

### **Phase 2: Linear Reference Processing (Loop)**
*Performs the following for BOTH **GRCh38** and **CHM13**:*

* **Step 30: Surjection & Header Fix**
    * **Surjection:** Projects the graph alignment (.gam) onto the linear reference path.
    * **Sanitization:** Uses `awk` (no Python dependency) to clean non-standard contig names (e.g., `GRCh38#0#chr1` → `chr1`) and fix sequence lengths in the SAM header to match the specific reference assembly.
    * **Output:** Sorted, indexed, "DeepVariant-ready" BAM.
* **Step 40: Small Variant Calling**
    * Runs **Google DeepVariant** (WGS model) to call SNPs and small Indels.
* **Step 45: Specialized Filtering**
    * Splits the DeepVariant VCF into:
        1.  **High-Confidence SNPs:** `QUAL>=30 && DP>=10`.
        2.  **Indels:** Normalized and filtered (`QUAL>=15`).
* **Step 60: CNV Calling**
    * Runs **CNVkit** in `wgs` mode.
    * **Auto-Sex Detection:** Automatically configures flags (`--male` or default female) based on user config.
* **Step 65: CNV QC and Filtering**
    * Filters CNV VCFs to PASS calls and generates stats.

### **Phase 3: Structural Variants & Reporting**
* **Step 50: SV Calling**
    * Packs the GAM alignment and calls Structural Variants directly from the graph topology using `vg call`.
* **Step 55: SV QC and Filtering**
    * Filters SV VCFs to high-confidence calls (QUAL >= 30) and generates stats.
* **Step 90: MultiQC**
    * Aggregates FastQC, Samtools, BCFTools, and CNVkit reports into a single interactive HTML summary.

## 📋 Requirements
* **OS:** Linux (Ubuntu/Debian recommended or WSL2).
* **Docker:** Must be installed and running.
* **Hardware:** Recommended 12+ Threads, 48GB+ RAM.

### Docker Images Used
| Tool | Image |
| :--- | :--- |
| **vg** | `quay.io/vgteam/vg:v1.70.0` |
| **DeepVariant** | `google/deepvariant:1.10.0-beta` |
| **Samtools** | `staphb/samtools:1.19` |
| **BCFtools** | `quay.io/biocontainers/bcftools:1.19` |
| **CNVkit** | `etal/cnvkit:0.9.10` |
| **MultiQC** | `ewels/multiqc:v1.19` |

## 💡 Quick Tip: Resuming & Re-running
The pipeline creates a file named `completed_steps_SAMPLE.txt` in your working directory.
* **To Resume:** Just run the script again. It skips steps listed in this file.
* **To Re-run a Step:** Delete the specific line from the text file.
    * *Example:* To re-do CNV calling for GRCh38, delete the line `60_GRCh38` and re-run the script.

## 📜 Usage
1.  Edit the **CONFIGURATION** section at the top of `vg_unified.sh` (Sample ID, Sex, Paths).
2.  Run the script:
    ```bash
    bash vg_unified_v1.0.sh
    ```
---

## Step Identifiers & Execution Control

The pipeline utilizes a checkpoint system stored in `completed_steps_${SAMPLE}.txt`. To **exclude** a step from a run or skip steps that have already finished, ensure the corresponding **Step ID** listed below is present in that file.

### 1. Global Processing Steps

These steps run once per sample, regardless of how many reference genomes (GRCh38/CHM13) are being targeted.

| Step ID | Description | Tool |
| --- | --- | --- |
| `05_FastQC` | Generates quality control reports for raw FASTQ reads. | FastQC |
| `10_Mapping` | Align raw reads to the pangenome graph (Output: GAM). | vg giraffe |
| `15_GAM_QC` | Calculates alignment statistics for the graph mapping. | vg stats |
| `50_SV_PACK` | Generates the coverage "pack" file required for SV calling. | vg pack |
| `90_MultiQC` | Aggregates all stats and logs into a final HTML report. | MultiQC |

### 2. Reference-Specific Steps

These steps run inside the reference loop. Replace `{REF}` with `GRCh38` or `CHM13` (e.g., `40_GRCh38`).

| Step ID | Description | Tool |
| --- | --- | --- |
| `30_{REF}` | Surjects graph alignments to linear BAM and fixes headers. | vg surject |
| `40_{REF}` | Performs small variant calling (SNPs and Indels). | DeepVariant |
| `45_{REF}` | Splits and filters VCFs into PASS-only SNPs and Indels. | bcftools |
| `50_SV_CALL_{REF}` | Identifies variants (including SVs) from the graph. | vg call |
| `55_SV_Filter_{REF}` | Filters structural variants by size ( 50bp). | bcftools |
| `60_{REF}` | Identifies Copy Number Variants and exports to VCF. | CNVkit |
| `65_CNV_Filter_{REF}` | Filters CNV calls for `PASS` status. | bcftools |

---

### How to skip a step manually

If a step fails or you wish to bypass a specific module, append the ID to your completion file:

```bash
echo "60_GRCh38" >> completed_steps_ME.txt

```
flowchart TD
    %% Global Inputs
    Reads[("Raw Reads (FASTQ)")]
    Graph[("Pangenome Graph (GBZ)")]

    %% Pre-Processing
    subgraph Global_Steps ["Global Processing"]
        direction TB
        Step05("05_FastQC")
        Step10("10_Mapping")
        Step15("15_GAM_QC")
        Step50P("50_SV_PACK")
    end

    Reads --> Step05
    Reads --> Step10
    Graph --> Step10
    Step10 --> Step15
    Step10 -- "GAM File" --> Step50P

    %% Loop Logic
    subgraph Ref_Loop ["Reference Loop (GRCh38 / CHM13)"]
        direction TB
        
        %% Linear Path
        Step30("30_{REF} (Surjection)")
        
        %% DeepVariant Path
        Step40("40_{REF} (DeepVariant)")
        Step45("45_{REF} (Filtering)")
        
        %% CNV Path
        Step60("60_{REF} (CNVkit)")
        Step65("65_CNV_Filter_{REF}")

        %% SV Path
        Step50C("50_SV_CALL_{REF}")
        Step55("55_SV_Filter_{REF}")

        %% Connections
        Step30 -- "BAM File" --> Step40
        Step30 -- "BAM File" --> Step60
        
        Step40 --> Step45
        Step60 --> Step65
        
        Step50P -- "Pack File" --> Step50C
        Step50C --> Step55
    end

    %% Final Output
    Step90("90_MultiQC")

    %% Aggregation
    Step10 -- "GAM" --> Step30
    Step05 & Step15 --> Step90
    Step45 & Step65 & Step55 --> Step90

    %% Styling
    classDef core fill:#f9f,stroke:#333,stroke-width:2px;
    classDef qc fill:#bbf,stroke:#333,stroke-width:1px;
    class Step10,Step30,Step40,Step50P,Step50C,Step60 core;
    class Step05,Step15,Step90,Step45,Step55,Step65 qc;
