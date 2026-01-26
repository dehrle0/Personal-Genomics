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
    bash vg_unified.sh
    ```
